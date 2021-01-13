// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "assert_macro.h"
#include "stream_func.h"
#include "movable.h"
#include "exceptions.h"
#include "quaternion.h"
#include "iowrapper.h"
#include "tokenizer.h"
#include "glossary.h"
#include "modulo.h"
#include "random.h"
#include "space.h"
#include "simul.h"

extern Random RNG;

/** The default implementation is invalid */
void Movable::translate(Vector const&)
{
    ABORT_NOW("Movable::translate() called for immobile object");
}

/**
The default implementation:
@code
 translate( w - position() );
@endcode
can be redefined in derived class for efficiency.
*/
void Movable::setPosition(Vector const& vec)
{
    assert_true( mobile() );
    translate( vec - position() );
}

/*
 if mobile() is true,
 the Object is translated by `[ rot * Object::position() - Object::position() ]`
*/
void Movable::rotate(Rotation const& rot)
{
    if ( mobile() )
    {
        Vector pos = position();
        translate(rot*pos-pos);
    }
}

/**
The default implementation:
@code
 Vector G = position();
 translate( -G );
 rotate( T );
 translate(  G ); 
@endcode
can be redefined in derived class for efficiency.
*/
void Movable::revolve(Rotation const& T)
{
    Vector G = position();
    translate( -G );
    rotate( T );
    translate(  G );
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 There are different ways to specify a position:
 
 keyword & parameters |   Position (X, Y, Z)
 ---------------------|---------------------------------------------------------
 `A B C`              | The specified vector (A,B,C)
 `inside`             | A random position inside the current Space
 `edge E`             | At distance E from the edge of the current Space
 `surface E`          | On the surface of the current Space\n By projecting a point at distance E from the surface.
 `line L T`           | Selected randomly with -L/2 < X < L/2; norm(Y,Z) < T
 `sphere R T`         | At distance R +/- T/2 from the origin\n `R-T/2 < norm(X,Y,Z) < R+T/2`
 `ball R`             | At distance R at most from the origin\n `norm(X,Y,Z) < R`
 `disc R T`           | in 2D, a disc in the XY-plane \n in 3D, a disc in the XY-plane of thickness T in Z
 `discXZ R T`         | Disc in the XZ-plane of radius R, thickness T
 `discYZ R T`         | Disc in the YZ-plane of radius R, thickness T
 `equator R T`        | At distance R from the origin, and T from the XY plane:\n `norm(X,Y) < R` `norm(Z) < T`
 `circle R T`         | Circle of radius R and thickness T \n At distance T from the circle of radius R
 `cylinder W R`       | Cylinder of axis X, W=thickness in X, R=radius in YZ
 `ellipse A B C`      | Inside the ellipse or ellipsoid of main axes 2A, 2B and 2C
 `arc L Theta`        | A piece of circle of length L and covering an angle Theta
 `stripe L R`         | Random vector with L < X < R
 `square R`           | Random vector with -R < X < R; -R < Y < R; -R < Z < R;
 `rectangle A B C`    | Random vector with -A < X < A; -B < Y < B; -C < Z < C;
 `gradient S E R`     | Linear density gradient, from 0 at X=S to 1 at X=E, within a cylinder of radius R

 Each primitive describes a certain area in Space, and in most cases the returned position is
 chosen randomly inside this area following a uniform probability.
 */

Vector Movable::readPrimitive(std::istream& is, const Space* spc)
{
    char c = Tokenizer::get_space(is, false);

    if ( c == EOF )
        return Vector(0,0,0);

    if ( isalpha(c) )
    {
        std::string tok = Tokenizer::get_symbol(is);
        
        if ( spc )
        {
            if ( tok == "inside" || tok == "random" )
                return spc->randomPlace();
            
            if ( tok == "edge" )
            {
                real R = 0;
                is >> R;
                if ( R < REAL_EPSILON )
                    throw InvalidParameter("you must specify a distance R > 0 in `edge R`");
                return spc->randomPlaceNearEdge(R);
            }
            
            if ( tok == "surface" )
            {
                real e = 0.1;
                is >> e;
                return spc->randomPlaceOnEdge(e);
            }

            if ( tok == "outside_sphere" )
            {
                real R = 0;
                is >> R;
                if ( R < 0 )
                    throw InvalidParameter("you must specify a radius R >= 0 in `outside_sphere R`");
                Vector P;
                do
                    P = spc->randomPlace();
                while ( P.norm() < R );
                return P;
            }
            
            if (tok == "outside_ZYring")
            {
                real R = 0,x = 0;
                is>> R >> x;
                if ( R < 0 )
                    throw InvalidParameter("you must specify a radius R >= 0 in `outside_ZYring x R`");
                Vector P;
                do
                {
                    P = spc->randomPlace();
                    P.XX = 0;
                }
                while ( P.norm() < R );
                P.XX = x;
                return P;
            }
            
            if ( tok == "stripe" )
            {
                real s = -0.5, e = 0.5;
                is >> s >> e;
                Vector inf, sup;
                spc->boundaries(inf, sup);
                Vector pos = inf + (sup-inf).e_mul(Vector::prand());
                pos.XX = RNG.real_uniform(s, e);
                return pos;
            }
        }
        
        if ( tok == "sphere" )
        {
            real R = -1, T = 0;
            is >> R >> T;
            if ( R < 0 )
                throw InvalidParameter("you must specify a radius R >= 0 in `sphere R`");
            if ( T < 0 )
                throw InvalidParameter("the thickness T must be >= 0 in `sphere R T`");
            return Vector::randU(R) + Vector::randU(T*0.5);
        }
        
        if ( tok == "equator" )
        {
            real R = 0, T = 0;
            is >> R >> T;
            if ( R < 0 )
                throw InvalidParameter("you must specify a radius R >= 0 in `equator R T`");
            if ( T < 0 )
                throw InvalidParameter("the thickness T must be >= 0 in `equator R T`");
            Vector2 vec2 = Vector2::randU();
            return Vector(R*vec2.XX, R*vec2.YY, T*RNG.sreal_half());
        }
       
        if ( tok == "cylinder" )
        {
            real L = -1, R = -1;
            is >> L >> R;
            if ( L < 0 )
                throw InvalidParameter("you must specify a length L >= 0 in `cylinder L R`");
            if ( R < 0 )
                throw InvalidParameter("you must specify a radius R >= 0 in `cylinder L R`");
            Vector2 YZ = Vector2::randB(R);
            return Vector(L*RNG.sreal_half(), YZ.XX, YZ.YY);
        }
        
        if ( tok == "circle" )
        {
            real R = -1, T = 0;
            is >> R >> T;
            if ( R < 0 )
                throw InvalidParameter("you must specify a radius R >= 0 in `circle R T`");
            if ( T < 0 )
                throw InvalidParameter("the thickness T must be >= 0 in `circle R T`");
#if ( DIM == 3 )
            Vector2 XY = Vector2::randU(R);
            return Vector3(XY.XX, XY.YY, 0) + (0.5*T) * Vector3::randU();
#endif
            return Vector::randU(R) + Vector::randU(T*0.5);
        }
        
        if ( tok == "ball" )
        {
            real R = -1;
            is >> R;
            if ( R < 0 )
                throw InvalidParameter("you must specify a radius R >= 0 in `ball R`");
            return Vector::randB(R);
        }
        
        if ( tok == "disc" || tok == "discXY" )
        {
            real R = -1, T = 0;
            is >> R >> T;
            if ( R < 0 )
                throw InvalidParameter("you must specify a radius R >= 0 in `disc R`");
#if ( DIM == 3 )
            //in 3D, a disc in the XY-plane of thickness T in Z-direction
            if ( T < 0 )
                throw InvalidParameter("the thickness T must be >= 0 in `disc R T`");
            Vector2 V = Vector2::randB(R);
            return Vector(V.XX, V.YY, T*RNG.sreal_half());
#endif
            //in 2D, a disc in the XY-plane
            return Vector::randB(R);
        }
        
        if ( tok == "discXZ"  )
        {
            real R = -1, T = 0;
            is >> R >> T;
            if ( R < 0 )
                throw InvalidParameter("you must specify a radius R >= 0 in `discXZ R`");
            if ( T < 0 )
                throw InvalidParameter("the thickness T must be >= 0 in `discXZ R T`");
            Vector2 V = Vector2::randB(R);
            return Vector(V.XX, T*RNG.sreal_half(), V.YY);
        }
        
        if ( tok == "discYZ"  )
        {
            real R = -1, T = 0;
            is >> R >> T;
            if ( R < 0 )
                throw InvalidParameter("you must specify a radius R >= 0 in `discYZ R`");
            if ( T < 0 )
                throw InvalidParameter("the thickness T must be >= 0 in `discYZ R T`");
            Vector2 V = Vector2::randB(R);
            return Vector(T*RNG.sreal_half(), V.XX, V.YY);
        }
        
        if ( tok == "ellipse" )
        {
            real x = 1, y = 1, z = 0;
            is >> x >> y >> z;
            return Vector(x,y,z).e_mul(Vector::randB());
        }
        
        if ( tok == "ellipse_surface" )
        {
            real x = 1, y = 1, z = 0;
            is >> x >> y >> z;
            return Vector(x,y,z).e_mul(Vector::randU());
        }
        
        if ( tok == "line" )
        {
            real L = -1, T = 0;
            is >> L >> T;
            if ( L < 0 )
                throw InvalidParameter("you must specify a length L >= 0 in `line L`");
            if ( T < 0 )
                throw InvalidParameter("the thickness T must be >= 0 in `line L T`");
#if ( DIM == 3 )
            Vector2 V = Vector2::randB(T);
            return Vector(L*RNG.sreal_half(), V.XX, V.YY);
#endif
            return Vector(L*RNG.sreal_half(), T*RNG.sreal_half(), 0);
        }
        
        if ( tok == "arc" )
        {
            real L = -1, A = 1.57;
            is >> L >> A;
            
            if ( L <= 0 )
                throw InvalidParameter("you must specify a length L >= 0 in `arc L`");
            
            real x = 0, y = 0;
            if ( A == 0 ) {
                x = 0;
                y = L * RNG.sreal_half();
            }
            else {
                real R = L / A;
                real angle = A * RNG.sreal_half();
                x = R * cos(angle) - R; // origin centered on arc
                y = R * sin(angle);
            }
#if ( DIM < 3 )
            return Vector(x, y, 0);
#else
            // in 3D this distributes on a portion of the sphere, with a gradient
            real a = RNG.sreal() * M_PI;
            return Vector(x, y*cos(a), y*sin(a));
#endif
        }
        
        if ( tok == "square" )
        {
            real x = 1;
            is >> x;
            return Vector::srand(x);
        }
        
        if ( tok == "rectangle" )
        {
            real x = 0, y = 0, z = 0;
            is >> x >> y >> z;
            return Vector(x,y,z).e_mul(Vector::srand());
        }
        
        if ( tok == "gradient" )
        {
            real s = -10, e = 10, R = 10;
            is >> s >> e >> R;
            real x;
            do {
                x = RNG.preal();
            } while ( RNG.preal() > x );
            return Vector(s+(e-s)*x, R*RNG.sreal(), R*RNG.sreal());
        }
        
#if ( 1 )
        /// A contribution from Beat Rupp
        if ( tok == "segment" || tok == "newsegment" )
        {
            real bending = 0, length = 0, thickness = 0, rotation = 0;
            is >> bending >> length >> thickness >> rotation;
            real x=0, y=0;
            
            // straight
            if ( bending == 0 ) {
                x = thickness * RNG.sreal_half();
                y = length * RNG.preal();
            } else {
                real radius = length / (bending * M_PI);
                real radiusInner = radius - thickness/2.0;
                real theta = fabs( length / radius );
                real angle = RNG.preal() * theta;
                // substract R to have the arc start from 0,0:
                x = (radiusInner + thickness * RNG.preal()) * cos( angle ) - radius;
                y = (radiusInner + thickness * RNG.preal()) * sin( angle );            
            }
            
            real c = cos(rotation);
            real s = sin(rotation);        
            
            // rotate
            return Vector(c*x + s*y , -s*x + c*y, 0 );        
        }
#endif
        if ( tok == "center" || tok == "origin" )
            return Vector(0,0,0);
            
        throw InvalidParameter("Unknown position `"+tok+"'");
    }
    
    // expect a vector to be specified:
    real x = 0, y = 0, z = 0;
    is >> x >> y >> z;
    return Vector(x,y,z);
}


//------------------------------------------------------------------------------
/**
 A position is defined with a SHAPE followed by a number of TRANSFORMATION.
 
 TRANSFORMATION         |  Result
 -----------------------|-------------------------------------------------------
 `at X Y Z`             | Translate by specified vector (X,Y,Z)
 `add SHAPE`            | Translate by a vector chosen according to SHAPE
 `align VECTOR`         | Rotate to align parallel with specified vector
 `turn ROTATION`        | Apply specified rotation
 `blur REAL`            | Add centered Gaussian noise of variance REAL
 `to X Y Z`             | Interpolate with the previously specified position
 `or POSITION`          | flip randomly between two specified positions
 
 A vector is set according to SHAPE, and the transformations are applied one after 
 the other, in the order in which they were given.\n

 Examples:
 @code
   position = 1 0 0
   position = circle 3 at 1 0
   position = square 3 align 1 1 0 at 1 1
 @endcode
 */ 
Vector Movable::readPosition(std::istream& is, const Space* spc)
{
    std::string tok;
    Vector pos(0,0,0);
    std::streampos isp = 0;
    
    try
    {
        if ( is.fail() )
            return pos;
        
        isp = is.tellg();
        pos = readPrimitive(is, spc);
        is.clear();
        
        while ( !is.eof() )
        {
            isp = is.tellg();
            tok = Tokenizer::get_symbol(is);

            if ( !is.good() )
                break;

            if ( tok.size() == 0 || !isalpha(tok[0]) )
                throw InvalidParameter("keyword expected here");
            
            // Translation is specified with 'at' or 'move'
            if ( tok == "at"  ||  tok == "move" )
            {
                Vector vec(0,0,0);
                is >> vec;
                pos = pos + vec;
            }
            // Convolve with shape
            else if ( tok == "add" )
            {
                Vector vec = readPrimitive(is, spc);
                pos = pos + vec;
            }
            // Alignment with a vector is specified with 'align'
            else if ( tok == "align" )
            {
                Vector vec = readDirection(is, pos, spc);
                Rotation rot = Rotation::rotationToVector(vec, RNG);
                pos = rot * pos;
            }
            // Rotation is specified with 'turn'
            else if ( tok == "turn" )
            {
                Rotation rot = readRotation(is, pos, spc);
                pos = rot * pos;
            }
            // Gaussian noise specified with 'blur'
            else if ( tok == "blur" )
            {
                real blur = 0;
                is >> blur;
                pos += Vector::randG(blur);
            }
            // This is handled in another function, and we just rewind
            else if ( tok == "to" )
            {
                is.seekg(isp);
                return pos;
            }
            // Alternative specified with 'or'
            else if ( tok == "or" )
            {
                if ( RNG.flip() )
                    return readPosition(is, spc);
            }
            else if ( tok == "and" )
            {
                std::string c = Tokenizer::get_symbol(is);
                std::string o = Tokenizer::get_symbol(is);
                std::string p = Tokenizer::get_symbol(is);
                ///\todo: specify a croping plane, eg:  "circle 4 and X > 0"
                throw InvalidParameter("unimplemented feature `and'");
            }
            else
            {
                /*
                We need to work around a bug in the stream extraction operator,
                which eats extra characters ('a','n','e','E') if doubles are read
                19.10.2015
                */
                is.clear();
                is.seekg(isp);
                is.seekg(-1, std::ios_base::cur);
                char c = is.peek();
                if ( c=='a' || c=='b' )
                    continue;
                
                throw InvalidParameter("unknown transformation `"+tok+"'");
            }
        }
    }
    catch ( InvalidParameter& e )
    {
        e << "\n" << StreamFunc::get_line(is, isp, PREF);
        throw;
    }
    return pos;
}


//------------------------------------------------------------------------------
/**
 Reads a direction which is a unit vector (norm = 1):
 
 Keyword                                       |  Resulting Vector
 ----------------------------------------------|------------------------------------------------------------
 `REAL REAL REAL`                              | the vector of norm 1 co-aligned with given vector
 `parallel REAL REAL REAL`                     | one of the two vectors of norm 1 parallel with given vector
 `orthogonal REAL REAL REAL`                   | a vector of norm 1 perpendicular to the given vector
 `horizontal` \n `parallel X`                  | (+1,0,0) or (-1,0,0), randomly chosen with equal chance
 `vertical`\n `parallel Y`                     | (0,+1,0) or (0,-1,0), randomly chosen with equal chance
 `parallel Z`                                  | (0,0,+1) or (0,0,-1), randomly chosen with equal chance
 `parallel XY`\n `parallel XZ`\n `parallel YZ` | A random vector in the specified plane
 `radial`                                      | directed from the origin to the current point
 `circular`                                    | perpendicular to axis joining the current point to the origin
 `or DIRECTION`                                | flip randomly between two specified directions

 
 If a Space is defined, one may also use:
 
 Keyword         |   Resulting Vector
 ----------------|----------------------------------------------------
 `tangent`       | parallel to the surface of the Space
 `normal`        | perpendicular to the surface
 `centrifuge`    | normal to the surface, directed outward
 `centripete`    | normal to the surface, directed inward


 Note: when the rotation is not uniquely determined in 3D (eg. `horizontal`), 
 cytosim will pick uniformly among all the possible rotations that fulfill the requirements.
 */


Vector Movable::readDirection(std::istream& is, Vector const& pos, const Space* spc)
{
    char c = Tokenizer::get_space(is, false);
    
    if ( c == EOF )
        return Vector::randU();

    if ( isalpha(c) )
    {
        const std::string tok = Tokenizer::get_symbol(is);
        
        if ( tok == "random" )
            return Vector::randU();

        if ( tok == "parallel" )
        {
            char c = Tokenizer::get_space(is, false);
            
            // an axis or a plane can be specified:
            if ( c == 'X' || c == 'Y' || c == 'Z' )
            {
                std::string k = Tokenizer::get_symbol(is);
                
                if ( k == "X" )
                    return Vector(RNG.sflip(), 0, 0);
                if ( k == "Y" )
                    return Vector(0, RNG.sflip(), 0);
                if ( k == "Z" && DIM == 3 )
                    return Vector(0, 0, RNG.sflip());
                if ( k == "XY" )
                {
#if ( DIM < 3 )
                    return Vector::randU();
#else
                    Vector2 h = Vector2::randU();
                    return Vector(h.XX, h.YY, 0);
#endif
                }
#if ( DIM == 3 )
                if ( k == "XZ" && DIM == 3 )
                {
                    Vector2 h = Vector2::randU();
                    return Vector(h.XX, 0, h.YY);
                }
                if ( k == "YZ" && DIM == 3 )
                {
                    Vector2 h = Vector2::randU();
                    return Vector(0, h.XX, h.YY);
                }
#endif
                throw InvalidParameter("Unexpected keyword `"+k+"' after `parallel`");
            }
            
            Vector vec;
            if ( is >> vec )
                return vec.normalized();
            throw InvalidParameter("expected vector after `parallel`");
        }
        
        if ( tok == "orthogonal" )
        {
            Vector vec;
            if ( is >> vec )
                return vec.randOrthoU(1);
            throw InvalidParameter("expected vector after `orthogonal`");
        }

        if ( tok == "horizontal" )
            return Vector(RNG.sflip(), 0, 0);
        
        if ( tok == "vertical" )
            return Vector(0, RNG.sflip(), 0);
        
        if ( tok == "radial" )
            return pos.normalized();

        if ( tok == "circular" )
            return pos.randOrthoU(1);
        
#if ( DIM >= 2 )
        if ( tok == "orthoradial" )
        {
            // in the XY plane orthogonal to the position
            real x = pos.XX;
            real y = pos.YY;
            real s = RNG.sflip() / sqrt( x*x + y*y );
            return Vector(-s*y, s*x, 0);
        }
#endif
      
        if ( spc )
        {
            if ( tok == "tangent" )
                return spc->normalToEdge(pos).randOrthoU(1);
           
#if ( DIM == 3 )
            Vector dir(0,0,1);
#elif ( DIM == 2 )
            real dir = 1;
#endif
            
#if ( DIM > 1 )
            if ( tok == "clockwise" )
                return cross(dir, spc->normalToEdge(pos));
            
            if ( tok == "anticlockwise" )
                return -cross(dir, spc->normalToEdge(pos));
#endif
            
            if ( tok == "normal" )
                return RNG.sflip() * spc->normalToEdge(pos);
            
            if ( tok == "centrifuge" )
                return -spc->normalToEdge(pos);
            
            if ( tok == "centripete" )
                return spc->normalToEdge(pos);
        }
        
        throw InvalidParameter("Unknown direction `"+tok+"'");
    }
    
    // accept a Vector:
    Vector vec(0,0,0);
    if ( is >> vec )
    {
        real n = vec.norm();
        if ( n < REAL_EPSILON )
            throw InvalidParameter("vector is singular (norm is too small)");
        return vec / n;
    }
    throw InvalidParameter("expected vector specifying a `direction`");
}


/**
 The initial orientation of objects is defined by a rotation, which can be
 specified as follows:
 
 Keyword                   |   Rotation / Result
 --------------------------|--------------------------------------------------------------------
 `random`                  | A rotation selected uniformly among all possible rotations
 `identity`                | The object is not rotated
 `angle A B C`             | As specified by 3 Euler angles in radians (in 2D, only A is needed)
 `degree A B C`            | As specified by 3 Euler angles in degrees (in 2D, only A is needed)
 `quat q0 q1 q2 q3`        | As specified by the Quaternion (q0, q1, q2, q3)
 DIRECTION                 | see @ref Movable::readDirection
 DIRECTION or DIRECTION    | flip randomly between two specified directions
 
 In the last case, a rotation will be built that transforms (1, 0, 0) into the given vector,
 after normalization. In 3D, this does not define the rotation uniquely (eg. `horizontal`),
 and cytosim will randomly pick one of the possible rotations that fulfill the requirements,
 with equal probability for all.
*/

Rotation Movable::readRotation(std::istream& is, Vector const& pos, const Space* spc)
{
    char c = Tokenizer::get_space(is, false);
    
    if ( c == EOF )
        return Rotation::randomRotation(RNG);

    std::streampos isp = is.tellg();
    if ( isalpha(c) )
    {
        std::string tok = Tokenizer::get_symbol(is);
        
        if ( tok == "random" )
            return Rotation::randomRotation(RNG);
        else if ( tok == "identity" || tok == "none" )
        {
            return Rotation::one();
        }
        else if ( tok == "angle" )
        {
            Torque a;
            is >> a;
#if ( 1 )
            std::streampos ist = is.tellg();
            std::string unit;
            is >> unit;
            if ( unit == "degree" )
                a *= M_PI/180.0;
            else
                is.seekg(ist);
#endif
            return Rotation::rotationFromAngles(a);
        }
        else if ( tok == "degree" )
        {
            Torque a;
            is >> a;
            a *= M_PI/180.0;
            return Rotation::rotationFromAngles(a);
        }
#if ( DIM == 3 )
        else if ( tok == "quat" )
        {
            Quaternion<real> quat;
            is >> quat;
            quat.normalize();
            Rotation rot;
            quat.setMatrix3(rot.data());
            return rot;
        }
#endif
        
        is.clear();
        is.seekg(isp);
    }

    // The last option is to specity a vector:
    Vector vec = readDirection(is, pos, spc);
    
    isp = is.tellg();
    if ( "or" == Tokenizer::get_symbol(is) )
    {
        if ( RNG.flip() )
            vec = readDirection(is, pos, spc);
    }
    else
    {
        is.seekg(isp);
    }
    
    /*
     A single Vector does not uniquely define a rotation in 3D:
     return a random rotation that is picked uniformly among
     the all possible rotations transforming (1,0,0) in vec.
     */
    return Rotation::rotationToVector(vec, RNG);
}


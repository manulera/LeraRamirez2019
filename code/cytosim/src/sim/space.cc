// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space.h"
#include "space_prop.h"
#include "exceptions.h"
#include "point_interpolated.h"
#include "messages.h"
#include "iowrapper.h"
#include "meca.h"
extern Random RNG;


Space::Space(const SpaceProp* p) 
: prop(p)
{
    assert_true(prop);
    
    for ( unsigned int d = 0; d < DMAX; ++d )
        setLength(d, 0);
    
    readLengths(p->dimensions);
}


Space::~Space()
{
    //std::clog << "~Space(" << prop->name() << ")\n";
    prop = 0;
}


void Space::readLengths(const std::string& str)
{
    //std::clog << "Space::readLengths(" << str << ")\n";
    unsigned d = 0;
    char const* ptr = str.c_str();
    char const* end = ptr + str.size();
    char * nxt = const_cast<char*>(ptr);

    while ( ptr < end )
    {
        real s = strtod(ptr, &nxt);
        if ( nxt == ptr )
            break;
        
        setLength(d++, s);
        ptr = nxt;
    }
    if ( d )
        resize();

    if ( nxt < end )
        std::cerr << "Warning: ignored trailing `" << nxt << "' in Space:geometry\n";
}


/**
 Checks that the number of specified dimensions is greater or equal to `required`,
 and if `positive == true`, also check that they are non-negative values.
 */
void Space::checkLengths(unsigned required, int strict) const
{
    for ( unsigned d = 0; d < required; ++d )
    {
        if ( strict && length(d) <= 0 )
            throw InvalidParameter(prop->name()+":dimension[",d,"] must be > 0");
        if ( length(d) < 0 )
            throw InvalidParameter(prop->name()+":dimension[",d,"] must be >= 0");
    }
}


void Space::setLength(unsigned int d, const real v)
{
    //std::clog << "Space `" << prop->name() << "' length[" << d << "] = " << v << "\n";
    if ( d < DMAX )
    {
        mLength[d]    = v;
        mLength2[d]   = 2*v;
        mLengthSqr[d] = v*v;
    }
    else
        throw InvalidParameter("exceeded space:dimension hard-coded limits");
}


void Space::resize(unsigned int d, const real v)
{
    setLength(d, v);
    //std::cerr << " dim[" << d << "] = " << mLength[d] << std::endl;
    resize();
}


//------------------------------------------------------------------------------
#pragma mark - Random Places

/**
 Provide a uniform random distribution in the volume by Monte-Carlo.
 
 Algorithm: throw a point in the rectangular volume provided by extension()
 until inside() returns true.
*/
Vector Space::randomPlace() const
{
    unsigned long nb_trials = 1<<13;
    Vector res, inf, sup;
    boundaries(inf, sup);
    Vector dif = sup - inf;
    
    unsigned long ouf = 0;
    do {
        
        res = inf + dif.e_mul(Vector::prand());

        if ( ++ouf > nb_trials )
        {
            MSG.warning("placement failed in Space::randomPlace()\n");
            return Vector(0,0,0);
            //throw InvalidParameter("placement failed in Space::randomPlace()");
        }
        
    } while ( ! inside(res) );
    
    return res;
}


/**
 Return a `point` for which:
 - inside(point) = true
 - inside(point, radius) = false
 */
Vector Space::randomPlaceNearEdge(real rad, unsigned long nb_trials) const
{
    unsigned long ouf = 0;
    Vector res;
    do {
        res = randomPlace();
        assert_true( inside(res) );
        if ( ++ouf > nb_trials )
            throw InvalidParameter("placement failed in Space::randomPlaceNearEdge()");
    } while ( allInside(res, rad) );
    return res;
}


/**
 Return a random point on the edge, using this method:
 - toss a random point `pos` within the range extended by `rad`.
 - project `pos` on the edge
 - return projection if the distance to `pos` is less than `rad`
 .
 */
Vector Space::randomPlaceOnEdge(real rad, unsigned long nb_trials) const
{
    if ( rad <= 0 )
        throw InvalidParameter("Space:distance_to_edge must be > 0");

    unsigned long ouf = 0;
    real d, rr = rad * rad;
    Vector pos, res, inf, dif;
    
    boundaries(inf, dif);
    inf -= Vector(rad, rad, rad);
    dif += Vector(rad, rad, rad) - inf;
    
    do {
        pos = inf + dif.e_mul(Vector::prand());
        project(pos, res);
        d = ( pos - res ).normSqr();
        if ( ++ouf > nb_trials )
            throw InvalidParameter("placement failed in Space::randomPlaceOnEdge()");
    } while ( d > rr );
    
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Inside/Outside

/**
 A bead is entirely inside if:
 - its center is inside,
 - the minimal distance (center-to-edge) is greater than the radius
 .
 */
bool Space::allInside(const real center[], const real rad) const
{
    if ( ! inside(center) )
        return false;

    return ( distanceToEdgeSqr(center) >= rad * rad );
}

/**
 A bead is entirely outside if:
 - its center is outside,
 - the minimal distance (center-to-edge) is greater than the radius
 .
 
 Attention: this is not equivalent to !allInside(center, radius)
 */
bool Space::allOutside(const real center[], const real rad) const
{
    if ( inside(center) )
        return false;
        
    return ( distanceToEdgeSqr(center) >= rad * rad );
}

//------------------------------------------------------------------------------
#pragma mark - Project

/**
this code is equivalent to SpaceInflate::project(), with a negative radius
 */
void Space::project(const real pos[], real prj[], const real rad) const
{
    if ( rad < 0 )
        ABORT_NOW("radius should not be negative");

    project(pos, prj);

    
    ///\todo problem in project() with radius if point is exactly on the box (n==0)
    //if (n==0) we do not know the orthogonal direction to follow. We should
    //take another point near by, and project from there.
    

    real n = 0, x, pw[DIM];
    for ( int d = 0; d < DIM; ++d )
    {
        x = pos[d] - prj[d];
        pw[d] = x;
        n += x * x;
    }
            
    if ( n > 0 )
        n = ( inside(pos) ? +rad : -rad ) / sqrt(n);
    else {
        throw Exception("in project(..., radius): the point is on the edge");
        //printf("point % .3f % .3f % .3f :", pos[0], pos[1], pos[2]);
        //printf("inside = %i :", inside(point));
        //printf("proj  % .3f % .3f % .3f\n", prj[0], prj[1], prj[2]);
    }
    
    for ( int d=0; d<DIM; ++d )
        prj[d] += n * pw[d];
}

//------------------------------------------------------------------------------
#pragma mark - Misc


real Space::max_extension() const
{
    Vector inf, sup;
    boundaries(inf, sup);
    return std::max(inf.norm_inf(), sup.norm_inf());
}

/**
 The volume is estimated with a simple monte-carlo approach:
 - throw points in the rectangular volume provided by boundaries()
 - count how many are inside the volume with inside()
 .
 Then
 @code
 volume ~ ( number-of-points-inside / number-of-point ) * volume-of-rectangle
 @endcode
 */
real Space::estimateVolume(unsigned long cnt, bool verbose) const
{
    Vector inf, sup;
    boundaries(inf, sup);
    Vector dif = sup - inf;
    
    unsigned long in = 0;
    for ( unsigned long i = 0; i < cnt; ++i )
    {
        Vector pos;
        pos.XX = inf.XX + dif.XX * RNG.preal();
#if ( DIM > 1 )
        pos.YY = inf.YY + dif.YY * RNG.preal();
#endif
#if ( DIM > 2 )
        pos.ZZ = inf.ZZ + dif.ZZ * RNG.preal();
#endif
        in += inside(pos);
    }
    
    real vol = in / real(cnt);
    for ( int d = 0; d < DIM; ++d )
        vol *= dif[d];
    
    if ( verbose )
        printf("Monte-Carlo estimated volume of `%s` is %.6f +/- %.6f\n",
               prop->geometry.c_str(), vol, vol/sqrt(in));

    return vol;
}


/**
 This uses Space::project to reflect `w` on the edge of the Space,
 until the result eventually falls inside.
 
 In most geometries, this works well, but if the distance from the point
 to the edge is very large compared to the width of the space, the number
 of iterations may be large.
*/
void Space::bounce(Vector& pos) const
{
    Vector p;
    
    // bounce once on the edge, and return if inside
    project(pos, p);
    pos = p + p - pos;
    
    if ( !inside(pos) )
    {
        // bounce on the edge, and return if inside
        int cnt = 0;
        do {
            project(pos, p);
            pos = p + p - pos;
            if ( inside(pos) )
                return;
        } while ( ++cnt < 8 );

        std::cerr << "Warning: Space:bounce once failed: reduce time_step?\n";

        // Place point on edge, if iterations have failed to bring it inside:
        project(pos, p);
        pos = p;
    }
}


/** 
 `normalToEdge(const Vector&)` uses an iterative method to find
 the normal to the edge, using Space::project().
 
 If you know for certain that `point[]` is far from the edge,
 the normal can be more directly obtained from the projection:
 @code 
 project(point, proj);
 normal = ( proj - point ).normalized()
 @endcode
 
 */
Vector Space::normalToEdge(const Vector& pos) const
{
    const real goal = 10000*REAL_EPSILON*REAL_EPSILON;
    
    Vector P, M, prj, res;
    project(pos, prj);

    real H = 1;
    for ( unsigned i = 0; i < 12; ++i )
    {
        H /= 2;
        for ( unsigned j = 0; j < 16; ++j )
        {
            //start from a random vector:
            res = Vector::randU(H);
            
            for ( unsigned n = 0; n < 32; ++n )
            {
                project(prj+res, P);
                project(prj-res, M);
                
                // refine the estimate:
                Vector ref = 0.5 * ( M - P );
                res += ref;
                
                // check convergence:
                if ( ref.normSqr() < goal )
                {
                    if ( 2 * res.norm() < H )
                        res.normalize(H);
                    else
                    {
                        if ( inside(prj+res) )
                            return res.normalized(-1);
                        else
                            return res.normalized();
                    }
                }
            }
        }
    }
    
    printf("warning: convergence failure in normalToEdge()\n");
    printf("         error = %e at height = %e\n", P.distance(prj), H);
    if ( inside(prj+res) )
        return res.normalized(-1);
    else
        return res.normalized();
}


//------------------------------------------------------------------------------
real  Space::distanceToEdgeSqr(Vector const& pos) const
{
    real prj[DIM];
    
    project(pos, prj);
    
    real res = (pos[0]-prj[0]) * (pos[0]-prj[0]);
    for ( int dd=1; dd < DIM; ++dd )
        res += (pos[dd]-prj[dd]) * (pos[dd]-prj[dd]);
    
    return res;
}


real  Space::signedDistanceToEdge(Vector const& pos) const
{
    if ( inside(pos) )
        return -distanceToEdge(pos);
    else
        return +distanceToEdge(pos);
}



//------------------------------------------------------------------------------
#pragma mark - Interactions

/**
 Call the appropriate interaction from `meca`, to force `pe` to be on the edge of the Space.
 
 This implementation uses `pos` to find the local normal to the edge of the Space.
 and then calls Meca::addPlaneClamp, with the approprimate aguments.
 This generates a friction-less potential centered on the edge.
 */

void Space::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    Vector prj;
    project(pos, prj);
    Vector dir = pos - prj;
    real n = dir.normSqr();
    if ( n > 0 )
        meca.addPlaneClamp( pe, prj, dir, stiff/n );
}

/**
 Call the appropriate interaction from `meca`, to confine `pe`, which is at position `pos`.
 
 The default implementation projects `pos`,
 to calculate the direction of the normal to the edge of the Space,
 and then calls Meca::addPlaneClamp, with the approprimate aguments.
 This generates a friction-less potential centered on the edge.
 */

void Space::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    Vector prj;
    project(pos, prj, rad);
    Vector dir = pos - prj;
    real n = dir.normSqr();
    if ( n > 0 )
        meca.addPlaneClamp( pe, prj, dir, stiff/n );
}

#if ( 0 )

/**
 This calls Space::setInteraction(pos, PointExact, meca, stiff) twice,
 to generate a force on `pi` (which is at position `pos`) toward the surface.
 */
void Space::setInteraction(Vector const& pos, PointInterpolated const& pi, Meca & meca, real stiff) const
{
    setInteraction(pos, pi.exact1(), meca, pi.coef2()*stiff);
    setInteraction(pos, pi.exact2(), meca, pi.coef1()*stiff);
}


/**
 This will add a force component if:
 - ( conf == inside ) && ! Space::inside(pos)
 - ( conf == outside ) &&  Space::inside(pos)
 - ( conf == surface )
 .
 */
void Space::setInteraction(PointInterpolated const& pi, Meca & meca, real stiff, Confinement conf) const
{
    if ( conf == CONFINE_ON )
    {
        setInteraction(pi.pos(), pi, meca, stiff);
    }
    else if ( conf == CONFINE_INSIDE )
    {
        Vector pos = pi.pos();
        if ( ! inside(pos) )
            setInteraction(pos, pi, meca, stiff);
    }
    else if ( conf == CONFINE_OUTSIDE )
    {
        Vector pos = pi.pos();
        if ( inside(pos) )
            setInteraction(pos, pi, meca, stiff);
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark - IO


void Space::write(Outputter& out) const
{
    out.put_line(prop->shape+" ", ' ');
    out.writeUInt16(DMAX);
    for ( unsigned d = 0; d < DMAX; ++d )
        out.writeFloat(mLength[d]);
}


void Space::read(Inputter& in, Simul&, Tag)
{
#ifdef BACKWARD_COMPATIBILITY

    if ( in.formatID() < 35 )
        return;
    
    if ( in.formatID() < 36 )
    {
        for ( unsigned d = 0; d < 3; ++d )
            setLength(d, in.readFloat());
        resize();
        return;
    }
    
    if ( in.formatID() < 41 )
    {
        unsigned n = in.readUInt8();
        for ( unsigned d = 0; d < n; ++d )
            setLength(d, in.readFloat());
        resize();
        return;
    }
    
#endif
    
    // read the 'shape' stored as a space-terminated string
    int c = in.get_char();
    if ( c != ' ' )
        in.unget(c);
    
    std::string str;
    in.get_line(str, ' ');
    
    // check that this matches current Space:
    if ( str.compare(0, prop->shape.size(), prop->shape) )
    {
        std::cerr << PREF << "Space `" << prop->name() << "':";
        std::cerr << " file: " << str << " property: " << prop->shape << '\n';
        throw InvalidIO("shape missmatch");
    }
    
    unsigned n;
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 43 )
        n = in.readUInt8();
    else
#endif
        n = in.readUInt16();
    
    for ( unsigned d = 0; d < n; ++d )
        setLength(d, in.readFloat());
    
    resize();
}

//------------------------------------------------------------------------------
#pragma mark - Display


#ifdef DISPLAY
#include "gle.h"


void Space::displaySection(const int dim, const real pos, const real step) const
{
    Vector inf, sup;
    boundaries(inf, sup);
    Vector q, p( pos, pos, pos );
    int xx = ( dim + 1 ) % DIM;
    int yy = ( xx + 1 ) % DIM;
    
    real xs = sup[xx];
    real ys = sup[yy];
    real inc = step * ( xs > ys ? xs : ys );

    glBegin(GL_LINE_LOOP);
    p[yy] = ys;
    for ( real a = -xs; a < xs; a += inc )
    {
        p[xx] = a;
        project(p, q);
        gle::gleVertex(q);
    };
    p[xx] = xs;
    for ( real a = -ys; a < ys; a += inc )
    {
        p[yy] = -a;
        project(p, q);
        gle::gleVertex(q);
    };
    p[yy] = -ys;
    for ( real a = -xs; a < xs; a += inc )
    {
        p[xx] = -a;
        project(p, q);
        gle::gleVertex(q);
    };
    p[xx] = -xs;
    for ( real a = -ys; a < ys; a += inc )
    {
        p[yy] = a;
        project(p, q);
        gle::gleVertex(q);
    };
    glEnd();
}

#else

void Space::displaySection(const int dim, const real pos, const real step) const
{
    //you will get this output if objects for play was not compiled properly:
    //DISPLAY should be defined on the compiler command, with: -DDISPLAY
    printf("dummy Space::displaySection()");
}

#endif

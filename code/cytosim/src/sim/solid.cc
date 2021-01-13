// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "assert_macro.h"
#include "solid.h"
#include "solid_prop.h"
#include "exceptions.h"
#include "hand_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"
#include "simul.h"
#include "space.h"
#include "wrist.h"


#if ( DIM == 3 )
#   include "quaternion.h"
#   include "matrix3.h"
#   include "clapack.h"
#endif

extern Random RNG;



//------------------------------------------------------------------------------
#pragma mark -


void Solid::step()
{
}


void Solid::setInteractions(Meca & meca) const
{
#ifdef NEW_RADIAL_FLOW
    PRINT_ONCE("NEW_RADIAL_FLOW enabled: Solids converge to the same point\n");
    /// Special code for Maria Burdyniuk
    real now = objset()->simul.time();
    if ( prop->flow_time[0] > now )
    {
        PointExact pt(this,0);
        Vector dir = prop->flow_center - pt.pos();
        real s = dragCoefficient() / ( prop->flow_time[1] - now );
        meca.addForce(pt, dir * s);
    }
#endif
    
    switch ( prop->confine )
    {
        case CONFINE_OFF:
            break;
            
        case CONFINE_INSIDE:
        {
            Space const* spc = prop->confine_space_ptr;
            
            for ( unsigned pp = 0; pp < nbPoints(); ++pp )
            {
                const real rad = soRadius[pp];
                // confine all massive points:
                if ( rad > 0 )
                {
                    Vector pos = posP(pp);
                    if ( ! spc->inside(pos) )
                        spc->setInteraction(pos, PointExact(this, pp), meca, prop->confine_stiffness);
                }
            }
        } break;
            
        case CONFINE_OUTSIDE:
        {
            Space const* spc = prop->confine_space_ptr;
            
            for ( unsigned pp = 0; pp < nbPoints(); ++pp )
            {
                const real rad = soRadius[pp];
                // confine all massive points:
                if ( rad > 0 )
                {
                    Vector pos = posP(pp);
                    if ( spc->inside(pos) )
                        spc->setInteraction(pos, PointExact(this, pp), meca, prop->confine_stiffness);
                }
            }
        } break;
            
        case CONFINE_ALL_INSIDE:
        {
            Space const* spc = prop->confine_space_ptr;
            
            for ( unsigned pp = 0; pp < nbPoints(); ++pp )
            {
                const real rad = soRadius[pp];
                // confine all massive points:
                if ( rad > 0 )
                {
                    Vector pos = posP(pp);
                    if ( ! spc->allInside(pos, rad) )
                        spc->setInteraction(pos, PointExact(this, pp), rad, meca, prop->confine_stiffness);
                }
            }
        } break;
            
        case CONFINE_ON:
        {
            Space const* spc = prop->confine_space_ptr;
            
            for ( unsigned pp = 0; pp < nbPoints(); ++pp )
            {
                // only confine massive points:
                if ( soRadius[pp] > 0 )
                    spc->setInteraction(posP(pp), PointExact(this, pp), meca, prop->confine_stiffness);
            }
        } break;
        
        case CONFINE_POINT0:
        {
            Space const* spc = prop->confine_space_ptr;
            // confine Point 0:
            spc->setInteraction(posP(0), PointExact(this, 0), meca, prop->confine_stiffness);
        } break;

        default:
            throw InvalidParameter("Invalid solid::confine");
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/// set 'mat' of order DIM to `val * I`
void set_diagonal_matrix(real * mat, const real val)
{
#if ( DIM == 1 )
    mat[0] = val;
#elif ( DIM == 2 )
    mat[0] = val;
    mat[1] = 0.0;
    mat[2] = 0.0;
    mat[3] = val;
#else
    mat[0] = val;
    mat[1] = 0.0;
    mat[2] = 0.0;
    mat[3] = 0.0;
    mat[4] = val;
    mat[5] = 0.0;
    mat[6] = 0.0;
    mat[7] = 0.0;
    mat[8] = val;
#endif
}



void Solid::reset()
{
    soRadius    = 0;
    soShape     = 0;
    soShapeSize = 0;
    soDrag      = 0;
    set_diagonal_matrix(soMom, 1.0);
    soReshapeTimer = RNG.pint(7);
}


Solid::Solid (SolidProp const* p)
: prop(p)
{
    reset();
}


Solid::Solid(const Solid & o)
: PointSet(o)
{
    reset();
    prop = o.prop;
    allocatePoints(o.nbPoints());
    for ( unsigned p = 0; p < nbPoints(); ++p )
        soRadius[p] = o.soRadius[p];
    fixShape();
}


Solid & Solid::operator =(const Solid & o)
{
    reset();
    prop = o.prop;
    PointSet::operator=(o);
    for ( unsigned p = 0; p < nbPoints(); ++p )
        soRadius[p] = o.soRadius[p];
    fixShape();
    return *this;
}


Solid::~Solid()
{
    deallocatePoints();
    prop = 0;
}


/**
 This calls PointSet::allocatePoints().
 If PointSet::allocatePoints() allocated memory, it will return the 
 size of the new array, and in that case, the same size is allocated for other arrays.
 */
unsigned Solid::allocatePoints(const unsigned nbp)
{
    unsigned ms = PointSet::allocatePoints(nbp);
    if ( ms )
    {
        //std::clog << "Solid::allocatePoints " << ms << std::endl;
        
        // allocate a new array of the right size:
        real  *  soShape_new = new real[DIM*ms];
        real  *  soRadius_new = new real[ms];
        
        //set the radii to zero (no drag) by default:
        for ( unsigned p = 0; p < ms; ++p )
            soRadius_new[p] = 0;
        
        // copy the current values in the new array:
        if ( soShape )
        {
            for ( unsigned p = 0; p < nbPoints(); ++p )
            {
                soRadius_new[p] = soRadius[p];
                for ( int d = 0; d < DIM; ++d )
                    soShape_new[DIM*p+d] = soShape[DIM*p+d];
            }
            // delete the 'current' array:
            delete[] soShape;
            delete[] soRadius;
        }
        // the 'new' array becomes the 'current' one:
        soShape = soShape_new;
        soRadius = soRadius_new;
    }
    return ms;
}


void Solid::deallocatePoints()
{
    PointSet::deallocatePoints();
    if ( soRadius )
    {
        delete[] soRadius;
        soRadius = 0;
    }
    if ( soShape )
    {
        delete[] soShape;
        soShape = 0;
    }
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 @ingroup NewObject
 
 There are different ways to specify the number and positions of points in a Solid:
 
 @code
 new solid NAME
 {
   point0 = [INTEGER,] POSITION, RADIUS [, SINGLE_SPEC]
   point1 = [INTEGER,] POSITION, RADIUS [, SINGLE_SPEC]
   point2 = [INTEGER,] POSITION, RADIUS [, SINGLE_SPEC]
   etc.
 }
 @endcode
 
 each `point#` specifies a number of points to be added.
 The first parameter (`INTEGER`) specifies the number of points.
 The second argument (`POSITION`) specifies their position with respect to the center.
 The keywords are the same as for other position in cytosim (see examples below).
 The last argument (`RADIUS`) specifies the radius of the bead attached at this point,
 and it can be zero.
 
 Examples:
 
 @code
 new solid blob
 {
   point0 = center, 1.0
   point1 = 10, sphere 1, 0, grafted
   ...
 }
 @endcode
 
 `POSITION` can be a `VECTOR`, or the usual keywords:
 - `center`
 - `ball RADIUS`
 - `sphere RADIUS`
 - `equator RADIUS`
 .
 
 Another way to specify points of a Solid:

 @code
 new solid NAME
 {
   sphere0 = POSITION, RADIUS [, SINGLE_SPEC]
   sphere1 = POSITION, RADIUS [, SINGLE_SPEC]
   etc.
 }
 @endcode
 
 each `sphere#` specifies one sphere to be added.
 The first argument (`POSITION`) specifies the position with respect to the center.
 The keywords are the same as for other position in cytosim (see examples below).
 The second argument (`RADIUS`) specifies the radius of the bead attached at this point,
 and it should not be zero.
 

 <h3> Add Singles to a Solid </h3>
 
 The parameter 'attach' can be used to add Single to the points of a Solid:
 
 @code
 new solid NAME
 {
   point0   = ... , SINGLE_SPEC
   sphere0  = ... , SINGLE_SPEC
   etc.
   attach   = SINGLE_SPEC [, SINGLE_SPEC] ...
   attach0  = SINGLE_SPEC [, SINGLE_SPEC] ...
   attach1  = SINGLE_SPEC [, SINGLE_SPEC] ...
   etc.
 }
 @endcode
 
 Where `SINGLE_SPEC` is string containing at most 3 words: `[INTEGER] NAME [each]`,
 where the `INTEGER` specifies the number of Singles, `NAME` specifies their name,
 and the optional word `each` species that the command applies to every point.
 
 The command `attach` applies to all the points of the Solid, while `attach0`,
 `attach1`, etc. apply to the points specified by `point0`, `point1`, etc. only.
 With `attach`, the Singles are distributed randomly on all the points,
 and if `each` is specified, the specification is repeated for each point.
 
 For example if `grafted` is the name of a Single, one can use:

 @code
 new solid NAME
 {
    attach0 = 1 grafted each
    attach1 = 10 grafted
 }
 @endcode
 */

ObjectList Solid::build(Glossary & opt, Simul& simul)
{
    ObjectList res;
    std::string str;
    unsigned inp = 0, inx = 0, nbp = 1;

    // interpret each instruction as a command to add points:
    std::string var = "point0";
    while ( opt.has_key(var) )
    {
        inx = 0;
        nbp = 1;
        // optionally specify a number of points
        if ( opt.is_number(var) == 2 && opt.set(nbp, var) )
            ++inx;
        
        if ( nbp > 0 )
        {
            // get sphere radius:
            real sr = 0;
            opt.set(sr, var, inx+1);
            
            if ( sr < 0 )
                throw InvalidParameter("the radius of solid:sphere must be >= 0");

            unsigned fip = nbPoints();
            // add 'nbp' points:
            for ( unsigned n = 0; n < nbp; ++n )
            {
                // get position:
                Vector vec(0,0,0);
                std::istringstream iss(opt.value(var, inx));
                vec = Movable::readPosition(iss, 0);
                addSphere(vec, sr);
            }
            
            // attach Single to this set of points:
            inx += 2;
            while ( opt.set(str, var, inx++) )
                res.append(simul.singles.makeWrists(this, fip, nbp, str));
            
            // attach Single to this set of points:
            inx = 0;
            var = "attach" + sMath::repr(inp);
            while ( opt.set(str, var, inx++) )
                res.append(simul.singles.makeWrists(this, fip, nbp, str));
        }
        
        var = "point" + sMath::repr(++inp);
    }
    
    // interpret each instruction as a command to add spheres:
    var = "sphere0";
    while ( opt.has_key(var) )
    {
        inx = 0;

        // get sphere radius:
        real sr = 0;
        opt.set(sr, var, inx+1);
        
        if ( sr <= 0 )
            throw InvalidParameter("the radius of sphere specified in solid must be > 0");

        // get position:
        Vector vec(0,0,0);
        std::istringstream iss(opt.value(var, inx));
        vec = Movable::readPosition(iss, 0);
        
        // add a bead with a local coordinate system
        unsigned ref = addSphere(vec, sr);
        addTriad(sr);
        
        // attach Single on the surface of this sphere:
        inx += 2;
        while ( opt.set(str, var, inx++) )
        {
            unsigned num = 1;
            std::istringstream iss(str);
            iss >> num;
            
            if ( iss.fail() )
            {
                num = 1;
                iss.clear();
            }
           
            iss >> str;
            
            SingleProp * sip = simul.findProperty<SingleProp*>("single", str);
            if ( sip == 0 )
                throw InvalidParameter("could not find fiber:attach single `"+str+"'");

            // add Wrists anchored on the local coordinate system:
            for ( unsigned i = 0; i < num; ++i )
            {
                Vector vec = Vector::randU();
                res.push_back(new Wrist(sip, this, ref, vec));
            }
        }
        var = "sphere" + sMath::repr(++inp);
    }

    
    // attach Singles to be distributed over all the points:
    inx = 0;
    while ( opt.set(str, "attach", inx++) )
        res.append(simul.singles.makeWrists(this, 0, nbPoints(), str));

    // final verification of the number of points:
    nbp = 0;
    if ( opt.set(nbp, "nb_points")  &&  nbp != nbPoints() )
    {
        throw InvalidParameter("could not find the number of points specified in solid:nb_points");
    }
    
    //std::cerr << *this << std::endl;
    return res;
}


unsigned Solid::addSphere(Vector const& vec, real rad)
{
    if ( rad < 0 )
        throw InvalidParameter("solid:sphere's radius should be >= 0");

    unsigned inx = addPoint(vec);
    soRadius[inx] = rad;
    //std::clog << "addSphere(" << vec << ", " << rad << ") for " << reference() << " index " << inx << "\n";
    return inx;
}


unsigned Solid::addTriad(real arm)
{
    if ( nbPoints() < 1 )
        throw InvalidParameter("cannot add Triad to solid without point");
    
    if ( arm <= 0 )
        throw InvalidParameter("solid::arm should be > 0");
    
    unsigned inx = lastPoint();

    //std::clog << "Solid::addTriad(" << arm << ") at index " << inx << "\n";
    Vector vec = posPoint(inx);
    
    if ( DIM > 0 ) addPoint(vec+Vector(arm,0,0));
    if ( DIM > 1 ) addPoint(vec+Vector(0,arm,0));
    if ( DIM > 2 ) addPoint(vec+Vector(0,0,arm));
    
    return inx;
}


void Solid::radius(const unsigned indx, real rad)
{
    assert_true( indx < nbPoints() );
    if ( rad < 0 )
        throw InvalidParameter("solid:radius must be positive");
    soRadius[indx] = rad;
}


Vector Solid::centroid() const
{
    if ( nbPoints() == 0 )
        ABORT_NOW("cannot calculate centroid of a Solid without point");
    
    if ( nbPoints() == 1 )
        return posP(0);
    
    Vector res(0,0,0);
    real sum = 0;
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
    {
        if ( soRadius[pp] > 0 )
        {
            res += soRadius[pp] * posP(pp);
            sum += soRadius[pp];
        }
    }
    if ( sum < REAL_EPSILON )
        ABORT_NOW("cannot calculate centroid of a Solid without drag sphere");
    
    res /= sum;
    return res;
}



/**
 fixShape() copies the current shape in the array soShape[],
 and calculates the moment of inertia of the ensemble of points.
 The reference soShape[] is used by 'reshape()', and 'rescale()'.
 */
void Solid::fixShape()
{
    if ( nbPoints() == 0 )
        throw InvalidParameter("Solid has no points!");
    
    //std::clog << "Fixing Solid " << reference() << " with " << nbPoints() << " points\n";
    
    Vector avg, sec;
    calculateMomentum(avg, sec, true);
    
    // store momentum of the current shape:
    soShapeSqr = sec.e_sum();
    
    //we store the current points:
    soShapeSize = nbPoints();
    //set reference to current shape translated by -G (center) :
    for ( unsigned nn = 0; nn < DIM*soShapeSize; nn += DIM )
    {
        for ( int dd = 0; dd < DIM; ++dd )
            soShape[nn+dd] = psPos[nn+dd] - avg[dd];
    }
    
    setDragCoefficient();
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 The function rescale the reference shape soShape[], that was specified last time fixShape() was called.
 If axis==-1 (default), then all dimensions are scaled uniformly.
 The next call to reshape() will then apply the new reference to the current shape.
 */
void Solid::scaleShape(const real sx, const real sy, const real sz)
{
    //scale in only in the specified dimension
    for ( unsigned pp = 0; pp < DIM * soShapeSize; pp += DIM)
    {
        soShape[pp  ] *= sx;
#if ( DIM > 1 )
        soShape[pp+1] *= sy;
#endif
#if ( DIM > 2 )
        soShape[pp+2] *= sz;
#endif
    }
    
    
    //recalculate the momentum needed in rescale():
    soShapeSqr = 0;
    for ( unsigned pp = 0; pp < DIM * soShapeSize; ++pp )
        soShapeSqr += soShape[pp] * soShape[pp];
    
    setDragCoefficient();
}


/**
 Rescale the current cloud of points around its center of gravity,
 to recover the same 'size' as the reference soShape[]. 
 Size is measured as sum( ( x - g )^2 ).
 */
void Solid::rescale()
{
    Vector avg, sec;
    calculateMomentum(avg, sec, true);
    
    // calculate the momentum of the current shape:
    real sz = sec.e_sum();
    
    if ( sz > 0 )
    {
        // calculate the scaling factor to restore the size to 'soShapeSqr':
        real scale = sqrt( soShapeSqr / sz );
    
        // scale the shape around the center of gravity:
        for ( unsigned p = 0; p < DIM * nbPoints(); p += DIM )
        {
            for ( unsigned d = 0; d < DIM; ++d )
                psPos[p+d] = scale * ( psPos[p+d] - avg[d] ) + avg[d];
        }
    }
}



/**
 reshapeReally() finds the best isometric transformation = rotation + translation
 to bring the reference (soShape[]) onto the current shape (PointSet::psPos[]),
 and then replaces psPos[] by the transformed soShape[]. 
 This restores the shape of the cloud of point which is stored in soShape[],
 into the current position and orientation of the object.
 The best translation is the ones that conserves the center of gravity,
 The best rotation is obtained differently in 2D and 3D, and is unique.

 @todo: store the rotation and translation calculated by reshapeReally()
*/

#if ( DIM == 1 )

void Solid::reshape()
{    
    //we check that the number of points is the same as when fixShape() was called.
    if ( soShapeSize != nbPoints() )
        ABORT_NOW("mismatch with current number of points: forgot to call fixShape()?");
         
    real cc = 0, a = 0;
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
    {
        a  += psPos[pp] * soShape[pp];
        cc += psPos[pp];
    }
    
    cc /= real( nbPoints() );
    real s = a / fabs(a);
    
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
        psPos[pp] = s * soShape[pp] + cc;
}

#elif ( DIM == 2 )

void Solid::reshape()
{    
    // the number of points should be the same as when fixShape() was called.
    if ( soShapeSize != nbPoints() )
        ABORT_NOW("mismatch with current number of points: forgot to call fixShape()?");
    
    Vector avg = PointSet::position();
    
    /*
     The best rotation is obtained by simple math on the cross products
     and vector products of soShape[] and psPos[]: (see it on paper)
    */
    
    real a = 0, b = 0;
    
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
    {
        a += psPos[DIM*pp] * soShape[DIM*pp  ] + psPos[DIM*pp+1] * soShape[DIM*pp+1];
        b += soShape[DIM*pp] * psPos[DIM*pp+1] - soShape[DIM*pp+1] * psPos[DIM*pp  ];
    }
    
    real n = sqrt( a*a + b*b );
    
    // cosine and sinus of the rotation:
    real c = 1, s = 0;
    if ( n > REAL_EPSILON ) {
        c = a / n;
        s = b / n;
    }
    
    //printf(" n %8.3f, c %8.3f, s %8.3f norm = %8.3f\n", n, c, s, c*c + s*s);
    
    // apply transformation = rotation + translation:
    
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
    {
        psPos[DIM*pp  ] = c * soShape[DIM*pp] - s * soShape[DIM*pp+1] + avg.XX;
        psPos[DIM*pp+1] = s * soShape[DIM*pp] + c * soShape[DIM*pp+1] + avg.YY;
    }
}

#elif ( DIM == 3 )

void Solid::reshape()
{
    // the number of points should be the same as when fixShape() was called.
    if ( soShapeSize != nbPoints() )
        ABORT_NOW("mismatch with current number of points: fixShape() was not called?");
    
    /*
     We follow the procedure described by Berthold K.P. Horn in
     "Closed-form solution of absolute orientation using unit quaternions"
     Journal of the optical society of America A, Vol 4, Page 629, April 1987
    */
    
    Vector avg = PointSet::position();
    
    real S[3*3] = { 0 };
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
    {
        for ( unsigned dd = 0; dd < DIM; ++dd )
            for ( unsigned ee = 0; ee < DIM; ++ee )
                S[dd+3*ee] += soShape[DIM*pp+dd] * psPos[DIM*pp+ee];
    }
    
    //the scaling is arbitrary here, but it keeps the magnitude of the matrix:
    real scale = 1.0 / ( fabs(S[0])+fabs(S[4])+fabs(S[8]) );
    
    real N[4*4];
    
    N[0+4*0] = scale * ( S[0+3*0] + S[1+3*1] + S[2+3*2] );
    N[0+4*1] = scale * ( S[1+3*2] - S[2+3*1] );
    N[0+4*2] = scale * ( S[2+3*0] - S[0+3*2] );
    N[0+4*3] = scale * ( S[0+3*1] - S[1+3*0] );
    N[1+4*1] = scale * ( S[0+3*0] - S[1+3*1] - S[2+3*2] );
    N[1+4*2] = scale * ( S[0+3*1] + S[1+3*0] );
    N[1+4*3] = scale * ( S[2+3*0] + S[0+3*2] );
    N[2+4*2] = scale * ( S[1+3*1] - S[0+3*0] - S[2+3*2] );
    N[2+4*3] = scale * ( S[1+3*2] + S[2+3*1] );
    N[3+4*3] = scale * ( S[2+3*2] - S[1+3*1] - S[0+3*0] );
    
    /* 
     Use lapack to find the largest Eigenvalue, and associated Eigenvector,
     which is the quaternion corresponding to the best rotation
     */
    
    int nbvalues;
    real eValue[4];
    Quaternion<real> quat;
    real work[8*4];
    int iwork[5*4];
    int ifail[4];
    
    int info = 0;
    lapack_xsyevx('V','I','U', 4, N, 4, 0, 0, 4, 4, REAL_EPSILON,
                  &nbvalues, eValue, quat, 4, work, 8*4, iwork, ifail, &info);
    
    //MSG("optimal LWORK = %i\n", work[0] );
    //MSG("eigenvalue %6.2f,", eValue[0]);
    //quat.println();
    
    if ( info == 0 )
    {
        //get the rotation matrix corresponding to the quaternion:
        quat.setMatrix3(S);
    
        //apply the transformation = rotation + translation:
        for ( unsigned pp = 0; pp < nbPoints(); ++pp )
        {
            psPos[DIM*pp  ] = avg.XX + S[0]*soShape[DIM*pp]+ S[3]*soShape[DIM*pp+1] + S[6]*soShape[DIM*pp+2];
            psPos[DIM*pp+1] = avg.YY + S[1]*soShape[DIM*pp]+ S[4]*soShape[DIM*pp+1] + S[7]*soShape[DIM*pp+2];
            psPos[DIM*pp+2] = avg.ZZ + S[2]*soShape[DIM*pp]+ S[5]*soShape[DIM*pp+1] + S[8]*soShape[DIM*pp+2];
        }
    }
    else
    {
        //apply translation:
        for ( unsigned pp = 0; pp < nbPoints(); ++pp )
        {
            psPos[DIM*pp  ] = avg.XX + soShape[DIM*pp  ];
            psPos[DIM*pp+1] = avg.YY + soShape[DIM*pp+1];
            psPos[DIM*pp+2] = avg.ZZ + soShape[DIM*pp+2];
        }
        
        printf("Solid::reshapeReally(): lapack_xsyevx() failed with code %i\n", info);
    }
}
#endif


/**
 
 getPoints() calls rescale() often and reshapeReally() occasionally, because
 - reshapeReally() corrects for all kind of numerical drift but is CPU expensive
 - rescale() corrects for 2d order numerical drift, which are dominant.
 .
 
 The calls for different solids are shifted by using the identity() of each Solid.
 */
void Solid::getPoints(const real * x)
{
    PointSet::getPoints(x);
    
    // for one point, nothing should be done
    if ( nbPoints() < 2 )
        return;
    
    if ( ++soReshapeTimer > 7 )
    {
        reshape();
        soReshapeTimer = 0;
    }
    else
        rescale();
}


//------------------------------------------------------------------------------
#pragma mark -


/**
 returns 6 * M_PI * viscosity * sum(radius);
 */
real Solid::dragCoefficient() const
{
    real sumR = 0;
    
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
        sumR += soRadius[pp];
    
    return 6 * M_PI * prop->viscosity * sumR;
}


/**
 Stokes relations:
 Translation:
   muT = 6 * M_PI * viscosity * radius;
   d(position)/dt = force / muT
 Rotation:
   muR = 8 * M_PI * viscosity * radius^3
   d(angle)/dt = force-torque / muR
 */
void Solid::setDragCoefficient()
{
    real sumR = 0;            //the total drag coef.
    soDragRot = 0;            //the total rotational drag coef.  
    soCenter.zero();          //the centroid of the points weighted by their drag coefficients
#if ( DIM == 2 )
    real roti = 0;            //in 2D, the total rotational inertia
#endif
    
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
    {
        real R = soRadius[pp];
        if ( R > 0 )
        {
            sumR      += R;
            soDragRot += R * R * R;
            soCenter  += R * posP(pp);
#if ( DIM == 2 )
            roti      += R * posP(pp).normSqr();
#endif
        }
    }
    
    soCenter  /= sumR;
    soDrag     = sumR * 6 * M_PI;
    soDragRot *= 8 * M_PI;
    
    //std::clog << "Solid " << reference() << " drag  " << soDrag << "\n";
    if ( soDrag < REAL_EPSILON )
        throw InvalidParameter("The Solid's drag coefficient is null");
    
#if ( DIM == 2 )
    real m = soDragRot + 6 * M_PI * roti - soDrag * soCenter.normSqr();
    //sanity checks:
    if ( m < REAL_EPSILON )
        throw InvalidParameter("ill-formed Solid has zero rotational drag");
    soMom[0] = 1.0 / m;
#endif
}



/**
 setDragCoefficient() is called by fixShape(), and it is not necessary to
 call it here again.
*/
void Solid::prepareMecable()
{
    //setDragCoefficient();
    
    makeProjection();
}


real Solid::addBrownianForces(real* rhs, real const* rnd, real sc) const
{    
    // Brownian amplitude
    const real drag = prop->viscosity * soDrag;
    real b = sqrt( 2 * sc * drag / nbPoints() );

    for ( unsigned jj = 0; jj < DIM*nbPoints(); ++jj )
        rhs[jj] += b * rnd[jj];
    
    return b / drag;
}

#pragma mark -


#if ( DIM == 1 )

/**
 The projection in 1D is just summing all the forces,
 and distributing equally to all the points:
*/
void Solid::makeProjection()
{
}

void Solid::setSpeedsFromForces(const real* X, const real sc, real* Y, bool) const
{
    real T = 0;
    for ( unsigned p = 0; p < nbPoints(); ++p )
        T += X[p];
    
    T *= sc / ( prop->viscosity * soDrag );
    
    for ( unsigned p = 0; p < nbPoints(); ++p )
        Y[p] = T;
}

#elif ( DIM == 2 )


/**
 Recalculate soCenter, and the rotational moment of inertia.
 */
void Solid::makeProjection()
{
    soCenter = centroid();
    
#if ( 0 )
    /*
     In 2D the rotational moment of inertia is a scalar that is invariant
     by rotation, and it is not normally necessary to recalculate it here
     */
    
    soCenter.zero();
    real roti = 0;
    real sumR = 0;
    
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
    {
        real R = soRadius[pp];
        if ( R > 0 )
        {
            sumR     += R;
            soCenter += R * posP(pp);
            roti     += R * posP(pp).normSqr();
        }
    }
    
    soCenter /= sumR;
    
    real m = soDragRot + 6 * M_PI * roti - soDrag * soCenter.normSqr();
 
    std::clog << "Solid2D " << std::setw(8) << reference() << "  err_roti  " << fabs(soMom[0]-1.0/m) << "\n";

    if ( m < REAL_EPSILON )
        throw InvalidParameter("ill-formed Solid has zero rotational drag");
    
    soMom[0] = 1.0 / m;
#endif
}


void Solid::setSpeedsFromForces(const real* X, const real sc, real* Y, bool) const
{
    real  TX = 0, TY = 0;  //Translation
    real  R  = 0;          //Infinitesimal Rotation (a vector in Z)
    
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
    {
        TX += X[pp*DIM  ];
        TY += X[pp*DIM+1];
        R  += psPos[pp*DIM] * X[pp*DIM+1] - psPos[pp*DIM+1] * X[pp*DIM];
    }
    
    const real alpha = sc * soMom[0] / prop->viscosity;
    const real beta = sc / ( prop->viscosity * soDrag );

    R = alpha * ( R + cross(Vector(TX,TY),soCenter) );
    Vector T = beta * Vector(TX,TY) + cross(soCenter,R);
    
    for ( unsigned p = 0; p < nbPoints(); ++p )
    {
        Y[p*DIM  ] = T.XX - R * psPos[p*DIM+1];
        Y[p*DIM+1] = T.YY + R * psPos[p*DIM  ];
    }
}


#elif ( DIM == 3 )

/**
 To project in 3D, we calculate the resulting tensor by summing all
 the forces on all points, reducing it at the center of gravity.
 From this, we can deduce the forces compatible with solid motion,
 which is a combination of translation and rotation.
 */
void Solid::makeProjection()
{
    soCenter = centroid();
    
    if ( nbPoints() == 1 )
    {
        set_diagonal_matrix(soMom, 1.0 / soDragRot);
        return;
    }

    ///\todo: from reshape, we know the rotation matrix from the stored shape
    //to the current shape. We could use it to transform the inertia matrix
    unsigned cnt = 0;

    real m0=0, m3=0, m6=0, m4=0, m7=0, m8=0;
    
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
    {
        if ( soRadius[pp] > 0 )
        {
            ++cnt;
            const Vector pos = posP(pp);
            real px = soRadius[pp] * pos.XX;
            real py = soRadius[pp] * pos.YY;
            real pz = soRadius[pp] * pos.ZZ;
            m0 += px * pos.XX;
            m3 += px * pos.YY;
            m6 += px * pos.ZZ;
            m4 += py * pos.YY;
            m7 += py * pos.ZZ;
            m8 += pz * pos.ZZ;
        }
    }
    
    if ( cnt == 1 )
    {
        set_diagonal_matrix(soMom, 1.0 / soDragRot);
        return;
    }
    
    // scale to get the correct mobility:
    const real alpha = 6 * M_PI;
    
    // calculate the diagonal term of the matrix:
    const real diag = soDragRot - soDrag*soCenter.normSqr();
    
    // finally set the matrix in front of R in setSpeedsFromForces()
    soMom[0+DIM*0] = diag + alpha * (m4+m8) + soDrag * soCenter[0] * soCenter[0];
    soMom[0+DIM*1] =      - alpha *  m3     + soDrag * soCenter[0] * soCenter[1];
    soMom[0+DIM*2] =      - alpha *  m6     + soDrag * soCenter[0] * soCenter[2];
    soMom[1+DIM*1] = diag + alpha * (m0+m8) + soDrag * soCenter[1] * soCenter[1];
    soMom[1+DIM*2] =      - alpha *  m7     + soDrag * soCenter[1] * soCenter[2];
    soMom[2+DIM*2] = diag + alpha * (m0+m4) + soDrag * soCenter[2] * soCenter[2];
    
    // The matrix should be symmetric positive definite,
    // and we can invert it using the cholesky factorization:
    int info = 0;
    
    lapack_xpotf2('U', DIM, soMom, DIM, &info);
    
    if ( info )
        ABORT_NOW("failed to factorize Solid momentum matrix");
    
    // invert matrix:
    lapack_xpotri('U', DIM, soMom, DIM, &info);
    
    if ( info )
        ABORT_NOW("failed to invert Solid momentum matrix");
    
    // make matrix symmetric:
    soMom[1+DIM*0] = soMom[0+DIM*1];
    soMom[2+DIM*0] = soMom[0+DIM*2];
    soMom[2+DIM*1] = soMom[1+DIM*2];
    
#if ( 0 )
    std::clog << "Solid " << reference() << " " << cnt << "\n";
    Matrix3(soMom).write(std::clog);
#endif
}


/**
 This calculated Y <- P * X, where
 P is the projection associated with the constraints of motion without
 deformation (solid object)
 
 We calculate the total force and momentum in zero, and distribute
 it according to solid motion mechanics.
*/
void Solid::setSpeedsFromForces(const real* X, const real sc, real* Y, bool) const
{    
    real TX=0, TY=0, TZ=0;    //Translation
    real RX=0, RY=0, RZ=0;    //Rotation
    
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
    {
        TX += X[pp*DIM  ];
        TY += X[pp*DIM+1];
        TZ += X[pp*DIM+2];
        RX += psPos[pp*DIM+1] * X[pp*DIM+2] - psPos[pp*DIM+2] * X[pp*DIM+1];
        RY += psPos[pp*DIM+2] * X[pp*DIM  ] - psPos[pp*DIM  ] * X[pp*DIM+2];
        RZ += psPos[pp*DIM  ] * X[pp*DIM+1] - psPos[pp*DIM+1] * X[pp*DIM  ];
    }
    
    Vector V = Vector(RX,RY,RZ) + cross( Vector(TX,TY,TZ), soCenter );
    Vector R;
    
    const real alpha = sc / prop->viscosity;
    const real beta = sc / ( prop->viscosity * soDrag );

    // R = alpha * ( soMom * V ):
    R.XX = alpha * ( soMom[0] * V.XX + soMom[3] * V.YY + soMom[6] * V.ZZ );
    R.YY = alpha * ( soMom[1] * V.XX + soMom[4] * V.YY + soMom[7] * V.ZZ );
    R.ZZ = alpha * ( soMom[2] * V.XX + soMom[5] * V.YY + soMom[8] * V.ZZ );
    
    Vector T = beta * Vector(TX,TY,TZ) + cross(soCenter, R);
    
    for ( unsigned pp = 0; pp < nbPoints(); ++pp )
    {
        Y[pp*DIM  ] = T.XX + R.YY * psPos[pp*DIM+2] - R.ZZ * psPos[pp*DIM+1];
        Y[pp*DIM+1] = T.YY + R.ZZ * psPos[pp*DIM  ] - R.XX * psPos[pp*DIM+2];
        Y[pp*DIM+2] = T.ZZ + R.XX * psPos[pp*DIM+1] - R.YY * psPos[pp*DIM  ];
    }
}

#endif


//------------------------------------------------------------------------------
#pragma mark -


void Solid::write(Outputter& out) const
{
    out.writeUInt16(nbPoints());
    for ( unsigned pp = 0; pp < nbPoints() ; ++pp )
    {
        out.writeFloatVector(psPos + DIM * pp, DIM, '\n');
        out.writeSoftSpace(2);
        out.writeFloat(soRadius[pp]);
    }
}


void Solid::read(Inputter & in, Simul&, Tag)
{
    try {
                     
        unsigned nbp = in.readUInt16();
        setNbPoints(nbp);
        for ( unsigned pp = 0; pp < nbp ; ++pp )
        {
            in.readFloatVector(psPos+DIM*pp, DIM);
            soRadius[pp] = in.readFloat();
        }
        
    }
    catch( Exception & e ) {
        
        e << ", in Solid::read()";
        clearPoints();
        throw;
        
    }
    
    fixShape();
}


void Solid::write(std::ostream& os, bool write_shape) const
{
    std::streamsize p = os.precision();
    os.precision(3);
    os << "new solid " << reference() << '\n';
    os << "{\n";
    os << " nb_points = " << nbPoints() << '\n';
    for ( unsigned n = 0; n < nbPoints() ; ++n )
    {
        os << " point" << n << " = ";
        if ( write_shape )
            os << std::setw(8) << std::fixed << Vector(soShape+DIM*n);
        else
            os << std::setw(8) << std::fixed << Vector(psPos+DIM*n);
        if ( radius(n) > 0 )
            os << ", " << radius(n);
        os << '\n';
    }
    os << "}" << '\n';
    os.precision(p);
}


std::ostream& operator << (std::ostream& os, Solid const& obj)
{
    obj.write(os, false);
    return os;
}


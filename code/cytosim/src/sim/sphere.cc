// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "sim.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "messages.h"
#include "glossary.h"
#include "point_exact.h"
#include "sphere_prop.h"
#include "object_set.h"
#include "space_prop.h"
#include "space.h"
#include "sphere.h"
#include "wrist.h"
#include "meca.h"
#include "modulo.h"
#include "simul.h"

extern Random RNG;

//------------------- construction and destruction ---------------------------

/**
 The Sphere is returned with no points
 */
Sphere::Sphere(SphereProp const* p)
: spRadius(0), spDrag(0), spDragRot(0), spAllocated(0), spProj(0), prop(p)
{
}


/*
 This will create the center point
 */
Sphere::Sphere(SphereProp const* p, real rad)
: spRadius(rad), spDrag(0), spDragRot(0), spAllocated(0), spProj(0), prop(p)
{
    if ( prop == 0 )
        throw InvalidParameter("Sphere:prop should be specified");
    
    if ( rad <= 0 )
        throw InvalidParameter("sphere:radius should be > 0");
    
    // center point
    assert_true( nbPoints() == 0 );
    addPoint( Vector(0,0,0) );
    
    // reference points to track the orientation of the sphere
    if ( DIM >= 2 )
        addPoint( Vector(spRadius,0,0) );
    if ( DIM == 3 ) {
        addPoint( Vector(0,spRadius,0) );
        addPoint( Vector(0,0,spRadius) );
    }
    
    // this only needs to be called once:
    setDragCoefficient();
}



Sphere::Sphere(const Sphere & o)
: PointSet(o)
{
    prop     = o.prop;
    spRadius = o.spRadius;
    setDragCoefficient();
}


Sphere & Sphere::operator =(const Sphere & o)
{
    prop     = o.prop;
    spRadius = o.spRadius;
    setDragCoefficient();
    return *this;
}


Sphere::~Sphere()
{
    //free memory
    if ( spProj ) delete[] spProj;
    prop = 0;
}

//------------------------------------------------------------------------------
#pragma mark -


/*
 if PointSet::allocatePoints() allocated memory, it will return the
 size of the new array, and we allocate the same size for other arrays.
 */
unsigned Sphere::allocatePoints(const unsigned nbp)
{
    unsigned ms = PointSet::allocatePoints(nbp);
    if ( ms )
    {
        //std::clog << "Sphere::allocatePoints " << ms << std::endl;
        allocateProjection(ms);
    }
    return ms;
}


/*
 here 'cp' is the vector from the center to the point to be added,
 in other words, the position of the point in the local reference frame.
 */
unsigned Sphere::addSurfacePoint(Vector const& cp)
{
    return addPoint(posP(0)+cp.normalized(spRadius));
}


/**
 @ingroup NewObject
 
 Specify radius and number of surface points of a Sphere:
 @code
 new sphere NAME
 {
    radius = REAL
    point0 = INTEGER, POSITION [, SINGLE_SPEC]
 }
 @endcode
 
 The `INTEGER` specifies the number of points created, and `POSITION` can be a
 `VECTOR`, or the string 'surface'.  Multiple `SINGLE_SPEC` can be specified.
 
 <h3> Add Singles to a Sphere </h3>
 
 The parameter 'attach' can be used to add Single to the points of a Solid:
 
 @code
 new sphere NAME
 {
    radius   = ...
    point0   = ...
    etc.
    attach   = SINGLE_SPEC [, SINGLE_SPEC] ...
    attach0  = SINGLE_SPEC [, SINGLE_SPEC] ...
    etc.
 }
 @endcode
 
 Where `SINGLE_SPEC` is string containing at most 3 words: `[INTEGER] NAME [each]`,
 where the `INTEGER` specifies the number of Singles, `NAME` specifies their name,
 and the optional word `each` species that the command applies to every point.
 
 The command `attach` applies to all the points of the Solid, while `attach0`,
 `attach1`, etc. apply to the points specified by `point0`, `point1`, etc.
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
ObjectList Sphere::build(Glossary & opt, Simul& simul)
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
        if ( opt.is_number(var) == 2 && opt.set(nbp, var) )
            ++inx;
        
        if ( nbp > 0 )
        {
            unsigned fip = nbPoints();
            // add 'nbp' points:
            for ( unsigned n = 0; n < nbp; ++n )
            {
                Vector vec(0,0,0);
                str = opt.value(var, inx);
                if ( str == "surface" )
                    vec = Vector::randU(radius());
                else
                {
                    std::istringstream iss(str);
                    vec = Movable::readPosition(iss, 0);
                    if ( 8 * vec.norm() < spRadius )
                        throw InvalidParameter(var+" cannot be brought to the Sphere surface");
                }
                addSurfacePoint(vec);
            }
            
            // attach Single to this set of points:
            ++inx;
            while ( opt.set(str, var, inx++) )
                res.append(simul.singles.makeWrists(this, fip, nbp, str));
            
            // attach Single to this set of points:
            inx = 0;
            var = "attach" + sMath::repr(inp);
            while ( opt.set(str, var, inx++) )
                res.append(simul.singles.makeWrists(this, fip, nbp, str));
        }
        
        // set next keyword:
        var = "point" + sMath::repr(++inp);
    }
    
    
    // attach Singles distributed over the surface points:
    inx = 0;
    while ( opt.set(str, "attach", inx++) )
        res.append(simul.singles.makeWrists(this, nbRefPts, nbSurfacePoints(), str));

    
    // final verification of the number of points:
    nbp = 0;
    if ( opt.set(nbp, "nb_points")  &&  nbp != nbPoints() )
    {
        throw InvalidParameter("could not find the number of points specified in solid:nb_points");
    }
    
    //std::cerr << *this << std::endl;
    return res;
}


//------------------------------------------------------------------------------
void Sphere::setInteractions(Meca & meca) const
{
    switch ( prop->confine )
    {
        case CONFINE_OFF:
            break;

        case CONFINE_INSIDE:
        {
            Space const* spc = prop->confine_space_ptr;
            
            Vector cen(psPos);
            if ( ! spc->inside(cen) )
                spc->setInteraction(cen, PointExact(this, 0), meca, prop->confine_stiffness);
        } break;
        
        case CONFINE_ALL_INSIDE:
        {
            const Space* spc = prop->confine_space_ptr;
            
            Vector cen(psPos);
            if ( ! spc->allInside(cen, spRadius) )
                spc->setInteraction(cen, PointExact(this, 0), spRadius, meca, prop->confine_stiffness);
        } break;
        
        case CONFINE_ON:
        {
            const Space* spc = prop->confine_space_ptr;
            spc->setInteraction(posP(0), PointExact(this, 0), meca, prop->confine_stiffness);
        }
            
        default:
            throw InvalidParameter("Invalid sphere::confine");            
    }
}


void Sphere::resize(const real R)
{
    //std::clog << "Sphere::resize " << R << std::endl;
    if ( R > 0 )
    {
        spRadius = R;
        reshape();
        //recalculate drag:
        setDragCoefficient();
    }
}

/**
 the mobility is that of a sphere in an infinite fluid:
 Stokes law:
 
 mu_translation = 6 * PI * viscosity * radius
 dposition/dt   = mu_trans * force
 
 mu_rotation = 8 * PI * viscosity * radius^3
 dangle/dt   = mu_rotation * torque
 */
void Sphere::setDragCoefficientStokes()
{
    assert_true( spRadius > 0 );
    
    const real rad = spRadius;
    
    //hydrodynamic not corrected: infinite fluid is assumed
    spDrag    = 6 * M_PI * prop->viscosity * rad;
    spDragRot = 8 * M_PI * prop->viscosity * rad * rad * rad;

    //MSG("Sphere of radius %.3f has mobility %.2e\n", spRadius, spDrag);
}


/**
 Expect higher friction due to flow around the sphere in a narrow tube.
 This is only valid if (r -a)/a << 1, where r = radius of the tube, and
 a = radius of the sphere.
 
 The formula are taken from:
 <em>The Motion of a Closely-Fitting Sphere in a Fluid-Filled Tube</em>\n
 <b>P. Bungay and H. Brenner, Int. J. Multiphase Flow</b>\n
 Vol 1, pp. 25-56, 1973 (see 3.6, 4.68a and 5.11)
 */
void Sphere::setDragCoefficientPiston()
{
    assert_true( spRadius > 0 );
    assert_true( prop->confine_space_ptr );
    
    const real rad = spRadius;
    real cell_radius = prop->confine_space_ptr->length(1);
    real eps  = ( cell_radius - rad ) / rad;
    
    if ( eps <= 0 )
        throw InvalidParameter("Error: piston formula invalid if sphere is larger than the cell");

    if ( eps > 1 )
        throw InvalidParameter("Error: piston formula invalid if sphere and cylinder do not fit");

    spDrag    = 9*M_PI*M_PI * prop->viscosity * rad / ( 4*sqrt(sMath::power(eps,5)/2) );
    spDragRot = 2*M_PI*M_PI * prop->viscosity * rad * rad * rad / sqrt(eps/2);
        
    //report the reduced mobility of the sphere:
    //MSG("Sphere of radius %.3f has drag coefficient %.2e, due to piston effect\n", spRadius, spDrag);
}


void Sphere::setDragCoefficient()
{
    setDragCoefficientStokes();

    if ( prop->piston_effect )
    {
        if ( prop->confine_space_ptr )
            setDragCoefficientPiston();
        else
            MSG("Piston effect ignored because space is undefined\n");
    }
}


#pragma mark -

void Sphere::prepareMecable()
{
    // setDragCoefficient() was already called by the constructor
    //setDragCoefficient();
    
    assert_true( spDrag > 0 );
    assert_true( spDragRot > 0 );
    
    makeProjection();
}

//------------------------------------------------------------------------------

real Sphere::addBrownianForces(real* rhs, real const* rnd, real sc) const
{
    real bT = sqrt( 2 * sc * spDrag );
    real bS = prop->point_mobility > 0 ? sqrt( 2 * sc / prop->point_mobility ) : 0;

    Vector F(0, 0, 0);
    Torque T(nullTorque);

    real cx = psPos[0];
    real cy = psPos[1];
    real cz = psPos[2];

    /*
     Add random forces to the surface points, and calculate the resulting force
     and momentum in F and T. They will be subtracted from the reference points.
     */
    for ( unsigned dp = DIM*nbRefPts; dp < DIM*nbPoints(); dp+=DIM )
    {
        Vector fp = bS * Vector(rnd+dp);
        
        F += fp;
        
        rhs[dp  ] += fp.XX;
        
#if   ( DIM == 2 )
        rhs[dp+1] += fp.YY;
        T += cross(Vector(psPos[dp]-cx, psPos[dp+1]-cy), fp);
#elif ( DIM == 3 )
        rhs[dp+1] += fp.YY;
        rhs[dp+2] += fp.ZZ;
        T += cross(Vector(psPos[dp]-cx, psPos[dp+1]-cy, psPos[dp+2]-cz), fp);
#endif
    }

    /*
     The Torque is distributed to the surface points.
     In 2D, there is one point, and the coefficient is therefore 1.
     in 3D, there are 3 points, but always one is parallel to the axis of the torque,
     and the decomposition over these 3 points gives a factor 2.
     */
    T /= - ( DIM - 1 ) * spRadius * spRadius;
    Vector R = cross(Vector(cx,cy,cz), T);

    for ( unsigned dp = DIM; dp < DIM*nbRefPts; dp+=DIM )
    {
        Vector fp = bT * Vector(rnd+dp);
#if   ( DIM == 2 )
        rhs[dp]   += R.XX - T * psPos[dp+1] + fp.XX;
        rhs[dp+1] += R.YY + T * psPos[dp  ] + fp.YY;
        F += fp + cross(T, Vector(psPos[dp]-cx, psPos[dp+1]-cy));
#elif ( DIM == 3 )
        rhs[dp  ] += R.XX + T.YY * psPos[dp+2] - T.ZZ * psPos[dp+1] + fp.XX;
        rhs[dp+1] += R.YY + T.ZZ * psPos[dp  ] - T.XX * psPos[dp+2] + fp.YY;
        rhs[dp+2] += R.ZZ + T.XX * psPos[dp+1] - T.YY * psPos[dp  ] + fp.ZZ;
        F += fp + cross(T, Vector(psPos[dp]-cx, psPos[dp+1]-cy, psPos[dp+2]-cz));
#endif
    }
    
    // center of the sphere:
#if   ( DIM == 2 )
    rhs[0] -= F.XX + bT * rnd[0];
    rhs[1] -= F.YY + bT * rnd[1];
#elif ( DIM == 3 )
    rhs[0] -= F.XX + bT * rnd[0];
    rhs[1] -= F.YY + bT * rnd[1];
    rhs[2] -= F.ZZ + bT * rnd[2];
#endif

    return std::max(bT/spDrag, bS*prop->point_mobility);
}


void Sphere::orthogonalizeRef(unsigned i)
{
#if ( DIM == 3 )
    const unsigned ix = 1 + i;
    const unsigned iy = 1 + (i+1)%3;
    const unsigned iz = 1 + (i+2)%3;
    
    Vector cen(psPos);
    assert_true( nbPoints() >= nbRefPts );
    
    // reduce to the center of mass an normalize
    Vector tmpX =   posP(ix) - cen;
    Vector tmpY =   posP(iy) - cen;
    Vector tmpZ = ( posP(iz) - cen ).normalized();
    
    // make tmpY orthogonal to tmpZ, and normalized
    tmpY -= (tmpZ*tmpY) * tmpZ;
    tmpY.normalize();
    
    // make tmpX orthogonal to tmpZ and tmpY
    tmpX -= (tmpZ*tmpX) * tmpZ + (tmpY*tmpX) * tmpY;
    tmpX.normalize();
    
    // store corrected vectors back into the array
    ( cen + spRadius * tmpX ).put(psPos+DIM*ix);
    ( cen + spRadius * tmpY ).put(psPos+DIM*iy);
    ( cen + spRadius * tmpZ ).put(psPos+DIM*iz);
#endif
}


/**
 we get rid of finite-step errors but conserve the shape
 by projecting back onto the sphere,
 without changing the position of point zero (the center)
*/
void Sphere::reshape()
{
    assert_true( nbPoints() > 0 );
    assert_true( spRadius > 0 );
    Vector axis;
    Vector cen(psPos);
    
    for ( unsigned j = 1; j < nbPoints(); ++j )
    {
        axis = ( posP(j) - cen ).normalized(spRadius);
        setPoint(j, cen + axis);
    }
    
#if ( DIM == 3 )
    orthogonalizeRef(RNG.pint(3));
#endif
}


//------------------------------------------------------------------------------
//------------------- methods for the projection -------------------------------
#pragma mark -


void Sphere::allocateProjection(const unsigned int nbp)
{
    //std::clog << "Sphere::allocateProjection(" << nbp << ")" << std::endl;
    if ( spAllocated < nbp )
    {        
        if ( spProj ) delete[] spProj;
        
        spAllocated = nbp;
        spProj = new real[DIM*spAllocated];
    }
}


#if (DIM == 1)

//this is unsafe, don't use the sphere in 1D!
void Sphere::makeProjection() { ABORT_NOW("Sphere is not implemented in 1D"); }
void Sphere::addSurfaceSpeedsFromForces(const real*, real, real*) const {}
void Sphere::setSpeedsFromForces(const real* X, real, real* Y, bool) const {}

#elif ( DIM == 2 || DIM == 3 )

/**
 prepare variables for the projection 
 */
void Sphere::makeProjection()
{
    //allocate more memory if needed
    allocateProjection(nbPoints());
    assert_true( spAllocated >= nbPoints() );
    assert_true( nbPoints() >= nbRefPts );

    //preparations for the motion of the Surfaces:
    //the reference points will be omitted!
        
    //calculate the nonzero components of J, omitting factors of 2
    real curv = 1.0 / spRadius;
    for ( unsigned dp = DIM*nbRefPts; dp < DIM*nbPoints(); dp+=DIM)
    {
        spProj[dp  ] = curv * ( psPos[dp  ] - psPos[0] );
#if ( DIM >= 2 )
        spProj[dp+1] = curv * ( psPos[dp+1] - psPos[1] );
#endif
#if ( DIM == 3 )
        spProj[dp+2] = curv * ( psPos[dp+2] - psPos[2] );
#endif
    }
}


void Sphere::setSphereSpeedsFromForces(const real* X, const real sc, real* Y) const
{
    // total force:
    Vector F(0,0,0);
    
    // total torque:
#if   ( DIM == 2 )
    real T = 0;
#elif ( DIM == 3 )
    Vector T(0,0,0);
#endif
    
    for ( unsigned dp = 0; dp < DIM*nbPoints(); dp+=DIM )
    {
        F.XX +=  X[dp  ];
        F.YY +=  X[dp+1];
#if   ( DIM == 2 )
        T    += psPos[dp] * X[dp+1] - psPos[dp+1] * X[dp];
#elif ( DIM == 3 )
        F.ZZ +=  X[dp+2];
        T.XX += psPos[dp+1] * X[dp+2] - psPos[dp+2] * X[dp+1];
        T.YY += psPos[dp+2] * X[dp  ] - psPos[dp  ] * X[dp+2];
        T.ZZ += psPos[dp  ] * X[dp+1] - psPos[dp+1] * X[dp  ];
#endif
    }
    
    Vector cen(psPos);

    T -= cross(cen, F);       // reduce the torque to the center of mass
    T *= sc/spDragRot;          // multiply by the mobility and maybe time_step
    F  = F*(sc/spDrag) + cross(cen, T);
    
    for ( unsigned dp = 0; dp < DIM*nbPoints(); dp+=DIM )
    {
#if   ( DIM == 2 )
        Y[dp]   = F.XX - T * psPos[dp+1];
        Y[dp+1] = F.YY + T * psPos[dp];
#elif ( DIM == 3 )
        Y[dp  ] = F.XX + T.YY * psPos[dp+2] - T.ZZ * psPos[dp+1];
        Y[dp+1] = F.YY + T.ZZ * psPos[dp  ] - T.XX * psPos[dp+2];
        Y[dp+2] = F.ZZ + T.XX * psPos[dp+1] - T.YY * psPos[dp  ];
#endif
    }
}


void Sphere::addSurfaceSpeedsFromForces(const real* X, real sc, real* Y) const
{
    //scale by point mobility:
    sc *= prop->point_mobility;
    
    // no surface-motion on the center and reference points
    assert_true( nbPoints() >= nbRefPts );
    
    // skip the reference-points:
    for ( unsigned dp = DIM*nbRefPts; dp < DIM*nbPoints(); dp+=DIM )
    {
        real a = spProj[dp] * X[dp];
        
        for ( int d = 1; d < DIM; ++d )
            a += spProj[dp+d] * X[dp+d];
        
        for ( int d = 0; d < DIM; ++d )
            Y[dp+d] += sc * ( X[dp+d] - a * spProj[dp+d] );
    }
}


void Sphere::setSpeedsFromForces(const real* X, const real sc, real* Y, bool) const
{
    // set motions from the rigid body
    setSphereSpeedsFromForces( X, sc, Y );
    //blas_xzero(DIM*nbPoints(), Y);
    
    // add contribution from the Surface motion
    addSurfaceSpeedsFromForces( X, sc, Y );
}
#endif



//------------------------------------------------------------------------------
#pragma mark -

void Sphere::write(Outputter& out) const
{
    out.writeFloat(radius());
    PointSet::write(out);
}


void Sphere::read(Inputter & in, Simul& sim, Tag tag)
{
    try {
        real rad;
#ifdef BACKWARD_COMPATIBILITY
        if ( in.formatID() < 36 )
            rad = radius();
        else
#endif
        rad = in.readFloat();
        PointSet::read(in, sim, tag);
        resize(rad);
    }
    catch( Exception & e ) {
        e << ", in Sphere::read()";
        clearPoints();
        throw;
    }
}



std::ostream& operator << (std::ostream& os, Sphere const& obj)
{
    os << "new sphere " << obj.reference() << '\n';
    os << "{\n";
    os << " nb_points = " << obj.nbPoints() << '\n';
    for ( unsigned pp = 0; pp < obj.nbPoints() ; ++pp )
        os << " point" << pp << " = " << obj.posP(pp) << '\n';
    os << "}" << '\n';
    return os;
}


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "messages.h"
#include "aster.h"
#include "solid.h"
#include "solid_prop.h"
#include "fiber_prop.h"
#include "pointsonsphere.h"
#include "point_exact.h"
#include "point_interpolated.h"
#include "glossary.h"
#include "simul.h"
#include "meca.h"

extern Random RNG;


void Aster::step()
{
    assert_true( linked() );
    assert_true( asLinks.size()+1 == nbOrganized() );
    
    Simul & sim = objset()->simul;

    // nucleation:
    for ( unsigned ii = 0; ii < asLinks.size(); ++ii )
    {
        if ( 0 == fiber(ii) &&  RNG.test(prop->nucleation_rate_prob) )
        {
            Glossary opt;
            sim.add(makeFiber(ii, opt, sim));
        }
    }
}

#if   ( DIM == 1 )
#    define INTERLINK interLink2
#elif ( DIM == 2 )
#    define INTERLINK interLink3
#else
#    define INTERLINK interLink4
#endif

/*
 Note on possible optimization:
 The coefficients of the interpolations to the Solid points are constant in time,
 and so we could simply set a matrix once, and keep it over time.
 Specifically, we would introduce a new matrix in Meca, `mK` and set it only once.
 We can then include these additional terms directly in Meca::addLinearForces():
 @code
 Y = Y + ( mB + mC + mK ) * X
 @endcode
 */
void Aster::setInteractions(Meca & meca) const
{
    assert_true( linked() );
    assert_true( asLinks.size()+1 == nbOrganized() );

    Solid const* sol = solid();
    
    if ( sol == 0 )
        return;
    
    unsigned pts[] = { 0, 1, 2, 3 };

    for ( unsigned n = 0 ; n < asLinks.size(); ++n )
    {
        Fiber * fib = fiber(n);

        if ( fib )
        {
            AsterLink const& link = asLinks[n];
            unsigned off = sol->matIndex() + link.ref;
            
#ifdef BACKWARD_COMPATIBILITY
            if ( link.alt > 0 )
            {
                meca.interLink(PointExact(sol, link.ref), fib->exactEnd(prop->focus), prop->stiffness[0]);
               if ( fib->length() > link.len )
                {
                    meca.interLink(fib->interpolate(link.len, prop->focus), PointExact(sol, link.alt), prop->stiffness[1]);
                }
                else
                {
                    FiberEnd tip = ( prop->focus == PLUS_END ? MINUS_END : PLUS_END );
                    // link the opposite end to an interpolation of the two solid-points:
                    real c = fib->length() / link.len;
                    meca.interLink(PointInterpolated(sol, link.ref, link.alt, c), fib->exactEnd(tip), prop->stiffness[1]);
                }
                continue;
            }
#endif
            if ( link.degenerate )
                meca.interLink(fib->exactEnd(prop->focus), PointExact(sol, link.ref), prop->stiffness[0]);
            else
                meca.INTERLINK(fib->exactEnd(prop->focus), off, pts, link.coef1, prop->stiffness[0]);
            
            real len = link.len;
            
            if ( fib->length() >= len )
            {
                if ( len > 0 )
                    meca.INTERLINK(fib->interpolate(len, prop->focus), off, pts, link.coef2, prop->stiffness[1]);
                else
                    meca.INTERLINK(fib->exactEnd(prop->focus), off, pts, link.coef2, prop->stiffness[1]);
            }
            else
            {
                // link the opposite fiber end to a new interpolation:
                FiberEnd end = ( prop->focus == PLUS_END ? MINUS_END : PLUS_END );
                real c = fib->length() / len;
                real u = 1.0 - c;
                real coef[4];
                for ( int d = 0; d < 4; ++d )
                    coef[d] = u * link.coef1[d] + c * link.coef2[d];
                meca.INTERLINK(fib->exactEnd(end), off, pts, coef, prop->stiffness[1]);
            }
        }
    }
}


Aster::~Aster()
{
    //MSG(31, "destroying %c%lu\n", TAG, identity() );
    prop = 0;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 @defgroup NewAster How to create an Aster
 @ingroup NewObject
 
 By default the aster creates a radial distribution of fiber,
 and only the radius need to be specified:
 
 @code
 new aster NAME
 {
   nb_fibers = INTEGER
   radius = OUTER_RADIUS, INNER_RADIUS
   ...
 }
 @endcode
 
 The configuration of the Aster can also be customized by providing
 points on which the fibers are attached.
 
 Syntax:
 @code
 new aster NAME
 {
   nb_fibers = INTEGER
   anchor0 = POINT, POINT
   anchor1 = POINT, POINT
   ...
 }
 @endcode
 
 Each POINT can be specified in two ways:
 - as a VECTOR, to specify the position of a new point,
 - as `point0`, `point1`, etc. to refer to an already defined point.
 .
 A point with `index=0` is always added at the center of the Aster.
 
 Example:
 
 @code
 new aster centrosome
 {
   nb_fibers = 3
   anchor0 = point0, 1  0 0
   anchor1 = point0, 0  1 0
   anchor2 = point0, 0 -1 0
 }
 @endcode
 
 It is also possible to add point to the Solid, before adding the anchors.
 The syntax is the same as for customizing a solid (@ref Solid::build):
 
 @code
 new aster centrosome
 {
   point0 = 0 -0.2 0, 0.2
   point1 = 0  0   0, 0.2
   point2 = 0 +0.2 0, 0.2

   nb_fibers = 3
   anchor0 = point0,  0.5 -0.2 0
   anchor1 = point1,  0.5  0   0
   anchor2 = point2,  0.5 +0.2 0
 }
 @endcode
 
 */
ObjectList Aster::build(Glossary& opt, Simul& simul)
{
    assert_true(prop);
    assert_true(nbOrganized()==0);

    opt.set(asRadius, "radius");
        
    if ( asRadius <= 0 )
        throw InvalidParameter("aster:radius must be specified and > 0");

    unsigned origin = 0;
    ObjectList res = makeSolid(opt, simul, origin);
    
    if ( nbOrganized() != 1 )
        throw InvalidParameter("could not make aster:solid");

    placeFibers(opt, origin);
    
    //solid()->write(std::clog);
    
    for ( unsigned n = 0; n < asLinks.size(); ++n )
        res.append(makeFiber(n, opt, simul));
    
    return res;
}


ObjectList Aster::makeFiber(unsigned inx, Glossary& opt, Simul& simul)
{
    ObjectList res = simul.fibers.newObjects(prop->fibers, opt);
    
    if ( res.empty() )
        throw InvalidParameter("could not create aster:fiber");

    Fiber * fib = Fiber::toFiber(res[0]);

    if ( !fib )
        throw InvalidParameter("unexpected object returned by fibers.newObjects()");

    grasp(fib, inx+1);

    Vector pos = posLink1(inx);
    
    Vector dir = posLink2(inx) - pos;
    real n = dir.normSqr();
    if ( n > REAL_EPSILON )
    {
        if ( prop->focus == PLUS_END )
            dir /= -sqrt(n);
        else
            dir /= sqrt(n);
    }
    else
        dir = Vector::randU();
    
    ObjectSet::rotateObjects(res, Rotation::rotationToVector(dir, RNG));
    ObjectSet::translateObjects(res, pos - fib->posEnd(prop->focus));
    
    return res;
}


/**
 Anchor Fiber stored at index `inx` between positions A and B
 Here the positions A and B are given in a local reference frame, 
 made of unit vectors which have a length 'asRadius'.
 Thus everything will be scaled by 'asRadius' eventually.
 */
void Aster::anchorFiber(unsigned inx, Vector const& A, Vector const& B, unsigned ref)
{
    assert_true( inx < asLinks.size() );
    asLinks[inx].set(A, B);
    asLinks[inx].len *= asRadius;
    asLinks[inx].ref = ref;
    //asLinks[inx].dump(std::clog);
}


ObjectList Aster::makeSolid(Glossary& opt, Simul& simul, unsigned& origin)
{
    ObjectList res;
    Solid * sol = 0;
    
    // find the Solid specified:
    std::string spec;
    if ( opt.set(spec, "solid") )
    {
        ObjectList objs = simul.solids.findObjects(spec);
        
        if ( objs.size() != 1 )
            throw InvalidParameter("could not find aster:solid `"+spec+"'");
        
        if ( objs[0]->tag() != Solid::TAG )
            throw InvalidParameter("object `"+spec+"' is not a Solid");
        
        // we found a Solid
        sol = static_cast<Solid*>(objs[0]);
        
        if ( sol->nbBuddies() > 0 )
            throw InvalidParameter("aster:solid `"+spec+"' is already engaged and thus unsuitable");
        
        // add local coordinate system:
        origin = sol->addTriad(asRadius);
    }
    
#ifdef BACKWARD_COMPATIBILITY
    // make a solid from scratch
    if ( sol == 0 )
    {
        Property * p = simul.properties.find_or_die("solid", prop->solid);
        sol = new Solid(static_cast<SolidProp*>(p));
        res.append(sol->build(opt, simul));
        res.push_back(sol);
        // add a massive bead if needed:
        if ( sol->dragCoefficient() < REAL_EPSILON )
            sol->addSphere(Vector(0,0,0), asRadius);
        // add local coordinate system:
        origin = sol->addTriad(asRadius);
        //std::clog << "Aster::makeSolid() created solid " << sol->reference() << "\n";
    }
#endif

    sol->fixShape();

    // check that there is at least one point:
    if ( sol->dragCoefficient() < REAL_EPSILON )
        throw InvalidParameter("Aster's drag coefficient is null: please specify 'point0=center, RADIUS'");
    
    grasp(sol, 0);

    return res;
}


/**
 Generate a random distribution of points on the unit disc,
 with the distance between two points never below `sep`.
 */
unsigned tossPointsDisc(unsigned nbp, Vector2 pts[], real sep, unsigned max_nb_trials)
{
    const real ss = sep * sep;
    unsigned ouf = 0;
    unsigned n = 0;
    
    while ( n < nbp )
    {
    toss:
        if ( ++ouf > max_nb_trials )
            break;
        
        Vector2 xy = Vector2::randB();

        for ( int i = 0; i < n; ++i )
            if ( xy.distanceSqr(pts[i]) < ss )
                goto toss;
        
        pts[n++] = xy;
        ouf = 0;
    }
    return n;
}


/**
 Generate a random distribution of points on the unit circle,
 with the distance between two points never below `sep`.
 */
unsigned tossPointsCap(unsigned nbp, Vector pts[], real cap, real sep, unsigned max_nb_trials)
{
    const real ss = sep * sep;
    unsigned ouf = 0;
    unsigned n = 0;
    
    while ( n < nbp )
    {
    toss:
        if ( ++ouf > max_nb_trials )
            break;
        
        real a = M_PI * RNG.sreal();
        real u = 1.0 - cap * RNG.preal();
        real v = sqrt( 1.0 - u * u );
        Vector pos(u, v*cos(a), v*sin(a));
        
        for ( int i = 0; i < n; ++i )
            if ( pos.distanceSqr(pts[i]) < ss )
                goto toss;
        
        pts[n++] = pos;
        ouf = 0;
    }
    return n;
}

/**
 Generate a random distribution of points on the unit circle,
 with the distance between two points never below `sep`.
 */
unsigned tossPoints(unsigned nbp, Vector pts[], real sep, unsigned max_nb_trials)
{
    const real ss = sep * sep;
    unsigned ouf = 0;
    unsigned n = 0;
    
    while ( n < nbp )
    {
    toss:
        if ( ++ouf > max_nb_trials )
            break;
        
        Vector pos = Vector::randU();
        
        for ( int i = 0; i < n; ++i )
            if ( pos.distanceSqr(pts[i]) < ss )
                goto toss;
        
        pts[n++] = pos;
        ouf = 0;
    }
    return n;
}


/**
 One can specify the `radius` of the aster, and `nb_fibers`.
 
 The aster 'type' can be:
 - `astral` fiberd are anchored at random positions near the center, pointing outward
 - `radial` fibers are anchored always at the same distance from the center, pointing radially
 - `regular` fibers are anchored regularly over the surface and point radially
 - `angular` only
 .
 */
void Aster::placeFibers(Glossary & opt, unsigned origin)
{
    unsigned nbf = 0;
    opt.set(nbf, "nb_fibers");

    real dis = 0;
    if ( opt.set(dis, "radius", 1) && dis > asRadius )
        throw InvalidParameter("aster:radius[1] must be smaller than aster:radius[0]");

    const real alpha = dis / asRadius;
    
    unsigned type = 0;
    opt.set(type, "type", KeyList<unsigned>("radial", 0, "astral", 1, "regular", 2, "angular", 3, "disc", 4));
    
    if ( type == 4 )
    {
        // This is a special case for Yeast spindles
        // use a separation of 40 nm by default, corresponding to Microtubules
        real sep = 0.04;
        opt.set(sep, "fiber_separation");
        Vector2 * pts = new Vector2[nbf];
        unsigned ouf = 0;
        unsigned cnt = 0;
        do {
            cnt = tossPointsDisc(nbf, pts, sep/asRadius, 1<<10);
        } while ( cnt < nbf && ++ouf < 1<<10 );
        //std::clog << "toss(" << nbf << ") placed " << cnt << "\n";
        asLinks.resize(cnt);
        real d = dis * 0.5;
        for ( unsigned n = 0; n < cnt; ++n )
        {
            // orient anchors by default along the X-axis:
            real y = pts[n].YY;
#if ( DIM == 2 )
            anchorFiber(n, Vector2(-d,y), Vector2(+d,y), origin);
#elif ( DIM == 3 )
            real x = pts[n].XX;
            anchorFiber(n, Vector3(-d,x,y), Vector3(+d,x,y), origin);
#endif
        }
        delete[] pts;
    }
    else if ( type == 3 )
    {
        /*
         For type 'angular' all fibers are restricted within an specified solid angle,
         and their orientation is radial
         by GAELLE LETORT, 14.03.2017
         */
        real angle = M_PI;
        opt.set(angle, "angle") || opt.set(angle, "aster_angle");
        if ( angle < REAL_EPSILON )
            throw InvalidParameter("aster::angle must be > 0");
        asLinks.resize(nbf);
#if ( DIM == 1 )
        // No effect of angle in 1D, same as default
        for ( unsigned n = 0; n < nbf; ++n )
        {
            Vector D(n%2?1:-1, 0);
            anchorFiber(n, Vector(0, 0), D, origin);
        }
#elif ( DIM == 2 )
        real delta = 2 * angle / real(nbf);
        // points are evenly distributed from -aster_angle to aster_angle
        real ang = -angle;
        for ( unsigned n = 0; n < nbf; ++n )
        {
            Vector P(cos(ang), sin(ang));
            anchorFiber(n, alpha*P, P, origin);
            ang += delta;
        }
#else
        real cap = 1.0 - cos(angle);
        real sep, sep0 = sqrt( 2 * M_PI * cap / nbf );
        Vector * pts = new Vector[nbf];
        unsigned ouf = 0;
        unsigned cnt = 0;
        do {
            sep = 512 * sep0 / real(ouf+512);
            cnt = tossPointsCap(nbf, pts, cap, sep, 1024);
            //std::clog << "toss(" << nbf << ") placed " << cnt << " with sep = " << sep << "\n";
        } while ( cnt < nbf && ++ouf < 1024 );
        //std::clog << "toss(" << nbf << ") placed " << cnt << " with sep = " << sep << "\n";
        asLinks.resize(cnt);
        for ( unsigned n = 0; n < cnt; ++n )
            anchorFiber(n, alpha*pts[n], pts[n], origin);
        delete[] pts;
#endif
    }
    else if ( type == 2 )
    {
        /*
         For type 'regular' we put fibers regularly on the surface,
         */
        asLinks.resize(nbf);
#if ( DIM == 1 )
        for ( unsigned n = 0; n < nbf; ++n )
        {
            Vector D(n%2?1:-1, 0);
            anchorFiber(n, Vector(0, 0), D, origin);
        }
#elif ( DIM == 2 )
        real ang = 0, delta = 2 * M_PI / real(nbf);
        for ( unsigned n = 0; n < nbf; ++n )
        {
            Vector P(cos(ang), sin(ang));
            anchorFiber(n, alpha*P, P, origin);
            ang += delta;
        }
#else
        //we use PointOnSphere to distribute points on the sphere
        PointsOnSphere sphere(RNG, nbf);
        Vector P;
        for ( unsigned n = 0; n < nbf; ++n )
        {
            sphere.copyPoint(P, n);
            anchorFiber(n, alpha*P, P, origin);
        }
#endif
    }
    else if ( type == 1 )
    {
        /*
         For type 'astral' we put fibers randomly on the surface,
         with a constrain based on the scalar product: position*direction > 0
         */
        if ( dis <= 0 )
            throw InvalidParameter("aster:radius[1] must be specified and >= 0");
        asLinks.resize(nbf);
        for ( unsigned n = 0; n < nbf; ++n )
        {
            Vector P = Vector::randB();
            Vector D = Vector::randU();
            while ( D * P < 0 )
                D = Vector::randU();
            anchorFiber(n, P-alpha*D, P, origin);
        }
    }
    else if ( type == 0 )
    {
        /*
         For type 'radial' we put fibers randomly on the surface,
          and set their direction as purely radial.
         */
        // use a separation of 40 nm by default, corresponding to Microtubules
        real sep = 0.04;
        opt.set(sep, "fiber_separation");
        
        Vector * pts = new Vector[nbf];
        unsigned ouf = 0;
        unsigned cnt = 0;
        do {
            cnt = tossPoints(nbf, pts, sep/asRadius, 1<<10);
        } while ( cnt < nbf && ++ouf < 1<<10 );
        //std::clog << "toss(" << nbf << ") placed " << cnt << "\n";
        asLinks.resize(cnt);
        for ( unsigned n = 0; n < cnt; ++n )
            anchorFiber(n, alpha*pts[n], pts[n], origin);
        delete[] pts;
    }
    else
        throw InvalidParameter("unknown aster::type");
}


//------------------------------------------------------------------------------
#pragma mark -

void Aster::write(Outputter& out) const
{
    Organizer::write(out);
    
    out.writeSoftNewline();
    for ( unsigned ii = 0; ii < asLinks.size(); ++ii )
    {
        out.writeSoftNewline();
        asLinks[ii].write(out);
    }
}


void Aster::read(Inputter & in, Simul& sim, Tag tag)
{
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 40 )
        in.readUInt16();
#endif
    
    Organizer::read(in, sim, tag);
    
    assert_true( nbOrganized() > 0 );
    assert_true( organized(0)->tag() == Solid::TAG );
    
    Solid * sol = solid();
    
    asRadius = ( sol->posPoint(0) - sol->posPoint(1) ).norm();
    
    try
    {
        unsigned nc = nbOrganized() - 1;
        asLinks.resize(nc);
        for ( unsigned inx = 0; inx < nc; ++inx )
        {
#ifdef BACKWARD_COMPATIBILITY
            if ( in.formatID() < 47 )
            {
                asLinks[inx].reset();
                unsigned a = in.readUInt16();
                unsigned b = in.readUInt16();
                if ( a >= sol->nbPoints() || b >= sol->nbPoints() )
                    throw InvalidIO("invalid point index");
                asLinks[inx].ref = a;
                asLinks[inx].coef1[0] = 1.0;
                asLinks[inx].alt = b;
                asLinks[inx].len = (sol->posPoint(a)-sol->posPoint(b)).norm();
                continue;
            }
#endif
            asLinks[inx].read(in);
            asLinks[inx].len *= asRadius;
        }
        
        if ( nc > 0 )
        {
            unsigned ref = asLinks[0].ref;
            asRadius = ( sol->posPoint(ref) - sol->posPoint(ref) ).norm();
        }
    }
    catch( Exception & e )
    {
        e << ", in Aster::read()";
        throw;
    }
}


//------------------------------------------------------------------------------
#pragma mark -

Vector Aster::posLink1(unsigned inx) const
{
    Solid const* sol = solid();
    real const* coef = asLinks[inx].coef1;
    const unsigned ref = asLinks[inx].ref;
    
#ifdef BACKWARD_COMPATIBILITY
    if ( asLinks[inx].alt > 0 )
        return sol->posPoint(ref);
#endif

    Vector res = coef[0] * sol->posPoint(ref);
    for ( int i = 1; i <= DIM; ++i )
        res += coef[i] * sol->posPoint(i+ref);
    
    return res;
}

Vector Aster::posLink2(unsigned inx) const
{
    Solid const* sol = solid();
    real const* coef = asLinks[inx].coef2;
    const unsigned ref = asLinks[inx].ref;
    
#ifdef BACKWARD_COMPATIBILITY
    if ( asLinks[inx].alt > 0 )
        return sol->posPoint(asLinks[inx].alt);
#endif

    Vector res = coef[0] * sol->posPoint(ref);
    for ( int i = 1; i <= DIM; ++i )
        res += coef[i] * sol->posPoint(i+ref);
    
    return res;
}

Vector Aster::posFiber2(unsigned inx) const
{
    Fiber const* fib = fiber(inx);
    real len = asLinks[inx].len;
    
    if ( fib->length() >= len )
    {
        if ( len > 0 )
            return fib->pos(len, prop->focus);
        else
            return fib->posEnd(prop->focus);
    }
    else
    {
        // link the opposite end to an interpolation of the two solid-points:
        return fib->posEnd( prop->focus == PLUS_END ? MINUS_END : PLUS_END );
    }
}

/**
 This sets the ends of the link number `inx`
 or returns zero if the link does not exist
 */
unsigned Aster::getLink(unsigned inx, Vector& pos1, Vector& pos2) const
{
    unsigned n = inx / 2;
    
    if ( n < asLinks.size() )
    {
        Fiber const* fib = fiber(n);
        
        if ( inx & 1 )
        {
            pos1 = posLink1(n);
            if ( fib )
                pos2 = fib->posEnd(prop->focus);
            else
                pos2 = pos1;
        }
        else
        {
            if ( fib )
            {
                real len = asLinks[n].len;
                if ( fib->length() >= len )
                {
                    pos1 = posLink2(n);
                }
                else
                {
                    // interpolate between the two solid-points:
                    real c = fib->length() / len;
                    pos1 = ( 1.0 - c ) * posLink1(n) + c * posLink2(n);
                }
                pos2 = posFiber2(n);
            }
            else
            {
                pos1 = posLink2(n);
                pos2 = pos1;
            }
        }
        
        return 1;
    }
    
    return 0;
}



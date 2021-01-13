// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_force.h"
#include "exceptions.h"
#include "point_exact.h"
#include "glossary.h"
#include "meca.h"


SpaceForce::SpaceForce(const SpaceProp* p, Glossary& opt)
: Space(p)
{
    force.zero();
    center.zero();
    stiffness = 0;
    
    opt.set(force, "force");
    opt.set(center, "center");
    opt.set(stiffness, "stiffness");
}


void SpaceForce::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-1, -1, -1);
    sup.set( 1,  1,  1);
}


real SpaceForce::volume() const
{
    throw InvalidParameter("invalid use of space `force'");
    return -1;
}


void SpaceForce::project( const real w[], real p[] ) const
{
    throw InvalidParameter("Invalid use of space `force'");
}


void SpaceForce::setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff) const
{
    throw InvalidParameter("Invalid use of space `force'");
}


void SpaceForce::setInteraction(Vector const& pos, PointExact const&, real rad, Meca &, real stiff) const
{
    throw InvalidParameter("Invalid use of space `force'");
}


void SpaceForce::setInteractions(Meca & meca, FiberSet const&) const
{
    real *const base = meca.base();

    if ( stiffness > 0 )
    {
        Vector sc = stiffness * center;
        const unsigned nbp = meca.nbPoints();
        // generate an isotropic squeezing force:
        for ( int p = 0; p < nbp; ++p )
        {
            meca.mB(p,p)  -= stiffness;
            sc.add_to(base+DIM*p);
        }
    }
    else
    {
        const unsigned nbu = meca.size();
        for ( unsigned u = 0; u < nbu; u += DIM )
            force.add_to(base+u);
    }
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpaceForce::display() const
{
    Vector vec(length(0), length(1), length(2));
    real d = 0.1 * vec.norm();
    gleArrow(Vector(0,0,0), vec, d);
    return true;
}

#else

bool SpaceForce::display() const
{
    return false;
}


#endif



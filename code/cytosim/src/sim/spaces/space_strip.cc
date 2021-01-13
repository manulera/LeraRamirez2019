// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_strip.h"
#include "exceptions.h"
#include "point_exact.h"
#include "meca.h"


SpaceStrip::SpaceStrip(const SpaceProp* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("strip is not usable in 1D");
}


void SpaceStrip::setModulo(Modulo * mod) const
{
    for ( int d = 0; d < DIM-1; ++d )
        mod->enable(d, length(d));
}


void SpaceStrip::resize()
{
    checkLengths(DIM-1, true);
    
    if ( length(DIM-1) < 0 )
        throw InvalidParameter("strip:dimension[DIM-1] must be >= 0");
}


void SpaceStrip::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length(0),-length(1),-length(2));
    sup.set( length(0), length(1), length(2));
}

//------------------------------------------------------------------------------
#pragma mark -


#if (DIM == 1)

real SpaceStrip::volume() const
{
    return length2(0);
}

bool  SpaceStrip::inside( const real point[] ) const
{
    if ( point[0] >  length(0) ) return false;
    if ( point[0] < -length(0) ) return false;
    return true;
}


void SpaceStrip::project( const real point[], real proj[] ) const
{
    if ( point[0] > 0 )
        proj[0] =  length(0);
    else
        proj[0] = -length(0);
}

#endif


//------------------------------------------------------------------------------

#if (DIM == 2)

real SpaceStrip::volume() const
{
    return length2(0) * length2(1);
}

bool  SpaceStrip::inside( const real point[] ) const
{
    if ( point[1] >  length(1) ) return false;
    if ( point[1] < -length(1) ) return false;
    return true;
}


void SpaceStrip::project( const real point[], real proj[] ) const
{
    proj[0] = point[0];
    
    if ( point[1] > 0 )
        proj[1] =  length(1);
    else
        proj[1] = -length(1);
}

#endif

//------------------------------------------------------------------------------

#if (DIM == 3)

real SpaceStrip::volume() const
{
    return length2(0) * length2(1) * length2(2);
}

bool  SpaceStrip::inside( const real point[] ) const
{
    if ( point[2] >  length(2) ) return false;
    if ( point[2] < -length(2) ) return false;
    return true;
}

void SpaceStrip::project( const real point[], real proj[] ) const
{
    proj[0] = point[0];
    proj[1] = point[1];
    
    if ( point[2] > 0 )
        proj[2] =  length(2);
    else
        proj[2] = -length(2);
}

#endif

//------------------------------------------------------------------------------

void SpaceStrip::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    unsigned inx = DIM-1 + DIM * pe.matIndex();
    
    meca.mC(inx, inx) -= stiff;

#if ( DIM == 2 )
    if ( pos.YY > 0 )
        meca.base(inx) += stiff * length(1);
    else
        meca.base(inx) -= stiff * length(1);
#elif ( DIM == 3 )
    if ( pos.ZZ > 0 )
        meca.base(inx) += stiff * length(2);
    else
        meca.base(inx) -= stiff * length(2);
#endif
}


void SpaceStrip::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    unsigned inx = DIM-1 + DIM * pe.matIndex();
    
    meca.mC(inx, inx) -= stiff;

#if ( DIM == 2 )
    if ( pos.YY > 0 )
        meca.base(inx) += stiff * ( length(1) - rad );
    else
        meca.base(inx) -= stiff * ( length(1) - rad );
#elif ( DIM == 3 )
    if ( pos.ZZ > 0 )
        meca.base(inx) += stiff * ( length(2) - rad );
    else
        meca.base(inx) -= stiff * ( length(2) - rad );
#endif
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------
#pragma mark -

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpaceStrip::display() const
{
    const real X = length(0);
    const real Y = ( DIM > 1 ) ? length(1) : 1;
    const real Z = ( DIM > 2 ) ? length(2) : 0;
    
#if ( DIM > 2 )
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex( -X,  Y, -Z );
    gleVertex(  X,  Y, -Z );
    gleVertex( -X, -Y, -Z );
    gleVertex(  X, -Y, -Z );
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex( -X,  Y, Z );
    gleVertex( -X, -Y, Z );
    gleVertex(  X,  Y, Z );
    gleVertex(  X, -Y, Z );
    glEnd();
#else
    glBegin(GL_LINES);
    gleVertex( -X,  Y, Z );
    gleVertex(  X,  Y, Z );
    gleVertex( -X, -Y, Z );
    gleVertex(  X, -Y, Z );
    glEnd();
#endif

    glLineStipple(2, 0x0303);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
#if ( DIM > 2 )
    gleVertex(  X,  Y,  Z );
    gleVertex(  X,  Y, -Z );
    gleVertex(  X, -Y,  Z );
    gleVertex(  X, -Y, -Z );
    gleVertex( -X,  Y,  Z );
    gleVertex( -X,  Y, -Z );
    gleVertex( -X, -Y,  Z );
    gleVertex( -X, -Y, -Z );
#else
    gleVertex(  X,  Y, Z );
    gleVertex(  X, -Y, Z );
    gleVertex( -X,  Y, Z );
    gleVertex( -X, -Y, Z );
#endif
    glEnd();
    glDisable(GL_LINE_STIPPLE);

    return true;
}

#else

bool SpaceStrip::display() const
{
    return false;
}

#endif


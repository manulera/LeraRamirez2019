// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_periodic.h"
#include "exceptions.h"


SpacePeriodic::SpacePeriodic(const SpaceProp* p)
: Space(p)
{
}


void SpacePeriodic::setModulo(Modulo * mod) const
{
    for ( int d = 0; d < DIM; ++d )
        mod->enable(d, length(d));
}


void SpacePeriodic::resize()
{
    checkLengths(DIM, true);
}


void SpacePeriodic::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length(0),-length(1),-length(2));
    sup.set( length(0), length(1), length(2));
}

//------------------------------------------------------------------------------
#pragma mark -

#if (DIM == 1)

real SpacePeriodic::volume() const
{
    return length2(0);
}

bool  SpacePeriodic::inside( const real point[] ) const
{
    return true;
}


void SpacePeriodic::project( const real point[], real proj[] ) const
{
    throw InvalidParameter("A periodic space has no edge!");
}

#endif


//------------------------------------------------------------------------------

#if (DIM == 2)

real SpacePeriodic::volume() const
{
    return length2(0) * length2(1);
}

bool  SpacePeriodic::inside( const real point[] ) const
{
    return true;
}


void SpacePeriodic::project( const real point[], real proj[] ) const
{
    throw InvalidParameter("A periodic space has no edge!");
}

#endif

//------------------------------------------------------------------------------

#if (DIM == 3)

real SpacePeriodic::volume() const
{
    return length2(0) * length2(1) * length2(2);
}

bool  SpacePeriodic::inside( const real point[] ) const
{
    return true;
}

void SpacePeriodic::project( const real point[], real proj[] ) const
{
    throw InvalidParameter("A periodic space has no edge!");
}

#endif

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------
#pragma mark -

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpacePeriodic::display() const
{
    const real X = length(0);
    const real Y = ( DIM > 1 ) ? length(1) : 10;
    const real Z = ( DIM > 2 ) ? length(2) : 0;
    
    glLineStipple(2, 0x0303);
    glEnable(GL_LINE_STIPPLE);

#if ( DIM == 1 )
    glBegin(GL_LINES);
    gleVertex(  X, -Y, 0 );
    gleVertex(  X,  Y, 0 );
    gleVertex( -X,  Y, 0 );
    gleVertex( -X, -Y, 0 );
    glEnd();    
#endif
    
#if ( DIM > 1 )
    glBegin(GL_LINE_LOOP);
    gleVertex(  X,  Y, Z );
    gleVertex(  X, -Y, Z );
    gleVertex( -X, -Y, Z );
    gleVertex( -X,  Y, Z );
    glEnd();
#endif

#if ( DIM > 2 )
    glBegin(GL_LINE_LOOP);
    gleVertex(  X,  Y, -Z );
    gleVertex(  X, -Y, -Z );
    gleVertex( -X, -Y, -Z );
    gleVertex( -X,  Y, -Z );
    glEnd();
#endif

    glDisable(GL_LINE_STIPPLE);
    return true;
}

#else

bool SpacePeriodic::display() const
{
    return false;
}

#endif


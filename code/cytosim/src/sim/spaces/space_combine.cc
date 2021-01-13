// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_combine.h"
#include "exceptions.h"


SpaceCombine::SpaceCombine(Space * big, Space * small)
: Space(big->prop)
{
    if ( !big )
        throw InvalidParameter("space:combine invoked with void argument");
    if ( !small )
        throw InvalidParameter("space:combine invoked with void argument");

    outer = big;
    inner = small;
}


SpaceCombine::~SpaceCombine()
{
    if ( inner )
    {
        delete( inner );
        inner = 0;
    }
    if ( outer )
    {
        delete( outer );
        outer = 0;
    }
}


//------------------------------------------------------------------------------
bool  SpaceCombine::inside( const real point[] ) const
{
    return outer->inside(point) && ! inner->inside(point);
}



//------------------------------------------------------------------------------
void SpaceCombine::project(const real point[], real proj[]) const
{
    real proj2[DIM];
    
    inner->project( point, proj );
    outer->project( point, proj2 );
    
    // calculate the distance to the two projections:
    real pp = 0;
    for ( int dd = 0; dd < DIM; ++dd )
    {
        pp += ( point[dd] -  proj[dd] ) * ( point[dd] -  proj[dd] )
            - ( point[dd] - proj2[dd] ) * ( point[dd] - proj2[dd] );
    }
    
    // keep the shortest distance:
    if ( pp > 0 )
    {
        for ( int dd = 0; dd < DIM; ++dd )
            proj[dd] = proj2[dd];
    }
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"

bool SpaceCombine::display() const
{
    glLineStipple(1, 0x000F);
    glEnable(GL_LINE_STIPPLE);
    inner->display();
    glDisable(GL_LINE_STIPPLE);
    
    outer->display();
    return true;
}

#else

bool SpaceCombine::display() const
{
    return false;
}

#endif

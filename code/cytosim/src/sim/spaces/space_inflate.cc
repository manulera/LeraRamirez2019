// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_inflate.h"
#include "exceptions.h"


SpaceInflate::SpaceInflate(const SpaceProp* sp, const Space * spc, real len)
: Space(sp), radius(mLength[0]), radiusSqr(mLengthSqr[0])
{
    if ( DIM <= 1 )
        throw InvalidParameter("space:inflation is only valid in DIM=2 or 3");
    
    mSpace = spc;
    Space::resize(0, len);
}


SpaceInflate::~SpaceInflate()
{
    if ( mSpace )
    {
        delete( mSpace );
        mSpace = 0;
    }
}

//------------------------------------------------------------------------------

void SpaceInflate::boundaries(Vector& inf, Vector& sup) const
{
    mSpace->boundaries(inf, sup);
    inf -= Vector(radius,radius,radius);
    sup += Vector(radius,radius,radius);
}


bool SpaceInflate::inside(const real point[]) const
{
    bool in = mSpace->inside(point);
    
    if ( radius > 0  &&  in )
        return true;

    if ( radius < 0  && !in )
        return false;
    
    real proj[DIM];
    mSpace->project(point, proj);
        
    real n = 0;
    for ( int d = 0; d < DIM; ++d )
        n += (point[d]-proj[d])*(point[d]-proj[d]);
        
    if ( radius > 0 )
        return ( n <= radiusSqr );
    else
        return ( n >= radiusSqr );
}





//------------------------------------------------------------------------------
void SpaceInflate::project(const real point[], real proj[]) const
{
    ABORT_NOW("unfinished SpaceInflate");
#if ( 0 )
    switch ( mSpace->project(point, proj) )
    {
        case 0:
            return 0;
            
        case INTER_POINT:
            if ( radius ) {
                *radius = radius;
                return INTER_ARC;
            }
            //if (radius==0) we fall on 1
            
        case INTER_ARC:
            if ( radius ) {
                *radius += radius;
                return INTER_ARC;
            }
            //if (radius==0) we fall on 1
            
        case 1: {
            
            real n = 0, pw[DIM];
            for ( int d = 0; d < DIM; ++d ) {
                pw[d] = point[d] - proj[d];
                n += pw[d] * pw[d];
            }
            
            if ( n > 0 )
                n = ( mSpace->inside(point) ? -1 : +1 ) * radius / sqrt(n);
            else {
                throw Exception("insufficient space in SpaceInflate");
            }
            
            for ( int d=0; d<DIM; ++d )
                proj[d] += n * pw[d];
            
        } return 1;
    }
#endif
}

//------------------------------------------------------------------------------

void SpaceInflate::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    ABORT_NOW("unfinished SpaceInflate");
}


void SpaceInflate::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    ABORT_NOW("unfinished SpaceInflate");
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpaceInflate::display() const
{
    glLineStipple(1, 0x000F);
    glEnable(GL_LINE_STIPPLE);
    mSpace->display();
    glDisable(GL_LINE_STIPPLE);
    return false;
}

#else

bool SpaceInflate::display() const
{
    return false;
}


#endif


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_rotate.h"
#include "exceptions.h"


SpaceRotate::SpaceRotate(Space * spc)
: Space(spc->prop)
{
    if ( DIM <= 1 )
        throw InvalidParameter("space:rotate is only valid in DIM=2 or 3");
    if ( !spc )
        throw InvalidParameter("space:rotate invoked with void argument");

    mSpace = spc;
}


SpaceRotate::~SpaceRotate()
{
    if ( mSpace )
    {
        delete( mSpace );
        mSpace = 0;
    }
}


void SpaceRotate::forward(const real src[], real dest[]) const
{
    dest[0] = -src[2];
    dest[1] =  src[1];
    dest[2] =  src[0];
}


void SpaceRotate::backward(const real src[], real dest[]) const
{
    dest[0] =  src[2];
    dest[1] =  src[1];
    dest[2] = -src[0];
}


void SpaceRotate::boundaries(Vector& inf, Vector& sup) const
{
    mSpace->boundaries(inf, sup);
    Vector tmp = inf;
    forward(tmp, inf);
    tmp = sup;
    forward(tmp, sup);
}


bool SpaceRotate::inside( const real point[] ) const
{
    real rot[DIM];
    backward(point, rot);
    return mSpace->inside(rot);
}


void SpaceRotate::project(const real point[], real proj[]) const
{
    real rot[DIM], pro[DIM];
    backward(point, rot);
    mSpace->project(rot, pro);
    forward(pro, proj);
}


/* It would be necessary to create a fake Meca, to swap the indices between X and Z */
void SpaceRotate::setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff) const
{
    ABORT_NOW("unfinished SpaceRotate");
}


/* It would be necessary to create a fake Meca, to swap the indices between X and Z */
void SpaceRotate::setInteraction(Vector const& pos, PointExact const&, real rad, Meca &, real stiff) const
{
    ABORT_NOW("unfinished SpaceRotate");
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"

bool SpaceRotate::display() const
{
    glRotated(90, 0,  1, 0);
    mSpace->display();
    glRotated(90, 0, -1, 0);
    return true;
}

#else

bool SpaceRotate::display() const
{
    return false;
}

#endif

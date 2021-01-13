// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "assert_macro.h"
#include "point_exact.h"
#include "point_interpolated.h"
#include "mecable.h"
#include "iowrapper.h"
#include "simul.h"


void PointExact::read(Inputter & in, Simul& sim)
{
    Tag tag = 0;
    mPS = Simul::toMecable(sim.readReference(in, tag));
    if ( mPS )
        mPoint = in.readUInt16();
    else
        mPoint = 0;
}


void PointExact::write(Outputter& out) const
{
    out.writeSoftSpace();
    if ( mPS ) {
        mPS->writeReference(out);
        out.writeUInt16(mPoint);
    }
    else {
        Object::writeNullReference(out);
    }
}


bool PointExact::overlapping(const PointExact & p) const
{
    return ( mPS == p.mPS  &&  mPoint == p.mPoint );
}


bool PointExact::near(const PointExact & p) const
{
    return ( mPS == p.mPS  &&
            ( mPoint == p.mPoint || mPoint == p.mPoint+1 || mPoint+1 == p.mPoint ));
}


std::ostream& operator << (std::ostream& os, PointExact const& p)
{
    if ( p.mecable() )
        os << "(" << p.mecable()->reference() << " pt " << p.point() << ")";
    else
        os << "(null)";

    return os;
}

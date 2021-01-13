// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "point_interpolated.h"
#include "point_exact.h"
#include "point_set.h"
#include "fiber_locus.h"


PointInterpolated::PointInterpolated(FiberLocus const& loc, real abs)
{
    mPS     = loc.fiber();
    mPoint1 = loc.point();
    mPoint2 = loc.point()+1;
    mCoef   = abs / loc.len();
}


//------------------------------------------------------------------------------


bool PointInterpolated::overlapping(const PointExact & p) const
{
    return ( mPS==p.mPS  &&
             ( mPoint1==p.mPoint  ||  mPoint2==p.mPoint ));
}


bool PointInterpolated::overlapping(const PointInterpolated & p) const
{
    return ( mPS==p.mPS  &&
             ( mPoint1==p.mPoint1  ||  mPoint1==p.mPoint2 ||
               mPoint2==p.mPoint1  ||  mPoint2==p.mPoint2 ));
}


std::ostream& operator << (std::ostream& os, PointInterpolated const& p)
{
    if ( p.mecable() )
        os << "(" << p.mecable()->reference() << " int " << p.point1() << " " << p.point2() << " " << p.coef1() << ")";
    else
        os << "(null)";
    
    return os;
}

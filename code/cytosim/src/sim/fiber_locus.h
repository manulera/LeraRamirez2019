// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FIBER_LOCUS_H
#define FIBER_LOCUS_H


#include "real.h"
#include "vector.h"
#include "fiber.h"

class PointInterpolated;

/// Represents the segment between two consecutive points of a Fiber
/** 
 FiberLocus is used to refer to the entire segment of a Fiber.

 It is used to calculate the distance to this segment,
 or the intersection of the segment with a plane.
 
 @todo The logical name for this object would be FiberSegment
 */
class FiberLocus 
{
private:
    
    /// Fiber to which the segment belongs to
    Fiber const*   mFib;
    
    /// index of segment
    unsigned int   mPoint;
    
public:
    
    /// construct without initialization
    FiberLocus() {}
    
    /// constructor
    FiberLocus(Fiber const* f, int r) : mFib(f), mPoint(r) {}

    /// the Fiber
    Fiber const* fiber()       const { return mFib; }
    
    /// index of segment
    unsigned int point()       const { return mPoint; }
    
    /// abscissa at start of segment (i.e. corresponding to point())
    real         abscissa1()   const { return mFib->abscissaPoint(mPoint); }
    
    /// abscissa of second point
    real         abscissa2()   const { return mFib->abscissaPoint(mPoint+1); }

    /// the length of the segment
    real         len()         const { return mFib->segmentation(); }
    
    /// position of first point
    Vector       pos1()        const { return mFib->posP(mPoint); }
    
    /// position of second point
    Vector       pos2()        const { return mFib->posP(mPoint+1); }

    /// interpolated position, where c is in [0, 1]
    Vector       pos(real c)   const { return mFib->interpolatePoints(mPoint, mPoint+1, c); }
    
    /// that is [ pos2() + pos1() ] / 2
    Vector       center()      const { return mFib->interpolatePoints(mPoint, mPoint+1, 0.5); }

    /// that is pos2() - pos1()
    Vector       diff()        const { return mFib->diffPoints(mPoint); }

    /// that is ( pos2() - pos1() ).normalized()
    Vector       dir()         const { return mFib->dirPoint(mPoint); }
    
    /// PointExact corresponding to first point
    PointExact   exact1()      const { return PointExact(mFib, mPoint); }
    
    /// PointExact corresponding to second point
    PointExact   exact2()      const { return PointExact(mFib, mPoint+1); }
    
    /// true if the segment is the first of the Fiber
    bool         isFirst()     const { return ( mPoint == 0 ); }

    /// true if the segment is not the first of the Fiber
    bool         notFirst()    const { return ( mPoint > 0 ); }
    
    /// true if the segment is the last of the fiber
    bool         isLast()      const { return ( mPoint+2 == mFib->nbPoints() ); }
    
    /// true if the segment is not the last of the fiber
    bool         notLast()     const { return ( mPoint+2 < mFib->nbPoints() ); }

    
    /// calculate the projection of `w` on the line supporting the segment
    real         projectPoint0(Vector const& w, real& dist) const;

    /// calculate the projection of `w` on the line supporting the segment
    real         projectPoint(Vector const& w, real& dist) const;

    /// faster projectionPoint, but incompatible with periodic boundary conditions
    real         projectPointF(const real[], real& dist) const;

    /// calculates the closest distance between two segments
    int          shortestDistance(FiberLocus const& that, real& a, real& b, real& dis) const;

};

/// print for debug purpose
std::ostream& operator << (std::ostream&, FiberLocus const&);


#endif


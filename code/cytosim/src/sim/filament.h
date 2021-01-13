// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FILAMENT_H
#define FILAMENT_H

// Option to segment Fibers according to their curvature (experimental)
//#define CURVATURE_DEPENDENT_SEGMENTATION

#include "sim.h"
#include "vector.h"
#include "common.h"
#include "point_set.h"
#include "point_interpolated.h"


class PointExact;


/// Mecable with linear geometry
/**
 This PointSet describes a thin flexible filament that is longitudinally incompressible.
 The curvilinear length of the filament can be changed by growP(), growM(), cutP() and cutM().
 
 \par Number of points:
 
 The best number of points to describe a Filament is automatically calculated:
 It is the integer `number_of_points` that minimizes:
 @code
 fabs( length() / number_of_points - FiberProp::segmentation )
 @endcode
 
 where segmentation is a parameter of the fiber class.
 All the segments in a fiber all have the same length
 @code
 Filament::segmentation() = length() / ( number_of_points - 1 )
 @endcode

 Note that Filament::segmentation() is not always equal to FiberProp::segmentation.
 If the fibers have various length, their segmentation() will be different,
 even though they all share the same value of FiberProp::segmentation.

 See related functions: length(), nbPoints() and segmentation().
 
 \par Longitudinal incompressibility:
 
 Successive model-points are kept at a constant distance via constrained dynamics:
 @code
 ( posP(N+1)-posP(N) ).norm() == Filament::segmentation()
 @endcode
 
 \par Origin:
 
 An abscissa is a curvilinear distance taken along the Fiber,
 and the Filament provides an origin to make this independent of the model-points. 
 Thus even if the fiber lengthen from its ends, a position described by an abscissa will
 stay associated with the same local lattice site.
 
 Functions are provided in Filament to convert abscissa measured from different references,
 and to obtain positions of the fiber for a given abcissa.

 \par Derived classes:
 
 The class FiberBinder keeps track of its position using an abscissa from the origin,
 and all Hand objects are built from this class.
 The class Fiber keeps track of the FiberBinder that are attached to itself.
*/
class Filament : public PointSet
{
    /// the ideal number of points for ratio = length / segmentation
    static unsigned bestNumberOfPoints(real ratio);

    /// calculate length of given string of points
    static real trueLength(const real* pts, unsigned n_pts);
    
private:
        
    /// actual section length: distance between consecutive points
    real         fnCut;
    
    /// target segmentation length (equal to parameter 'fiber:segmentation')
    real         fnSegmentation;
  
#ifdef CURVATURE_DEPENDENT_SEGMENTATION
    /** number of time steps between each attempt to remove/add a point
    also number of states over which curvature information is averaged */
    static const unsigned RECUT_PERIOD = 3;
    
    /// error due to the cutting at different steps
    real         fnCutError;
    
    /// index into the fnCutError[]
    unsigned int fnCutErrorIndex;
#endif
    
    /// abscissa of the minus-end (equal to zero initially)
    real         fnAbscissa;
    
    
    /// oldest method to restore the distance between successive model-points
    static void  reshape_sure(unsigned, real*, real cut);
    
    /// guess the movements needed to restore the distance between two points
    static void  reshape_guess(real*, unsigned, const Vector*, real, real*, real*);

    /// apply the forces movements needed to the distance between two points
    static void  reshape_apply(unsigned, const real*, real*, const real*, const Vector*);
    
    /// restore the distance between two points
    static void  reshape_two(const real*, real*, real cut);

    /// iterative method to restore the distance between successive model-points
    static int   reshape_it(unsigned, const real*, real*, real cut);

protected:
    
    /// flag to update
    bool         needUpdate;
    
    /// callback to signal that update is needed, to be called after a change in length
    virtual void postUpdate() { needUpdate = true; }
    
public:
    
    /// vector orthogonal to backbone at the origin
    mutable Vector normal;
    
    /// Constructor
    Filament();
    
    /// Destructor
    ~Filament() {}

    //---------------------

    /// set position of MINUS_END and direction (length and Nb of points are not modified)
    /** dir does not need to be normalized */
    void         setStraight(Vector const& pos, Vector const& dir);

    /// set position of 'ref' and direction of Fiber
    void         setStraight(Vector const& pos, Vector const& dir, FiberEnd ref);

    /// set position of 'ref', direction and length of Fiber
    void         setStraight(Vector const& pos, Vector const& dir, real len, FiberEnd ref);
    
    /// import shape from the given array of size DIM*n_pts, and create a shape with `np` points
    void         setShape(const real pts[], unsigned n_pts, unsigned np);
    
    /// change the current segmentation to force `length()==len` (normally not needed)
    void         imposeLength(real len) { fnCut = len / ( nbPoints() - 1 ); }

    //---------------------
    
    /// exact representation of given end
    PointExact         exactEnd(FiberEnd) const;

    /// interpolation representing MINUS_END
    PointInterpolated  interpolateEndM() const { return PointInterpolated(this, 0, 1, 0); }

    /// interpolation representing PLUS_END
    PointInterpolated  interpolateEndP() const { return PointInterpolated(this, nbPoints()-2, nbPoints()-1, 1); }

    /// interpolation representing a given end
    PointInterpolated  interpolateEnd(FiberEnd) const;

    /// interpolation representing the mid-point between the two ends
    PointInterpolated  interpolateCenter() const;

    /// interpolation of the site specified by its distance from the ORIGIN
    PointInterpolated  interpolate(real ab) const;
    
    /// interpolation of the site specified from the MINUS_END
    PointInterpolated  interpolateM(real ab) const;
    
    /// interpolation of a site specified by its distance from a FiberEnd
    PointInterpolated  interpolate(real ab, FiberEnd from) const;
    
    //---------------------
    
    /// the total length of the Fiber, estimated from the segmentation and number of segment
    real         length()                const { return nbSegments() * fnCut; }
    
    /// the sum of the distance between model points (used for debugging)
    real         trueLength()            const { return trueLength(psPos, nbPoints()); }
    
    /// true if ( abscissaM() <= a <= abscissaP() )
    bool         betweenMP(const real a) const { return abscissaM() <= a && a <= abscissaP(); }
    
    /// true if abscissa is smaller than abscissa of PLUS_END
    bool         belowP(const real a)    const { return a <= abscissaP(); }
    
    /// true if abscissa is greater than abscissa of MINUS_END
    bool         aboveM(const real a)    const { return abscissaM() <= a; }
    
    /// calculate the domain in which ab is located (near a FiberEnd, or central)
    FiberEnd     whichEndDomain(real a, real lambda) const;

    //---------------------
    
    /// signed distance from ORIGIN to MINUS_END (abscissa of MINUS_END)
    real         abscissaM()             const { return fnAbscissa; }
    
    /// abscissa of center, midway between MINUS_END and PLUS_END
    real         abscissaC()             const { return fnAbscissa + 0.5 * length(); }

    /// signed distance from ORIGIN to PLUS_END (abscissa of PLUS_END)
    real         abscissaP()             const { return fnAbscissa + length(); }
    
    /// signed distance from ORIGIN to model point specified with index (or intermediate position)
    real         abscissaPoint(const real n) const { return fnAbscissa + fnCut * n; }

    /// signed distance from the ORIGIN to the specified FiberEnd
    real         abscissaEnd(FiberEnd end)  const;
    
    /// converts abscissa from the specified FiberEnd, to abscissa from the ORIGIN
    real         abscissaFrom(real ab, FiberEnd ref) const;

    //---------------------

#if ( DIM == 1 )
    /// position of a point specified by abscissa from the MINUS_END
    Vector       posM(real ab) const
    {
        return psPos[1] > psPos[0] ? Vector(psPos[0]+ab, 0) : Vector(psPos[0]-ab, 0);
    }
    
    /// position of a point specified by abscissa from the ORIGIN
    Vector       pos(real ab) const
    {
        return psPos[1] > psPos[0] ? Vector(psPos[0]+ab-fnAbscissa, 0) : Vector(psPos[0]-ab+fnAbscissa, 0);
    }
#else
    /// position of a point specified by abscissa from the MINUS_END
    Vector       posM(real ab) const;
    
    /// position of a point specified by abscissa from the ORIGIN
    Vector       pos(real ab) const { return posM(ab-fnAbscissa); }
#endif

    /// position of a point specified by abscissa `ab` from reference `ref`
    Vector       pos(real ab, FiberEnd ref) const;

    /// position of the point taken mid-way along the curve
    Vector       posMiddle() const { return posM(0.5*length()); }
    
    /// position of a FiberEnd
    Vector       posEnd(FiberEnd which) const;
    
    /// position of MINUS_END
    Vector       posEndM() const { return Vector(psPos); }

    /// position of PLUS_END
    Vector       posEndP() const { return Vector(psPos+DIM*(nbPoints()-1)); }

    //---------------------
    
    /// vector between two consecutive points `p` and `p+1` (alias to diffPoints())
    Vector       diffP(unsigned p) const { return diffPoints(p); }

#if ( 1 )
    /// normalized tangent vector to the fiber within segment [p, p+1]
    /** We divide by fnCut, which should be the distance between points */
    Vector       dirPoint(unsigned p)  const { return diffPoints(p) / fnCut; }
#else
    /// normalized tangent vector to the fiber within segment [p, p+1]
    /** We normalize the difference between points */
    Vector       dirPoint(unsigned p)  const { return diffPoints(p).normalized(); }
#endif

        
    /// normalized tangent vector to the fiber at given abscissa from the origin
    Vector       dir(real ab) const;

    /// normalized tangent vector to the fiber at given abscissa from given reference
    Vector       dir(real ab, FiberEnd from) const;
    
    /// normalized tangent vector to the fiber at given end
    Vector       dirEnd(FiberEnd which) const;
    
    /// normalized tangent vector to the fiber at MINUS_END
    Vector       dirEndM() const { return dirPoint(0); }
    
    /// normalized tangent vector to the fiber at PLUS_END
    Vector       dirEndP() const { return dirPoint(lastSegment()); }

    /// dot-product (force at the end of the Fiber).(direction of Fiber growth)
    real         projectedForceEnd(FiberEnd which) const;

    /// force on the PLUS_END projected on the direction of elongation
    real         projectedForceEndM() const { return -netForce(0) * dirPoint(0); }

    /// force on the PLUS_END projected on the direction of elongation
    real         projectedForceEndP() const { return netForce(lastPoint()) * dirPoint(lastSegment()); }
    
    /// angle of segment `p` in the XY-plane
    real         angleXY(unsigned p) const;
    
    //--------------------- Segmentation / discrete representation
    
    /// set desired segmentation (the length of the segments might be different)
    void         segmentation(real c) { assert_true(c>0); fnSegmentation = c; }
    
    /// the current segment length (distance between successive model-points)
    real         segmentation() const { return fnCut; }
    
    /// returns segmentation() ^ 3
    real         segmentationCub() const { return fnCut*fnCut*fnCut; }
    
    /// recalculate fiber to have `np` model points
    void         resegment(unsigned np);
    
    /// automatically select the number of points if needed, and resegment the fiber
    void         adjustSegmentation();
    
    /// restore the distance between successive model-points
    void         reshape();
    
    /// change position
    void         getPoints(const real * x);
    
    /// invert polarity
    void         flip();
    
    //--------------------- Info
    
    /// calculate the minimum and maximum segment length
    void         minMaxSegments(real&, real&) const;

    /// calculate average and variance of the segment length
    void         infoSegments(real&, real&) const;

    /// curvature calculated at joint `p`, where `0 < p < nbPoints()-1`
    real         curvature(unsigned p) const;
    
    /// normalized energy associated with bending
    real         bendingEnergy0() const;

    /// the cosine of the maximum segment angle: indicate the errors due to curvature
    real         minCosinus() const;
    
    /// number of joints at which ( cosine(angle) < threshold )
    unsigned     nbKinks(real threshold = 0) const;
    
    /// calculate intersection between segment `s` and the plane defined by <em> n.pos + a = 0 </em>
    real         planarIntersect(unsigned s, Vector const& n, const real a) const;

    //--------------------- Growing/Shrinking
    
    /// merge two fibers by attaching `fib` at the PLUS_END of `this`
    void         join(Filament const* fib);

    /// increase/decrease length of Fiber by `dlen`, at the MINUS_END
    void         growM(real dlen);
    
    /// add a segment of length segmentation() at the MINUS_END
    void         addSegmentM();
    
    /// remove a portion of length `dlen` including the MINUS_END
    void         cutM(real dlen);
    
    /// increase/decrease length of Fiber by `dlen`, at the PLUS_END
    void         growP(real dlen);
    
    /// add a segment of length segmentation() at the PLUS_END
    void         addSegmentP();
    
    /// remove a portion of length `dlen` including the PLUS_END
    void         cutP(real dlen);
    
    /// grow at specified end (PLUS_END or MINUS_END)
    void         grow(FiberEnd end, real dlen);
    
    /// shorten or lengthen Fiber without changing the position of `ref`
    void         adjustLength(real len, FiberEnd ref);

    /// Discard model points in [ 0, P-1 ] and keep [ P, end ]
    void         truncateM(unsigned p);

    /// Keep model points [ 0, P ] and discard the others
    void         truncateP(unsigned p);

    //---------------------
    
    /// check the length of the segments
    void         checkLength() const;
    
    /// dump for debugging
    void         dump(std::ostream&) const;
    
    /// write to Outputter
    void         write(Outputter&) const;
    
    /// read from Inputter
    void         read(Inputter&, Simul&, Tag);
    
    
    

};


#endif

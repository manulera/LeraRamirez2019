// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef POINT_SET_H
#define POINT_SET_H

#include "assert_macro.h"
#include "dim.h"
#include "movable.h"
#include "object.h"
#include "mecable.h"
#include "isometry.h"

class PointInterpolated;
class Space;
class Modulo;
class Simul;

/// Array of points describing a physical object
/**
 This implements the interface defined by Mecable,
 and defines an Object with a variable number of points.
 */
class PointSet : public Mecable
{
    
private:
    
    /// Allocation size of arrays psPos[] and psFor[]
    unsigned int  psAllocated;
    
    /// Constructor base
    void          psConstructor();

protected:
        
    /// psPos[] of size DIM*psAllocated contains DIM*pSize point-coordinates
    real*         psPos;
    
    /// psFor[] contains force-coordinates and is allocated in Meca
    const real*   psFor;
    
    //--------------------------------------------------------------------------
public:
        
    /// Constructor
    PointSet();
    
    /// Copy constructor
    PointSet(const PointSet &);
    
    /// Assignement operator
    PointSet & operator =(const PointSet &);
    
    /// Destructor
    virtual ~PointSet()  { deallocatePoints(); }
    
    
    /// allocate memory to store 'nbp' points
    virtual unsigned allocatePoints(unsigned nbp);
    
    /// free memory allocated by allocatePoints()
    virtual void  deallocatePoints();
    
    /// Set the number of points in the array
    void          setNbPoints(const unsigned n)  { allocatePoints(n); pSize = n; }

    //--------------------------------------------------------------------------

    /// Position of point 'p' of the object (this is the non-virtual equivalent of posPoint())
    Vector        posP(const unsigned p) const { assert_true(psPos && p<pSize); return Vector(psPos+DIM*p); }
    
    ///  Position of point 'p' of the object
    Vector        posPoint(unsigned p)   const { assert_true(psPos && p<pSize); return Vector(psPos+DIM*p); }

    /// Address of point `p`
    const real*   data()                 const { return psPos; }
    
    /// Address of point `p`
    const real*   addrPoint(const unsigned p) const { return psPos + DIM*p; }

    /// Add a point, returning the array index that was used
    unsigned int  addPoint(Vector const& w);
    
    /// Remove `nbp` points starting from index `inx`
    void          removePoints(unsigned inx, unsigned nbp);
    
    /// Shift `nbp` points starting from index `inx`
    void          shiftPoints(unsigned inx, unsigned nbp);
    
    /// Remove all points with indices [ 0, p-1 ], keep [ p, nbPoints() ]
    virtual void  truncateM(unsigned int p);
    
    /// Keep points [ 0, p ], remove other points
    virtual void  truncateP(unsigned int p);
    
    /// Remove all points
    void          clearPoints()  { pSize = 0; }
    
    /// Set all coordinates to zero (nicer for debug/testing)
    void          resetPoints();
    
    /// Add random noise uniformly to all coordinate (used for testing purposes)
    void          addNoise(real amount);
    
    /// copy current coordinates to argument
    virtual void  putPoints(real *) const;
    
    /// replace current coordinates by provided ones
    virtual void  getPoints(const real *);
    
    //--------------------------------------------------------------------------
    
    /// return force-vector on point `p` calculated at previous step by Meca
    Vector        netForce(const unsigned p) const;
    
    /// replace current forces by the ones provided
    virtual void  getForces(const real * ptr) { psFor = ptr; }
    
    //--------------------------------------------------------------------------
    //These functions are defined here to enable inlining, which may be faster
    
    /// Set position of point p to w
    void   setPoint(unsigned P, Vector const& w)
    {
        assert_true( P < pSize );
        w.put(psPos+DIM*P);
    }
    
    /// Shift point at index p by vector w
    void   movePoint(const unsigned P, Vector const& w)
    {
        assert_true( P < pSize );
        w.add_to(psPos+DIM*P);
    }
    
    
    /// Difference of two points = (P+1) - P
    static inline Vector diffPoints(const real* src, const unsigned P)
    {
#if ( DIM == 1 )
        return Vector(src[P+1]-src[P], 0);
#elif ( DIM == 2 )
        return Vector(src[2*P+2]-src[2*P], src[2*P+3]-src[2*P+1]);
#elif ( DIM == 3 )
        return Vector(src[3*P+3]-src[3*P], src[3*P+4]-src[3*P+1], src[3*P+5]-src[3*P+2]);
#else
        return Vector(src[DIM*P+DIM]-src[DIM*P], src[DIM*P+DIM+1]-src[DIM*P+1], src[DIM*P+DIM+2]-src[DIM*P+2]);
#endif
    }

    /// Difference of two points = Q - P
    static inline Vector diffPoints(const real* src, const unsigned P, const unsigned Q)
    {
#if ( DIM == 1 )
        return Vector(src[P]-src[Q], 0);
#elif ( DIM == 2 )
        return Vector(src[2*Q]-src[2*P], src[2*Q+1]-src[2*P+1]);
#elif ( DIM == 3 )
        return Vector(src[3*Q]-src[3*P], src[3*Q+1]-src[3*P+1], src[3*Q+2]-src[3*P+2]);
#else
        return Vector(src[DIM*Q]-src[DIM*P], src[DIM*Q+1]-src[DIM*P+1], src[DIM*Q+2]-src[DIM*P+2]);
#endif
    }

    /// Difference of two consecutive points: (P+1) - (P)
    Vector diffPoints(const unsigned P) const
    {
        assert_true( P+1 < pSize );
        return diffPoints(psPos, P);
    }

    /// Difference of two points = Q - P = vector PQ
    Vector diffPoints(const unsigned P, const unsigned Q) const
    {
        assert_true( P < pSize );
        assert_true( Q < pSize );
        return diffPoints(psPos, P, Q);
    }

    /// Calculate intermediate position = P + a ( Q - P )
    Vector interpolatePoints(const unsigned P, const unsigned Q, const real a) const
    {
        assert_true( P < pSize );
        assert_true( Q < pSize );
#if ( DIM == 1 )
        return Vector(psPos[P]+a*(psPos[Q]-psPos[P]), 0);
#elif ( DIM == 2 )
        return Vector(psPos[2*P]+a*(psPos[2*Q]-psPos[2*P]), psPos[2*P+1]+a*(psPos[2*Q+1]-psPos[2*P+1]));
#elif ( DIM == 3 )
        return Vector(psPos[3*P]+a*(psPos[3*Q]-psPos[3*P]), psPos[3*P+1]+a*(psPos[3*Q+1]-psPos[3*P+1]), psPos[3*P+2]+a*(psPos[3*Q+2]-psPos[3*P+2]));
#else
        const real * pp = psPos + DIM*P;
        const real * qq = psPos + DIM*Q;
        return Vector(pp[0]+a*(qq[0]-pp[0]), p[1]+a*(qq[1]-pp[1]), p[2]+a*(qq[2]-pp[2]));
#endif
    }
    
    /// calculate first and second momentum of point coordinates
    void          calculateMomentum(Vector&, Vector&, bool sub);
    
    //--------------------------------------------------------------------------
    //                      Position-related functions
    //--------------------------------------------------------------------------
    
    /// Position of center of gravity
    virtual Vector position() const;
    
    /// true if object accepts translations
    virtual bool  mobile() const { return true; }

    /// Translate object (moves all the points by w)
    virtual void  translate(Vector const&);
    
    /// Rotate object by given rotation
    virtual void  rotate(Rotation const&);
    
    /// Modulo around the first point
    virtual void  foldPosition(Modulo const*);
    
    /// true is all points are inside Space
    bool          allInside(Space const*) const;
    
    //--------------------------------------------------------------------------
    
    /// Write to file
    void          write(Outputter&) const;

    /// write as text
    void          write(std::ostream&) const;

    /// Read from file
    void          read(Inputter&, Simul&, Tag);
    
};



/// output operator:
std::ostream& operator << (std::ostream& os, PointSet const&);


#endif

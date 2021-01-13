// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FIBER_BINDER_H
#define FIBER_BINDER_H

#include "assert_macro.h"
#include "point_interpolated.h"
#include "fiber.h"
#include "sim.h"


/// a location on a Fiber represented by its abscissa from the Fiber's origin
/**
 A FiberBinder has a pointer to a Fiber, which is:
 - zero if the Binder is not attached,
 - the corresponding Fiber when the Binder is attached.
 .
 
 When the FiberBinder is attached, its location on the Fiber
 is stored as a curvilinear abscissa (fbAbs) taken along the fiber, 
 from a fixed reference on the fiber.
 The abscissa is independent from the Fiber's model points.
 
 //@todo FiberBinder should include the functionalities of Digit
*/
class FiberBinder: public Node
{
private:
    
    /// the corresponding interpolation, which is kept up-to-date
    /**
     normally, ( inter.mecable() == fbFiber ) and ( abscissaInter() == fbAbs )
     */
    PointInterpolated inter;
    
    /// the abscissa of the interpolated point, which should be equal to `abscissa()`
    real abscissaInter() const { return fbFiber->abscissaPoint(inter.point1()+inter.coef1()); }

protected:
    
    /// the Fiber on which it is attached, or 0 if not attached
    Fiber*     fbFiber;
    
    /// the abscissa from the origin of the Fiber
    real       fbAbs;
    
public:
    
    /// construct as unattached
    FiberBinder() : fbFiber(0), fbAbs(0) {}
    
    /// construct at the given distance from the origin
    FiberBinder(Fiber* f, real a);

    //--------------------------------------------------------------------------
    
    /// set as bound at position `a` on Fiber `f`
    void         locate(Fiber* f, real a);

    /// set as bound at position represented by FiberBinder
    void         locate(FiberBinder const& a) { locate(a.fbFiber, a.fbAbs); }
    
    /// set as unbound
    void         delocate();

    /// move to a different fiber, at same abscissa
    virtual void relocate(Fiber* f);
    
    /// move to a different abscissa on the same fiber
    virtual void relocate(real a);

    /// move to a different fiber, at given position
    virtual void relocate(Fiber* f, real a);
    
    /// relocate to MINUS_END on the current fiber
    virtual void moveToEndM();
    
    /// relocate to PLUS_END on the current fiber
    virtual void moveToEndP();
    
    /// relocate to the specified tip of the current fiber
    virtual void moveToEnd(FiberEnd);

    //--------------------------------------------------------------------------
    
    /// true if not attached
    bool         unattached()    const { return fbFiber == 0; }

    /// true if attached
    bool         attached()      const { return fbFiber != 0; }
    
    /// Fiber to which this is attached, or zero if not attached
    Fiber*       fiber()         const { return fbFiber; }
    
    /// position in space (using current interpolation)
    Vector       pos()           const { assert_false(bad()); return inter.pos(); }
    
    /// position (always works)
    Vector       posHand()       const { return fbFiber->pos(fbAbs); }
    
    /// direction of Fiber obtained by normalization
    Vector       dir()           const { return inter.dir(); }
    
    /// the direction of the Fiber at the point of attachment
    Vector       dirFiber()      const { return fbFiber->dirPoint(inter.point1()); }
    
    /// the abscissa, from the origin of the Fiber
    real         abscissa()      const { return fbAbs; }

    /// abscissa, counted from the MINUS_END
    real         abscissaFromM() const { return fbAbs - fbFiber->abscissaM(); }

    /// abscissa, counted from the PLUS_END in the reverse direction
    real         abscissaFromP() const { return fbFiber->abscissaP() - fbAbs; }

    /// abscissa, counted from the specified FiberEnd (in reversed direction for the PLUS_END)
    real         abscissaFrom(FiberEnd from) const;
            
    /// nearest end to the current attachment point
    FiberEnd     nearestEnd() const;
    
    /// true if abscissa is below abscissaP
    bool         belowP()        const { return fbFiber->belowP(fbAbs); }
    
    /// true if abscissa is above abscissaM
    bool         aboveM()        const { return fbFiber->aboveM(fbAbs); }
    
    /// true if abscissa is within the fiber boundaries
    bool         betweenMP()     const { return fbFiber->betweenMP(fbAbs); }
    
    //--------------------------------------------------------------------------
    
    /// the interpolation
    const PointInterpolated& interpolation() const { assert_false(bad()); return inter; }
    
    /// set a valid PointInterpolated
    void         updateBinder()                    { inter = fbFiber->interpolate(fbAbs); }
    
    //--------------------------------------------------------------------------
    
    /// a static_cast<> of Node::next()
    FiberBinder *  next()  const  { return static_cast<FiberBinder*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    FiberBinder *  prev()  const  { return static_cast<FiberBinder*>(nPrev); }
    
    //--------------------------------------------------------------------------
    
    /// file output
    virtual void read(Inputter&, Simul&);
    
    /// file input
    virtual void write(Outputter&) const;
 
    //--------------------------------------------------------------------------
    
    /// check that fbAbs is within Fiber::abscissaM() and Fiber::abscissaP()
    void         checkAbscissa() const;
    
    /// check validity of the interpolation (debuging purposes)
    int          bad() const;
    
    /// Added by Manu: useful to access lattice values
    virtual int          lattice_site() const {return 0;};
    
    virtual int          lattice_val_add() const {return 0;};
};

/// output operator for debugging purpose
std::ostream& operator << (std::ostream&, FiberBinder const&);


#endif


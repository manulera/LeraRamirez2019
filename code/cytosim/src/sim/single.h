// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SINGLE_H
#define SINGLE_H

#include "dim.h"
#include "vector.h"
#include "movable.h"
#include "object.h"
#include "hand_monitor.h"
#include "point_exact.h"
#include "single_prop.h"
#include "hand.h"

class Meca;
class Modulo;
class Glossary;
class FiberGrid;
class Fiber;
class PointDisp;


/// A point-like object containing one Hand.
/**
 A Single contains one pointer to Hand, and consequently
 inherit the 2 possible states: `attached` or `free`.
 
 By default:
 - Free Single are diffusing, and try to bind to nearby Fibers,
 - Attached Singles are moving along the Fiber to which their Hand is attached.
 .
 
 However, two derived classes change this behavior:
 -# a Picket is fixed in position and do not diffuse,
 -# a Wrist is attached to one model point of a Mecable.
 .
 
 Attached Wrist and Picket exert a force on the Fiber to which the Hand is attached.
 For WristLong and PicketLong, this force can have a non-zero resting length.
 For these class in which the Hand can be under tension, `hasForce()` returns true.

 Wrist and Picket can be distinguished with Single::base():
 - for Single and Picket, this returns zero,
 - for Wrist, this returns the Mecable on which the Wrist is attached.
 .

 @ingroup SingleGroup
 */

class Single : public HandMonitor, public Object
{
private:
    
    /// specialization of HandMonitor
    void      afterAttachment(Hand const*);
    /// specialization of HandMonitor
    void      beforeDetachment(Hand const*);
    /// specialization of HandMonitor
    Vector    otherPosition(Hand const*) const { return posFoot(); }
    /// = identity() of the Object on which a Wrist is attached, or Single::identity()
    ObjectID  nucleatorID() const { if ( base() ) return base()->identity(); return Object::identity(); }
    /// specialization of HandMonitor
    real      interactionLength() const;
    
protected:
    
    /// the position of the foot
    Vector        sPos;

    /// the motor domain
    Hand *        sHand;
    
#ifdef TRAP_SINGLES
    /// the position of the trapped hand
    Vector        trapPos;
#endif
public:
    
    /// property
    SingleProp const* prop;

    /// constructor at specified position
    Single(SingleProp const*, Vector const& = Vector(0,0,0));

    /// destructor
    virtual ~Single();
    
    //--------------------------------------------------------------------------
    
    /// a reference to the Hand
    Hand const*  hand()                          const  { return sHand; }
    
    /// sHand->attached()
    bool    attached()                           const  { return sHand->attached(); }
    
    /// Fiber to which this is attached
    Fiber*  fiber()                              const  { return sHand->fiber(); }
    
    /// attachment position of Hand along fiber (call is invalid if Hand is not attached)
    real    abscissa()                           const  { return sHand->abscissa(); }
    
    /// position of the Hand (call is invalid if Hand is not attached)
    Vector  posHand()                            const  { return sHand->pos(); }
    
    /// direction of Fiber at attachment point (call is invalid if Hand is not attached)
    Vector  dirFiber()                           const  { return sHand->dirFiber(); }
    
    /// attach Hand at the given site
    void    attach(FiberBinder const& fb)               { sHand->attach(fb); }
    
    /// attach Hand at the given abscissa
    void    attachTo(Fiber * f, real ab)                { sHand->attachTo(f, ab); }

    /// attach Hand at specified position
    void    attachTo(Fiber * f, real ab, FiberEnd from) { sHand->attachTo(f, ab, from); }
    
    /// attach Hand at given Fiber end
    void    attachToEnd(Fiber * f, FiberEnd end)        { sHand->attachToEnd(f, end); }

    /// move Hand at given end
    void    moveToEnd(FiberEnd end)                     { sHand->moveToEnd(end); }
    
    /// detach
    void    detach()                                    { sHand->detach(); }

    //--------------------------------------------------------------------------
    
    ///return the position in space of the object
    virtual Vector  position() const;
    
    /// Single can be translated only if it is not attached
    virtual bool    mobile()                     const  { return !sHand->attached(); }
    
    /// translate object's position by the given vector
    virtual void    translate(Vector const& w)          { sPos += w; }
    
    /// move object to specified position
    virtual void    setPosition(Vector const& w)        { sPos = w; }

    /// modulo the position of the grafted
    virtual void    foldPosition(const Modulo * s);
    
    /// set the position randomly inside prop->confine_space
    void            randomizePosition();

    //--------------------------------------------------------------------------
    
    /// the position of the anchoring point
    virtual Vector  posFoot()                    const  { return sPos; }
    
    /// position on the side of fiber used for sideInteractions
    virtual Vector  posSide()                    const  { return sHand->pos(); }
    
    /// the Mecable to which this is anchored, or zero
    virtual Mecable const* base()                const  { return 0; }
    
    /// true if Single creates an interaction
    virtual bool    hasForce() const                    { return false; }

    /// force = stiffness * ( position_anchor - position_hand ), or zero for a diffusible Single
    virtual Vector  force()                      const  { return Vector(0,0,0); }

    /// Monte-Carlo step if the Hand is not attached
    virtual void    stepF(const FiberGrid&);
    
    /// Monte-Carlo step if the Hand is attached
    virtual void    stepA();
    
    /// add interactions to the Meca
    virtual void    setInteractions(Meca &) const;
    
    //--------------------------------------------------------------------------
    
    /// a static_cast<> of Node::next()
    Single*         next()   const  { return static_cast<Single*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Single*         prev()   const  { return static_cast<Single*>(nPrev); }

    //--------------------------------------------------------------------------

    /// a unique character identifying the class
    static const Tag TAG = 's';
    
    /// return unique character identifying the class
    virtual Tag     tag() const { return TAG; }
    
    /// return Object Property
    Property const* property() const { return prop; }
    
    /// read from file
    virtual void    read(Inputter&, Simul&, Tag);
    
    /// write to file
    virtual void    write(Outputter&) const;
    
    /// return PointDisp of Hand
    PointDisp *     disp() const { return sHand->prop->disp; }

#ifdef TRAP_SINGLES

    /// return whether it is trapped or not
    bool trapped() const {return sHand->trapped();};
    
#if (TRAP_SINGLES==1)
    /// Virtual HandMonitor functions
    Hand * trappedHand() const;
    
    /// release the single
    void   untrap( bool first_call);
    
    Vector trappedHandPos() const;
    void   afterTrapping();

    Hand * trappedOtherHand() const;
    Vector trappedOtherHandPos() const;
#endif
#if (TRAP_SINGLES==2)
    /// Release the single
    void   untrap();
    
    /// Trap and include in the right list
    void   trap(HandMonitor * h);
#endif

    void   stepTrappedF(const FiberGrid& grid, Vector const & pos);
    void   stepTrappedA(Vector force);
    void   stepTrappedF_AA();
    
    /// Write the reference (to be able to call it from HandMOnitor)
    void writeReference(Outputter& out) const;
    
    ///
    void getHands(Hand * &h1, Hand * &h2) const {h1 = sHand; h2 =0;};
#endif
    
};


#endif

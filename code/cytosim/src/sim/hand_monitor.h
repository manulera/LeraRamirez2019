// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef HAND_MONITOR
#define HAND_MONITOR

#include "real.h"
#include "vector.h"
#include "inventoried.h"
#include "iowrapper.h"
#include <deque>


class Hand;
class Simul;
class FiberBinder;
class FiberGrid;
class Fiber;
#define TRAP_SINGLES 2

/// base class to monitor and control Hand's actions
/**
 The HandMonitor is a mechanism for a Hand to access data from the container
 class (Single or Couple) to which it belongs. This is the only way for a Hand 
 to know features at this higher level.
 */
class HandMonitor
{

public:
    
    virtual ~HandMonitor() {}
    
    /// Returning `false` prevents the attachment (this is called before every attempt)
    virtual bool allowAttachment(const FiberBinder & site) { return true; }
    
    /// called after attachement
    virtual void afterAttachment(Hand const*) {}
    
    /// called before detachment
    virtual void beforeDetachment(Hand const*) {}
    
    /// identity() for associated object
    virtual ObjectID nucleatorID() const { return 0; }
    
    /// return the Hand that is not the argument, in a Couple
    virtual Hand * otherHand(Hand const*) const { return 0; }

    /// return the position of the Hand that is not the argument, for a Couple
    /** If the hand is not part of a Couple, this returns Vector(0,0,0) */
    virtual Vector otherPosition(Hand const*) const { return Vector(0,0,0); }

    /// return the direction of the Fiber for the Hand that is not the argument, in a Couple
    /** If the hand is not part of a Couple, this returns a random unit vector */
    virtual Vector otherDirection(Hand const*) const { return Vector::randU(); }

    /// resting length of the interaction
    virtual real   interactionLength() const { return 0; }
    
#ifdef TRAP_SINGLES
    
    virtual bool   trapped() const {return false;};
    
    virtual void   stepTrappedA(Vector force){};
    
    virtual void   stepTrappedF(const FiberGrid& grid, Vector const & pos){};
    
    virtual void   stepTrappedF_AA(){};
    
#if (TRAP_SINGLES==1)
    
    virtual void   afterTrapping(){};
    /// Allow hand_monitors to bind to each other
    virtual void   trap( Hand * h){};
    /// return the position of the hand engaged in the trap in THIS HandMonitor
    virtual Vector trappedHandPos() const{ return 0; };

    virtual void   untrap(bool first_call){};
    /// return the hand engaged in the trap in THE PARTNER HandMonitor
    virtual Hand * trappedOtherHand() const{ return 0; };
    
    /// return the position of the hand engaged in the trap in THE PARTNER HandMonitor
    virtual Vector trappedOtherHandPos() const{ return 0; };

    /// return the hand engaged in the trap in THIS HandMonitor.
    virtual Hand * trappedHand() const { return 0; }
#elif (TRAP_SINGLES==2)
    /// Allow hand_monitors to bind to each other
    virtual void   trap( HandMonitor * h){};
    
    virtual void   untrap(){};
    
    virtual Vector3 trapCenter(real exp_shift1,real exp_shift2,real exp_shift3) const {return Vector3(0,0,0);}
    virtual real trap_rate() const {return 0;};

    virtual real untrap_rate() const {return 0;};
    
    virtual int partnerLattice(Hand * h, Fiber * f) const {return 0;};
    virtual int checkPartnerLattice() const {return 0;};
#endif

    
    /// Write the reference (used by trapping haMons)
    virtual void writeReference(Outputter& out) const {};
    
    HandMonitor * trapped_haMon;
    
    /// Get the hands
    virtual void getHands(Hand * &h1, Hand * &h2) const {};
#endif
};


#endif

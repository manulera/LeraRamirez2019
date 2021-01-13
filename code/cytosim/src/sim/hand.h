// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef HAND_H
#define HAND_H

#include "fiber_binder.h"

class HandMonitor;
class FiberGrid;
class FiberProp;
class HandProp;
class Simul;
/// Simulates the stochastic binding/unbinding of a FiberBinder
/**
 A Hand is always part of a larger construct, for example Single or Couple.
 
 Hand is the parent to many class that implement different fiber-related activities.
 Hand provides binding/unbinding capacity to these derived classes.
 
 Attachment occurs with constant rate @ref HandPar "attach_rate" to any fiber located
 at distance  @ref HandPar "attach_range" or less.
 If attachment occurs, it happens on the closest point of the fiber,
 which is either the projection of the current position on the fiber axis, 
 or one of the fiber end.
 
 You can restrict binding to happen only at the ends by setting `bind_only_end`,
 and the associated cutoff distance `bind_end_range`.

 Detachment increases exponentially with force:
 @code
 off_rate = unbinding_rate * exp( force.norm() / unbinding_force )
 @endcode

 See @ref HandPar
 @ingroup HandGroup
 */
class Hand : public FiberBinder
{
    friend class Couple;
    friend class Flipper;
#ifdef TRAP_SINGLES
    friend class Trapper;
#endif
    
private:
    
    /// disabled default constructor
    Hand();
    
protected:

    /// the monitor associated with this Hand
    HandMonitor *  haMonitor;
    /// Gillespie normalized time for attachment (must be set at detachment)
    real           nextAttach;
    
    /// Gillespie normalized time for detachment (must be set at attachment)
    real           nextDetach;
    
    /// reset Gillespie's counters
    void           resetTimers();

    /// test for detachment with rate prop->unbinding_rate
    bool           testDetachment();
    
    /// test for detachment with Kramers theory
    bool           testKramersDetachment(real force);
    
public:
    
    /// Property is constant, so we do not need to make it private
    HandProp const* prop;

    /// Property
    HandProp const* property() const { return prop; }
    
    /// constructor
    /**
     To create a new Hand, you need to have a HandProp.
     HandProp is parent to several classes, exactly mirroring the hierarchy of Hands.
     
     The correct derived class can be created by its associated HandProp:
     @code
     HandProp * hp = HandProp::newProperty(name, opt);
     Hand * h = hp->newHand(this);
     @endcode
     */
    Hand(HandProp const*, HandMonitor* h);

    /// destructor
    virtual ~Hand();

    /// return other Hand if part of a Couple, and zero otherwise
    Hand *         otherHand() const;
    
    /// tell if attachment at given site is permitted
    virtual int   attachmentAllowed(FiberBinder& site) const;
    
    /// Do something in case the hand can change to allow attachemnt (the attachmentAllowed is a constant)
    virtual bool   attachmentSecondTry(FiberBinder& site) {return false;};
    
    /// attach the hand at the position described by site
    virtual void   attach(FiberBinder const& site);
    
    /// detach
    virtual void   detach();
    
    /// move along the Fiber to specified abscissa, or detach
    virtual void   moveTo(real abs);
    
    /// move along the Fiber by the abscissa offset dabs, or detach
    void           moveBy(real dabs) { moveTo(fbAbs+dabs); }

    /// simulate when the Hand is not attached
    virtual void   stepUnattached(const FiberGrid&, Vector const & pos);

    /// simulate when the Hand is attached but not under load
    virtual void   stepUnloaded();

    /// simulate when the Hand is attached and under load
    virtual void   stepLoaded(Vector const & force);
    
    /// check abscissa against fiber edge, and calls handle functions if necessary.
    void           checkFiberRange();

    /// this is called when disassembly occured PLUS_END
    virtual void   handleDisassemblyM();
    
    /// this is called when the attachment point is below the MINUS_END
    virtual void   handleDisassemblyP();

    /// attach at abscissa of given Fiber (calls attach(FiberBinder))
    void           attachTo(Fiber *, real ab);
    
    /// attach at specified distance `ab` from FiberEnd (calls attach(FiberBinder))
    void           attachTo(Fiber *, real ab, FiberEnd from);
    
    /// attach at the given end of Fiber (calls attach(FiberBinder))
    void           attachToEnd(Fiber *, FiberEnd end);
    
    /// file input
    void           read(Inputter&, Simul&);
    
    /// file output
    void           write(Outputter&) const;
    
    
    
    /// Added by Manu
    virtual real  propensUnloaded(){return 0;};
    
    virtual real  propensLoaded(Vector const& force){return 0;};
    
    virtual int  gillestepUnloaded(){return 0;};
    
    virtual int  gillestepLoaded(Vector const &){return 0;};
    
    virtual bool needs_sweep(){return false;};
    
    virtual void sweep(){};
    
    virtual bool first_sweep(){return 0;};

    virtual bool find_lattice_neighbour(unsigned int const v) const {return 0;} ;
#ifdef MULTI_LATTICE
    virtual void set_lat_id(unsigned int const){};
    virtual unsigned int get_fiberlat_id() const{return 0;};
    virtual void random_multi_lattice(){return;};
#endif
    
#ifdef TRAP_SINGLES

    bool trapped();
    
    void stepTrapped(Vector force);
    
    virtual void stepUnattachedTrappedAA();
#if (TRAP_SINGLES==1)
    
    // Set trapped_haMon to h->haMonitor
    virtual void trap(Hand * h);
    // Set trapped_haMon to zero
    virtual void untrap();
    
    Vector trappedOtherHandPos();
    Hand * trapped_hand() const;
#endif
#if (TRAP_SINGLES==2)
    
    HandMonitor * trappedHaMon();
    
#endif
#endif
};

#endif


#ifndef TRAPPER_H
#define TRAPPER_H

#include "couple.h"

class TrapperProp;

class Trapper : public Couple
{
public:
    /// property
    TrapperProp const* prop;
    
    /// create following the specifications in the CoupleProp
    Trapper(TrapperProp const*, Vector const & w = Vector(0,0,0));
    
    /// destructor TODO: THIS HAS TO BE DEFINED PROPERLY
    virtual ~Trapper();
    
    /// Only to add the flipping step (only when double bound)
    void stepAA();
    
    /// Monitor the detachment
    void afterAttachment(Hand const* h);
    
    void beforeDetachment(Hand const* h);
    
    // We have to store and read the property of the hand each time point
    void write(Outputter& out) const;
    
    void read(Inputter & in, Simul& sim, Tag tag);
    
    /// The usual gillespie counters
    real next_trap;
    
    real next_untrap;
    
    /// Returns the force vector from the hand in the couple to the trapped single
    Vector forceTrap();
    
    // The previous one is a general function, these ones are knowing that the trapped hand is either the 1 or 2 (fewer fucntion calls involved)
    Vector forceTrap1();
    Vector forceTrap2();
    
    void step_untrapped();
    
    void step_trapped(Vector force);
    
    /// See if in the lattice site, with a different lat_val we find a "target" for trapping. If it is found, verify whether it is not trapped, and then trap it.
    bool try2Trap(Hand * cHand);
    
    /// Right now there seems no need to override this function, but maybe eventually?
    void setInteractions(Meca &) const;
    
    void untrap(bool first_call);
    
    bool trapped() const {return (cHand1->trapped()||cHand2->trapped());}
    
    /// Virtual HandMonitor functions
    Hand * trappedHand() const;
    Vector trappedHandPos() const;
    Hand * trappedOtherHand() const;
    Vector trappedOtherHandPos() const;
    
    /// Write the reference (to be able to call it from HandMOnitor)
    void writeReference(Outputter& out) const;
    
};
#endif

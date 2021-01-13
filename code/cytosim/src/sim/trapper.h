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
    
    void stepFA(const FiberGrid& grid);
    
    void stepAF(const FiberGrid& grid);
    
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
    
    Vector force() const;
    Vector forceTrap1() const;
    Vector forceTrap2() const;
    
    /// The stiffness between the elements changes depending on whether all three elements are bound or only two.
    real get_stiffness() const;
    
    /// The trapping/untrapping is controlled from the trapper. The attachment once it is trapped is controlled from Single->stepTrappedF
    void stepUntrapped(Hand * cHand);
    void stepTrappedA(Vector force);
    void stepTrappedF(const FiberGrid& grid, Vector const & pos);
    void stepTrappedF_AA();
    
    /// See if in the lattice site, with a different lat_val we find a "target" for trapping. If it is found, verify whether it is not trapped, and then trap it.
    bool try2Trap(Hand * cHand);
    
    void setInteractions(Meca &) const;
    void setInteractionsAF(Meca & meca) const;
    void setInteractionsFA(Meca & meca) const;
    
    void trap(HandMonitor * h);
    void untrap();
    
    bool trapped() const {return trapped_haMon;}
    
    /// Write the reference (to be able to call it from HandMOnitor)
    void writeReference(Outputter& out) const;
    
    /// The point where the virtual trilink is
    Vector3 trapCenter(real exp_shift1,real exp_shift2,real exp_shift3) const;
    
    /// Pointer to the trapped hand
    Hand * trappedHand() const;
    
    void checkAll() const;
    
    real trap_rate() const;
    
    real untrap_rate() const;
    
#ifdef MULTI_LATTICE
    int partnerLattice(Hand * h, Fiber * f) const;
    int checkPartnerLattice() const;
    int checkPartnerLattice(Hand * h) const;
#endif
};
#endif

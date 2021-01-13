// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SWAPPER_H
#define SWAPPER_H

#include "couple.h"

class SwapperProp;

class Swapper : public Couple
{
public:
    /// property
    SwapperProp const* prop;
    // A flag to tell whether a call to attach() or detach() should update the
    // nodelists of couples. If the detachment results from a swap of hands, then
    // it should not.
    bool    has_swapped;
    
protected:

    real    next_swap;
public:
    
    /// create following the specifications in the CoupleProp
    Swapper(SwapperProp const*, Vector const & w = Vector(0,0,0));
    
    /// destructor TODO: THIS HAS TO BE DEFINED PROPERLY
    virtual ~Swapper();
    
    /// Swap the hand
    void swap();
    
    /// Only to add the flipping step (only when double bound)
    void stepAA();
    

    /// Monitor the detachment
    void afterAttachment(Hand const* h);
    
    void beforeDetachment(Hand const* h);
    
    // We have to store and read the property of the hand each time point
    void write(Outputter& out) const;
    
    void read(Inputter & in, Simul& sim, Tag tag);
    
    bool allowAttachment(const FiberBinder & fb);

};
#endif


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef FLIPPER_H
#define FLIPPER_H

#include "couple.h"

class FlipperProp;

class Flipper : public Couple
{
public:
    /// property
    FlipperProp const* prop;
    // A flag to tell whether a call to attach() or detach() should update the
    // nodelists of couples. If the detachment results from a swap of hands, then
    // it should not.
    bool    has_flipped;

    // Rate of flipping that switches along with the hand
    real    flipRate1;
    real    flipRate2;

    // Next flipping (the usual counter)
    real    nextFlip1;
    real    nextFlip2;

protected:
    
//        Hand    * base1;
    Hand    * altHand1;
//        Hand    * base2;
    Hand    * altHand2;

public:
    
    /// create following the specifications in the CoupleProp
    Flipper(FlipperProp const*, Vector const & w = Vector(0,0,0));
    
    /// destructor TODO: THIS HAS TO BE DEFINED PROPERLY
    virtual ~Flipper();

    /// Only to add the flipping step (only when double bound)
    void stepAA();

    /// Swap the hands
    virtual void swapHands1();

    virtual void swapHands2();

    /// Monitor the detachment
    void afterAttachment(Hand const* h);

    void beforeDetachment(Hand const* h);

    // We have to store and read the property of the hand each time point
    void write(Outputter& out) const;

    void read(Inputter & in, Simul& sim, Tag tag);

    virtual bool detachment_after_stepAA(Hand * & cHand_main,Hand * & altHand_main, HandProp const* prop_main, void (Flipper::*swapHandsMain)(), Hand * & cHand_other,Hand * & altHand_other, HandProp const* prop_other, void (Flipper::*swapHandsOther)());
    
    bool hand1Flipped();
        
    bool hand2Flipped();
};
#endif


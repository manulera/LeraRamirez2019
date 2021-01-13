// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "flipper.h"
#include "flipper_prop.h"
#include "exceptions.h"
#include "random.h"
#include "modulo.h"
#include "meca.h"
#include "couple_set.h"
#include "simul.h"
#include "sim.h"

extern Random RNG;
extern Modulo const* modulo;

//------------------------------------------------------------------------------
Flipper::Flipper(FlipperProp const* p, Vector const& w)
: Couple(p, w), prop(p)
{
    has_flipped = false;
    if (p->hand1_flips)
        altHand1 = prop->hand1_alter_prop->newHand(this);
    if (p->hand2_flips)
        altHand2 = prop->hand2_alter_prop->newHand(this);

    flipRate1 = prop->flip_rate1_for_dt;
    flipRate2 = prop->flip_rate2_for_dt;
    
}


Flipper::~Flipper()
{
    prop = 0;
}

//void attach_without_monitor(Fiber * f, real abs)
//{
//    FiberBinder(f, abs);
//}

void Flipper::swapHands1()
{
    if (cHand1->attached())
    {
        has_flipped = true;
        real abs1 = cHand1->abscissa();
        Fiber * f1 = cHand1->fiber();
        
        altHand1->attachTo(f1, abs1);
        altHand1->resetTimers();
        cHand1->detach();
        has_flipped = false;
    }
    
    std::swap(cHand1, altHand1);
    bool unflipped = ( cHand1->prop == Couple::prop->hand1_prop );
    flipRate1 = unflipped ? prop->flip_rate1_for_dt : prop->flip_rate1_bak_dt;
    nextFlip1 = RNG.exponential();
}

void Flipper::swapHands2()
{
    if (cHand2->attached())
    {
        has_flipped = true;
        real abs1 = cHand2->abscissa();
        Fiber * f1 = cHand2->fiber();
        
        altHand2->attachTo(f1, abs1);
        altHand2->resetTimers();
        cHand2->detach();
        has_flipped = false;
    }
    
    std::swap(cHand2, altHand2);
    bool unflipped = ( cHand2->prop == Couple::prop->hand2_prop );
    flipRate2 = unflipped ? prop->flip_rate2_for_dt : prop->flip_rate2_bak_dt;
    nextFlip2 = RNG.exponential();
}

bool Flipper::detachment_after_stepAA(Hand * & cHand_main,Hand * & altHand_main, HandProp const* prop_main, void (Flipper::*swapHandsMain)(), Hand * & cHand_other,Hand * & altHand_other, HandProp const* prop_other, void (Flipper::*swapHandsOther)())
{
    bool main_detached = !cHand_main->attached();
    // If the hand we are looking at was detached
    if (main_detached)
    {
        // And it was in the alter mode
        if (cHand_main->prop!=prop_main)
            // Set it to the base (the base should be detached, since it is set to detach by SwapHands())
            (this->*swapHandsMain)();
        // If the other one was in the alter mode
        if (cHand_other->prop!=prop_other)
        {
            (this->*swapHandsOther)();
        }
    }
    return main_detached;
}

void Flipper::stepAA()
{
    Couple::stepAA();
    // There should not be unattached hands in the form of "alter", since the swap step is only called
    // in this function. For the same reason, here is the place to check whether both hands are attached
    // and if not revert them to the base state. Not further checks should be required.
    
    // This is a bit ugly but better than duplicating code I guess
    if ( detachment_after_stepAA(cHand1, altHand1, prop->hand1_prop,& Flipper::swapHands1, cHand2, altHand2, prop->hand2_prop, & Flipper::swapHands2)||
         detachment_after_stepAA(cHand2, altHand2, prop->hand2_prop,& Flipper::swapHands2, cHand1, altHand1, prop->hand1_prop, & Flipper::swapHands1))
        return;
    
    // Testing the swap of hand 1
    if (prop->hand1_flips)
    {
        nextFlip1-=flipRate1;
        if (nextFlip1<0)
            swapHands1();
    }
    if (prop->hand2_flips)
    {
        nextFlip2-=flipRate2;
        if (nextFlip2<0)
            swapHands2();
    }
}

void Flipper::afterAttachment(Hand const* h)
{
    // This should only be done if the attachment is from the solution, not when switching hands
    if (!has_flipped) {
        Couple::afterAttachment(h);
        if (linking()) {
            if (prop->hand1_flips)
                nextFlip1 = RNG.exponential();
            if (prop->hand2_flips)
                nextFlip2 = RNG.exponential();
        }
    }
}
void Flipper::beforeDetachment(Hand const* h)
{
    // This should only be done if the attachment is from the solution, not when switching hands
    if (!has_flipped) {
        Couple::beforeDetachment(h);
    }
}


void Flipper::write(Outputter& out) const
{
    Couple::write(out);
    // Property of hand 1
    out.writeUInt16(cHand1->prop->index());
    // Property of hand 2
    out.writeUInt16(cHand2->prop->index());
}

void Flipper::read(Inputter & in, Simul& sim, Tag tag)
{
    Couple::read(in, sim, tag);
    try {
        cHand1->prop = sim.findProperty<HandProp*>("hand",in.readUInt16());
        cHand2->prop = sim.findProperty<HandProp*>("hand",in.readUInt16());
    }
    catch( Exception & e ) {
        e << ", in Flipper::read()";
        throw;
    }
}
bool Flipper::hand1Flipped()
{
    return cHand1->prop!=prop->hand1_prop;
}
bool Flipper::hand2Flipped()
{
    return cHand2->prop!=prop->hand2_prop;
}




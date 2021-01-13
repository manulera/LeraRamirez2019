// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "swapper.h"
#include "swapper_prop.h"
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
Swapper::Swapper(SwapperProp const* p, Vector const& w)
: Couple(p, w), prop(p)
{
    has_swapped = false;
    next_swap = 0;
}


Swapper::~Swapper()
{
    prop = 0;
}


void Swapper::swap()
{
    FiberBinder fb1 = FiberBinder(cHand1->fiber(), cHand1->abscissa());
    FiberBinder fb2 = FiberBinder(cHand2->fiber(), cHand2->abscissa());
    // Flag needs to be set before testing attachment allowance
    has_swapped = true;
    if (cHand1->attachmentAllowed(fb2) &&
        cHand2->attachmentAllowed(fb1))
    {
        cHand1->detach();
        cHand2->detach();
        cHand1->attach(fb2);
        cHand2->attach(fb1);
    }
    has_swapped = false;

}

void Swapper::stepAA()
{
    next_swap-=prop->swap_rate_dt;
    if (next_swap<0)
    {
        
        swap();
        next_swap = RNG.exponential();
    }

    // Do it before so that it does not detach
    Couple::stepAA();
    
}

void Swapper::afterAttachment(Hand const* h)
{
    // This should only be done if the attachment is from the solution, not when switching hands
    if (!has_swapped) {
        Couple::afterAttachment(h);
        if (linking()) {
            next_swap = RNG.exponential();
        }
    }
}
void Swapper::beforeDetachment(Hand const* h)
{
    // This should only be done if the attachment is from the solution, not when switching hands
    if (!has_swapped) {
        Couple::beforeDetachment(h);
    }
}


void Swapper::write(Outputter& out) const
{
    Couple::write(out);
    // Property of hand 1
    out.writeUInt16(cHand1->prop->index());
    // Property of hand 2
    out.writeUInt16(cHand2->prop->index());
}

void Swapper::read(Inputter & in, Simul& sim, Tag tag)
{
    Couple::read(in, sim, tag);
    try {
        cHand1->prop = sim.findProperty<HandProp*>("hand",in.readUInt16());
        cHand2->prop = sim.findProperty<HandProp*>("hand",in.readUInt16());
    }
    catch( Exception & e ) {
        e << ", in Swapper::read()";
        throw;
    }
}

bool Swapper::allowAttachment(const FiberBinder & fb)
{
    if (has_swapped)
        return true;
    else
        return Couple::allowAttachment(fb);
}

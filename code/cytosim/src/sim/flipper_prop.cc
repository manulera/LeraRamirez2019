#include "dim.h"
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"
#include "common.h"
#include "simul_prop.h"
#include "flipper_prop.h"
#include "flipper.h"

void FlipperProp::clear()
{
    CoupleProp::clear();
    flip_rate1_for    = 0;
    flip_rate1_bak    = 0;
    flip_rate2_for    = 0;
    flip_rate2_bak    = 0;

    hand1_alter_prop   = 0;
    hand2_alter_prop   = 0;
    
    hand1_alter       = "";
    hand2_alter       = "";
    
    hand1_flips       = false;
    hand2_flips       = false;
    
}

void FlipperProp::read(Glossary& glos)
{
    CoupleProp::read(glos);
    

    glos.set(flip_rate1_for,           "hand1",1);
    if (flip_rate1_for>0)
    {
        glos.set(hand1_alter,          "hand1",2);
        glos.set(flip_rate1_bak,          "hand1",3);
    }
    

    glos.set(flip_rate2_for,           "hand2",1);
    if (flip_rate2_for>0)
    {
        glos.set(hand2_alter,          "hand2",2);
        glos.set(flip_rate2_bak,          "hand2",3);
    }
    
}

void FlipperProp::complete(Simul const* sim)
{
    CoupleProp::complete(sim);
    hand1_flips = flip_rate1_for>0;
    hand2_flips = flip_rate2_for>0;
    
    if (hand1_flips)
    {
        flip_rate1_for_dt = flip_rate1_for * sim->prop->time_step/sim->prop->handmonitor_pace;
        flip_rate1_bak_dt = flip_rate1_bak * sim->prop->time_step/sim->prop->handmonitor_pace;
        if ( hand1_alter.empty() )
            throw InvalidParameter("must define the flip for hand1");
        hand1_alter_prop = static_cast<HandProp*>(sim->properties.find_or_die("hand", hand1_alter));
    }

    if (hand2_flips)
    {
        flip_rate2_for_dt = flip_rate2_for * sim->prop->time_step/sim->prop->handmonitor_pace;
        flip_rate2_bak_dt = flip_rate2_bak * sim->prop->time_step/sim->prop->handmonitor_pace;
        if ( hand2_alter.empty() )
            throw InvalidParameter("must define the flip for hand2");
        hand2_alter_prop = static_cast<HandProp*>(sim->properties.find_or_die("hand", hand2_alter));
    }
}


Couple * FlipperProp::newCouple(Glossary*) const
{
    //std::clog << "CrosslinkProp::newCouple" << std::endl;
    if ( length > 0 )
        throw InvalidParameter("Flipper is not implemented for couples of length different from zero, if you want to make that, create a flipper_long class\n");
    else
        return new Flipper(this);
}


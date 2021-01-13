// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "couple_prop.h"
#include "couple.h"
#include "couple_long.h"
#include "hand_prop.h"
#include "glossary.h"
#include "property_list.h"
#include "simul_prop.h"
#include <cmath>


/**
 This returns a new Couple if ( prop::length <= 0 ),
 or a CoupleLong if ( prop::length > 0 ).
 */
Couple * CoupleProp::newCouple(Glossary*) const
{
    //std::clog << "CoupleProp::newCouple" << std::endl;
    if ( length > 0 )
        return new CoupleLong(this);
    
    return new Couple(this);
}


void CoupleProp::clear()
{
    hand1             = "";
    hand2             = "";
    hand1_prop        = 0;
    hand2_prop        = 0;
    stiffness         = -1;
    length            = 0;
    diffusion         = 0;
    fast_diffusion    = false;
    trans_activated   = 0;
    stiff             = true;
    specificity       = BIND_ALWAYS;
    activity          = "diffuse";
    
    confine           = CONFINE_INSIDE;
    //confine_stiffness = 0;
    confine_space     = "first";
    confine_space_ptr = 0;
    
    is_gillespie      = false;

}


void CoupleProp::read(Glossary& glos)
{
    glos.set(hand1,           "hand1");
    glos.set(hand2,           "hand2");
    glos.set(stiffness,       "stiffness");
    glos.set(length,          "length");
    
    if ( glos.value_is("diffusion", 0, "fast") )
        fast_diffusion = 1;
    else
        glos.set(diffusion,       "diffusion");
    glos.set(fast_diffusion,  "fast_diffusion");
    
    glos.set(trans_activated, "trans_activated");
    glos.set(stiff,           "stiff");
    glos.set(specificity,     "specificity",
             KeyList<Specificity>("off",          BIND_ALWAYS,
                                  "none",         BIND_ALWAYS,
                                  "orthogonal",   BIND_ORTHOGONAL,
                                  "parallel",     BIND_PARALLEL,
                                  "not_parallel", BIND_NOT_PARALLEL,
                                  "antiparallel", BIND_ANTIPARALLEL,
                                  "not_antiparallel", BIND_NOT_ANTIPARALLEL));
    glos.set(activity,        "activity");
    
    glos.set(confine,         "confine",
             KeyList<Confinement>("off",     CONFINE_OFF,
                                  "on",      CONFINE_ON,
#ifdef BACKWARD_COMPATIBILITY
                                  "none",    CONFINE_OFF,
                                  "surface", CONFINE_ON,
#endif
                                  "inside",  CONFINE_INSIDE));
    
    //glos.set(confine_stiffness, "confine", 1);
    glos.set(confine_space,   "confine", 2);

    //glos.set(confine_stiffness, "confine_stiffness");
    glos.set(confine_space,  "confine_space");
    
#ifdef BACKWARD_COMPATIBILITY
    if ( confine_space == "current" )
        confine_space = "last";
#endif
}


void CoupleProp::complete(Simul const* sim)
{
    confine_space_ptr = sim->findSpace(confine_space);

    if ( confine_space_ptr )
        confine_space = confine_space_ptr->property()->name();

    if ( sim->prop->strict && confine != CONFINE_OFF )
    {
        if ( !confine_space_ptr )
            throw InvalidParameter(name()+":confine_space `"+confine_space+"' was not found");
    }
    
    if ( length < 0 )
        throw InvalidParameter("couple:length must be >= 0");
    
    if ( diffusion < 0 )
        throw InvalidParameter("couple:diffusion must be >= 0");

    /**
     We want for one degree of freedom to fulfill `var(dx) = 2 D dt`
     And we use: dx = diffusion_dt * RNG.sreal()
     Since `sreal()` is uniformly distributed, its variance is 1/3,
     and we need `diffusion_dt^2 = 6 D dt`
     */
    diffusion_dt = sqrt( 6.0 * diffusion * sim->prop->time_step/sim->prop->handmonitor_pace );

    if ( stiffness < 0 )
        throw InvalidParameter("couple:stiffness must be specified and >= 0");
    
    if ( hand1.empty() )
        throw InvalidParameter("couple:hand1 must be defined");
    hand1_prop = static_cast<HandProp*>(sim->properties.find_or_die("hand", hand1));
   
    if ( hand2.empty() )
        throw InvalidParameter("couple:hand2 must be defined");
    hand2_prop = static_cast<HandProp*>(sim->properties.find_or_die("hand", hand2));
    
    if ( sim->prop->strict )
    {
        hand1_prop->checkStiffness(stiffness, length, 2, sim->prop->kT);
        /*
         If the length of a Couple (L) is longer than the attachment range of its hands,
         a Couple would place a pair of Fibers at a distance L, thus preventing further
         Couples from linking these two Fibers.
         In most cases, this is not desirable and physically inconsistent.
         */
        if ( length > hand1_prop->binding_range )
            MSG.warning("hand1:binding_range should be >= couple:length for "+name());

        if ( hand2_prop != hand1_prop )
        {
            hand2_prop->checkStiffness(stiffness, length, 2, sim->prop->kT);
        
            if ( length > hand2_prop->binding_range )
                MSG.warning("hand2:binding_range should be >= couple:length for "+name());
        }
    }
    //Added by manu
    if (hand1_prop->is_gillespie||hand2_prop->is_gillespie)
        is_gillespie=true;
    
    hand1_prop->complete_from_couple(sim,stiffness);
    hand2_prop->complete_from_couple(sim,stiffness);
}



void CoupleProp::write_values(std::ostream & os) const
{
    write_value(os, "hand1",           hand1);
    write_value(os, "hand2",           hand2);
    write_value(os, "stiffness",       stiffness);
    write_value(os, "length",          length);
    write_value(os, "diffusion",       diffusion);
    write_value(os, "fast_diffusion",  fast_diffusion);
    write_value(os, "trans_activated", trans_activated);
    write_value(os, "stiff",           stiff);
    write_value(os, "specificity",     specificity);
    write_value(os, "confine",         confine, "", confine_space);
    write_value(os, "activity",        activity);
    write_value(os, "is_gillespie",    is_gillespie);
}


real CoupleProp::spaceVolume() const
{
    if ( confine_space_ptr == 0 )
        throw InvalidParameter("no couple:confinement defined");
    
    real volume = confine_space_ptr->volume();
    
    if ( volume <= 0 )
        throw InvalidParameter("couple:confinement has null volume");
    
    return volume;
}



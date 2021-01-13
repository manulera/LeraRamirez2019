// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "single_prop.h"
#include "glossary.h"
#include "messages.h"

#include "property_list.h"
#include "simul_prop.h"
#include "hand_prop.h"
#include "single.h"
#include "wrist.h"
#include "wrist_long.h"
#include "picket.h"
#include "picket_long.h"

/**
 @defgroup SingleGroup Single and Derived Activities
 @ingroup ObjectGroup
 @ingroup NewObject
 @brief A Single contains one Hand, and can thus bind to one Fiber.

 List of classes accessible by specifying single:activity:
 
 `activity`    |   Class           | Description
 --------------|-------------------|--------------------------------------------
 `diffuse`     | Single            | a single Hand that is mobile (default)
 `fixed`       | Picket PicketLong | a single Hand anchored at a fixed position

 Another class Wrist is used automatically to anchor a Single to a Mecable.
 
 Example:
 @code
 set single grafted
 {
   hand = kinesin
   stiffness = 100
   activity = fixed
 } 
 @endcode
 */
Single * SingleProp::newSingle() const
{
    //std::clog << "SingleProp::newSingle" << std::endl;
    if ( activity == "fixed" )
    {
        if ( length > 0 )
            return new PicketLong(this);
        else
            return new Picket(this);
    }
    else if ( activity == "diffuse" )
    {
        return new Single(this);
    }
    else 
        throw InvalidParameter("unknown Single activity `"+activity+"'");
    return new Single(this);
}


/**
 The Wrist requires a anchor point to be created
 */
Wrist * SingleProp::newWrist(Mecable const* mec, const unsigned point) const
{
    //std::clog << "SingleProp::newWrist" << std::endl;
    if ( length > 0 )
        return new WristLong(this, mec, point);
    else
        return new Wrist(this, mec, point);
}

//------------------------------------------------------------------------------
#pragma mark -

void SingleProp::clear()
{
    hand              = "";
    hand_prop         = 0;
    stiffness         = 0;
    length            = 0;
    diffusion         = 0;
    fast_diffusion    = false;
    activity          = "diffuse";
    
    confine           = CONFINE_INSIDE;
    //confine_stiffness = 0;
    confine_space     = "first";
    confine_space_ptr = 0;
}


void SingleProp::read(Glossary& glos)
{
    glos.set(hand,           "hand");
    glos.set(stiffness,      "stiffness");
    glos.set(length,         "length");
    if ( glos.value_is("diffusion", 0, "fast") )
        fast_diffusion = 1;
    else
        glos.set(diffusion,  "diffusion");
    glos.set(fast_diffusion, "fast_diffusion");
    glos.set(activity,       "activity");

    glos.set(confine,        "confine",
             KeyList<Confinement>("off",     CONFINE_OFF,
                                  "on",      CONFINE_ON,
#ifdef BACKWARD_COMPATIBILITY
                                  "none",    CONFINE_OFF,
                                  "surface", CONFINE_ON,
#endif
                                  "inside",  CONFINE_INSIDE));
    //glos.set(confine_stiffness,  "confine", 1);
    glos.set(confine_space,  "confine", 2);

#ifdef BACKWARD_COMPATIBILITY
    if ( confine_space == "current" )
        confine_space = "last";
#endif
}


void SingleProp::complete(Simul const* sim)
{
    confine_space_ptr = sim->findSpace(confine_space);
    
    if ( confine_space_ptr )
        confine_space = confine_space_ptr->property()->name();

    if ( sim->prop->strict && confine != CONFINE_OFF )
    {
        if ( !confine_space_ptr )
            throw InvalidParameter(name()+":confine_space `"+confine_space+"' was not found");
    }

    if ( hand.empty() )
        throw InvalidParameter("single:hand must be defined");
    hand_prop = static_cast<HandProp*>(sim->properties.find_or_die("hand", hand));
    
    if ( hand_prop == 0 )
        throw InvalidParameter("unknown single:hand '"+hand+"'");

    if ( diffusion < 0 )
        throw InvalidParameter("single:diffusion must be >= 0");

    /**
     We want for one degree of freedom to fulfill `var(dx) = 2 D dt`
     And we use: dx = diffusion_dt * RNG.sreal()
     Since `sreal()` is uniformly distributed, its variance is 1/3,
     and we need `diffusion_dt^2 = 6 D dt`
     */
    diffusion_dt = sqrt( 6.0 * diffusion * sim->prop->time_step/sim->prop->handmonitor_pace );
    
    if ( stiffness < 0 )
        throw InvalidParameter("single:stiffness must be >= 0");

    if ( length < 0 )
        throw InvalidParameter("single:length must be >= 0");

    if ( stiffness > 0  &&  sim->prop->strict )
    {
        hand_prop->checkStiffness(stiffness, length, 1, sim->prop->kT);
        
        /*
         If the length of a Single (L) is longer than the attachment range of its hands,
         a Couple would place a pair of Fibers at a distance L, thus preventing further
         Singles from linking these two Fibers.
         In most cases, this is not desirable and physically inconsistent.
         */
        if ( length > hand_prop->binding_range )
            throw InvalidParameter("hand:binding_range must be >= single:length");
            //MSG.warning("Attachment prevented because single:length > hand:binding_range\n");
    }
}


//------------------------------------------------------------------------------

void SingleProp::write_values(std::ostream & os) const
{
    write_value(os, "hand",           hand);
    write_value(os, "stiffness",      stiffness);
    write_value(os, "length",         length);
    write_value(os, "diffusion",      diffusion);
    write_value(os, "fast_diffusion", fast_diffusion);
    write_value(os, "confine",        confine, "", confine_space);
    write_value(os, "activity",       activity);
}


//------------------------------------------------------------------------------

real SingleProp::spaceVolume() const
{
    if ( confine_space_ptr == 0 )
        throw InvalidParameter("no single:confinement defined");
    
    real volume = confine_space_ptr->volume();
    
    if ( volume <= 0 )
        throw InvalidParameter("single:confinement has null volume");
    
    return volume;
}

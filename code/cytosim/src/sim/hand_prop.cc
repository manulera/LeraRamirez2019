// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "hand_prop.h"
#include "smath.h"
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"
#include "common.h"
#include "sim.h"
#include "property_list.h"
#include "simul_prop.h"
#include "hand.h"
#include "hand_monitor.h"

#include "motor_prop.h"
#include "walker_prop.h"
#include "ase_walker_prop.h"
#include "slider_prop.h"
#include "nucleator_prop.h"
#include "regulator_prop.h"
#include "rescuer_prop.h"
#include "tracker_prop.h"
#include "cutter_prop.h"
#include "chewer_prop.h"
#include "mighty_prop.h"
#include "actor_prop.h"

#ifdef NEW_HANDS
#include "kinesin_prop.h"
#include "dynein_prop.h"
#include "myosin_prop.h"
#endif


/**
 @defgroup HandGroup Hand and Derived Activities
 @ingroup ObjectGroup
 @ingroup NewObject
 @brief A Hand can bind to a Fiber, and derived class can do more things.

 A plain Hand can only bind and unbind from a Fiber.
 Derived classes are available that implement more complex functionalities,
 for example molecular motors or severing enzymes.
 
 List of classes accessible by specifying `hand:activity`:
 
 @ref HandGroup
 
 `activity`    |   Class       | Parameters         |  Property
 --------------|---------------|--------------------|---------------
 `bind`        | Hand          | @ref HandPar       | HandProp
 `move`        | Motor         | @ref MotorPar      | MotorProp
 `nucleate`    | Nucleator     | @ref NucleatorPar  | NucleatorProp
 `slide`       | Slider        | @ref SliderPar     | SliderProp
 `track`       | Tracker       | @ref TrackerPar    | TrackerProp
 `rescue`      | Rescuer       | @ref RescuerPar    | RescuerProp
 `regulate`    | Regulator     | @ref RegulatorPar  | RegulatorProp
 `cut`         | Cutter        | @ref CutterPar     | CutterProp
 `chew`        | Chewer        | @ref ChewerPar     | ChewerProp
 `mighty`      | Mighty        | @ref MightyPar     | MightyProp
 `act`         | Actor         | @ref ActorPar      | ActorProp
 
 <h2>Digital Hands:</h2>
 
 `activity`    |   Class       | Parameters         |  Property
 --------------|---------------|--------------------|---------------
 `digit`       | Digit         | @ref DigitPar      | DigitProp
 `walk`        | Walker        | @ref WalkerPar     | WalkerProp
 `kinesin`*    | Kinesin       | @ref KinesinPar    | KinesinProp
 `dynein`*     | Dynein        | @ref DyneinPar     | DyneinProp
 `myosin`*     | Myosin        | @ref MyosinPar     | MyosinProp
 
 * Unfinished classes.
 
 Example:
 @code
 set hand motor
 {
   binding = 10, 0.05
   unbinding = 0.2, 3
 
   activity = move
   max_speed = 1
   stall_force = 5
 } 
 @endcode
 */
HandProp * HandProp::newProperty(const std::string& name, Glossary& glos)
{
    HandProp * hp = 0;
    
    std::string a;
    if ( glos.peek(a, "activity") )
    {
        
        if ( a == "move" || a == "motor" )
            hp = new MotorProp(name);
        else if ( a == "digit" )
            hp = new DigitProp(name);
        else if ( a == "walk" )
            hp = new WalkerProp(name);
        //ADDED
        else if ( a == "ase_walk" )
            hp = new AseWalkerProp(name);
        else if ( a == "slide" )
            hp = new SliderProp(name);
#ifdef NEW_HANDS
        else if ( a == "kinesin" )
            hp = new KinesinProp(name);
        else if ( a == "dynein" )
            hp = new DyneinProp(name);
        else if ( a == "myosin" )
            hp = new MyosinProp(name);
#endif
        else if ( a == "nucleate" )
            hp = new NucleatorProp(name);
        else if ( a == "regulate" )
            hp = new RegulatorProp(name);
        else if ( a == "track" )
            hp = new TrackerProp(name);
        else if ( a == "rescue" )
            hp = new RescuerProp(name);
        else if ( a == "cut" )
            hp = new CutterProp(name);
        else if ( a == "chew" )
            hp = new ChewerProp(name);
        else if ( a == "mighty" )
            hp = new MightyProp(name);
        else if ( a == "act" )
            hp = new ActorProp(name);
        else if ( a == "bind" || a == "none" )
            hp = new HandProp(name);
        else
        {
#if ( 1 )
            throw InvalidParameter("unknown hand:activity `"+a+"'");
#else
            // try to proceed with an incorrect Object:
            hp = new HandProp(name);
            std::cerr << "WARNING: unknown hand:activity `" << a << "'" << std::endl;
#endif
        }
    }
    else
        hp = new HandProp(name);
    
    return hp;
}


Hand * HandProp::newHand(HandMonitor* h) const
{
    return new Hand(this, h);
}


//------------------------------------------------------------------------------
void HandProp::clear()
{
    binding_rate       = 0;
    binding_range      = 0;
    binding_key        = (~0);  //all bits at 1
    
    unbinding_rate     = 0;
    unbinding_force    = INFINITY;
    unbinding_force_inv = 0;

    bind_also_ends     = false;
    bind_only_end      = NO_END;
    bind_end_range     = 0;
	bind_only_free_end = false;
    hold_growing_end   = false;
    hold_shrinking_end = false;
    
    activity           = "bind";
    display            = "";
    display_fresh      = false;
    is_gillespie       = false;
    lat_val            = 0;
    overlap_affinity   = 0;
}


void HandProp::read(Glossary& glos)
{
    glos.set(binding_rate,       "binding_rate");
    glos.set(binding_range,      "binding_range");
    glos.set(binding_key,        "binding_key");
    //alternative syntax:
    glos.set(binding_rate,       "binding", 0);
    glos.set(binding_range,      "binding", 1);
    glos.set(binding_key,        "binding", 2);
    
    glos.set(unbinding_rate,     "unbinding_rate");
    glos.set(unbinding_force,    "unbinding_force");
    //alternative syntax:
    glos.set(unbinding_rate,     "unbinding", 0);
    glos.set(unbinding_force,    "unbinding", 1);
    
    glos.set(bind_also_ends, "bind_also_ends") || glos.set(bind_also_ends, "bind_also_end");

	glos.set(bind_only_end,      "bind_only_end",
             KeyList<FiberEnd>("none",        NO_END,
                               "plus_end",    PLUS_END,
                               "minus_end",   MINUS_END,
                               "both_ends",   BOTH_ENDS));
    glos.set(bind_end_range,     "bind_only_end", 1);
    glos.set(bind_end_range,     "bind_end_range");

#ifdef BACKWARD_COMPATIBILITY
    glos.set(bind_only_end,      "bind_end",
             KeyList<FiberEnd>("none",        NO_END,
                               "plus_end",    PLUS_END,
                               "minus_end",   MINUS_END,
                               "both_ends",   BOTH_ENDS));
    glos.set(bind_end_range,     "bind_end", 1);
#endif
    
    
    glos.set(hold_growing_end,   "hold_growing_end");
    glos.set(hold_shrinking_end, "hold_shrinking_end");
    glos.set(bind_only_free_end, "bind_only_free_end");
    
    glos.set(activity,           "activity");
    if ( glos.set(display, "display") )
        display_fresh = true;
    
#ifdef BACKWARD_COMPATIBILITY
    if ( glos.set(hold_growing_end, "hold_growing_ends") )
        MSG.warning("hand:hold_growing_ends was renamed hold_growing_end\n");
#endif
    
    //Added by manu
    glos.set(is_gillespie,           "is_gillespie");
    glos.set(overlap_affinity, "overlap_affinity");
}


void HandProp::complete(Simul const* sim)
{    
    if ( sim==0 || sim->prop->time_step < REAL_EPSILON )
        throw InvalidParameter("simul:time_step is not defined");
    
    binding_range_sqr = binding_range * binding_range;
    binding_rate_dt   = binding_rate * sim->prop->time_step/sim->prop->handmonitor_pace;
    unbinding_rate_dt = unbinding_rate * sim->prop->time_step/sim->prop->handmonitor_pace;
    
    binding_rate_dt_4 = 4 * binding_rate_dt;
    
    if ( binding_range < 0 )
        throw InvalidParameter(name()+":binding_range must be >= 0");
    
    if ( binding_rate < 0 )
        throw InvalidParameter(name()+":binding_rate must be positive");
    
    if ( unbinding_rate < 0 )
        throw InvalidParameter(name()+":unbinding_rate must be positive");
    
    if ( sim && sim->prop->strict )
    {
        if ( binding_rate_dt_4 > sim->prop->acceptable_rate )
            MSG.warning(name()+":binding_rate is too high: decrease time_step\n");
    
        if ( unbinding_rate_dt > sim->prop->acceptable_rate )
            MSG.warning(name()+":unbinding_rate is too high: decrease time_step\n");
    }
    
#ifdef BACKWARD_COMPATIBILITY
    if ( unbinding_force == 0 )
    {
        MSG.warning("assuming that hand:unbinding_force=+inf, because the set value was zero\n");
        unbinding_force = INFINITY;
    }
#endif
    
    if ( unbinding_force <= 0 )
        throw InvalidParameter(name()+":unbinding_force must be > 0");

    // this should be zero if 'unbinding_force = inf':
    unbinding_force_inv = 1.0 / unbinding_force;

    /*
     The exponential term in Kramer's theory can easily become numerically "infinite",
     and since `zero * infinite` is undefined, we disable here
     the exponential argument if the unbinding rate is null:
     */
    if ( unbinding_rate == 0 )
        unbinding_force_inv = 0;

    //std::clog << name() << " unbinding_force_inv = " << unbinding_force_inv << std::endl;
}



/**
 Compare the energy in a link when it binds at its maximum distance,
 with the Thermal energy
 
 @todo the warning may not be relevant for long Links
 */
void HandProp::checkStiffness(real stiff, real len, real mul, real kT) const
{
    real dis = binding_range - len;
    real en = ( stiff * dis * dis ) / kT;
    
    if ( en > 10.0 )
    {
        std::ostringstream oss;
        oss << "hand `" << name() << "' overcomes high energy when binding:\n";
        oss << PREF << "stiffness * binding_range^2 = " << en << " kT\n";
        //oss << PREF << "you could decrease stiffness or binding_range" << std::endl;
        MSG.warning(oss.str());
    }
    
    
    real ap = exp( stiff * dis * unbinding_force_inv );
    
    if ( ap > 10.0 )
    {
        std::ostringstream oss;
        oss << "hand `" << name() << "' may unbind just after binding:\n";
        oss << PREF << "exp( stiffness * binding_range / unbinding_force ) = " << ap << "\n";
        //oss << PREF << "you might want decrease stiffness or binding_range" << std::endl;
        MSG.warning(oss.str());
    }
}


void HandProp::write_values(std::ostream & os) const
{
    write_value(os, "binding",            binding_rate, binding_range);
    write_value(os, "binding_key",        binding_key);
    write_value(os, "unbinding",          unbinding_rate, unbinding_force);

    write_value(os, "bind_also_ends",     bind_also_ends);
    write_value(os, "hold_growing_end",   hold_growing_end);
    write_value(os, "hold_shrinking_end", hold_shrinking_end);
	write_value(os, "bind_only_end",      bind_only_end, bind_end_range);
	write_value(os, "bind_only_free_end", bind_only_free_end);
	
    write_value(os, "display",            "("+display+")");
    write_value(os, "activity",           activity);
    write_value(os, "is_gillespie",       is_gillespie);
}


/**
 Estimate attachment propensity per unit length of fiber
 */
real HandProp::bindingSection(bool with_time_step) const
{
    real rate = ( with_time_step ?  binding_rate_dt : binding_rate );

#if ( DIM == 2 )
    return 2 * binding_range * rate;
#elif ( DIM == 3 )
    return M_PI * binding_range * binding_range * rate;
#else
    return 0;
#endif
}



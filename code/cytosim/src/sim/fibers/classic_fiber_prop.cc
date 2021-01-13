// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cmath>
#include "sim.h"
#include "messages.h"
#include "classic_fiber_prop.h"
#include "classic_fiber.h"
#include "property_list.h"
#include "simul_prop.h"
#include "exceptions.h"
#include "glossary.h"


Fiber* ClassicFiberProp::newFiber() const
{
    return new ClassicFiber(this);
}


void ClassicFiberProp::clear()
{
    FiberProp::clear();
    
    for ( int i = 0; i < 2; ++i )
    {
        growing_speed[i]            = 0;
        growing_off_speed[i]        = 0;
        growing_force[i]            = INFINITY;
        shrinking_speed[i]          = 0;
        catastrophe_rate[i]         = 0;
        catastrophe_rate_stalled[i] = 0;
        catastrophe_coef[i]         = 0;
        rescue_rate[i]              = 0;
        rebirth_rate[i]             = 0;
    }
    
#ifdef NEW_LENGTH_DEPENDENT_CATASTROPHE
    catastrophe_length  = 0;
#endif
#ifdef NEW_CATASTROPHE_OUTSIDE
    catastrophe_outside = false;
#endif
    persistent          = false;
}


void ClassicFiberProp::read(Glossary& glos)
{
    FiberProp::read(glos);
    
    glos.set(growing_speed,            2, "growing_speed");
    glos.set(growing_off_speed,        2, "growing_off_speed");
    glos.set(growing_force,            2, "growing_force");
    glos.set(shrinking_speed,          2, "shrinking_speed");
    glos.set(catastrophe_rate,         2, "catastrophe_rate");
    glos.set(catastrophe_rate_stalled, 2, "catastrophe_rate_stalled");
    glos.set(rescue_rate,              2, "rescue_rate");
    glos.set(rebirth_rate,             2, "rebirth_rate");
    glos.set(persistent,                  "persistent");

#ifdef NEW_LENGTH_DEPENDENT_CATASTROPHE
    glos.set(catastrophe_length,       "catastrophe_length");
#endif
#ifdef NEW_CATASTROPHE_OUTSIDE
    glos.set(catastrophe_outside,      "catastrophe_outside");
#endif


#ifdef BACKWARD_COMPATIBILITY
    
    if ( glos.set(growing_force[0], "dynamic_force") )
        MSG.warning("fiber:dynamic_force was renamed growing_force\n");

    // change made 13 October 2014
    if ( glos.nb_values("growing_force") == 1 )
    {
        real x = 0;
        if ( glos.peek(x, "catastrophe_rate", 1) && x > 0 )
        {
            MSG.warning("catastrophe_rate[1] was renamed catastrophe_rate_stalled[0]\n");
            catastrophe_rate_stalled[0] = x;
            catastrophe_rate[1] = 0;
        }
        if ( glos.peek(x, "growing_speed", 1) && x != 0 )
            throw InvalidParameter("fiber:growing_speed[1] was renamed growing_off_speed[0]\n");
    }

    int f = 0;
    if ( glos.set(f, "fate", KeyList<int>("none", 0, "destroy", 1, "rescue", 2)))
    {
        MSG.warning("fiber:fate is deprecated: use `persistent` and `rebirth_rate`\n");
        persistent = ( f != 1 );
        rebirth_rate[0] = ( f == 2 ? INFINITY : 0 );
    }
    
    bool ds;
    if ( glos.set(ds, "delete_stub") )
    {
        MSG.warning("please use `persistent` instead of `delete_stub`\n");
        persistent = !ds;
    }

#endif
}


void ClassicFiberProp::complete(Simul const* sim)
{
    FiberProp::complete(sim);
    
    for ( int i = 0; i < 2; ++i )
    {
        if ( growing_speed[i] < 0 )
            throw InvalidParameter("fiber:growing_speed should be >= 0");
        growing_speed_dt[i] = growing_speed[i] * sim->prop->time_step;

        if ( growing_off_speed[i] > 0 )
            throw InvalidParameter("growing_off_speed should be <= 0");

        if ( growing_speed[i] + growing_off_speed[i] < 0 )
            throw InvalidParameter("fiber:growing_speed+growing_off_speed should be >= 0");
        growing_off_speed_dt[i] = growing_off_speed[i] * sim->prop->time_step;

        if ( growing_force[i] <= 0 )
            throw InvalidParameter("fiber:growing_force should be > 0");

        if ( shrinking_speed[i] > 0 )
            throw InvalidParameter("fiber:shrinking_speed should be <= 0");
        shrinking_speed_dt[i]  = shrinking_speed[i] * sim->prop->time_step;

        if ( catastrophe_rate[i] < 0 )
            throw InvalidParameter("fiber:catastrophe_rate should be >= 0");
        catastrophe_rate_dt[i] = catastrophe_rate[i] * sim->prop->time_step;

        if ( catastrophe_rate_stalled[i] < 0 )
            throw InvalidParameter("fiber:catastrophe_rate_stalled should be >= 0");
        catastrophe_rate_stalled_dt[i] = catastrophe_rate_stalled[i] * sim->prop->time_step;

        catastrophe_coef[i] = 0;
        
        if ( catastrophe_rate_stalled[i] > 0 && catastrophe_rate[i] > 0 )
        {
            if ( catastrophe_rate[i] > catastrophe_rate_stalled[i] )
                throw InvalidParameter("fiber:catastrophe_rate_stalled must be greater than catastrophe_rate");

            if ( sim->prop->strict  &&  growing_speed[i] + growing_off_speed[i] <= 0 )
                MSG.warning("fiber:growing_speed + growing_off_speed <= 0\n");
                
            catastrophe_coef[i] = ( catastrophe_rate_stalled[i]/catastrophe_rate[i] - 1.0 )
            / ( ( growing_speed[i] + growing_off_speed[i] ) * sim->prop->time_step );
            
            if ( catastrophe_coef[i] < 0 )
                throw InvalidParameter("inconsistent fiber:dynamic parameters");
        }
        
        if ( rescue_rate[i] < 0 )
            throw InvalidParameter("fiber:rescue_rate should be >= 0");
        rescue_prob[i] = 1 - exp( -rescue_rate[i] * sim->prop->time_step );

        if ( rebirth_rate[i] < 0 )
            throw InvalidParameter("fiber:rebirth_rate should be >= 0");
        rebirth_prob[i] = 1 - exp( -rebirth_rate[i] * sim->prop->time_step );
    }
}


void ClassicFiberProp::write_values(std::ostream & os) const
{
    FiberProp::write_values(os);
    
    write_value(os, "growing_speed",            growing_speed, 2);
    write_value(os, "growing_off_speed",        growing_off_speed, 2);
    write_value(os, "growing_force",            growing_force, 2);
    write_value(os, "shrinking_speed",          shrinking_speed, 2);
    write_value(os, "catastrophe_rate",         catastrophe_rate, 2);
    write_value(os, "catastrophe_rate_stalled", catastrophe_rate_stalled, 2);
    write_value(os, "rescue_rate",              rescue_rate, 2);
    write_value(os, "rebirth_rate",             rebirth_rate, 2);
    write_value(os, "persistent",               persistent);
#ifdef NEW_LENGTH_DEPENDENT_CATASTROPHE
    write_value(os, "catastrophe_length",       catastrophe_length);
#endif
#ifdef NEW_CATASTROPHE_OUTSIDE
    write_value(os, "catastrophe_outside",      catastrophe_outside);
#endif
}


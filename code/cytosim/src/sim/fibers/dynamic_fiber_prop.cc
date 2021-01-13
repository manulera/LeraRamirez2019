// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cmath>
#include "sim.h"
#include "dynamic_fiber_prop.h"
#include "dynamic_fiber.h"
#include "property_list.h"
#include "simul_prop.h"
#include "exceptions.h"
#include "glossary.h"


Fiber* DynamicFiberProp::newFiber() const
{
    return new DynamicFiber(this);
}


void DynamicFiberProp::clear()
{
    FiberProp::clear();
    
    unit_length    = 0.008;
    persistent     = false;
    zone_radius    = INFINITY;
    zone_space     = "";
    
    for ( int i = 0; i < 2; ++i )
    {
        growing_speed[i]        = 0;
        growing_off_speed[i]    = 0;
        growing_force[i]        = INFINITY;
        hydrolysis_rate[i]      = 0;
        zone_hydrolysis_rate[i] = 0;
        shrinking_speed[i]      = 0;
        rebirth_rate[i]         = 0;
    }
}


void DynamicFiberProp::read(Glossary& glos)
{
    FiberProp::read(glos);
    
    glos.set(unit_length,             "unit_length");
    glos.set(growing_speed,        2, "growing_speed");
    glos.set(growing_off_speed,    2, "growing_off_speed");
    glos.set(growing_force,        2, "growing_force");
    glos.set(hydrolysis_rate,      2, "hydrolysis_rate");
    glos.set(shrinking_speed,      2, "shrinking_speed");
    glos.set(rebirth_rate,         2, "rebirth_rate");
    glos.set(persistent,              "persistent");
    
    glos.set(zone_space,              "zone_space");
    glos.set(zone_radius,             "zone_radius");
    glos.set(zone_hydrolysis_rate, 2, "zone_hydrolysis_rate");

#ifdef BACKWARD_COMPATIBILITY
    
    if ( glos.set(growing_force[0], "dynamic_force") )
        MSG.warning("fiber:dynamic_force was renamed growing_force\n");
    
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


void DynamicFiberProp::complete(Simul const* sim)
{
    FiberProp::complete(sim);

    for ( int i = 0; i < 2; ++i )
    {
        if ( growing_force[i] <= 0 )
            throw InvalidParameter("fiber:growing_force should be > 0");
        
        if ( growing_speed[i] < 0 )
            throw InvalidParameter("fiber:growing_speed should be >= 0");
        growing_rate_dt[i] = sim->prop->time_step * fabs(growing_speed[i]) / unit_length;

        if ( growing_off_speed[i] > 0 )
            throw InvalidParameter("growing_off_speed should be <= 0");
        growing_off_rate_dt[i] = sim->prop->time_step * growing_off_speed[i] / unit_length;

        if ( hydrolysis_rate[i] < 0 )
            throw InvalidParameter("fiber:hydrolysis_rate should be >= 0");
        hydrolysis_rate_2dt[i] = 2 * sim->prop->time_step * hydrolysis_rate[i];
        
        
        zone_space_ptr = sim->findSpace(zone_space);
        
        if ( zone_radius < 0 )
            throw InvalidParameter("fiber:zone_radius should be >= 0");
        if ( zone_hydrolysis_rate[i] < 0 )
            throw InvalidParameter("fiber:zone_hydrolysis_rate should be >= 0");
        zone_hydrolysis_rate_2dt[i] = 2 * sim->prop->time_step * zone_hydrolysis_rate[i];
        
        if ( shrinking_speed[i] > 0 )
            throw InvalidParameter("fiber:shrinking_speed should be <= 0");
        shrinking_rate_dt[i] = sim->prop->time_step * fabs(shrinking_speed[i]) / unit_length;
        
        if ( rebirth_rate[i] < 0 )
            throw InvalidParameter("fiber:rebirth_rate should be >= 0");
        rebirth_prob[i] = 1 - exp( -rebirth_rate[i] * sim->prop->time_step );
    }

    if ( min_length <= 0 )
        min_length = 3 * unit_length;
    
    zone_radius_sqr = zone_radius * zone_radius;
    
    if ( sim->prop->verbose && sim->prop->strict )
    {
        /*
         Using formula from:
         A theory of microtubule catastrophes and their regulation</b>\n
         Brun L, Rupp B, Ward J, Nedelec F\n
         PNAS 106 (50) 21173-21178; 2009\n
         */
        const real h = hydrolysis_rate[0];
        real g = growing_speed[0] / unit_length;
        real ctime = ( 7*h*h + 12*g*h + 3*g*g ) / ( 3*h*h * ( 2*h + 3*g ) );
        //const real ctime = g / ( 3*h*h );  // that is only true if g >> h
        std::clog << std::setprecision(5);
        std::clog << "  DynamicFiber h " << h << " g " << g << " :";
        std::clog << " catastrophe_time " << ctime << "  rate " << 1/ctime;
        std::clog << " length " << g*unit_length*ctime << '\n';
        /*
        ctime = g / ( 3*h*h );  // that is only true if g >> h
        std::clog << "  DynamicFiber h " << h << " g " << g << " :";
        std::clog << " catastrophe_time " << ctime << "  rate " << 1/ctime;
        std::clog << " length " << g*unit_length*ctime << '\n';
         */
        g *= 0.5;
        ctime = ( 7*h*h + 12*g*h + 3*g*g ) / ( 3*h*h * ( 2*h + 3*g ) );
        std::clog << "  DynamicFiber h " << h << " g " << g << " :";
        std::clog << " catastrophe_time " << ctime << "  rate " << 1/ctime;
        std::clog << " length " << g*unit_length*ctime << '\n';
    }
}


void DynamicFiberProp::write_values(std::ostream & os) const
{
    FiberProp::write_values(os);
    
    write_value(os, "unit_length",          unit_length);
    write_value(os, "growing_speed",        growing_speed, 2);
    write_value(os, "growing_off_speed",    growing_off_speed, 2);
    write_value(os, "growing_force",        growing_force, 2);
    write_value(os, "hydrolysis_rate",      hydrolysis_rate, 2);
    write_value(os, "shrinking_speed",      shrinking_speed, 2);
    write_value(os, "rebirth_rate",         rebirth_rate, 2);
    write_value(os, "persistent",           persistent);

    write_value(os, "zone_space",           zone_space);
    write_value(os, "zone_radius",          zone_radius);
    write_value(os, "zone_hydrolysis_rate", zone_hydrolysis_rate, 2);
}


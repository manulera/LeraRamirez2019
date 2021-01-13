// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "fiber_prop.h"
#include "sim.h"
#include "smath.h"
#include <cmath>
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"
#include "property_list.h"
#include "single_prop.h"
#include "simul_prop.h"
#include "simul.h"
#include "fiber.h"

/**
 This is virtualized to return a derived Fiber if appropriate
 */
Fiber* FiberProp::newFiber() const
{
    return new Fiber(this);
}


/**
 @addtogroup FiberGroup
 @{
 <hr>
 
 When creating a new Fiber, you may specify:
 - the initial length,
 - the initial state of the PLUS_END and MINUS_END,
 - if the position refers to the center or to the tip of the fiber
 - the shape, using a set of points
 .
 
 Syntax:
 
 @code
 new fiber ...
 {
   length = REAL, LENGTH_MODIFIER
   end_state = PLUS_END_STATE, MINUS_END_STATE
   reference = REFERENCE
 }
 @endcode
 
 The optional LENGTH_MODIFIER can be:
 - `exponential`,
 - REAL
 .
 This introduces variability, without changing the mean length.
 The second form generates a flat distribution of width 2*LENGTH_MODIFIER.
 
 The initial states PLUS_END_STATE and MINUS_END_STATE can be:
 - 0 = white
 - 1 = green
 - 4 = red
 .
 
 Optional reference specificiation:
 - center [default]
 - plus_end
 - minus_end
 .
 
 To specify the shape of a Fiber directly, use:
 @code
 new fiber ...
 {
   shape = POSITION, POSITION, ...
 }
 @endcode
 
 Examples:
 
 @code
 new fiber ...
 {
   length = 1
   plus_end_state = 1
   minus_end_state = 0
 }
 @endcode

 which is equivalent to:
 @code
 new fiber ...
 {
   length = 1
   end_state = green, white
 }
 @endcode
 
 @code
 new fiber ...
 {
   position = 0 0 0
   orientation = 1 0 0
   shape = -4 -3 0, -3 0 0, -1 2 0, 1  3 0
 }
 @endcode

 @}
 */
Fiber* FiberProp::newFiber(Glossary& opt) const
{
    Fiber * fib = newFiber();
    real len = length;
    
    /* 
     initial length and reference point for placement can be specified in 'opt'
     */
#ifdef BACKWARD_COMPATIBILITY
    opt.set(len, "initial_length");
#endif
    opt.set(len, "length") || opt.set(len, "fiber_length");

    // exponential distribution:
    if ( opt.value_is("length", 1, "exponential") )
    {
        len *= RNG.exponential();
        if ( len < min_length )
            len = min_length;
        if ( len > max_length )
            len = max_length;
    }
    else
    {
        // add variability without changing mean:
        real var = 0;
        if ( opt.set(var, "length", 1) || opt.set(var, "fiber_length",1) )
        {
            len += var * RNG.sreal();
            if ( len < min_length )
                len = min_length;
        }
    }
    
#if ( 1 )
    // specify the model-points directly:
    if ( opt.has_key("points") )
    {
        unsigned nbp = opt.nb_values("points");
        fib->setNbPoints(nbp);
        
        Vector vec(0,0,0);
        for ( unsigned p = 0; p < nbp; ++p )
        {
            if ( ! opt.set(vec, "points", p) )
                throw InvalidParameter("fiber:points must be a list of comma-separated vectors");
            fib->setPoint(p, vec);
        }
        fib->imposeLength(len);
    }
    else
#endif
    if ( opt.has_key("shape") )
    {
        unsigned nbp = opt.nb_values("shape");
        
        if ( nbp < 2 )
            throw InvalidParameter("fiber:shape must be a list of comma-separated vectors");

        real* tmp = new real[DIM*nbp];

        Vector vec(0,0,0);
        for ( unsigned p = 0; p < nbp; ++p )
        {
            if ( ! opt.set(vec, "shape", p) )
                throw InvalidParameter("fiber:shape must be a list of comma-separated vectors");
            vec.put(tmp+DIM*p);
        }
        fib->setShape(tmp, nbp, 0);
        if ( fib->nbPoints() < 2 )
            throw InvalidParameter("the vectors specified in fiber:shape must be different");
        fib->reshape();
        delete[] tmp;
    }
    else
    {        
        FiberEnd ref = CENTER;

        opt.set(ref, "reference",
                KeyList<FiberEnd>("plus_end", PLUS_END, "minus_end", MINUS_END, "center", CENTER));
        
        // initialize points:
        fib->setStraight(Vector(0,0,0), Vector(1,0,0), len, ref);
    }
    
    // possible dynamic states of the ends
    KeyList<int> keys("white",     STATE_WHITE,
                      "green",     STATE_GREEN,
                      "yellow",    STATE_YELLOW,
                      "orange",    STATE_ORANGE,
                      "red",       STATE_RED,
                      "static",    STATE_WHITE,
                      "growing",   STATE_GREEN,
                      "shrinking", STATE_RED);
    
    int s = STATE_WHITE;
    
    if ( opt.set(s, "plus_end_state") || opt.set(s, "end_state", keys) )
        fib->setDynamicStateP(s);
    
    if ( opt.set(s, "minus_end_state") || opt.set(s, "end_state", keys, 1) )
        fib->setDynamicStateM(s);

    fib->update();

    return fib;
}


//------------------------------------------------------------------------------
void FiberProp::clear()
{
    rigidity            = -1;
    segmentation        = 1;
    length              = 1;
    min_length          = 0;
    max_length          = INFINITY;
    total_polymer       = INFINITY;

    viscosity           = -1;
    hydrodynamic_radius[0] = 0.0125;  // radius of a Microtubule
    hydrodynamic_radius[1] = 10;
    surface_effect      = false;
    cylinder_height     = 0;
    
    binding_key         = (~0);  //all bits at 1

    lattice             = 0;
    lattice_unit        = 0;
    lattice_cut_fiber   = 0;
    lattice_flux_speed  = 0;
    lattice_binding_rate = 0;
    lattice_unbinding_rate = 0;
    
    confine             = CONFINE_OFF;
    confine_stiffness   = -1;
    confine_space       = "first";
    confine_space_ptr   = 0;
    
    steric              = 0;
    steric_radius       = 0;
    steric_range        = 0;
    
    field               = "none";
    field_ptr           = 0;
    
    glue                = 0;
    glue_single         = "none";
    glue_prop           = 0;
    glue_set            = 0;
    
    colinear_force      = 0;
    max_chewing_speed   = 0;
    
    activity            = "none";
    display             = "";
    display_fresh       = false;
    
    total_length        = 0;
    free_polymer        = 1;
    time_step           = 0;
    
#ifdef BACKWARD_COMPATIBILITY
    squeeze             = 0;
    squeeze_force       = 0;
    squeeze_range       = 1;
#endif
    
    // ADDED BY MANU
    sweep_digits        = false;
}

//------------------------------------------------------------------------------
void FiberProp::read(Glossary& glos)
{
    glos.set(rigidity,          "rigidity");
    glos.set(segmentation,      "segmentation");
    glos.set(length,            "length");
    glos.set(min_length,        "min_length");
    glos.set(max_length,        "max_length");
    glos.set(total_polymer,     "total_polymer");

    glos.set(viscosity,         "viscosity");
    glos.set(hydrodynamic_radius, 2, "hydrodynamic_radius");
    glos.set(surface_effect,    "surface_effect");
    glos.set(cylinder_height,   "surface_effect", 1);
    
    glos.set(binding_key,       "binding_key");
    
    glos.set(lattice,           "lattice");
    glos.set(lattice_unit,      "lattice", 1);
    
    glos.set(lattice_unit,      "lattice_unit");
    glos.set(lattice_cut_fiber, "lattice_cut_fiber");
    glos.set(lattice_flux_speed,"lattice_flux_speed");
    glos.set(lattice_binding_rate, "lattice_binding_rate");
    glos.set(lattice_unbinding_rate, "lattice_unbinding_rate");
    
    
    glos.set(confine,           "confine",
             KeyList<Confinement>("off",       CONFINE_OFF,
                                  "on",        CONFINE_ON,
                                  "inside",    CONFINE_INSIDE,
                                  "outside",   CONFINE_OUTSIDE,
#ifdef BACKWARD_COMPATIBILITY
                                  "none",      CONFINE_OFF,
                                  "surface",   CONFINE_ON,
#endif
                                  "plus_end",  CONFINE_PLUS_END,
                                  "minus_end", CONFINE_MINUS_END,
                                  "both_ends", CONFINE_BOTH_ENDS,
                                  "plus_out",  CONFINE_PLUS_OUT));
    glos.set(confine_stiffness, "confine", 1);
    glos.set(confine_space,     "confine", 2);

    glos.set(confine_stiffness, "confine_stiffness");
    glos.set(confine_space,     "confine_space");

#ifdef BACKWARD_COMPATIBILITY
    if ( confine_space == "current" )
        confine_space = "last";

    glos.set(confine,           "confined",
             KeyList<Confinement>("none",      CONFINE_OFF,
                                  "inside",    CONFINE_INSIDE,
                                  "outside",   CONFINE_OUTSIDE,
                                  "surface",   CONFINE_ON,
                                  "minus_end", CONFINE_MINUS_END,
                                  "plus_end",  CONFINE_PLUS_END));
    glos.set(confine_stiffness, "confined", 1);

    glos.set(squeeze,           "squeeze");
    glos.set(squeeze_force,     "squeeze", 1);
    glos.set(squeeze_range,     "squeeze", 2);
#endif
    
    glos.set(steric,            "steric");
    glos.set(steric_radius,     "steric", 1);
    glos.set(steric_range,      "steric", 2);
    glos.set(steric_radius,     "steric_radius");
    glos.set(steric_range,      "steric_range");
    
    glos.set(field,             "field");
    glos.set(glue,              "glue");
    glos.set(glue_single,       "glue", 1);
    
    glos.set(colinear_force,    "colinear_force");
    glos.set(max_chewing_speed, "max_chewing_speed");
    
    
    glos.set(activity, "activity");
    if ( glos.set(display, "display") )
        display_fresh = true;
    
    //ADDED BY MANU
    glos.set(sweep_digits, "sweep_digits");
    
}


void FiberProp::complete(Simul const* sim)
{
    time_step = sim->prop->time_step;
    
    if ( viscosity < 0 )
        viscosity = sim->prop->viscosity;
    
    if ( viscosity < 0 )
        throw InvalidParameter("fiber:viscosity or simul:viscosity should be defined");
    
    confine_space_ptr = sim->findSpace(confine_space);
    
    if ( confine_space_ptr )
        confine_space = confine_space_ptr->property()->name();

    if ( sim->prop->strict && confine != CONFINE_OFF )
    {
        if ( !confine_space_ptr )
            throw InvalidParameter(name()+":confine_space `"+confine_space+"' was not found");
    
        if ( confine_stiffness < 0 )
            throw InvalidParameter(name()+":confine_stiffness must be specified and >= 0");
    }

    if ( length <= 0 )
        throw InvalidParameter("fiber:length should be > 0");

    if ( min_length < 0 )
        throw InvalidParameter("fiber:min_length should be >= 0");

    if ( max_length < 0 )
        throw InvalidParameter("fiber:max_length should be >= 0");

#ifdef BACKWARD_COMPATIBILITY
    if ( total_polymer == 0 )
        total_polymer = INFINITY;
#endif
    if ( total_polymer <= 0 )
        throw InvalidParameter("fiber:total_polymer should be > 0 (you can specify 'inf')");

    if ( glue )
    {
        glue_set = const_cast<SingleSet*>(&sim->singles);
        if ( sim->prop->strict )
            glue_prop = static_cast<SingleProp*>(sim->properties.find_or_die("single", glue_single));
    }
    
    if ( field != "none" )
    {
        Property * fp = sim->properties.find("field", field);
        field_ptr = static_cast<Field*>(sim->fields.ObjectSet::first(fp));
    }
    
    if ( lattice && sim->prop->strict )
    {
        if ( lattice_unit <= 0 )
            throw InvalidParameter("fiber:lattice_unit (known as fiber:lattice[1]) must be specified and > 0");

        if ( lattice_flux_speed || lattice_binding_rate || lattice_unbinding_rate )
        {
            if ( field.empty() )
                throw InvalidParameter("fiber:lattice features require fiber:field to be specified");

            if ( field_ptr == 0 )
                throw InvalidParameter("fiber:field not found");
        }
    }
    lattice_unbinding_prob = 1 - exp( -lattice_unbinding_rate * sim->prop->time_step );
    
    if ( rigidity < 0 )
        throw InvalidParameter("fiber:rigidity must be specified and >= 0");
    
    if ( segmentation <= 0 )
        throw InvalidParameter("fiber:segmentation must be > 0");
 
#if ( 1 )
    if ( sim )
    {
        // Adjust the segmentation of all Fibers with this FiberProp:
        for ( Fiber* fib = sim->fibers.first(); fib; fib=fib->next() )
        {
            if ( fib->property() == this  &&  fib->segmentation() != segmentation )
            {
                fib->segmentation(segmentation);
                fib->update();
            }
        }
    }
#endif
    
    if ( steric && steric_radius <= 0 )
        throw InvalidParameter("fiber:steric[1] (radius) must be specified and > 0");
    
    if ( hydrodynamic_radius[0] <= 0 )
        throw InvalidParameter("fiber:hydrodynamic_radius[0] must be > 0");
    
    if ( hydrodynamic_radius[1] <= 0 )
        throw InvalidParameter("fiber:hydrodynamic_radius[1] must be > 0");
    
    if ( max_chewing_speed < 0 )
        throw InvalidParameter("fiber:max_chewing_speed must be >= 0");
    max_chewing_speed_dt = max_chewing_speed * sim->prop->time_step;
    
#if ( 0 )
    //print some information on the 'stiffness' of the matrix
    Fiber fib(this);
    fib.setStraight(Vector(0,0,0), Vector(1,0,0), 10, CENTER);
    
    fib.setDragCoefficient();
    real mob_dt = sim->prop->time_step * fib.nbPoints() / fib.dragCoefficient();
    
    real stiffness = 100;
    real coef1 = mob_dt * stiffness;

    MSG(5, "Numerical hardness (stiffness=%.1f): %7.2f\n", stiffness, coef1);

    real rod   = segmentation;
    real coef2 = mob_dt * rigidity / ( rod * rod * rod );
    
    MSG(5, "Numerical hardness (rigidity=%.1f): %7.2f\n", rigidity, coef2);
#endif
}


//------------------------------------------------------------------------------

void FiberProp::write_values(std::ostream & os) const
{
    write_value(os, "rigidity",            rigidity);
    write_value(os, "segmentation",        segmentation);
    write_value(os, "length",              length);
    write_value(os, "min_length",          min_length);
    write_value(os, "max_length",          max_length);
    write_value(os, "total_polymer",       total_polymer);
    write_value(os, "viscosity",           viscosity);
    write_value(os, "hydrodynamic_radius", hydrodynamic_radius, 2);
    write_value(os, "surface_effect",      surface_effect, cylinder_height);
#ifdef BACKWARD_COMPATIBILITY
    write_value(os, "squeeze",             squeeze, squeeze_force, squeeze_range);
#endif
    write_value(os, "binding_key",         binding_key);
    write_value(os, "lattice",             lattice, lattice_unit);
    write_value(os, "lattice_cut_fiber",   lattice_cut_fiber);
    write_value(os, "lattice_flux_speed",  lattice_flux_speed);
    write_value(os, "lattice_binding_rate", lattice_binding_rate);
    write_value(os, "lattice_unbinding_rate", lattice_unbinding_rate);
    write_value(os, "confine",             confine, confine_stiffness, confine_space);
    write_value(os, "steric",              steric, steric_radius, steric_range);
    write_value(os, "field",               field);
    write_value(os, "glue",                glue, glue_single);
    write_value(os, "colinear_force",      colinear_force);
    write_value(os, "max_chewing_speed",   max_chewing_speed);
    write_value(os, "activity",            activity);
    write_value(os, "display",             "("+display+")");
}


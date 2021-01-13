// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "simul_prop.h"
#include "assert_macro.h"
#include "space_prop.h"
#include "space_set.h"
#include "simul.h"
#include "exceptions.h"
#include "messages.h"
#include "glossary.h"
#include "property_list.h"
#include "filepath.h"
#include "random.h"

extern Random RNG;
extern bool functionKey[];

//------------------------------------------------------------------------------
void SimulProp::clear()
{
    viscosity         = 1;
#ifdef NEW_CYTOPLASMIC_FLOW
    flow.zero();
#endif
    time              = 0;
    time_step         = 0;
    kT                = 0.0042;
    tolerance         = 0.05;
    acceptable_rate   = 0.5;
    precondition      = 1;
    random_seed       = 0;
    steric            = 0;
    
    steric_stiffness_push[0] = 100;
    steric_stiffness_pull[0] = 100;
    steric_stiffness_push[1] = 100;
    steric_stiffness_pull[1] = 100;

    steric_max_range  = -1;
    binding_grid_step = -1;
    
    strict            = 0;
    verbose           = 0;

    config            = "config.cym";
    trajectory_file   = "objects.cmo";
    property_file     = "properties.cmo";
    clear_trajectory  = true;
    
    display           = "";
    display_fresh     = false;
    
    watch_dist_solids =     0;
    watch_dist_solids_t =   0;
    watch_mt_nb =           0;
    watch_mt_nb_t =         0;
    watch_AAcouples =       0;
    watch_AAcouples_name = "";
    watch_AAcouples_t =     0;
    handmonitor_pace      =     1;
    
}


void SimulProp::read(Glossary& glos)
{
    // a dimensionality can be specified to stop the program from running
    unsigned d = DIM;
    if ( glos.set(d, "dimensionality") || glos.set(d, "dim") )
    {
        if ( d != DIM )
            throw InvalidParameter("requested dimensionality is not fulfilled");
    }
    
    glos.set(viscosity,         "viscosity");
#ifdef NEW_CYTOPLASMIC_FLOW
    glos.set(flow,              "flow");
#endif
    glos.set(time,              "time");
    glos.set(time_step,         "time_step");
    glos.set(kT,                "kT");

    glos.set(tolerance,         "tolerance");
    glos.set(acceptable_rate,   "acceptable_rate");
    glos.set(precondition,      "precondition");
    
    glos.set(steric,                   "steric", KeyList<int>("off", 0, "on", 1));
    glos.set(steric_stiffness_push[0], "steric", 1);
    glos.set(steric_stiffness_pull[0], "steric", 2);
    glos.set(steric_stiffness_push, 2, "steric_stiffness_push");
    glos.set(steric_stiffness_pull, 2, "steric_stiffness_pull");
    glos.set(steric_max_range,         "steric_max_range");

    glos.set(binding_grid_step, "binding_grid_step");
    
    // these parameters are not written:
    glos.set(strict,            "strict");
    if ( glos.set(verbose, "verbose") )
        MSG.setVerbose(verbose);
    
    glos.set(functionKey, 17,   "function_key");
    
    // names of files and path:
    glos.set(config,            "config");
    glos.set(config,            ".cytosim");
    glos.set(config,            ".cym");
    
    glos.set(property_file,     "property_file");
    glos.set(property_file,     "property");
    
#ifdef BACKWARD_COMPATIBILITY
    glos.set(trajectory_file,   "object_file");
    bool a = false;
    if ( glos.set(a, "append_file") )
        clear_trajectory = !a;
#endif

    glos.set(trajectory_file,   "trajectory_file");
    glos.set(trajectory_file,   "trajectory");
    glos.set(trajectory_file,   ".cmo");
    glos.set(clear_trajectory,  "clear_trajectory");
    glos.set(random_seed,  "random_seed");
    
    if ( glos.set(display, "display") )
        display_fresh = true;
    
    
    // ADDED BY MAANU: SPINDLE WATCH VARIABLES: there should be some checks in here, but
    // its ok by now (forinstance if name of the property exists at all, etc)
    
    glos.set(watch_dist_solids, "watch_dist_solids", 0);
    glos.set(watch_dist_solids_t, "watch_dist_solids", 1);
    
    glos.set(watch_mt_nb, "watch_mt_nb", 0);
    glos.set(watch_mt_nb_t, "watch_mt_nb", 1);
    
    glos.set(watch_AAcouples, "watch_AAcouples",0);
    glos.set(watch_AAcouples_t, "watch_AAcouples",1);
    glos.set(watch_AAcouples_name, "watch_AAcouples",2);
    glos.set(handmonitor_pace, "handmonitor_pace");
}


/**
 If the Global parameters have changed, we update all derived parameters.
 This makes it possible to change the time-step in the middle of a config file.
 */
void SimulProp::complete(Simul const* sim)
{
    if ( !RNG.seeded() )
    {
        // initialize random number generator
        if ( random_seed )
            RNG.seed(random_seed);
        else
        {
            random_seed = RNG.seedTimer();
            MSG(5, "Cytosim: time-generated random seed 0x%lx\n", random_seed);
        }
    }
    
    if ( strict )
    {
        if ( time_step <= 0 )
            throw InvalidParameter("simul:time_step must be specified and > 0");
        
        if ( kT < 0 )
            throw InvalidParameter("simul:kT must be > 0");
        
        if ( kT == 0 && tolerance > 0.001 )
            throw InvalidParameter("if simul:kT==0, simul:tolerance must be set and <= 0.001");
    }
    /*
     If the Global parameters have changed, we update all derived parameters.
     To avoid an infinite recurence, the main SimulProp (*this) was
     not included in the PropertyList Simul::properties;
     */
    if ( sim )
        sim->properties.complete(sim);
}

//------------------------------------------------------------------------------

void SimulProp::write_values(std::ostream & os) const
{
    write_value(os, "time",            time);
    write_value(os, "time_step",       time_step);
    write_value(os, "kT",              kT);
    write_value(os, "viscosity",       viscosity);
#ifdef NEW_CYTOPLASMIC_FLOW
    write_value(os, "flow", flow);
#endif
    os << std::endl;
    write_value(os, "tolerance",       tolerance);
    write_value(os, "acceptable_rate", acceptable_rate);
    write_value(os, "precondition",    precondition);
    write_value(os, "random_seed",     random_seed);
    os << std::endl;
    write_value(os, "steric", steric, steric_stiffness_push[0], steric_stiffness_pull[0]);
    write_value(os, "steric_max_range",  steric_max_range);
    write_value(os, "binding_grid_step", binding_grid_step);
    write_value(os, "verbose", verbose);
    os << std::endl;

    write_value(os, "display", "("+display+")");
}


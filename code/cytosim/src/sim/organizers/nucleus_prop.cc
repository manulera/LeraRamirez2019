// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "property_list.h"
#include "nucleus_prop.h"
#include "sphere_prop.h"
#include "bundle_prop.h"
#include "property_list.h"
#include "simul_prop.h"
#include "glossary.h"
#include "smath.h"


//------------------------------------------------------------------------------
void NucleusProp::clear()
{
    sphere          = "";
    sphere_prop     = 0;
    fibers          = "";
    bundles         = "";
    bundle_prop     = 0;
    nb_fibers       = 0;
    nb_bundles      = 0;
    stiffness       = -1;
    focus           = MINUS_END;
    nucleation_rate = 0;
}


void NucleusProp::read(Glossary& glos)
{
    glos.set(sphere,          "sphere");
    glos.set(fibers,          "fibers");
    glos.set(nb_fibers,       "nb_fibers");
    glos.set(stiffness,       "stiffness");
    glos.set(focus,           "focus", KeyList<FiberEnd>("plus_end", PLUS_END, "minus_end", MINUS_END));
    glos.set(nucleation_rate, "nucleation_rate");
    glos.set(bundles,         "bundles");
    glos.set(nb_bundles,      "nb_bundles");
}


void NucleusProp::complete(Simul const* sim)
{
    if ( nb_fibers > 0 )
    {
        if ( fibers.empty() )
            throw InvalidParameter("nucleus:fibers must be specified");
        
        sim->properties.find_or_die("fiber", fibers);
    }
    
    if ( sphere.empty() )
        throw InvalidParameter("nucleus:sphere (a sphere) must be specified");
    
    sphere_prop = static_cast<SphereProp*>(sim->properties.find_or_die("sphere", sphere));
    
    if ( nb_bundles > 0 )
        bundle_prop = static_cast<BundleProp*>(sim->properties.find_or_die("bundle", bundles));
        
    if ( stiffness < 0 )
        throw InvalidParameter("nucleus:stiffness must be specified and >= 0");

    if ( nucleation_rate < 0 )
        throw InvalidParameter("nucleus:nucleation_rate must be positive");

    nucleation_rate_prob = 1 - exp( -nucleation_rate*sim->prop->time_step );
}



//------------------------------------------------------------------------------

void NucleusProp::write_values(std::ostream & os) const
{
    write_value(os, "sphere",          sphere);
    write_value(os, "fibers",          fibers);
    write_value(os, "nb_fibers",       nb_fibers);
    write_value(os, "stiffness",       stiffness);
    write_value(os, "focus",           focus);
    write_value(os, "bundles",         bundles);
    write_value(os, "nb_bundles",      nb_bundles);
    write_value(os, "nucleation_rate", nucleation_rate);
}


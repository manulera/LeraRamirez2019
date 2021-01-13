// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "aster_prop.h"
#include "property_list.h"
#include "simul_prop.h"
#include "solid_prop.h"
#include "fiber_prop.h"
#include "glossary.h"
#include "smath.h"

//------------------------------------------------------------------------------
void AsterProp::clear()
{
    solid           = "";
    fibers          = "";
    stiffness[0]    = -1;
    stiffness[1]    = -1;
    focus           = MINUS_END;
    nucleation_rate = 0;
}


void AsterProp::read(Glossary& glos)
{
    glos.set(solid,           "solid");
    glos.set(fibers,          "fibers");
    glos.set(stiffness, 2,    "stiffness");
    glos.set(focus,           "focus", KeyList<FiberEnd>("plus_end", PLUS_END, "minus_end", MINUS_END));
    glos.set(nucleation_rate, "nucleation_rate");
}


void AsterProp::complete(Simul const* sim)
{
    if ( fibers.empty() )
        throw InvalidParameter("aster:fibers must be specified");
    
    sim->properties.find_or_die("fiber", fibers);
    
    if ( stiffness[0] < 0 )
        throw InvalidParameter("aster:stiffness[0] must be specified and >= 0");
    
    if ( stiffness[1] < 0 )
        throw InvalidParameter("aster:stiffness[1] must be specified and >= 0");

    if ( nucleation_rate < 0 )
        throw InvalidParameter("aster:nucleation_rate must be >= 0");
    
    nucleation_rate_prob = 1 - exp( -nucleation_rate*sim->prop->time_step );
}



//------------------------------------------------------------------------------

void AsterProp::write_values(std::ostream & os) const
{
    write_value(os, "solid",           solid);
    write_value(os, "fibers",          fibers);
    write_value(os, "stiffness",       stiffness[0], stiffness[1]);
    write_value(os, "focus",           focus);
    write_value(os, "nucleation_rate", nucleation_rate);
}


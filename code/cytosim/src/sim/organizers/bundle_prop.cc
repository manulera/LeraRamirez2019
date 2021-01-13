// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "bundle_prop.h"
#include "property_list.h"
#include "glossary.h"
#include "smath.h"
#include "simul.h"


//------------------------------------------------------------------------------
void BundleProp::clear()
{
    fibers     = "";
    nb_fibers  = 0;
    stiffness  = -1;
    overlap    = -1;
    focus      = MINUS_END;
    nucleate   = true;
}


void BundleProp::read(Glossary& glos)
{
    glos.set(fibers,       "fibers");
    glos.set(nb_fibers,    "nb_fibers");
    glos.set(stiffness,    "stiffness");
    glos.set(overlap,      "overlap");
    glos.set(focus,        "focus", KeyList<FiberEnd>("plus_end", PLUS_END, "minus_end", MINUS_END));
    glos.set(nucleate,     "nucleate");
}


void BundleProp::complete(Simul const* sim)
{
    if ( fibers.empty() )
        throw InvalidParameter("bundle:fibers must be specified");
    
    sim->properties.find_or_die("fiber", fibers);

    if ( nb_fibers <= 0 )
        throw InvalidParameter("bundle:nb_fibers must be specified and > 0");
    
    if ( overlap < 0 )
        throw InvalidParameter("bundle:overlap must be specified and >= 0");
    
    if ( stiffness < 0 )
        throw InvalidParameter("bundle:stiffness must be specified and >= 0");
}



//------------------------------------------------------------------------------

void BundleProp::write_values(std::ostream & os) const
{
    write_value(os, "fibers",    fibers);
    write_value(os, "nb_fibers", nb_fibers);
    write_value(os, "stiffness", stiffness);
    write_value(os, "overlap",   overlap);
    write_value(os, "focus",     focus);
    write_value(os, "nucleate",  nucleate);
}


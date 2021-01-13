// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "glapp_prop.h"
#include "glossary.h"


//------------------------------------------------------------------------------
void glAppProp::clear()
{
    full_screen    = 0;
    buffered       = 0;

    show_message   = 0;
    message        = "Please, visit www.cytosim.org";
}


void glAppProp::read(Glossary& glos)
{
    glos.set(full_screen,     "full_screen");
    glos.set(show_message,    "show_message");
    glos.set(buffered,        "buffered");
}
        
//------------------------------------------------------------------------------

void glAppProp::write_values(std::ostream & os) const
{
    write_value(os, "full_screen",    full_screen);
    write_value(os, "buffered",       buffered);
}



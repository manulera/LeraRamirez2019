// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "display_prop.h"
#include "glossary.h"
#include "dim.h"

//------------------------------------------------------------------------------
void DisplayProp::clear()
{
    style          = 1;
    tile           = 0;
    fold           = 1;
    
    couple_flip    = 1;
    couple_select  = 7;
    single_select  = 3;
    
    point_value    = 0;
    point_size     = 5;
    link_size      = 4;
    line_width     = 2;
}

//------------------------------------------------------------------------------
void DisplayProp::read(Glossary& glos)
{
    glos.set(style,         "style");
    glos.set(tile,          "tile");
    glos.set(fold,          "fold");
    glos.set(fold,          "tile", 1);

    glos.set(tile,          "periodic");
    glos.set(fold,          "periodic", 1);

    glos.set(tile,          "tiled");
    glos.set(fold,          "tiled", 1);
    
    glos.set(couple_flip,   "couple_flip");
    glos.set(couple_select, "couple_select");
    glos.set(single_select, "single_select");
    
    glos.set(point_value,   "point_value");
    glos.set(point_size,    "point_size");
    // unless specified, `link_size` will be equal to `line_width`:
    if ( glos.set(line_width, "line_width") )
        link_size = line_width;
    glos.set(link_size,     "link_size");
}


//------------------------------------------------------------------------------

void DisplayProp::write_values(std::ostream & os) const
{
    write_value(os, "style",         style);
    write_value(os, "tile",          tile, fold);
    write_value(os, "couple_flip",   couple_flip);
    write_value(os, "couple_select", couple_select);
    write_value(os, "single_select", single_select);
    write_value(os, "point_value",   point_value);
    write_value(os, "point_size",    point_size);
    write_value(os, "link_size",     link_size);
    write_value(os, "line_width",    line_width);
}



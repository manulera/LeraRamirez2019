// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "fiber_disp.h"
#include "glossary.h"
#include "filament.h"
#include "sim.h"


//------------------------------------------------------------------------------

void FiberDisp::clear()
{
    style            = 0;
    visible          = 1;
    color            = 0xFFFFFFFF;
    colorT           = 0xFFFFFF00;
    inner_color      = 0x777777FF;
    coloring         = 0;
    color_left       = 0xCCCCCCFF;
    color_right      = 0x00CC00FF;
    
    line_style       = 1;
    line_width       = 2;
    line_caps        = 1;
    
    point_style      = 0;
    point_size       = 5;
    point_interval   = 1;

    end_style[0]     = 0;
    end_style[1]     = 0;

    end_size[0]      = 6;
    end_size[1]      = 6;
    
    end_length[0]    = 0;
    end_length[1]    = 0;
    
    end_color[0]     = 0xFFFFFFFF;  // white
    end_color[1]     = 0x00FF00FF;  // green
    end_color[2]     = 0xFFFF00FF;  // yellow
    end_color[3]     = 0xFF7538FF;  // orange
    end_color[4]     = 0xFF0000FF;  // red
    
    lattice_style    = 0;
    lattice_scale    = 1;
    
    label_style      = 0;    
    speckle_style    = 0;
    speckle_interval = 1;
    
    exclude          = 0;
    right_direction.set(1,0,0);
    
    mask             = 0;
    mask_bitfield    = 0;

    tension          = 10;
    forces           = 0;
    forces_color     = 0xFF0000FF;
    
    explode          = 0;
    explode_range    = 0;
    show_average     = 0;
}


//------------------------------------------------------------------------------
void FiberDisp::read(Glossary& glos)
{
    glos.set(style,            "style", KeyList<int>("line", 0, "actin", 1, "microtubule", 2));
    glos.set(visible,          "visible");
    glos.set(color,            "color");
    glos.set(colorT,           "color", 1);
    glos.set(inner_color,      "color", 2);
    glos.set(inner_color,      "inner_color");
    glos.set(coloring,         "coloring");
    glos.set(color_left,       "color_left");
    glos.set(color_right,      "color_right");
    
#ifdef BACKWARD_COMPATIBILITY
    glos.set(line_width,       "line");
    glos.set(line_style,       "line", 1);
#endif
    glos.set(line_width,       "lines");
    glos.set(line_style,       "lines", KeyList<int>("none", 0, "line", 1, "force", 2, "curvature", 3, "angle", 4), 1);
    glos.set(line_caps,        "lines", 2);
    glos.set(line_style,       "line_style", KeyList<int>("none", 0, "line", 1, "force", 2, "curvature", 3, "angle", 4));
    glos.set(line_width,       "line_width");
    glos.set(line_width,       "width");
    glos.set(line_caps,        "line_caps");

#ifdef BACKWARD_COMPATIBILITY
    glos.set(point_size,       "point");
    glos.set(point_style,      "point", 1);
#endif
    glos.set(point_size,       "points");
    glos.set(point_style,      "points", 1);
    glos.set(point_interval,   "points", 2);
    
    glos.set(point_interval,   "point_interval");

    if ( point_interval <= 0 )
        point_interval = 1;

    glos.set(point_style,      "point_style");
    glos.set(point_size,       "point_size");
    glos.set(point_size,       "size");

    if ( glos.set(end_size[0], "plus_end") )
        end_style[0] = 2;
    glos.set(end_style[0],     "plus_end", 1);
    
    if ( glos.set(end_size[1], "minus_end") )
        end_style[1] = 3;
    glos.set(end_style[1],     "minus_end", 1);
    
    glos.set(end_style, 2,     "end_style");
    glos.set(end_size,  2,     "end_size");
    glos.set(end_length, 2,    "end_length");
    glos.set(end_color, 5,     "end_color");
    
#ifdef BACKWARD_COMPATIBILITY
    glos.set(end_length, 2,    "end_section");
    glos.set(lattice_style,    "show_lattice");
    glos.set(lattice_scale,    "lattice_max");
#endif
    
    glos.set(lattice_scale,    "lattice_scale");
    glos.set(lattice_style,    "lattice");
    glos.set(lattice_scale,    "lattice", 1);

    glos.set(label_style,      "labels");
    glos.set(label_style,      "label_style");
    
    glos.set(speckle_style,    "speckles", KeyList<int>("off", 0, "random", 1, "regular", 2));
    glos.set(speckle_interval, "speckles", 1);

    glos.set(speckle_style,    "speckle_style");
    glos.set(speckle_interval, "interval");
    glos.set(speckle_style,    "speckles");
    glos.set(speckle_interval, "speckles", 1);

    if ( speckle_interval <= 0 )
        speckle_interval = 1;

    glos.set(exclude,          "exclude");
    glos.set(right_direction,  "exclude", 1);
    glos.set(right_direction,  "right_direction");
    
    if ( glos.set(mask, "mask") )
        mask_bitfield = RNG.pint_bits(mask);
    glos.set(mask_bitfield, "mask", 1);

    glos.set(tension,          "tension");
    glos.set(forces,           "forces");
    glos.set(forces_color,     "forces", 1);
    
    glos.set(explode,          "explode");
    glos.set(explode_range,    "explode", 1);
    
#ifdef BACKWARD_COMPATIBILITY
    if ( glos.set(explode_range, "display_shift") )
        explode = 1;
#endif
    
    glos.set(show_average,     "show_average");
}


//------------------------------------------------------------------------------

void FiberDisp::write_values(std::ostream & os) const
{
    write_value(os, "style",        style);
    write_value(os, "visible",      visible);
    if ( colorT.visible() )
        write_value(os, "color",    color, colorT);
    else
        write_value(os, "color",    color);
    write_value(os, "inner_color",  inner_color);
    write_value(os, "coloring",     coloring);
    write_value(os, "color_left",   color_left);
    write_value(os, "color_right",  color_right);
    
    write_value(os, "points",       point_size, point_style, point_interval);
    write_value(os, "lines",        line_width, line_style, line_caps);
    write_value(os, "plus_end",     end_size[0], end_style[0]);
    write_value(os, "minus_end",    end_size[1], end_style[1]);
    write_value(os, "end_length",   end_length, 2);
    write_value(os, "end_color",    end_color, 5);
 
    write_value(os, "lattice",      lattice_style, lattice_scale);
    write_value(os, "labels",       label_style);
    write_value(os, "speckles",     speckle_style, speckle_interval);
    write_value(os, "exclude",      exclude, right_direction);
    write_value(os, "mask",         mask, mask_bitfield);

    write_value(os, "tension",      tension);
    write_value(os, "forces",       forces);
    write_value(os, "explode",      explode, explode_range);
    write_value(os, "show_average", show_average);
}


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef DISPLAY_PROP_H
#define DISPLAY_PROP_H

#include "real.h"
#include "property.h"
#include "gle_color.h"
#include "inventoried.h"
class View;


/// Property for Play
class DisplayProp : public Property
{

public:
    
    /**
     @defgroup DisplayPar Display parameters: World
     @ingroup DisplayParameters
     @{
     */
    
    /// style of display { 1, 2, 3 }
    /**
     3 styles are implemented:
     - style 1 used OpenGL lines and points. It is suitable for 2D work.
     - style 2 is a faster display, also suitable for 2D.
     - style 3 draw real tubes and uses OpenGL lighting for rendering. It is nice for 3D.
     .
     */
    unsigned       style;

    /// if true, repeat the display for periodic boundary conditions
    int            tile;
    
    /// if true, translate objects to place them in the root cell for periodic boundary conditions
    int            fold;
    
    /// default diameter of points
    real           point_size;
    
    /// default width of hookean links
    real           link_size;

    /// default width of lines
    real           line_width;
    
    /// if set > 0, this defines the unit size used for `point_size` and `line_width` 
    /**
     Set this parameter to specify the fiber radius and point size in real units.

     `point_size` and `line_width` are normally set in pixels, but if `point_value` is set,
     then the specifications are understood in multiples of `point_value`,
     which itself is given simulation units (aka. real distance).
     
     For example, if you set `line_width=5` and `point_value=0.005`,
     the fibers will be displayed with a diameter of 0.025.
     
     <em> default = 0 </em>
     */
    real           point_value;
    
    /// if `true`, unattached Couples are display randomly with one or the other Hand (default=false)
    unsigned       couple_flip;
    
    /// selection bitfield for Couples
    unsigned       couple_select;
    
    /// selection bitfield for Singles
    unsigned       single_select;
    
    /// @}

public:
    
    /// constructor
    DisplayProp(const std::string& n) : Property(n) { clear(); }
    
    /// destructor
    ~DisplayProp() { }
    
    /// identifies the property
    std::string category() const { return "simul:display"; }
        
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// return a carbon copy of object
    Property* clone() const { return new DisplayProp(*this); }
    
    /// write all values
    void write_values(std::ostream &) const;

};


#endif



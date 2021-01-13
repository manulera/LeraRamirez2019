// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef FIBER_DISP_H
#define FIBER_DISP_H

#include "real.h"
#include "assert_macro.h"
#include "gle_color.h"
#include "property.h"
#include "vector.h"

class Glossary;

/// Display parameters for a class of Fiber
/**
 Holds the display attributes for a certain class of Fiber.
 
 There is one FiberDisp for each FiberProp.
 */
class FiberDisp : public Property
{
public:

    /// possible values for fiber:coloring
    enum ColoringModes {
        COLORING_OFF,
        COLORING_RANDOM,
        COLORING_DIRECTION,
        COLORING_MARK,
        COLORING_FLECK,
        COLORING_AGE
    };
    
public:
    
    /**
     @defgroup FiberDispPar Display parameters: Fibers
     @ingroup DisplayParameters
     @{
     */
    
    /// general rendering style
    /**
     Possible values of `style`:
     - 0 or 'line'        : line or cylinders for style=3 (this is the default),
     - 1 or 'actin'       : actin-like rendering using beads for monomers,
     - 2 or 'microtubule' : microtubule-like rendering using beads for monomers.
     .
     */
    int          style;
    
    /// visibility flag
    int          visible;
    
    /// color of fiber
    gle_color    color;
    
    /// second color (set as color[1])
    gle_color    colorT;

    /// color of inner surfaces of cylinder in 3D display
    gle_color    inner_color;

    /// if true, vary the colors used to display the fibers
    /**
     Values for `coloring`:
     - 0 : no coloring,
     - 1 : color fibers according to ID-number,
     - 2 : color fibers depending on direction, relative to `right_direction`,
     - 3 : color fibers depending on the mark,
     - 4 : color fibers by connectivity.
     .
     */
    int          coloring;
    
    /// color for fibers parallel to `right_direction` (default is green)
    gle_color    color_right;

    /// color for fibers antiparallel  to `right_direction` (default is white)
    gle_color    color_left;

    /// width of lines (also known as `line[0]` or `width`)
    real         line_width;

    /// style for lines (also known as `line[1]`)
    /**
     Possible values of `line_style`:
     - 0 : hide,
     - 1 : plain lines,
     - 2 : color-based display  of longitudinal tensions,
     - 3 : color set by the local curvature,
     - 4 : color set by the angular orientation of each segment relative to the X-axis
     .
     */
    int          line_style;
    
    /// if true, close the end of the fiber (valid only for style==3)
    /**
     Possible values of `line_caps`:
     - 0: leave fibers open (unfinished),
     - 1: use a disc to make a flat end,
     - 2: use a hemisphere to make a round end.
     This is only valid for style==3
     */
    int          line_caps;
    
    /// diameter of points (also known as `point[0]` or `size`)
    /**
     `point_size` and `line_width` are normally set in pixels, 
     but if `display`:point_value is set, their value is understood 
     in multiples of `point_value`, which itself is a distance.
     
     For example, if you set line_width=2.5 and point_value=0.01,
     the fibers will be displayed with a diameter of 0.025.
     */
    real         point_size;
    
    /// style for display of points (also known as `point[1]`)
    /**
     Possible point_style:
     - 0 : hide,
     - 1 : show model points,
     - 2 : show arrowheads along fiber, separated by `point_interval`,
     - 3 : show middle point of each fiber
     - 5 : actin-like rendering
     - 7 : microtubule-like rendering
     .
     */
    int          point_style;
    
    /// distance between arrows for `point_style=2` (also known as `point[2]`)
    real         point_interval;
    

    /// style of fiber tips for { PLUS_END, MINUS_END }
    /**
     end_style[0] determines the style of the PLUS_END,
     and end_style[1] the style of the MINUS_END.
     
     Possible end_style:
     - 0 : hide,
     - 1 : display a disc/sphere,
     - 2 : display a cone,
     - 3 : display a disc,
     - 4 : draw arrowhead,
     - 5 : draw arrowhead in the inverted direction (for actin)
     .
     */
    int          end_style[2];
    
    /// size of fiber tips for { PLUS_END, MINUS_END }
    /**
     You can also specify:
     @code
     plus_end  = SIZE, STYLE
     minus_end = SIZE, STYLE
     @endcode
     */
    real         end_size[2];
    
    /// length of a section displayed near the fiber tips
    /**
     if `end_length[0] > 0`, a section near the PLUS_END is drawn with the color of the PLUS_END.
     if `end_length[1] > 0`, a section near the MINUS_END is drawn with the color of the MINUS_END.
     */
    real         end_length[2];
    
    /// colors of the different FiberTip states
    /**
     This determines the set of color that are used to display the fiber tips,
     according to their assembly state (Fiber::dynamicState):
     - static ends (dynamic-state 0) use end_color[0],
     - growing end (dynamic-state 1), use end_color[1],
     - shrinking end (dynamic-state 4), use end_color[4]
     .
     */
    gle_color    end_color[5];
    
    
    /// if true, specify the style for displaying lattice content (also known as `lattice[0]`)
    int          lattice_style;
    
    /// defines the range of colors when displaying the lattice (also known as `lattice[1]`)
    real         lattice_scale;

    
    /// style of labels
    /**
     Possible `label_style`:
     - 0 : hide,
     - 1 : index of model points
     - 2 : abscissa along fiber
     .
     */
    int          label_style;
    
    
    /// style for speckle display (also know as `speckles`)
    /**
     Possible values for `speckle_style`:
     - 0 : hide,
     - 1 : random speckles, separated on average by `interval`,
     - 2 : regular speckes, separated by `interval`.
     .
     */
    int          speckle_style;
    
    /// distance between speckles (also known as `speckles[1]`)
    real         speckle_interval;
    
    
    /// a bit-field to hide certain categories of fibers
    /**
     Possible values for `exclude`:
     - 0 : all fibers are displayed,
     - 1 : show only right-pointing fibers,
     - 2 : show only left-pointing fibers,
     - 4 : show only counter-clockwise fibers,
     - 8 : show only clockwise fibers.
     .
     
     You may also address each bit directly, knowning that:
     - 1st bit true: hide left-pointing fibers
     - 2nd bit true: hide right-pointing fibers
     - 3rd bit true: hide clockwise fibers
     - 4th bit true: hide counter-clockwise fibers
     .
     */
    int          exclude;
    
    /// the direction used for hiding left- or right-pointing fibers, etc. (known as `exclude[1]`)
    Vector       right_direction;
    
    
    
    /// number of bits equal to `1` in the mask_bitfield
    /**
     This parameter can be used to hide a fraction of the fiber.
     Each fiber will be visible with a probability `1/2^mask`.
     `mask_bitfield' is set randomly with `mask` bits set to 1, 
     When the parameter is read.
     */
    unsigned int mask;
    
    /// selection bitfield used to hide some fibers (known as `mask[1]`)
    /**
     `mask_bitfield` is a 32-bitfield that is compared to the signature of the object,
     itself a random bitfield. The Object is hidden if the result is non-zero.
     So if `mask_bitfield` has many 1s, fewer filaments will be visible.
     Note that `mask_bitfield' is set randomly if `mask` is given.
     */
    unsigned int mask_bitfield;
    
    
    /// conversion coefficient from tension to color, for `line_style==2`
    /**
     The values of `tension` determines how longitudinal tensions are displayed:
     - tension < 0 : compressive forces are highlighted,
     - tension > 0 : tensile forces are highlighted.
     .
     A longitudinal tension equal to `tension` will be displayed with a blue tint,
     while `3*tension` will be displayed red.
     Lower `tension` values will yield brighter colors for the same force.
     */
    real         tension;
    
    /// ( if > 0 ) display the net forces FP acting on model points
    /**
     The force is displayed as segments of length `forces*PF`.
     A color can be specified as forces[1]
     */
    real         forces;

    
    /// the 'explosion' effect shift the fibers in space
    /**
     This can be useful to visualize dense regions,
     but is only implemented for style=2
     */
    int          explode;
    
    /// amount of lateral shift to separate fibers when display is exploded (known as `explode[1]`)
    real         explode_range;
    
    
    /// if true, display the average fiber
    /**
     The 'average fiber' is calculated from the centroid of the fiber tips,
     and the centroid of the polymer mass.
     It is useful to evaluate the amount of order in the network.
     */
    int          show_average;

    /// @}
    
    /// this color is specified as forces[1]
    gle_color    forces_color;

public:
    
    /// constructor
    FiberDisp(const std::string& n) : Property(n)  { clear(); }
    
    /// destructor
    ~FiberDisp() { }
    
    /// identifies the property
    std::string category() const { return "fiber:display"; }

    /// clear to default values
    void clear();
    
    /// set from glossary
    void read(Glossary&);
    
    /// return a carbon copy of object
    Property* clone() const { return new FiberDisp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};


#endif


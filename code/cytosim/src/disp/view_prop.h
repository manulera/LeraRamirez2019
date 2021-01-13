// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef VIEW_PROP_H
#define VIEW_PROP_H

#include "real.h"
#include "vector3.h"
#include "quaternion.h"
#include "property.h"
#include "gle_color.h"

///properties needed to define a view
class ViewProp : public Property
{
protected:
    
    /// text displayed near top right corner of window
    mutable std::string message_top;
    
    /// text displayed near bottom left corner of window
    mutable std::string message_bot;
    
    /// number of OpenGL clipping planes
    const static unsigned int NB_CLIP_PLANES = 4;

public:
    
    /**
     @defgroup ViewPar Display Parameters: View
     @ingroup DisplayParameters
     @{
     */
    
    /// color of background
    gle_color        back_color;
    
    /// color used to highlight objects
    gle_color        front_color;
    
    /// color used to paint the inner surface of objects
    gle_color        inner_color;

    /// string at start of `message` (if `none` is specified, no message is shown)
    std::string      label;
    
    /// flag to enable OpenGL depth buffer (default=1)
    /**
     This is useful for 3D rendering.
     http://en.wikipedia.org/wiki/Z-buffering
     */
    int              depth_test;
    
    /// flag to enable perspective view in 3D
    /**
     By default, cytosim uses a orthographic projection to view the 3D space,
     but it will use a 3D perspective if 'perspective==true'.
     This is only meaningful in 3D mode.
     */
    int              perspective;
    
    /// flag to enable native device resolution on mac osx
    /**
     This works only if you use Renaud Blanch's modified GLUT
     http://iihm.imag.fr/blanch/software/glut-macosx
     */
    int              retina;
    
    /// flag to enable OpenGL stencil buffer (default=0)
    int              stencil;
    
    /// flag to perform depth-clamp (default=false)
    /** http://www.opengl.org/registry/specs/NV/depth_clamp.txt */
    int              depth_clamp;
    
    /// if > 0, enables OpenGL full scene anti-aliasing (default=0)
    /**
     This defines the number of samples used to build an image.
     Higher values result in nicer (but slower) display.
     http://en.wikipedia.org/wiki/Multisample_anti-aliasing
     Many graphic cards only support 8 samples max, so try 4 or 8.
     */
    int              multisample;

    /// zoom factor = ratio between visible area and `view_size`
    real             zoom;
    
    /// size of area visible in the window, in sim-units (default=10)
    real             view_size;

    /// enables the display area to be set from the size of the simulation space
    /**
     If ( `auto_scale` > 0 ), `view_size` is set automatically to match the simulation space.
     This is on by default.
     */
    unsigned int     auto_scale;
    
    /// the point that is in the center of the window in real-world coordinates
    Vector3          focus;

    /// additional translation used by autoTrack
    Vector3          focus_shift;
    
    /// orientation of display
    Quaternion<real> rotation;
    
    
    /// modifies the display to show only the front, the back or a slice of the world
    /**
     possible values are:
     - `off`    (0)
     - `front`  (1)
     - `back`   (2)
     - `slice`  (3)
     .
     */
    unsigned int     slice;
    
    /// enables auto_translation, auto_zoom or auto_rotation
    /**
     if ( `traveling` > 0 ), this sets the interval of time in milli-seconds
     at which the model-view transformation will be updated, by applying
     `auto_translation`, `auto_zoom` and `auto_rotation`.
     */
    unsigned int     traveling;
    
    /// translation speed of display (known as `traveling[1]`)
    /**
     This is a speed in distance / wall-time
     */
    Vector3          auto_translation;
    
    /// rotational speed of display (known as `traveling[2]`)
    /**
     This is a speed in quaternion / wall-time
     */
    Quaternion<real> auto_rotation;
    
    /// zooming speed of display (known as `traveling[3]`)
    /**
     This is a speed in zoom-unit / wall-time
     - auto_zoom > 0 : zoom closer
     - auto_zoom < 0 : zoom away
     .
     */
    real             auto_zoom;
    
    /// automatically adjust view to keep fibers in window
    /**
     Possible values:
     - 0 : off
     - 1 : translate to track the center of gravity of the cloud of fiber-points
     - 2 : rotate to align the principal direction of the fiber
     - 3 : translate and rotate ( 1 and 2 are combined )
     - 4 : rotate to align two principal directions
     - 5 : translate and rotate ( 1 and 4 are combined )
     .
     The translation defined by focus is applied after this adjustment.
     */
    unsigned int     track_fibers;
    
    /// position of window on screen (top-left corner, in pixels)
    int              window_position[2];
    
    /// desired size of window in pixels (also known as `size`)
    int              window_size[2];
    
    /// size of scale-bar in sim-world units (set as `scale_bar[0]`)
    real             scale_bar_size;
    
    /// color of scale-bar (set as `scale_bar[1]`)
    gle_color        scale_bar_color;
    
    /// display flag for scale-bar (set as `scale_bar[2]`)
    unsigned int     scale_bar_mode;

    /// display flag for displaying X-Y-Z axes
    unsigned int     show_axes;
    
    /// length of axes (set a `show_axes[1]`, default=1)
    real             axes_size;

    /// on/off flags for clipping (defined as `clip_plane?`)
    /**
     Up to 4 clipping planes can be defined: clip_plane0 to clip_plane3
     
     Syntax:
     @code
        clip_plane? = BOOL, VECTOR, REAL
     @endcode
     The Boolean enables the clipping plane.
     The plane is specified by a normal vector `n` (VECTOR) and a scalar `a` (REAL).
     The visible half-space corresponds to <em> n.pos + a > 0 </em>
     
     Example:
     To define a slice perpendicular to the X-axis of width 2: 
     @code
     set simul:display *
     {
        clip_plane1 = 1,  1 0 0, 1
        clip_plane2 = 1, -1 0 0, 1
     }
     @endcode
     */
    unsigned int     clip_plane_mode[NB_CLIP_PLANES];

    /// direction perpendicular to clipping plane (defined as `clip_plane?[1]`)
    Vector3          clip_plane_vector[NB_CLIP_PLANES];
    
    /// scalar offset defining the equation of the clipping plane (defined as `clip_plane?[2]`)
    real             clip_plane_scalar[NB_CLIP_PLANES];

    
    /// characteristics of OpenGL fog (also known as `fog[0]`)
    int              fog_type;
    
    /// density of fog (also known as `fog[1]`)
    real             fog_param;
    
    /// color of fog (also known as `fog[2]`)
    gle_color        fog_color;
 
    /// @}
    
public:
   
    /// constructor
    ViewProp(const std::string& n) : Property(n)  { clear(); }
    
    /// destructor
    ~ViewProp()  { }
    
    /// identifies the property
    std::string category() const { return "view"; }
     
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// return a carbon copy of object
    Property* clone() const { return new ViewProp(*this); }

    /// write all values
    void write_values(std::ostream &) const;

};

#endif

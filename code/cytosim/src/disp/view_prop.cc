// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "view_prop.h"
#include "glossary.h"
#include "smath.h"

//------------------------------------------------------------------------------
void ViewProp::clear()
{
    back_color   = 0x000000FF;
    front_color  = 0xFFFFFFFF;
    inner_color  = 0x444444FF;
    label        = "Cytosim";
    message_top  = "";
    message_bot  = "";
    depth_test   = 1;
    stencil      = 0;
    multisample  = 0;
    depth_clamp  = 0;
    retina       = 0;
    perspective  = 0;

    zoom         = 1;
    view_size    = 10;
    auto_scale   = 1;
    
    focus.zero();
    focus_shift.zero();
    rotation.set(1,0,0,0);
    
    slice        = 0;
    traveling    = 0;
    auto_zoom    = 0;
    auto_rotation.set(1,0,0,0);
    auto_translation.zero();
    
    window_size[0]     = 768;
    window_size[1]     = 768;
    window_position[0] = 0;
    window_position[1] = 50;
    
    scale_bar_size     = 10;
    scale_bar_color    = 0xFFFF88AA;
    scale_bar_mode     = 0;
    show_axes          = 0;
    axes_size          = 1;
    
    for ( unsigned k = 0; k < NB_CLIP_PLANES; ++k )
    {
        clip_plane_mode[k] = 0;
        clip_plane_vector[k].set(1,0,0);
        clip_plane_scalar[k] = 0;
    }
    
    track_fibers       = 0;
    
    fog_type           = 0;
    fog_param          = 1;
    fog_color          = 0x000000FF;
}


void ViewProp::read(Glossary& glos)
{
    glos.set(back_color,  "background_color");
    if ( glos.set(back_color,  "back_color") )
    {
        fog_color = back_color;
        front_color = back_color.inverted();
    }
    glos.set(front_color, "front_color");
    glos.set(inner_color, "inner_color");

    glos.set(label,       "label");
    glos.set(depth_test,  "depth_test");
    glos.set(stencil,     "stencil");
    glos.set(multisample, "samples");
    glos.set(multisample, "multisample");
    glos.set(depth_clamp, "depth_clamp");
    glos.set(retina,      "retina");
    glos.set(perspective, "perspective");

// BACKWARD_COMPATIBILITY
    glos.set(multisample, "gl_samples");
    glos.set(multisample, "multisample");
    
    glos.set(zoom,        "zoom");
    if ( glos.set(view_size, "view_size") )
        auto_scale = 0;
    glos.set(auto_scale,  "auto_scale");
    glos.set(focus,       "focus");
    glos.set(rotation,    "rotation");
    
    // normalize quaternion:
    if ( rotation.norm() > 0.001 )
        rotation.normalize();
    else
        rotation.set(1,0,0,0);

    glos.set(slice, "slice",
             KeyList<unsigned>("off",   0,
                               "front", 1,
                               "back",  2,
                               "slice", 3));
    
    glos.set(traveling,        "traveling");
    glos.set(auto_translation, "traveling", 1);
    glos.set(auto_rotation,    "traveling", 2);
    glos.set(auto_zoom,        "traveling", 3);
    
    glos.set(track_fibers,     "track_fibers");

    // normalize quaternion:
    if ( auto_rotation.norm() > 0.001 )
        auto_rotation.normalize();
    else
        auto_rotation.set(1,0,0,0);
    
    glos.set(window_position, 2,  "window_position");
    
    // A square window is made if only one value is given.
    if ( glos.set(window_size, 2, "window_size") == 1 )
        window_size[1] = window_size[0];

    // A square window is made if only one value is given.
    if ( glos.set(window_size, 2, "image_size") == 1 )
        window_size[1] = window_size[0];

    // 'size' is an alias to set the size of the window.
    if ( glos.set(window_size, 2, "size") == 1 )
        window_size[1] = window_size[0];
    /*
     window_size can be changed here, but we cannot resize
     the window, since we do not have access to GLUT 
     */
    //std::clog << this << " window " << window_size[0] << "x" << window_size[1] << std::endl;
    
    glos.set(scale_bar_size,  "scale_bar");
    glos.set(scale_bar_color, "scale_bar", 1);
    glos.set(scale_bar_mode,  "scale_bar", 2);
    glos.set(scale_bar_color, "scale_bar_color");

    glos.set(show_axes,       "show_axes");
    glos.set(axes_size,       "show_axes", 1);
    glos.set(axes_size,       "axes_size");

    for ( unsigned k = 0; k < NB_CLIP_PLANES; ++k )
    {
        std::string var = "clip_plane" + sMath::repr(k);
        glos.set(clip_plane_mode[k],   var);
        glos.set(clip_plane_vector[k], var, 1);
        glos.set(clip_plane_scalar[k], var, 2);
    }
    
    KeyList<GLint> fogValues("off",          0,
                             "linear",       1,
                             "exponential",  2,
                             "exponential2", 3);

    glos.set(fog_type,     "fog_type", fogValues);
    glos.set(fog_param,    "fog_param");
    glos.set(fog_color,    "fog_color");
 
    glos.set(fog_type,     "fog", fogValues);
    glos.set(fog_param,    "fog", 1);
    glos.set(fog_color,    "fog", 2);
}

//------------------------------------------------------------------------------

void ViewProp::write_values(std::ostream & os) const
{
    write_value(os, "back_color",    back_color);
    write_value(os, "front_color",   front_color);
    write_value(os, "inner_color",   inner_color);
    write_value(os, "label",         label);
    write_value(os, "depth_test",    depth_test);
    write_value(os, "stencil",       stencil);
    write_value(os, "samples",       multisample);
    write_value(os, "depth_clamp",   depth_clamp);
    write_value(os, "perspective",   perspective);
    write_value(os, "retina",        retina);
    write_value(os, "zoom",          zoom);
    write_value(os, "focus",         focus+focus_shift);
    write_value(os, "rotation",      rotation);
    write_value(os, "slice",         slice);
    write_value(os, "traveling",     traveling, auto_translation, auto_rotation, auto_zoom);
    write_value(os, "window_size",   window_size, 2);
    //write_value(os, "window_position", window_position, 2);

    write_value(os, "view_size",     view_size);
    write_value(os, "scale_bar",     scale_bar_size, scale_bar_color, scale_bar_mode);
    write_value(os, "show_axes",     show_axes, axes_size);

    for ( unsigned k = 0; k < NB_CLIP_PLANES; ++k )
    {
        std::string var = "clip_plane" + sMath::repr(k);
        write_value(os, var, clip_plane_mode[k], clip_plane_vector[k], clip_plane_scalar[k]);
    }
    
    write_value(os, "track_fibers",  track_fibers);
    write_value(os, "auto_scale",    auto_scale);
    write_value(os, "fog",           fog_type, fog_param, fog_color);
}



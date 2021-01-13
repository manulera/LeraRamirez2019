// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_TORUS_H
#define SPACE_TORUS_H

#include "space.h"

///a torus of constant diameter centered on the origin
/**
 Space `torus` is defined by two parameters: 
 @code
    torus radius width
 @endcode
 
 With:
 - `radius` = the main radius of the torus
 - `width`  = the diameter of the torus in its crosssections.
 .
 
 */
class SpaceTorus : public Space
{
private:
    
    /// main radius
    real  bRadius;
    
    /// thickness
    real  bWidth, bWidthSqr;
    
    /// project on the backbone
    void project0(const real point[], real proj[]) const;
    
public:
        
    /// constructor
    SpaceTorus(const SpaceProp* p);
        
    /// this is called if any length has been changed
    void        resize();
    
    /// bounding box
    void        boundaries(Vector& inf, Vector& sup) const;

    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(const real point[]) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    void        project(const real point[], real proj[]) const;
    
    /// OpenGL display function, return true is display was done
    bool        display() const;
    
};

#endif

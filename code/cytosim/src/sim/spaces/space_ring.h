// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_RING_H
#define SPACE_RING_H

#include "space.h"

/// a cylinder of axis X, but without caps
/**
 Space `ring' is radial symmetric along the X-axis.
 The crosssection in the YZ plane is a disc.
 The end discs are not part of the surface.
 project() will always project on the curvy surface of the cylinder.

 @code
    ring length radius
 @endcode

 With:
 - length = half-length of the cylinder along X
 - radius = radius of the cylinder
 .

 @ingroup SpaceGroup
 */
class SpaceRing : public Space
{    
    /// apply a force directed towards the edge of the Space
    static void setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff, real len, real rad);

private:
    
    /// half the length of the cylinder (alias to mLength[0])
    real &      length;
    
    /// the radius of the ring (alias to mLength[1])
    real &      radius;
    
    /// the square of the radius (alias to mLengthSqr[1])
    real &      radiusSqr;
    
public:
        
    ///creator
    SpaceRing(const SpaceProp*);
    
    /// check number and validity of specified lengths
    void        resize() { Space::checkLengths(2, true); }
    
    /// bounding box
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(const real point[]) const;
    
    /// true if the bead is inside the Space
    bool        allInside(const real point[], real rad) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    void        project(const real point[], real proj[]) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, real rad, Meca &, real stiff) const;

    
    /// OpenGL display function, return true is display was done
    bool        display() const;
    
};

#endif


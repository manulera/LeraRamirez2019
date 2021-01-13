// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_STRIP_H
#define SPACE_STRIP_H

#include "space.h"
#include "modulo.h"

///a rectangular Space with partial periodic boundary conditions
/**
 Space `periodic` implements periodic boundary condition in all but the last dimension.
 The volume only has edge in the last dimension, and otherwise wraps on itself.
 The last dimension is Y in 2D and Z in 3D.
 
 @code
    strip sizeX sizeY sizeZ
 @endcode
 
 With:
 - sizeX = half-width along X
 - sizeY = half-width along Y
 - sizeZ = half-width along Z
 .
 
 */
class SpaceStrip : public Space
{
public:
    
    /// creator
    SpaceStrip(const SpaceProp*);

    /// check number and validity of specified lengths
    void       resize();
    
    /// initialize Modulo Object
    void       setModulo(Modulo*) const;
    
    /// bounding box
    void       boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real       volume()           const;
    
    /// true if the point is inside the Space
    bool       inside(const real point[]) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    void       project(const real point[], real proj[]) const;

    
    /// apply a force directed towards the edge of the Space
    void       setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void       setInteraction(Vector const& pos, PointExact const&, real rad, Meca &, real stiff) const;
    
    
    /// OpenGL display function, return true is display was done
    bool       display() const;
    
};

#endif


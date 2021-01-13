// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_PERIODIC_H
#define SPACE_PERIODIC_H

#include "space.h"
#include "modulo.h"

/// a rectangular Space with periodic boundary conditions
/**
 Space `periodic` implements periodic boundary condition in all dimensions.
 The volume has no edge and wraps on itself.
 
 @code
    periodic sizeX sizeY sizeZ
 @endcode
 
 With:
 - sizeX = half-width along X
 - sizeY = half-width along Y
 - sizeZ = half-width along Z
 .
 
 */
class SpacePeriodic : public Space
{
    
public:
    
    /// creator
    SpacePeriodic(const SpaceProp*);

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
    
    /// OpenGL display function, return true is display was done
    bool       display() const;
    
};

#endif


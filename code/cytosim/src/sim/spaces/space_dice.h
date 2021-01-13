// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_DICE_H
#define SPACE_DICE_H

#include "space.h"

/// A rectangle ( or a cube ) with rounded edges. 
/**
 Space `dice' is a cube with smooth edges.

 It is build by expanding an inner supporting cube of size `sizeX-radius`, etc.
 Mathematically, a point is inside the `dice' if it is at most at distance
 `radius` from the inner supporting cube.
 The dice is included in the rectangular space of sizeX, sizeY, etc.

 @code 
    dice sizeX sizeY sizeZ radius
 @endcode

 With:
 - sizeX = half-width along X
 - sizeY = half-width along Y
 - sizeZ = half-width along Z
 - radius = rounding radius of edges
 .

 Note: Dice::setInteraction() relies on project(), and numerical instabilities
 may arise in particular if `radius << size`, because determining a tangent
 plane becomes imprecise.
*/
class SpaceDice : public Space
{
public:
    
    /// the radius by which the corners are smoothed (alias to mLength[3])
    real &      radius;
    
    /// the square of the radius
    real &      radiusSqr;
    
public:
        
    /// constructor
    SpaceDice(const SpaceProp*);
    
    /// check number and validity of specified lengths
    void        resize() { Space::checkLengths(4, false); }
    
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

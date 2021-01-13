// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_CAPSULE_H
#define SPACE_CAPSULE_H

#include "space.h"

/// a spherocylinder (cylinder capped with hemispheres)
/**
 Space `capsule' is cylinder ending with hemispheres (a spherocylinder)

 @code
    capsule length radius 
 @endcode
 
 With: 
 - length: half the length of the central cylinder
 - radius: the radius of the hemisphere and central cylinder
 .

 @ingroup SpaceGroup
 */
class SpaceCapsule : public Space
{    
    /// apply a force directed towards the edge of the Space
    static void setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff, real len, real rad);

private:
    
    /// half the length of the central cylinder (alias to mLength[0])
    real &      length;
    
    /// the radius of the hemisphere (alias to mLength[1])
    real &      radius;
    
    /// the square of the radius (alias to mLengthSqr[1])
    real &      radiusSqr;

public:
        
    /// creator
    SpaceCapsule(const SpaceProp*);
    
    /// check number and validity of specified lengths
    void        resize() { Space::checkLengths(2, false); }

    /// bounding box
    void        boundaries(Vector& inf, Vector& sup) const;

    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(const real point[]) const;
    
    /// true if the bead is inside the Space
    bool        allInside(const real point[], real rad) const;
    
    /// a random position inside the volume
    Vector      randomPlace() const;

    /// set `proj` as the point on the edge that is closest to `point`
    void        project(const real point[], real proj[]) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, real rad, Meca &, real stiff) const;
    
    /// openGL display function, return true is display was done
    bool        display() const;
    
};

#endif

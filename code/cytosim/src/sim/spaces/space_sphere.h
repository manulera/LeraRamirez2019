// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_SPHERE_H
#define SPACE_SPHERE_H

#include "space.h"

/// sphere centered at the origin.
/**
 Space `sphere` is a sphere centered around the origin
 
 @code 
    sphere radius
 @endcode
 
With:
 - radius = radius of the sphere
 .
 
 @ingroup SpaceGroup
 */

class SpaceSphere : public Space
{
protected:
    
    /// the radius of the sphere (set by mLength[0])
    real & radius;
    
    /// square of the radius (set by mLengthSqr[0])
    real & radiusSqr;
    
public:
    
    /// constructor
    SpaceSphere(const SpaceProp*);

    /// check number and validity of specified lengths
    void        resize() { Space::checkLengths(1, false); }

    /// bounding box
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(const real point[]) const;
    
    /// a random position inside the volume
    Vector      randomPlace() const { return Vector::randB(radius); }
    
    /// direct normal direction calculation
    Vector      normalToEdge(Vector const& pos) const { return pos.normalized(); }

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


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_FORCE_H
#define SPACE_FORCE_H

#include "space.h"

/// generates a force field
/**
 Space `force` generates a force field that applies to every objects.
 It is sufficient to create a space, for the force to be active.
 This space cannot be used for confinement.
 
 @code
 set space myspace
 {
    shape = force
 }
 
 new space myspace
 {
    force = 1 0 0
 }
 @endcode

 
 @ingroup SpaceGroup
 */
class SpaceForce : public Space
{
    /// stiffness of interaction
    real        stiffness;
    
    /// center
    Vector      center;
    
    /// force applied in every point
    Vector      force;
    
public:
    
    ///creator
    SpaceForce(const SpaceProp*, Glossary&);
    
    /// check number and validity of specified lengths
    void        resize() { Space::checkLengths(DIM, false); }
    
    /// bounding box
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(const real point[]) const { return true; }

    /// true if a sphere (center w, radius) fits in the space, edges included
    bool        allInside(const real point[], real rad) const { return true; }
    
    /// true if a sphere (center w[], radius) is entirely outside
    bool        allOutside(const real point[], real rad) const { return false; }
    
    /// set `proj` as the point on the edge that is closest to `point`
    void        project(const real point[], real proj[]) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, real rad, Meca &, real stiff) const;

    /// apply force to all objects in Meca
    void        setInteractions(Meca&, FiberSet const&) const;
    
    /// OpenGL display function, return true is display was done
    bool        display() const;
    
};

#endif



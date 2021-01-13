// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_DISC_H
#define SPACE_DISC_H

#include "space.h"

/// A disc centered at the origin, with variable radius.
/**
 Space `disc` is a disc centered around the origin
 
 @code 
    disc radius
 @endcode
 
 Forces registered with 'setInteractions' are added, and used to update the
 radius of the Space. How fast the radius changes is set by the value 'mobility'
 in SpaceProp.
 
 @ingroup SpaceGroup
 
 FJN, Strasbourg 29.01.2017
 */

class SpaceDisc : public Space
{
private:
    
    /// the radius of the sphere (set by mLength[0])
    real & radius;
    
    /// radial force
    mutable real   rForce;
    
public:
    
    /// constructor
    SpaceDisc(const SpaceProp*);

    /// check number and validity of specified lengths
    void        resize() { }

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

    
    /// add interactions to a Meca
    void	    setInteractions(Meca &, FiberSet const&) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, real rad, Meca &, real stiff) const;
    
    ///	the step function can change the radius
    void        step();

    
    /// OpenGL display function, return true is display was done
    bool        display() const;
    
};

#endif


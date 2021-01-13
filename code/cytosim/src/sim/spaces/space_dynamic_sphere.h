// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_DYNAMIC_SPHERE_H
#define SPACE_DYNAMIC_SPHERE_H

#include "space_sphere.h"

/// A disc centered at the origin, with variable radius.
/**
 Space `dynamic_sphere` is a disc or a sphere centered around the origin
 
 @code 
    disc radius
 @endcode
 
 Forces registered with 'setInteractions' are added, and used to update the
 radius of the Space. How fast the radius changes is set by the value 'mobility'
 in SpaceProp.
 
 @ingroup SpaceGroup
 
 FJN, Strasbourg 29.01.2017
 */

class SpaceDynamicSphere : public SpaceSphere
{
private:
    
    /// radial force
    mutable real   rForce;
    
public:
    
    /// constructor
    SpaceDynamicSphere(const SpaceProp*);

    /// check number and validity of specified lengths
    void        resize() { }
    
    /// add interactions to a Meca
    void	    setInteractions(Meca &, FiberSet const&) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, real rad, Meca &, real stiff) const;
    
    ///	the step function can change the radius
    void        step();

};

#endif


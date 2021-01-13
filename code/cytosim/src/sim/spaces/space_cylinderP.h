// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_CYLINDERP_H
#define SPACE_CYLINDERP_H

#include "space.h"
#include "modulo.h"

///a cylinder of axis X that is periodic along X
/**
 Space `cylinderP' is a cylinder with periodic boundary conditions
 along the X-axis. It has no ends and loops on itself like a torus,
 but without the curvature.

 @code
    cylinderP length radius
 @endcode

 With:
 - length = half-length of the cylinder along X
 - radius = radius of the cylinder
 .
 

 @ingroup SpaceGroup
 */
class SpaceCylinderP : public Space
{
private:
    
    /// half the length of the central cylinder (alias to mLength[0])
    real &      length;
    
    /// the radius of the cylinder (alias to mLength[1])
    real &      radius;
    
    /// the square of the radius (alias to mLengthSqr[1])
    real &      radiusSqr;

public:
        
    ///creator
    SpaceCylinderP(const SpaceProp*);
    
    /// check number and validity of specified lengths
    void        resize();
    
    /// initialize Modulo Object
    void        setModulo(Modulo*) const;
    
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

    
    /// OpenGL display function, return true is display was done
    bool        display() const;
    
};

#endif


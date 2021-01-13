// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_CYLINDERZ_H
#define SPACE_CYLINDERZ_H

#include "space.h"

///a cylinder of axis Z
/**
 Space `cylinderZ' is radial symmetric along the Z-axis.
 The crosssection in the XY plane is a disc.

 @code
    cylinderZ radius bottom top
 @endcode
 
 With:
 - radius = radius of cylinder
 - bottom = smallest Z
 - top = highest Z
 .

 @ingroup SpaceGroup
 */
class SpaceCylinderZ : public Space
{    
    /// apply a force directed towards the edge of the Space
    static void setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff, real, real, real);

private:
    
    /// the radius of the cylinder (alias to mLength[0])
    real &      radius;
    
    /// the square of the radius (alias to mLengthSqr[0])
    real &      radiusSqr;
    
    /// position in Z of the bottom limit (alias to mLength[1])
    real &      bottom;
    
    /// position in Z of the top limit (alias to mLength[2])
    real &      top;
    
public:
        
    ///creator
    SpaceCylinderZ(const SpaceProp*);
    
    /// check number and validity of specified lengths
    void        resize();
    
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


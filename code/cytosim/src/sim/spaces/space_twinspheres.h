// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_TWINSPHERES_H
#define SPACE_TWINSPHERES_H

#include "space.h"
class FiberSet;

/// UNFINISHED the union of two spheres placed with the common part in the center.
/**
 Space `twin_spheres` is the union of two spheres centered on the X-axis
 
 (DO NOT USE: this class is UNFINISHED)
 @code 
    sphere radiusL radiusR overlap
 @endcode
 
With:
 - radiusL = radius of the left sphere
 - radiusR = radius of the right sphere
 - overlap = amount of overlap between the two spheres
 .
 
 @ingroup SpaceGroup
 @todo: Volume calculation
 @todo: Interaction
 */

class SpaceTwinSpheres : public Space
{
private:
    
    /// the radius of the spheres (alias to mLength[1])
    real &      radiusL, & radiusLSqr;
    
    /// square of the radius (alias to mLengthSqr[1])
    real &      radiusR, & radiusRSqr;
    
    /// amount of overlap
    real &      overlap;
    
    /// position of the center of two spheres, radius of neck
    real        cenL, cenR, neck;
   
    /// project on sphere
    static real projectS(real rad, real cen, const real point[], real proj[]);

    /// project on disc
    static real projectC(real rad, const real point[], real proj[]);

public:
    
    /// constructor
    SpaceTwinSpheres(const SpaceProp*);
    
    /// check number and validity of specified lengths
    void        resize();
    
    /// bounding box
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// direct normal direction calculation
    Vector      normalToEdge(Vector const& pos) const { return pos.normalized(); }

    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(const real point[]) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    void        project(const real point[], real proj[]) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, real rad, Meca &, real stiff) const;
    
    /// add interactions to a Meca
    void        setInteractions(Meca &, FiberSet const&) const;

    /// OpenGL display function, return true is display was done
    bool        display() const;
    
};

#endif


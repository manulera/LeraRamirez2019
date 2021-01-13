// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_POLYGONZ_H
#define SPACE_POLYGONZ_H

#include "space.h"
#include "polygon.h"

/// an axisymmetric volume obtained by rotating a polygon around the Z axis
/**
 This is only valid in 3D.
 The volume is built by rotating a closed 2D polygon around the Z axis.
 
 The coordinates of the 2D polygon (X Z) are read from a file.
 The offset `shift` is added to the X-coordinate before the polygon is rotated around Z.
 Volume is estimated by Monte-Carlo, and takes an instant.

 @code
    polygonZ file_name shift_x shift_y
 @endcode

 @ingroup SpaceGroup
*/
class SpacePolygonZ : public Space
{
private:
    
    ///pointer to the points defining the polygon in 2D
    Polygon           mPoly;
        
    ///pre-calculated bounding box since this is called often
    Vector            mInf, mSup;
    
    /// Volume calculated from polygon
    real              mVolume;


public:
        
    ///creator
    SpacePolygonZ(const SpaceProp *, Glossary&);
    
    ///destructor
    ~SpacePolygonZ();
    
    /// bounding box
    void        boundaries(Vector& inf, Vector& sup) const { inf=mInf; sup=mSup; }
    
    /// the volume inside
    real        volume() const { return mVolume; }
    
    /// true if the point is inside the Space
    bool        inside(const real point[]) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    void        project(const real point[], real proj[]) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, PointExact const&, real rad, Meca &, real stiff) const;
    
    /// add interactions between fibers and reentrant corners
    void        setInteractions(Meca &, FiberSet const&) const;

    /// estimate Volume using a crude Monte-Carlo method with `cnt` calls to Space::inside()
    real        estimateVolumeZ(unsigned long cnt) const;

    /// update since length have changed
    void        resize();
    
    /// OpenGL display function
    void        displayZ(bool show_rings) const;

    /// OpenGL display function, return true is display was done
    bool        display() const;
    
};

#endif


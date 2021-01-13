// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_POLYGON_H
#define SPACE_POLYGON_H

#include "space.h"
#include "polygon.h"

/// a polygonal convex region in space
/**
 Space `polygon` implements a polygon. It works best for convex polygon.
 In 3D, and additional HEIGHT can be specified to describe a generalized 
 cylinder of axis Z, that has the 2D polygon as cross-section.
 
 The coordinates of the polygon are read from a file.

 @code
    polygon file_name HEIGHT
 @endcode

 @ingroup SpaceGroup
 @todo add SpacePolygon::setInteraction() for re-entrant corners
*/
class SpacePolygon : public Space
{
private:
    
    /// The 2D polygon object
    Polygon           mPoly;
        
    /// pre-calculated bounding box since this is called often
    Vector            mInf, mSup;
    
    /// Volume calculated from polygon
    real              mVolume;
    
    /// half the total height (set by mLength[0])
    real              mHeight;


public:
        
    ///creator
    SpacePolygon(const SpaceProp *, Glossary&);
    
    ///destructor
    ~SpacePolygon();
    
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

    /// update since length have changed
    void        resize();
    
    /// OpenGL display function, return true is display was done
    bool        display() const;
    
};

#endif


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_INFLATE_H
#define SPACE_INFLATE_H

#include "space.h"

///a Space inflated in all directions, by a uniform length (for convex shapes)
/**
 Space `inflate` is built from another Space by inflation or deflation
 
 @code
 new space ...
 {
   geometry = SPACE
   inflate = RADIUS
 }
 @endcode
 
 if `inflate > 0` the space is extended by `RADIUS`,
 and if `inflate < 0` the Space is reduced.

 The method works if the original shape is convex, and otherwise may fail.
 Inflation tends to work better than deflation.
*/

class SpaceInflate : public Space
{
private:
    
    ///original Space
    Space const* mSpace;

    /// the amount by which the corners are smoothed (alias to mLength[0])
    real &      radius;
    
    /// the square of the radius (alias to mLengthSqr[0])
    real &      radiusSqr;
    
    /// Volume calculated by Monte-Carlo
    real        mcVolume;

public:
        
    ///creator
    SpaceInflate(const SpaceProp*, const Space * space, real);
    
    ///destructor
    ~SpaceInflate();
    
    /// bounding box
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real        volume() const { assert_true(mcVolume>0); return mcVolume; }
    
    /// this is called if any length has been changed
    void        resize() { mcVolume = estimateVolume(1<<20); }

    /// true if the point is inside the Space
    bool        inside(const real point[]) const;
    
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

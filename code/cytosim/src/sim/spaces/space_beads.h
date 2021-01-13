// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_BEADS_H
#define SPACE_BEADS_H

#include "space.h"
#include "array.h"


class Bead;


/// a volume defined as the union of Beads
/** 
 Space `bead' is a volume defined as the union of Beads
 
 @code 
    beads BEAD_NAME
 @endcode
 
 The implementation is quite incomplete, as only inside() works.
 */

class SpaceBeads : public Space
{

    typedef Array<Bead*> BeadList;
    
    BeadList        mBeads;
    
    real            bbMin[3];
    real            bbMax[3];

    /// set bounding box
    void        setBoundaries();
    
public:

    /// constructor
    SpaceBeads(SpaceProp const*);
    
    /// length have changed
    void        resize();
    
    /// bounding box
    void        boundaries(Vector& inf, Vector& sup) const;
    
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

    /// find the beads
    void        step();
    
    /// OpenGL display function, return true is display was done
    bool        display() const;
    
};

#endif


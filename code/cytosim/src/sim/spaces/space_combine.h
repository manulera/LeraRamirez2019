// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_COMBINE_H
#define SPACE_COMBINE_H

#include "space.h"

/// An outer Space, from which an inner- Space is removed
/**
 SpaceCombine is a volume defined from other Space

 */
class SpaceCombine : public Space
{
private:
    
    ///Space in which objects are contained
    Space * outer;
    
    ///Space in which objects are excluded
    Space * inner;
    
public:
        
    ///creator
    SpaceCombine(Space * big, Space * small);
    
    ///destructor
    ~SpaceCombine();
    
    /// bounding box
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// true if the point is inside the Space
    bool        inside(const real point[]) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    void        project(const real point[], real proj[]) const;
    
    /// on the edge of the innerspace
    Vector      randomPlaceNearEdge(real rad, unsigned long nb_trials) const { return inner->randomPlaceNearEdge(rad, nb_trials); }

    /// OpenGL display function, return true is display was done
    bool        display() const;
    
};

#endif

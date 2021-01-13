// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_H
#define SPACE_H

#include <string>

#include "sim.h"
#include "real.h"
#include "node.h"
#include "vector.h"
#include "object.h"
#include "common.h"
#include "modulo.h"
#include "space_prop.h"


class PointExact;
class PointInterpolated;
class FiberSet;
class Modulo;
class Simul;
class Meca;

//------------------------------------------------------------------------------

/// Defines the spatial constrains in cytosim
/**
Confined Space needs to define two important functions:\n
 - inside(x), which tells if a position is inside the space or outside,
 - project(x,p), which projects x perpendicularly on the edge of the space.
 .
The edges are considered to be inside.
*/
class Space : public Object
{
public:
    
    /// max number of dimensions
    static const unsigned DMAX = 8;

protected:
    
    /// dimensions that define the geometry
    real         mLength[DMAX];
    
    /// double of each dimension
    real         mLength2[DMAX];
    
    /// square of each dimension
    real         mLengthSqr[DMAX];
    
    /// set the length at index `d`, but also mLength2[d] and mLengthSqr[d]
    void setLength(unsigned d, real v);
    
public:
    
    /// parameters
    const SpaceProp* prop;
    
    /// constructor
    Space(const SpaceProp* p);
    
    /// destructor
    virtual ~Space();
    
    //------------------------------ BASIC -------------------------------------
    
    /// return dimension `d`
    real length(unsigned int d) const { assert_true(d<DMAX); return mLength[d]; }
 
    /// return double dimension `d`
    real length2(unsigned int d) const { assert_true(d<DMAX); return mLength2[d]; }

    /// return squared dimension `d`
    real lengthSqr(unsigned int d) const { assert_true(d<DMAX); return mLengthSqr[d]; }
    
    /// read dimensions from a string
    void readLengths(const std::string&);
    
    /// check that all `required` lengths are positive
    void checkLengths(unsigned required, int strict) const;
    
    /// change dimension `d` to `v`, and update derived variables by calling resize()
    void resize(unsigned d, real v);
    
    /// this is called if any length has been changed
    virtual void resize() {}

    /// initialize Modulo if this Space has some periodic dimensions
    virtual void setModulo(Modulo * m) const { m->disable(); }
    
    //------------------------------ OBJECT ------------------------------------
    
    /// the volume inside in 3D, or the surface area in 2D
    virtual real   volume() const = 0;

    /// return the bounds for the coordinates of the points inside the Space
    /**
     set inf as [ min(X), min(Y), min(Z) ]
     and sup as [ max(X), max(Y), max(Z) ]
     for any point (X, Y, Z) contained inside the Space.
     
     It thus defines a cuboid aligned with the main axes, and containing the entire volume.
     */
    virtual void   boundaries(Vector& inf, Vector& sup) const = 0;

    /// return the maximum absolute value of any coordinate
    real           max_extension() const;
    
    /// true if `point` is inside or on the edge of this Space
    virtual bool   inside(const real point[]) const = 0;
    
    /// set `proj` as the point on the edge that is closest to `point`
    /*
     If the edge is a smooth surface, this should correspond to the usual orthogonal projection.
     */
    virtual void   project(const real point[], real proj[]) const = 0;
    
    /// apply a force directed towards the edge of this Space, for a point located at `pos`
    virtual void   setInteraction(Vector const& pos, PointExact const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of this Space deflated by `radius`
    virtual void   setInteraction(Vector const& pos, PointExact const&, real rad, Meca &, real stiff) const;
    
#if ( 0 )
    /// apply a force directed towards the edge of this Space
    virtual void   setInteraction(Vector const&, PointInterpolated const&, Meca &, real stiff) const;

    /// apply a force directed towards the edge of this Space
    virtual void   setInteraction(PointInterpolated const&, Meca &, real stiff, Confinement conf) const;
#endif
    
    /// true if all points of the sphere (`center`, `radius`) are inside this Space
    virtual bool   allInside(const real center[], real rad) const;
    
    /// true if no point of the sphere (`center`, `radius`) is inside this Space
    virtual bool   allOutside(const real center[], real rad) const;
    
    //---------------------------- DERIVED -------------------------------------
    
    /// true if `point` is outside this Space ( defined as !inside(point) )
    bool           outside(const real point[])  const { return ! inside(point); }
    
    /// project `point` on this Space deflated by `radius`, putting the result in `proj`
    void           project(const real point[], real proj[], real rad) const;
    
    /// return the projection of `point` on edge of this Space
    Vector         project(real point[]) const { Vector P; project(point, P); return P; }
    
    
    /// the square of the distance to the edge of this Space
    real           distanceToEdgeSqr(Vector const&) const;
    
    /// the distance to the edge, always positive
    real           distanceToEdge(Vector const& pos) const { return sqrt(distanceToEdgeSqr(pos)); }
    
    /// the distance to the edge, positive if `point` is outside, and negative if inside
    real           signedDistanceToEdge(Vector const&) const;
    
    /// bring a position back inside, as if it bounced off the walls of the Space 
    void           bounce(Vector&) const;
    
    
    /// a Vector perpendicular to the space edge at `point`, directed towards the outside
    virtual Vector normalToEdge(Vector const& pos) const;
    
    /// a random position inside the volume, uniformly distributed in the volume
    virtual Vector randomPlace() const;

    /// a random position located inside and at most at distance `radius` from the edge
    virtual Vector randomPlaceNearEdge(real rad, unsigned long nb_trials = 10000) const;
    
    /// a random position located on the edge
    Vector         randomPlaceOnEdge(real rad, unsigned long nb_trials = 10000) const;
    
    /// estimate Volume using a crude Monte-Carlo method with `cnt` calls to Space::inside()
    real           estimateVolume(unsigned long cnt, bool verbose=false) const;

    //------------------------------ SIMULATION ---------------------------------
    
    /// one Monte-Carlo simulation step
    virtual void   step() {}
    
    /// add interactions to a Meca
    virtual void   setInteractions(Meca &, FiberSet const&) const {}

    //------------------------------ READ/WRITE --------------------------------
    
    /// a unique character identifying the class
    static const Tag TAG = 'e';
    
    /// return unique character identifying the class
    Tag           tag() const { return TAG; }
    
    /// return Object Property
    Property const* property() const { return prop; }
    
    /// read from file
    virtual void  read(Inputter&, Simul&, Tag);
    
    /// write to file
    virtual void  write(Outputter&) const;
    
    /// a static_cast<> of Node::next()
    Space*        next()  const  { return static_cast<Space*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Space*        prev()  const  { return static_cast<Space*>(nPrev); }
    
    //------------------------------ DISPLAY ----------------------------------
    
    /// a shape-specific openGL display function, return true is display was done
    /**
     In 2D, this should cover the inside area using polygons primities.
     in 3D, this should draw the surface of the space, using polygon primities.
     */
    virtual bool  display()  const { return false; }

    /// display the outline of a section of the box
    void          displaySection(int dim, real pos, real step) const;

};

#endif


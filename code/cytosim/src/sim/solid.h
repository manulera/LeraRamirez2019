// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SOLID_H
#define SOLID_H

#include "dim.h"
#include "array.h"
#include "object.h"
#include "point_set.h"
#include "solid_prop.h"

class Meca;
class Glossary;

/// Undeformable cloud of points
/**
 This is a PointSet behaving like a undeformable cloud of points.
 Each point can have its own radius and together they define the viscous drag
 of the Solid in the medium.
 
 \par Geometry:
 
 The ensemble can rotate and translate like a rigid body under external forces,
 but the relative configuration of the points in space is fixed: 
 the distance between any two points is constant.  

 A snapshot of the current geometry is saved in soShape[] by fixShape().
 This configuration is reapplied to the current points by reshape(). 
 reshape() however find the best isometric transformation of soShape[] 
 into the current configuration to maintain the current position and the
 current orientation of the object.
 
 \par Viscous Drag:
 
 The distance between the points, and their radii define a total drag
 coefficient according to Stokes's law applied to individual spheres.
 Points that have a radius = 0 do not induce viscous drag.
 The hydrodynamic interactions between the beads in the ensemble,
 and more advanced hydrodynamic effects are neglected.
 The drag coefficent for translation is simply the sum of Stoke's law, 
 for all points that have a radius > 0.
 The rotational drag coefficient involves the second momentum of the configuration.
 
 \par Related classes:
 
 Solid is an extension of Bead. 
 A Solid with only one point is equivalent to a Bead, but slower to simulate.
*/
class Solid : public PointSet
{
private:
    
    /// the mean of the the points weighted by their drag coefficients
    Vector         soCenter;
    
    /// the dimension used in Stoke's law to calculate mobility of each point
    real     *     soRadius;
    
    /// array to store the reference shape of the solid, as coordinates
    real     *     soShape;
    
    /// the number of points when fixShape() was last called, used for verifications.
    unsigned int   soShapeSize;
    
    /// a counter used in reshape()
    unsigned int   soReshapeTimer;
    
    /// the reduced total (all points summed) drag coefficient for translation
    real           soDrag;
    
    /// the reduced total drag coefficient for rotation
    real           soDragRot;
    
    /// matrix containing the reduced momentum of inertia
    real           soMom[DIM*DIM];
    
    /// second momentum of the reference shape
    real           soShapeSqr;
    
    /// reset private variables
    void           reset();
    
public:
    
    /// Property
    SolidProp const* prop;
    
    /// allocate memory to hold 'size' points
    virtual unsigned allocatePoints(unsigned size);
    
    /// free all memory allocated by allocatePoints()
    void        deallocatePoints();
    
    /// initialize according to options given in Glossary
    ObjectList  build(Glossary&, Simul&);
    
    /// create following the specifications in the SolidProp
    Solid(SolidProp const*);
    
    /// Copy constructor
    Solid(const Solid&);
    
    /// Assignement operator
    Solid& operator =(const Solid&);

    /// destructor
    virtual    ~Solid();
    
    //--------------------------------------------------------------------------
    
    /// sets the mobility
    void        setDragCoefficient();

    /// total translation drag-coefficient (force = drag * speed)
    real        dragCoefficient() const;
    
    /// prepare for Meca
    void        prepareMecable();

    /// prepare for constrained projection
    void        makeProjection();
    
    /// calculates the speed of points in Y, for the forces given in X
    void        setSpeedsFromForces(const real* X, real, real* Y, bool) const;
    
    /// add contribution of Brownian forces
    real        addBrownianForces(real* rhs, real const* rnd, real sc) const;
    
    /// monte-carlo step
    void        step();
    
    //--------------------------------------------------------------------------

    /// set the reference shape as a copy of the current one
    void        fixShape();
    
    /// scale the reference shape
    void        scaleShape(real, real, real);
    
    
    /// scale current shape to match the reference set in fixShape()
    void        rescale();
    
    /// restore the reference shape in the place and orientation of the current one
    void        reshape();
    
    /// set position
    void        getPoints(const real * x);

    //--------------------------------------------------------------------------

    /// add the interactions due to confinement
    void        setInteractions(Meca &) const;
    
    /// the radius of the sphere attached at point pp
    real        radius(const unsigned int pp) const { return soRadius[pp]; }
    
    /// set the radius of the sphere which is attached at point index pos
    void        radius(unsigned int pp, real radius);
        
    /// mean of all spheres weighted with their drag coefficients (or equivalently radius)
    Vector      centroid() const;
    
    /// Position of center of gravity
    Vector      position() const { return centroid(); }

    /// add a new point with a sphere
    unsigned    addSphere(Vector const& w, real radius);
    
    /// add DIM points to build a local coordinate system around the last point
    unsigned    addTriad(real arm_length);
    
    //---------------------------- next / prev ---------------------------------
    
    /// a static_cast<> of Node::next()
    Solid *     next()  const { return static_cast<Solid*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Solid *     prev()  const { return static_cast<Solid*>(nPrev); }
    
    //------------------------------ read/write --------------------------------

    /// a unique character identifying the class
    static const Tag TAG = 'd';
    
    /// return unique character identifying the class
    Tag         tag() const { return TAG; }
    
    /// return Object Property
    Property const* property() const { return prop; }
    
    ///read from file
    void        read(Inputter&, Simul&, Tag);
    
    ///write to file
    void        write(Outputter&) const;

    ///write to ostream
    void        write(std::ostream&, bool write_shape = false) const;
};

/// output operator:
std::ostream& operator << (std::ostream& os, Solid const&);

#endif

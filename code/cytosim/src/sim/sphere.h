// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPHERE_H
#define SPHERE_H

#include "dim.h"
#include "array.h"
#include "object.h"
#include "point_set.h"
#include "sphere_prop.h"

class Meca;
class Wrist;
class GraftProp;
class Glossary;

/// Spherical object with a viscous surface
/** 
 A PointSet representing a spherical object using:
 - a radius,
 - the position of the center (point index 0),
 - fixed points on the surface to keep track of the orientation,
 - mobile points on the surface.
 .

 nbRefPts 'fixed' points provide a reference frame for the sphere:
 nbRefPts=2 in 2D and  nbRefPts=4 in 3D.

 The sphere can move as a solid body by rotation and translation.
 In addition, the surface-points can move on the surface. This motion includes
 diffusion and force-induced drag and is characterized by a mobility scalar. 
 Finally, a mobile point can also carry a Single.

 This class was started by Dietrich Foethke in 2005 to simulate the nucleus of S.pombe.
 Related classes are Bead and Solid.
*/
class Sphere : public PointSet
{
public:
    
    /// number of reference points, including center: 1, 2, 4 for DIM = 1, 2 and 3
    static const unsigned nbRefPts = DIM+(DIM==3);

    /// number of reference points
    static unsigned  nbRefPoints() { return nbRefPts; }

private:
    
    /// radius
    real             spRadius;
    
    /// drag coefficients for translation and rotation
    real             spDrag, spDragRot;
    
    
    /// size of allocated projection memory
    unsigned int     spAllocated;
        
    /// radial vectors used for projecting the forces perpendicular to constraints
    real*            spProj;
        
    //--------------------------------------------------------------------------
    
    /// move the reference points such as to restore a orthogonal reference
    void             orthogonalizeRef(unsigned i);
    
    /// allocate memory for the projection matrices
    void             allocateProjection(unsigned int nb_points);
    
    /// set speed of points in Y (rigid body motion), for the forces given in X, scaled by sc
    void             setSphereSpeedsFromForces(const real* X, real sc, real* Y) const;
    
    /// add surface motions to Y, for the forces given in X, scaled by sc
    void             addSurfaceSpeedsFromForces(const real* X, real sc, real* Y) const;
    
    ///sets the number of points in the Fiber
    virtual unsigned allocatePoints(unsigned nbp);

public:
    
    /// Property
    SphereProp const* prop;
       
    //------------------- construction and destruction -------------------------
    /// create but do not initialize
    Sphere(SphereProp const* p);

    /// create and initialize following the specifications of SphereProp
    Sphere(SphereProp const* p, real radius);
    
    /// Copy constructor
    Sphere(const Sphere&);
    
    /// Assignement operator
    Sphere& operator =(const Sphere&);

    /// destructor
    virtual      ~Sphere();
    
    //-------------------------------- info ------------------------------------
    
    /// calculate mobility with piston effect
    void         setDragCoefficientPiston();
    
    /// calculate mobility with piston effect
    void         setDragCoefficientStokes();
    
    /// calculate mobility
    void         setDragCoefficient();
    
    /// total drag-coefficient of object (force = drag * speed)
    real         dragCoefficient() const { return spDrag; }
    
    /// calculate mobility and diffusion constant
    void         prepareMecable();

    /// returns position of center of gravity (the center of the sphere)
    Vector       position()        const { return posP(0); }
    
    /// radius of the sphere
    real         radius()          const { return spRadius; }

    /// change radius
    void         resize(real);
    
    /// add the interactions due to confinement
    void         setInteractions(Meca &) const;
    
    //------------------- technical functions and mathematics ------------------
        
    /// add contribution of Brownian forces
    real         addBrownianForces(real* rhs, real const* rnd, real sc) const;

    /// bring all surface points at distance spRadius from center, by moving them radially
    void         reshape();
    
    /// set position
    void         getPoints(const real * x) { PointSet::getPoints(x); reshape(); }

    /// normalize point and add center
    unsigned     addSurfacePoint(Vector const&);
    
    /// number of points on the surface
    unsigned     nbSurfacePoints() const { return nbPoints() - nbRefPts; }
    
    /// initialize according to options given in Glossary
    ObjectList   build(Glossary&, Simul&);

    //------------------- methods for the projection ---------------------------
    
    /// prepare for constrained projection
    void         makeProjection();
    
    /// calculate speed of points in Y, for the forces given in X, scaled by sc
    void         setSpeedsFromForces(const real* X, real, real* Y, bool) const;
    
    //--------------------------------------------------------------------------
    
    /// a static_cast<> of Node::next()
    Sphere * next() const { return static_cast<Sphere*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Sphere * prev() const { return static_cast<Sphere*>(nPrev); }
    
    //------------------------------ read/write --------------------------------
    
    /// a unique character identifying the class
    static const Tag TAG = 'o';
    
    /// return unique character identifying the class
    Tag          tag() const { return TAG; }
    
    /// return Object Property
    Property const* property() const { return prop; }
    
    /// write to file
    void         write(Outputter&) const;
    
    /// read from file
    void         read(Inputter&, Simul&, Tag);
};


/// output operator:
std::ostream& operator << (std::ostream& os, Sphere const&);

#endif

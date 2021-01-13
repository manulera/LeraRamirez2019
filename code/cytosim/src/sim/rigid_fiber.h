// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef RIGID_FIBER_H
#define RIGID_FIBER_H

#include "filament.h"


/**
 Enable this option to build the projection matrix explicitly.
 Alternatively, the projection is calculated directly using vectors only.
 Having two methods is useful for cross-validation, but the matrix version is SLOWER
 
 Conclusion : do not enable this normally
*/
//#define PROJECT_WITH_MATRIX


class Matrix;

/// incompressible Filament with bending elasticity
/**
 Implements the methods of a Mecable for the Filament:
 
 -# setSpeedsFromForces() includes longitudinal incompressibility,
 which means keeping successive points equidistants:
 norm( point(p+1) - point(p) ) = segmentation()
 
 -# addRigidity() implements bending elasticity.
 .
*/
class RigidFiber : public Filament
{
private:
    
    /// allocation size for projection
    size_t      rfAllocated;
    
    /// Lagrange multipliers associated with longitudinal imcompressibility
    real   *    rfLag;
    
    /// stored normalized differences of successive model-points (array of size DIM*nbSegments)
    real   *    rfDiff;
    
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
    /// local filament direction used to calculate anisotropic drag
    real   *    rfDir;
#endif
    
    /// memory allocated to hold nbPoints() values (used as temporary variables)
    real   *    rfTMP, * rfVTP;
    
#ifdef PROJECT_WITH_MATRIX
    
    /* variables used for projecting with a matrix ( rigid_fiber_projectmat.cc ) */
    
    /// projection matrix
    real   *    mtP;
    
    /// differential of projection matrix
    real   *    mtDiffP;
    
    /// intermediate of calculus
    real   *    mtJJtiJ;
    
#else
    
    /* variables used for projecting without an explicit matrix ( rigid_fiber_project.cc ) */
    
    /// J*J', a nbSegments^2 matrix. We store the diagonal and one off-diagonal
    real   *    mtJJt, * mtJJtU;
    
#endif
    
    /// vector for the projection correction of size nbSegments
    real   *    mtJJtiJforce;
    
protected:
    
    /// mobility of the points (all points have the same drag coefficient)
    real        rfDragPoint;
    
    /// rigidity scaling factor used in addRigidity()
    real        rfRigidity;
    
    /// calculate the normalized difference of successive model-points in rfDiff[]
    void        storeDirections();

private:
    
    /// reset the memory pointers for the projection
    void        buildProjection();
    
    /// allocate memory for the projection
    void        allocateProjection(unsigned int nb_points);
    
    /// free the memory for the projection
    void        destroyProjection();

public:
    
    /// Constructor
    RigidFiber();
    
    /// copy constructor
    RigidFiber(RigidFiber const&);
    
    /// copy assignment
    RigidFiber& operator=(RigidFiber const&);
    
    /// Destructor
    virtual    ~RigidFiber();
    
    ///sets the number of points in the Fiber
    virtual unsigned allocatePoints(unsigned nbp);
    
    
    /// compute Lagrange multiplier corresponding to the longitudinal tensions in the segments
    void        computeTensions(const real* force);
    
    /// longitudinal force along segment `p`
    /**
     Tensions are calculated as the Lagrange multipliers associated with the
     constrains of conserved segments lengths.
     The tension is:
     - positive when the fiber is under tension (is being pulled)
     - negative when the fiber is under compression (is being pushed)
     .
     */
    real        tension(unsigned p) const { assert_true(p+1<nbPoints()); return rfLag[p]; }
    
    /// total drag-coefficient of object (force = drag * speed)
    real        dragCoefficient() const { return  nbPoints() * rfDragPoint; }
    
    /// drag coefficient of one point
    real        dragPoint() const { return rfDragPoint; }
    
    //--------------------- Projection  / Dynamics
    
    /// prepare for projection
    void        makeProjection();
    
    /// apply the projection: Y <- s * P * X
    /** work should be of size nbSegments() at least */
    void        projectForces(const real* X, real s, real* Y, real* work) const;
    
    
    /// prepare the correction to the projection
    void        makeProjectionDiff(const real* ) const;
    
    /// add the contribution from the projection correction
    void        addProjectionDiff(const real*, real*) const;
    
    
    /// add displacements due to the Brownian motion to rhs[]
    real        addBrownianForces(real* rhs, real const* rnd, real sc) const;
    
    /// calculate the speeds from the forces, including projection
    void        setSpeedsFromForces(const real* X, real, real* Y, bool) const;
    
    //--------------------- Rigidity

    /// add the rigidity force corresponding to configuration X into vector Y
    void        addRigidity(const real* X, real* Y) const;
    
    /// add rigidity terms to a symmetric matrix
    void        addRigidityMatrix(Matrix &, int offset, int dim) const;
    
    /// add rigidity terms to upper side of matrix
    void        addRigidityUpper(real*) const;
    
    /// add terms to loop the fiber on itself
    void        loopRigidity(const real* X, real* Y) const;

};


#endif

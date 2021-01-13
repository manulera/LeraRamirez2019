// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


#ifndef MECABLE_H
#define MECABLE_H

#include "dim.h"
#include "object.h"
#include "matrix.h"
#include "buddy.h"
#include "sim.h"

class Meca;

/// Can be simulated using a Meca.
/**
 A Mecable is an Object made of points that can can be simulated in Meca.
 
 Mecable defines an interface that is implemented in Bead and PointSet.
 Mecable is also a Buddy, and can thus be part of an Organizer.
 
 */
class Mecable : public Object, public Buddy
{
protected:
    
    /// Number of points in the set
    unsigned  pSize;
    
private:
    
    /// index in the matrices and vectors using in Meca
    Matrix::index_type mIndex;
    
    /// matrix block used for preconditionning
    real *    pBlock;
    
    /// allocated size of pBlock
    unsigned  pBlockAllocated;
    
    /// actual size of pBlock
    unsigned  pBlockSize;
    
    /// flag that block is used for preconditionning
    bool      pBlockUse;
    
    /// Disabled copy constructor (@todo: write copy constructor)
    Mecable(Mecable const&);
    
    /// Disabled copy assignment (@todo: write copy assignement)
    Mecable& operator=(Mecable const&);
    
public:
        
    /// The constructor resets the pointers
    Mecable();
    
    /// Destructor de-allocates memory
    virtual ~Mecable();
    
    //--------------------------------------------------------------------------
    
    /// Number of points
    unsigned nbPoints()     const { return pSize; }
    
    /// Index of the last point = nbPoints() - 1
    unsigned lastPoint()    const { return pSize - 1; }
    
    /// Number of segments = nbPoints() - 1
    unsigned nbSegments()   const { return pSize - 1; }
    
    /// Index of the last segment = nbPoints() - 2
    unsigned lastSegment()  const { return pSize - 2; }

    /// return position of point `p`
    virtual Vector posPoint(unsigned p) const = 0;
    
    /// copy current coordinates to provided array
    virtual void   putPoints(real[]) const = 0;

    /// replace current coordinates of points by values from the provided array
    virtual void   getPoints(const real*) = 0;
    
    //--------------------------------------------------------------------------
    
    /// return force-vector on point `p` calculated at previous step by Meca
    virtual Vector netForce(unsigned p) const = 0;
    
    /// replace current forces by the ones provided as argument
    virtual void   getForces(const real* force) = 0;
    
    /// compute Lagrange multiplier corresponding to mechanical constraints
    virtual void   computeTensions(const real* force) {}

    //--------------------------------------------------------------------------
    
    /// Store the index where coordinates are located in Meca
    void           matIndex(Matrix::index_type inx)  { mIndex = inx; }
    
    /// Index in mB of the first point. the index in the vectors is DIM*matIndex()
    /** X1 is stored at DIM*matIndex(), Y1 at DIM*matIndex()+1, Z1 at DIM*matIndex()+2
     then X2, Y2, Z2...
     */
    Matrix::index_type matIndex()      const { return mIndex; }
    
    /// Ensure that block can hold a `n x n` full matrix, return true if new memory was allocated
    bool          allocateBlock(unsigned n);
   
    /// true if preconditionner block is ready
    bool          useBlock()           const { return pBlockUse; }
    
    /// change preconditionning flag
    void          useBlock(bool b)           { pBlockUse = b; }

    /// return allocated block
    real *        block()              const { return pBlock; }
    
    /// return size of block in previous iteration
    unsigned      blockSize()          const { return pBlockSize; }
    
    //--------------------------------------------------------------------------
    /// Calculate the mobility coefficient
    virtual void  setDragCoefficient() = 0;
    
    /// The total drag coefficient of the object ( force = drag * speed )
    virtual real  dragCoefficient() const = 0;

    
    /// prepare the Mecable to solve the mechanics in Meca::solve()
    /**
     This should prepare necessary variables to solve the system:
     - set rigidity coefficients, for addRigidity() to work properly
     - set drag mobility, for setSpeedsFromForces() to work,
     - set matrix/variables necessary for constrained dynamics
     .
     */
    virtual void  prepareMecable() = 0;
        
    /// Add Brownian noise terms to a force vector (sc = kT / dt)
    virtual real  addBrownianForces(real* rhs, real const* rnd, real sc) const { return INFINITY; }
    
    //--------------------------------------------------------------------------
    
    /// Add rigidity terms Y <- Y + Rigidity * X
    /**
        Rigidity can be any force acting internally to the objects
     for example, the bending rigidity of Fibers.
     This version is used to calculate the Matrix * Vector in Meca.
     */
    virtual void  addRigidity(const real* X, real* Y) const {}
    
    /// Add rigidity matrix elements (which should be symmetric) to provided matrix
    /**
       The function should add terms to the upper part of matrix `mat`, at indices starting from `offset`.
     It should fill at maximum the upper part of the diagonal block corresponding to indices [offset, offset+dim*nbPoints()].
     It should be consistent with addRigidity(), adding exactly the same terms.
     */
    virtual void  addRigidityMatrix(Matrix &, int offset, int dim) const {}

    /// Fill upper diagonal of `mat` with matrix elements
    /**
     The function should add terms to the upper part of matrix `mat`.
     The array `mat` should be square of size `DIM*nbPoints()`.
     This version is used to build the preconditionner in Meca.
     It should be consistent with addRigidity(), adding exactly the same terms.
     */
    virtual void  addRigidityUpper(real * mat) const {}

    /// Calculate speeds from given forces
    /**
        The function should perform Y <- sc * mobility * X.
        - X and Y are vectors of size DIM*nbPoints().
        - sc is a provided scalar, and mobility is known by the object. 
        - if rhs==true, the call was made with X containing the true force in the system.
     */
    virtual void  setSpeedsFromForces(const real* X, real sc, real* Y, bool rhs=false) const = 0;
    
    //--------------------------------------------------------------------------

    /// set the terms obtained from the linearization of the Projection operator, from the given forces
    /** This is enabled by a keyword ADD_PROJECTION_DIFF in meca.cc */
    virtual void  makeProjectionDiff(const real* force) const {}
    
    /// add terms from projection correction terms: Y <- Y + P' * X;
    /** This is enabled by a keyword ADD_PROJECTION_DIFF in meca.cc */
    virtual void  addProjectionDiff(const real* X, real* Y) const {}
    
    //--------------------------------------------------------------------------
    /// add the interactions (for example due to confinements)
    virtual void  setInteractions(Meca &) const {}
    
};

#endif

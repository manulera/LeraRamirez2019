// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MECA_H
#define MECA_H

#include "array.h"
#include "vector.h"
#include "matsym.h"
#include "matsparse.h"
#include "matsparsesym.h"
#include "matsparsesym1.h"
#include "matsparsesym2.h"
#include "matsparsesym4.h"
#include "matsparsesymblk.h"
//#include "matsparsebandsym.h"

class Mecable;
class PointExact;
class PointInterpolated;
class SimulProp;
class Modulo;


/// A class to calculate the motion of objects in Cytosim
/**
Meca solves the motion of objects defined by points (i.e. Mecable),
using an equation that includes terms for each interaction between Objects,
and also forces that are internal to an object, for instance bending elasticity
for Fibers, and external forces such as confinements.
The equation is formulated using linear-algebra:
 
 @code
 d vPTS/dt = mobility * mP * ( Force + mdiffP * vPTS )
 @endcode
 
 with
 
 @code
 Force = vBAS + ( mB + mC + mR ) * vPTS
 @endcode

 The equation is solved for a small increment of time `time_step`, in the presence
 of Brownian motion, and at low Reynolds number, ie. a regime in which inertial
 forces that are proportional to mass are negligible.
 
 The equation contains `DIM * nbPts` degrees of freedom, where `nbPts` is the
 total number of points in the system. It contains vectors and matrices.
 The different  terms of the equation are:
 
 - Vector vPTS containing all the Mecable coordinates (x, y, z):
   Fiber, Sphere, Solid and other PointSet. 
 
 - Vector vBAS is of same size as vPTS, and includes the constant part obtained by
   linearization of the forces. It includes for instance the positions of Single,
   calibrated random forces simulating Brownian motion, and also offsets for periodic
   boundary conditions.
 
 - Matrix mB is the isotropic part obtained after linearization of the forces.
   It operates similarly and independently on the different dimension X, Y and Z.
   mB is square of size nbPts, symmetric and sparse.
 
 - Matrix mC is the non-isotropic part obtained after linearization of the forces.
   mC is square of size DIM*nbPts, symmetric and sparse.
 .
 
 Typically, mB and mC will inherit the stiffness coefficients of the interactions, 
 while vBAS will get forces (stiffness * position). They are set by the member functions
 interLink(), interLongLink(), interSideLink(), interSlidingLink(), etc.

 - mR add the bending elasticity for RigidFiber, or other internal forces.
   mR is symmetric of size DIM*nbPts, diagonal by blocks, each block corresponding to a Fiber. 
 
 - mP applies the projection due to constrained dynamics.
   For RigidFiber, this maintains the distance between neighboring points (longitudinal incompressibility). 
   mP is symmetric of size DIM*nbPts, diagonal by blocks, each block corresponding to a Fiber. 
   mP is not actually calculated as a matrix:
   its application on each block is done by Mecable::setSpeedsFromForces()
 
 - mdiffP is a term coming from the derivative of the projection P.
   It can provide better numerical stability in some situations where the filament are stretched.
   You can however undefine ADD_PROJECTION_DIFF in meca.cc to remove mdiffP.
 .
 
 
 Note: All Links have no effect if the given PointExacts or PointInterpolated have a point
 in common, because the matrix elements would not be calcuated correctly in that case.
 Generally, such interactions are anyway not desirable. It would correspond for 
 example to a link between two point of the same segment, without effect since the segment is straight,
 or between two successive segments on the same Fiber, which at best would fold it in a non-physical way.

 */

class Meca
{
public:
 
    /// same as Matrix::index_type
    typedef Matrix::index_type index_type;
    
    /// enables graphical display of all interactions
    bool displayInteractions;

private:
    
    /// local copy of the SimulProp::time_step
    real            time_step;
    
    /// list of Mecable containing points to simulate
    Array<Mecable*> objs;
    
    /// total number of points in the system
    unsigned int    nbPts;
    
    /// size of the currently allocated memory
    unsigned int    allocated;
    
    /// max block size
    unsigned int    largestBlock;

    //--------------------------------------------------------------------------
    // Vectors of size DIM * nbPts
    
    real*  vPTS;         ///< position of the points
    real*  vSOL;         ///< positions after the dynamics has been solved
    real*  vBAS;         ///< base points of forces
    real*  vRND;         ///< vector of Gaussian random numbers
    real*  vRHS;         ///< right hand side of the dynamic system
    real*  vFOR;         ///< the calculated forces, with Brownian components
    real*  vTMP;         ///< intermediate of calculus

    //--------------------------------------------------------------------------
    
    /// true if the matrix mC is non-zero
    bool   use_mC;
    
public:
    /// isotropic symmetric part of the dynamic; square matrix of size `nbPts`
    /** 
        For interactions which have identical coefficients on the X, Y, Z subspaces
    */
    MatrixSparseSymmetric1  mB;

    
    /// non-isotropic symmetric part of the dynamic; square matrix of size `DIM*nbPts`
    /** 
        For interactions which have different coefficients on the X, Y, Z subspaces,
        or which create interactions between two different subspaces.
    */
    MatrixSparseSymmetricBlock  mC;
    
    /// base for force
    real*   base()              { return vBAS; }

    /// base for force
    real&   base(index_type ix) { return vBAS[ix]; }
    
    /// position of a point
    //Vector  pos(index_type ix) const { return Vector(vPTS+DIM*ix); }
    
private:
    
    /// allocate memory
    void  allocate(unsigned);
    
    /// release memory
    void  release();
    
    /// prepare matrices for 'solve'
    void  prepareMatrices();

    /// add the linear part of forces:  Y <- Y + ( mB + mC ) * X
    void  addLinearForces(const real* X, real* Y) const;
    
    /// add forces due to bending elasticity
    void  addRigidity(const real* X, real* Y) const;
    
    /// extract the matrix diagonal block corresponding to a Mecable
    void  extractBlock(real* res, const Mecable*, real*, real*) const;
    
    /// slowly extract the matrix diagonal block corresponding to a Mecable
    void  extractBlockS(real* res, const Mecable*) const;

    /// allocate memory, compute preconditionner and return true if completed
    int   computePreconditionner(int mode);
    
    /// compute preconditionner using the provided temporary memory
    int   computePreconditionner(Mecable*, int mode, int*, real*, real*, int);
    
public:
    

    /// constructor
    Meca();
    
    /// destructor
    ~Meca() { release(); }
    
    /// Clear list of Mecable
    void  clear();
    
    /// Add a Mecable to the list of objects to be simulated
    void  add(Mecable* ob);
    
    /// Number of Mecable
    unsigned nbMecables() const { return objs.size(); }
    
    /// number of points in the system
    unsigned nbPoints() const { return nbPts; }
    
    /// Implementation of LinearSolvers::LinearOperator
    unsigned size() const { return DIM * nbPts; }
    
    /// calculate Y = M*X, where M is the matrix associated with the system
    void multiply(const real* X, real* Y) const;

    /// apply preconditionner: Y = P*X
    void precondition(const real* X, real* Y) const;
    
    //--------------------------- FORCE ELEMENTS -------------------------------
    
    /// Add a constant force at a PointExact
    void  addForce(index_type, Vector const& force);

    /// Add a constant force at a PointExact
    void  addForce(PointExact const&, Vector const& force);
    
    /// Add a constant force at a PointInterpolated
    void  addForce(PointInterpolated const&, Vector const& force);
    
    /// Add a torque to a segment
    void  addTorque(PointInterpolated const&, Torque const& torque);
    
    /// Add a torque to constrain the segment in direction `dir`
    void  addTorqueClamp(PointInterpolated const&, Vector const& dir, real weight);
    
    /// Add an explicit torque to constrain the segment to be parallel
    void  interTorque(PointInterpolated const&, PointInterpolated const&, real weight);

    /// Add an explicit torque to constrain the segment to an angle defined by (sinus, cosinus)
    void  interTorque(PointInterpolated const&, PointInterpolated const&, real cosinus, real sinus, real weight);
    
#if (DIM == 2)
    void interTorque2D(PointInterpolated const&, PointInterpolated const&, real cosinus, real sinus, real torque_weight);
#endif

    /// Force of stiffness weight with fixed position g
    void  addPointClamp(PointExact const&, Vector, real weight);
    
    /// Force of stiffness weight with fixed position g
    void  addPointClamp(PointInterpolated const&, Vector, real weight);
    
    /// Force of stiffness `weight` centered on `g`, of length `len`
    void  addLongPointClamp(Vector const& pos, PointExact const&, Vector const& center, real len, real weight);
    
    /// Force of stiffness `weight` centered on `g`, of length `len`
    void  addLongPointClamp(PointExact const&, Vector const& center, real len, real weight);
    
    /// Force of stiffness `weight` centered on `g`, of length `len`
    void  addLongPointClamp(PointInterpolated const&, Vector const& center, real len, real weight);
    
    /// Force of stiffness `weight` centered on `g`, of length `len`, in YZ plane
    void  addLongPointClampXY(PointExact const&, real len, real weight);
    
    /// Force of stiffness `weight` centered on `g`, of length `len`, in YZ plane
    void  addLongPointClampYZ(PointExact const&, real len, real weight);
    
#if ( DIM == 2 )
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  addSidePointClamp2D(PointInterpolated const&, Vector const&, real arm, real weight);
#elif ( DIM == 3 )
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  addSidePointClamp3D(PointInterpolated const&, Vector const&, Vector const& arm, real weight);
#endif
    /// Force of stiffness weight with fixed position g, on the side
    void  addSidePointClamp(PointInterpolated const&, Vector const&, real len, real weight);
    
    /// Force of stiffness weight with a line defined by `g` and its tangent `dir`
    void  addLineClamp(PointExact const&, Vector const& g, Vector const& dir, real weight);
    
    /// Force of stiffness weight with a plane defined by `g` and its tangent `dir`
    void  addLineClamp(PointInterpolated const&, Vector const& g, Vector const& dir, real weight);
    
    /// Force of stiffness weight with a plane defined by `g` and its normal `dir`
    void  addPlaneClamp(PointExact const&, Vector const& g, Vector const& dir, real weight);
    
    /// Force of stiffness weight with a plane defined by `g` and its normal `dir`
    void  addPlaneClamp(PointInterpolated const&, Vector const& g, Vector const& dir, real weight);

    //------------ ZERO-RESTING LENGTH ELEMENTS LINKING POINTS -----------------
    
    /// zero-resting length force of stiffness weight
    void  interLink(PointExact const&, PointExact const&, real weight);
    
    /// zero-resting length force of stiffness weight
    void  interLink(PointInterpolated const&, PointExact const&, real weight);
    
    /// zero-resting length force of stiffness weight
    void  interLink(PointInterpolated const&, PointInterpolated const&, real weight);
    
    /// zero-resting length force of stiffness weight between points
    void  interLink2(PointExact const&, unsigned off, const unsigned[], const real[], real weight);
    
    /// zero-resting length force of stiffness weight between points
    void  interLink3(PointExact const&, unsigned off, const unsigned[], const real[], real weight);

    /// zero-resting length force of stiffness weight between points
    void  interLink4(PointExact const&, unsigned off, const unsigned[], const real[], real weight);
    
    /// zero-resting length force of stiffness weight between points
    void  interLink2(PointInterpolated const&, unsigned off, const unsigned[], const real[], real weight);
    
    /// zero-resting length force of stiffness weight between points
    void  interLink3(PointInterpolated const&, unsigned off, const unsigned[], const real[], real weight);

    /// zero-resting length force of stiffness weight between points
    void  interLink4(PointInterpolated const&, unsigned off, const unsigned[], const real[], real weight);

    //----------------------- ELEMENTS LINKING POINTS --------------------------

    /// Force of stiffness weight and resting length `len`
    void  interLongLink(PointExact const&, PointExact const&, real len, real weight);
    
    /// Force of stiffness weight and resting length len
    void  interLongLink(PointInterpolated const&, PointExact const&, real len, real weight);
    
    /// Force of stiffness weight and resting length len
    void  interLongLink(PointInterpolated const&, PointInterpolated const&, real len, real weight);

#if ( DIM == 2 )
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideLink2D(PointInterpolated const&, PointExact const&, real arm, real weight);
#elif ( DIM == 3 )
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideLink3D(PointInterpolated const&, PointExact const&, Vector const& arm, real weight);
    
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideLinkS(PointInterpolated const&, PointExact const&, Vector const& arm, real len, real weight);
#endif
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideLink(PointInterpolated const&, PointExact const&, real len, real weight);

    
#if ( DIM == 2 )
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideLink2D(PointInterpolated const&, PointInterpolated const&, real arm, real weight);
#elif ( DIM == 3 )
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideLinkS(PointInterpolated const&, PointInterpolated const&, Vector const& arm, real len, real weight);
#endif
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideLink(PointInterpolated const&, PointInterpolated const&, real len, real weight);

#if ( DIM == 2 )
    /// Force of stiffness weight and resting length len, on the sides of both fibers
    void  interSideSideLink2D(PointInterpolated const&, PointInterpolated const&, real len, real weight, int side1, int side2);
#endif
    /// Force of stiffness weight and resting length len, on the sides of both fibers
    void  interSideSideLink(PointInterpolated const&, PointInterpolated const&, real len, real weight);

    /// Force of stiffness weight and resting length len, which can slide on first fiber
    void  interSlidingLink(PointInterpolated const&, PointExact const&, real weight);
    
    /// Force of stiffness weight and resting length len, which can slide on first fiber
    void  interSlidingLink(PointInterpolated const&, PointInterpolated const&, real weight);

    
#if ( DIM == 2 )
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideSlidingLink2D(PointInterpolated const&, PointExact const&, real arm, real weight);

    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideSlidingLinkS(PointInterpolated const&, PointExact const&, real arm, real weight);
#elif ( DIM == 3 )
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideSlidingLinkS(PointInterpolated const&, PointExact const&, Vector const& dir, real arm, real weight);
#endif
    /// Force of stiffness weight and resting length len, on the side of first point, which can slide
    void  interSideSlidingLink(PointInterpolated const&, PointExact const&, real len, real weight);
    
    
#if ( DIM == 2 )
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideSlidingLink2D(PointInterpolated const&, PointInterpolated const&, real arm, real weight);
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideSlidingLinkS(PointInterpolated const&, PointInterpolated const&, real arm, real weight);
#elif ( DIM == 3 )
    /// Force of stiffness weight and resting length len, on the side of first fiber
    void  interSideSlidingLinkS(PointInterpolated const&, PointInterpolated const&, Vector const& dir, real arm, real weight);
#endif
    
    /// Force of stiffness weight and resting length len, on the side of first point, which can slide
    void  interSideSlidingLink(PointInterpolated const&, PointInterpolated const&, real len, real weight);
    
    
    /// Create a 3-way link with same weight
    void  interSimpleTriLink(PointInterpolated const& pt1, PointInterpolated const& pt2,PointInterpolated const& pt3, real w);
    
    /// Create a link with different weights
    void  interTriLink(PointInterpolated const& pt1, const real w1, PointInterpolated const& pt2, const real w2, PointInterpolated const& pt3, const real w3);

    
    /// Linearized Coulomb repulsive force (experimental)
    void  interCoulomb(PointExact const&, PointExact const&, real weight);
    
    
    /// Allocate the memory necessary to solve(). This must be called after the last add()
    void  prepare(SimulProp const*);
    
    /// Calculate motion of all Mecables in the system
    unsigned solve(SimulProp const*, int precondition);
    
    /// calculate Forces on Mecables and Lagrange multipliers for Fiber, without thermal motion
    void  computeForces();
    
    
    
    /// Extract the complete dynamic matrix in a standard array
    void  getSystem(int order, real * matrix) const;
    
    /// Save complete matrix in binary format
    void  dumpSystem(FILE *) const;
    
    /// Save elasticity matrix in binary format
    void  dumpElasticity(FILE *) const;
    
    /// Save projection matrix in binary format
    void  dumpProjection(FILE *) const;
    
    /// Save preconditionner in binary format
    void  dumpPreconditionner(FILE *) const;
    
    /// Save drag coefficients associated with each degree of freedom in binary format
    void  dumpDrag(FILE *) const;
    
    /// Save the object ID associated with each degree of freedom
    void  dumpObjectID(FILE *) const;
    
    /// Output vectors and matrices in various files (for debugging)
    void  dump() const;
 
    /// Output vectors and matrices in various files (for debugging)
    void  dumpSparse();
    
};

#endif


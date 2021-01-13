// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/*
 Conjugate gradient and related iterative methods
 to solve linear systems: http://www.netlib.org/templates
*/
 
#include "real.h"

/// Conjugate-Gradient with C-like interface (equivalent to LinearSolver)
/**
 The only difference between SolverC and the class Solver is calling semantic.
 The calculations are the same.
 */
namespace SolverC
{
    
    /// Conjugate Gradient, no Preconditionning
    void CG(int size, const real* rhs, real* solution,
            void (*matVect)( const real*, real*),
            int& max_iter, real& max_residual);
    
    
    /// Conjugate Gradient, with Preconditionning
    void CGP(int size, const real* rhs, real* solution,
             void (*matVect)( const real*, real* ),
             void (*precond)( const real*, real* ),
             int& max_iter, real& max_residual);
    
    
    /// Bi-Conjugate Gradient
    void BCG(int size, const real* rhs, real* solution,
             void (*matVect)( const real*, real* ),
             void (*matVectTrans)( const real*, real* ),
             int& max_iter, real& max_residual);
    
    
    /// Bi-Conjugate Gradient Stabilized
    int BCGS(int size, const real* rhs, real* solution,
             void (*matVect)( const real*, real* ),
             int& nb_iterations, real& max_residual);
    
    
    /// Bi-Conjugate Gradient Stabilized with Preconditionning
    int BCGSP(int size, const real* rhs, real* solution,
              void (*matVect)( const real*, real* ),
              void (*precond)( const real*, real* ),
              int& nb_iterations, real& max_residual);
    
    
    /// Memory allocation function
    void allocate(int size,
                  real**v1=0, real**v2=0, real**v3=0, real**v4=0,
                  real**v5=0, real**v6=0, real**v7=0, real**v8=0);
    
    /// Memory release function
    void release(real**v1=0, real**v2=0, real**v3=0, real**v4=0,
                 real**v5=0, real**v6=0, real**v7=0, real**v8=0);
    
}


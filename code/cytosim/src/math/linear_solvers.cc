// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/*
 Conjugate gradient and related iterative methods
 to solve linear systems: http://www.netlib.org/templates
*/
 
#include "linear_solvers.h"
#include <iostream>
#include "cblas.h"
#include "real.h"



inline double dot(const int size, const real* x, const real* y)
{
#ifdef REAL_IS_FLOAT
    return blas_dsdot(size, x, 1, y, 1);
#else
    return blas_xdot(size, x, 1, y, 1);
#endif
}


/**
 Here is defined which norm is used to measure the residual
 */
bool LinearSolvers::Monitor::finished(const unsigned size, const real * x)
{
#if ( 0 )
    // this is the standard Euclidian norm
    mResid = blas_xnrm2(size, x, 1);
#else
    // use the 'infinite' norm
    mResid = blas_xnrm8(size, x);
#endif
    
#if ( 1 )
    if ( mIter > mIterOld+31 )
    {
        if ( mResid >= mResidOld )
        {
            std::cerr << "Warning: slow convergence (time_step may be too big)";
            std::cerr << " residuals: " << std::scientific << mResidOld << " at iteration " << mIterOld;
            std::cerr << ", " << std::scientific << mResid << " at " << mIter << std::endl;
        }
        mIterOld = mIter;
        mResidOld = mResid;
    }
#endif
    
#if ( 0 )
    if ( mIter > 64 )
        std::cerr << "Solver step " << mIter << " residual " << mResid << std::endl;
#endif
    
    if ( mResid != mResid )
    {
        std::cerr << "Solver diverged: step " << mIter << " resid = not-a-number" << std::endl;
        mResid = INFINITY;
        return true;
    }
    
    if ( mIter > mIterMax )
        return true;
    
    return ( mResid < mResidMax );
}


/**
 By default, the preconditionning operator is the identity
 */
void LinearSolvers::LinearOperator::precondition(const real* X, real* Y) const
{
    blas_xcopy(size(), X, 1, Y, 1);
}




#pragma mark - Iterative Methods

/**
 Conjugate Gradient, no Preconditionning
 */

void LinearSolvers::CG(const LinearOperator& mat, const real* rhs, real* x, Monitor & monitor, Allocator & allocator)
{
    const unsigned size = mat.size();
    allocator.allocate(size, 4);
    real * d = allocator.bind(0);
    real * s = allocator.bind(1);
    real * r = allocator.bind(2);
    real * q = allocator.bind(3);
    
    double alpha, beta, dold, dnew;
    
    blas_xcopy( size, rhs, 1, r, 1 );
    mat.multiply( x, s );
    blas_xaxpy( size, -1, s, 1, r, 1);            //   r <- rhs - A * x
    blas_xcopy( size, r, 1, d, 1 );               //   d <- r 
    dnew = dot( size, r, r );
        
    while ( ! monitor.finished(size, r) )
    {
        mat.multiply( d, q );                     //   q = A * d
        
        alpha = dnew / dot(size, d, q);
        blas_xaxpy( size,  alpha, d, 1, x, 1 );   //   x += alpha * d
        blas_xaxpy( size, -alpha, q, 1, r, 1 );   //   r -= alpha * q
        
        dold = dnew;
        dnew = dot(size, r, r);
        beta = dnew / dold;
        blas_xscal( size, beta, d, 1 );
        blas_xaxpy( size, 1, r, 1, d, 1 );        //   d = beta * d + r
        
        ++monitor;
    }
    
    allocator.relax();
}



/**
Conjugate Gradient, with Preconditioning
*/

void LinearSolvers::CGP(const LinearOperator& mat, const real* rhs, real* x, Monitor & monitor, Allocator & allocator)
{    
    const unsigned size = mat.size();
    allocator.allocate(size, 4);
    real * d  = allocator.bind(0);
    real * s  = allocator.bind(1);
    real * r  = allocator.bind(2);
    real * q  = allocator.bind(3);

    double alpha, beta, dold, dnew;
    
    blas_xcopy( size, rhs, 1, r, 1 );
    mat.multiply( x, s );
    blas_xaxpy( size, -1, s, 1, r, 1);             //   r = rhs - M * x
    
    mat.precondition( r, d );                      //   d <- inv(M) * r
    
    dnew = dot(size, r, d);
    
    while ( ! monitor.finished(size, r) )
    {
        mat.multiply( d, q );                      //   q = M * d
        
        alpha = dnew / dot(size, d, q);
        blas_xaxpy( size,  alpha, d, 1, x, 1 );    //   x += alpha * d
        blas_xaxpy( size, -alpha, q, 1, r, 1 );    //   r -= alpha * q
        
        mat.precondition( r, s );                  //   s = inv(M) * r;
        
        dold = dnew;
        dnew = dot(size, r, s);
        beta = dnew / dold;
        blas_xscal( size, beta, d, 1 );
        blas_xaxpy( size, 1, s, 1, d, 1 );         //   d = beta * d + s
        
        ++monitor;
    }
    
    allocator.relax();
}





/**
 Bi-Conjugate Gradient
*/

void LinearSolvers::BCG(const LinearOperator& mat, const real* rhs, real* x, Monitor & monitor, Allocator & allocator)
{
    const unsigned size = mat.size();
    allocator.allocate(size, 6);
    real * r  = allocator.bind(0);
    real * rb = allocator.bind(1);
    real * p  = allocator.bind(2);
    real * pb = allocator.bind(3);
    real * q  = allocator.bind(4);
    real * qb = allocator.bind(5);
    
    double alpha, beta, dold, dnew;
    
    blas_xcopy( size, rhs, 1, r, 1 );
    mat.multiply( x, rb );
    blas_xaxpy( size, -1, rb, 1, r, 1);            //   r = rhs - A * x
    
    blas_xcopy( size, r, 1, p, 1 );
    blas_xcopy( size, r, 1, rb, 1 );
    blas_xcopy( size, r, 1, pb, 1 );
    
    dnew = dot(size, rb, r);
    
    while ( ! monitor.finished(size, r) )
    {
        mat.multiply( p, q );                      //   q = A * p
        mat.transMultiply( pb, qb );               //   qb = A' * pb
        
        alpha = dnew / dot(size, pb, q);
        blas_xaxpy( size,  alpha, p, 1, x, 1 );    //   x  += alpha * p
        blas_xaxpy( size, -alpha, q, 1, r, 1 );    //   r  -= alpha * q
        blas_xaxpy( size, -alpha, qb, 1, rb, 1 );  //   rb -= alpha * qb
        
        dold = dnew;
        dnew = dot(size, r, rb);
        beta = dnew / dold;
        blas_xscal( size, beta, p, 1 );
        blas_xaxpy( size, 1, r, 1, p, 1 );         //   p  = beta * p  + r
        blas_xscal( size, beta, pb, 1 );
        blas_xaxpy( size, 1, rb, 1, pb, 1 );       //   pb = beta * pb + rb
    
        ++monitor;
    }
    
    allocator.relax();
}



/**
 Bi-Conjugate Gradient Stabilized
*/

void LinearSolvers::BCGS(const LinearOperator& mat, const real* rhs, real* x, Monitor & monitor, Allocator & allocator)
{
    const unsigned size = mat.size();
    allocator.allocate(size, 5);
    real * r      = allocator.bind(0);
    real * rtilde = allocator.bind(1);
    real * p      = allocator.bind(2);
    real * t      = allocator.bind(3);
    real * v      = allocator.bind(4);

    double rho_1 = 1, rho_2, alpha = 0, beta = 0, omega = 1;
    
    blas_xcopy( size, rhs, 1, r, 1 );
    mat.multiply( x, rtilde );
    blas_xaxpy( size, -1.0, rtilde, 1, r, 1);       // r = rhs - A * x
    blas_xcopy( size, r, 1, rtilde, 1 );
    
    while ( !monitor.finished(size, r) )
    {
        rho_2 = rho_1;
        rho_1 = dot(size, rtilde, r);
        
        if ( rho_1 == 0.0 )
        {
            monitor.finish(2, size, r);
            break;
        }
        
        beta = ( rho_1 / rho_2 ) * ( alpha / omega );
        if ( beta == 0 )
        {
            // p = r;
            blas_xcopy(size, r, 1, p, 1 );
        }
        else {
            // p = r + beta * ( p - omega * v )
            blas_xaxpy(size, -omega, v, 1, p, 1);
#ifdef __INTEL_MKL__
            blas_xaxpby(size, 1.0, r, 1, beta, p, 1);
#else
            blas_xscal(size, beta, p, 1);
            blas_xaxpy(size, 1.0, r, 1, p, 1);
#endif
        }
        
        mat.multiply( p, v );                     // v = A * p;
        alpha = rho_1 / dot(size, rtilde, v);
        
        blas_xaxpy(size, -alpha, v, 1, r, 1);     // r = r - alpha * v;
        blas_xaxpy(size,  alpha, p, 1, x, 1);     // x = x + alpha * p;
        
        mat.multiply( r, t );                     // t = A * r;
        omega = dot(size, t, r) / dot(size, t, t);
        
        if ( omega == 0.0 )
        {
            monitor.finish(3, size, r);
            break;
        }
        
        blas_xaxpy(size,  omega, r, 1, x, 1);     // x = x + omega * r;
        blas_xaxpy(size, -omega, t, 1, r, 1);     // r = r - omega * t;
        
        ++monitor;
    }

    allocator.relax();
}


/**
 Bi-Conjugate Gradient Stabilized with Preconditionning
*/

void LinearSolvers::BCGSP(const LinearOperator& mat, const real* rhs, real* x, Monitor & monitor, Allocator & allocator)
{
    const unsigned size = mat.size();
    allocator.allocate(size, 7);
    real * r      = allocator.bind(0);
    real * rtilde = allocator.bind(1);
    real * p      = allocator.bind(2);
    real * t      = allocator.bind(3);
    real * v      = allocator.bind(4);
    real * phat   = allocator.bind(5);
    real * shat   = allocator.bind(6);

    double rho_1 = 1, rho_2, alpha = 0, beta = 0, omega = 1.0, delta;
    
    blas_xcopy( size, rhs, 1, r, 1 );
    mat.multiply( x, rtilde );
    blas_xaxpy( size, -1.0, rtilde, 1, r, 1);         // r = rhs - A * x
    blas_xcopy( size, r, 1, rtilde, 1 );              // r_tilde = r
    
    while ( ! monitor.finished(size, r) )
    {
        rho_2 = rho_1;
        rho_1 = dot(size, rtilde, r);
        
        if ( rho_1 == 0.0 )
        {
            monitor.finish(2, size, r);
            break;
        }
        
        beta = ( rho_1 / rho_2 ) * ( alpha / omega );
        if ( beta == 0.0 )
        {
            // p = r;
            blas_xcopy(size, r, 1, p, 1 );
        }
        else {
            // p = r + beta * ( p - omega * v )
            blas_xaxpy(size, -omega, v, 1, p, 1);
#ifdef __INTEL_MKL__
            blas_xaxpby(size, 1.0, r, 1, beta, p, 1);
#else
            blas_xscal(size, beta, p, 1);
            blas_xaxpy(size, 1.0, r, 1, p, 1);
#endif
        }
        
        mat.precondition( p, phat );                // phat = PC * p;
        mat.multiply( phat, v );                    // v = M * phat;
        
        delta = dot(size, rtilde, v);
        if ( delta == 0.0 )
        {
            monitor.finish(4, size, r);
            break;
        }
        
        alpha = rho_1 / delta;
        blas_xaxpy(size, -alpha, v, 1, r, 1);       // r = r - alpha * v;
        blas_xaxpy(size,  alpha, phat, 1, x, 1);    // x = x + alpha * phat;
        
        if ( monitor.finished(size, r) )
            break;

        mat.precondition( r, shat );                // shat = PC * r
        mat.multiply( shat, t );                    // t = M * shat
        
        omega = dot(size, t, r) / dot(size, t, t);
        
        if ( omega == 0.0 )
        {
            monitor.finish(3, size, r);
            break;
        }
        
        blas_xaxpy(size,  omega, shat, 1, x, 1);    // x = x + omega * shat
        blas_xaxpy(size, -omega, t, 1, r, 1);       // r = r - omega * t
        
        ++monitor;
    }
    
    allocator.relax();
}


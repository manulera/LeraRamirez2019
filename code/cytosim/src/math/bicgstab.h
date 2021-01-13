// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef BICGSTAB_H
#define BICGSTAB_H

#include <iostream>
#include <cmath>

#include "real.h"
#include "cblas.h"

inline double dot(const int size, const real* x, const real* y)
{
#ifdef REAL_IS_FLOAT
    return blas_dsdot(size, x, 1, y, 1);
#else
    return blas_xdot(size, x, 1, y, 1);
#endif
}

/// Templated iterative solvers for systems of linear equations
/**
 The linear system (and the preconditionner) is defined by class LinearOperator:    
 @code
    class LinearOperator
    {
    public:
        /// size of the matrix
        unsigned int size() const;
        
        /// apply operator to a vector
        void multiply(const real*, real*) const;
        
        /// apply transposed operator to vector
        void transMultiply(const real*, real*) const;
        
        /// apply preconditionning 
        void precondition(const real*, real*) const;
    };
 @endcode
 
 The iterative solver is followed by class Monitor,
 where the desired convergence criteria can be specified.
 Monitor will also keep track of iteration counts.
 An suitable implementation of Monitor is given here.
 
 F. Nedelec, 27.03.2012 - 21.02.2013
*/
namespace LinearSolvers
{
    
    /// records the number of iterations, and the convergence
    class Monitor
    {
    private:
        
        int      mFlag;
        
        unsigned mIter,  mIterMax, mIterOld;
        
        real     mResid, mResidMax, mResidOld;
        
    public:
        
        /// set the maximum number of iterations, and the residual threshold
        Monitor(unsigned i, real r) { reset(); mIterMax = i; mResidMax = r; }
        
        /// destructor
        virtual ~Monitor() {};
        
        /// reset interation count and achieved residual
        void reset() { mFlag = 0; mIter = 0; mResid = INFINITY; mIterOld = 0; mResidOld = INFINITY; }
        
        /// increment iteration count
        void operator ++() { ++mIter; }
        
        /// the termination code
        int flag()       const { return mFlag; }
        
        /// iteration count
        int iterations() const { return mIter; }
        
        /// last achieved residual
        real residual()  const { return mResid; }
        
        /// true if achieve residual < residual threshold
        bool converged() const { return mResid < mResidMax; }
        
        /// calculate residual from `x` and return true if threshold is achieved
        virtual bool finished(unsigned size, const real* x);

        /// calculate residual from `x` and set flag to `f`
        void finish(int f, unsigned size, const real* x) { mFlag = f; finished(size, x); }
    };
    
    
    /// allocates vectors of real
    class Allocator
    {
    private:
        
        /// size of the vector to be allocated
        size_t siz;
        
        /// number of vectors allocated
        size_t alc;
        
        /// memory
        real * mem;
        
    public:
        
        /// initialize
        Allocator()  { siz = 0; alc = 0; mem = 0; }
        
        /// calls release()
        ~Allocator() { release(); }
        
        /// allocate n vectors of size `s`
        void allocate(size_t s, unsigned int n)
        {
            // Keep memory aligned to 64 bytes:
            const size_t chunk = 64 / sizeof(real);
            siz = ( s + chunk - 1 ) & ~( chunk -1 );
            size_t a = siz * n;
            if ( a > alc )
            {
                if ( mem )
                    delete[] mem;
                mem = new real[a];
                alc = a;
                //std::clog << "Allocator::allocate "<<a<<std::endl;
            }
        }
        
        /// release memory
        void release()
        {
            if ( mem )
            {
                delete[] mem;
                //std::clog << "Allocator::release "<<std::endl;
            }
            mem = 0;
            alc = 0;
            siz = 0;
        }
        
        /// called to declare that memory is not needed anymore
        void relax()
        {
            //release();
        }
        
        /// return the memory allocated for i-th vector
        real * bind(unsigned int i)
        {
            if ( mem == 0 )
                return 0;
            if ( (i+1)*siz > alc )
                return 0;
            //std::clog << "Allocator::bind " << i << " " << mem+i*siz << std::endl;
            return mem + i * siz;
        }
    };
    
    
    /// Bi-Conjugate Gradient Stabilized without Preconditionning
    template < typename LinearOperator, typename Monitor, typename Allocator >
    void BCGS(const LinearOperator& mat, const real* rhs, real* x, Monitor& monitor, Allocator& allocator)
    {
        double rho = 1.0, rho_old = 1.0, alpha = 0.0, beta = 0.0, omega = 1.0;
        
        const unsigned int size = mat.size();
        allocator.allocate(size, 5);
        real * r  = allocator.bind(0);
        real * r0 = allocator.bind(1);
        real * p  = allocator.bind(2);
        real * t  = allocator.bind(3);
        real * v  = allocator.bind(4);
        
        blas_xcopy(size, rhs, 1, r, 1);
        if ( dot(size, x, x) > 0 )
        {
            mat.multiply(x, r0);
            blas_xaxpy(size, -1.0, r0, 1, r, 1);    // r = rhs - A * x
        }
        blas_xcopy(size, r, 1, r0, 1);              // r0 = r
        
        rho = dot(size, r, r);
        blas_xcopy(size, r, 1, p, 1);

        if ( monitor.finished(size, r) )
            return;

        goto start;
        
        while ( ! monitor.finished(size, r) )
        {
            rho_old = rho;
            rho = dot(size, r0, r);
            
            if ( rho == 0.0 )
            {
#if ( 1 )
                /* The residual vector became nearly orthogonal to the
                 arbitrarily chosen direction r0, and we restart with a new r0 */
                blas_xcopy(size, rhs, 1, r, 1);
                mat.multiply(x, r0);
                blas_xaxpy(size, -1.0, r0, 1, r, 1);    // r = rhs - A * x
                blas_xcopy(size, r, 1, r0, 1);          // r0 = r
                rho = dot(size, r0, r0);
#else
                monitor.finish(2, size, r);
                break;
#endif
            }
            
            beta = ( rho / rho_old ) * ( alpha / omega );
            // p = r + beta * ( p - omega * v )
            blas_xaxpy(size, -omega, v, 1, p, 1);
#ifdef __INTEL_MKL__
            blas_xaxpby(size, 1.0, r, 1, beta, p, 1);
#else
            blas_xscal(size, beta, p, 1);
            blas_xaxpy(size, 1.0, r, 1, p, 1);
#endif
        start:
            
            mat.multiply( p, v );                     // v = A * p;
            alpha = rho / dot(size, r0, v);
            
            blas_xaxpy(size, -alpha, v, 1, r, 1);     // r = r - alpha * v;
            blas_xaxpy(size,  alpha, p, 1, x, 1);     // x = x + alpha * p;
            
            //if ( monitor.finished(size, r) )
            //    break;
            
            mat.multiply( r, t );                     // t = A * r;
            
            double tdt = dot(size, t, t);
            
            if ( tdt > 0.0 )
            {
                omega = dot(size, t, r) / tdt;
                
                if ( omega == 0.0 )
                {
                    monitor.finish(3, size, r);
                    break;
                }
                
                blas_xaxpy(size,  omega, r, 1, x, 1);    // x = x + omega * r
                blas_xaxpy(size, -omega, t, 1, r, 1);    // r = r - omega * t
            }
            else
                omega = 0.0;
            
            ++monitor;
         }
        
        allocator.relax();
    }
    
    
    /// Bi-Conjugate Gradient Stabilized with Preconditionning
    template < typename LinearOperator, typename Monitor, typename Allocator >
    void BCGSP(const LinearOperator& mat, const real* rhs, real* x, Monitor& monitor, Allocator& allocator)
    {
        double rho = 1.0, rho_old = 1.0, alpha = 0.0, beta = 0.0, omega = 1.0, delta;
        
        const unsigned int size = mat.size();
        allocator.allocate(size, 7);
        real * r    = allocator.bind(0);
        real * r0   = allocator.bind(1);
        real * p    = allocator.bind(2);
        real * t    = allocator.bind(3);
        real * v    = allocator.bind(4);
        real * phat = allocator.bind(5);
        real * shat = allocator.bind(6);
        
        blas_xcopy(size, rhs, 1, r, 1);
        if ( dot(size, x, x) > 0 )
        {
            mat.multiply(x, r0);
            blas_xaxpy(size, -1.0, r0, 1, r, 1);    // r = rhs - A * x
        }
        blas_xcopy(size, r, 1, r0, 1);              // r0 = r
        
        rho = dot(size, r, r);
        blas_xcopy(size, r, 1, p, 1);

        if ( monitor.finished(size, r) )
            return;

        goto start;

        while ( ! monitor.finished(size, r) )
        {
            rho_old = rho;
            rho = dot(size, r0, r);
            
            if ( rho == 0.0 )
            {
#if ( 1 )
                /* The residual vector became nearly orthogonal to the
                 arbitrarily chosen direction r0, and we restart with a new r0 */
                blas_xcopy(size, rhs, 1, r, 1);
                mat.multiply(x, r0);
                blas_xaxpy(size, -1.0, r0, 1, r, 1);    // r = rhs - A * x
                blas_xcopy(size, r, 1, r0, 1);          // r0 = r
                rho = dot(size, r0, r0);
#else
                monitor.finish(2, size, r);
                break;
#endif
            }
            
            beta = ( rho / rho_old ) * ( alpha / omega );
            // p = r + beta * ( p - omega * v )
            blas_xaxpy(size, -omega, v, 1, p, 1);
#ifdef __INTEL_MKL__
            blas_xaxpby(size, 1.0, r, 1, beta, p, 1);
#else
            blas_xscal(size, beta, p, 1);
            blas_xaxpy(size, 1.0, r, 1, p, 1);
#endif
        start:
            
            mat.precondition( p, phat );                // phat = PC * p;
            mat.multiply( phat, v );                    // v = M * phat;
            
            delta = dot(size, r0, v);
            if ( delta == 0.0 )
            {
                monitor.finish(4, size, r);
                break;
            }
            
            alpha = rho / delta;
            blas_xaxpy(size, -alpha, v, 1, r, 1);       // r = r - alpha * v;
            blas_xaxpy(size,  alpha, phat, 1, x, 1);    // x = x + alpha * phat;

            mat.precondition( r, shat );                // shat = PC * r
            mat.multiply( shat, t );                    // t = M * shat
            
            double tdt = dot(size, t, t);
            
            if ( tdt > 0.0 )
            {
                omega = dot(size, t, r) / tdt;
            
                if ( omega == 0.0 )
                {
                    monitor.finish(3, size, r);
                    break;
                }
                
                blas_xaxpy(size,  omega, shat, 1, x, 1);    // x = x + omega * shat
                blas_xaxpy(size, -omega,    t, 1, r, 1);    // r = r - omega * t
            }
            else
                omega = 0.0;
  
            ++monitor;
        }
        
        allocator.relax();
    }
};

#endif


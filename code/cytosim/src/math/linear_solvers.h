// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


#ifndef LINEAR_SOLVERS_H
#define LINEAR_SOLVERS_H

#include <cstdlib>
#include <cmath>

#include "real.h"

/// Iterative solvers for systems of linear equations
/**
 The linear system is defined by class LinearOperator,
 which should be derived to implement its functions.
 
 The iterative solver is monitored by class Monitor,
 where the desired convergence criteria can be specified,
 and that will also keep track of iteration counts.
 
 F. Nedelec, 22.03.2012
*/
namespace LinearSolvers
{
    
    /// defines the functions that characterize the linear transformation
    class LinearOperator
    {
    public:
        /// size of the matrix
        virtual unsigned size() const = 0;
        
        /// apply operator to a vector
        virtual void multiply(const real*, real*) const = 0;
        
        /// apply transposed operator to vector
        virtual void transMultiply(const real*, real*) const {}
        
        /// apply preconditionning 
        virtual void precondition(const real*, real*) const;
    };
    
    
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
        void allocate(size_t s, unsigned n)
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
        real * bind(unsigned i)
        {
            if ( mem == 0 )
                return 0;
            if ( (i+1)*siz > alc )
                return 0;
            //std::clog << "Allocator::bind " << i << " " << mem+i*siz << std::endl;
            return mem + i * siz;
        }
    };

    //----------- iterative methods to solve a linear system:
    
    /// Conjugate Gradient
    void CG(const LinearOperator& mat, const real* rhs, real* solution, Monitor &, Allocator &);
    
    /// Conjugate Gradient, with Preconditionning
    void CGP(const LinearOperator& mat, const real* rhs, real* solution, Monitor &, Allocator &);
    
    /// Bi-Conjugate Gradient
    void BCG(const LinearOperator& mat, const real* rhs, real* solution, Monitor &, Allocator &);
    
    /// Bi-Conjugate Gradient Stabilized
    void BCGS(const LinearOperator& mat, const real* rhs, real* solution, Monitor &, Allocator &);
    
    /// Bi-Conjugate Gradient Stabilized with Preconditionning
    void BCGSP(const LinearOperator& mat, const real* rhs, real* solution, Monitor &, Allocator &);
        
};

#endif


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/*
 Conjugate gradient and related iterative methods
 to solve linear systems: http://www.netlib.org/templates
*/

#include "conjgradient.h"
#include "cblas.h"

/*
 The linear system is defined by functions provided as arguments:
 int size
 matVect(const real* X, real* Y)
 precond(const real* X, real* Y)
 */

#ifndef NORM
   #define NORM blas_xnrm2
#endif


inline double dot(const int size, const real* x, const real* y)
{
#ifdef REAL_IS_FLOAT
    return blas_dsdot(size, x, 1, y, 1);
#else
    return blas_xdot(size, x, 1, y, 1);
#endif
}

//-------all purpose allocation function:

void SolverC::allocate(int size,
                       real** vec1, real** vec2, real** vec3, real** vec4,
                       real** vec5, real** vec6, real** vec7, real** vec8)
{
    real* * vec[] = { vec1, vec2, vec3, vec4, vec5, vec6, vec7, vec8 };
    
    for ( int ii = 0; ii < 8; ++ii )
    {
        if ( vec[ii] )
        {
            if ( *vec[ii] )
                delete( *vec[ii] );
            *vec[ii] = new real[size];
            //blas_xzero(size, *vec[ii]);
        }
    }
}

void SolverC::release(real** vec1, real** vec2, real** vec3, real** vec4,
                      real** vec5, real** vec6, real** vec7, real** vec8)
{
    real* * vec[] = { vec1, vec2, vec3, vec4, vec5, vec6, vec7, vec8 };
    
    for ( int ii = 0; ii < 8; ++ii )
    {
        if ( vec[ii] && *vec[ii] )
            delete( *vec[ii] );
    }
}



//=============================================================================
//              Conjugate Gradient, no Preconditionning
//=============================================================================


void SolverC::CG(int size, const real* rhs, real* x,
                 void (*matVect)( const real*, real* ),
                 int& max_itr, real& max_res)
{
    real* d = 0, *s=0, *r=0, *q=0;
    
    allocate( size, &d, &s, &r, &q );
    
    double alpha, beta, dold, dnew;
    
    blas_xcopy( size, rhs, 1, r, 1 );
    matVect( x, s );
    blas_xaxpy( size, -1, s, 1, r, 1);            //   r <- rhs - A * x
    blas_xcopy( size, r, 1, d, 1 );               //   d <- r 
    dnew = dot( size, r, r );
    
    real res = 0;
    int  itr = 0;
    
    for ( ; itr <= max_itr; ++itr )
    {
        matVect( d, q );                          //   q = A * d
        
        alpha = dnew / dot(size, d, q);
        blas_xaxpy( size,  alpha, d, 1, x, 1 );   //   x += alpha * d
        blas_xaxpy( size, -alpha, q, 1, r, 1 );   //   r -= alpha * q
        
        dold = dnew;
        dnew = dot(size, r, r);
        beta = dnew / dold;
        blas_xscal( size, beta, d, 1 );
        blas_xaxpy( size, 1, r, 1, d, 1 );        //   d = beta * d + r
        
        real res = NORM( size, r, 1 );
        if ( res < max_res )
            break;
    }
    
    max_res = res;
    max_itr = itr;
    release(&d, &s, &r, &q);
}



//=============================================================================
//              Conjugate Gradient, with Preconditioning
//=============================================================================



void SolverC::CGP(int size, const real* rhs, real* x, 
                  void (*matVect)( const real*, real* ),
                  void (*precond)( const real*, real* ),
                  int& max_itr, real& max_res)
{
    real* d = 0, *s=0, *r=0, *q=0;
    allocate( size, &d, &s, &r, &q );

    double alpha, beta, dold, dnew;
    
    blas_xcopy( size, rhs, 1, r, 1 );
    matVect( x, s );
    blas_xaxpy( size, -1, s, 1, r, 1);             //   r = rhs - M * x
    
    precond( r, d );                               //   d <- inv(M) * r
    
    dnew = dot(size, r, d);
    
    real res = 0;
    int  itr = 0;
    
    for ( ; itr <= max_itr; ++itr )
    {
        matVect( d, q );                           //   q = M * d
        
        alpha = dnew / dot(size, d, q);
        blas_xaxpy( size,  alpha, d, 1, x, 1 );    //   x += alpha * d
        blas_xaxpy( size, -alpha, q, 1, r, 1 );    //   r -= alpha * q
        
        precond( r, s );                           //   s = inv(M) * r;
        
        dold = dnew;
        dnew = dot(size, r, s);
        beta = dnew / dold;
        blas_xscal( size, beta, d, 1 );
        blas_xaxpy( size, 1, s, 1, d, 1 );         //   d = beta * d + s
        
        real res = NORM( size, r, 1 );
        if ( res < max_res )
            break;
    }

    max_res = res;
    max_itr = itr;
    release(&d, &s, &r, &q);
}



//=============================================================================
//                      Bi-Conjugate Gradient
//=============================================================================



void SolverC::BCG(int size, const real* rhs, real* x, 
                   void (*matVect)( const real*, real* ),
                   void (*matVectTrans)( const real*, real* ),
                   int& max_itr, real& max_res)
{
    real* r=0, *rb=0, *p=0, *pb=0, *q=0, *qb=0;
    allocate( size, &r, &rb, &p, &pb, &q, &qb );
    
    real alpha, beta, dold, dnew;
    
    blas_xcopy( size, rhs, 1, r, 1 );
    matVect( x, rb );
    blas_xaxpy( size, -1, rb, 1, r, 1);            //   r = rhs - A * x
    
    blas_xcopy( size, r, 1, p, 1 );
    blas_xcopy( size, r, 1, rb, 1 );
    blas_xcopy( size, r, 1, pb, 1 );
    
    dnew = dot(size, rb, r);
    
    real res = 0;
    int  itr = 0;
    
    for ( ; itr <= max_itr; ++itr )
    {
        matVect( p, q );                           //   q = A * p
        matVectTrans( pb, qb );                    //   qb = A' * pb
        
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
        
        res = NORM( size, r, 1 );
        if ( res < max_res )
            break;
    }

    max_res = res;
    max_itr = itr;
    release(&r, &rb, &p, &pb, &q, &qb );
}



//=============================================================================
//                 Bi-Conjugate Gradient Stabilized
//=============================================================================



int SolverC::BCGS(int size, const real* rhs, real* x,
                  void (*matVect)( const real*, real* ),
                  int& max_itr, real& max_res)
{
    double rho_1, rho_2 = 1, alpha = 0, beta, omega = 1.0;
    real* r=0, *rtilde=0, *p=0, *t=0,  *v=0;
    
    allocate( size, &r, &rtilde, &p, &t, &v );
    
    blas_xcopy( size, rhs, 1, r, 1 );
    matVect( x, rtilde );
    blas_xaxpy( size, -1.0, rtilde, 1, r, 1);       // r = rhs - A * x
    blas_xcopy( size, r, 1, rtilde, 1 );
    
    real res = 0;
    int  ret = 1;
    int  itr = 0;
    
    for ( ; itr <= max_itr; ++itr )
    {
        rho_1 = dot(size, rtilde, r);
        
        if ( rho_1 == 0 )
        {
            res = NORM(size, r, 1);
            ret = 2;
            break;
        }
        
        if ( itr == 1 )
            blas_xcopy(size, r, 1, p, 1 );        // p = r;
        else {
            beta = (rho_1/rho_2) * (alpha/omega);
            blas_xaxpy(size, -omega, v, 1, p, 1);
            blas_xscal(size, beta, p, 1);
            blas_xaxpy(size, 1.0, r, 1, p, 1);    // p = r + beta*(p-omega*v);
        }
        
        matVect( p, v );                          // v = A * p;
        alpha = rho_1 / dot(size, rtilde, v);
        blas_xaxpy(size, -alpha, v, 1, r, 1);     // r = r - alpha * v;
        blas_xaxpy(size,  alpha, p, 1, x, 1);     // x = x + alpha * p;
        
        res = NORM(size, r, 1);
        if ( res < max_res )
        {
            ret = 0;
            break;
        }
        
        matVect( r, t );                          // t = A * s;
        
        omega = dot(size, t, r) / dot(size, t, t);
        blas_xaxpy(size,  omega, r, 1, x, 1);     // x = x + omega * r;
        blas_xaxpy(size, -omega, t, 1, r, 1);     // r = r - omega * t;
        
        res = NORM(size, r, 1);
        
        if ( res < max_res )
        {
            ret = 0;
            break;
        }
        
        if ( omega == 0 )
        {
            ret = 3;
            break;
        }
        
        rho_2 = rho_1;
    }

    max_res = res;
    max_itr = itr;
    release(&r, &rtilde, &p, &t, &v);
    return ret;
}



//=============================================================================
//        Bi-Conjugate Gradient Stabilized with Preconditionning
//=============================================================================


int SolverC::BCGSP(int size, const real* rhs, real* x, 
                   void (*matVect)( const real*, real* ),
                   void (*precond)( const real*, real* ),
                   int& max_itr, real& max_res)
{
    double rho_1, rho_2 = 1, alpha = 0, beta, omega = 1.0, delta;
    real* r=0, *rtilde=0, *p=0, *phat=0, *shat=0, *t=0, *v=0;
    
    allocate( size, &r, &rtilde, &p, &t, &v, &phat, &shat );
    
    blas_xcopy( size, rhs, 1, r, 1 );
    matVect( x, rtilde );
    blas_xaxpy( size, -1.0, rtilde, 1, r, 1);         // r = rhs - A * x
    
    blas_xcopy( size, r, 1, rtilde, 1 );              // r_tilde = r
     
    real res = 0;
    int  ret = 1;
    int  itr = 0;
    
    for ( ; itr <= max_itr; ++itr )
    {
        rho_1 = dot(size, rtilde, r);
        
        if ( rho_1 == 0 )
        {
            res = NORM(size, r, 1);
            ret = 2;
            break;
        }
        
        if ( itr == 1 )
            blas_xcopy(size, r, 1, p, 1 );          // p = r;
        else {
            beta = (rho_1/rho_2) * (alpha/omega);
            //we should test here the value of beta, which is scalar
            blas_xaxpy(size, -omega, v, 1, p, 1);
            blas_xscal(size, beta, p, 1);
            blas_xaxpy(size, 1.0, r, 1, p, 1);      // p = r + beta*(p-omega*v);
        }
        
        precond( p, phat );                         // phat = inv(M) * p;
        matVect( phat, v );                         // v = M * phat;
        
        //added test for failure detected by D. Foethke, Jan 2005
        delta = dot(size, rtilde, v);
        if ( delta == 0 )
        {
            res = NORM(size, r, 1);
            ret = 4;
            break;
        }
        alpha = rho_1 / delta;
        blas_xaxpy(size, -alpha, v, 1, r, 1);       // r = r - alpha * v;
        blas_xaxpy(size,  alpha, phat, 1, x, 1);    // x = x + alpha * phat;
        
        res = NORM(size, r, 1);
        if ( res < max_res )
        {
            ret = 0;
            break;
        }
        
        precond( r, shat );                         // shat = inv(M) * r
        matVect( shat, t );                         // t = M * shat
        
        omega = dot(size, t, r) / dot(size, t, t);
        blas_xaxpy(size,  omega, shat, 1, x, 1);    // x = x + omega * shat
        blas_xaxpy(size, -omega, t, 1, r, 1);       // r = r - omega * t
        
        res = NORM(size, r, 1);
        if ( res < max_res )
        {
            ret = 0;
            break;
        }
        if ( omega == 0 )
        {
            ret = 3;
            break;
        }
        rho_2 = rho_1;
    }
    
    max_res = res;
    max_itr = itr;
    release(&r, &rtilde, &p, &t, &v, &phat, &shat);
    return ret;
}


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


#include "exceptions.h"
#include "vecprint.h"


/**
 This option replaces lapack_xpttrf() and lapack_xptts2() by custom functions
 that use multiplications instead of divisions.
 */
//#define NEW_CUSTOM_DPTTS2

/**
 This is a C-translation of the reference implementation of LAPACK dpttrf()
 (that is Thomas' algorithm to factorize a tridiagonal matrix)
*/
void custom_dpttrf(int size, real* D, real* E, int* info)
{
    for ( int n = 0; n < size-1; ++n )
    {
        if ( D[n] < 0 )
        {
            *info = n;
            return;
        }
        real en = E[n];
        E[n] = en / D[n];
        D[n+1] = D[n+1] - en * E[n];
    }
}

/**
 This is a C-translation of the reference implementation of LAPACK dptts2()
 
 *     Solve A * X = B using the factorization A = L*D*L**T,
 *     overwriting each right hand side vector with its solution.
 *
 DO I = 2, N
 B( I ) = B( I ) - B( I-1 )*E( I-1 )
 CONTINUE
 B( N ) = B( N ) / D( N )
 DO I = N - 1, 1, -1
 B( I ) = B( I ) / D( I ) - B( I+1 )*E( I )
 CONTINUE
 */
void custom_dptts2(int size, const real* D, const real* E, real * B)
{
    for ( int n = 1; n < size; ++n )
        B[n] = B[n] - B[n-1] * E[n-1];
    
    B[size-1] = B[size-1] / D[size-1];
    
    for ( int n = size-2; n >= 0; --n )
        B[n] = B[n] / D[n] - B[n+1] * E[n];
}


/**
 Custom version
 */
void custom_xpttrf(int size, real* D, real* E, int* info)
{
#if ( 0 )
    lapack_xpttrf(size, D, E, info);
    if ( *info == 0 )
    {
        //invert diagonal terms:
        for ( int n = 0; n < size; ++n )
            D[n] = 1.0 / D[n];
    }
#else
    for ( int n = 0; n < size-1; ++n )
    {
        real e = E[n];
        if ( D[n] < 0 )
        {
            *info = n;
            return;
        }
        D[n] = 1.0 / D[n];
        E[n] = e * D[n];
        D[n+1] = D[n+1] - e * E[n];
    }
    if ( D[size-1] < 0 )
    {
        *info = size-1;
        return;
    }
    D[size-1] = 1.0 / D[size-1];
#endif
}


/**
 This works like rptts2, but it assumes that the elements of D have been inverted
 */
void custom_xptts2(int size, const real* D, const real* E, real * B)
{
    for ( int n = 1; n < size; ++n )
        B[n] = B[n] - B[n-1] * E[n-1];
    
    B[size-1] = B[size-1] * D[size-1];
    
    for ( int n = size-2; n >= 0; --n )
        B[n] = B[n] * D[n] - B[n+1] * E[n];
}


/**
 Custom factorization that uses different operations
 */
void custom_xpttrf_alt(int size, real* D, real* E, int* info)
{
    D[0] = 1.0 / D[0];
    
    for ( int n = 1; n < size; ++n )
        D[n] = 1.0 / ( D[n] - E[n-1] * D[n-1] * E[n-1] );
}


/**
 Custom version without divisions
 */
void custom_xptts2_alt(int size, real const* D, real const* E, real * B)
{
    B[0] = D[0] * B[0];
    
    for ( int n = 1; n < size; ++n )
        B[n] = D[n] * ( B[n] - B[n-1] * E[n-1] );
    
    for ( int n = size-2; n >= 0; --n )
        B[n] = B[n] - D[n] * E[n] * B[n+1];
}


//------------------------------------------------------------------------------
#pragma mark -


void RigidFiber::buildProjection()
{
    //reset all variables for the projections:
    rfAllocated       = 0;
    mtJJt             = 0;
    mtJJtiJforce      = 0;
}


void RigidFiber::allocateProjection(const unsigned int nbp)
{
    if ( rfAllocated < nbp )
    {
        //std::clog << reference() << "allocateProjection(" << nbp << ")\n";
        if ( mtJJt )
            delete[] mtJJt;
        
        // Keep memory aligned to 64 bytes:
        const size_t chunk = 64 / sizeof(real);
        // make a multiple of chunk to align memory:
        rfAllocated  = ( nbp + chunk - 1 ) & ~( chunk -1 );
        
        mtJJt        = new real[3*rfAllocated];
        mtJJtU       = mtJJt + rfAllocated;
        mtJJtiJforce = mtJJt + rfAllocated*2;
    }
}


void RigidFiber::destroyProjection()
{
    //std::clog << reference() << "destroyProjection\n";
    if ( mtJJt )        delete[] mtJJt;
    mtJJt        = 0;
    mtJJtU       = 0;
    mtJJtiJforce = 0;
}


void RigidFiber::makeProjection()
{
    assert_true( rfAllocated >= nbPoints() );
    
    //set the diagonal and off-diagonal of J*J'
    const unsigned nbu = nbPoints() - 2;
    const real*const diff = rfDiff;
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
   real b = 1;
#endif

    for ( unsigned jj = 0; jj < nbu; ++jj )
    {
        const real* X = diff + DIM * jj;
#if ( DIM == 2 )
        real xn = X[0]*X[2] + X[1]*X[3];
#else
        real xn = X[0]*X[3] + X[1]*X[4] + X[2]*X[5];
#endif
        
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
        real a = 0.25 * ( 1 + xn ) * ( 1 + xn );
        mtJJt[jj]  = 2 + a + b;
        mtJJtU[jj] = -xn - a;
        b = a;
#elif ( DIM == 2 )
        mtJJt[jj]  = 2 * ( X[0]*X[0] + X[1]*X[1] );
        // the diagonal term should be equal to 2, since diff[] vectors are normalized
        //mtJJt[jj]  = 2;
        mtJJtU[jj] = -xn;
#else
        mtJJt[jj]  = 2 * ( X[0]*X[0] + X[1]*X[1] + X[2]*X[2] );
        // the diagonal term should be equal to 2, since diff[] vectors are normalized
        //mtJJt[jj]  = 2;
        mtJJtU[jj] = -xn;
#endif
    }
    
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
    mtJJt[nbu] = 3 + b;
#elif ( DIM == 2 )
    const real* X = diff + DIM*nbu;
    mtJJt[nbu] = 2 * ( X[0]*X[0] + X[1]*X[1] );
    // the diagonal term should be equal to 2, since diff[] vectors are normalized
    //mtJJt[nbu] = 2;
#else
    const real* X = diff + DIM*nbu;
    mtJJt[nbu] = 2 * ( X[0]*X[0] + X[1]*X[1] + X[2]*X[2] );
    // the diagonal term should be equal to 2, since diff[] vectors are normalized
    //mtJJt[nbu] = 2;
#endif
    
    //std::clog << "D="; VecPrint::vecPrint(std::clog, nbu+1, mtJJt, 4);
    //std::clog << "E="; VecPrint::vecPrint(std::clog, nbu, mtJJtU, 4);
    
    int info = 0;
#ifdef NEW_CUSTOM_DPTTS2
    custom_xpttrf(nbu+1, mtJJt, mtJJtU, &info);
#else
    lapack_xpttrf(nbu+1, mtJJt, mtJJtU, &info);
#endif
    if ( info )
    {
        std::clog << " info = " << info << std::endl;
        std::clog << "D="; VecPrint::vecPrint(std::clog, info, mtJJt, 3);
        std::clog << "E="; VecPrint::vecPrint(std::clog, info-1, mtJJtU, 3);
        //std::clog << "X="; VecPrint::vecPrint(std::clog, DIM*(nbu+2), psPos);
        throw Exception("could not build Fiber's projection matrix");
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Perform first calculation needed by projectForces:
 tmp <- J * X
 */
void projectForcesU_(unsigned nbs, const real* diff, const real* X, real* mul)
{
#pragma simd
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        const unsigned kk = DIM * jj;
        const real *const x = X + kk;
        mul[jj] = diff[kk  ] * ( x[DIM  ] - x[0] )
                + diff[kk+1] * ( x[DIM+1] - x[1] )
#if ( DIM > 2 )
                + diff[kk+2] * ( x[DIM+2] - x[2] )
#endif
        ;
    }
}

/**
 Perform second calculation needed by projectForces:
 Y <- s * ( X + Jt * tmp )
 */
void projectForcesD_(unsigned nbs, const real* diff, const real s, const real* X, const real* mul, real* Y)
{
    const unsigned ee = DIM * nbs;
    for ( unsigned d = 0; d < DIM; ++d )
    {
        Y[   d] = s * ( X[   d] + diff[       d] * mul[    0] );
        Y[ee+d] = s * ( X[ee+d] - diff[ee-DIM+d] * mul[nbs-1] );
    }
    
    for ( unsigned jj = 1; jj < nbs; ++jj )
    {
        const unsigned kk = DIM*jj;
        Y[kk  ] = s * ( X[kk  ] + diff[kk  ] * mul[jj] - diff[kk-DIM  ] * mul[jj-1] );
        Y[kk+1] = s * ( X[kk+1] + diff[kk+1] * mul[jj] - diff[kk-DIM+1] * mul[jj-1] );
#if ( DIM > 2 )
        Y[kk+2] = s * ( X[kk+2] + diff[kk+2] * mul[jj] - diff[kk-DIM+2] * mul[jj-1] );
#endif
    }
}

#if defined __SSE3__ && ( DIM == 2 )

#include <pmmintrin.h>

typedef __m128d vec2;
#define SSE(x) _mm_##x##_pd

/**
 Perform first calculation needed by projectForces:
 */
inline void projectForcesU_SSE(unsigned nbs, const real* diff, const real* X, real* mul)
{
    const real* pM = diff;
    const real* pX = X;
    real *pT = mul;
    
    if ( nbs & 1 )
    {
        *pT++ = diff[0] * ( X[DIM] - X[0] ) + diff[1] * ( X[DIM+1] - X[1] );
        pX += DIM;
        pM += DIM;
    }
    
    vec2 y, x = SSE(load)(pX);
    
    //we calculate the terms two by two, with the vectorized operations
    real *const end = mul + nbs;
    while ( pT < end )
    {
        y = SSE(load)(pX+2);
        pX += 4;
        vec2 a = SSE(mul)(SSE(sub)(y, x), SSE(load)(pM));
        x = SSE(load)(pX);
        vec2 b = SSE(mul)(SSE(sub)(x, y), SSE(load)(pM+2));
        pM += 4;
        SSE(storeu)(pT, SSE(hadd)(a, b));
        pT += 2;
    }
}

/**
 Perform second calculation needed by projectForces:
 */
inline void projectForcesD_SSE(unsigned nbs, const real* diff, const real s, const real* X, const real* mul, real* Y)
{
    real *pY = Y;
    const real* pX = X;
    const real* pD = diff;
    
    vec2 x = SSE(load)(X);
    vec2 ss = SSE(set1)(s);
    
    real const* pM = mul;
    real const*const end = mul + nbs;
    while ( pM < end )
    {
        pX += DIM;
        vec2 y = SSE(mul)(SSE(load)(pD), SSE(load1)(pM));
        ++pM;
        pD += DIM;
        SSE(store)(pY, SSE(mul)(ss, SSE(add)(x, y)));
        pY += DIM;
        x = SSE(sub)(SSE(load)(pX), y);
    }
    SSE(store)(pY, SSE(mul)(ss, x));
}

#endif

#if defined __AVX__ && ( DIM == 2 )

#include <immintrin.h>
typedef __m256d vec4;
#define AVX(x) _mm256_##x##_pd

/**
 Perform first calculation needed by projectForces:
 tmp <- J * X
 F. Nedelec, 9.12.2016
 */
inline void projectForcesU_AVX(unsigned nbs, const real* diff, const real* X, real* tmp)
{
    const real* pM = diff;
    const real* pX = X;
    
    real *pT = tmp;
    
    if ( nbs & 1 )
    {
        *pT++ = diff[0] * ( X[DIM] - X[0] ) + diff[1] * ( X[DIM+1] - X[1] );
        pX += DIM;
        pM += DIM;
    }
    
    vec2 y, x = SSE(load)(pX);
    
    //we calculate the terms 4 by 4
    real *const end4 = tmp + nbs % 4;
    while ( pT < end4 )
    {
        y = SSE(load)(pX+2);
        pX += 4;
        vec2 a = SSE(mul)(SSE(sub)(y, x), SSE(load)(pM));
        x = SSE(load)(pX);
        vec2 b = SSE(mul)(SSE(sub)(x, y), SSE(load)(pM+2));
        pM += 4;
        SSE(storeu)(pT, SSE(hadd)(a, b));
        pT += 2;
    }
    
    real *const end8 = tmp + nbs % 8;
    //we now calculate the terms 8 by 8
    while ( pT < end8 )
    {
        vec4 a = AVX(mul)(AVX(sub)(AVX(load)(pX+2), AVX(load)(pX  )), AVX(load)(pM  ));
        vec4 b = AVX(mul)(AVX(sub)(AVX(load)(pX+6), AVX(load)(pX+4)), AVX(load)(pM+4));
        pM += 8;
        pX += 8;
#ifdef __AVX2__
        AVX(store)(pT, AVX(permute4x64)(AVX(hadd)(a, b),0xD8));
#else
        AVX(store)(pT, AVX(hadd)(AVX(permute2f128)(a,b,0x20), AVX(permute2f128)(a,b,0x31)));
#endif
        pT += 4;
    }
    
    real *const end = tmp + nbs;
#pragma ivdep
    while ( pT < end )
    {
        vec4 a = AVX(mul)(AVX(sub)(AVX(load)(pX+2 ), AVX(load)(pX   )), AVX(load)(pM   ));
        vec4 b = AVX(mul)(AVX(sub)(AVX(load)(pX+6 ), AVX(load)(pX+4 )), AVX(load)(pM+4 ));
        vec4 c = AVX(mul)(AVX(sub)(AVX(load)(pX+10), AVX(load)(pX+8 )), AVX(load)(pM+8 ));
        vec4 d = AVX(mul)(AVX(sub)(AVX(load)(pX+14), AVX(load)(pX+12)), AVX(load)(pM+12));
        pM += 16;
        pX += 16;
#ifdef __AVX2__
        AVX(store)(pT  , AVX(permute4x64)(AVX(hadd)(a, b),0xD8));
        AVX(store)(pT+4, AVX(permute4x64)(AVX(hadd)(c, d),0xD8));
#else
        AVX(store)(pT  , AVX(hadd)(AVX(permute2f128)(a,b,0x20), AVX(permute2f128)(a,b,0x31)));
        AVX(store)(pT+4, AVX(hadd)(AVX(permute2f128)(c,d,0x20), AVX(permute2f128)(c,d,0x31)));
#endif
        pT += 8;
    }
}

/**
 Perform second calculation needed by projectForces:
 Y <- s * ( X + Jt * tmp )
 F. Nedelec, 9.12.2016
 */
inline void projectForcesD_AVX(unsigned nbs, const real* diff, const real s, const real* X, const real* mul, real* Y)
{
    real *pY = Y;
    const real* pX = X;
    const real* pD = diff;
    
    vec4 cc = AVX(setzero)();
    vec4 ss = AVX(set1)(s);
    
    bool odd = nbs % 2;
    real const* pM = mul;
    real const*const end = mul + nbs - odd;
    
    while ( pM < end )
    {
        vec4 t = AVX(broadcast)((__m128d*)pM);
        vec4 x = AVX(load)(pX);
        pM += 2;
        vec4 m = AVX(shuffle)(t, t, 0b1100);
        vec4 d = AVX(mul)(m, AVX(load)(pD));
        pD += 4;
        vec4 n = AVX(permute2f128)(cc,d,0x21);
        cc = d;
        vec4 z = AVX(add)(x, AVX(sub)(d, n));
        pX += 4;
        AVX(store)(pY, AVX(mul)(ss, z));
        pY += 4;
    }
    
    vec2 c = _mm256_castpd256_pd128(AVX(permute2f128)(cc,cc,0x81));
    
    if ( odd )
    {
        assert( pM + 1 == mul + nbs );
        vec2 m = SSE(load1)(pM);
        vec2 x = SSE(mul)(m, SSE(load)(pD));
        vec2 z = SSE(add)(SSE(load)(pX), SSE(sub)(x, c));
        SSE(store)(pY, SSE(mul)(SSE(set1)(s), z));
        c = x;
        pY += 2;
        pX += 2;
    }
    
    vec2 z = SSE(sub)(SSE(load)(pX), c);
    SSE(store)(pY, SSE(mul)(SSE(set1)(s), z));
    assert( pY == Y + DIM * nbs );
    assert( pX == X + DIM * nbs );
}

#endif



/**
 Perform first calculation needed by projectForces:
 */
inline void projectForcesU_PTR(unsigned nbs, const real* diff, const real* X, real* mul)
{
    const real * pX = X + DIM;
    const real * pM = diff;
    real x3, x0 = X[0];
    real x4, x1 = X[1];
#if ( DIM >= 3 )
    real x5, x2 = X[2];
#endif
    real *const end = mul + nbs;
    
    //normally optimized version
    for ( real* pT = mul; pT < end; ++pT )
    {
        x3 = pX[0];
        x4 = pX[1];
#if ( DIM == 2 )
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1);
#elif ( DIM == 3 )
        x5 = pX[2];
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1) + pM[2] * (x5 - x2);
        x2 = x5;
#endif
        pX += DIM;
        pM += DIM;
        x0 = x3;
        x1 = x4;
    }
}



/**
 Perform first calculation needed by projectForces:
 */
inline void projectForcesU_PTR2(unsigned nbs, const real* diff, const real* X, real* mul)
{
    const real * pX = X + DIM;
    const real * pM = diff;
    real x3, x0 = X[0];
    real x4, x1 = X[1];
#if ( DIM >= 3 )
    real x5, x2 = X[2];
#endif
    real *const end = mul + nbs;
    
    //further optimization with manual loop-unrolling
    real* pT = mul;
    if ( nbs & 1 )
    {
        x3 = pX[0];
        x4 = pX[1];
#if ( DIM == 2 )
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1);
#elif ( DIM == 3 )
        x5 = pX[2];
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1) + pM[2] * (x5 - x2);
        x2 = x5;
#endif
        ++pT;
        pX += DIM;
        pM += DIM;
        x0 = x3;
        x1 = x4;
    }
    
    while ( pT < end )
    {
        x3 = pX[0];
        x4 = pX[1];
#if ( DIM == 2 )
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1);
#elif ( DIM == 3 )
        x5 = pX[2];
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1) + pM[2] * (x5 - x2);
#endif
        
#if ( DIM == 2 )
        x0 = pX[2];
        x1 = pX[3];
        pT[1] = pM[2] * (x0 - x3) + pM[3] * (x1 - x4);
#elif ( DIM == 3 )
        x0 = pX[3];
        x1 = pX[4];
        x2 = pX[5];
        pT[1] = pM[3] * (x0 - x3) + pM[4] * (x1 - x4) + pM[5] * (x2 - x5);
#endif
        
        pT += 2;
        pX += 2*DIM;
        pM += 2*DIM;
    }
    assert_true( pT == end );
}

/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD_PTR(unsigned nbs, const real* diff, const real s, const real* X, const real* mul, real* Y)
{
    // Y <- X + Jt * tmp :
    real x0 = X[0];
    real x1 = X[1];
#if ( DIM > 2 )
    real x2 = X[2];
#endif
    
    const real* pX = X+DIM;
    const real* pM = diff;
    real *pY = Y;
    real const*const end = mul+nbs;
    for ( real const* pT = mul; pT < end; ++pT )
    {
        real y0 = *pT * pM[0];
        real y1 = *pT * pM[1];
#if ( DIM > 2 )
        real y2 = *pT * pM[2];
#endif
        pM  += DIM;
        pY[0]  = s * ( x0 + y0 );
        pY[1]  = s * ( x1 + y1 );
#if ( DIM > 2 )
        pY[2]  = s * ( x2 + y2 );
#endif
        pY  += DIM;
        x0     = pX[0] - y0;
        x1     = pX[1] - y1;
#if ( DIM > 2 )
        x2     = pX[2] - y2;
#endif
        pX  += DIM;
    }
    pY[0] = s * x0;
    pY[1] = s * x1;
#if ( DIM > 2 )
    pY[2] = s * x2;
#endif
}



//------------------------------------------------------------------------------
#pragma mark -


#if (DIM==2) && !defined(REAL_IS_FLOAT)
#  if defined(__AVX2__)
#    warning "Using AVX2 implementation"
#    define projectForcesU projectForcesU_AVX
#    define projectForcesD projectForcesD_AVX
#  elif defined(__SSE3__)
#    warning "Using SSE3 implementation"
#    define projectForcesU projectForcesU_SSE
#    define projectForcesD projectForcesD_SSE
#  endif
#else
#  define projectForcesU projectForcesU_
#  define projectForcesD projectForcesD_
#endif

/**
 Using the Intel SSE code if available
 */
void RigidFiber::projectForces(const real* X, const real s, real* Y, real* tmp) const
{
    const unsigned nbs = nbSegments();
    //printf("X  "); VecPrint::vecPrint(std::clog, DIM*nbPoints(), X);

#ifdef NEW_ANISOTROPIC_FIBER_DRAG
    
    scaleTangentially(nbPoints(), X, rfDir, rfVTP);
    projectForcesU(nbs, rfDiff, rfVTP, tmp);

#else
    
    projectForcesU(nbs, rfDiff, X, tmp);
    
#endif
    
    // tmp <- inv( J * Jt ) * tmp to find the multipliers
#ifdef NEW_CUSTOM_DPTTS2
    custom_xptts2(nbs, mtJJt, mtJJtU, tmp);
#else
    lapack_xptts2(nbs, 1, mtJJt, mtJJtU, tmp, nbs);
#endif

    projectForcesD(nbs, rfDiff, s, X, tmp, Y);
    //printf("Y  "); VecPrint::vecPrint(std::clog, DIM*nbPoints(), Y);
}


void RigidFiber::computeTensions(const real* force)
{
    const unsigned nbs = nbSegments();
    
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
    
    scaleTangentially(nbPoints(), force, rfDir, rfVTP);
    projectForcesU(nbs, rfDiff, rfVTP, rfLag);
    
#else

    projectForcesU(nbs, rfDiff, force, rfLag);
    
#endif
    
    // tmp <- inv( J * Jt ) * tmp to find the multipliers
#ifdef NEW_CUSTOM_DPTTS2
    custom_xptts2(nbs, mtJJt, mtJJtU, rfLag);
#else
    lapack_xptts2(nbs, 1, mtJJt, mtJJtU, rfLag, nbs);
#endif
}



//------------------------------------------------------------------------------
#pragma mark - Projection DIFF


void RigidFiber::makeProjectionDiff( const real* force ) const
{
    const unsigned nbs = nbSegments();
    assert_true( nbs > 0 );
    
#if ( 0 )
    // Calculates the Lagrange multipliers associated with the constraints

    //mtJJtiJforce <- J * force
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        unsigned kk = DIM*jj;
        mtJJtiJforce[jj] = rfDiff[kk] * ( force[kk+DIM] - force[kk] );
        for ( unsigned d = 1; d < DIM; ++d )
            mtJJtiJforce[jj] += rfDiff[kk+d] * ( force[kk+DIM+d] - force[kk+d] );
    }
    
    // mtJJtiJforce <- inv( J * Jt ) * J * force
#ifdef NEW_CUSTOM_DPTTS2
    custom_xptts2(nbs, mtJJt, mtJJtU, mtJJtiJforce);
#else
    lapack_xptts2(nbs, 1, mtJJt, mtJJtU, mtJJtiJforce, nbs);
#endif

#endif
#if ( 0 )
    // verify that the alternative calculations gives identical results:
    real n = max_diff(nbs, mtJJtiJforce, rfLag);
    if ( n > 1e-6 )
    {
        std::clog << "Error= \n" << n << "\n";
        std::clog << "Lagrange: "; VecPrint::vecPrint(std::clog, sMath::min(12u,nbs), mtJJtiJforce);
        std::clog << "Multipl.: "; VecPrint::vecPrint(std::clog, sMath::min(12u,nbs), rfLag);
        std::clog << "\n";
    }
#endif
    
    //----- we remove compressive forces ( negative Lagrange-multipliers )
    const real sc = 1.0 / segmentation();
    
#pragma ivdep
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        if ( rfLag[jj] < 0 )
            mtJJtiJforce[jj] = 0;
        else
            mtJJtiJforce[jj] = rfLag[jj] * sc;
    }
    //std::clog << "proj_diff:"; VecPrint::vecPrint(std::clog, sMath::min(12u,nbs), mtJJtiJforce);
}


//------------------------------------------------------------------------------

//straightforward implementation:
inline void add_projection(const unsigned nbs, const real* mul, const real* X, real* Y)
{
#pragma simd
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        if ( mul[jj] )
        {
            for ( unsigned d = 0; d < DIM; ++d )
            {
                real w = mul[jj] * ( X[DIM*jj+DIM+d] - X[DIM*jj+d] );
                Y[DIM*jj    +d] += w;
                Y[DIM*jj+DIM+d] -= w;
            }
        }
    }
}

inline void add_projectionF(const unsigned nbs, const real* mul, const real* X, real* Y)
{    
#pragma ivdep
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        const real LM = mul[jj];
        if ( LM )
        {
            const unsigned ll = DIM * jj;
            const unsigned kk = DIM * jj + DIM;
            
            real w0 = LM * ( X[kk  ] - X[ll  ] );
            real w1 = LM * ( X[kk+1] - X[ll+1] );
#if ( DIM == 3 )
            real w2 = LM * ( X[kk+2] - X[ll+2] );
#endif
            
            Y[ll  ] += w0;
            Y[ll+1] += w1;
#if ( DIM == 3 )
            Y[ll+2] += w2;
#endif
            Y[kk  ] -= w0;
            Y[kk+1] -= w1;
#if ( DIM == 3 )
            Y[kk+2] -= w2;
#endif
        }
    }
}


#if defined __AVX__ && ( DIM == 2 )

#include <immintrin.h>

typedef __m128d vec2;
#define SSE(x) _mm_##x##_pd

typedef __m256d vec4;
#define AVX(x) _mm256_##x##_pd


inline void add_projectionAVX(const unsigned nbs, const real* mul, const real* X, real* Y)
{
    real * pY = Y;
    real const* pX = X;
    real const* pM = mul;
    
    if ( nbs % 2 )
    {
        vec2 m = SSE(load1)(pM);
        ++pM;
        vec2 s = SSE(mul)(SSE(sub)(SSE(load)(pX+DIM), SSE(load)(pX)), m);
        pX += DIM;
        SSE(store)(pY    , SSE(add)(SSE(load)(pY    ), s));
        SSE(store)(pY+DIM, SSE(sub)(SSE(load)(pY+DIM), s));
        pY += DIM;
    }
    
    real const*const end = mul + nbs;
    while ( pM < end )
    {
        vec4 a = AVX(broadcast)((__m128d*)pM);
        vec4 m = AVX(shuffle)(a, a, 0b1100);

        pM += DIM;
        vec4 s = AVX(mul)(m, AVX(sub)(AVX(load)(pX+2), AVX(load)(pX)));
        pX += 2*DIM;
        
        AVX(store)(pY  , AVX(add)(AVX(load)(pY  ), s));
        AVX(store)(pY+2, AVX(sub)(AVX(load)(pY+2), s));
        pY += 2*DIM;
    }
    assert_true(pM==end);
}

#endif



void RigidFiber::addProjectionDiff(const real* X, real* Y) const
{
#if ( 0 )
    // debug code to compare with default implementation
    unsigned nbp = nbPoints()*DIM;
    real * vec = new real[nbp];
    blas_xcopy(nbp, Y, 1, vec, 1);
    add_projection(nbSegments(), mtJJtiJforce, X, vec);
#endif

    
    //add_projection(nbSegments(), mtJJtiJforce, X, Y);
#if defined __AVX__ && ( DIM == 2 )
    add_projectionAVX(nbSegments(), mtJJtiJforce, X, Y);
#else
    add_projectionF(nbSegments(), mtJJtiJforce, X, Y);
#endif
    
    
#if ( 0 )
    // debug code to compare with default implementation
    real n = max_diff(nbp, Y, vec);
    if ( n > 1e-6 )
    {
        std::clog << "Error " << n << "\n";
        VecPrint::vecPrint(std::clog, sMath::min(16u,nbp), vec);
        VecPrint::vecPrint(std::clog, sMath::min(16u,nbp), Y);
    }
    delete vec;
#endif
}



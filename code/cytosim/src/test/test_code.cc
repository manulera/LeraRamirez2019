// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <sys/time.h>

#define DIM 2

#include "real.h"
#include "smath.h"
#include "tictoc.h"
#include "random.h"
#include "vecprint.h"
#include <cassert>

#include <emmintrin.h>
typedef __m128d vec2;
#define SSE(x) _mm_##x##_pd


#include <immintrin.h>
typedef __m256d vec4;
#define AVX(x) _mm256_##x##_pd


//------------------------------------------------------------------------------

extern Random RNG;

const real scalar = 1.0;

const unsigned NBS = 31;
const unsigned SIZE = DIM * ( NBS + 1 );
const unsigned ALOC = DIM * ( NBS + 3 );

real x[ALOC], y[ALOC], z[ALOC];
real diff[ALOC], pos[ALOC], res[ALOC];


void zero(int s, real *vec)
{
    for ( int i=0; i<s; ++i )
        vec[i] = 0;
}

void init()
{
    for ( unsigned ii=0; ii<ALOC; ++ii )
    {
        x[ii]    = RNG.sreal();
        y[ii]    = RNG.sreal();
        z[ii]    = RNG.sreal();
        diff[ii] = 1+RNG.preal();
        pos[ii]  = ii+0.5*RNG.preal();
        res[ii]  = 0;
    }
}

//------------------------------------------------------------------------------
#pragma mark - RIGIDITY

/*
 This is the simple implementation
 */
void add_rigidity1(const unsigned nbt, const real* X, const real rigid, real* Y)
{
#pragma ivdep
    for ( unsigned jj = 0; jj < nbt; ++jj )
    {
        real f = rigid * ( X[jj] - 2*X[jj+DIM] + X[jj+DIM*2] );
        Y[jj      ] -=   f;
        Y[jj+DIM  ] += f+f;
        Y[jj+DIM*2] -=   f;
    }
}

/*
 In this version the loop is unrolled, pointers are used
 and further optimization are made by replacing
 ( a0 -2*a1 + a2 ) by (a2-a1)-(a1-a0).
 */
void add_rigidity3(const unsigned nbt, const real* X, const real rigid, real* Y)
{
    const real * xn = X + DIM;
    
    real x0 = xn[0];
    real x1 = xn[1];
#if ( DIM == 3 )
    real x2 = xn[2];
#endif
    
    real d0 = x0 - X[0];
    real d1 = x1 - X[1];
#if ( DIM == 3 )
    real d2 = x2 - X[2];
#endif
    
    real df0, of0 = 0, odf0 = 0;
    real df1, of1 = 0, odf1 = 0;
#if ( DIM == 3 )
    real df2, of2 = 0, odf2 = 0;
#endif
    
    xn += DIM;
    
    real * yp = Y;
    real *const end = Y + nbt;
    while ( yp < end )
    {
        real e0 = *xn - x0;
        x0 = *xn;
        ++xn;
        real f0 = rigid * ( e0 - d0 );
        d0      = e0;
        df0     = f0 - of0;
        of0     = f0;
        *yp    += odf0 - df0;
        odf0    = df0;
        ++yp;
        
        real e1 = *xn - x1;
        x1 = *xn;
        ++xn;
        real f1 = rigid * ( e1 - d1 );
        d1      = e1;
        df1     = f1 - of1;
        of1     = f1;
        *yp    += odf1 - df1;
        odf1    = df1;
        ++yp;
        
#if ( DIM == 3 )
        real e2 = *xn - x2;
        x2 = *xn;
        ++xn;
        real f2 = rigid * ( e2 - d2 );
        d2      = e2;
        df2     = f2 - of2;
        of2     = f2;
        *yp    += odf2 - df2;
        odf2    = df2;
        ++yp;
#endif
    }
    
    yp[0]   += df0 + of0;
    yp[1]   += df1 + of1;
#if ( DIM == 3 )
    yp[2]   += df2 + of2;
#endif
    
    yp += DIM;
    
    yp[0] -= of0;
    yp[1] -= of1;
#if ( DIM == 3 )
    yp[2] -= of2;
#endif
}

void add_rigiditySSE(const unsigned nbt, const real* X, const real rigid, real* Y)
{
    vec2 R = SSE(set1)(rigid);
    
    real * yp = Y;
    const real * xn = X + 2*DIM;
    real *const end = Y + nbt;
    
    vec2 x   = SSE(load)(X+DIM);
    vec2 d   = SSE(sub)(x, SSE(load)(X));
    vec2 df  = SSE(setzero)();
    vec2 of  = SSE(setzero)();
    vec2 odf = SSE(setzero)();
    
    while ( yp < end )
    {
        vec2 z = SSE(load)(xn);
        xn += DIM;
        vec2 e = SSE(sub)(z, x);
        x = z;
        
        vec2 f = SSE(mul)(R, SSE(sub)(e, d));
        d  = e;
        df = SSE(sub)(f, of);
        of = f;
        SSE(store)(yp, SSE(add)(SSE(load)(yp), SSE(sub)(odf, df)));
        yp += DIM;
        odf = df;
    }
    
    SSE(store)(yp, SSE(add)(SSE(load)(yp), SSE(add)(df, of)));
    yp += DIM;
    SSE(store)(yp, SSE(sub)(SSE(load)(yp), of));
}


void add_rigidityAVX(const unsigned nbt, const real* X, const real rigid, real* Y)
{
    /**
     THIS DOES NOT PRODUCE THE CORRECT RESULT!!!
     */
    //std::cerr << "INCORRECT CODE\n";
    vec4 R = AVX(set1)(rigid);
    
    real * yp = Y;
    const real * xn = X + 2;
    real *const end = Y + nbt;
    
    vec4 x   = AVX(load)(X+2);
    vec4 d   = AVX(sub)(x, AVX(load)(X));
    vec4 df  = AVX(setzero)();
    vec4 of  = AVX(setzero)();
    vec4 odf = AVX(setzero)();
    
    while ( yp < end )
    {
        vec4 z = AVX(load)(xn);
        xn += 4;
        vec4 e = AVX(sub)(z, x);
        x = z;
        
        vec4 f = AVX(mul)(R, AVX(sub)(e, d));
        d  = e;
        df = AVX(sub)(f, of);
        of = f;
        AVX(store)(yp, AVX(add)(AVX(load)(yp), AVX(sub)(odf, df)));
        yp += 4;
        odf = df;
    }
    
    AVX(store)(yp, AVX(add)(AVX(load)(yp), AVX(add)(df, of)));
    yp += 4;
    AVX(store)(yp, AVX(sub)(AVX(load)(yp), of));
}



void add_rigidityF(const unsigned nbt, const real* X, const real rigid, real* Y)
{
    real const* E = X + nbt + DIM;
    
    for ( int d = 0; d < DIM; ++d )
    {
        Y[d]         += rigid * ( X[d+DIM] + X[d+DIM] - X[d] - X[d+2*DIM] );
        Y[nbt+DIM+d] += rigid * ( E[d-DIM] + E[d-DIM] - E[d] - E[d-2*DIM] );
    }
    
    if ( nbt == DIM )
    {
        for ( int d = 0; d < DIM; ++d )
            Y[d+DIM] += 2 * rigid * ( X[d] - X[d+DIM] - X[d+DIM] + X[d+DIM*2] );
    }
    else
    {
        for ( unsigned ii = DIM*2; ii < nbt; ++ii )
            Y[ii] += rigid * ( - X[ii-DIM*2] - 6*X[ii] + 4*( X[ii-DIM] + X[ii+DIM] ) - X[ii+DIM*2] );
        
        for ( int d = 0; d < DIM; ++d )
        {
            Y[d+DIM] += rigid * ( X[d] + X[d] - 5*X[d+DIM] + 4*X[d+DIM*2] - X[d+DIM*3] );
            Y[nbt+d] += rigid * ( E[d] + E[d] - 5*E[d-DIM] + 4*E[d-DIM*2] - E[d-DIM*3] );
        }
    }
}




inline void testR(unsigned cnt, void (*func)(const unsigned, const real*, real, real*), char const* str)
{
    unsigned nbt = DIM * ( NBS - 1 );
    TicToc::tic();
    for ( unsigned int ii=0; ii<cnt; ++ii )
    {
        func(nbt, y, scalar, x);
        func(nbt, x, scalar, z);
        func(nbt, z, scalar, x);
    }
    TicToc::toc(str);
    zero(ALOC, x);
    func(nbt, pos, scalar, x);
    VecPrint::vecPrint(std::cout, sMath::min(16u,SIZE), x);
}


void testRigidity(unsigned cnt)
{
    testR(cnt, add_rigidity1,    "add_rigidity1   ");
    testR(cnt, add_rigidity3,    "add_rigidity3   ");
    testR(cnt, add_rigidityF,    "add_rigidityF   ");
    testR(cnt, add_rigiditySSE,  "add_rigiditySSE ");
    testR(cnt, add_rigidityAVX,  "add_rigidityAVX ");
}


//------------------------------------------------------------------------------
#pragma mark - PROJECT

void projectForcesU_(unsigned nbs, const real* diff, const real* X, real* mul)
{
#pragma simd
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        const unsigned kk = DIM * jj;
        mul[jj] = diff[kk  ] * ( X[kk+DIM  ] - X[kk  ] )
                + diff[kk+1] * ( X[kk+DIM+1] - X[kk+1] )
#if ( DIM > 2 )
                + diff[kk+2] * ( X[kk+DIM+2] - X[kk+2] )
#endif
        ;
    }
}


void projectForcesU__(unsigned nbs, const real* diff, const real* X, real* mul)
{
#pragma ivdep
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        const unsigned kk = DIM * jj;
        const real* const x = X + DIM * jj;
        mul[jj] = diff[kk  ] * ( x[DIM  ] - x[0] )
                + diff[kk+1] * ( x[DIM+1] - x[1] )
#if ( DIM > 2 )
                + diff[kk+2] * ( x[DIM+2] - x[2] )
#endif
        ;
    }
}


#ifdef __SSE3__

/**
 Perform first calculation needed by projectForces:
 tmp <- J * X
 */
inline void projectForcesU_SSE(unsigned nbs, const real* diff, const real* X, real* tmp)
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

    //we calculate the terms two by two, with the vectorized operations
    real *const end = tmp + nbs;
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

#endif

#ifdef __AVX__

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
        vec4 c = AVX(permute2f128)(a,b,0x20);
        vec4 d = AVX(permute2f128)(a,b,0x31);
        AVX(store)(pT, AVX(hadd)(c, d));
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

#endif


//------------------------------------------------------------------------------
#pragma - projectForcesD

/**
 Perform second calculation needed by projectForces:
 Y <- s * ( X + Jt * tmp )
 */
void projectForcesD_(unsigned nbs, const real* diff, const real s, const real* X, const real* lag, real* Y)
{
    const unsigned ee = DIM * nbs;
    
    // end points are special cases:
    for ( unsigned d = 0; d < DIM; ++d )
    {
        Y[   d] = s * ( X[   d] + diff[       d] * lag[    0] );
        Y[ee+d] = s * ( X[ee+d] - diff[ee-DIM+d] * lag[nbs-1] );
    }

#pragma simd
    for ( unsigned jj = 1; jj < nbs; ++jj )
    {
        const unsigned kk = DIM * jj;
        Y[kk  ] = s * ( X[kk  ] + diff[kk  ] * lag[jj] - diff[kk-DIM  ] * lag[jj-1] );
        Y[kk+1] = s * ( X[kk+1] + diff[kk+1] * lag[jj] - diff[kk-DIM+1] * lag[jj-1] );
#if ( DIM > 2 )
        Y[kk+2] = s * ( X[kk+2] + diff[kk+2] * lag[jj] - diff[kk-DIM+2] * lag[jj-1] );
#endif
    }
}



/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD__(unsigned nbs, const real* diff, const real s, const real* X, const real* lag, real* Y)
{
    real a0 = diff[0] * lag[0];
    real a1 = diff[1] * lag[0];
#if ( DIM > 2 )
    real a2 = diff[2] * lag[0];
#endif
    
    Y[0] = s * ( X[0] + a0 );
    Y[1] = s * ( X[1] + a1 );
#if ( DIM > 2 )
    Y[2] = s * ( X[2] + a2 );
#endif
    
    for ( unsigned jj = 1; jj < nbs; ++jj )
    {
        const unsigned kk = DIM * jj;
        real b0 = diff[kk  ] * lag[jj];
        real b1 = diff[kk+1] * lag[jj];
        
        Y[kk  ] = s * ( X[kk  ] + b0 - a0 );
        Y[kk+1] = s * ( X[kk+1] + b1 - a1 );
        a0 = b0;
        a1 = b1;
#if ( DIM > 2 )
        real b2 = diff[kk+2] * lag[jj];
        Y[kk+2] = s * ( X[kk+2] + b2 - a2 );
        a2 = b2;
#endif
    }
    
    const unsigned ee = DIM * nbs;
    Y[ee  ] = s * ( X[ee  ] - a0 );
    Y[ee+1] = s * ( X[ee+1] - a1 );
#if ( DIM > 2 )
    Y[ee+2] = s * ( X[ee+2] - a2 );
#endif
}


/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD___(unsigned nbs, const real* diff, const real s, const real* X, const real* lag, real* Y)
{
    real a0 = X[0];
    real a1 = X[1];
#if ( DIM > 2 )
    real a2 = X[2];
#endif
    
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        const unsigned kk = DIM * jj;
        real b0 = diff[kk  ] * lag[jj];
        real b1 = diff[kk+1] * lag[jj];
#if ( DIM > 2 )
        real b2 = diff[kk+2] * lag[jj];
#endif

        Y[kk  ] = s * ( a0 + b0 );
        Y[kk+1] = s * ( a1 + b1 );
#if ( DIM > 2 )
        Y[kk+2] = s * ( a2 + b2 );
#endif
        
        a0 = X[DIM+kk  ] - b0;
        a1 = X[DIM+kk+1] - b1;
#if ( DIM > 2 )
        a2 = X[DIM+kk+2] - b2;
#endif
    }
    
    const unsigned ee = DIM * nbs;
    Y[ee  ] = s * a0;
    Y[ee+1] = s * a1;
#if ( DIM > 2 )
    Y[ee+2] = s * ( X[ee+2] - a2 );
#endif
}


/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD_PTR(unsigned nbs, const real* diff, const real s, const real* X, const real* lag, real* Y)
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
    real const*const end = lag + nbs;
    for ( real const* pT = lag; pT < end; ++pT )
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

#if defined __SSE3__ && ( DIM == 2 )

/**
 Perform second calculation needed by projectForces:
 */
inline void projectForcesD_SSE(unsigned nbs, const real* diff, const real s, const real* X, const real* mul, real* Y)
{
    real *pY = Y;
    const real* pX = X;
    const real* pD = diff;
    
    vec2 cc = SSE(load)(X);
    vec2 ss = SSE(set1)(s);
    
    real const* pM = mul;
    real const*const end = mul + nbs;
    while ( pM < end )
    {
        pX += DIM;
        vec2 x = SSE(load)(pX);
        vec2 d = SSE(mul)(SSE(load)(pD), SSE(load1)(pM));
        ++pM;
        pD += DIM;
        SSE(store)(pY, SSE(mul)(ss, SSE(add)(cc, d)));
        pY += DIM;
        cc = SSE(sub)(x, d);
    }
    SSE(store)(pY, SSE(mul)(ss, cc));
}

#endif

#if defined __AVX__ && ( DIM == 2 )


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
        vec4 m = AVX(shuffle)(t, t, 0xC);
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




inline void testU(unsigned cnt, void (*func)(unsigned, const real*, const real*, real*), char const* str)
{
    TicToc::tic();
    for ( unsigned ii=0; ii<cnt; ++ii )
    {
        func(NBS, x, y, z);
        func(NBS, z, x, y);
        func(NBS, y, z, x);
    }
    TicToc::toc(str);
    zero(ALOC, res);
    func(NBS, diff, pos, res);
    VecPrint::vecPrint(std::cout, sMath::min(12u,SIZE), res);
}


inline void testD(unsigned cnt, void (*func)(unsigned, const real*, real, const real*, const real*, real*), char const* str)
{
    TicToc::tic();
    for ( unsigned ii=0; ii<cnt; ++ii )
    {
        func(NBS, diff, 1.0, x, y, z);
        func(NBS, diff, 0.5, y, z, x);
        func(NBS, diff, 2.0, z, x, y);
    }
    TicToc::toc(str);
    zero(ALOC, x);
    func(NBS, diff, 1.0, pos, res, x);
    VecPrint::vecPrint(std::cout, sMath::min(20u,SIZE+1), x);
}


void testProjection(unsigned cnt)
{
    testU(cnt, projectForcesU_,    "projectForcesU_   ");
    testU(cnt, projectForcesU__,   "projectForcesU__  ");
    testU(cnt, projectForcesU_SSE, "projectForcesU_SSE");
#if defined __AVX__ && ( DIM == 2 )
    testU(cnt, projectForcesU_AVX, "projectForcesU_AVX");
#endif
    
    testD(cnt, projectForcesD_,    "projectForcesD_   ");
    testD(cnt, projectForcesD__,   "projectForcesD__  ");
    testD(cnt, projectForcesD___,  "projectForcesD___ ");
    //testD(cnt, projectForcesD_PTR, "projectForcesD_PTR");
    testD(cnt, projectForcesD_SSE, "projectForcesD_SSE");
#if defined __AVX__ && ( DIM == 2 )
    testD(cnt, projectForcesD_AVX, "projectForcesD_AVX");
#endif
}




int main(int argc, char* argv[])
{
    //re-seed the random number generator:
    RNG.seedTimer();
    
    init();
    testRigidity(1<<20);
    
    init();
    testProjection(1<<22);
    
    return EXIT_SUCCESS;
}

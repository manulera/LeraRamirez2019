// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "sim.h"
#include "rigid_fiber.h"
#include "cblas.h"
#include "clapack.h"
#include "matrix.h"
#include "smath.h"
#include "random.h"
//#include "vecprint.h"

extern Random RNG;

//------------------------------------------------------------------------------
RigidFiber::RigidFiber()
{
    buildProjection();
    rfDragPoint = 0;
    rfRigidity  = 0;
    rfDiff      = 0;
    rfLag       = 0;
    rfTMP       = 0;
    rfVTP       = 0;
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
    rfDir       = 0;
#endif
}


RigidFiber::~RigidFiber()
{
    destroyProjection();
    if ( rfDiff )
    {
        delete(rfDiff);
        rfDiff = 0;
        rfLag  = 0;
        rfTMP  = 0;
    }
}


//------------------------------------------------------------------------------
RigidFiber::RigidFiber(RigidFiber const&)
{
    ABORT_NOW("unfinished: cannot copy a Fiber");
}


RigidFiber& RigidFiber::operator=(RigidFiber const&)
{
    ABORT_NOW("unfinished: cannot copy a Fiber");
}


//------------------------------------------------------------------------------
unsigned RigidFiber::allocatePoints(const unsigned nbp)
{
    unsigned ms = PointSet::allocatePoints(nbp);
    /*
     if PointSet::allocatePoints() allocated memory, it will return the 
     size of the new array, and we allocate the same size for other arrays.
     */
    if ( ms )
    {
        //std::clog << "RigidFiber::allocatePoints " << ms << std::endl;
        allocateProjection(ms);
        
        // allocate memory:
        if ( rfDiff )
            delete(rfDiff);
        
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
        rfDiff = new real[ms*(4*DIM+1)];
        rfLag  = rfDiff + ms*DIM;
        rfTMP  = rfLag + ms;
        rfDir  = rfTMP + ms*DIM;
        rfVTP  = rfDir + ms*DIM;
#else
        rfDiff = new real[ms*(2*DIM+1)];
        rfLag  = rfDiff + ms*DIM;
        rfTMP  = rfLag + ms;
#endif
        
        // reset Lagrange multipliers
        for ( unsigned p = 0; p < ms; ++p )
            rfLag[p] = 0;
    }
    return ms;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 The argument should be: sc = kT / dt;
 */
real RigidFiber::addBrownianForces(real* rhs, real const* rnd, real sc) const
{
    real b = sqrt( 2 * sc * rfDragPoint );

    for ( unsigned jj = 0; jj < DIM*nbPoints(); ++jj )
        rhs[jj] += b * rnd[jj];
    
    return b / rfDragPoint;
}


//------------------------------------------------------------------------------

/**
 Calculate the normalized difference between successive model point of the fiber:
 @code
 for ( int n = 0; n < DIM*lastPoint(); ++n )
     rfDiff[n] = ( psPos[n+DIM] - psPos[n] ) / segmentation();
 @endcode
 */

void RigidFiber::storeDirections()
{
#if ( 1 )
    // assume here that successive points are correctly separated
    const real sc  = 1.0 / segmentation();
    const unsigned end = DIM * lastPoint();
    for ( unsigned p = 0; p < end; ++p )
        rfDiff[p] = sc * ( psPos[p+DIM] - psPos[p] );
#else
    for ( unsigned p = 0; p < lastPoint(); ++p )
        diffPoints(p).normalized().put(rfDiff+DIM*p);
#endif
    
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
    /*
     Calculate the average direction at each model point of the fiber:
     - for extremities, the direction of the segment is used.
     - for intermediate points, the directions of the two flanking segments are averaged
     .
     
     The result is stored in rfDir[].
     Note: rfDir[] is calculated from rfDiff[]
     */
    
    for ( unsigned d = 0; d < DIM; ++d )
    {
        rfDir[d]     = rfDiff[d];
        rfDir[d+end] = rfDiff[d+end-DIM];
    }
    
    for ( unsigned p = DIM ; p < end; ++p )
        rfDir[p] = 0.5 * ( rfDiff[p-DIM] + rfDiff[p] );

    //VecPrint::vecPrint(std::clog, last+DIM, rfDir);
#endif
}


#ifdef NEW_ANISOTROPIC_FIBER_DRAG

/**
 This will perform:
 @code
 Y = X + (TT') X
 @endcode
 
 Where T is the local direction of the fiber.
 This is used to multiply the tangential component of X by a factor 2,
 without changing the orthogonal components
 */
void scaleTangentially(unsigned nbp, const real* X, const real* dir, real * Y)
{
    const unsigned last = DIM * nbp;
    for ( unsigned dp = 0; dp < last; dp += DIM )
    {
#if ( DIM == 2 )
        real s = X[dp]*dir[dp] + X[dp+1]*dir[dp+1];
        Y[dp  ] = X[dp  ] + s * dir[dp  ];
        Y[dp+1] = X[dp+1] + s * dir[dp+1];
#elif ( DIM == 3 )
        real s = X[dp]*dir[dp] + X[dp+1]*dir[dp+1] + X[dp+2]*dir[dp+2];
        Y[dp  ] = X[dp  ] + s * dir[dp  ];
        Y[dp+1] = X[dp+1] + s * dir[dp+1];
        Y[dp+2] = X[dp+2] + s * dir[dp+2];
#endif
    }
}
#endif


/**
 If `rhs == true`, then the local array `rfLag[]` is used,
 thus saving the intermediary calculation of the Lagrange multiplier,
 that correspond to the longitudinal tension in the fiber for later use.
*/
void RigidFiber::setSpeedsFromForces(const real* X, const real sc, real* Y, bool rhs) const
{
    assert_true( X != Y );

    if ( rhs )
        projectForces(X, sc/rfDragPoint, Y, rfLag);
    else
        projectForces(X, sc/rfDragPoint, Y, rfTMP);
    
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
    scaleTangentially(nbPoints(), Y, rfDir, Y);
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Project

#if ( DIM > 1 )

#  ifdef PROJECT_WITH_MATRIX
#     include "rigid_fiber_projectmat.cc"
#     warning "projection matrices are built explicitly"
#  else
#     include "rigid_fiber_project.cc"
#  endif

#else

void RigidFiber::buildProjection()   {}  //DIM == 1
void RigidFiber::makeProjection()    {}  //DIM == 1
void RigidFiber::destroyProjection() {}  //DIM == 1
void RigidFiber::allocateProjection(unsigned int) {}  //DIM == 1

void RigidFiber::projectForces(const real* X, real s, real* Y, real*) const
{
    real sum = X[0];
    for ( unsigned int ii = 1; ii < nbPoints(); ++ii )
        sum += X[ii];
    
    sum = s * sum / (real) nbPoints();
    for ( unsigned int ii = 0; ii < nbPoints(); ++ii )
        Y[ii] = sum;
}

void RigidFiber::computeTensions(const real*) {} //DIM == 1
void RigidFiber::makeProjectionDiff(const real*) const {} //DIM == 1
void RigidFiber::addProjectionDiff(const real*, real*) const {} //DIM == 1

#endif


//-----------------------------------------------------------------------
#pragma mark -

#if ( 0 )

/**
 only the upper terms are set
 */
void RigidFiber::addRigidityMatrix(Matrix & mat, const int offset, const int dim) const
{
    const real R = rfRigidity;
    for ( unsigned ii = 0; ii < nbPoints() - 2 ; ++ii )
    {
        mat(offset+dim* ii   , offset+dim* ii   ) -= R;
        mat(offset+dim* ii   , offset+dim*(ii+1)) += R * 2;
        mat(offset+dim* ii   , offset+dim*(ii+2)) -= R;
        mat(offset+dim*(ii+1), offset+dim*(ii+1)) -= R * 4;
        mat(offset+dim*(ii+1), offset+dim*(ii+2)) += R * 2;
        mat(offset+dim*(ii+2), offset+dim*(ii+2)) -= R;
    }
}

#else

/**
 Set elements of matrix `mat` corresponding to the elastic terms of the Fiber.
 The dimension of the matrix must be `dim * this->nbPoints()`
 Only the upper diagonal terms corresponding to the first subspace are set
 */
void RigidFiber::addRigidityMatrix(Matrix & mat, const int s, const int dim) const
{
    const real R = rfRigidity;
    int Z = nbPoints();
    if ( Z < 3 ) return;
    
    int ddd = 2*dim;
    int e = s + dim * ( Z - 2 );

    mat(s    , s    ) -= R;
    mat(s    , s+dim) += R * 2;
    mat(s    , s+ddd) -= R;
    
    mat(e    , e+dim) += R * 2;
    mat(e+dim, e+dim) -= R;

    if ( 3 < Z )
    {
        mat(s+dim, s+ddd) += R * 4;
        mat(s+dim, s+dim) -= R * 5;
        mat(e    , e    ) -= R * 5;
        mat(s+dim, s+dim+ddd) -= R;
    }
    else
    {
        mat(s+dim, s+dim) -= R * 4;
    }
    
    for ( int n = s+ddd; n < e ; n += dim )
    {
        mat(n, n    ) -= R * 6;
        mat(n, n+dim) += R * 4;
        mat(n, n+ddd) -= R;
    }
}

#endif

/**
 Set elements of matrix `mat` corresponding to the elastic terms of the Fiber.
 The array `mat` must be square of dimension `dim * this->nbPoints()`
 Only the upper diagonal terms corresponding to the first subspace are set
 */
void RigidFiber::addRigidityUpper(real * mat) const
{
    const real R = rfRigidity;
    int Z = nbPoints();
    if ( Z < 3 ) return;
    
    int LDD = Z*DIM;
    int ddd = 2*DIM;
    int e = DIM * ( Z - 2 );
    
    mat[0                  ] -= R;
    mat[        LDD*DIM    ] += R * 2;
    mat[        LDD*ddd    ] -= R;
    
    mat[e     + LDD*(e+DIM)] += R * 2;
    mat[e+DIM + LDD*(e+DIM)] -= R;
    
    if ( 3 < Z )
    {
        mat[DIM + LDD*ddd      ] += R * 4;
        mat[DIM + LDD*DIM      ] -= R * 5;
        mat[e   + LDD*e        ] -= R * 5;
        mat[DIM + LDD*(DIM+ddd)] -= R;
    }
    else
    {
        mat[DIM + LDD*DIM] -= R * 4;
    }
    
    for ( int n = ddd; n < e ; n += DIM )
    {
        mat[n + LDD*(n    )] -= R * 6;
        mat[n + LDD*(n+DIM)] += R * 4;
        mat[n + LDD*(n+ddd)] -= R;
    }
}

//------------------------------------------------------------------------------

/*
 This is the reference implementation
 */
inline void add_rigidity1(const unsigned nbt, const real* X, const real rigid, real* Y)
{
    for ( unsigned jj = 0; jj < nbt; ++jj )
    {
        real f = rigid * ( X[jj] - X[jj+DIM] - X[jj+DIM] + X[jj+DIM*2] );
        Y[jj      ] -=   f;
        Y[jj+DIM  ] += f+f;
        Y[jj+DIM*2] -=   f;
    }
}

/*
 In this version the loop is unrolled
 */
inline void add_rigidity2(const unsigned nbt, const real* X, const real rigid, real* Y)
{    
    real dx = X[DIM  ] - X[0];
    real dy = X[DIM+1] - X[1];
#if ( DIM > 2 )
    real dz = X[DIM+2] - X[2];
#endif
    
    real * yv = Y;
    const real*const end = X + nbt + DIM;
    
    const real* xv = X+DIM;
    while ( xv < end )
    {
        real d0 = xv[DIM] - xv[0];
        real f0 = rigid * ( d0 - dx );
        dx = d0;
        yv[0    ] -=    f0;
        yv[DIM  ] += f0+f0;
        yv[DIM*2] -=    f0;
        ++yv;
        ++xv;
        
        real d1 = xv[DIM] - xv[0];
        real f1 = rigid * ( d1 - dy );
        dy = d1;
        yv[0    ] -=    f1;
        yv[DIM  ] += f1+f1;
        yv[DIM*2] -=    f1;
        ++yv;
        ++xv;
        
#if ( DIM > 2 )
        real d2 = xv[DIM] - xv[0];
        real f2 = rigid * ( d2 - dz );
        dz = d2;
        Y[0    ] -=    f2;
        Y[DIM  ] += f2*f2;
        Y[DIM*2] -=    f2;
        ++yv;
        ++xv;
#endif
    }
}


/*
 In this version the loop is unrolled, pointers are used
 and further optimization are made by replacing
 ( a0 -2*a1 + a2 ) by (a2-a1)-(a1-a0).
 */
inline void add_rigidity3(const unsigned nbt, const real* X, const real rigid, real* Y)
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



#if defined __SSE3__ && ( DIM == 2 )

#include <pmmintrin.h>

/*
 In this version the loop is unrolled, pointers are used
 and further optimization are made by calculating
 (a2-a1)-(a1-a0) instead of ( a0 -2*a1 + a2 ).
 
 Fast version with SSE 128bit vector-arithmetics
 */

#include <emmintrin.h>
typedef __m128d vec2;
#define SSE(x) _mm_##x##_pd

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

#endif


void add_rigidityF(const unsigned nbt, const real* X, const real rigid, real* Y)
{
    real const* E = X + nbt + DIM;
    
    for ( unsigned d = 0; d < DIM; ++d )
    {
        Y[d]         += rigid * ( X[d+DIM] + X[d+DIM] - X[d] - X[d+2*DIM] );
        Y[nbt+DIM+d] += rigid * ( E[d-DIM] + E[d-DIM] - E[d] - E[d-2*DIM] );
    }
    
    if ( nbt == DIM )
    {
        for ( unsigned d = 0; d < DIM; ++d )
            Y[d+DIM] += 2 * rigid * ( X[d] - X[d+DIM] - X[d+DIM] + X[d+DIM*2] );
    }
    else
    {
        for ( unsigned ii = DIM*2; ii < nbt; ++ii )
            Y[ii] += rigid * ( - X[ii-DIM*2] - 6*X[ii] + 4*( X[ii-DIM] + X[ii+DIM] ) - X[ii+DIM*2] );
        
        for ( unsigned d = 0; d < DIM; ++d )
        {
            Y[d+DIM] += rigid * ( X[d] + X[d] - 5*X[d+DIM] + 4*X[d+DIM*2] - X[d+DIM*3] );
            Y[nbt+d] += rigid * ( E[d] + E[d] - 5*E[d-DIM] + 4*E[d-DIM*2] - E[d-DIM*3] );
        }
    }
}

//------------------------------------------------------------------------------

/**
 calculate the second-differential of points,
 scale by the rigidity term, and add to vector Y
*/
void RigidFiber::addRigidity(const real* X, real* Y) const
{
    int nbt = nbPoints() - 2;
    if ( nbt > 0 )
    {
        /// the manual SSE3 code is only valid with double precision and if DIM == 2:
#if (DIM==2) && defined(__SSE3__) && !defined(REAL_IS_FLOAT)
        #warning "Using SSE3 implementation"
        add_rigiditySSE(DIM*nbt, X, rfRigidity, Y);
#else
        add_rigidity3(DIM*nbt, X, rfRigidity, Y);
#endif
    
#if ( 0 )
        if ( prop->loop )
            loopRigidity(X,Y);
#endif
    }
}



/**
 Add rigidity terms between the last and first points, to loop the fiber onto itself.
 Done with Serge DMITRIEFF, 2015
 */
void RigidFiber::loopRigidity(const real* X, real* Y) const
{
    const unsigned lpx = DIM * lastPoint();
    if ( nbPoints() < 4 )
        return;
    
#if ( 1 )
    /*
     Better method done with Serge DMITRIEFF:
     link first and last point in the same way as all other points
     this makes the entire fiber mechanically symmetric and all points equivalent
    */
    for ( int d = 0; d < DIM; ++ d )
    {
        // 3-point rigidity term on [last, 0, 1]
        real f = rfRigidity * ( X[d+lpx] - 2*X[d] + X[d+DIM] );
        Y[d+lpx] -=   f;
        Y[d    ] += f+f;
        Y[d+DIM] -=   f;
        
        // 3-point rigidity term on [last-1, last, 0]
        real g = rfRigidity * ( X[d+lpx-DIM] - 2*X[d+lpx] + X[d] );
        Y[d+lpx-DIM] -=   g;
        Y[d+lpx    ] += g+g;
        Y[d        ] -=   g;
    }
    
#elif ( 1 )
    
    // Make links to overlap first and last points
    for ( int d = 0; d < DIM; ++ d )
    {
        // 3-point rigidity term on [last-1, 0, 1]
        real f = rfRigidity * ( X[d+lpx-DIM] - 2*X[d] + X[d+DIM] );
        Y[d+lpx-DIM] -=   f;
        Y[d        ] += f+f;
        Y[d    +DIM] -=   f;
        
        // 3-point rigidity term on [last-1, last, DIM]
        real g = rfRigidity * ( X[d+lpx-DIM] - 2*X[d+lpx] + X[d+DIM] );
        Y[d+lpx-DIM] -=   g;
        Y[d+lpx    ] += g+g;
        Y[d    +DIM] -=   g;
    }
    
#else
    
    for ( int d = 0; d < DIM; ++ d )
    {
        // 4-point rigidity term between [last-1, last] and [0 and 1]
        real f = rfRigidity * ( X[d+lpx-DIM] - X[d+lpx] - X[d] + X[d+DIM] );
        Y[d+lpx-DIM] -= f;
        Y[d+lpx    ] += f;
        Y[d        ] += f;
        Y[d    +DIM] -= f;
        
        // link last and first points, with an arbitrary (high) stiffness
        real g = 1000 * rfRigidity * ( X[d+lpx] - X[d] );
        Y[d+lpx] += g;
        Y[d    ] -= g;
    }

#endif
}


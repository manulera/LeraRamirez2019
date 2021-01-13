// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matsparsesym5.h"
#include "cblas.h"
#include "smath.h"

#include <iomanip>
#include <sstream>

#if defined(__SSE3__) &&  !defined(REAL_IS_FLOAT)
#   define MATRIX5_USES_INTEL_SSE3
#endif

//------------------------------------------------------------------------------
MatrixSparseSymmetric5::MatrixSparseSymmetric5()
: mxSize(size_)
{
    mxSize      = 0;
    mxAllocated = 0;
    
    col     = 0;
    colSize = 0;
    colMax  = 0;
    diag    = 0;
    
#ifdef MATRIX5_OPTIMIZED_MULTIPLY
    nCol    = 0;
    nmax    = 0;
    ija     = 0;
    sa      = 0;
#endif
}


//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::allocate( const unsigned int alc )
{
    if ( alc > mxAllocated )
    {
        Element ** col_new     = new Element*[alc];
        unsigned int * colSize_new  = new unsigned int[alc];
        unsigned int * colMax_new  = new unsigned int[alc];
        real    *  diag_new    = new real[alc];
        
        unsigned int ii = 0;
        if ( col )
        {
            for ( ; ii < mxAllocated; ++ii )
            {
                col_new[ii]     = col[ii];
                colSize_new[ii] = colSize[ii];
                colMax_new[ii]  = colMax[ii];
                diag_new[ii]    = diag[ii];
            }
            delete[] col;
            delete[] colSize;
            delete[] colMax;
            delete[] diag;
        }
        
        for ( ; ii < alc; ++ii )
        {
            col_new[ii]     = 0;
            colSize_new[ii] = 0;
            colMax_new[ii]  = 0;
            diag_new[ii]    = 0;
        }
        
        col       = col_new;
        colSize   = colSize_new;
        colMax    = colMax_new;
        diag      = diag_new;
        mxAllocated = alc;
        
#ifdef MATRIX5_OPTIMIZED_MULTIPLY
        if ( nCol ) delete[] nCol;
        nCol = new index_type[mxAllocated+1];
        for ( unsigned int jj = 0; jj <= mxSize; ++jj )
            nCol[jj] = jj;
#endif
    }
}

//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::deallocate()
{
    if ( col )
    {
        for ( unsigned int ii = 0; ii < mxAllocated; ++ii )
        {
            if ( col[ii] )
            {
                delete[] col[ii];
            }
        }
        delete[] col;       col     = 0;
        delete[] colSize;   colSize = 0;
        delete[] colMax;    colMax  = 0;
        delete[] diag;      diag    = 0;
    }
    mxAllocated = 0;
}

//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::allocateColumn( const index_type jj, unsigned int sz )
{
    assert_true( jj < mxSize );
    assert_true( sz > 0 );
    //printf("new S-COL %i %i\n", jj, sz );
    
    if ( sz > colMax[jj] )
    {
        const unsigned chunk = 4;
        sz = ( sz + chunk - 1 ) & ~( chunk -1 );
 
        Element * col_new = new Element[sz];
        
        if ( col[jj] )
        {
            
            //copy what is there
            for ( unsigned int ii = 0; ii < colMax[jj]; ++ii )
                col_new[ii] = col[jj][ii];
            delete[] col[jj];
            
        }
        col[jj]    = col_new;
        colMax[jj] = sz;
    }
}

//------------------------------------------------------------------------------
real& MatrixSparseSymmetric5::operator()( index_type ii, index_type jj )
{
    //this allocate the position if necessary
    assert_true( ii < mxSize );
    assert_true( jj < mxSize );
    
    //check for diagonal element
    if ( ii == jj )
        return diag[ii];
    
    // swap to get ii <= jj (address upper triangle)
    if ( jj < ii )
        std::swap(ii, jj);
    
    //check if the column is empty:
    if ( colSize[jj] == 0 )
    {
        allocateColumn( jj, 1 );
        colSize[jj]     = 1;
        col[jj][0].reset(ii);
        return col[jj][0].val;
    }
    
    //we keep the Elements ordered in the column by increasing indices
    Element * e    = col[jj];
    Element * last = col[jj] + colSize[jj];
    while ( e < last  &&  e->line < ii )
        ++e;
    if ( e < last  &&  e->line == ii )
        return e->val;
    
    //we will have to create/allocate a new Element, before the current one
    if ( colMax[jj] <= colSize[jj] )
    {
        index_type indx = e - col[jj];
        allocateColumn( jj, colSize[jj]+1 );
        assert_true( indx < colMax[jj] );
        e    = col[jj] + indx;
        last = col[jj] + colSize[jj];
    }
    
    if ( e == last )
    {
        //add the requested term at the end of the list:
        assert_true( colMax[jj] > colSize[jj] );
        e->reset(ii);
        ++colSize[jj];
        assert_true( e == col[jj] + colSize[jj] -1 );
        return e->val;
    }
    else {
        Element * p;
        //shift the list by one spot
        do {
            p = last - 1;
            last->line = p->line;
            last->val  = p->val;
            last = p;
        } while ( e < last );
        e->reset(ii);
        ++colSize[jj];
        return e->val;
    }
}

//------------------------------------------------------------------------------
real* MatrixSparseSymmetric5::addr( index_type ii, index_type jj ) const
{
    if ( ii == jj )
        return diag+ii;
    
    // swap to get ii <= jj (address upper triangle)
    if ( jj < ii )
        std::swap(ii, jj);
    
    for ( unsigned int kk = 0; kk < colSize[jj]; ++kk )
        if ( col[jj][kk].line == ii )
            return &( col[jj][kk].val );
    return 0;
}


//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::makeZero()
{
    for ( unsigned int ii = 0; ii < mxSize; ++ii )
    {
        colSize[ii] = 0;
        diag[ii] = 0;
    }
}


//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::scale( const real a )
{
    for ( unsigned int ii = 0; ii < mxSize; ++ii )
    {
        diag[ii] *= a;
        for ( unsigned int jj = 0; jj < colSize[ii]; ++jj )
            col[ii][jj].val *= a;
    }
}


void MatrixSparseSymmetric5::addTriangularBlock(real* mat, const unsigned ldd,
                                                const index_type si, const unsigned nb,
                                                const unsigned dim) const
{
    index_type up = si + nb;
    assert_true( up <= mxSize );

    for ( index_type jj = 0; jj < nb; ++jj )
    {
        mat[dim*(jj+ldd*jj)] += diag[jj+si];
        for ( unsigned int kk = 0; kk < colSize[jj+si]; ++kk )
        {
            index_type ii = col[jj+si][kk].line;
            if ( si <= ii && ii < up )
            {
                ii -= si;
                assert_true( ii <= jj );
                mat[dim*(ii+ldd*jj)] += col[jj+si][kk].val;
                //printf("Sp %4i %4i % .4f\n", ii, jj, a );
            }
        }
    }
}


//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::addDiagonalBlock(real* mat, unsigned ldd, index_type si, unsigned nb) const
{
    assert_true( si + nb <= mxSize );
    
    for ( index_type jj = 0; jj < nb; ++jj )
    {
        mat[jj+ldd*jj] += diag[jj+si];
        for ( unsigned int kk = 0 ; kk < colSize[jj+si] ; ++kk )
        {
            index_type ii = col[jj+si][kk].line;
            if ( si <= ii )
            {
                ii -= si;
                if ( ii < nb )
                {
                    assert_true( ii <= jj );
                    mat[ii+ldd*jj] += col[jj+si][kk].val;
                    if ( ii != jj )
                        mat[jj+ldd*ii] += col[jj+si][kk].val;
                    //printf("Sp %4i %4i % .4f\n", ii, jj, a );
                }
            }
        }
    }
}


//------------------------------------------------------------------------------
int MatrixSparseSymmetric5::bad() const
{
    if ( mxSize <= 0 ) return 1;
    for ( unsigned int jj = 0; jj < mxSize; ++jj )
    {
        for ( unsigned int kk = 0 ; kk < colSize[jj] ; ++kk )
        {
            if ( col[jj][kk].line >= mxSize ) return 2;
            if ( col[jj][kk].line <= jj )   return 3;
        }
    }
    return 0;
}


//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::printSparse(std::ostream& os) const
{
    os.precision(8);
    for ( unsigned int jj = 0; jj < mxSize; ++jj )
    {
        os << jj << " " << jj << " ";
        os << diag[jj] << std::endl;
        for ( unsigned int kk = 0 ; kk < colSize[jj] ; ++kk )
        {
            os << col[jj][kk].line << " " << jj << " ";
            os << col[jj][kk].val << std::endl;
        }
    }
}


//------------------------------------------------------------------------------
bool MatrixSparseSymmetric5::nonZero() const
{
    //check for any non-zero sparse term:
    for ( unsigned int jj = 0; jj < mxSize; ++jj )
    {
        if ( diag[jj] )
            return true;
        for ( unsigned int kk = 0 ; kk < colSize[jj] ; ++kk )
            if ( col[jj][kk].val )
                return true;
    }
    //if here, the matrix is empty
    return false;
}

//------------------------------------------------------------------------------
unsigned int MatrixSparseSymmetric5::nbElements() const
{
    //all allocated elements are counted, even if zero
    unsigned int cnt = 0;
    for ( unsigned int jj = 0; jj < mxSize; ++jj )
        cnt += colSize[jj];
    return cnt;
}

//------------------------------------------------------------------------------
std::string MatrixSparseSymmetric5::what() const
{
    std::ostringstream msg;
#ifdef MATRIX_USES_INTEL_SSE3
    msg << "mSS5i (" << nbElements() << ")";
#else
    msg << "mSS5 (" << nbElements() << ")";
#endif
    return msg.str();
}


#ifndef MATRIX5_OPTIMIZED_MULTIPLY

//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::prepareForMultiply()
{
}


//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::vecMulAdd( const real* X, real* Y ) const
{
    for ( unsigned int jj = 0; jj < mxSize; ++jj )
    {
        Y[jj] += diag[jj] * X[jj];
        for ( int kk = 0 ; kk < colSize[jj] ; ++kk )
        {
            const index_type ii = col[jj][kk].line;
            const real a = col[jj][kk].val;
            Y[ii] += a * X[jj];
            if ( ii != jj )
                Y[jj] += a * X[ii];
        }
    }
}

//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::vecMulAddIso2D( const real* X, real* Y ) const
{
    for ( unsigned int jj = 0; jj < mxSize; ++jj )
    {
        const index_type Djj = 2 * jj;
        Y[Djj  ] += diag[jj] * X[Djj  ];
        Y[Djj+1] += diag[jj] * X[Djj+1];
        for ( unsigned int kk = 0 ; kk < colSize[jj] ; ++kk )
        {
            const index_type Dii = 2 * col[jj][kk].line;
            const real  a = col[jj][kk].val;
            Y[Dii  ] += a * X[Djj  ];
            Y[Dii+1] += a * X[Djj+1];
            if ( Dii != Djj )
            {
                Y[Djj  ] += a * X[Dii  ];
                Y[Djj+1] += a * X[Dii+1];
            }
        }
    }
}

//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::vecMulAddIso3D( const real* X, real* Y ) const
{
    for ( unsigned int jj = 0; jj < mxSize; ++jj )
    {
        const index_type Djj = 3 * jj;
        Y[Djj  ] += diag[jj] * X[Djj  ];
        Y[Djj+1] += diag[jj] * X[Djj+1];
        Y[Djj+2] += diag[jj] * X[Djj+2];
        for ( unsigned int kk = 0 ; kk < colSize[jj] ; ++kk )
        {
            const index_type Dii = 3 * col[jj][kk].line;
            const real  a = col[jj][kk].val;
            Y[Dii  ] += a * X[Djj  ];
            Y[Dii+1] += a * X[Djj+1];
            Y[Dii+2] += a * X[Djj+2];
            if ( Dii != Djj )
            {
                Y[Djj  ] += a * X[Dii  ];
                Y[Djj+1] += a * X[Dii+1];
                Y[Djj+2] += a * X[Dii+2];
            }
        }
    }
}

#else

//==============================================================================
//                           MATRIX_OPTIMIZED_MULTIPLY
//==============================================================================
//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::prepareForMultiply()
{
    
    //set array nCol[]
    int inx = mxSize;
    nCol[mxSize] = mxSize;
    for ( unsigned int jj = mxSize-1; jj > 0; --jj )
    {
        if ( colSize[jj] > 0  || diag[jj] != 0 )  inx = jj;
        nCol[jj] = inx;
    }
    nCol[0] = 0;
    
    //count number of non-zero elements, including diagonal
    index_type nbe = 1+mxSize;
    for ( unsigned int jj = 0; jj < mxSize; ++jj )
        nbe += colSize[jj];
    
    //we use a sparse matrix indexed storage from Numerical Recipes
    //allocate arrays ( ija = indices, sa = values )
    if ( nbe > nmax )
    {
        if ( ija )  delete[] ija;
        if ( sa )   delete[] sa;
        
        nmax  = nbe + mxSize;
        ija   = new index_type[nmax];
        sa    = new real[nmax];
    }
    
    //create the sparse representation, zero-based indices
    ija[0] = mxSize+1;
    index_type kk = mxSize;
    for ( unsigned int jj = 0; jj < mxSize; ++jj )
    {
        // diagonal term first:
        sa[jj]  = diag[jj];
        for ( unsigned int cc = 0; cc < colSize[jj]; ++cc )
        {
            ++kk;
            if ( kk >= nbe ) ABORT_NOW("internal out of range error");
            sa[kk]  = col[jj][cc].val;
            ija[kk] = col[jj][cc].line;
            //printf("col %i: %i %f\n", jj, ija[kk], sa[kk]);
        }
        ija[jj+1] = kk+1;
    }
    assert_true( kk+1 == nbe );
}


//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::vecMulAdd( const real* X, real* Y ) const
{
    for ( index_type jj = nCol[0]; jj < mxSize; jj = nCol[jj+1] )
    {
        real X0 = X[jj];
        real Y0 = Y[jj] + sa[jj] * X0;
        const index_type end = ija[jj+1];
        for ( index_type kk = ija[jj]; kk < end; ++kk )
        {
            index_type ii = ija[kk];
            real a =  sa[kk];
            assert_true( jj != ii );
            Y0    += a * X[ii];
            Y[ii] += a * X0;
        }
        Y[jj] = Y0;
    }
}

//------------------------------------------------------------------------------

#ifdef MATRIX_USES_INTEL_SSE3

#include <pmmintrin.h>
#warning "Using SSE3 implementation"

void MatrixSparseSymmetric5::vecMulAddIso2D( const real* X, real* Y ) const
{
    __m128d x, y, a, t;
    for ( index_type jj = nCol[0]; jj < mxSize; jj = nCol[jj+1] )
    {
        const index_type Djj = 2 * jj;
        x = _mm_load_pd(X+Djj);
        a = _mm_loaddup_pd(sa+jj);
        y = _mm_add_pd(_mm_load_pd(Y+Djj), _mm_mul_pd(a, x));
        const index_type end = ija[jj+1];
        for ( index_type kk = ija[jj]; kk < end; ++kk )
        {
            const index_type Dii = 2 * ija[kk];
            a = _mm_loaddup_pd(sa+kk);
            t = _mm_add_pd(_mm_load_pd(Y+Dii), _mm_mul_pd(x, a));
            y = _mm_add_pd(y, _mm_mul_pd(_mm_load_pd(X+Dii), a));
            _mm_store_pd(Y+Dii, t);
        }
        _mm_store_pd(Y+Djj, y);
    }
}
#else
void MatrixSparseSymmetric5::vecMulAddIso2D( const real* X, real* Y ) const
{
    for ( index_type jj = nCol[0]; jj < mxSize; jj = nCol[jj+1] )
    {
        index_type Djj = 2 * jj;
        real X0 = X[Djj  ];
        real X1 = X[Djj+1];
        real Y0 = Y[Djj  ] + sa[jj] * X0;
        real Y1 = Y[Djj+1] + sa[jj] * X1;
        const index_type end = ija[jj+1];
        for ( index_type kk = ija[jj]; kk < end; ++kk )
        {
            index_type Dii = 2 * ija[kk];
            assert_true( Djj != Dii );
            real  a = sa[kk];
            Y0       += a * X[Dii  ];
            Y1       += a * X[Dii+1];
            Y[Dii  ] += a * X0;
            Y[Dii+1] += a * X1;
        }
        Y[Djj  ] = Y0;
        Y[Djj+1] = Y1;
    }
}
#endif

//------------------------------------------------------------------------------
void MatrixSparseSymmetric5::vecMulAddIso3D( const real* X, real* Y ) const
{
    for ( index_type jj = nCol[0]; jj < mxSize; jj = nCol[jj+1] )
    {
        index_type Djj = 3 * jj;
        real X0 = X[Djj  ];
        real X1 = X[Djj+1];
        real X2 = X[Djj+2];
        real Y0 = Y[Djj  ] + sa[jj] * X0;
        real Y1 = Y[Djj+1] + sa[jj] * X1;
        real Y2 = Y[Djj+2] + sa[jj] * X2;
        const index_type end = ija[jj+1];
        for ( index_type kk = ija[jj]; kk < end; ++kk )
        {
            index_type Dii = 3 * ija[kk];
            assert_true( Djj != Dii );
            real  a = sa[kk];
            Y0       += a * X[Dii  ];
            Y1       += a * X[Dii+1];
            Y2       += a * X[Dii+2];
            Y[Dii  ] += a * X0;
            Y[Dii+1] += a * X1;
            Y[Dii+2] += a * X2;
        }
        Y[Djj  ] = Y0;
        Y[Djj+1] = Y1;
        Y[Djj+2] = Y2;
    }
}

#endif


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matsparsesym1.h"
#include "cblas.h"
#include "smath.h"

#include <iomanip>
#include <sstream>

#if defined(__SSE3__) &&  !defined(REAL_IS_FLOAT)
#   define MATRIX1_USES_INTEL_SIMD
#endif


MatrixSparseSymmetric1::MatrixSparseSymmetric1()
{
    allocated_ = 0;
    col_       = 0;
    col_size_  = 0;
    col_max_   = 0;
    
#ifdef MATRIX1_OPTIMIZED_MULTIPLY
    nmax_      = 0;
    ija_       = 0;
    sa_        = 0;
    col_next_  = new index_type[1];
    col_next_[0] = 0;
#endif
}



void MatrixSparseSymmetric1::allocate(unsigned alc)
{
    if ( alc > allocated_ )
    {
        /*
         'chunk' can be increased to gain performance:
         more memory will be used, but reallocation will be less frequent
         */
        const unsigned chunk = 64;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );

        //fprintf(stderr, "MSS1 allocate matrix %u\n", alc);
        Element ** col_new      = new Element*[alc];
        unsigned * col_size_new = new unsigned[alc];
        unsigned * col_max_new  = new unsigned[alc];
        
        index_type ii = 0;
        if ( col_ )
        {
            for ( ; ii < allocated_; ++ii )
            {
                col_new[ii]      = col_[ii];
                col_size_new[ii] = col_size_[ii];
                col_max_new[ii]  = col_max_[ii];
            }
            delete[] col_;
            delete[] col_size_;
            delete[] col_max_;
        }
        
        col_       = col_new;
        col_size_  = col_size_new;
        col_max_   = col_max_new;
        allocated_ = alc;

        for ( ; ii < alc; ++ii )
        {
            col_[ii]      = 0;
            col_size_[ii] = 0;
            col_max_[ii]  = 0;
        }
        
#ifdef MATRIX1_OPTIMIZED_MULTIPLY
        if ( col_next_ )
            delete[] col_next_;
        col_next_ = new index_type[allocated_+1];
#endif
    }
}


void MatrixSparseSymmetric1::deallocate()
{
    if ( col_ )
    {
        for ( index_type ii = 0; ii < allocated_; ++ii )
        {
            if ( col_[ii] )
                delete[] col_[ii];
        }
        delete[] col_;       col_      = 0;
        delete[] col_size_;  col_size_ = 0;
        delete[] col_max_;   col_max_  = 0;
#ifdef MATRIX1_OPTIMIZED_MULTIPLY
        delete[] col_next_;  col_next_ = 0;
#endif
    }
    allocated_ = 0;
}


/// copy `cnt` elements from `src` to `dst`
void copy(unsigned cnt, MatrixSparseSymmetric1::Element * src, MatrixSparseSymmetric1::Element * dst)
{
    for ( unsigned ii = 0; ii < cnt; ++ii )
        dst[ii] = src[ii];
}

/// move `cnt` elements to next index, starting at vec[0]
void shift(unsigned cnt, MatrixSparseSymmetric1::Element * vec)
{
    for ( unsigned ii = cnt; ii > 0; --ii )
        vec[ii] = vec[ii-1];
}


void MatrixSparseSymmetric1::allocateColumn(const index_type jj, unsigned alc)
{
    assert_true( jj < size_ );
    
    if ( alc > col_max_[jj] )
    {
        //fprintf(stderr, "MSS1 allocate column %i size %u\n", jj, alc);
        const unsigned chunk = 16;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );
        Element * col_new = new Element[alc];
        
        if ( col_[jj] )
        {
            //copy over previous column elements
            copy(col_size_[jj], col_[jj], col_new);
            
            //release old memory
            delete[] col_[jj];
        }
        col_[jj]     = col_new;
        col_max_[jj] = alc;
    }
}


real& MatrixSparseSymmetric1::diagonal(index_type ix)
{
    assert_true( ix < size_ );
    
    Element * col;
    
    if ( col_size_[ix] == 0 )
    {
        allocateColumn(ix, 1);
        col = col_[ix];
        //diagonal term always first:
        col->reset(ix);
        col_size_[ix] = 1;
    }
    else
    {
        col = col_[ix];
        assert_true( col->inx == ix );
    }
    
    return col->val;
}

/**
 This allocate to be able to hold the matrix element if necessary
*/
real& MatrixSparseSymmetric1::operator()(index_type ii, index_type jj)
{
    assert_true( ii < size_ );
    assert_true( jj < size_ );
    //fprintf(stderr, "MSS1( %6i %6i )\n", ii, jj);
    
    Element * col;

    // swap to get ii > jj (address lower triangle)
    if ( ii < jj )
        std::swap(ii, jj);
    else if ( ii == jj )
    {
        // return diagonal element
        if ( col_size_[jj] <= 0 )
        {
            allocateColumn(jj, 1);
            col = col_[jj];
            // put diagonal term always first:
            col->reset(jj);
            col_size_[jj] = 1;
            return col->val;
        }
        else
        {
            col = col_[jj];
            assert_true( col->inx == jj );
        }
        return col->val;
    }
 
    //check if the column is empty:
    if ( col_size_[jj] < 2 )
    {
        allocateColumn(jj, 2);
        col = col_[jj];
        if ( col_size_[jj] == 0 )
        {
            // put diagonal term always first:
            col->reset(jj);
        }
        //add the requested term:
        col[1].reset(ii);
        col_size_[jj] = 2;
        return col[1].val;
    }
    
    col = col_[jj];
    Element * e = col + 1;
    Element * lst = col + col_size_[jj] - 1;
    
    //search, knowing that elements are kept ordered in the column:
    while ( e->inx < ii )
    {
        if ( ++e > lst )
        {
            // add element last
            unsigned n = col_size_[jj];
            if ( n >= col_max_[jj] )
            {
                allocateColumn(jj, n+1);
                col = col_[jj];
            }
            ++col_size_[jj];
            col[n].reset(ii);
            return col[n].val;
        }
    }
    
    if ( e->inx == ii )
        return e->val;
    
    unsigned n = e - col;

    assert_true( col[n].inx > ii );
    
    // allocate space for new Element if necessary:
    if ( col_size_[jj] >= col_max_[jj] )
    {
        const unsigned chunk = 16;
        unsigned alc = ( col_max_[jj]+ 1 + chunk - 1 ) & ~( chunk -1 );
        Element * col_new = new Element[alc];
        copy(n, col, col_new);
        col_new[n].reset(ii);
        copy(col_size_[jj]-n, col+n, col_new+n+1);
        ++col_size_[jj];
        col_[jj] = col_new;
        col_max_[jj] = alc;
        delete[] col;
        return col_new[n].val;
    }
    
    shift(col_size_[jj]-n, col+n);
    ++col_size_[jj];
    
    // add the requested term
    col[n].reset(ii);
    
    //printColumn(jj);
    return col[n].val;
}


real* MatrixSparseSymmetric1::addr(index_type ii, index_type jj) const
{
    // swap to get ii <= jj (address lower triangle)
    if ( ii < jj )
        std::swap(ii, jj);
    
    for ( unsigned kk = 0; kk < col_size_[jj]; ++kk )
        if ( col_[jj][kk].inx == ii )
            return &( col_[jj][kk].val );
    
    return 0;
}


//------------------------------------------------------------------------------
#pragma mark -

void MatrixSparseSymmetric1::makeZero()
{
    for ( index_type jj = 0; jj < size_; ++jj )
        col_size_[jj] = 0;
}


bool MatrixSparseSymmetric1::nonZero() const
{
    //check for any non-zero sparse term:
    for ( unsigned jj = 0; jj < size_; ++jj )
        for ( unsigned kk = 0 ; kk < col_size_[jj] ; ++kk )
            if ( col_[jj][kk].val )
                return true;
    
    //if here, the matrix is empty
    return false;
}


void MatrixSparseSymmetric1::scale(const real alpha)
{
    for ( index_type jj = 0; jj < size_; ++jj )
        for ( index_type n = 0; n < col_size_[jj]; ++n )
            col_[jj][n].val *= alpha;
}


void MatrixSparseSymmetric1::addTriangularBlock(real* mat, const unsigned ldd,
                                                const index_type si,
                                                const unsigned nb,
                                                const unsigned dim) const
{
    index_type up = si + nb;
    assert_true( up <= size_ );
    
    for ( index_type jj = si; jj < up; ++jj )
    {
        for ( unsigned n = 0; n < col_size_[jj]; ++n )
        {
            index_type ii = col_[jj][n].inx;
            // assuming lower triangle is stored:
            assert_true( ii >= jj );
            if ( ii < up )
            {
                //printf("MSS1 %4i %4i % .4f\n", ii, jj, a);
                mat[dim*( jj-si + ldd * (ii-si) )] += col_[jj][n].val;
            }
        }
    }
}


void MatrixSparseSymmetric1::addDiagonalBlock(real* mat, unsigned ldd,
                                              const index_type si,
                                              const unsigned nb) const
{
    index_type up = si + nb;
    assert_true( up <= size_ );
    
    for ( index_type jj = si; jj < up; ++jj )
    {
        for ( unsigned n = 0; n < col_size_[jj]; ++n )
        {
            index_type ii = col_[jj][n].inx;
            // assuming lower triangle is stored:
            assert_true( ii >= jj );
            if ( ii < up )
            {
                //printf("MSS1 %4i %4i % .4f\n", ii, jj, a);
                mat[jj-si+ldd*(ii-si)] += col_[jj][n].val;
                if ( jj != ii )
                    mat[ii-si+ldd*(jj-si)] += col_[jj][n].val;
            }
        }
    }
}



int MatrixSparseSymmetric1::bad() const
{
    if ( size_ <= 0 ) return 1;
    for ( unsigned jj = 0; jj < size_; ++jj )
    {
        for ( unsigned kk = 0 ; kk < col_size_[jj] ; ++kk )
        {
            if ( col_[jj][kk].inx >= size_ ) return 2;
            if ( col_[jj][kk].inx <= jj )   return 3;
        }
    }
    return 0;
}


unsigned MatrixSparseSymmetric1::nbElements() const
{
    //all allocated elements are counted, even if zero
    unsigned cnt = 0;
    for ( unsigned jj = 0; jj < size_; ++jj )
        cnt += col_size_[jj];
    return cnt;
}


std::string MatrixSparseSymmetric1::what() const
{
    std::ostringstream msg;
#ifdef MATRIX1_USES_INTEL_SIMD
    msg << "mSS1i (" << nbElements() << ")";
#else
    msg << "mSS1 (" << nbElements() << ")";
#endif
    return msg.str();
}


void MatrixSparseSymmetric1::printSparse(std::ostream & os) const
{
    std::streamsize p = os.precision();
    os.precision(8);
    for ( unsigned jj = 0; jj < size_; ++jj )
    {
        for ( unsigned n = 0 ; n < col_size_[jj] ; ++n )
        {
            os << col_[jj][n].inx << " " << jj << " ";
            os << col_[jj][n].val << "\n";
        }
    }
    os.precision(p);
}


void MatrixSparseSymmetric1::printColumns(std::ostream & os)
{
    os << "MSS1 size " << size_ << ":";
    for ( int jj = 0; jj < size_; ++jj )
    {
        os << "\n   " << jj << "   " << col_size_[jj];
#if MATRIX1_OPTIMIZED_MULTIPY
        os << " " << col_next_[jj];
#endif
    }
    os << std::endl;
}


void MatrixSparseSymmetric1::printColumn(std::ostream & os, const index_type jj)
{
    Element const* col = col_[jj];
    os << "MSS1 col " << jj << ":";
    for ( unsigned n = 0; n < col_size_[jj]; ++n )
    {
        os << "\n" << col[n].inx << " :";
        os << " " << col[n].val;
    }
    os << std::endl;
}

//------------------------------------------------------------------------------
#ifndef MATRIX1_OPTIMIZED_MULTIPLY

void MatrixSparseSymmetric1::prepareForMultiply()
{
}


void MatrixSparseSymmetric1::vecMulAdd(const real* X, real* Y) const
{
    for ( index_type jj = 0; jj < size_; ++jj )
    {
        if ( col_size_[jj] > 0 )
        {
            const real X0 = X[jj];
            assert_true( col_[jj][0].inx == jj );
            real Y0 = Y[jj] + col_[jj][0].val * X0;
            for ( unsigned kk = 1 ; kk < col_size_[jj] ; ++kk )
            {
                const index_type ii = col_[jj][kk].inx;
                const real a = col_[jj][kk].val;
                Y[ii] += a * X0;
                assert_true( ii != jj );
                Y0 += a * X[ii];
            }
            Y[jj] = Y0;
        }
    }
}


void MatrixSparseSymmetric1::vecMulAddIso2D(const real* X, real* Y) const
{
    for ( index_type jj = 0; jj < size_; ++jj )
    {
        if ( col_size_[jj] > 0 )
        {
            const index_type Djj = 2 * jj;
            const real X0 = X[Djj  ];
            const real X1 = X[Djj+1];
            assert_true( col_[jj][0].inx == jj );
            real Y0 = Y[Djj  ] + col_[jj][0].val * X0;
            real Y1 = Y[Djj+1] + col_[jj][0].val * X1;
            for ( unsigned n = 1 ; n < col_size_[jj] ; ++n )
            {
                const index_type Dii = 2 * col_[jj][n].inx;
                const real  a = col_[jj][n].val;
                Y[Dii  ] += a * X0;
                Y[Dii+1] += a * X1;
                assert_true( Dii != Djj );
                Y0 += a * X[Dii  ];
                Y1 += a * X[Dii+1];
            }
            Y[Djj  ] = Y0;
            Y[Djj+1] = Y1;
        }
    }
}


void MatrixSparseSymmetric1::vecMulAddIso3D(const real* X, real* Y) const
{
    for ( index_type jj = 0; jj < size_; ++jj )
    {
        if ( col_size_[jj] > 0 )
        {
            const index_type Djj = 3 * jj;
            const real X0 = X[Djj  ];
            const real X1 = X[Djj+1];
            const real X2 = X[Djj+2];
            assert_true( col_[jj][0].inx == jj );
            real Y0 = Y[Djj  ] + col_[jj][0].val * X0;
            real Y1 = Y[Djj+1] + col_[jj][0].val * X1;
            real Y2 = Y[Djj+2] + col_[jj][0].val * X2;
            for ( unsigned kk = 1 ; kk < col_size_[jj] ; ++kk )
            {
                const index_type Dii = 3 * col_[jj][kk].inx;
                const real  a = col_[jj][kk].val;
                Y[Dii  ] += a * X0;
                Y[Dii+1] += a * X1;
                Y[Dii+2] += a * X2;
                assert_true( Dii != Djj );
                Y0 += a * X[Dii  ];
                Y1 += a * X[Dii+1];
                Y2 += a * X[Dii+2];
            }
            Y[Djj  ] = Y0;
            Y[Djj+1] = Y1;
            Y[Djj+2] = Y2;
        }
    }
}

#else  // MATRIX1_OPTIMIZED_MULTIPLY defined below


void MatrixSparseSymmetric1::setNextColumn()
{
    col_next_[size_] = size_;

    if ( size_ > 0 )
    {
        unsigned inx = size_;
        unsigned nxt = size_;
        while ( --inx > 0 )
        {
            if ( col_size_[inx] > 0 )
                nxt = inx;
            col_next_[inx] = nxt;
        }
        if ( col_size_[0] > 0 )
            col_next_[0] = 0;
        else
            col_next_[0] = nxt;
    }
}


void MatrixSparseSymmetric1::prepareForMultiply()
{
    assert_true( size_ <= allocated_ );
    
    setNextColumn();
    
    //count number of non-zero elements, including diagonal
    unsigned nbe = 1;
    for ( unsigned jj = 0; jj < size_; ++jj )
    {
        if ( col_size_[jj] > 0 )
            nbe += col_size_[jj];
        else
            nbe ++;
    }
    
    //allocate classical sparse matrix storage (Numerical Recipes)
    if ( nbe > nmax_ )
    {
        if ( ija_ ) delete[] ija_;
        if ( sa_ )  delete[] sa_;
        
        nmax_  = nbe + size_;
        ija_   = new index_type[nmax_];
        sa_    = new real[nmax_];
    }
    
    //create the sparse representation, described in numerical-recipe
    //indices start at zero, unlike in numerical recipe
    ija_[0] = size_+1;
    index_type kk = size_;
    for ( unsigned jj = 0; jj < size_; ++jj )
    {
        if ( col_size_[jj] > 0 )
        {
            // diagonal term first:
            assert_true( col_[jj][0].inx == jj );
            sa_[jj] = col_[jj][0].val;
            // other elements:
            for ( unsigned cc = 1; cc < col_size_[jj]; ++cc )
            {
                ++kk;
                assert_true( kk < nbe );
                sa_[kk]  = col_[jj][cc].val;
                ija_[kk] = col_[jj][cc].inx;
            }
        }
        else {
            sa_[jj] = 0;
        }
        ija_[jj+1] = kk+1;
    }
    if ( kk+1 != nbe ) ABORT_NOW("internal error");
    
    //printSparse(std::clog);
    //printSparseArray(std::clog);
}


void MatrixSparseSymmetric1::printSparseArray(std::ostream& os) const
{
    unsigned end = ija_[size_];

    os << "ija ";
    for ( index_type n = 0; n < end; ++n )
        os << " " << std::setw(6) << ija_[n];
    os << "\n";

    std::streamsize p = os.precision();
    os.precision(2);
    os << "sa  ";
    for ( index_type n = 0; n < end; ++n )
        os << " " << std::setw(6) << sa_[n];
    os << "\n";
    os.precision(p);
}


//------------------------------------------------------------------------------
#pragma mark -

void MatrixSparseSymmetric1::vecMulAdd(const real* X, real* Y) const
{
    for ( index_type jj = col_next_[0]; jj < size_; jj = col_next_[jj+1] )
    {
        assert_true( col_size_[jj] > 0 );
        assert_true( col_[jj][0].inx == jj );
        real X0 = X[jj];
        real Y0 = Y[jj] + sa_[jj] * X0;
        const index_type end = ija_[jj+1];
        for ( index_type kk = ija_[jj]; kk < end; ++kk )
        {
            real a = sa_[kk];
            index_type ii = ija_[kk];
            Y[ii] += a * X0;
            Y0    += a * X[ii];
        }
        Y[jj] = Y0;
    }
}


#ifdef MATRIX1_USES_INTEL_SIMD

#include <immintrin.h>

#if defined __AVX2__ && defined __ICC
#   define USE_SIMD_FMA
#endif

#ifdef USE_SIMD_FMA
#   warning "Using SSE3 implementation with Fused Multiply Add instructions"
#else
#   warning "Using SSE3 implementation"
#endif

#define SSE(x) _mm_##x##_pd
typedef __m128d vec2;


void MatrixSparseSymmetric1::vecMulAddIso2D(const real* X, real* Y) const
{
    for ( index_type jj = col_next_[0]; jj < size_; jj = col_next_[jj+1] )
    {
        assert_true( col_size_[jj] > 0 );
        assert_true( col_[jj][0].inx == jj );
        vec2 xx = SSE(load)(X+2*jj);
        vec2 aa = SSE(load1)(sa_+jj);
#ifdef USE_SIMD_FMA
        vec2 yy = SSE(fmadd)(aa, xx, SSE(load)(Y+2*jj));
#else
        vec2 yy = SSE(add)(SSE(load)(Y+2*jj), SSE(mul)(aa, xx));
#endif
        //vec2 yy = aa * xx + SSE(load)(Y+2*jj);

        const index_type end = ija_[jj+1];
        for ( index_type kk = ija_[jj]; kk < end; ++kk )
        {
            aa = SSE(load1)(sa_+kk);
#ifdef USE_SIMD_FMA
            vec2 tt = SSE(fmadd)(xx, aa, SSE(load)(Y+2*ija_[kk]));
            yy = SSE(fmadd)(SSE(load)(X+2*ija_[kk]), aa, yy);
#else
            vec2 tt = SSE(add)(SSE(load)(Y+2*ija_[kk]), SSE(mul)(xx, aa));
            yy = SSE(add)(yy, SSE(mul)(SSE(load)(X+2*ija_[kk]), aa));
#endif
            SSE(store)(Y+2*ija_[kk], tt);
        }
        SSE(store)(Y+2*jj, yy);
    }
}

#else

void MatrixSparseSymmetric1::vecMulAddIso2D(const real* X, real* Y) const
{    
    for ( index_type jj = col_next_[0]; jj < size_; jj = col_next_[jj+1] )
    {
        assert_true( col_size_[jj] > 0 );
        assert_true( col_[jj][0].inx == jj );
        index_type Djj = 2 * jj;
        real X0 = X[Djj  ];
        real X1 = X[Djj+1];
        real Y0 = Y[Djj  ] + sa_[jj] * X0;
        real Y1 = Y[Djj+1] + sa_[jj] * X1;
        const index_type end = ija_[jj+1];
        for ( index_type kk = ija_[jj]; kk < end; ++kk )
        {
            index_type Dii = 2 * ija_[kk];
            assert_true( Djj != Dii );
            real a = sa_[kk];
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


void MatrixSparseSymmetric1::vecMulAddIso3D(const real* X, real* Y) const
{
    //for ( index_type jj = 0; jj < size_; ++jj )
    for ( index_type jj = col_next_[0]; jj < size_; jj = col_next_[jj+1] )
    {
        assert_true( col_size_[jj] > 0 );
        assert_true( col_[jj][0].inx == jj );
        index_type Djj = 3 * jj;
        real X0 = X[Djj  ];
        real X1 = X[Djj+1];
        real X2 = X[Djj+2];
        real Y0 = Y[Djj  ] + sa_[jj] * X0;
        real Y1 = Y[Djj+1] + sa_[jj] * X1;
        real Y2 = Y[Djj+2] + sa_[jj] * X2;
        const index_type next = ija_[jj+1];
        for ( index_type kk = ija_[jj]; kk < next; ++kk )
        {
            index_type Dii = 3 * ija_[kk];
            assert_true( Djj != Dii );
            real a = sa_[kk];
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


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


#include "matsparsesymblk.h"
#include <sstream>


MatrixSparseSymmetricBlock::MatrixSparseSymmetricBlock()
{
    allocated_ = 0;
    col_       = 0;
    col_max_   = 0;
    col_size_  = 0;
    
    col_next_  = new index_type[1];
    col_next_[0] = 0;
}



void MatrixSparseSymmetricBlock::allocate(unsigned alc)
{
    if ( alc % BSZ )
        ABORT_NOW("the matrix dimension in not a multiple of its block size");
    
    if ( alc > allocated_ )
    {
        /*
         'chunk' can be increased to gain performance:
          more memory will be used, but reallocation will be less frequent
        */
        const unsigned chunk = 64;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );

        //fprintf(stderr, "MSSB allocate matrix %u\n", alc);
        Element ** col_new      = new Element*[alc+1];
        unsigned * col_max_new  = new unsigned[alc+1];
        unsigned * col_size_new = new unsigned[alc+1];
       
        index_type ii = 0;
        if ( col_ )
        {
            for ( ; ii <= allocated_; ++ii )
            {
                col_new[ii]      = col_[ii];
                col_max_new[ii]  = col_max_[ii];
                col_size_new[ii] = col_size_[ii];
            }
            delete[] col_;
            delete[] col_max_;
            delete[] col_size_;
       }
        
        col_       = col_new;
        col_max_   = col_max_new;
        col_size_  = col_size_new;
        allocated_ = alc;

        for ( ; ii <= alc; ++ii )
        {
            col_[ii]      = 0;
            col_max_[ii]  = 0;
            col_size_[ii] = 0;
        }

        delete[] col_next_;
        col_next_ = new index_type[allocated_+1];
    }
}


void MatrixSparseSymmetricBlock::deallocate()
{
    if ( col_ )
    {
        for ( index_type ii = 0; ii <= allocated_; ++ii )
        {
            if ( col_[ii] )
                delete[] col_[ii];
        }
        delete[] col_;       col_      = 0;
        delete[] col_max_;   col_max_  = 0;
        delete[] col_size_;  col_size_ = 0;
        delete[] col_next_;  col_next_ = 0;
    }
    allocated_ = 0;
}


/// copy `cnt` elements from `src` to `dst`
void copy(unsigned cnt, MatrixSparseSymmetricBlock::Element * src, MatrixSparseSymmetricBlock::Element * dst)
{
    for ( unsigned ii = 0; ii < cnt; ++ii )
        dst[ii] = src[ii];
}


void MatrixSparseSymmetricBlock::allocateColumn(const index_type jj, unsigned alc)
{
    assert_true( jj < size_ );
    
    if ( alc > col_max_[jj] )
    {
        //fprintf(stderr, "MSS1 allocate column %i size %u\n", jj, alc);
        /*
         'chunk' can be increased, to gain performance:
         more memory will be used, but reallocation will be less frequent
         */
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


/**
 This allocate to be able to hold the matrix element if necessary
*/
MatrixSparseSymmetricBlock::Element& MatrixSparseSymmetricBlock::block(index_type ii, index_type jj)
{
    assert_true( ii < size_ );
    assert_true( jj < size_ );
    assert_true( ii % BSZ == 0 );
    assert_true( jj % BSZ == 0 );
    
    //fprintf(stderr, "MSS( %6i %6i )\n", ii, jj);
    
    // swap to get ii > jj (address lower triangle)
    if ( ii < jj )
        std::swap(ii, jj);
    
    Element * col = col_[jj];
    
    if ( col_size_[jj] > 0 )
    {
        Element * e = col;
        Element * lst = col + col_size_[jj] - 1;
        
        //check all elements in the column:
        while ( e <= lst )
        {
            if ( e->inx == ii )
                return *e;
            ++e;
        }
    }
    else
    {
        allocateColumn(jj, 2);
        col = col_[jj];
        // put diagonal term always first:
        col->reset(jj);
        if ( ii == jj )
        {
            col_size_[jj] = 1;
            return col[0];
        }
        //add the requested term:
        col[1].reset(ii);
        col_size_[jj] = 2;
        return col[1];
    }
    
    // add the requested term at the end:
    unsigned n = col_size_[jj];
    
    // allocate space for new Element if necessary:
    if ( n >= col_max_[jj] )
    {
        allocateColumn(jj, n+1);
        col = col_[jj];
    }
    
    assert_true( n < col_max_[jj] );
    col[n].reset(ii);
    ++col_size_[jj];
    
    //printColumn(jj);
    return col[n];
}


real& MatrixSparseSymmetricBlock::operator()(index_type ii, index_type jj)
{
    index_type iir = ii % BSZ;
    index_type jjr = jj % BSZ;

    // swap to address lower triangle
    if ( ii < jj )
    {
        Element & e = block(jj-jjr, ii-iir);
        return *e.addr(jjr, iir);
    }
    else
    {
        Element & e = block(ii-iir, jj-jjr);
        return *e.addr(iir, jjr);
    }
}


real* MatrixSparseSymmetricBlock::addr(index_type iid, index_type jjd) const
{
    index_type ii, jj, iir, jjr;
    // swap to address lower triangle
    if ( iid < jjd )
    {
        iir = jjd % BSZ;
        ii  = jjd - iir;
        jjr = iid % BSZ;
        jj  = iid - jjr;
    }
    else
    {
        iir = iid % BSZ;
        ii  = iid - iir;
        jjr = jjd % BSZ;
        jj  = jjd - jjr;
    }
    
    for ( unsigned n = 0; n < col_size_[jj]; ++n )
        if ( col_[jj][n].inx == ii )
            return col_[jj][n].addr(iir, jjr);
    
    return 0;
}


//------------------------------------------------------------------------------
#pragma mark -

void MatrixSparseSymmetricBlock::makeZero()
{
    for ( index_type jj = 0; jj < size_; ++jj )
        col_size_[jj] = 0;

    // add a terminal value for AVX implementation of mulVec:
    if ( size_ > 0 )
    {
        col_size_[size_] = 0;
        if ( !col_[size_] )
        {
            col_[size_] = new Element[1];
            col_max_[size_] = 1;
        }
        col_[size_][0].reset(0);
    }
}


bool MatrixSparseSymmetricBlock::nonZero() const
{
    //check for any non-zero sparse term:
    for ( unsigned jj = 0; jj < size_; ++jj )
        for ( unsigned n = 0 ; n < col_size_[jj] ; ++n )
        {
            real * M = col_[jj][n].val;
            for ( int u = 0; u < BSZ*BSZ; ++u )
                if ( M[u] )
                    return true;
        }
    
    //if here, the matrix is empty
    return false;
}


void MatrixSparseSymmetricBlock::scale(const real alpha)
{
    for ( index_type jj = 0; jj < size_; ++jj )
        for ( index_type n = 0; n < col_size_[jj]; ++n )
            col_[jj][n].scale(alpha);
}


/// add all elements of block 'S' to 'M'
void addBlockF(real * M, unsigned ldd, real const* S)
{
    M[0    ] += S[0];
#if ( BSZ == 2 )
    M[1    ] += S[1];
    M[  ldd] += S[2];
    M[1+ldd] += S[3];
#elif ( BSZ == 3 )
    M[1      ] += S[1];
    M[2      ] += S[2];
    M[  ldd  ] += S[3];
    M[1+ldd  ] += S[4];
    M[2+ldd  ] += S[5];
    M[  ldd*2] += S[6];
    M[1+ldd*2] += S[7];
    M[2+ldd*2] += S[8];
#endif
}

/// add lower elements of block 'S' to upper triangle of 'M'
void addBlockU(real * M, unsigned ldd, real const* S)
{
    M[0    ] += S[0];
#if ( BSZ == 2 )
    M[  ldd] += S[1];
    M[1+ldd] += S[3];
    assert_true( S[2] == 0 );
#elif ( BSZ == 3 )
    M[  ldd  ] += S[1];
    M[1+ldd  ] += S[4];
    M[  ldd*2] += S[2];
    M[1+ldd*2] += S[5];
    M[2+ldd*2] += S[8];
    assert_true( S[3] == 0 );
    assert_true( S[6] == 0 );
    assert_true( S[7] == 0 );
#endif
}

/// add lower elements of block 'S' to both upper and lower triangles of 'M'
void addBlockS(real * M, unsigned ldd, real const* S)
{
    M[0    ] += S[0];
#if ( BSZ == 2 )
    M[1    ] += S[1];
    M[  ldd] += S[1];
    M[1+ldd] += S[3];
    assert_true( S[2] == 0 );
#elif ( BSZ == 3 )
    M[1      ] += S[1];
    M[2      ] += S[2];
    M[  ldd  ] += S[1];
    M[1+ldd  ] += S[4];
    M[2+ldd  ] += S[5];
    M[  ldd*2] += S[2];
    M[1+ldd*2] += S[5];
    M[2+ldd*2] += S[8];
    assert_true( S[3] == 0 );
    assert_true( S[6] == 0 );
    assert_true( S[7] == 0 );
#endif
}

/// add all elements of block 'S' to 'M', with transposition
void addBlockT(real * M, unsigned ldd, real const* S)
{
    M[0    ] += S[0];
#if ( BSZ == 2 )
    M[1    ] += S[2];
    M[  ldd] += S[1];
    M[1+ldd] += S[3];
#elif ( BSZ == 3 )
    M[1      ] += S[3];
    M[2      ] += S[6];
    M[  ldd  ] += S[1];
    M[1+ldd  ] += S[4];
    M[2+ldd  ] += S[7];
    M[  ldd*2] += S[2];
    M[1+ldd*2] += S[5];
    M[2+ldd*2] += S[8];
#endif
}


void MatrixSparseSymmetricBlock::addTriangularBlock(real* mat, const unsigned ldd,
                                                const index_type si,
                                                const unsigned nb,
                                                const unsigned dim) const
{
    if ( si % BSZ )  ABORT_NOW("incompatible index");
    if ( nb % BSZ )  ABORT_NOW("incompatible size");

    index_type up = si + nb;
    index_type off = si + ldd * si;
    assert_true( up <= size_ );
    
    for ( index_type jj = si; jj < up; ++jj )
    {
        if ( col_size_[jj] > 0 )
        {
            assert_true(col_[jj][0].inx == jj);
            addBlockU(mat + ( jj + ldd*jj ) - off, ldd, col_[jj][0].val);
            for ( unsigned n = 1; n < col_size_[jj]; ++n )
            {
                index_type ii = col_[jj][n].inx;
                // assuming lower triangle is stored:
                assert_true( ii > jj );
                if ( ii < up )
                {
                    //printf("MSSB %4i %4i % .4f\n", ii, jj, a);
                    addBlockT(mat + ( jj + ldd*ii ) - off, ldd, col_[jj][n].val);
                }
            }
        }
    }
}


void MatrixSparseSymmetricBlock::addDiagonalBlock(real* mat, unsigned ldd,
                                              const index_type si,
                                              const unsigned nb) const
{
    if ( si % BSZ )  ABORT_NOW("incompatible index");
    if ( nb % BSZ )  ABORT_NOW("incompatible size");
    
    index_type up = si + nb;
    index_type off = si + ldd * si;
    assert_true( up <= size_ );
    
    for ( index_type jj = si; jj < up; ++jj )
    {
        if ( col_size_[jj] > 0 )
        {
            assert_true(col_[jj][0].inx == jj);
            addBlockS(mat+( jj + ldd*jj )-off, ldd, col_[jj][0].val);
            for ( unsigned n = 1; n < col_size_[jj]; ++n )
            {
                index_type ii = col_[jj][n].inx;
                // assuming lower triangle is stored:
                assert_true( ii > jj );
                if ( ii < up )
                {
                    //printf("MSSB %4i %4i % .4f\n", ii, jj, a);
                    addBlockF(mat + ( ii + ldd*jj ) - off, ldd, col_[jj][n].val);
                    addBlockT(mat + ( jj + ldd*ii ) - off, ldd, col_[jj][n].val);
                }
            }
        }
    }
}


int MatrixSparseSymmetricBlock::bad() const
{
    if ( size_ <= 0 ) return 1;
    for ( unsigned jj = 0; jj < size_; ++jj )
    {
        for ( unsigned n = 0 ; n < col_size_[jj] ; ++n )
        {
            if ( col_[jj][n].inx >= size_ ) return 2;
            if ( col_[jj][n].inx <= jj )    return 3;
        }
    }
    return 0;
}


unsigned MatrixSparseSymmetricBlock::nbElements() const
{
    //all allocated elements are counted, even if zero
    unsigned cnt = 0;
    for ( unsigned jj = 0; jj < size_; jj+=BSZ )
        cnt += col_size_[jj];
    return cnt;
}


//------------------------------------------------------------------------------
#pragma mark -


std::string MatrixSparseSymmetricBlock::what() const
{
    std::ostringstream msg;
#if defined(__AVX2__) &&  !defined(REAL_IS_FLOAT)
    msg << "mSSB* (" << BSZ*BSZ << "x" << nbElements() << ")";
#elif defined(__SSE3__) &&  !defined(REAL_IS_FLOAT)
    msg << "mSSB+ (" << BSZ*BSZ << "x" << nbElements() << ")";
#else
    msg << "mSSB (" << BSZ*BSZ << "x" << nbElements() << ")";
#endif
    return msg.str();
}



void MatrixSparseSymmetricBlock::printSparse(std::ostream & os) const
{
    std::streamsize p = os.precision();
    os.precision(8);
    for ( unsigned jj = 0; jj < size_; ++jj )
    {
        for ( unsigned n = 0 ; n < col_size_[jj] ; ++n )
        {
            index_type ii = col_[jj][n].inx;
            real * M = col_[jj][n].val;
            if ( ii == jj )
            {
                for ( int x = 0; x < BSZ; ++x )
                for ( int y = x; y < BSZ; ++y )
                    os << ii+y << " " << jj+x << " " << M[y+BSZ*x] << std::endl;
            }
            else
            {
                for ( int x = 0; x < BSZ; ++x )
                for ( int y = 0; y < BSZ; ++y )
                    os << ii+y << " " << jj+x << " " << M[y+BSZ*x] << std::endl;
            }
        }
    }
    os.precision(p);
}


void MatrixSparseSymmetricBlock::printColumns(std::ostream & os)
{
    os << "MSS3 size " << size_ << ":";
    for ( int jj = 0; jj < size_; ++jj )
    {
        os << "\n   " << jj << "   " << col_size_[jj];
        os << " " << col_next_[jj];
    }
    os << std::endl;
}


void MatrixSparseSymmetricBlock::printColumn(std::ostream & os, const index_type jj)
{
    Element const* col = col_[jj];
    os << "MSS3 col " << jj << ":";
    for ( unsigned n = 0; n < col_size_[jj]; ++n )
    {
        os << "\n" << col[n].inx << " :";
        for ( unsigned u = 0; u < BSZ*BSZ; ++u )
            os << " " << col[n].val[u];
    }
    os << std::endl;
}


//------------------------------------------------------------------------------
#pragma mark -

void MatrixSparseSymmetricBlock::prepareForMultiply()
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
    
    //printColumns(std::clog);
}




#if ( BSZ == 1 )

void MatrixSparseSymmetricBlock::vecMulAdd(const real* X, real* Y) const
{
    for ( index_type jj = col_next_[0]; jj < size_; jj = col_next_[jj+1] )
    {
        assert_true( jj < size_ );
        assert_true( col_size_[jj] > 0 );
        const real X0 = X[jj];
        assert_true( col_[jj][0].inx == jj );
        real * M = col_[jj][0].val;
        real Y0 = Y[jj] + M[0] * X0;
        for ( unsigned n = 1 ; n < col_size_[jj] ; ++n )
        {
            const index_type ii = col_[jj][n].inx;
            assert_true( ii < size_ );
            assert_true( ii != jj );
            real * M = col_[jj][n].val;
            Y[ii] += M[0] * X0;
            Y0 += M[0] * X[ii];
        }
        Y[jj] = Y0;
    }
}

#elif ( BSZ == 2 )

#if defined(__AVX2__) &&  !defined(REAL_IS_FLOAT)

#include <immintrin.h>

#if defined __AVX2__ && defined __ICC
#   define USE_SIMD_FMA
#endif

#ifdef USE_SIMD_FMA
#   warning "Using AVX2 implementation with Fused Multiply Add instructions"
#else
#   warning "Using AVX2 implementation"
#endif

typedef __m128d vec2;
#define SSE(x) _mm_##x##_pd
typedef __m256d vec4;
#define AVX(x) _mm256_##x##_pd

void MatrixSparseSymmetricBlock::vecMulAdd(const real* X, real* Y) const
{
    vec2 d01, d23;
    
    assert_true( col_next_[0] <= size_ );
    {
        //load next 2x2 matrix diagonal element:
        real * m = col_[col_next_[0]][0].val;
        d01 = SSE(load)(m  );
        d23 = SSE(load)(m+2);
    }
    
    //std::clog << "MatrixSparseSymmetricBlock using Intel SIMD\n";
    index_type nn;
    for ( index_type jj = col_next_[0]; jj < size_; jj = nn )
    {
        nn = col_next_[jj+1];
        assert_true( col_size_[jj] > 0 );
        vec2 x0 = SSE(load1)(X+jj  );
        vec2 x1 = SSE(load1)(X+jj+1);

        assert_true( col_[jj][0].inx == jj );
        
        vec2 m01 = d01;
        vec2 m23 = d23;
        vec2 m13 = SSE(unpackhi)(m01, m23);

        // load first 2x2 matrix element into 2 vectors:
        real * n = col_[jj][1].val;
        vec2 n01 = SSE(load)(n  );
        vec2 n23 = SSE(load)(n+2);

        assert_true( nn <= size_ );
        {
            assert_true( col_max_[nn] > 0 );
            // load diagonal 2x2 matrix element from next column:
            real * nd = col_[nn][0].val;
            d01 = SSE(load)(nd  );
            d23 = SSE(load)(nd+2);
        }
        
        //multiply with the symmetrized block:
        //real Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1;
        //real Y1 = Y[jj+1] + M[1] * X0 + M[3] * X1;
#ifdef USE_SIMD_FMA
        //std::clog << "MatrixSparseSymmetricBlock using FMA instructions\n";
        vec2 uu = SSE(fmadd)(m01, x0, SSE(load)(Y+jj));
        vec2 yy = SSE(fmadd)(m13, x1, uu);
#else
        vec2 uu = SSE(add)(SSE(mul)(m01, x0), SSE(load)(Y+jj));
        vec2 yy = SSE(add)(SSE(mul)(m13, x1), uu);
#endif
        
        const unsigned sup = col_size_[jj]-1;
        for ( unsigned n = 1 ; n < sup ; ++n )
        {
            const index_type ii = col_[jj][n].inx;
            assert_true(ii != jj);
            
            //
            vec2 m01 = n01;
            vec2 m23 = n23;
            {
                // load next 2x2 matrix element into 2 vectors:
                real * m = col_[jj][n+1].val;
                n01 = SSE(loadu)(m  );
                n23 = SSE(loadu)(m+2);
            }
            
            // multiply with the full block:
            //Y[ii  ] += M[0] * X0 + M[2] * X1;
            //Y[ii+1] += M[1] * X0 + M[3] * X1;
#ifdef USE_SIMD_FMA
            vec2 uu = SSE(fmadd)(m01, x0, SSE(load)(Y+ii));
            vec2 tt = SSE(fmadd)(m23, x1, uu);
#else
            vec2 uu = SSE(add)(SSE(mul)(m01, x0), SSE(load)(Y+ii));
            vec2 tt = SSE(add)(SSE(mul)(m23, x1), uu);
#endif
            SSE(store)(Y+ii, tt);
            
            // multiply with the transposed block:
            //Y0 += M[0] * X[ii] + M[1] * X[ii+1];
            //Y1 += M[2] * X[ii] + M[3] * X[ii+1];
#ifdef USE_SIMD_FMA
            vec2 m02 = SSE(unpacklo)(m01, m23);
            vec2 m13 = SSE(unpackhi)(m01, m23);
            vec2 vv = SSE(fmadd)(m02, SSE(load1)(X+ii), yy);
            yy = SSE(fmadd)(m13, SSE(load1)(X+ii+1), vv);
#else
            vec2 xx = SSE(load)(X+ii);
            vec2 vv = SSE(hadd)(SSE(mul)(m01, xx), SSE(mul)(m23, xx));
            yy = SSE(add)(vv, yy);
#endif
        }
        
        if ( sup > 0 )
        {
            const index_type ii = col_[jj][sup].inx;
            // multiply with the full block:
            //Y[ii  ] += M[0] * X0 + M[2] * X1;
            //Y[ii+1] += M[1] * X0 + M[3] * X1;
#ifdef USE_SIMD_FMA
            uu = SSE(fmadd)(n01, x0, SSE(load)(Y+ii));
            vec2 tt = SSE(fmadd)(n23, x1, uu);
#else
            uu = SSE(add)(SSE(mul)(n01, x0), SSE(load)(Y+ii));
            vec2 tt = SSE(add)(SSE(mul)(n23, x1), uu);
#endif
            SSE(store)(Y+ii, tt);
            
            // multiply with the transposed block:
            //Y0 += M[0] * X[ii] + M[1] * X[ii+1];
            //Y1 += M[2] * X[ii] + M[3] * X[ii+1];
#ifdef USE_SIMD_FMA
            vec2 m02 = SSE(unpacklo)(n01, n23);
            m13 = SSE(unpackhi)(n01, n23);
            vec2 vv = SSE(fmadd)(m02, SSE(load1)(X+ii), yy);
            yy = SSE(fmadd)(m13, SSE(load1)(X+ii+1), vv);
#else
            vec2 xx = SSE(load)(X+ii);
            vec2 vv = SSE(hadd)(SSE(mul)(n01, xx), SSE(mul)(n23, xx));
            yy = SSE(add)(vv, yy);
#endif
        }
        //Y[jj  ] = Y0;
        //Y[jj+1] = Y1;
        SSE(store)(Y+jj, yy);
    }
}

#elif defined(__SSE3__) &&  !defined(REAL_IS_FLOAT)

#include <immintrin.h>
#warning "Using SSE3 implementation"

typedef __m128d vec2;
#define SSE(x) _mm_##x##_pd

void MatrixSparseSymmetricBlock::vecMulAdd(const real* X, real* Y) const
{
    //std::clog << "MatrixSparseSymmetricBlock using Intel SIMD\n";
    for ( index_type jj = col_next_[0]; jj < size_; jj = col_next_[jj+1] )
    {
        assert_true( col_size_[jj] > 0 );
        vec2 x0 = SSE(load1)(X+jj  );
        vec2 x1 = SSE(load1)(X+jj+1);
        
        assert_true( col_[jj][0].inx == jj );
        
        // load 2x2 matrix element into 2 vectors:
        real * m = col_[jj][0].val;
        vec2 m01 = SSE(load)(m  );
        vec2 m23 = SSE(load)(m+2);
        vec2 m13 = SSE(unpackhi)(m01, m23);
        
        //multiply with the symmetrized block:
        //real Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1;
        //real Y1 = Y[jj+1] + M[1] * X0 + M[3] * X1;
        vec2 uu = SSE(add)(SSE(mul)(m01, x0), SSE(load)(Y+jj));
        vec2 yy = SSE(add)(SSE(mul)(m13, x1), uu);
        
        for ( unsigned n = 1 ; n < col_size_[jj]; ++n )
        {
            const index_type ii = col_[jj][n].inx;
            assert_true(ii != jj);
            
            // load 2x2 matrix element into 2 vectors:
            real * m = col_[jj][n].val;
            m01 = SSE(loadu)(m  );
            m23 = SSE(loadu)(m+2);
            
            // multiply with the full block:
            //Y[ii  ] += M[0] * X0 + M[2] * X1;
            //Y[ii+1] += M[1] * X0 + M[3] * X1;
            vec2 uu = SSE(add)(SSE(mul)(m01, x0), SSE(load)(Y+ii));
            vec2 tt = SSE(add)(SSE(mul)(m23, x1), uu);
            SSE(store)(Y+ii, tt);
            
            // multiply with the transposed block:
            //Y0 += M[0] * X[ii] + M[1] * X[ii+1];
            //Y1 += M[2] * X[ii] + M[3] * X[ii+1];
            vec2 xx = SSE(load)(X+ii);
            vec2 vv = SSE(hadd)(SSE(mul)(m01, xx), SSE(mul)(m23, xx));
            yy = SSE(add)(vv, yy);
        }
        //Y[jj  ] = Y0;
        //Y[jj+1] = Y1;
        SSE(store)(Y+jj, yy);
    }
}

#else

void MatrixSparseSymmetricBlock::vecMulAdd(const real* X, real* Y) const
{
    for ( index_type jj = col_next_[0]; jj < size_; jj = col_next_[jj+1] )
    {
        assert_true( col_size_[jj] > 0 );
        assert_true( jj < size_ );
        const real X0 = X[jj  ];
        const real X1 = X[jj+1];
        assert_true( col_[jj][0].inx == jj );
        real * M = col_[jj][0].val;
        assert_true( M[2] == 0 );
        // multiply with the symmetrized block:
        real Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1;
        real Y1 = Y[jj+1] + M[1] * X0 + M[3] * X1;
        for ( unsigned n = 1 ; n < col_size_[jj] ; ++n )
        {
            const index_type ii = col_[jj][n].inx;
            assert_true( ii < size_ );
            assert_true( ii != jj );
            real * M = col_[jj][n].val;
            // multiply with the full block:
            Y[ii  ] += M[0] * X0 + M[2] * X1;
            Y[ii+1] += M[1] * X0 + M[3] * X1;
            // multiply with the transposed block:
            Y0 += M[0] * X[ii] + M[1] * X[ii+1];
            Y1 += M[2] * X[ii] + M[3] * X[ii+1];
        }
        Y[jj  ] = Y0;
        Y[jj+1] = Y1;
    }
}

#endif

#else  // BSZ == 3

///@todo Optimize with INTEL AVX2
void MatrixSparseSymmetricBlock::vecMulAdd(const real* X, real* Y) const
{
    for ( index_type jj = col_next_[0]; jj < size_; jj = col_next_[jj+1] )
    {
        assert_true( jj < size_ );
        assert_true( col_size_[jj] > 0 );
        const real X0 = X[jj  ];
        const real X1 = X[jj+1];
        const real X2 = X[jj+2];
        assert_true( col_[jj][0].inx == jj );
        real * M = col_[jj][0].val;
        assert_true( M[3] == 0 );
        assert_true( M[6] == 0 );
        assert_true( M[7] == 0 );
        // multiply with the symmetrized block:
        real Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1 + M[2] * X2;
        real Y1 = Y[jj+1] + M[1] * X0 + M[4] * X1 + M[5] * X2;
        real Y2 = Y[jj+2] + M[2] * X0 + M[5] * X1 + M[8] * X2;
        
        for ( unsigned n = 1 ; n < col_size_[jj] ; ++n )
        {
            const index_type ii = col_[jj][n].inx;
            assert_true( ii < size_ );
            assert_true( ii != jj );
            real * M = col_[jj][n].val;
            // multiply with the full block:
            Y[ii  ] +=  M[0] * X0 + M[3] * X1 + M[6] * X2;
            Y[ii+1] +=  M[1] * X0 + M[4] * X1 + M[7] * X2;
            Y[ii+2] +=  M[2] * X0 + M[5] * X1 + M[8] * X2;
            // multiply with the transposed block:
            Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
            Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
            Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        }
        Y[jj  ] = Y0;
        Y[jj+1] = Y1;
        Y[jj+2] = Y2;
    }
}

#endif


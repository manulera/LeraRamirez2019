// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matsparsesym4.h"
#include "cblas.h"
#include "smath.h"

#include <iomanip>
#include <sstream>

#define MATRIX4_OPTIMIZED_MULTIPLY

#if defined(__SSE3__) &&  !defined(REAL_IS_FLOAT)
#   define MATRIX_USES_INTEL_SSE3
#endif


MatrixSparseSymmetric4::MatrixSparseSymmetric4()
: mxSize(size_)
{
    mxAllocated = 0;
    col        = 0;
    colSize    = 0;
    colMax     = 0;
    nCol       = new index_type[1];
    nCol[0]    = 0;
}


void MatrixSparseSymmetric4::allocate( const unsigned int alc )
{
    if ( alc > mxAllocated )
    {
        Element ** col_new         = new Element*[alc];
        unsigned int * colSize_new = new unsigned int[alc];
        unsigned int * colMax_new  = new unsigned int[alc];
        
        unsigned int ii = 0;
        if ( col )
        {
            for ( ; ii < mxAllocated; ++ii )
            {
                col_new[ii]      = col[ii];
                colSize_new[ii]  = colSize[ii];
                colMax_new[ii]   = colMax[ii];
            }
            delete[] col;
            delete[] colSize;
            delete[] colMax;
        }
        
        for ( ; ii < alc; ++ii )
        {
            col_new[ii]      = 0;
            colSize_new[ii]  = 0;
            colMax_new[ii]   = 0;
        }
        
        col         = col_new;
        colSize     = colSize_new;
        colMax      = colMax_new;
        mxAllocated = alc;

        if ( nCol )
            delete[] nCol;
        nCol = new index_type[mxAllocated+1];
    }
}


void MatrixSparseSymmetric4::deallocate()
{
    if ( col )
    {
        for ( unsigned int ii = 0; ii < mxAllocated; ++ii )
            if ( col[ii] )
            {
                delete[] col[ii];
            };
        delete[] col;       col     = 0;
        delete[] colSize;   colSize = 0;
        delete[] colMax;    colMax  = 0;
        delete[] nCol;      nCol    = 0;
    }
    mxAllocated = 0;
}


void MatrixSparseSymmetric4::allocateColumn( const index_type jj, unsigned int sz )
{
    assert_true( jj < mxSize );
    assert_true( sz > 0 );
    //printf("new S-COL %i %i\n", jj, sz );
    
    if ( sz > colMax[jj] )
    {
        const unsigned chunk = 4;
        sz = ( sz + chunk - 1 ) & ~( chunk -1 );

        Element * col_new  = new Element[sz];
        
        if ( col[jj] )
        {
            
            //copy what is there
            for ( unsigned int ii = 0; ii < colMax[jj]; ++ii )
            {
                col_new[ii] = col[jj][ii];
            }
            delete[] col[jj];
            
        }
        col[jj]  = col_new;
        colMax[jj] = sz;
    }
}


real& MatrixSparseSymmetric4::operator()( index_type ii, index_type jj )
{
    //this allocate the position if necessary
    assert_true( ii < mxSize );
    assert_true( jj < mxSize );
    
    Element * e, * last;
    
    // swap to get ii <= jj (address upper triangle)
    if ( jj < ii )
        std::swap(ii, jj);
    
    //check if the column is empty:
    if ( colSize[jj] == 0 )
    {
        allocateColumn( jj, 2 );
        e = col[jj];
        
        //always include the diagonal term first:
        e->indx = jj;
        e->val  = 0.;
        colSize[jj] = 1;
        
        if ( ii != jj )
        {
            //add the requested term:
            ++e;
            e->indx = ii;
            e->val  = 0.;
            colSize[jj] = 2;
        }
        return e->val;
    }
    
    e = col[jj];
    assert_true( e->indx == jj ); //the first term should be the diagonal
    
    //check if diagonal term is requested
    if ( ii == jj )
        return e->val;
    
    ///\todo optimize this search in MatrixSparse::operator()
    last = e + colSize[jj];
    ++e; //we have checked the first term already, which is the diagonal
    while ( e < last )
    {
        if ( e->indx == ii )
            return e->val;
        ++e;
    }
    
    //we will have to create/allocate a new Element
    if ( colMax[jj] <= colSize[jj] )
    {
        allocateColumn( jj, colSize[jj]+1 );
        e = col[jj] + colSize[jj];
    }
    
    assert_true( colMax[jj] > colSize[jj] );
    
    //add the requested term:
    e->indx = ii;
    e->val  = 0.;
    ++colSize[jj];
    return e->val;
}


real* MatrixSparseSymmetric4::addr( index_type ii, index_type jj ) const
{
    // swap to get ii <= jj (address upper triangle)
    if ( jj < ii )
        std::swap(ii, jj);
    
    for ( unsigned int kk = 0; kk < colSize[jj]; ++kk )
        if ( col[jj][kk].indx == ii )
            return &( col[jj][kk].val );
    return 0;
}


void MatrixSparseSymmetric4::makeZero()
{
    for ( index_type ii = 0; ii < mxSize; ++ii )
        colSize[ii] = 0;
}


void MatrixSparseSymmetric4::scale( real a )
{
    for ( index_type ii = 0; ii < mxSize; ++ii )
        for ( unsigned int jj = 0; jj < colSize[ii]; ++jj )
            col[ii][jj].val *= a;
}


void MatrixSparseSymmetric4::addTriangularBlock(real* mat, const unsigned ldd,
                                                const index_type si,
                                                const unsigned nb,
                                                const unsigned dim) const
{
    index_type up = si + nb;
    assert_true( up <= size_ );
    
    for ( index_type jj = si; jj < up; ++jj )
    {
        for ( unsigned n = 0; n < colSize[jj]; ++n )
        {
            index_type ii = col[jj][n].indx;
            if ( si <= ii && ii < up )
            {
                //printf("MSS1 %4i %4i % .4f\n", ii, jj, a);
                mat[dim*( ii-si + ldd * (jj-si) )] += col[jj][n].val;
            }
        }
    }
}


void MatrixSparseSymmetric4::addDiagonalBlock(real* mat, unsigned ldd,
                                              const index_type si,
                                              const unsigned nb) const
{
    index_type up = si + nb;
    assert_true( up <= size_ );
    
    for ( index_type jj = si; jj < up; ++jj )
    {
        for ( unsigned n = 0; n < colSize[jj]; ++n )
        {
            index_type ii = col[jj][n].indx;
            if ( si <= ii && ii < up )
            {
                //printf("MSS1 %4i %4i % .4f\n", ii, jj, a);
                mat[ii-si+ldd*(jj-si)] += col[jj][n].val;
                if ( jj != ii )
                    mat[jj-si+ldd*(ii-si)] += col[jj][n].val;
            }
        }
    }
}

//------------------------------------------------------------------------------

int MatrixSparseSymmetric4::bad() const
{
    if ( mxSize <= 0 ) return 1;
    for ( index_type jj = 0; jj < mxSize; ++jj )
    {
        for ( unsigned int kk = 0 ; kk < colSize[jj] ; ++kk )
        {
            if ( col[jj][kk].indx >= mxSize ) return 2;
            if ( col[jj][kk].indx <= jj ) return 3;
        }
    }
    return 0;
}


void MatrixSparseSymmetric4::printSparse(std::ostream& os) const
{
    os.precision(8);
    for ( index_type jj = 0; jj < mxSize; ++jj )
    {
        for ( index_type kk = 0 ; kk < colSize[jj] ; ++kk )
        {
            os << col[jj][kk].indx << " " << jj << " ";
            os << col[jj][kk].val << std::endl;
        }
    }
}


bool MatrixSparseSymmetric4::nonZero() const
{
    //check for any non-zero sparse term:
    for ( index_type jj = 0; jj < mxSize; ++jj )
        for ( index_type kk = 0 ; kk < colSize[jj] ; ++kk )
            if ( col[jj][kk].val )
                return true;
    
    //if here, the matrix is empty
    return false;
}


unsigned int MatrixSparseSymmetric4::nbElements() const
{
    //all allocated elements are counted, even if val == zero
    unsigned int cnt = 0;
    for ( index_type jj = 0; jj < mxSize; ++jj )
        cnt += colSize[jj];
    return cnt;
}


std::string MatrixSparseSymmetric4::what() const
{
    std::ostringstream msg;
#ifdef MATRIX_USES_INTEL_SSE3
    msg << "mSS4i (" << nbElements() << ")";
#else
    msg << "mSS4 (" << nbElements() << ")";
#endif
    return msg.str();
}


int MatrixSparseSymmetric4::compare( const void * a, const void * b )
{
    return ((Element*) a)->indx - ((Element*)b)->indx;
}


void MatrixSparseSymmetric4::prepareForMultiply()
{
    assert_true( mxSize <= mxAllocated );

    //update nCol[], a pointer to the next non-empty column:
    index_type inx = mxSize;
    nCol[mxSize] = mxSize;
    for ( int jj = mxSize-1; jj >= 0; --jj )
    {
        if ( colSize[jj] > 0 )
            inx = jj;
        nCol[jj] = inx;
    }
    
    //the second optimization is not so clear, in terms of CPU gain.
#if ( 0 )
    {
        //we reorder the columns to have ordered memory access
        for ( index_type jj = nCol[0]; jj < mxSize; jj = nCol[jj+1] )
        {
            if ( colSize[jj] > 1 )
            {
                //printf("\nColumn %i:\n", jj);
                //for ( int kk=0; kk<colSize[jj]; ++kk ) printf("%i:% .2f ", col[jj][kk].indx, col[jj][kk].val);
                
                //optimize by putting the diagonal term first
                //the diagonal term is usually not very far, so we do a simple search
                for ( unsigned int kk = 0; kk < colSize[jj]; ++kk )
                {
                    if ( col[jj][kk].indx == jj )
                    {
                        if ( kk > 0 )
                        {
                            Element tmp = col[jj][0];
                            col[jj][0]  = col[jj][kk];
                            col[jj][kk] = tmp;
                        }
                        break;
                    }
                }
                assert_true( sizeof(**col) == sizeof(Element) );
                qsort(col[jj]+1, colSize[jj]-1, sizeof(**col), &compare);
                
                //printf("\n");printf("\nColumn %i:\n", jj);
                //for ( int kk=0; kk<colSize[jj]; ++kk ) printf("%i:% .2f ", col[jj][kk].indx, col[jj][kk].val);
            }
        }
    }
#endif
}


//------------------------------------------------------------------------------
#ifndef MATRIX4_OPTIMIZED_MULTIPLY


void MatrixSparseSymmetric4::vecMulAdd( const real* X, real* Y ) const
{
    for ( index_type jj = 0; jj < mxSize; ++jj )
    {
        for ( unsigned int kk = 0 ; kk < colSize[jj] ; ++kk )
        {
            const int ii = col[jj][kk].line;
            const real a = col[jj][kk].val;
            Y[ii] += a * X[jj];
            if ( ii != jj )
                Y[jj] += a * X[ii];
        }
    }
}

void MatrixSparseSymmetric4::vecMulAddIso2D( const real* X, real* Y ) const
{
    for ( index_type int jj = 0; jj < mxSize; ++jj )
    {
        const index_type Djj = 2 * jj;
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

void MatrixSparseSymmetric4::vecMulAddIso3D( const real* X, real* Y ) const
{
    for ( index_type jj = 0; jj < mxSize; ++jj )
    {
        const index_type Djj = 3 * jj;
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

#else  // MATRIX4_OPTIMIZED_MULTIPLY


void MatrixSparseSymmetric4::vecMulAdd( const real* X, real* Y ) const
{
    const Element * elmt, * last;
    
    //we only visit the non-empty columns
    for ( index_type jj = nCol[0]; jj < mxSize; jj = nCol[jj+1] )
    {
        //assert_true( colSize[jj] > 0 );
        elmt = col[jj];
        //the diagonal term should be first
        //assert_true( elmt->indx == jj );
        if ( colSize[jj] == 1 )
        {
            Y[jj] += elmt->val * X[jj];
        }
        else {
            last = elmt + colSize[jj];
            real X1 = X[jj];
            real Y1 = Y[jj];
            real dY1 = elmt->val * X[jj];
            ++elmt;
            
#if ( 0 ) //0 = optimization: manual loop unrolling
            //this loop can be unrolled
            for ( ; elmt < last; ++elmt )
            {
                Element   E = *elmt;
                //assert_true( E.indx != jj );
                Y[E.indx]  +=  E.val * X1;
                dY1        +=  E.val * X[E.indx];
            }
#else
            //this is uses loop-unrolling for (maybe) faster code
            const int roll = 2;
            
            //printf("jj=%i, %i\n", jj, colSize[jj]-1);
            int start = (colSize[jj]-1) % roll;
            for ( int ii=0; ii < start; ++ii, ++elmt )
            {
                Element   E = *elmt;
                //assert_true( E.indx != jj );
                Y[E.indx]  +=  E.val * X1;
                dY1        +=  E.val * X[E.indx];
            }
            
            for ( ; elmt < last; elmt += roll )
            {
                //this is an unrolled loop:
                Element  E0 = elmt[0];
                Element  E1 = elmt[1];
                //assert_true( E0.indx != jj );
                //assert_true( E1.indx != jj );
                Y[E0.indx] += E0.val * X1;
                Y[E1.indx] += E1.val * X1;
                dY1        += E0.val * X[E0.indx] + E1.val * X[E1.indx];
            }
            //assert_true( elmt == last );
#endif
            Y[jj] = Y1 + dY1;
        }
    }
}

#ifdef MATRIX_USES_INTEL_SSE3

#include <pmmintrin.h>
#warning "Using SSE3 implementation"

void MatrixSparseSymmetric4::vecMulAddIso2D( const real* X, real* Y ) const
{
    //we only visit the non-empty columns
    for ( index_type jj = nCol[0]; jj < mxSize; jj = nCol[jj+1] )
    {
        int ll = 2 * jj;
        int siz = colSize[jj];
        
        Element * e = col[jj];
        __m128d x = _mm_load_pd(X+ll);
        
        //the diagonal term should be first
        __m128d t, a = _mm_loaddup_pd(&(e->val));
        ++e;
        
        __m128d y = _mm_load_pd(Y+ll);
        __m128d dy = _mm_mul_pd(a, x);
        
        for ( int kk = 1;  kk < siz; ++kk, ++e )
        {
            a = _mm_loaddup_pd(&(e->val));
            int ii  = 2 * e->indx;
            
            t = _mm_add_pd(_mm_load_pd(Y+ii), _mm_mul_pd(a, x));
            dy = _mm_add_pd(dy, _mm_mul_pd(a, _mm_load_pd(X+ii)));
            _mm_store_pd(Y+ii, t);
        }
        _mm_store_pd(Y+ll, _mm_add_pd(y, dy));
    }
}

#else

void MatrixSparseSymmetric4::vecMulAddIso2D( const real* X, real* Y ) const
{
    int kk, ii, siz;
    real val, X1, X2, Y1, Y2, P1, P2;
    Element * e;
    
    //we only visit the non-empty columns
    for ( index_type jj = nCol[0]; jj < mxSize; jj = nCol[jj+1] )
    {
        int ll = 2 * jj;
        siz = colSize[jj];
        assert_true( siz > 0 );
        
        e = col[jj];
        
        X1 = X[ll  ];
        X2 = X[ll+1];
        
        //the diagonal term should be first
        assert_true( e->indx == jj );
        val = e->val;
        ++e;
        kk = 1;
        Y1 = Y[ll  ] + val * X1;
        Y2 = Y[ll+1] + val * X2;
        
        while ( kk < siz )
        {
            val = e->val;
            ii  = 2 * e->indx;
            
            assert_true( ii != ll );
            
            P1 = val * X1;
            P2 = val * X2;
            
            Y1 += val * X[ii  ];
            Y2 += val * X[ii+1];
            
            ++kk;
            ++e;
            
            Y[ii  ] += P1;
            Y[ii+1] += P2;
        }
        
        Y[ll  ] = Y1;
        Y[ll+1] = Y2;
    }
}

#endif

//------------------------------------------------------------------------------
void MatrixSparseSymmetric4::vecMulAddIso3D( const real* X, real* Y ) const
{
    int kk, ii, siz;
    real val, X1, X2, X3, Y1, Y2, Y3;
    real P1, P2, P3;
    Element * e;
    
    for ( index_type jj = nCol[0]; jj < mxSize; jj = nCol[jj+1] )
    {
        int ll = 3 * jj;
        siz = colSize[jj];
        assert_true( siz > 0 );
        
        e = col[jj];
        
        X1 = X[ll  ];
        X2 = X[ll+1];
        X3 = X[ll+2];
        
        //the diagonal term should be first
        assert_true( e->indx == jj );
        val = e->val;
        ++e;
        kk = 1;
        Y1 = Y[ll  ] + val * X1;
        Y2 = Y[ll+1] + val * X2;
        Y3 = Y[ll+2] + val * X3;
        
        while ( kk < siz )
        {
            val = e->val;
            ii  = 3 * e->indx;
            
            assert_true( ii != ll );
            
            P1 = val * X1;
            P2 = val * X2;
            P3 = val * X3;
            
            Y1 += val * X[ii  ];
            Y2 += val * X[ii+1];
            Y3 += val * X[ii+2];
            
            ++kk;
            ++e;
            
            Y[ii  ] += P1;
            Y[ii+1] += P2;
            Y[ii+2] += P3;
        }
        
        Y[ll  ] = Y1;
        Y[ll+1] = Y2;
        Y[ll+2] = Y3;
    }
}

#endif


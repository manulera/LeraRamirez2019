// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matsparsebandsym.h"
#include "cblas.h"
#include "smath.h"
#include <iomanip>
#include <sstream>

#define FASTER_CODE


const int AVAILABLE_CELL = -1;
const int LAST_IN_COLUMN = -2;


//------------------------------------------------------------------------------
MatrixSparseBandSymmetric::MatrixSparseBandSymmetric()
{
    mxSize    = 0;
    mxAllocated = 0;
    diag      = 1;
    mxAllocated_diag = 0;
    val       = 0;
    mxCol     = 0;
    mxRow     = 0;
}


//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::allocate( const unsigned int sz, const unsigned int di )
{
    mxSize = sz;
    
    if ( mxSize > mxAllocated )
    {
        if ( val ) delete[] val;
        val = 0;
        
        real ** mxCol_new = new real*[mxSize];
        int  ** mxRow_new = new  int*[mxSize];
        
        unsigned int ii = 0;
        if ( mxCol )
        {
            for ( ; ii < mxAllocated; ++ii )
            {
                mxCol_new[ii] =  mxCol[ii];
                mxRow_new[ii] =  mxRow[ii];
            }
            delete[] mxCol;
            delete[] mxRow;
        }
        for ( ; ii < mxSize; ++ii )
        {
            mxCol_new[ii] = 0;
            mxRow_new[ii] = 0;
        }
        mxCol = mxCol_new;
        mxRow = mxRow_new;
        mxAllocated = mxSize;
    }
    
    assert_true( di < 20 );
    if ( ( val == 0 ) || ( ( di + 1 ) * mxSize > mxAllocated_diag ) )
    {
        if ( val )  delete[] val;
        mxAllocated_diag = ( di + 1 ) * mxAllocated;
        val = new real[mxAllocated_diag];
    }
    diag = di;
}


//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::deallocate()
{
    if ( val ) delete[] val;
    if ( mxCol )
    {
        for ( unsigned int ii = 0; ii < mxAllocated; ++ii )
            if ( mxCol[ii] )
            {
                delete[] mxCol[ii];
                delete[] mxRow[ii];
            }
        delete[] mxCol;
        delete[] mxRow;
    }
}


//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::allocateColumn( const unsigned int jj, unsigned int sz )
{
    assert_true( jj < mxSize );
    assert_true( sz > 0 );
    //printf("new S-COL %i %i\n", jj, sz );
    const unsigned chunk = 4;
    sz = ( sz + chunk - 1 ) & ~( chunk -1 );

    real* mxCol_new  = new real[sz];
    int*  mxRow_new  = new  int[sz];
    
    unsigned int ii = 0;
    if ( mxCol[jj] )
    {
        for ( ; mxRow[jj][ii] != LAST_IN_COLUMN ; ++ii )
        {
            mxCol_new[ii] =  mxCol[jj][ii];
            mxRow_new[ii] =  mxRow[jj][ii];
        }
        
        delete[] mxCol[jj];
        delete[] mxRow[jj];
    }
    for ( ; ii < sz-1 ; ++ii )
        mxRow_new[ii] = AVAILABLE_CELL;
    mxRow_new[sz-1] = LAST_IN_COLUMN;
    
    mxCol[jj]  = mxCol_new;
    mxRow[jj]  = mxRow_new;
}


//------------------------------------------------------------------------------
//allocate the position if necessary:
real& MatrixSparseBandSymmetric::operator()( index_type ii, index_type jj )
{
    assert_true( ii < mxSize );
    assert_true( jj < mxSize );
    
    // swap to get ii <= jj (address upper triangle)
    if ( jj < ii )
        std::swap(ii, jj);
    
    if ( jj <= ii+diag )
        return val[ ii+diag*(jj+1) ];
    
    int const* row = mxRow[jj];
    if ( row )
    {
        unsigned int kk = 0;
        for ( kk = 0; row[kk] >= 0; ++kk )
            if ( row[kk] == (int)ii )
                return mxCol[jj][kk];
        
        if ( mxRow[jj][kk] == LAST_IN_COLUMN )
            allocateColumn(jj, kk + 1);

        assert_true( mxCol[jj] != 0 );
        assert_true( mxRow[jj] != 0 );
        assert_true( mxRow[jj][kk] == AVAILABLE_CELL );
        mxRow[jj][kk] = ii;
        mxCol[jj][kk] = 0;
        return mxCol[jj][kk];
    }
    
    allocateColumn(jj, 1);
    assert_true( mxCol[jj] != 0 );
    assert_true( mxRow[jj] != 0 );
    assert_true( mxRow[jj][0] == AVAILABLE_CELL );
    mxRow[jj][0] = ii;
    mxCol[jj][0] = 0;
    return mxCol[jj][0];
}


//------------------------------------------------------------------------------
real* MatrixSparseBandSymmetric::addr( index_type ii, index_type jj) const
{
    assert_true( ii < mxSize );
    assert_true( jj < mxSize );
    
    // swap to get ii <= jj (address upper triangle)
    if ( jj < ii )
        std::swap(ii, jj);
    
    if ( jj <= ii+diag )
        return &( val[ii + diag * (jj+1) ] );
    
    int const* row = mxRow[jj];
    if ( row )
        for ( int kk = 0; row[kk] >= 0; ++kk )
            if ( row[kk] == (int)ii )
                return &( mxCol[jj][kk] );
    
    return 0;
}


//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::makeZero()
{
    for ( index_type ii = 0; ii < (diag+1)*mxSize; ++ii )
        val[ii] = 0;
    
    for ( index_type ii = 0; ii < mxSize; ++ii )
    {
        if ( mxRow[ii] )
            for ( index_type jj = 0; mxRow[ii][jj] >= 0; ++jj )
                mxRow[ii][jj] = AVAILABLE_CELL;
    }
}

//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::scale( const real a )
{
    for (index_type ii = 0; ii < (diag+1)*mxSize; ++ii )
        val[ii] *= a;
    
    for (index_type ii = 0; ii < mxSize; ++ii ) if ( mxRow[ii] )
        for ( int jj = 0; mxRow[ii][jj] >= 0; ++jj )
            mxCol[ii][jj] *= a;
}

//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::addDiagonalBlock(real* mat, const index_type ldd, const index_type x, const unsigned nb) const
{
    assert_true( x + nb <= mxSize );
    
    //set the diagonal, colum by column :
    //set the off-diagonals, line by line :
    for (index_type ii = 0; ii < nb ; ++ii )
    {
        mat[ii+ldd*ii] += val[ ii+x+diag*( ii+x+1 ) ];
        
        for (index_type jj = ii+1; (jj <= ii + diag) && (jj < nb) ; ++jj )
        {
            mat[ii+ldd*jj] += val[ ii+x+diag*( jj+x+1 ) ];
            mat[jj+ldd*ii] += val[ ii+x+diag*( jj+x+1 ) ];
        }
    }
    
    //set the sparse values :
    for (index_type jj = 0; jj < nb; ++jj )
    {
        int * row = mxRow[jj+x];
        if ( row != 0 )
        {
            real* col = mxCol[jj+x];
            for ( ; *row >= 0 ; ++row, ++col )
            {
                if ( *row > (int)x )
                {
                    index_type ii = *row - x;
                    if ( ii < nb )
                    {
                        assert_true( ii < jj );
                        mat[ii+ldd*jj] += * col;
                        mat[jj+ldd*ii] += * col;
                        //printf("Sp %4i %4i % .4f\n", ii, jj, a );
                    }
                }
            }
        }
    }
}


 void MatrixSparseBandSymmetric::addTriangularBlock(real* mat,
                                                    const index_type ldd,
                                                    const index_type si,
                                                    const unsigned nb,
                                                    const unsigned dim) const
{
    assert_true( si + nb <= mxSize );
    
    //set the diagonal, colum by column :
    //set the off-diagonals, line by line :
    for ( index_type ii = 0; ii < nb ; ++ii )
        for ( index_type jj = ii; (jj <= ii + diag) && (jj < nb); ++jj )
            mat[dim*ii+ldd*dim*jj] += val[ ii+si+diag*(jj+si+1) ];
    
    //set the sparse values :
    for ( index_type jj = 0; jj < nb; ++jj )
    {
        int * row = mxRow[jj+si];
        if ( row != 0 )
        {
            real* col = mxCol[jj+si];
            for ( ; *row >= 0 ; ++row, ++col )
            {
                if ( *row > (int)si )
                {
                    index_type ii = *row - si;
                    if ( ii < nb )
                    {
                        assert_true( ii < jj );
                        mat[dim*ii+ldd*dim*jj] += *col;
                        //printf("Sp %4i %4i % .4f\n", ii, jj, a );
                    }
                }
            }
        }
    }
}



//------------------------------------------------------------------------------
int MatrixSparseBandSymmetric::bad() const
{
    if ( mxSize <= 0 ) return 1;
    for ( index_type jj = 0; jj < mxSize; ++jj )
    {
        if ( mxRow[jj] )
        {
            for ( int ii = 0; mxRow[jj][ii] >= 0; ++ii )
                if ( mxRow[jj][ii] >= (int)mxSize ) return 2;
        }
    }
    return 0;
}


//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::printSparse(std::ostream& os) const
{
    os.precision(8);
    //print the diagonal, colum by column :
    for (index_type ii = 0; ii < mxSize ; ++ii )
    {
        for (index_type jj = ii; (jj <= ii + diag) && (jj < mxSize); ++jj )
        {
            os << ii << " " << jj << " ";
            os << val[ ii + diag * (jj+1) ] << std::endl;
        }
    }
    
    for (index_type jj = 0; jj < mxSize; ++jj ) if ( mxRow[jj] )
    {
        for (index_type ii = 0; mxRow[jj][ii] >= 0; ++ii )
        {
            assert_true( mxRow[jj][ii] < (int)mxSize );
            os << mxRow[jj][ii] << " " << jj << " ";
            os << mxCol[jj][ii] << std::endl;
        }
    }
}

//------------------------------------------------------------------------------
bool MatrixSparseBandSymmetric::nonZero() const
{
    //check the diagonal terms:
    for ( index_type ii = 0; ii < (diag+1)*mxSize; ++ii )
        if ( 0 != val[ii] )
            return true;
    
    //check the sparse terms:
    for ( index_type jj = 0; jj < mxSize; ++jj ) if ( mxRow[jj] )
        for ( index_type ii = 0; mxRow[jj][ii] >= 0; ++ii )
            if ( mxCol[jj][ii] != 0 )
                return true;
    
    //empty matrix:
    return false;
}

//------------------------------------------------------------------------------
unsigned int MatrixSparseBandSymmetric::nbElements() const
{
    //FIX: the diagonal elements are not counted here
    unsigned int cnt = 0;
    for ( index_type jj = 0; jj < mxSize; ++jj )
    {
        if ( mxRow[jj] )
            for ( int ii = 0; mxRow[jj][ii] >= 0; ++ii )
                if ( mxCol[jj][ii] != 0 )
                    ++cnt;
    }
    return cnt;
}

//------------------------------------------------------------------------------
std::string MatrixSparseBandSymmetric::what()  const
{
    std::ostringstream msg;
    msg << "SB " << diag << " (" << nbElements() << ")";
    return msg.str();
}


//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::vecMulAdd( const real* X, real* Y )  const
{
    
    blas_xsbmv('U', mxSize, diag, 1.0, val, diag+1, X, 1, 1.0, Y, 1);
    
    for ( index_type jj = 0; jj < mxSize; ++jj )
    {
        if ( mxRow[jj] )
        {
            index_type ii = 0;
            while (1)
            {
                if ( mxRow[jj][ii] < 0 ) break;
                unsigned int kk = mxRow[jj][ii];
                assert_true( jj < mxSize );
                assert_true( kk < mxSize );
                assert_true( kk != jj );
                Y[kk] += mxCol[jj][ii] * X[jj];
                Y[jj] += mxCol[jj][ii] * X[kk];
                ++ii;
            }
        }
    }
}

#ifndef FASTER_CODE

//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::vecMulAddIso2D( const real* X, real* Y ) const
{
    const int D = 2;
    blas_xsbmv('U', size, diag, 1.0, val, diag+1, X+0, D, 1.0, Y+0, D);
    blas_xsbmv('U', size, diag, 1.0, val, diag+1, X+1, D, 1.0, Y+1, D);
    
    int kk;
    for ( index_type jj = 0; jj < mxSize; ++jj )
    {
        if ( mxRow[jj] )
            for ( index_type ii = 0; ( kk = mxRow[jj][ii] ) >= 0; ++ii )
            {
                assert_true( jj < mxSize );
                assert_true( kk < mxSize );
                assert_true( kk != jj );
                blas_xaxpy( D, mxCol[jj][ii], X + D*jj, 1, Y + D*kk, 1 );
                blas_xaxpy( D, mxCol[jj][ii], X + D*kk, 1, Y + D*jj, 1 );
            }
    }
}

//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::vecMulAddIso3D( const real* X, real* Y ) const
{
    const int D = 3;
    blas_xsbmv('U', size, diag, 1.0, val, diag+1, X+0, D, 1.0, Y+0, D);
    blas_xsbmv('U', size, diag, 1.0, val, diag+1, X+1, D, 1.0, Y+1, D);
    blas_xsbmv('U', size, diag, 1.0, val, diag+1, X+2, D, 1.0, Y+2, D);
    
    int kk;
    for ( index_type jj = 0; jj < mxSize; ++jj )
    {
        if ( mxRow[jj] )
            for ( index_type ii = 0; ( kk = mxRow[jj][ii] ) >= 0; ++ii )
            {
                assert_true( jj < mxSize );
                assert_true( kk < mxSize );
                assert_true( kk != jj );
                blas_xaxpy( D, mxCol[jj][ii], X + D*jj, 1, Y + D*kk, 1 );
                blas_xaxpy( D, mxCol[jj][ii], X + D*kk, 1, Y + D*jj, 1 );
            }
    }
}

#else
//========================================================================
//                        FASTER_CODE
//========================================================================

//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::vecMulAddIso2D( const real* X, real* Y ) const
{
    const int D = 2;
    int * row;
    real* col;
    real* Y1, * Y2, X1, X2;
    
    blas_xsbmv('U', mxSize, diag, 1.0, val, diag+1, X+0, D, 1.0, Y+0, D);
    blas_xsbmv('U', mxSize, diag, 1.0, val, diag+1, X+1, D, 1.0, Y+1, D);
    
    for ( index_type jj = 0, ll = 0; jj < mxSize; ++jj, ll += D )
    {
        row = mxRow[jj];
        if ( row != 0 )
        {
            col = mxCol[jj];
            
            X1 = X[ll  ];
            X2 = X[ll+1];
            
            Y1 = &Y[ll  ];
            Y2 = &Y[ll+1];
            
            while (1)
            {
                if ( *row < 0 ) break;
                unsigned int kk = 2 * (*row);
                assert_true( kk < 2 * mxSize );
                assert_true( kk != ll );
                
                Y[kk  ] += (*col) * X1;
                Y[kk+1] += (*col) * X2;
                
                *Y1 += (*col) * X[kk  ];
                *Y2 += (*col) * X[kk+1];
                
                ++row;
                ++col;
            }
        }
    }
}

//------------------------------------------------------------------------------
void MatrixSparseBandSymmetric::vecMulAddIso3D( const real* X, real* Y ) const
{
    const int D = 3;
    int * row;
    real* col;
    real* Y1, * Y2, *Y3, X1, X2, X3;
    
    blas_xsbmv('U', mxSize, diag, 1.0, val, diag+1, X+0, D, 1.0, Y+0, D);
    blas_xsbmv('U', mxSize, diag, 1.0, val, diag+1, X+1, D, 1.0, Y+1, D);
    blas_xsbmv('U', mxSize, diag, 1.0, val, diag+1, X+2, D, 1.0, Y+2, D);
    
    for ( index_type jj = 0, ll = 0; jj < mxSize; ++jj, ll += D )
    {
        row = mxRow[jj];
        if ( row != 0 )
        {
            col = mxCol[jj];
            
            X1 = X[ll  ];
            X2 = X[ll+1];
            X3 = X[ll+2];
            
            Y1 = &Y[ll  ];
            Y2 = &Y[ll+1];
            Y3 = &Y[ll+2];
            
            while (1)
            {
                if ( *row < 0 ) break;
                unsigned int kk = D * (*row);
                assert_true( kk < D * mxSize );
                assert_true( kk != ll );
                
                Y[kk  ]   += (*col) * X1;
                Y[kk+1]   += (*col) * X2;
                Y[kk+2]   += (*col) * X3;
                
                *Y1    += (*col) * X[kk  ];
                *Y2    += (*col) * X[kk+1];
                *Y3    += (*col) * X[kk+2];
                
                ++row;
                ++col;
            }
        }
    }
}

#endif


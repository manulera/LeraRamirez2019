// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matsym.h"
#include "cblas.h"


//------------------------------------------------------------------------------
MatrixSymmetric::MatrixSymmetric()
{
    msAllocated = 0;
    val         = 0;
    in_charge = 1;
}

//------------------------------------------------------------------------------
void MatrixSymmetric::allocate(unsigned alc)
{
    if ( alc > msAllocated )
    {
        msAllocated = alc;
        if ( val )
            delete[] val;
        val = new real[alc*alc];
    }
}

//------------------------------------------------------------------------------
void MatrixSymmetric::deallocate()
{
    if ( in_charge )
    {
        if ( val ) delete[] val;
    }
    msAllocated = 0;
    val = 0;
}

//------------------------------------------------------------------------------
void MatrixSymmetric::makeZero()
{
    for ( index_type ii = 0; ii < size_ * size_; ++ii )
        val[ii] = 0;
}

//------------------------------------------------------------------------------
void MatrixSymmetric::scale( real alpha )
{
    for ( index_type ii = 0; ii < size_ * size_; ++ii )
        val[ii] *= alpha;
}

//------------------------------------------------------------------------------
real& MatrixSymmetric::operator()( index_type x, index_type y)
{
    assert_true( x < size_ );
    assert_true( y < size_ );
    if ( x < y )
        return val[ x + msLDD * y ];
    else
        return val[ y + msLDD * x ];
}

//------------------------------------------------------------------------------
real* MatrixSymmetric::addr( index_type x, index_type y) const
{
    assert_true( x < size_ );
    assert_true( y < size_ );
    if ( x < y )
        return val + x + msLDD * y;
    else
        return val + y + msLDD * x;
}

//------------------------------------------------------------------------------
bool MatrixSymmetric::nonZero() const
{
    return true;
}

//------------------------------------------------------------------------------
unsigned int MatrixSymmetric::nbElements() const
{
    return size_ * size_;
}

//------------------------------------------------------------------------------
std::string MatrixSymmetric::what() const
{
    return "full-symmetric";
}

//------------------------------------------------------------------------------
void MatrixSymmetric::vecMulAdd( const real* X, real* Y ) const
{
    blas_xsymv( 'U', size_, 1.0, val, size_, X, 1, 1.0, Y, 1 );
}

//------------------------------------------------------------------------------
void MatrixSymmetric::vecMulAddIso2D( const real* X, real* Y ) const
{
    blas_xsymv( 'U', size_, 1.0, val, size_, X+0, 2, 1.0, Y+0, 2 );
    blas_xsymv( 'U', size_, 1.0, val, size_, X+1, 2, 1.0, Y+1, 2 );
}

//------------------------------------------------------------------------------
void MatrixSymmetric::vecMulAddIso3D( const real* X, real* Y ) const
{
    blas_xsymv( 'U', size_, 1.0, val, size_, X+0, 3, 1.0, Y+0, 3 );
    blas_xsymv( 'U', size_, 1.0, val, size_, X+1, 3, 1.0, Y+1, 3 );
    blas_xsymv( 'U', size_, 1.0, val, size_, X+2, 3, 1.0, Y+2, 3 );
}



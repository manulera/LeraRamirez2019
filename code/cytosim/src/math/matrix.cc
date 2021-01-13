// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "assert_macro.h"
#include "matrix.h"
#include "cblas.h"
#include <iomanip>

//------------------------------------------------------------------------------
real Matrix::value(const index_type x, const index_type y) const
{
    real* v = addr( x, y );
    if ( v == 0 )
        return 0;
    else
        return *v;
}

//------------------------------------------------------------------------------
real Matrix::maxNorm() const
{
    const unsigned int sz = size();
    real result = 0;
    for ( unsigned int ii = 0; ii < sz; ++ii )
    {
        for ( unsigned int jj = 0; jj < sz; ++jj )
        {
            real* v = addr( ii, jj );
            if ( v  &&  ( *v > result ) )
                result = *v;
        }
    }
    return result;
}

//------------------------------------------------------------------------------
bool Matrix::nonZero() const
{
    const unsigned int sz = size();
    for ( unsigned int ii = 0; ii < sz; ++ii )
        for ( unsigned int jj = 0; jj < sz; ++jj )
            if ( 0 != value( ii, jj ) )
                return true;
    return false;
}

//------------------------------------------------------------------------------
unsigned int  Matrix::nbElements() const
{
    const unsigned int sz = size();
    unsigned int result = 0;
    for ( unsigned int ii = 0; ii < sz; ++ii )
        for ( unsigned int jj = 0; jj < sz; ++jj )
            result += ( 0 != value( ii, jj ) );
    return result;
}

//------------------------------------------------------------------------------
void Matrix::copyBlock(real* mat, unsigned ldd, index_type sx, unsigned nx, index_type sy, unsigned ny) const
{
    assert_true( sx + nx < size() );
    assert_true( sy + ny < size() );
    
    for ( unsigned ii = 0; ii < nx; ++ii )
        for ( unsigned jj = 0; jj < ny; ++jj )
            mat[ii + ldd * jj] = value( sx + ii, sy + jj );
}


void Matrix::addDiagonalBlock(real* mat, const index_type ldd, const index_type si, const unsigned nb) const
{
    assert_true( si + nb < size() );

    for ( unsigned ii = 0; ii < nb; ++ii )
        for ( unsigned jj = 0; jj < nb; ++jj )
            mat[ ii + ldd * jj ] += value( si + ii, si + jj );
}


void Matrix::addTriangularBlock(real* mat, const index_type ldd, const index_type si, const unsigned nb, const unsigned dim) const
{
    assert_true( si + nb < size() );

    for ( unsigned ii = 0; ii < nb; ++ii )
        for ( unsigned jj = ii; jj < nb; ++jj )
            mat[ dim*ii + ldd * dim*jj ] += value( si + ii, si + jj );
}

//------------------------------------------------------------------------------
void Matrix::vecMul( const real* X, real* Y ) const
{
    blas_xzero(size(), Y);
    vecMulAdd( X, Y );
}

//------------------------------------------------------------------------------
void Matrix::printFull(std::ostream & os) const
{
    char str[32];
    const unsigned sz = size();
    //printf("%i %i\n", size, size);
    for ( unsigned ii = 0; ii < sz; ++ii )
    {
        for ( unsigned jj = 0; jj < sz; ++jj )
        {
            real * a = addr(ii,jj);
            if ( a )
            {
                snprintf(str, sizeof(str), " %9.3f", *a);
                os << str;
            }
            else
                os << "       .  ";
        }
        os << std::endl;
    }
}

//------------------------------------------------------------------------------
void Matrix::printSparse(std::ostream & os) const
{
    const unsigned int sz = size();
    for ( unsigned int ii = 0; ii < sz; ++ii )
        for ( unsigned int jj = 0; jj < sz; ++jj )
            if ( addr( ii, jj ) )
                os << ii << " " << jj << " " << std::scientific << *addr(ii, jj) << std::endl;
}




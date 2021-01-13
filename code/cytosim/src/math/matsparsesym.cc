// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matsparsesym.h"
#include "cblas.h"
#include "smath.h"

#include <iomanip>
#include <sstream>


MatrixSparseSymmetric::MatrixSparseSymmetric()
{
    allocated_ = 0;
    col_       = 0;
    col_size_  = 0;
    col_max_   = 0;
}



void MatrixSparseSymmetric::allocate(unsigned alc)
{
    if ( alc > allocated_ )
    {
        const unsigned chunk = 8;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );

        //fprintf(stderr, "MSS allocate matrix %u\n", alc);
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
    }
}


void MatrixSparseSymmetric::deallocate()
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
    }
    allocated_ = 0;
}


/// copy `cnt` elements from `src` to `dst`
void copy(unsigned cnt, MatrixSparseSymmetric::Element * src, MatrixSparseSymmetric::Element * dst)
{
    for ( unsigned ii = 0; ii < cnt; ++ii )
        dst[ii] = src[ii];
}


void MatrixSparseSymmetric::allocateColumn(const index_type jj, unsigned alc)
{
    assert_true( jj < size_ );
    
    if ( alc > col_max_[jj] )
    {
        //fprintf(stderr, "MSS allocate column %i size %u\n", jj, alc);
        const unsigned chunk = 4;
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
real& MatrixSparseSymmetric::operator()(index_type ii, index_type jj)
{
    assert_true( ii < size_ );
    assert_true( jj < size_ );
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
                return e->val;
            ++e;
        }
    }
    
    // add the requested term at the end:
    unsigned n = col_size_[jj];

    // allocate space for new Element if necessary:
    if ( n >= col_max_[jj] )
    {
        allocateColumn(jj, n+1);
        col = col_[jj];
    }
    
    col[n].reset(ii);
    ++col_size_[jj];
    
    //printColumn(jj);
    return col[n].val;
}


real* MatrixSparseSymmetric::addr(index_type ii, index_type jj) const
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

void MatrixSparseSymmetric::makeZero()
{
    for ( index_type jj = 0; jj < size_; ++jj )
        col_size_[jj] = 0;
}


bool MatrixSparseSymmetric::nonZero() const
{
    //check for any non-zero sparse term:
    for ( unsigned jj = 0; jj < size_; ++jj )
        for ( unsigned kk = 0 ; kk < col_size_[jj] ; ++kk )
            if ( col_[jj][kk].val )
                return true;
    
    //if here, the matrix is empty
    return false;
}


void MatrixSparseSymmetric::scale(const real alpha)
{
    for ( index_type jj = 0; jj < size_; ++jj )
        for ( index_type n = 0; n < col_size_[jj]; ++n )
            col_[jj][n].val *= alpha;
}



void MatrixSparseSymmetric::addTriangularBlock(real* mat, const unsigned ldd,
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
            if ( si <= ii && ii < up )
            {
                if ( ii < jj )
                    mat[dim*(ii-si+ldd*(jj-si))] += col_[jj][n].val;
                else
                    mat[dim*(jj-si+ldd*(ii-si))] += col_[jj][n].val;
            }
        }
    }
}


void MatrixSparseSymmetric::addDiagonalBlock(real* mat, unsigned ldd,
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
            if ( si <= ii && ii < up )
            {
                //printf("MSS1 %4i %4i % .4f\n", ii, jj, a);
                mat[ii-si+ldd*(jj-si)] += col_[jj][n].val;
                if ( jj != ii )
                    mat[jj-si+ldd*(ii-si)] += col_[jj][n].val;
            }
        }
    }
}


int MatrixSparseSymmetric::bad() const
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


unsigned MatrixSparseSymmetric::nbElements() const
{
    //all allocated elements are counted, even if zero
    unsigned cnt = 0;
    for ( unsigned jj = 0; jj < size_; ++jj )
        cnt += col_size_[jj];
    return cnt;
}

//------------------------------------------------------------------------------
#pragma mark -

std::string MatrixSparseSymmetric::what() const
{
    std::ostringstream msg;
    msg << "mSS (" << nbElements() << ")";
    return msg.str();
}


void MatrixSparseSymmetric::printSparse(std::ostream & os) const
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


void MatrixSparseSymmetric::printColumns(std::ostream & os)
{
    os << "MSS size " << size_ << ":";
    for ( int jj = 0; jj < size_; ++jj )
    {
        os << "\n   " << jj << "   " << col_size_[jj];
    }
    os << std::endl;
}


void MatrixSparseSymmetric::printColumn(std::ostream & os, const index_type jj)
{
    Element const* col = col_[jj];
    os << "MSS col " << jj << ":";
    for ( unsigned n = 0; n < col_size_[jj]; ++n )
    {
        os << "\n" << col[n].inx << " :";
        os << " " << col[n].val;
    }
    os << std::endl;
}


//------------------------------------------------------------------------------
#pragma mark -

void MatrixSparseSymmetric::prepareForMultiply()
{
}


void MatrixSparseSymmetric::vecMulAdd(const real* X, real* Y) const
{
    for ( index_type jj = 0; jj < size_; ++jj )
    {
        for ( unsigned kk = 0 ; kk < col_size_[jj] ; ++kk )
        {
            const index_type ii = col_[jj][kk].inx;
            const real a = col_[jj][kk].val;
            Y[ii] += a * X[jj];
            if ( ii != jj )
                Y[jj] += a * X[ii];
        }
    }
}


void MatrixSparseSymmetric::vecMulAddIso2D(const real* X, real* Y) const
{
    for ( index_type jj = 0; jj < size_; ++jj )
    {
        const index_type Djj = 2 * jj;
        for ( unsigned kk = 0 ; kk < col_size_[jj] ; ++kk )
        {
            const index_type Dii = 2 * col_[jj][kk].inx;
            const real  a = col_[jj][kk].val;
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


void MatrixSparseSymmetric::vecMulAddIso3D(const real* X, real* Y) const
{
    for ( index_type jj = 0; jj < size_; ++jj )
    {
        const index_type Djj = 3 * jj;
        for ( unsigned kk = 0 ; kk < col_size_[jj] ; ++kk )
        {
            const index_type Dii = 3 * col_[jj][kk].inx;
            const real  a =     col_[jj][kk].val;
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


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matblock.h"
#include "cblas.h"

#include <iomanip>
#include <sstream>

//------------------------------------------------------------------------------
MatrixBlock::MatrixBlock()
{
    mxAllocated  = 0;
    block_     = 0;
    block_alc  = 0;
    block_nb   = 0;
    block_size = 0;
}

//------------------------------------------------------------------------------
void MatrixBlock::deallocate()
{
    if ( block_ == 0 ) return;
    
    for ( unsigned int ii=0; ii < block_nb; ++ii )
    {
        if ( block_[ii] )
            delete[] block_[ii];
    }
    
    delete[] block_;
    delete[] block_size;
    delete[] block_alc;
    mxAllocated = 0;
}

//------------------------------------------------------------------------------
void MatrixBlock::allocate( const unsigned int nbblock )
{
    assert_true( nbblock > 0 );
    block_nb = nbblock;
    if ( nbblock > mxAllocated )
    {
        const unsigned chunk = 4;
        unsigned nba = ( nbblock + chunk - 1 ) & ~( chunk -1 );

        //printf("new block-matrix sz %i\n", nbblock );
        real** block_new = new real*[nba];
        
        unsigned int * block_alc_new  = new unsigned int[nba];
        unsigned int * block_size_new = new unsigned int[nba];
        
        unsigned int ii = 0;
        
        if ( block_ )
        {
            for ( ii = 0; ii < mxAllocated; ++ii )
            {
                block_new[ii]            = block_[ii];
                block_alc_new[ii]        = block_alc[ii];
                block_size_new[ii]       = block_size[ii];
            }
            delete[] block_;
            delete[] block_alc;
            delete[] block_size;
         }
        
        for ( ; ii < nba; ++ii )
        {
            block_new[ii]            = 0;
            block_alc_new[ii]        = 0;
            block_size_new[ii]       = 0;
        }
        
        block_         = block_new;
        block_alc      = block_alc_new;
        block_size     = block_size_new;
        mxAllocated    = nba;
    }
}

//------------------------------------------------------------------------------
void MatrixBlock::setBlockSize( const unsigned int bb, unsigned int sz )
{
    assert_true( block_ != 0 );
    assert_true( bb < block_nb );
    
    block_size[bb] = sz;
    
    if ( block_alc[bb] < sz )
    {
        const unsigned chunk = 4;
        sz = ( sz + chunk - 1 ) & ~( chunk -1 );

        //printf("MatrixBlock::new block %i size %i\n", bb, sz );
        block_alc[bb] = sz;
        if ( block_[bb] )
            delete[] block_[bb];
        block_[bb] = new real[sz*sz];
    }
}

//------------------------------------------------------------------------------
unsigned int MatrixBlock::calculateSize()
{
    size_ = 0;
    for ( unsigned int ii = 0; ii < block_nb; ++ii )
        size_ += block_size[ii];
    return size_;
}



//------------------------------------------------------------------------------
real* MatrixBlock::addr(index_type ii, index_type jj) const
{
    unsigned int bx = 0;  //block index on x
    unsigned int by = 0;  //block index on y
    
    while ( ii >= block_size[bx] ) ii -= block_size[bx++];
    while ( jj >= block_size[by] ) jj -= block_size[by++];
    
    if ( bx != by ) // the element is not on a diagonal block
        return 0;
        
    if ( block_[bx] == 0 )  //the block in not allocated
        return 0;
    
    assert_true( ii < block_size[bx] );
    assert_true( jj < block_size[bx] );
    assert_true( block_size[bx] <= block_alc[bx] );
    
    return & block_[bx][ii + block_size[bx] * jj];
}

//------------------------------------------------------------------------------
real& MatrixBlock::operator()( index_type x, index_type y )
{
    return * addr( x, y );
}

//------------------------------------------------------------------------------
unsigned int MatrixBlock::nbElements() const
{
    unsigned int result=0;
    for ( unsigned int bb = 0; bb < block_nb; ++bb )
        result += block_size[bb] * block_size[bb];
    return result;
}

//------------------------------------------------------------------------------
unsigned int MatrixBlock::maxBlockSize()  const
{
    unsigned int result = 0;
    for ( unsigned int bb=0; bb < block_nb; ++bb )
        if ( result < block_size[bb] ) result = block_size[bb];
    return result;
}

//------------------------------------------------------------------------------
void MatrixBlock::setBlockToZero( const unsigned int bb )
{
    real* BS = block_[bb];
    if ( BS )
    {
        for ( unsigned int kk = 0; kk < block_size[bb] * block_size[bb]; ++kk )
            BS[kk] = 0;
    }
}

//------------------------------------------------------------------------------
void MatrixBlock::makeZero()
{
    for ( unsigned int ii = 0; ii < block_nb ; ++ii )
        setBlockToZero( ii );
}

//------------------------------------------------------------------------------
void MatrixBlock::scaleBlock( const unsigned int bb, const real a )
{
    real* BS = block_[bb];
    if ( BS ) {
        for ( unsigned int kk = 0; kk < block_size[bb] * block_size[bb]; ++kk )
            BS[kk] = a * BS[kk];
    }
}

//------------------------------------------------------------------------------
void MatrixBlock::scale(const real a)
{
    for ( unsigned int bb = 0; bb < block_nb ; ++bb )
        scaleBlock( bb, a );
}

//------------------------------------------------------------------------------
void MatrixBlock::vecMulAdd( const real* X, real* Y )  const
{
    index_type xx = 0;
    for ( unsigned int bb = 0; bb < block_nb; ++bb )
    {
        assert_true( block_[bb] );
        blas_xgemv('N', block_size[bb], block_size[bb], 1.0, 
                   block_[bb], block_size[bb], X+xx, 1, 1.0, Y+xx, 1);
        xx += block_size[bb];
    }
}

//------------------------------------------------------------------------------
void MatrixBlock::vecMul( const real* X, real* Y )  const
{
    index_type xx = 0;
    for ( unsigned int bb = 0; bb < block_nb; ++bb )
    {
        assert_true( block_[bb] );
        blas_xgemv('N', block_size[bb], block_size[bb], 1.0, 
                   block_[bb], block_size[bb], X+xx, 1, 0.0, Y+xx, 1);
        xx += block_size[bb];
    }
}

//------------------------------------------------------------------------------
real MatrixBlock::maxNorm() const
{
    real mx = 0;
    for ( unsigned int bb=0; bb < block_nb; ++bb )
    {
        real* X = block_[bb];
        index_type m = block_size[bb];
        for ( index_type ii = 0; ii < m*m; ++ii )
            if ( X[ii] > mx ) mx = X[ii];
    }
    return mx;
}


//------------------------------------------------------------------------------
/// printf debug function in sparse mode: i, j : value
void MatrixBlock::printSparse(std::ostream & os) const
{
    int offset=0;
    os.precision(8);
    for ( unsigned int bb=0; bb < block_nb; ++bb )
    {
        const unsigned int bsize = block_size[bb];
        for ( unsigned int ii=0; ii < bsize; ++ii )
            for ( unsigned int jj=0; jj < bsize; ++jj )
            {
                os << ii+offset << " " << jj+offset << " ";
                os << std::setw(16) << block_[bb][ii + bsize * jj] << std::endl;
            }
        offset += bsize;
    }
}

//------------------------------------------------------------------------------
std::string MatrixBlock::what() const
{
    std::ostringstream msg;
    msg << "mB " << block_nb;
    return msg.str();
}


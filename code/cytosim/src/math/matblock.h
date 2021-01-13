// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATBLOCK_H
#define MATBLOCK_H

#include "assert_macro.h"
#include "matrix.h"
#include <cstdio>

/// A non-symmetric real Matrix, diagonal by blocks.
class MatrixBlock : public Matrix
{    
private:

    /// size of memory which has been allocated
    unsigned int mxAllocated;
    
    /// number of blocks
    unsigned int block_nb;
    
    /// array of pointers to the blocks
    real** block_;
    
    /// array containing the size of each block
    unsigned int*   block_size;
    
    /// array specifying the allocated size of each block
    unsigned int*   block_alc;
        
public:

    /// default constructor
    MatrixBlock();
    
    /// the deallocation
    void deallocate();
    
    /// allocate the matrix to be able to hold nb_block (arg 1) blocks
    void allocate( unsigned int nb_block );
    
    /// allocate block b (arg 1) to be capable of holding (size*size) (arg 2)
    void setBlockSize( unsigned int b, unsigned int size);
    
    /// default destructor
    virtual ~MatrixBlock() { deallocate(); }
    
    
    /// return the address of first element in block ii
    real* block( const unsigned int ii ) const
    {
        assert_true( block_ );
        assert_true( ii < block_nb );
        assert_true( block_[ii] );
        assert_true( block_size[ii] <= block_alc[ii] );
        return block_[ii];
    }
    
    /// returns the number of blocks
    unsigned int nbBlocks() const 
    {
        return block_nb;
    }
    
    /// calculate and return the total size of the matrix
    unsigned int calculateSize();
        
    /// returns the size of block ii
    unsigned int blockSize( const unsigned int ii ) const 
    {
        assert_true( block_ );
        assert_true( ii < block_nb );
        return block_size[ii];
    }
    
    /// the address holding element (ii, jj)
    real* addr( index_type ii, index_type jj ) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    virtual real& operator()( index_type x, index_type y );
    
    /// reset all the values in block ii
    void setBlockToZero( unsigned int ii );
    
    /// reset the entire matrix
    void makeZero();
    
    /// scale block ii
    void scaleBlock( unsigned int ii, real a );
    
    /// scale the entire matrix
    void scale( real a );
    
    /// total number of elements allocated
    unsigned int  nbElements() const;
    
    /// size of the biggest block
    unsigned int  maxBlockSize() const;
    
    /// vector multiplication: Y <- M * X
    void vecMulAdd( const real* X, real* Y) const;
    
    /// vector multiplication: Y <- M * X
    void vecMul( const real* X, real* Y) const;
    
    /// isotropic vector multiplication: Y = Y + M * X, size(X) = size(Y) = 2 * size(M)
    virtual void vecMulAddIso2D( const real* X, real* Y ) const 
    {
        ABORT_NOW("unfinished");
    }
    
    /// isotropic vector multiplication: Y = Y + M * X, size(X) = size(Y) = 3 * size(M)
    virtual void vecMulAddIso3D( const real* X, real* Y ) const 
    {
        ABORT_NOW("unfinished");
    }
    
    /// maximum of the absolute value of all elements
    real maxNorm() const;
    
    /// printf debug function in sparse mode: i, j : value
    void printSparse(std::ostream &) const;
    
    /// returns a string which a description of the type of matrix
    std::string what() const;
};



#endif

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


#ifndef MATRIXSPARSEBANDSYM_H
#define MATRIXSPARSEBANDSYM_H

#include <cstdio>
#include "matrix.h"

/// a sparse real symmetric Matrix, with a variable number of full diagonals
class MatrixSparseBandSymmetric : public Matrix
{

private:
    
    /// size of the matrix
    unsigned int mxSize;
    
    /// size of memory which has been allocated
    unsigned int mxAllocated;
    
    /// size of memory allocate to hold the diagonals
    unsigned int mxAllocated_diag;
    
    // number of off-diagonals stored in val
    unsigned int diag;
    
    // the diagonal's values
    real * val;
    
    // the off-diagonal's values: normal sparse
    
    // array [size][?] holding the values for each column
    real ** mxCol;
    
    // array [size][?] holding the line index for each column
    int **  mxRow;
    
    void allocateColumn( unsigned int, unsigned int );
    
public:
    
    //size of (square) matrix
    unsigned int size() const { return mxSize; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSparseBandSymmetric();
    
    /// default destructor
    virtual ~MatrixSparseBandSymmetric()  { deallocate(); }
    
    /// set all the element to zero
    void makeZero();
    
    /// allocate the matrix, specifying the number of off diagonal
    void allocate( unsigned int sz, unsigned int diagonals );
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate( unsigned int sz ) 
    {
        allocate( sz, 1 );
    }
    
    /// returns the address of element at (x, y), no allocation is done
    real* addr( index_type x, index_type y ) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()( index_type x, index_type y );
    
    /// scale the matrix by a scalar factor
    void scale( real a );
    
    /// add  diagonal block ( x, x, x+sx, x+sx ) to `mat`
    void addDiagonalBlock(real* mat, index_type ldd, index_type si, unsigned nb) const;

    /// add upper triangular half of ( idx, idx, idx+siz, idx+siz ) to `mat`
    void addTriangularBlock(real* mat, unsigned ldd, index_type si, unsigned nb, unsigned dim) const;

    /// multiplication of a vector: Y = Y + M * X, dim(X) = dim(M)
    void vecMulAdd( const real* X, real* Y ) const;
    
    /// 2D isotropic multiplication of a vector: Y = Y + M * X
    void vecMulAddIso2D( const real* X, real* Y ) const;
    
    /// 3D isotropic multiplication of a vector: Y = Y + M * X
    void vecMulAddIso3D( const real* X, real* Y ) const;
    
    /// true if matrix is non-zero
    bool nonZero() const;
    
    /// number of element which are non-zero
    unsigned int  nbElements() const;
    
    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// printf debug function in sparse mode: i, j : value
    void printSparse(std::ostream&) const;
    
    /// debug function
    int bad() const;
    
};

#endif

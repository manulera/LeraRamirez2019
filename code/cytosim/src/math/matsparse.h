// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATSPARSE_H
#define MATSPARSE_H

#include <cstdio>
#include "matrix.h"

/// a real (non-symmetric) sparse Matrix
/**
 This class is not used currently
 */
class MatrixSparse : public Matrix
{
    
private:
    
    /// size of the matrix
    unsigned int mxSize;
    
    /// size of memory which has been allocated
    unsigned int mxAllocated;
    
    
    // array [ size ][ ? ] holding the values for each column
    real ** mxCol;
    
    // array [ size ][ ? ] holding the line index for each column
    int  ** mxRow;
    
    // allocate column to hold nb values
    void allocateColumn( index_type column_index, unsigned int nb_values );
    
public:
    
    //size of (square) matrix
    unsigned int size() const { return mxSize; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSparse();
    
    /// default destructor
    virtual ~MatrixSparse()  { deallocate(); }
    
    /// set all the element to zero
    void makeZero();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate( unsigned int sz );
        
    /// returns the address of element at (x, y), no allocation is done
    real* addr( index_type x, index_type y ) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()( index_type x, index_type y );
    
    /// scale the matrix by a scalar factor
    void scale( real a );
    
    /// add the diagonal block ( x, x, x+sx, x+sx ) from this matrix to M
    void addDiagonalBlock(real* mat, unsigned ldd, index_type si, unsigned nb) const;
    
    /// add this' data block ( idx, idx, idx+siz, idx+siz ) to upper triangular half of `mat`
    void addTriangularBlock(real* mat, index_type ldd, index_type si, unsigned nb, unsigned dim) const;
    
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
    void printSparse(std::ostream &) const;
    
    /// debug function
    int bad() const;
};


#endif

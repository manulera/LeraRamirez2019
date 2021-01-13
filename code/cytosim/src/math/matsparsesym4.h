// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATSPARSESYM4_H
#define MATSPARSESYM4_H

#include <cstdio>
#include "matrix.h"


///real symmetric sparse Matrix, faster equivalent to MatrixSparseSymmetric
/**
 MatrixSparseSymmetric4 uses a sparse storage, with arrays of elements for each column.
 There is also an array nCol[], to specify empty columns.
 The same storage is used for multiplication.
 */
class MatrixSparseSymmetric4 : public Matrix
{
    
private:
    
    ///Element describes an element in a sparse matrix
    // The elements are stored per columns, that is how the column index is known
    struct Element 
    {
        real val;    ///< The value of the element
        index_type  indx;   ///< The index of the line
    };
    
    ///accessory function for qsort:
    static int compare( const void *, const void * );
    
private:
        
    /// size of the matrix
    unsigned &   mxSize;
    
    /// size of memory which has been allocated
    unsigned int mxAllocated;
    
    /// array [size][indx] holding the Element for the 'indx' column
    Element ** col;
    
    ///number of elements in each column: colSize[indx] is the number of elements in column 'indx'
    unsigned int  * colSize;
    
    ///numer of element allocated in each column
    unsigned int  * colMax;
    
    ///nCol[ii] is the index of the first non-empty column of index >= ii
    unsigned int  * nCol;
    
    /// allocate column to hold nb values
    void allocateColumn( index_type column_index, unsigned int nb_values );
    
public:
    
    //size of (square) matrix
    unsigned int size() const { return mxSize; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSparseSymmetric4();
    
    /// default destructor
    virtual ~MatrixSparseSymmetric4()  { deallocate(); }
    
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
    
    /// add upper triangular half of 'this' block ( idx, idx, idx+siz, idx+siz ) to `mat`
    void addTriangularBlock(real* mat, index_type ldd, index_type si, unsigned nb, unsigned dim) const;
    
    ///optional optimization that may accelerate multiplications by a vector
    void prepareForMultiply();
    
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

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATSPARSESYM_H
#define MATSPARSESYM_H

#include <cstdio>
#include "matrix.h"

///real symmetric sparse Matrix
/**
 MatrixSparseSymmetric1 uses a sparse storage, with arrays of elements for each column.
 */
class MatrixSparseSymmetric : public Matrix
{
public:
    
    /// An element of the sparse matrix
    struct Element
    {
        real        val;   ///< The value of the element
        index_type  inx;   ///< The index of the line
        
        void reset(index_type i)
        {
            inx = i;
            val = 0.0;
        }
    };
    
private:
    
    /// amount of memory which has been allocated
    unsigned   allocated_;
    
    /// array col_[c][] holds Elements of column 'c'
    Element ** col_;
    
    /// col_size_[c] is the number of Elements in column 'c'
    unsigned * col_size_;
    
    /// col_max_[c] is the number of Elements allocated in column 'c'
    unsigned * col_max_;
    
    /// allocate column to hold specified number of values
    void allocateColumn(index_type col, unsigned nb);
    
public:

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSparseSymmetric();
    
    /// default destructor
    virtual ~MatrixSparseSymmetric()  { deallocate(); }
    
    /// set all the element to zero
    void makeZero();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(unsigned sz);
    
    /// returns the address of element at (x, y), no allocation is done
    real* addr(index_type x, index_type y) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(index_type x, index_type y);
    
    /// scale the matrix by a scalar factor
    void scale(real);
    
    /// add the diagonal block ( x, x, x+sx, x+sx ) from this matrix to M
    void addDiagonalBlock(real* mat, unsigned ldd, index_type si, unsigned nb) const;
    
    /// add upper triangular half of 'this' block ( idx, idx, idx+siz, idx+siz ) to `mat`
    void addTriangularBlock(real* mat, index_type ldd, index_type si, unsigned nb, unsigned dim) const;
    
    ///optional optimization that may accelerate multiplications by a vector
    void prepareForMultiply();
    
    /// multiplication of a vector: Y = Y + M * X with dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y) const;
    
    /// 2D isotropic multiplication of a vector: Y = Y + M * X with dim(X) = 2 * dim(M)
    void vecMulAddIso2D(const real* X, real* Y) const;
    
    /// 3D isotropic multiplication of a vector: Y = Y + M * X with dim(X) = 3 * dim(M)
    void vecMulAddIso3D(const real* X, real* Y) const;
    
    /// true if matrix is non-zero
    bool nonZero() const;
    
    /// number of element which are not null
    unsigned nbElements() const;
    
    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// printf debug function in sparse mode: i, j : value
    void printSparse(std::ostream &) const;
    
    /// print content of one column
    void printColumn(std::ostream &, index_type);
    
    /// print content of one column
    void printColumns(std::ostream &);

    /// debug function
    int bad() const;
};


#endif


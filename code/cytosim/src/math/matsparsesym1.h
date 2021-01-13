// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATSPARSESYM1_H
#define MATSPARSESYM1_H

#include <cstdio>
#include "matrix.h"

#define MATRIX1_OPTIMIZED_MULTIPLY

///real symmetric sparse Matrix, with optimized multiplication
/**
 MatrixSparseSymmetric1 uses a sparse storage, with arrays of elements for each column.
 For multiplication, it uses a another format, from Numerical Recipes.
 The conversion is done when prepareForMultiply() is called
 
 Elements are stored in order of increasing index in each column.
*/
class MatrixSparseSymmetric1 : public Matrix
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
    
#ifdef MATRIX1_OPTIMIZED_MULTIPLY
    
    /// col_next_[ii] is the index of the first non-empty column of index >= ii
    index_type * col_next_;
    
    /// update col_next_[], a pointer to the next non-empty column
    void setNextColumn();

    ///array of index for the optmized multiplication
    unsigned   nmax_;
    index_type * ija_;
    real       * sa_;

#endif

public:

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSparseSymmetric1();
    
    /// default destructor
    virtual ~MatrixSparseSymmetric1()  { deallocate(); }
    
    /// set all the element to zero
    void makeZero();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(unsigned sz);
    
    /// returns the address of element at (x, y), no allocation is done
    real* addr(index_type x, index_type y) const;

    /// set the diagonal term at given index
    real& diagonal(index_type ix);
    
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

    /// printf debug function in sparse mode: i, j : value
    void printSparseArray(std::ostream &) const;
    
    /// debug function
    int bad() const;
};


#endif


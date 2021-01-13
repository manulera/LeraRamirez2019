// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATSPARSESYM5_H
#define MATSPARSESYM5_H

#include <cstdio>
#include "matrix.h"

//#define MATRIX5_OPTIMIZED_MULTIPLY

///real symmetric sparse Matrix, with optimized multiplication
/** 
MatrixSparseSymmetric5 is a sparse matrix with a full storage for the diagonal.
The plan is to use a public sparse matrix package, for example:
 -Pardiso  http://www.pardiso-project.org
 -PETSc    http://www.mcs.anl.gov/petsc/petsc-as
Pardiso is included in Intel MKL and should be a good choice
*/
class MatrixSparseSymmetric5 : public Matrix
{
    
private:
    
    ///Element describes an element in a sparse matrix
    // The elements are stored per columns, that is how the column index is known
    struct Element 
    {
        real        val;    ///< The value of the element
        index_type  line;   ///< The index of the line
        
        void reset(index_type i)
        {
            line = i;
            val = 0.0;
        }
    };
    
private:
    /// size of matrix
    unsigned &   mxSize;
    
    /// amount of memory which has been allocated
    unsigned int mxAllocated;
    
    ///the values on the diagonal
    real  * diag;
    
    /// array col[c][] holds Elements of column 'c'
    Element ** col;
    
    /// colSize[c] is the number of Elements in column 'c'
    unsigned int  * colSize;
    
    /// colMax[c] number of Elements allocated in column 'c'
    unsigned int  * colMax;
    
    /// allocate column to hold specified number of values
    void allocateColumn( index_type column_index, unsigned int nb_values );
    
#ifdef MATRIX5_OPTIMIZED_MULTIPLY
    ///nCol[ii] is the index of the first non-empty column of index >= ii
    index_type * nCol;
    
    ///array of index for the optmized multiplication
    unsigned int  nmax;
    index_type  * ija;
    real        * sa;
#endif
    
public:

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSparseSymmetric5();
    
    /// default destructor
    virtual ~MatrixSparseSymmetric5()  { deallocate(); }
    
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


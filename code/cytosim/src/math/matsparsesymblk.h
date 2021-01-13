// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATSPARSESYM3_H
#define MATSPARSESYM3_H

#include <cstdio>
#include "matrix.h"
//#include "vecprint.h"


//The block size should match the dimensionaltiy of the simulation:
#ifndef BSZ
   #include "dim.h"
   #define BSZ DIM
#endif

///real symmetric sparse Matrix
/**
 MatrixSparseSymmetric3 uses a sparse storage, with arrays of elements for each column.
 Each element is a square block of size BSZ x BSZ.
 
 Elements are stored in random order in the column.
 
 F. Nedelec, 17--27 March 2017
 */
class MatrixSparseSymmetricBlock : public Matrix
{
public:
    
    /// A BSZ x BSZ block element of the sparse matrix
    class Element
    {
    public:
        
        /// values of the element
        real val[BSZ*BSZ];
        
        union {
            /// index of the line
            index_type  inx;
            
            /// padding to align to 16 bytes boundaries
            real pad[BSZ];
        };
        
        Element() {}
        
        void reset(index_type i)
        {
            inx = i;
            for ( unsigned u = 0; u < BSZ*BSZ; ++u )
                val[u] = 0.0;
        }
        
        void scale(const real alpha)
        {
            for ( unsigned u = 0; u < BSZ*BSZ; ++u )
                val[u] *= alpha;
        }
        
        void add_sym(const real* mat)
        {
            for ( int x = 0; x < BSZ; ++x )
            {
                val[x+BSZ*x] += mat[x+BSZ*x];
                for ( int y = x+1; y < BSZ; ++y )
                {
                    val[y+BSZ*x] += mat[y+BSZ*x];
                    val[x+BSZ*y] += mat[y+BSZ*x];
                }
            }
        }
        
        void add_sym(const real alpha, const real* mat)
        {
            for ( int x = 0; x < BSZ; ++x )
            {
                val[x+BSZ*x] += alpha * mat[x+BSZ*x];
                for ( int y = x+1; y < BSZ; ++y )
                {
                    val[y+BSZ*x] += alpha * mat[y+BSZ*x];
                    val[x+BSZ*y] += alpha * mat[y+BSZ*x];
                }
            }
        }
        
        void sub_lower(const real* mat)
        {
            for ( int x = 0; x < BSZ; ++x )
            for ( int y = x; y < BSZ; ++y )
                val[y+BSZ*x] -= mat[y+BSZ*x];
        }

        void add(const real alpha, const real* mat)
        {
            //VecPrint::matPrint(std::clog, BSZ, BSZ, mat, BSZ);
            for ( int u = 0; u < BSZ*BSZ; ++u )
                val[u] += alpha * mat[u];
        }
 
        void add_lower(const real* mat)
        {
            for ( int x = 0; x < BSZ; ++x )
                for ( int y = x; y < BSZ; ++y )
                    val[y+BSZ*x] += mat[y+BSZ*x];
        }

        void add_lower(const real alpha, const real* mat)
        {
            for ( int x = 0; x < BSZ; ++x )
            for ( int y = x; y < BSZ; ++y )
                val[y+BSZ*x] += alpha * mat[y+BSZ*x];
        }

        real* addr(int i, int j)
        {
            return val + ( i + BSZ*j );
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
    
    /// col_next_[ii] is the index of the first non-empty column of index >= ii
    index_type * col_next_;

public:

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSparseSymmetricBlock();
    
    /// default destructor
    virtual ~MatrixSparseSymmetricBlock()  { deallocate(); }
    
    /// set all the element to zero
    void makeZero();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(unsigned alc);
    
    /// returns element at indices
    Element& block(index_type x, index_type y);
    
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
    
    /// multiplication of a vector: Y = Y + M * X with dim(X) = dim(M)
    void vecMulAddSSE(const real* X, real* Y) const;
    
    /// 2D isotropic multiplication of a vector: Y = Y + M * X with dim(X) = 2 * dim(M)
    void vecMulAddIso2D(const real* X, real* Y) const {};
    
    /// 3D isotropic multiplication of a vector: Y = Y + M * X with dim(X) = 3 * dim(M)
    void vecMulAddIso3D(const real* X, real* Y) const {};

    /// true if matrix is non-zero
    bool nonZero() const;
    
    /// number of elements which are not null
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


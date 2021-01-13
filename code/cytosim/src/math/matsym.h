// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATSYM_H
#define MATSYM_H

#include <cstdio>
#include "matrix.h"

/// A real symmetric Matrix
/**
the full upper triangular is stored
 */
class MatrixSymmetric : public Matrix
{
private:
    
    /// leading dimension of array
    unsigned int msLDD;
    
    /// size of memory which has been allocated
    unsigned int msAllocated;
    
    // full upper triangle:
    real* val;
    
    // if 'false', destructor will not call delete[] val;
    int in_charge;
    
public:

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSymmetric();
    
    
    /// constructor from an existing array
    MatrixSymmetric(int s, real* array)
    : Matrix(s)
    {
        msLDD = s;
        val = array;
        in_charge = 0;
    }

    /// constructor from an existing array
    MatrixSymmetric(int s, real* array, int ldd)
    : Matrix(s)
    {
        msLDD = ldd;
        val = array;
        in_charge = 0;
    }
    
    /// default destructor
    virtual ~MatrixSymmetric()  { deallocate(); }
    
    /// set all the element to zero
    void makeZero();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(unsigned alc);
    
    /// returns the address of element at (x, y), no allocation is done
    real* addr(index_type x, index_type y) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(index_type i, index_type j);
    
    /// scale the matrix by a scalar factor
    void scale(real a);
    
    /// multiplication of a vector: Y = Y + M * X, dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y) const;
    
    /// 2D isotropic multiplication of a vector: Y = Y + M * X
    void vecMulAddIso2D(const real* X, real* Y) const;
    
    /// 3D isotropic multiplication of a vector: Y = Y + M * X
    void vecMulAddIso3D(const real* X, real* Y) const;
    
    /// true if matrix is non-zero
    bool nonZero() const;
    
    /// number of element which are non-zero
    unsigned int  nbElements() const;
    
    /// returns a string which a description of the type of matrix
    std::string what() const;
};

#endif

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATRIX_H
#define MATRIX_H

#include "assert_macro.h"
#include <iostream>
#include "real.h"

/// The interface for all the large matrices
class Matrix
{
public:
    
    /// type of an index into the matrix
    typedef unsigned index_type;
    
protected:

    /// size of matrix
    unsigned   size_;
    
private:
    
    /// Disabled copy constructor (@todo: write copy constructor)
    Matrix(Matrix const&);
    
    /// Disabled copy assignment (@todo: write copy assignement)
    Matrix& operator=(Matrix const&);

public:
    
    /// constructor
    Matrix() { size_ = 0; }

    /// constructor
    Matrix(unsigned s) { size_ = s; }
    
    /// empty destructor
    virtual ~Matrix() {}
    
    /// return the size of the matrix
    unsigned size() const { return size_; }
    
    /// change the size of the matrix
    void resize(unsigned s) { allocate(s); size_=s; }

    //----------------------------------------------------------------------
    
    /// allocate the matrix to hold ( sz * sz ), all values may be lost
    virtual void allocate(unsigned alc) = 0;
        
    /// returns the address of element at (x, y), no allocation is done
    virtual real*  addr(index_type x, index_type y) const = 0;
    
    /// returns the address of element at (x, y), allocating if necessary
    virtual real&  operator()(index_type x, index_type y) = 0;
    
    /// returns the value of element at (x, y) or zero if not allocated
    real value(index_type x, index_type y) const;
    
    //----------------------------------------------------------------------
    
    /// set all the elements to zero
    virtual void makeZero() = 0;
    
    /// scale the matrix by a scalar factor
    virtual void scale(real) = 0;
    
    /// copy the block ( x, y, x+sx, y+sy ) into `mat`
    void copyBlock(real* mat, unsigned ldd, index_type sx, unsigned nx, index_type sy, unsigned ny) const;
    
    /// add the block ( x, x, x+sx, x+sx ) from this matrix to `mat`
    virtual void addDiagonalBlock(real* mat, unsigned ldd, index_type si, unsigned nb) const;
    
    /// add upper triangular half of ( idx, idx, idx+siz, idx+siz ) to `mat`
    virtual void addTriangularBlock(real* mat, unsigned ldd, index_type si, unsigned nb, unsigned dim) const;
    
    //----------------------------------------------------------------------
    
    /// Optional optimization to accelerate multiplications below
    virtual void prepareForMultiply() {}
    
    /// Vector multiplication: Y <- Y + M * X, size(X) = size(Y) = size(M)
    virtual void vecMulAdd(const real* X, real* Y) const = 0;
    
    /// Vector multiplication: Y <- M * X, size(X) = size(Y) = size(M)
    virtual void vecMul(const real* X, real* Y) const;

    /// isotropic vector multiplication: Y = Y + M * X, size(X) = size(Y) = 2 * size(M)
    virtual void vecMulAddIso2D(const real*, real*) const = 0;
    
    /// isotropic vector multiplication: Y = Y + M * X, size(X) = size(Y) = 3 * size(M)
    virtual void vecMulAddIso3D(const real*, real*) const = 0;
    
    //----------------------------------------------------------------------
    
    /// maximum absolute value considering all the elements
    virtual real maxNorm() const;
    
    /// true if the matrix is non-zero
    virtual bool nonZero() const;
    
    /// number of element which are non-zero
    virtual unsigned nbElements() const;
    
    /// returns a string which a description of the type of matrix
    virtual std::string what() const = 0;
    
    /// printf debug function in sparse mode: i, j : value
    virtual void printSparse(std::ostream &) const;
    
    /// printf debug function in full lines, all columns
    virtual void printFull(std::ostream &) const;
    
};


#endif

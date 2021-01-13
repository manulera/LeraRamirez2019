// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATRIX1_H
#define MATRIX1_H

#include "matrixbase.h"
#include "vector1.h"

/// 1x1 matrix
class Matrix1 : public MatrixBase<1>
{
public:

    ///constructor
    Matrix1() {};
    
    /// automatic conversion from MatrixBase<1>
    Matrix1(MatrixBase<1> const& m) : MatrixBase<1>(m) {}

    
    
    /// extract column vector at index `jj`
    Vector1 getColumn(const unsigned) const
    {
        return Vector1(val[0], 0);
    }
    
    /// set column vector
    void setColumn(const unsigned, Vector1 const& vec)
    {
        val[0] = vec.XX;
    }
    
    /// extract line vector
    Vector1 getLine(const unsigned) const
    {
        return Vector1(val[0], 0);
    }
    
    /// set line vector
    void setLine(const unsigned, Vector1 const& vec)
    {
        val[0] = vec.XX;
    }

    /// set matrix from one vector
    void setColumns(Vector1 const& vec)
    {
        val[0] = vec.XX;
    }
    

    /// rotation angle
    real rotationAngle() const;
    
    /// rotation from angle
    static Matrix1 rotationFromAngles(real a);

    /// return a rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix1 rotationToVector(const Vector1& vec);
    
    /// return a random rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix1 rotationToVector(const Vector1& vec, Random&)
    {
        return rotationToVector(vec);
    }
    
    /// a random rotation chosen uniformly
    static Matrix1 randomRotation(Random&);

};


inline Vector1 operator * (const Matrix1 & a, const Vector1 & v)
{
    return Vector1(a.val[0] * v.XX, 0);
}

#endif

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATRIX2_H
#define MATRIX2_H


#include "matrixbase.h"
#include "vector2.h"

/// 2x2 matrix
class Matrix2 : public MatrixBase<2>
{
public:
    
    ///constructor
    Matrix2() {};

    /// automatic conversion from MatrixBase<2>
    Matrix2(MatrixBase<2> const& m) : MatrixBase<2>(m) { }
    
    
    /// extract column vector at index `n`
    Vector2 getColumn(const unsigned n) const
    {
        if ( n == 0 )
            return Vector2(val[0], val[1]);
        else
            return Vector2(val[2], val[3]);
    }
    
    /// set column vector at index `n`
    void setColumn(const unsigned n, Vector2 const& vec)
    {
        if ( n == 0 )
        {
            val[0] = vec.XX;
            val[1] = vec.YY;
        }
        else
        {
            val[2] = vec.XX;
            val[3] = vec.YY;
        }
    }
    
    /// extract line vector at index `n`
    Vector2 getLine(const unsigned n) const
    {
        if ( n == 0 )
            return Vector2(val[0], val[2]);
        else
            return Vector2(val[1], val[3]);
    }
    
    /// set line vector at index `n`
    void setLine(const unsigned n, Vector2 const& vec)
    {
        if ( n == 0 )
        {
            val[0] = vec.XX;
            val[2] = vec.YY;
        }
        else
        {
            val[1] = vec.XX;
            val[3] = vec.YY;
        }
    }
    
    /// set matrix from one vector
    void setColumns(Vector2 const& vec1, Vector2 const& vec2)
    {
        val[0] = vec1.XX;
        val[1] = vec1.YY;
        
        val[2] = vec2.XX;
        val[3] = vec2.YY;
    }
    
    /// return the determinant of the matrix
    real determinant() const;
    
    /// return the inverse of the matrix
    Matrix2 inverted() const;

    /// rotation angle
    real rotationAngle() const;
    
    /// return rotation of given rotation angle
    static Matrix2 rotationFromAngles(real a);
    
    /// a rotation around the Z axis, of given angle
    static Matrix2 rotationAroundAxis(real dir, real angle)
    {
        if ( dir > 0 )
            return rotationFromAngles(angle);
        else
            return rotationFromAngles(-angle);
    }
    
    /// return a rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix2 rotationToVector(const Vector2& vec);
    
    /// return a random rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix2 rotationToVector(const Vector2& vec, Random&)
    {
        return rotationToVector(vec);
    }
    
    /// a random rotation chosen uniformly
    static Matrix2 randomRotation(Random&);

};


inline Vector2 operator * (const Matrix2 & a, const Vector2 & v)
{
    return Vector2( a.val[0] * v.XX + a.val[2] * v.YY,
                    a.val[1] * v.XX + a.val[3] * v.YY );
}

#endif

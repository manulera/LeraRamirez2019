// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATRIX3_H
#define MATRIX3_H

#include "matrixbase.h"
#include "vector3.h"

/// 3x3 matrix
class Matrix3 : public MatrixBase<3>
{
public:
    
    ///constructor
    Matrix3() {};

    /// automatic conversion from MatrixBase<3>
    Matrix3(MatrixBase<3> const& m) : MatrixBase<3>(m) { }
    
    
    /// extract column vector at index `n`
    Vector3 getColumn(const unsigned n) const
    {
        assert_true(n < 3);
        real const* v = val + n*3;
        return Vector3(v[0], v[1], v[2]);
    }
    
    /// set column vector at index `n`
    void setColumn(const unsigned n, Vector3 const& vec)
    {
        assert_true(n < 3);
        real * v = val + n*3;
        v[0] = vec.XX;
        v[1] = vec.YY;
        v[2] = vec.ZZ;
    }
    
    /// extract line vector at index `n`
    Vector3 getLine(const unsigned n) const
    {
        assert_true(n < 3);
        real const* v = val + n;
        return Vector3(v[0], v[3], v[6]);
    }
    
    /// set line vector at index `n`
    void setLine(const unsigned n, Vector3 const& vec)
    {
        assert_true(n < 3);
        real * v = val + n;
        v[0] = vec.XX;
        v[3] = vec.YY;
        v[6] = vec.ZZ;
    }
    
    /// set matrix from one vector
    void setColumns(Vector3 const& vec1, Vector3 const& vec2, Vector3 const& vec3)
    {
        val[0] = vec1.XX;
        val[1] = vec1.YY;
        val[2] = vec1.ZZ;
        
        val[3] = vec2.XX;
        val[4] = vec2.YY;
        val[5] = vec2.ZZ;
        
        val[6] = vec3.XX;
        val[7] = vec3.YY;
        val[8] = vec3.ZZ;
    }
    
    /// return the determinant of the matrix
    real determinant() const;
    
    /// return the inverse of the matrix
    Matrix3 inverted() const;

    /// rotation angle
    real rotationAngle() const;

    /// calculate rotation angle and Euler angle of axis
    void getEulerAngles(real& angle, real&, real&) const;
    
    /// 
    static Matrix3 rotationAroundAxisEuler(const real a[3]);

    /// return rotation of angle a, around axis of azimuth b and elevation c
    static Matrix3 rotationFromAngles(const real a[3]);

    
    /// a rotation around the X axis of specified angle
    static Matrix3 rotationAroundX(real angle);
    
    /// a rotation around the Y axis of specified angle
    static Matrix3 rotationAroundY(real angle);
    
    /// a rotation around the Z axis of specified angle
    static Matrix3 rotationAroundZ(real angle);
    
    /// a rotation around one the axis X if `x==0`, Y if `x==1` or Z if `x==2`
    static Matrix3 rotationAroundPrincipalAxis(unsigned x, real angle);
    
    /// a rotation around the given axis
    static Matrix3 rotationAroundAxis(const Vector3&, real angle);    
    
    /// return a rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix3 rotationToVector(const Vector3&);
    
    /// return a random rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    /**
     In 3D, this rotation is chosen uniformly among all the rotation transforming (1,0,0) into `vec`.
     The function will fail if ( vec == 0 ).
     */
    static Matrix3 rotationToVector(const Vector3&, Random&);
    
    /// a random rotation chosen uniformly
    static Matrix3 randomRotation(Random&);
    
    /// a rotation of angle 'angle' around an axis chosen randomly
    static Matrix3 randomRotation(real angle, Random&);

};


inline Vector3 operator * (const Matrix3 & a, const Vector3 & v)
{
    return Vector3( a.val[0]*v.XX + a.val[3]*v.YY + a.val[6]*v.ZZ,
                    a.val[1]*v.XX + a.val[4]*v.YY + a.val[7]*v.ZZ,
                    a.val[2]*v.XX + a.val[5]*v.YY + a.val[8]*v.ZZ );
}


#endif

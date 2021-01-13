// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ISOMETRY_H
#define ISOMETRY_H

#include "dim.h"


#if (DIM==1)

   #include "matrix1.h"
   typedef Matrix1 MatrixD;

#elif (DIM==2)

   #include "matrix2.h"
   typedef Matrix2 MatrixD;

#elif (DIM==3)

   #include "matrix3.h"
   typedef Matrix3 MatrixD;

#else

    #error 'DIM' is not properly defined!

#endif


/// A Rotation is a matrix of dimension DIM x DIM
typedef MatrixD Rotation;


/// An affine transformation in space.
/**
 A Isometry contains a vector T and a rotation matrix M,
 and represents the affine transformation:
 @code
 X -> M.X + T
 @endcode
 */
class Isometry
{
public:
    
    /// translation component
    Vector  vec;
    
    /// rotation component
    MatrixD rot;
    
public:
    
    Isometry()
    {
        vec.zero();
        rot.makeOne();
    }

    Isometry(Vector const& v)
    {
        vec = v;
        rot.makeOne();
    }

    Isometry(Vector const& v, MatrixD const& r)
    {
        vec = v;
        rot = r;
    }

    void zero()
    {
        vec.zero();
        rot.makeOne();
    }

    Vector const& translation() const
    {
        return vec;
    }
    
    MatrixD const& rotation() const
    {
        return rot;
    }
    
    void translate(Vector const& v)
    {
        vec += v;
    }
   
    void rotate(MatrixD const& mat)
    {
        rot = mat * rot;
        vec = mat * vec;
    }

    void multiply(Isometry const& iso)
    {
        vec = iso.rot * vec + iso.vec;
        rot = iso.rot * rot;
    }
};


#endif

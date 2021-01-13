// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef VECTOR3_H
#define VECTOR3_H


#include "real.h"
#include "assert_macro.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

/// Vector3 is a vector with 3 `real` components.
/**
 Note: We assume that the coordinates XX, YY and ZZ are adjacent in memory,
 allowing easy conversion operators to and from C-array.
 Although this is not guaranteed by the C-standard, this is usually the case.
 */
class Vector3
{
    
public:
    
    /// dimensionality is 3
    static unsigned dimensionality() { return 3; }
    
    /// coordinates are public
    real XX, YY, ZZ;
    
    
    /// by default, coordinates are not initialized
    Vector3() {}
    
    /// construct from values
    Vector3(const real x, const real y, const real z) : XX(x), YY(y), ZZ(z) {}
    
    /// construct from address
    Vector3(const real v[]) : XX(v[0]), YY(v[1]), ZZ(v[2]) {}

    /// construct from array of size d
    Vector3(const real v[], int d)
    {
        if ( d > 2 )
        {
            XX = v[0];
            YY = v[1];
            ZZ = v[2];
        }
        else
        {
            if ( d > 1 )
            {
                XX = v[0];
                YY = v[1];
            }
            else
            {
                if ( d > 0 )
                    XX = v[0];
                else
                    XX = 0;
                YY = 0;
            }
            ZZ = 0;
        }
    }

    /// create new Vector with coordinates from the given array
    static Vector3 make(const real b[]) { return Vector3(b[0], b[1], b[2]); }

    /// destructor (is not-virtual: do not derive from this class)
    ~Vector3() {}
    
    
    /// address of coordinate array
    real * data()                { return &XX; }
    
    /// constant address of coordinate array
    real const* data()     const { return &XX; }
    
#if ( 1 )    
    /// implicit conversion to a modifiable real pointer
    operator real*()             { return &XX; }
    
    /// implicit conversion to a constant real pointer
    operator const real*() const { return &XX; }
#else
    /// value of a coordinate
    real operator[](unsigned ii) const
    {
        assert_true(ii<3);
        return (&XX)[ii];
    }
    
    /// modifiable access to individual coordinates
    real& operator[](int ii)
    {
        assert_true(ii<3);
        return (&XX)[ii];
    }
#endif
    
    /// copy coordinates from array of size d
    void get(const real v[], const int& d)
    {
        if ( d > 2 )
        {
            XX = v[0];
            YY = v[1];
            ZZ = v[2];
        }
        else
        {
            if ( d > 1 )
            {
                XX = v[0];
                YY = v[1];
            }
            else
            {
                if ( d > 0 )
                    XX = v[0];
                else
                    XX = 0;
                YY = 0;
            }
            ZZ = 0;
        }
    }
    
    /// replace coordinates by the ones provided
    void get(const float b[])
    {
        XX = b[0];
        YY = b[1];
        ZZ = b[2];
    }
    
    /// replace coordinates by the ones provided
    void get(const double b[])
    {
        XX = b[0];
        YY = b[1];
        ZZ = b[2];
    }
    
    /// copy coordinates to given array
    void put(float b[]) const
    {
        b[0] = XX;
        b[1] = YY;
        b[2] = ZZ;
    }
    
    /// copy coordinates to given array
    void put(double b[]) const
    {
        b[0] = XX;
        b[1] = YY;
        b[2] = ZZ;
    }
    
    /// add content to given address
    void add_to(real b[]) const
    {
        b[0] += XX;
        b[1] += YY;
        b[2] += ZZ;
    }
    
    /// add content scaled by `alpha` to given address
    void add_to(real alpha, real b[]) const
    {
        b[0] += alpha * XX;
        b[1] += alpha * YY;
        b[2] += alpha * ZZ;
    }

    /// add content `n` times to array `b` of size `ldd*n`
    void add_to(real b[], int n, int ldd) const
    {
        for ( int i = 0; i < n; ++i )
        {
            b[ldd*i  ] += XX;
            b[ldd*i+1] += YY;
            b[ldd*i+2] += ZZ;
        }
    }

    /// subtract to given address
    void sub_to(real b[]) const
    {
        b[0] -= XX;
        b[1] -= YY;
        b[2] -= ZZ;
    }
    
    /// add content scaled by `alpha` to given address
    void sub_to(real alpha, real b[]) const
    {
        b[0] -= alpha * XX;
        b[1] -= alpha * YY;
        b[2] -= alpha * ZZ;
    }
    
    /// set coordinates to zero
    void zero()
    {
        XX = 0;
        YY = 0;
        ZZ = 0;
    }

    /// change coordinates
    void set(const real x, const real y, const real z)
    {
        XX = x;
        YY = y;
        ZZ = z;
    }

    /// change signs of all coordinates
    void oppose()
    {
        XX = -XX;
        YY = -YY;
        ZZ = -ZZ;
    }
    
    //------------------------------------------------------------------
    
    /// the square of the standard norm
    real normSqr() const
    {
        return XX*XX + YY*YY + ZZ*ZZ;
    }

    
    /// the standard norm = sqrt(x^2+y^2+z^2)
    real norm() const
    {
        return sqrt(XX*XX + YY*YY + ZZ*ZZ);
    }

    /// the 2D norm = sqrt(x^2+y^2)
    real normXY() const
    {
        return sqrt(XX*XX + YY*YY);
    }

    /// the 2D norm = sqrt(x^2+z^2)
    real normXZ() const
    {
        return sqrt(XX*XX + ZZ*ZZ);
    }
    
    /// the 2D norm = sqrt(y^2+z^2)
    real normYZ() const
    {
        return sqrt(YY*YY + ZZ*ZZ);
    }

    /// square of the distance to other point == (a-this).normSqr()
    real distanceSqr(Vector3 const& a) const
    {
        return (a.XX-XX)*(a.XX-XX) + (a.YY-YY)*(a.YY-YY) + (a.ZZ-ZZ)*(a.ZZ-ZZ);
    }
    
    /// distance to other point == (a-this).norm()
    real distance(Vector3 const& a) const
    {
        return sqrt((a.XX-XX)*(a.XX-XX) + (a.YY-YY)*(a.YY-YY) + (a.ZZ-ZZ)*(a.ZZ-ZZ));
    }
    
    /// returns  min(x, y, z)
    real minimum() const
    {
        if ( XX > YY )
            return ( YY > ZZ ? ZZ : YY );
        else
            return ( ZZ > XX ? XX : ZZ );
    }
    
    /// returns  max(x, y, z)
    real maximum() const
    {
        if ( XX > YY )
            return ( XX > ZZ ? XX : ZZ );
        else
            return ( ZZ > YY ? ZZ : YY );
    }
    
    /// the infinite norm = max(|x|, |y|, |z|)
    real norm_inf() const
    {
        real x = fabs(XX);
        real y = fabs(YY);
        real z = fabs(ZZ);
        if ( x > y )
            return ( x > z ? x : z );
        else
            return ( z > y ? z : y );
    }
    
    /// true if no component is NaN
    bool valid() const
    {
        return ( XX == XX ) && ( YY == YY ) && ( ZZ == ZZ );
    }
    
    /// true if all components are zero
    bool null() const
    {
        return ( XX == 0 ) && ( YY == 0 ) && ( ZZ == 0 );
    }

    /// normalize to norm=n
    void normalize(const real n = 1.0)
    {
        real s = n / norm();
        XX *= s;
        YY *= s;
        ZZ *= s;
    }
    
    /// returns the colinear vector of norm=n
    const Vector3 normalized(const real n = 1.0) const
    {
        real sn = n / norm();
        return Vector3(XX*sn, YY*sn, ZZ*sn);
    }

    //------------------------------------------------------------------
    
    /// returns a perpendicular vector, of comparable but unspecified norm
    const Vector3 orthogonal() const
    {
        if ( fabs(XX) < fabs(YY) )
        {
            if ( fabs(XX) < fabs(ZZ) )
                return Vector3(0.0, -ZZ,  YY); //XX is the smallest
            else
                return Vector3( YY, -XX, 0.0); //ZZ is the smallest
        }
        else
        {
            if ( fabs(YY) < fabs(ZZ) )
                return Vector3(-ZZ, 0.0,  XX); //YY is the smallest
            else
                return Vector3( YY, -XX, 0.0); //ZZ is the smallest
        }
    }

    /// returns a perpendicular vector, of norm `n`
    const Vector3 orthogonal(const real n) const
    {
        if ( fabs(XX) < fabs(YY) )
        {
            if ( fabs(XX) < fabs(ZZ) )
            {
                // XX is the smallest component
                real s = n / sqrt(YY*YY+ZZ*ZZ);
                return Vector3(0.0, -s*ZZ, s*YY);
            }
            else
            {
                // ZZ is the smallest component
                real s = n / sqrt(XX*XX+YY*YY);
                return Vector3(s*YY, -s*XX, 0.0);
            }
        }
        else
        {
            if ( fabs(YY) < fabs(ZZ) )
            {
                // YY is the smallest component
                real s = n / sqrt(XX*XX+ZZ*ZZ);
                return Vector3(-s*ZZ, 0.0, s*XX);
            }
            else
            {
                // ZZ is the smallest component
                real s = n / sqrt(XX*XX+YY*YY);
                return Vector3(s*YY, -s*XX, 0.0);
            }
        }
    }
    
    /// returns a vector perpendicular to *this, close to `d`, and of norm = `n`
    const Vector3 orthogonal(Vector3 d, const real n)
    {
        d -= ( ( d * (*this) ) / normSqr() ) * (*this);
        return d.normalized(n);
    }

    /**
     Set vectors 'ex' and 'ey' to build an orthonormal basis (this, ex, ey),
     assuming that 'norm(*this) == 1'
     
     From `Building an Orthonormal Basis, Revisited`,
     Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
     */
    void orthonormal(Vector3& ex, Vector3& ey) const
    {
#ifndef NDEBUG
        if ( normSqr() < 0.9 )
        {
            // this should not happen...
            ex = orthogonal(1);
            ey = cross(*this, ex).normalized(1);
            std::clog << "rescued orthonormal(" << toString() << ")\n";
        }
#endif
#ifdef REAL_IS_FLOAT
        real s = copysignf(1.0f, ZZ);
#else
        real s = copysign(1.0, ZZ);
#endif
#if ( 1 )
        /// optimized version by Marc B. Reynolds
        const real a = YY / ( ZZ + s );
        const real b = YY * a;
        const real c = XX * a;
        // below normSqr(ex) = normSqr(this) + a*a*(normSqr(this)-s*s)
        ex.set(-ZZ - b, c, XX);
        ey.set(s * c, s * b - 1.0, s * YY);
        //if you do not mind an inverted basis, use  ey.set(c, b-s, YY);
#else
        /// original code from Duff et al.
        const real a = -1.0 / ( ZZ + s );
        const real b = XX * YY * a;
        // below normSqr(ex) = 1 + x*x*a*a*(normSqr(this)-s*s)
        ex.set(1.0 + s * XX * XX * a, s * b, -s * XX);
        ey.set(b, s + YY * YY * a, -YY);
#endif
    }
    
    
    /// return unit vector obtained by rotating `n` around `*this`, by angle defined by cosinus and sinus
    /**
     The result is a Vector orthogonal to *this
     */
    const Vector3 rotate(Vector3 const& n, real cs, real sn)
    {
        //Set two orthogonal vector to 'd' in arbitrary directions
        Vector3 e, f;
        orthonormal(e, f);
        // compute coordinates of n in reference frame (e, f):
        real ne = n * e;
        real nf = n * f;
        // normalization factor:
        real s = 1.0 / sqrt( ne * ne + nf * nf );
        // rotate 'n' by an angle defined by cs and sn:
        return (( ne * cs + nf * sn ) * s ) * e + (( nf * cs - ne * sn ) * s ) * f;
    }

    /// convert from cartesian to spherical coordinates ( r, theta, phi )
    const Vector3 spherical() const
    {
        return Vector3(sqrt(XX*XX+YY*YY+ZZ*ZZ),
                       atan2(YY, XX),
                       atan2(sqrt(XX*XX+YY*YY), ZZ));
    }
    
    /// convert from spherical to cartesian coordinates ( x, y, z )
    const Vector3 cartesian() const
    {
        return Vector3(XX*cos(YY)*sin(ZZ),
                       XX*sin(YY)*sin(ZZ),
                       XX*cos(ZZ));
    }
    
    //------------------------------------------------------------------
    
    /// addition of two vectors
    friend const Vector3 operator +(Vector3 const& a, Vector3 const& b)
    {
        return Vector3(a.XX+b.XX, a.YY+b.YY, a.ZZ+b.ZZ);
    }
    
    /// subtraction of two vectors
    friend const Vector3 operator -(Vector3 const& a, Vector3 const& b)
    {
        return Vector3(a.XX-b.XX, a.YY-b.YY, a.ZZ-b.ZZ);
    }
    
    /// unary + operator does nothing
    friend const Vector3 operator +(Vector3 const& b)
    {
        return b;
    }
    
    /// opposition of a vector
    friend const Vector3 operator -(Vector3 const& b)
    {
        return Vector3(-b.XX, -b.YY, -b.ZZ);
    }
    
    /// returns the element-by-element product
    const Vector3 e_mul(const real b[]) const
    {
        return Vector3(XX*b[0], YY*b[1], ZZ*b[2]);
    }

    /// returns the element-by-element division
    const Vector3 e_div(const real b[]) const
    {
        return Vector3(XX/b[0], YY/b[1], ZZ/b[2]);
    }
    
    /// returns a vector with each element squared
    const Vector3 e_squared() const
    {
        return Vector3(XX*XX, YY*YY, ZZ*ZZ);
    }
    
    /// returns sum of all coordinates
    real e_sum() const
    {
        return XX + YY + ZZ;
    }

    
    /// cross product of two vectors
    friend const Vector3 cross(Vector3 const& a, Vector3 const& b)
    {
        return Vector3(a.YY * b.ZZ - a.ZZ * b.YY,
                       a.ZZ * b.XX - a.XX * b.ZZ,
                       a.XX * b.YY - a.YY * b.XX);
    }

    
    /// scalar product of two vectors
    friend real operator *(Vector3 const& a, Vector3 const& b)
    {
        return a.XX*b.XX + a.YY*b.YY + a.ZZ*b.ZZ;
    }
    
    /// multiplication by scalar s
    friend const Vector3 operator *(Vector3 const& a, const real s)
    {
        return Vector3(s*a.XX, s*a.YY, s*a.ZZ);
    }
    
    /// mutiplication by scalar s
    friend const Vector3 operator *(const real s, Vector3 const& a)
    {
        return Vector3(s*a.XX, s*a.YY, s*a.ZZ);
    }
    
    /// division by scalar s
    friend const Vector3 operator /(Vector3 const& a, const real s)
    {
        return Vector3(a.XX/s, a.YY/s, a.ZZ/s);
    }
    
    /// addition of another vector b
    void operator +=(Vector3 const& b)
    {
        XX += b.XX;
        YY += b.YY;
        ZZ += b.ZZ;
    }
    
    /// subtraction of another vector b
    void operator -=(Vector3 const& b)
    {
        XX -= b.XX;
        YY -= b.YY;
        ZZ -= b.ZZ;
    }
    
    /// multiplication by a scalar
    void operator *=(const real b)
    {
        XX *= b;
        YY *= b;
        ZZ *= b;
    }
    
    /// division by a scalar
    void operator /=(const real b)
    {
        XX /= b;
        YY /= b;
        ZZ /= b;
    }
    
    //------------------------------------------------------------------
    
    /// equality test
    friend bool operator ==(Vector3 const& a, Vector3 const& b)
    {
        return ( a.XX==b.XX  &&  a.YY==b.YY  &&  a.ZZ==b.ZZ );
    }
    
    /// non-equality test
    friend bool operator !=(Vector3 const& a, Vector3 const& b)
    {
        return ( a.XX!=b.XX  ||  a.YY!=b.YY  ||  a.ZZ!=b.ZZ );
    }
    
    //------------------------------------------------------------------
    
    /// conversion to a string
    std::string toString() const
    {
        std::ostringstream oss;
        oss << XX << " " << YY << " " << ZZ;
        return oss.str();
    }
    
    /// conversion to ASCII string with given precision
    std::string toString(int w, int p) const
    {
        std::ostringstream oss;
        oss.precision(p);
        oss << std::setw(w) << XX << " ";
        oss << std::setw(w) << YY << " ";
        oss << std::setw(w) << ZZ;
        return oss.str();
    }

    /// print to a file
    void print(FILE * out = stdout) const
    {
        fprintf(out, "  %+9.3f %+9.3f %+9.3f", XX, YY, ZZ);
    }
    
    /// print to a file, surrounded by parenthesis
    void pprint(FILE * out = stdout) const
    {
        fprintf(out, "( %+9.3f %+9.3f %+9.3f )", XX, YY, ZZ);
    }

    /// print, followed by a new line
    void println(FILE * out = stdout) const
    {
        fprintf(out, "  %+9.3f %+9.3f %+9.3f\n", XX, YY, ZZ);
    }

    //------------------------------------------------------------------
    
    /// add a random component in [-s, s] to each coordinate
    void addRand(real s);
    
    
    /// a vector orthogonal to *this, with `norm == n, chosen randomly and uniformly
    const Vector3 randOrthoU(real n) const;
    
    /// a vector orthogonal to *this, with `norm <= n`, chosen randomly and uniformly
    const Vector3 randOrthoB(real n) const;
    
    
    /// Vector with random independent coordinates in [0,+1]
    static const Vector3 prand();
    
    /// Vector with random independent coordinates in [0,+n]
    static const Vector3 prand(real n);
    
    /// Vector with random independent coordinates in [-1,+1]
    static const Vector3 srand();
    
    /// Vector with random independent coordinates in [-n,+n]
    static const Vector3 srand(real n);
   
    
    /// random Vector of norm = 1; sampling is uniform
    static const Vector3 randU();
    
    /// return a random vector of norm = n; sampling is uniform
    static const Vector3 randU(real n);
    
    
    /// return a random vector of norm <= 1; sampling is uniform
    static const Vector3 randB();
    
    /// return a random vector of norm <= n; sampling is uniform
    static const Vector3 randB(real n);
    
    
    /// return a random vector with Normally distributed coordinates ~ N(0,n)
    static const Vector3 randG(real n);
    
};


//-------------------------- associated global functions -----------------------

/// stream input operator
std::istream & operator >> (std::istream&, Vector3&);

/// stream output operator
std::ostream & operator << (std::ostream&, Vector3 const&);

/// linear interpolation: returns a + x * b
const Vector3 interpolate(Vector3 const& a, real x, Vector3 const& b);

/// square of the distance between two points, equivalent to (a-b).normSqr()
real distanceSqr(Vector3 const& a, Vector3 const& b);

/// distance between two points, equivalent to (a-b).norm()
real distance(Vector3 const& a, Vector3 const& b);


#endif


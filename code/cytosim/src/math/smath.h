// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 Some basic mathematical functions
 Francois Nedelec, 
*/

#ifndef SMATH_H
#define SMATH_H

#include "real.h"
#include <cmath>
#include <sstream>
#include <stdint.h>
#include <iomanip>


#ifndef M_PI
/// Ratio of a circle's circumference to its diameter
const real M_PI=3.14159265358979323846264338327950288;
#endif

#ifndef M_E
const real M_E=2.7182818284590452354;
#endif


/// simple mathematical functions, mostly templated
namespace sMath
{
    // limit `x` inside [`a` , `b` ]:
    template <typename T>
    inline const T& constrain(T& x, const T& a, const T& b)
    {
        if ( x < a )
            return a;
        if ( x > b )
            return b;
        return x;
    }
    
    /// minimum of 2 values
    template <typename T> 
    inline const T& min(const T& a, const T& b)
    {
        return ( a > b ) ? b : a;
    }

    /// maximum of 2 values
    template <typename T> 
    inline const T& max(const T& a, const T& b)
    {
        return ( a > b ) ? a : b;
    }
    
    /// minimum of three arguments
    template <typename T>
    inline const T& min(const T& a, const T& b, const T& c)
    {
        if ( a > b )
            return ( b < c ) ? b : c;
        else
            return ( a < c ) ? a : c;
    }
    
    /// maximum of three arguments
    template <typename T>
    inline const T& max(const T& a, const T& b, const T& c)
    {
        if ( a > b )
            return ( a > c ) ? a : c;
        else
            return ( b > c ) ? b : c;
    }
    
    /// minimum of four arguments
    template <typename T>
    inline const T& min(const T& a, const T& b, const T& c, const T& d)
    {
        return min(min(a,b), min(c,d));
    }
    
    /// maximum of four arguments
    template <typename T>
    inline const T& max(const T& a, const T& b, const T& c, const T& d)
    {
        return max(max(a,b), max(c,d));
    }
    
    /// find which of three arguments is the smallest
    template <typename T>
    inline int arg_min(const T& a, const T& b, const T& c)
    {
        if ( a > b )
            return ( b > c ) ? 2 : 1;
        else
            return ( a > c ) ? 2 : 0;
    }
    
    /// find which of three arguments is the largest
    template <typename T>
    inline int arg_max(const T& a, const T& b, const T& c)
    {
        if ( a < b )
            return ( b < c ) ? 2 : 1;
        else
            return ( a < c ) ? 2 : 0;
    }

    /// find which of four arguments is the smallest
    template <typename T>
    inline int arg_min(const T& a, const T& b, const T& c, const T& d)
    {
        if ( a > b )
        {
            // consider ( b, c, d )
            if ( b > c )
                return ( c > d ) ? 3 : 2;
            else
                return ( b > d ) ? 3 : 1;
        }
        else
        {
            // consider ( a, c, d )
            if ( a > c )
                return ( c > d ) ? 3 : 2;
            else
                return ( a > d ) ? 3 : 0;
        }
    }
    
    /// find which of four arguments is the smallest
    template <typename T>
    inline int arg_max(const T& a, const T& b, const T& c, const T& d)
    {
        if ( a < b )
        {
            // consider ( b, c, d )
            if ( b < c )
                return ( c < d ) ? 3 : 2;
            else
                return ( b < d ) ? 3 : 1;
        }
        else
        {
            // consider ( a, c, d )
            if ( a < c )
                return ( c < d ) ? 3 : 2;
            else
                return ( a < d ) ? 3 : 0;
        }
    }


    /**
     Set vectors 'x' and 'y' to make an orthonormal basis (x, y, z)
     
     Building an Orthonormal Basis, Revisited
     Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
     */
    inline void orthonormal(const float z[3], float x[3], float y[3])
    {
        const float s = copysignf(1.0f, z[2]);
#if ( 1 )
        /// optimized version by Marc B. Reynolds
        const float a = z[1] / ( z[2] + s );
        const float b = z[1] * a;
        const float c = z[0] * a;
        x[0] = -z[2] - b;
        x[1] = c;
        x[2] = z[0];
        y[0] = s * c;
        y[1] = s * b - 1.0f;
        y[2] = s * z[1];
#else
        const float a = -1.0f / ( a[2] + s );
        const float b = a[0] * a[1] * a;
        x[0] = 1.0 + s * z[0] * z[0] * a;
        x[1] = s * b;
        x[2] = -s * z[0];
        y[0] = b;
        y[1] = s + z[1] * z[1] * a;
        y[2] = -z[1];
#endif
    }
    
    /**
     Set vectors 'x' and 'y' to make an orthonormal basis (x, y, z)
     
     Building an Orthonormal Basis, Revisited
     Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
     */
    inline void orthonormal(const double z[3], double x[3], double y[3])
    {
        const double s = copysign(1.0, z[2]);
#if ( 1 )
        /// optimized version by Marc B. Reynolds
        const double a = z[1] / ( z[2] + s );
        const double b = z[1] * a;
        const double c = z[0] * a;
        x[0] = -z[2] - b;
        x[1] = c;
        x[2] = z[0];
        y[0] = s * c;
        y[1] = s * b - 1.0;
        y[2] = s * z[1];
#else
        const double a = -1.0 / ( a[2] + s );
        const double b = a[0] * a[1] * a;
        x[0] = 1.0 + s * z[0] * z[0] * a;
        x[1] = s * b;
        x[2] = -s * z[0];
        y[0] = b;
        y[1] = s + z[1] * z[1] * a;
        y[2] = -z[1];
#endif
    }

    
#if ( 0 )
    /// the sign of a number
    template <typename T> 
    inline int sign(const T& a)
    {
        return a > T(0) ? 1 : ( a < T(0) ? -1 : 0 );
    }
    
    /// Composes a floating point value with the magnitude of x and the sign of y.
    template <typename T>
    inline T copysign(const T& x, const T& y)
    {
        if ( y > T(0) )
        {
            if ( x > T(0) )
                return x;
            else
                return -x;
        }
        else
        {
            if ( x > T(0) )
                return -x;
            else
                return x;
        }
    }
#endif
    
#if ( 0 )
    /// bit-hack absolute value
    inline float absf(float a)
    {
        union { uint32_t u; float f; } tmp;
        tmp.f = a;
        tmp.u &= 0x7FFFFFFFUL;
        return tmp.f;
    }
    
    /// bit-hack absolute value
    inline double absf(double a)
    {
        union { uint64_t u; double f; } tmp;
        tmp.f = a;
        tmp.u &= 0x7FFFFFFFFFFFFFFFULL;
        return tmp.f;
    }
#endif
    
    /// square of a number
    template <typename T> 
    inline T square(const T& a)
    {
        return a * a;
    }
    
    /// cube of a number
    template <typename T>
    inline T cube(const T& a)
    {
        return a * a * a;
    }
    
    /// power of `a` by positive integer exponent `n`
    /** This should be equivalent to std::pow(a, n) */
    template <typename T>
    inline T power_int(const T& a, unsigned n)
    {
        T x = a;
        T y = 1;
        while ( n )
        {
            if ( n & 1 )
                y = y * x;
            x = x * x;
            n = n >> 1;
        }
        return y;
    }
    
    
    ///power of `a` by integer exponent `n`
    template <typename T>
    inline T power(const T& a, const int n)
    {
        if ( n < 0 )
            return power_int(1.0/a, -n);
        return power_int(a, n);
    }
    
    
    ///power of `a` by integer exponent `n`
    template <typename T>
    T nextPowerOf2(T k)
    {
        if ( k & (k-1) )
        {
            do
                k &= k-1;
            while ( k & (k-1) );
            
            k <<= 1;
        }
        return k;
    }
    
    
    ///square of distance between two vectors in dimension `dim`
    template <typename T>
    inline T  distanceSqr(const T a[], const T b[], int dim)
    {
        T n = ( a[0] - b[0] ) * ( a[0] - b[0] );
        for( int d = 1; d < dim; ++d )
            n += ( a[d] - b[d] ) * ( a[d] - b[d] );
        return n;
    }
 
    ///usual distance between two vectors of dimension `dim`
    template <typename T>
    inline T  distance(const T a[], const T b[], int dim)
    {
        T n = ( a[0] - b[0] ) * ( a[0] - b[0] );
        for( int d = 1; d < dim; ++d )
            n += ( a[d] - b[d] ) * ( a[d] - b[d] );
        return sqrt(n);
    }
    
    /// return the usual base-10 representation of a number
    /** Note that with C++11, you can use std::to_string */
    template <typename T> 
    std::string repr(T const& x)
    {
        std::ostringstream oss;
        oss << x;
        return oss.str();
    }
    
    template <typename T>
    std::string repr(T const& x, unsigned width, unsigned precision)
    {
        std::ostringstream oss;
        oss << std::setw(width) << std::setprecision(precision) << std::fixed << x;
        return oss.str();
    }
    
    //------------------------------------------------------------------------------
#pragma mark -
    
    /// used for periodic boundary conditions:
    inline void fold(real& x, const real p)
    {
        while ( x >  p ) x -= p+p;
        while ( x < -p ) x += p+p;
    }

#ifdef WIN32
    
    //this is needed under windows:
    inline real remainder( const real a, const real b )
    {
        int p = (int)floor( 0.5 + a / b );
        if ( p )
            return a - p * b;
        else
            return a;
    }
    
    inline real round(real x)
    {
        if ( x < 0 )
            return -floor(0.5-x);
        else
            return  floor(0.5+x);
    }

#endif
    
    //------------------------------------------------------------------------------
#pragma mark -
    
    ///extract a 10-decimal digit form a number:
    /** 1st digit is really the first one, we do not start at zero */
    template <typename T> 
    inline int digit(T x, const int p)
    {
        for ( int q=1; q < p; ++q )
            x /= 10;
        return x % 10;
    }
    
    ///copy bytes
    inline void copyBytes( void * dest, const void * src, const unsigned cnt)
    {
        for ( unsigned ii=0; ii < cnt; ++ii )
            ((char*)dest)[ii] = ((char*)src)[ii];
    }
    

    //------------------------------------------------------------------------------
    
    /// return smallest power of 2 that is greater or equal to `x`
    inline unsigned next_power(unsigned x)
    {
        if ( x > 0 )
            return 0;
        --x;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        return x+1;
    }
    
    /// return smallest power of 2 that is greater or equal to `x`
    inline size_t next_power(size_t x)
    {
        if ( x == 0 )
            return 0;
        --x;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        x |= x >> 32;
        return x+1;
    }

    /// number of '1' bits in a 32-bits integer (Charlie Gordon & Don Clugston)
    /** Should use Intel SIMD instruction POPCNT */
    inline unsigned int count_bits(uint32_t v)
    {
        v = v - ((v >> 1) & 0x55555555);
        v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
        return (((v + (v >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
    }
    
    
    /// number of '1' bits, from: http://graphics.stanford.edu/~seander/bithacks.html
    /** Works up to 128 bits */
    template <typename T>
    unsigned int count_bits2(T v)
    {
        v = v - ((v >> 1) & (T)~(T)0/3);
        v = (v & (T)~(T)0/15*3) + ((v >> 2) & (T)~(T)0/15*3);
        v = (v + (v >> 4)) & (T)~(T)0/255*15;
        return (T)(v * ((T)~(T)0/255)) >> (sizeof(v) - 1) * 8;
    }
   
}


#endif //#ifdef SMATH_H

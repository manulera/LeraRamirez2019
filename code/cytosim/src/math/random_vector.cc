// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "vector.h"
#include "random.h"


/**
 Random vectors are generated using the global Random Generator `RNG`
 */
extern Random RNG;


//------------------------------------------------------------------------------
#pragma mark - 1D

const Vector1 Vector1::srand()        { return Vector1(  RNG.sreal(), 0); }
const Vector1 Vector1::srand(real n)  { return Vector1(n*RNG.sreal(), 0); }
const Vector1 Vector1::prand()        { return Vector1(  RNG.preal(), 0); }
const Vector1 Vector1::prand(real n)  { return Vector1(n*RNG.preal(), 0); }
const Vector1 Vector1::randU()        { return Vector1(  RNG.sflip(), 0); }
const Vector1 Vector1::randU(real n)  { return Vector1(n*RNG.sflip(), 0); }
void  Vector1::addRand(real n)        { XX += n*RNG.sreal(); }

const Vector1 Vector1::randB()        { return Vector1(   RNG.sreal(), 0); }
const Vector1 Vector1::randB(real n)  { return Vector1( n*RNG.sreal(), 0); }
const Vector1 Vector1::randG(real n)  { return Vector1( n*RNG.gauss(), 0);  }

const Vector1 Vector1::randOrthoU(real len) const { return Vector1(0, 0); }
const Vector1 Vector1::randOrthoB(real len) const { return Vector1(0, 0); }


//------------------------------------------------------------------------------
#pragma mark - 2D

const Vector2 Vector2::srand()        { return Vector2(  RNG.sreal(),   RNG.sreal()); }
const Vector2 Vector2::srand(real n)  { return Vector2(n*RNG.sreal(), n*RNG.sreal()); }
const Vector2 Vector2::prand()        { return Vector2(  RNG.preal(),   RNG.preal()); }
const Vector2 Vector2::prand(real n)  { return Vector2(n*RNG.preal(), n*RNG.preal()); }
const Vector2 Vector2::randG(real n)  { return Vector2(n*RNG.gauss(), n*RNG.gauss()); }
void  Vector2::addRand(real n)        { XX += n*RNG.sreal(); YY += n*RNG.sreal(); }

#if ( 0 )

const Vector2 Vector2::randU()
{
    real d, x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = x*x + y*y;
    } while ( d > 1.0  ||  d < 0.01 );
    d = sqrt( d );
    return Vector2( x/d, y/d );
}

const Vector2 Vector2::randU(const real n)
{
    real d, x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = x*x + y*y;
    } while ( d > 1.0  ||  d < 0.01 );
    d = n / sqrt( d );
    return Vector2( x*d, y*d );
}

#else

const Vector2 Vector2::randU()
{
    real d, x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = x*x + y*y;
    } while ( d > 1.0 );
    if ( d > 0 )
        return Vector2((x*x-y*y)/d, 2*x*y/d);
    else
        return Vector2(1, 0);
}

const Vector2 Vector2::randU(const real n)
{
    real d, x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = x*x + y*y;
    } while ( d > 1.0 );
    if ( d > 0 )
        return Vector2(n*(x*x-y*y)/d, n*2*x*y/d);
    else
        return Vector2(n, 0);
}

#endif


const Vector2 Vector2::randB()
{
    real x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
    } while ( x*x + y*y > 1.0 );
    return Vector2( x, y );
}


const Vector2 Vector2::randB(const real n)
{
    real x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
    } while ( x*x + y*y > 1.0 );
    return Vector2( x*n, y*n );
}


const Vector2 Vector2::randOrthoU(const real len) const
{
    real s = RNG.sflip() * len / sqrt( XX * XX + YY * YY );
    return Vector2(-s*YY, s*XX, 0);
}


const Vector2 Vector2::randOrthoB(const real len) const
{
    real s = RNG.sreal() * len / sqrt( XX * XX + YY * YY );
    return Vector2(-s*YY, s*XX, 0);
}

//------------------------------------------------------------------------------
#pragma mark - 3D

const Vector3 Vector3::srand()        { return Vector3(RNG.sreal(),     RNG.sreal(),   RNG.sreal()); }
const Vector3 Vector3::srand(real n)  { return Vector3(n*RNG.sreal(), n*RNG.sreal(), n*RNG.sreal()); }
const Vector3 Vector3::prand()        { return Vector3(RNG.preal(),     RNG.preal(),   RNG.preal()); }
const Vector3 Vector3::prand(real n)  { return Vector3(n*RNG.preal(), n*RNG.preal(), n*RNG.preal()); }
const Vector3 Vector3::randG(real n)  { return Vector3(n*RNG.gauss(), n*RNG.gauss(), n*RNG.gauss()); }
void  Vector3::addRand(real n)        { XX += n*RNG.sreal(); YY += n*RNG.sreal(); ZZ += n*RNG.sreal(); }


#if ( 0 )

/// hypercube rejection method
const Vector3 Vector3::randU()
{
    real x, y, z, d;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        z = RNG.sreal();
        d = x*x + y*y + z*z;
    } while ( d > 1.0  ||  d < 0.01 );
    d = 1.0 / sqrt( d );
    return Vector3( x*d, y*d, z*d );
}

/// hypercube rejection method
const Vector3 Vector3::randU(real n)
{
    real x, y, z, d;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        z = RNG.sreal();
        d = x*x + y*y + z*z;
    } while ( d > 1.0  ||  d < 0.01 );
    d = n / sqrt( d );
    return Vector3( x*d, y*d, z*d );
}

#elif ( 1 )

/**
 Derived from Marsaglia (1972)
 Allen & Tildesley "Computer Simulation of Liquids" Clarendon Pres, Oxford 1987
 http://mathworld.wolfram.com/SpherePointPicking.html
 This uses only 2 random-numbers!
*/
const Vector3 Vector3::randU()
{
    real x, y, d;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = 1.0 - x*x - y*y;
    } while ( d <= 0 );
    real h = 2 * sqrt(d);
    return Vector3( x*h, y*h, d+d-1.0 );
}

const Vector3 Vector3::randU(const real n)
{
    real x, y, d;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = 1.0 - x*x - y*y;
    } while ( d <= 0 );
    real h = ( n + n ) * sqrt(d);
    return Vector3( x*h, y*h, n*(d+d-1.0) );
}

#else


/**
 From Cook (1957)
 http://mathworld.wolfram.com/SpherePointPicking.html
 This uses 4 random-numbers, but avoids the square-root
 */
const Vector3 Vector3::randU()
{
    real x, y, z, t, d;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        z = RNG.sreal();
        t = RNG.sreal();
        d = x*x + y*y + z*z + t*t;
    } while ( d > 1.0 );
    return Vector3(2*(y*t+x*z)/d, 2*(z*t-x*y)/d, (x*x+t*t-y*y-z*z)/d);
}

const Vector3 Vector3::randU(const real n)
{
    real x, y, z, t, d;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        z = RNG.sreal();
        t = RNG.sreal();
        d = x*x + y*y + z*z + t*t;
    } while ( d > 1.0 );
    return Vector3(n*2*(y*t+x*z)/d, n*2*(z*t-x*y)/d, n*(x*x+t*t-y*y-z*z)/d);
}

#endif


const Vector3 Vector3::randB()
{
    real x, y, z;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        z = RNG.sreal();
    } while ( x*x + y*y + z*z > 1.0 );
    return Vector3( x, y, z );
}


const Vector3 Vector3::randB(const real n)
{
    real x, y, z;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        z = RNG.sreal();
    } while ( x*x + y*y + z*z > 1.0 );
    return Vector3( x*n, y*n, z*n );
}


#if ( 1 )

const Vector3 Vector3::randOrthoU(const real len) const
{
    real n = normSqr();
    if ( n > REAL_EPSILON )
    {
        Vector2 d = Vector2::randU();
        Vector3 x, y, z = *this / sqrt(n);
        z.orthonormal(x, y);
        return x * ( len * d.XX ) + y * ( len * d.YY );
    }
    return randU(len);
}

#else

/**
 This method is less efficient
 */
const Vector3 Vector3::randOrthoU(const real len) const
{
    Vector2 d = Vector2::randU();
    Vector3 b = orthogonal(1);
    Vector3 c = cross(*this, b).normalized();
    return b * ( len * d.XX ) + c * ( len * d.YY );
}

#endif

const Vector3 Vector3::randOrthoB(const real len) const
{
    real n = normSqr();
    if ( n > REAL_EPSILON )
    {
        Vector2 d = Vector2::randB();
        Vector3 x, y, z = *this / sqrt(n);
        z.orthonormal(x, y);
        return x * ( len * d.XX ) + y * ( len * d.YY );
    }
    return randU(len);
}


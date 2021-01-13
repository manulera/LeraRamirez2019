// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "random.h"

#include <climits>
#include <sys/time.h>
#include <cstring>
#include <ctime>


/// RNG = Random Number Generator
Random RNG;

//------------------------------------------------------------------------------
/**
 The generator is initialized with a zero state vector,
 and seed() or seedTimer() must be called before any
 random number can be produced.
 */
Random::Random()
: sfmt_buffer(&sfmt.state[0].u[0])
{
    if ( sizeof(uint32_t) != 4 )
    {
        fprintf(stderr, "Random can only work if sizeof(uint32_t) == 4\n");
        exit(1);
    }
    
    // initialize pointers:
    sfmt_ptr = sfmt_buffer;
    
    // zero state vector:
    memset(sfmt_buffer, 0, 4*SFMT_N32);
    gauss_ptr = gauss_buffer;
}


Random::~Random()
{
}


/**
 Get a uint32_t from t and c
 Better than uint32_t(x) in case x is floating point in [0,1]
 Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)
 */
uint32_t hash(time_t t, clock_t c)
{
    static uint32_t differ = 0;  // guarantee time-based seeds will change
    
    uint32_t h1 = 0;
    unsigned char* p = (unsigned char*) &t;
    for ( size_t i = 0; i < sizeof(t); ++i )
    {
        h1 *= UCHAR_MAX + 2U;
        h1 += p[i];
    }
    uint32_t h2 = 0;
    p = (unsigned char*) &c;
    for ( size_t j = 0; j < sizeof(c); ++j )
    {
        h2 *= UCHAR_MAX + 2U;
        h2 += p[j];
    }
    return ( h1 + differ++ ) ^ h2;
}


uint32_t Random::seedTimer()
{
    struct timeval now;
    gettimeofday(&now, 0);
    uint32_t s = hash( now.tv_sec, now.tv_usec );
    seed(s);
    return s;
}


bool Random::seeded() const
{
    for ( unsigned n = 0; n < SFMT_N32; ++n )
        if ( sfmt_buffer[n] )
            return true;
    
    return false;
}

//------------------------------------------------------------------------------
#pragma mark -

int Random::sign_exc(const real a)
{
    if ( a <= 0 )
    {
        if ( a < 0 )
            return -1;
        return RNG.sflip();
    }
    return 1;
}


float Random::pfloat()
{
    //This assumes IEEE Standard 754 Floating point numbers
    //32 bits: 1 for sign, 8 for exponents, 23 for fraction
    union { uint32_t i; float f; } tmp;
    tmp.i = URAND32();
    uint32_t E = 126;
    while (( tmp.i < BIT31 ) && ( E > 94 ))
    {
        tmp.i <<= 1; --E;
    }
    tmp.i = (( tmp.i << 1 ) >> 9 ) | ( E << 23 );
    return tmp.f;
}


float Random::sfloat()
{
    //This assumes IEEE Standard 754 Floating point numbers
    //32 bits: 1 for sign, 8 for exponents, 23 for fraction
    union { uint32_t i; float f; } tmp;
    tmp.i = URAND32();
    bool S = tmp.i & 1;
    uint32_t E = 126;
    while (( tmp.i < BIT31 ) && ( E > 94 ))
    {
        tmp.i <<= 1; --E;
    }
    tmp.i = (( tmp.i << 1 ) >> 9 ) | ( E << 23 );
    return S ? tmp.f : -tmp.f;
}

//------------------------------------------------------------------------------
#pragma mark -

double Random::pdouble()
{
    //This assumes IEEE Standard 754 Floating point numbers
    //64 bits: 1 for sign, 11 for exponents, 52 for Fraction
    //set all the 64 bits to random, using two random uint32_t:
    uint64_t l = URAND64();
    uint64_t E = 1022;
    while ( l < BIT63  &&  E > 959 )
    {
        l <<= 1; --E;
    }
    union { uint64_t i; double f; } tmp;
    tmp.i = ((l >> 11) & 0x000FFFFFFFFFFFFFULL) | ( E << 52 );
    return tmp.f;
}


double Random::sdouble()
{
    //This assumes IEEE Standard 754 Floating point numbers
    //64 bits: 1 for sign, 11 for exponents, 52 for Fraction
    //set all the 64 bits to random, using two random uint32_t:
    uint64_t l = URAND64();
    bool     S = l & 1;  //sign-bit
    uint64_t E = 1022;
    while ( l < BIT63  &&  E > 959 )
    {
        l <<= 1; --E;
    }
    union { uint64_t i; double f; } tmp;
    tmp.i = ((l >> 11) & 0x000FFFFFFFFFFFFFULL) | ( E << 52 );
    return S ? tmp.f : -tmp.f;
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 Fill array `vec[]` with Gaussian values ~ N(0,1).
 n_vec is the size of `vec` should be a multiple of 2.
 Return the number of values that were stored in `vec`
 */
size_t Random::gauss_refill(real vec[], const size_t n_vec)
{
    int32_t tmp[SFMT_N32];
    int32_t const*const t_end = tmp + SFMT_N32;
    int32_t *t = tmp;

    sfmt_fill_array32(&sfmt, (uint32_t*)tmp, SFMT_N32);
    
    real *v = vec;
    real *const v_end = vec + n_vec;

    while ( t < t_end )
    {
        real x = t[0] * 0x1p-31;
        real y = t[1] * 0x1p-31;
        real w = x * x + y * y;
        if ( w <= 1 )
        {
            w = sqrt( -2 * log(w) / w );
            v[0] = w * x;
            v[1] = w * y;
            v += 2;
            if ( v >= v_end )
                break;
        }
        t += 2;
    }
    return v - vec;
}

/**
 Fill array `gauss_buffer` with Gaussian values ~ N(0,1).
 Set `gauss_ptr` past the last position containing a valid number.
 The number of gaussian values set by this function is random,
 and it may even be zero.
 */
void Random::gauss_refill()
{
    int32_t tmp[SFMT_N32];
    int32_t const*const t_end = tmp + SFMT_N32;
    int32_t *t = tmp;

    sfmt_fill_array32(&sfmt, (uint32_t*)tmp, SFMT_N32);

    gauss_ptr = gauss_buffer;

    while ( t < t_end )
    {
        real x = t[0] * 0x1p-31;
        real y = t[1] * 0x1p-31;
        real w = x * x + y * y;
        if ( w <= 1 )
        {
            // could use fast reverse square root (_mm_rsqrt14_pd)
            w = sqrt( -2 * log(w) / w );
            gauss_ptr[0] = w * x;
            gauss_ptr[1] = w * y;
            gauss_ptr += 2;
        }
        t += 2;
    }
}


/**
 Signed real number, following a normal law N(0,v*v)
 using the polar rejection method (cf. Numerical Recipe)
 */
void Random::gauss_set(real & a, real & b, const real& v)
{
    real x, y, w;
    do {
        x = sreal();
        y = sreal();
        w = x * x + y * y;
    } while ( w >= 1.0 || w == 0 );
    /*
     formula below are only valid if ( w > 0 ),
     which may be false only with a minuscule probability
     */
    w = v * sqrt( -2 * log(w) / w );
    a = w * x;
    b = w * y;
}


#if ( 0 )

/**
 Fill `n` Gaussian values ~ N(0,1) in array `vec[]`.
 */
void Random::gauss_set(real vec[], size_t n_vec, const real& v = 1.0)
{
    unsigned u = n_vec % 8;
    unsigned w = u % 2;
    
    if ( w )
        vec[0] = v * gauss();
    
    for ( ; w < u; w += 2 )
        gauss_set(vec[w], vec[w+1], v);
    
    for ( ; u < n_vec; u += 8 )
    {
        gauss_set(vec[u  ], vec[u+1], v);
        gauss_set(vec[u+2], vec[u+3], v);
        gauss_set(vec[u+4], vec[u+5], v);
        gauss_set(vec[u+6], vec[u+7], v);
    }
}

#else

/**
 Fill `n` Gaussian values ~ N(0,1) in array `vec[]`.
 */
void Random::gauss_set(real vec[], size_t n_vec)
{
    real const* vec_end = vec + n_vec;

    do {
        // check if `gauss_buffer` provides enough values to complete `vec`
        if ( gauss_ptr - gauss_buffer > vec_end - vec )
        {
            while ( vec < vec_end )
            {
                --gauss_ptr;
                //assert(gauss_ptr >= gauss_buffer);
                *vec = *gauss_ptr;
                ++vec;
            }
            return;
        }
        
        // use all remaining `gauss_buffer`
        while ( --gauss_ptr >= gauss_buffer )
        {
            //assert(vec < vec_end);
            *vec = *gauss_ptr;
            ++vec;
        }
        
        gauss_refill();
        
    } while ( vec < vec_end );
}


/**
 Fill `n` Gaussian values ~ N(0,v) in array `vec[]`.
 */
void Random::gauss_set(real vec[], size_t n_vec, const real& v)
{
    real const* vec_end = vec + n_vec;
    
    do {
        // check if `gauss_buffer` provides enough values to complete `vec`
        if ( gauss_ptr - gauss_buffer > vec_end - vec )
        {
            while ( vec < vec_end )
            {
                --gauss_ptr;
                //assert(gauss_ptr >= gauss_buffer);
                *vec = v * (*gauss_ptr);
                ++vec;
            }
            return;
        }
        
        // use all remaining `gauss_buffer`
        while ( --gauss_ptr >= gauss_buffer )
        {
            //assert(vec < vec_end);
            *vec = v * (*gauss_ptr);
            ++vec;
        }
        
        gauss_refill();
        
    } while ( vec < vec_end );
}

#endif


/**
 this version uses cos() and sin() and is slower than gauss().
 const real PI = 3.14159265358979323846264338327950288;
 */
void Random::gauss_slow(real& x, real& y)
{
    real angle = real( URAND32() ) * 1.46291807926715968105133780430979e-9;
    //the constant is 2*pi/2^32
    real norm  = sqrt( -2 * log( preal_exc() ));
    x = norm * cos(angle);
    y = norm * sin(angle);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 integer in [0,n] for n < 2^32
 */
uint32_t Random::pint_slow(const uint32_t n)
{
    // Find which bits are used in n
    uint32_t used = n | ( n >> 1 );
    used |= (used >> 2);
    used |= (used >> 4);
    used |= (used >> 8);
    used |= (used >> 16);
    
    // Draw numbers until one is found in [0,n]
    uint32_t i;
    do
        i = URAND32() & used;  // toss unused bits to shorten search
    while ( i > n );
    return i;
}


/**
 returns a random integer with exactly `b` bits equal to `1`,
 but randomly positionned.
 */
uint32_t Random::pint_bits(int b)
{
    uint32_t n = 0;
    if ( b < 16 )
    {
        while ( b > 0 )
        {
            uint32_t x = 1 << ( URAND32() % 32 );
            if (!( n & x ))
            {
                n += x;
                --b;
            }
        }
    }
    else
    {
        n = ~0;
        while ( b < 32 )
        {
            uint32_t x = 1 << ( URAND32() % 32 );
            if ( n & x )
            {
                n -= x;
                ++b;
            }
        }
    }
    return n;
}


/**
 returns an integer in [0 n], with the ratios given in the array of ints
 */
uint32_t Random::pint_ratio(const uint32_t n, const int ratio[])
{
    int sum = 0;
    uint32_t ii;
    for ( ii = 0; ii < n; ++ii )
        sum += ratio[ii];
    // sum==0 may denotes a careless use of the function, with wrong arguments.
    // it might be safer to throw an exception
    if ( sum == 0 ) return 0;
    sum = (int) floor( preal() * sum );
    ii = 0;
    while ( sum >= ratio[ii] )
        sum -= ratio[ii++];
    return ii;
}

//------------------------------------------------------------------------------
#pragma mark -
/**
 Return Poisson distributed integer, with expectation=E  variance=E
 http://en.wikipedia.org/wiki/Poisson_distribution

 This routine is slow for large values of E.
 If E > 256, this returns a Gaussian distribution of parameter (E, E),
 which is a good approximation of the Poisson distribution
 
 Knuth D.E. The art of computer programming, Vol II: Seminumerical algorithms.

 This method fails for E > 700, in double precision
*/
uint32_t Random::poisson_knuth(const real E)
{
    if ( E > 256 )
        return static_cast<uint32_t>( gauss() * sqrt(E) + E );
    if ( E < 0 )
        return 0;
    real L = exp(-E);
    real p = preal();
    uint32_t k = 0;
    while ( p > L )
    {
        ++k;
        p *= preal();
    }
    return k;
}


/**
 Return Poisson distributed integer, with expectation=E  variance=E
 http://en.wikipedia.org/wiki/Poisson_distribution
 
 This routine is slow for large values of E.
 If E > 512, this returs a Gaussian distribution of parameter (E, E),
 which is a good approximation of the Poisson distribution.
 
 This method fails for E > 700, in double precision
 */
uint32_t Random::poisson(const real E)
{
    if ( E > 256 )
        return static_cast<uint32_t>( gauss() * sqrt(E) + E );
    if ( E < 0 )
        return 0;
    real p = exp(-E);
    real s = p;
    uint32_t k = 0;
    real u = preal();
    while ( u > s )
    {
        ++k;
        p *= E / k;
        s += p;
    }
    return k;
}


/**
 This is equivalent to calling poisson(exp(-E))
 The argument is EL = exp(-E)
 expectation=E  variance=E (see wikipedia, Poisson Distribution)
 */
uint32_t Random::poissonE(const real EL)
{
    real p = preal();
    uint32_t k = 0;
    while ( p > EL )
    {
        ++k;
        p *= preal();
    }
    return k;
}


uint32_t Random::geometric(const real P)
{
    if ( P < 0 )
        return 0;
    uint32_t pi = (uint32_t)( P * 0x1p32 );
    
    uint32_t s = 0;
    while ( URAND32() > pi )
        ++s;
    return s;
}


uint32_t Random::binomial(const int N, const real P)
{
    if ( P < 0 )
        return 0;
    uint32_t pi = (uint32_t)( P * 0x1p32 );
    
    uint32_t s = 0;
    for ( int x = 0; x < N; ++x )
        if ( URAND32() < pi )
            ++s;
    return s;
}


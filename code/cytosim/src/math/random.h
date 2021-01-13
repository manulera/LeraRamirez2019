// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef RANDOM_H
#define RANDOM_H

#include <stdint.h>
#include <cmath>
#include "real.h"

#define SFMT_MEXP 19937
#include "SFMT.h"


/// Random Number Generator
/**
 The generation of random bits is done with Mersenne Twister from U. of Hiroshima
 
 http://en.wikipedia.org/wiki/Mersenne_twister
 http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 
 This class provides convenient functions to generate floating points numbers,
 following different distributions, and other convenient functions.
*/
class Random
{
private:

    /// Mersenne Twister state variables
    sfmt_t    sfmt;
    
    /// alias to start of the Mersenne Twister data vector
    uint32_t *const sfmt_buffer;
    
    /// pointer to access current value
    uint32_t const* sfmt_ptr;
    
    /// extract next Mersenne Twister 32 random bits
    inline uint32_t URAND32()
    {
        if ( sfmt_ptr <= sfmt_buffer )
        {
            sfmt_gen_rand_all(&sfmt);
            sfmt_ptr = sfmt_buffer + SFMT_N32;
        }
        --sfmt_ptr;
        return *sfmt_ptr;
    }
    
    /// extract next Mersenne Twister 32 random bits
    inline int32_t RAND32()
    {
        return static_cast<int32_t>(URAND32());
    }
    
    inline uint64_t URAND64()
    {
        sfmt_ptr -= 2;
        if ( sfmt_ptr < sfmt_buffer )
        {
            sfmt_gen_rand_all(&sfmt);
            sfmt_ptr = sfmt_buffer + SFMT_N32 - 2;
        }
        return *reinterpret_cast<uint64_t const*>(sfmt_ptr);
    }
    
    inline int64_t RAND64()
    {
        return static_cast<int64_t>(URAND64());
    }

private:
    
    /// the last bit in a 32-bits integer
    static const uint32_t BIT31 = 1U << 31, FRAC32 = 0x7FFFFFU, EXPON32 = 127U << 23;
    
    /// exponent for a double precision float
    static const uint64_t BIT63 = 1ULL << 63, EXPON64 = 1023ULL << 52;
    
    /// used by gauss()
    real    gauss_buffer[SFMT_N32];
    
    real *  gauss_ptr;
    
    /// refill array `gauss_buffer` with normal law N(0,1). Set gauss_ptr
    void  gauss_refill();

public:
            
    /// Constructor sets the state vector to zero
    Random();
    
    /// destructor
    ~Random();
    
    uint32_t const* buffer() const { return sfmt_buffer; }
    
    void      refill()
    {
        sfmt_gen_rand_all(&sfmt);
        sfmt_ptr = sfmt_buffer + SFMT_N32;
    }
    
    /// true if state vector is not entirely zero
    bool      seeded() const;

    /// seed with integer
    void      seed(const uint32_t s)
    {
        sfmt_init_gen_rand(&sfmt, s);
        sfmt_gen_rand_all(&sfmt);
        sfmt_ptr = sfmt_buffer + SFMT_N32;
    }
    
    /// seed with time()
    uint32_t  seedTimer();

    /// unsigned integer in [0,2^32-1]
    uint32_t  pint()                     { return URAND32(); }
    
    /// unsigned integer in [0,n-1] for n < 2^32
    uint32_t  pint(const uint32_t n)     { return uint32_t(URAND32()*0x1p-32*n); }
    
    /// unsigned integer in [0,n] for n < 2^32
    uint32_t  pint_inc(const uint32_t n) { return uint32_t(URAND32()*0x1p-32*(n+1)); }
    
    /// integer in [0,n] for n < 2^32, (slow) integer based algorithm
    uint32_t  pint_slow(const uint32_t n);
    
    /// a random unsigned integer with exactly `b` bit equal to `1`
    uint32_t  pint_bits(int b);
    
    /// signed integer in [-2^31+1,2^31-1]; inlined for speed
    int32_t   sint()                     { return RAND32(); }
    
    /// integer in [-N, N], boundaries included
    int32_t   sint_inc(const int32_t n)  { return pint( 2*n+1 ) - n; }
    
    /// integer in [1-N, N-1], i.e. in ]-N,N[ with boundaries excluded
    int32_t   sint(const int32_t n)      { return pint( 2*n-1 ) - n + 1; }
    
    /// random integer in [low, high]  ( = low + pint(1+high-low) )
    int32_t   int_range(const int32_t low, const int32_t high)
    {
        if ( high >= low )
            return low + pint( 1 + high - low );
        else
            return high + pint( 1 + low - high );
    }

    /// integer in [0 N], with probabilities given in ratio[] of size N, with sum(ratio)>0
    uint32_t  pint_ratio(uint32_t n, const int ratio[]);

    /// integer k of probability distribution p(k,E) = exp(-E) * pow(E,k) / factorial(k)
    uint32_t  poisson(real E);
    
    /// integer k of probability distribution p(k,E) = EL * pow(E,k) / factorial(k)
    uint32_t  poissonE(real EL);
    
    /// integer k of probability distribution p(k,E) = exp(-E) * pow(E,k) / factorial(k)
    uint32_t  poisson_knuth(real E);

    /// number of successive unsuccessful trials, when success has probability p (result >= 0)
    uint32_t  geometric(real p);

    /// number of sucesses among n trials of probability p
    uint32_t  binomial(int n, real p);
    
    
#if ( 1 )  // 0 = enable crazy optimization
    /// returns true with probability (p), and false with probability (1-p)
    bool      test(const real& p)        { return ( URAND32() <  p * 0x1p32 ); }
    /// returns true with probability (1-p), and false with probability (p)
    bool      test_not(const real& p)    { return ( URAND32() >= p * 0x1p32 ); }
#else
    /// returns true with probability p, and false with probability (1-p)
    bool      test(const real& p)
    {
        //this bithack might be faster, using only 23 bits, but fails if p < ~1/2^23:
        union { uint32_t i; float f; } tmp;
        tmp.i = ( URAND32() >> 9 ) | EXPON32;
        return ( tmp.f < p + 1 );
    }
    bool      test_not(const real& p)
    {
        //this bithack might be faster, using only 23 bits, but fails if p < ~1/2^23:
        union { uint32_t i; float f; } tmp;
        tmp.i = ( URAND32() >> 9 ) | EXPON32;
        return ( tmp.f > p + 1 );
    }
#endif
    
    /// true with probability p / 2^32
    bool      test_uint(uint32_t p)      { return URAND32() < p; }
    
    /// True  or  False  with equal chance
    bool      flip()                     { return URAND32() & 1024U; } //any bit could be used
    
    /// returns -1  or  1 with equal chance
    int       sflip()                    { return URAND32() & 1024U ? -1 : 1; }

    /// True with probability 1/4
    bool      flip_4th()                 { return URAND32() < 1<<30; }
  
    /// return the sign of `a` if `a != 0` and -1 or +1 randomly, if `a == 0`
    int       sign_exc(const real a);
    
    /// fast (dirty) random float in [0,1[, requires IEEE Standard 754
    float     pfloat23()
    {
        //by setting random bits for the fraction-bits of a float IEEE 754, 
        //we get a random number between 1 and 2. We substract 1.0,
        //but that drops the lower bits, reducing the precision
        union { uint32_t i; float f; } tmp;
        tmp.i = EXPON32 | ( URAND32() >> 9 );
        return tmp.f - 1.0f;
    }
    
    /// random float in [0,1[, requires IEEE Standard 754 
    float     pfloat();
    
    /// random float in ]-1,1[, requires IEEE Standard 754
    float     sfloat();
    
    /// slow random double in [0,1[, using two uint32_t to set all the fraction bits, requires IEEE Standard 754
    double    pdouble();
    
    /// slow random double in ]-1,1[, using two uint32_t to set all the fraction bits, requires IEEE Standard 754
    double    sdouble();
    
    /// positive real number in [0,1[, zero included
    real      preal()                    { return URAND32() * 0x1p-32; }
    
    /// positive real number in [0,n[ = n * preal() : deprecated, use preal() * n
    //real      preal(const real n)        { return n * ( URAND32() * 0x1p-32 ); }
    
    /// signed real number in ]-1,1[, boundaries excluded
    real      sreal()                    { return RAND32() * 0x1p-31; }
    
    /// signed real number in ]-1/2, 1/2[, boundaries excluded
    real      sreal_half()               { return RAND32() * 0x1p-32; }

    /// non-zero real number in ]0,1]
    real      preal_exc()                { return URAND32() * 0x1p-32 + 0x1p-32; }
    
    
    /// non-zero real number in ]0,n]
    real      preal_exc(real n)          { return preal_exc() * n; }
    
    /// real number uniformly distributed in [a,b[
    real      real_uniform(real a, real b) { return a + preal() * ( b - a ); }
    
    
    /// random Gaussian number, following a normal law N(0,1)
    real      gauss()
    {
        while ( gauss_ptr <= gauss_buffer )
            gauss_refill();
        --gauss_ptr;
        return *gauss_ptr;
    }

    /// set two independent random numbers, both following a normal law N(0,v*v)
    void      gauss_set(real &, real &, const real& v);

    /// fill array `vec` with independent random numbers following normal law N(0,1).
    void      gauss_set(real vec[], size_t n);
    
    /// fill array `vec` with independent random numbers following normal law N(0,v*v).
    void      gauss_set(real vec[], size_t n, const real& v);
    
    /// fill array `vec` with normal law N(0,1). Return number of values set
    size_t    gauss_refill(real vec[], size_t n);

    /// signed real number, following a normal law N(0,1), slower algorithm
    void      gauss_slow(real &, real&);
    
    /// random in [0, inf[, with P(x) = exp(-x), mean = 1.0, variance = 1.0
    real      exponential() { return -log( URAND32() * 0x1p-32 + 0x1p-32 );  }
    
    /// exponentially distributed positive real, with P(x) = exp(-x/E) / E,  parameter E is 1/Rate
    real      exponential(const real E) { return -E * log( URAND32() * 0x1p-32 + 0x1p-32 );  }

    /// uniform choice among two given values
    template<typename T>
    T         choice(const T& x, const T& y)
    {
        if ( flip() )
            return x;
        else
            return y;
    }
    
    /// uniform choice within an array of `size` values
    template<typename T>
    T         choice(const T val[], int size)
    {
        return val[ pint(size) ];
    }
    
    /// uniform shuffling of array `T[]`.
    /** Algorithm from knuth's The Art of Programming, Vol 2 chp. 3.4.2 */
    template <typename T> 
    void mix(T val[], int size)
    {
        int  jj = size, kk;
        while ( jj > 1 )
        {
            kk = URAND32() % jj;
            --jj;
            T tmp   = val[jj];
            val[jj] = val[kk];
            val[kk] = tmp;
        }
    }
        
};



/**
 Linear congruential random number generator
 The coefficients are found in Numerical Recipes 3rd Ed. Chapter 7
 See also http://en.wikipedia.org/wiki/Linear_congruential_generator
 
 The low-order bits should never be relied on for any degree of randomness whatsoever. 
 To use the topmost bit:
 @code 
   if ( z & 0x80000000U )
 @endcode
 */

inline uint32_t lcrng1(uint32_t z) { return z * 2024337845U + 797082193U; }

inline uint32_t lcrng2(uint32_t z) { return z * 279470273U + 4294967291U; }

inline uint32_t lcrng3(uint32_t z) { return z * 1372383749U + 1289706101U; }

inline uint32_t lcrng4(uint32_t z) { return z * 1103515245U + 12345U; }

inline uint32_t lcrng3(uint32_t z, int n)
{
    uint32_t r = lcrng3(z);
    for ( int c = 0; c < n; ++c )
        r = lcrng3(r);
    return r;
}

#endif  //RANDOM_H

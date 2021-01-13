// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/*
 Tests Intel's Streaming SIMD
 F. Nedelec, HD, July 2013
 
 To compile: c++ -O4 tictoc.cc -mavx test.cc
 To generate assembly: c++ -S test.cc
 */

#include <cstdio>
#include "tictoc.h"
#include <immintrin.h>

typedef __m128d vec2;
#define SSE(x) _mm_##x##_pd

double sum(vec2 const& v)
{
    double * s = (double*)(&v);
    return s[0]+s[1];
}

void print(vec2 const& v, char const* x)
{
    double * s = (double*)(&v);
    printf("vec2 %s ( %5.2f %5.2f )\n", x, s[0], s[1]);
}


#ifdef __AVX__

typedef __m256  vecf;
#define AVS(x) _mm256_##x##_ps

typedef __m256d vec4;
#define AVX(x) _mm256_##x##_pd

double sum(vec4 const& v)
{
    double * s = (double*)(&v);
    return s[0]+s[1]+s[2]+s[3];
}

void print(vec4 const& v, char const* x)
{
    double * s = (double*)(&v);
    printf("vec4 %s ( %5.2f %5.2f %5.2f %5.2f )\n", x, s[0], s[1], s[2], s[3]);
}

void print(vecf const &v, char const* x)
{
    float * s = (float*)(&v);
    printf("vecf %s ( %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f )\n", x,
           s[7], s[6], s[5], s[4], s[3], s[2], s[1], s[0]);
}

#endif

//------------------------------------------------------------------------------

typedef double real;

const unsigned size = 1<<10;
real x[size], y[size];

void init()
{
    for ( unsigned ii=0; ii<size; ++ii )
    {
        x[ii] = 1.0/real(size-ii);
        y[ii] = real(size-ii);
    }
}


real scalar()
{
    real d = 0;
    for ( unsigned ii=0; ii<size; ++ii )
        d += x[ii] * y[ii];
    return d;
}

real vector2()
{
    vec2 s = SSE(setzero)();
    for ( unsigned ii=0; ii<size; ii+=2 )
        s = SSE(add)(s, SSE(mul)( SSE(load)(x+ii), SSE(load)(y+ii) ));
    _mm_empty();
    
    return sum(s);
}

#ifdef __AVX__


real vector4()
{
    vec4 s = AVX(setzero)();
    for ( unsigned ii=0; ii<size; ii+=4 )
        s = AVX(add)(s, AVX(mul)( AVX(load)(x+ii), AVX(load)(y+ii) ));
    _mm_empty();
    
    return sum(s);
}


real vectorU()
{
    vec4 v0 = AVX(setzero)();
    vec4 v1 = AVX(setzero)();
    vec4 v2 = AVX(setzero)();
    vec4 v3 = AVX(setzero)();
    
    for ( unsigned ii=0; ii<size; ii+=16 )
    {
        v0 = AVX(add)(v0, AVX(mul)( AVX(load)(x+ii   ), AVX(load)(y+ii   ) ));
        v1 = AVX(add)(v1, AVX(mul)( AVX(load)(x+ii+4 ), AVX(load)(y+ii+4 ) ));
        v2 = AVX(add)(v2, AVX(mul)( AVX(load)(x+ii+8 ), AVX(load)(y+ii+8 ) ));
        v3 = AVX(add)(v3, AVX(mul)( AVX(load)(x+ii+12), AVX(load)(y+ii+12) ));
    }
    
    vec4 s = AVX(add)(AVX(add)(v0, v1), AVX(add)(v2, v3));
    _mm_empty();
    
    return sum(s);
}

#else

#warning Unsupported Intel SIMD intruction set AVX
real vector4()
{
    return 0;
}

real vectorU()
{
    return 0;
}

#endif

void run(real (*func)(), const char name[])
{
    const int rep = 1<<14;
    real a = 0, b = 0, c = 0, d = 0;
    real e = 0, f = 0, g = 0, h = 0;
    fprintf(stderr, "%s:  ", name);
    TicToc::tic();
    for ( int ii=0; ii<rep; ++ii )
    {
        a = (*func)();
        b = (*func)();
        c = (*func)();
        d = (*func)();
        e = (*func)();
        f = (*func)();
        g = (*func)();
        h = (*func)();
    }
    
    real s = a + b + c + d + e + f + g + h;
    double ms = TicToc::toc();
    fprintf(stderr, " %f :  %.0f ms\n", s, ms);
}



void test_swap0()
{
    vec2 a = SSE(setr)(0, 1);
    vec2 b = SSE(setr)(2, 3);
    print(a, "a");
    print(b, "b");
    
    print(SSE(shuffle)(a,b,0b00), "0b00");
    print(SSE(shuffle)(a,b,0b01), "0b01");
    print(SSE(shuffle)(a,b,0b10), "0b10");
    print(SSE(shuffle)(a,b,0b11), "0b11");
    print(SSE(unpacklo)(a,b), "unpacklo");
    print(SSE(unpackhi)(a,b), "unpackhi");
}

#ifdef __AVX__

void test_swap1()
{
    vec4 a = AVX(setr)( 1, 2, 3, 4);
    vec4 b = AVX(setr)(-1,-2,-3,-4);
    print(a, "a");
    print(b, "b");
    
    vec4 c = AVX(permute)(a,0x05);
    vec4 d = AVX(permute)(b,0x05);
    print(c, "c");
    print(d, "d");
    
    c = AVX(permute2f128)(a,b,0x20);
    d = AVX(permute2f128)(a,b,0x31);
    print(c, "c");
    print(d, "d");
}


void test_swap2()
{
    vec4 a = AVX(setr)( 1, -1, 2, -2);
    vec4 b = AVX(setr)( 3, -3, 4, -4);
    print(a, "a");
    print(b, "b");
    
    vec4 c = AVX(permute2f128)(a,a,0x08);
    vec4 d = AVX(permute2f128)(a,b,0x21);
    vec4 e = AVX(permute2f128)(a,a,0x81);
    print(c, "c");
    print(d, "d");
    print(e, "e");
}


void test_swap4()
{
    //vec4 b = AVX(broadcast)((__m128d*)addr);
    vec4 b = AVX(setr)(1, 2, 3, 4);
    print(b, "source");

    print(AVX(shuffle)(b, b, 0xC), "shuffle 0xC");
    print(AVX(shuffle)(b, b, 0b1100), "shuffle 0b1100");
    print(AVX(unpacklo)(b, b), "unpacklo");
    print(AVX(unpackhi)(b, b), "unpackhi");
}

#endif

#ifdef __AVX__

#include "random.h"
extern Random RNG;


/**
 Calculate ~16 single precision random number, using SIMD AVX instructions
 The numbers are Gaussian distributed, have low precision, and can be NAN
 F. Nedelec 2/1/2017
 */
void fast_gaussian(float vf[], size_t size, const uint32_t vi[])
{
    vecf fac = AVS(set1)(0x1p-31);
    vecf two = AVS(set1)(-2.0);
    vecf one = AVS(set1)(1.0);

    for ( unsigned ii=0; ii<size; ii+=16 )
    {
        __m256i i = _mm256_load_si256((__m256i*)(vi+ii  ));
        __m256i j = _mm256_load_si256((__m256i*)(vi+ii+8));
        vecf a = AVS(mul)(fac, AVS(cvtepi32)(i));
        vecf b = AVS(mul)(fac, AVS(cvtepi32)(j));
        vecf n = AVS(add)(AVS(mul)(a,a), AVS(mul)(b,b));
        //print(a, "a");print(b, "b");
        //print(n, "n");
        vecf t = AVS(cmp)(n, one, 0x12);
        vecf s = AVS(div)(n, AVS(mul)(two, AVS(log)(n)));
        //print(t, 't');
        n = AVS(rsqrt)(s);
        a = AVS(mul)(n, a);
        b = AVS(mul)(n, b);
        //print(a, "a"); print(b, "b");
        AVS(store)(vf+ii  , a);
        AVS(store)(vf+ii+8, b);
    }
    _mm_empty();
}


// pack array by removing 'nan' values
float * remove_nans(float * s, float * e)
{
    if ( e <= s )
        return s;
    --e;
    while ( s < e )
    {
        // find the next `nan` going upward:
        while ( *s == *s )
        {
            ++s;
            if ( e <= s )
                return e + ( *e == *e );
        }
        // going downward, skip `nan` values:
        while ( *e != *e )
        {
            --e;
            if ( e <= s )
                return e;
        }
        // flip the two values:
        *s = *e;
        ++s;
        --e;
    }
    return e + ( *e == *e );
}


#endif

int main(int argc, char * argv[])
{
#ifdef __AVX__
    RNG.seedTimer();
    float vf[SFMT_N32] = { 0 };
    for ( int i = 0; i < 128; ++ i )
    {
        fast_gaussian(vf, SFMT_N32, RNG.buffer());
        float * end = remove_nans(vf, vf+SFMT_N32);
        for ( float * f = vf; f < end; ++f )
            printf("%10.6f\n", *f);
        RNG.refill();
    }
#endif
    if ( 1 )
    {
        test_swap0();
        //test_swap1();
        //test_swap2();
        //test_swap4();
    }
    if ( 0 )
    {
        init();
        run(scalar,  "scalar ");
        run(vector2, "vector2");
        run(vector4, "vector4");
        run(vectorU, "vectorU");
    }
    
    return 0;
}


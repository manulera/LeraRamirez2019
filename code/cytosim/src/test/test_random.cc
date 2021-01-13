// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "random.h"
#include <cstring>
#include "tictoc.h"

extern Random RNG;




void printBits(FILE * f, const void * v, const int size)
{
    for ( int ii=0; ii < size; ++ii )
    {
        char c = ((char*)v)[ii];
        for ( int jj=7; jj >=0; --jj )
            fprintf(f, "%d", ( c >> jj ) & 1 );
        fprintf(f, ".");
    }
    fprintf(f, "\n");
}



void speed_test()
{
    const unsigned int cnt = 1 << 30;
    TicToc::tic();
    uint32_t u = 10;
    for (uint32_t j=0; j<cnt; ++j)
    {
        u = RNG.pint(1024);
        RNG.pint(u);
    }
    TicToc::toc("int");
}


void test_int()
{
    int j;
    
    for (j=0; j<21; j++)
        printf(" %10u%s", RNG.pint(), (j%7)==6 ? "\n" : "");
    
    printf("\n");
    
    for (j=0; j<90; j++)
        printf(" %2u%s", RNG.pint(99), (j%30)==29 ? "\n" : "");
    
    printf("\n");
    
    for (j=0; j<90; j++)
        printf(" %2u%s", RNG.pint_slow(99), (j%30)==29 ? "\n" : "");
    
    printf("\n");
    
    for (j=0; j<100; j++)
        printf(" %3i%s", RNG.sint_inc(99), (j%20)==19 ? "\n" : "");
    
    printf("\n");
    
    for (j=0; j<42; j++)
        printf(" %10.7f%s", RNG.sreal(), (j%7)==6 ? "\n" : "");
    
    printf("\n");
    
    for (j=0; j<42; j++)
        printf(" %8f%s", RNG.preal(), (j%7)==6 ? "\n" : "");
    
    printf("\n");
}


void silly_test()
{
    const uint32_t up = 1 << 30;
    
    const uint32_t cnt = 1 << 24;
    uint32_t hit = 0;
    
    for (uint32_t j=0; j<cnt; ++j)
        hit += ( RNG.pint() < up );

    printf(" prob( pint() < 1^30 ) = %f\n", hit/(float)cnt);
}


float convertFix(uint32_t x)
{
    //This assumes IEEE Standard 754 Floating point numbers
    //32 bits: 1 for sign, 8 for exponents, 23 for fraction
    const uint32_t FRAC     = 0x7FFFFFU;
    const uint32_t EXPON    = 127 << 23;
    uint32_t result = EXPON | ( x & FRAC );
    return *(float*)&result - 1.0;
}



void testbits()
{
    const int SCALE=2;
    float x;
    for ( int ii=0; ii <= SCALE; ++ii )
    {
        x = ii / float(SCALE);
        printf(" %f :", x);
        printBits(stdout, &x, 4);
        // x = -ii / float(SCALE);
        // printf("%f :", x);
        // printBits(stdout, &x, 4);
    }
    
    double y;
    for ( int ii=0; ii <= 20; ++ii )
    {
        y = convertFix( RNG.pint() );
        printf(" %f :", y);
        printBits(stdout, &y,8);
    }
}


#define TEST test
void test_test( const real prob, const int MAX )
{
    int cnt = 0, a, b, c;
    for ( int jj=0; jj < MAX; ++jj )
    {
        a = RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob);
        b = RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob);
        c = RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob);
        cnt += a + b + c;
    }
    printf("prob = %f measured = %f cnt = %i\n", prob, cnt / double(12*MAX), cnt);
}

void test_RNG(const int MAX)
{
    for ( int jj=0; jj < MAX; ++jj )
    {
        RNG.preal();RNG.preal();RNG.preal();RNG.preal();RNG.preal();
        RNG.preal();RNG.preal();RNG.preal();RNG.preal();RNG.preal();
    }
}


void test_float()
{
    printf("pfloat:     ");
    float x;
    for ( int kk=0; kk < 10; ++kk )
    {
        x = RNG.pfloat();
        printf(" %+f", x);
    }
    printf("\n");
    printf("sfloat:     ");
    for ( int kk=0; kk < 10; ++kk )
    {
        x = RNG.sfloat();
        printf(" %+f", x);
    }
    printf("\n");
    
    double d;
    printf("pdouble:    ");
    for ( int kk=0; kk < 10; ++kk )
    {
        d = RNG.pdouble();
        printf(" %+f", d);
    }
    printf("\n");
    printf("sdouble:    ");
    for ( int kk=0; kk < 10; ++kk )
    {
        d = RNG.sdouble();
        printf(" %+f", d);
    }
    printf("\n");
}

//==========================================================================

void test_uniform()
{
    unsigned cnt = 1<<28;
    real avg = 0;
    real var = 0;
    for ( unsigned i = 0; i < cnt; ++i )
    {
        real x = RNG.sreal();
        real y = RNG.sreal();
        real z = RNG.sreal();
        real t = RNG.sreal();
        avg += x + y + z + t;
        var += x*x + y*y + z*z + t*t;
    }
    cnt *= 4;
    avg /= cnt;
    var = var/cnt - avg * avg;
    printf("UNIFORM      avg = %.12e   var = %.12e\n", avg, var);
}


void test_gauss()
{
    printf("Gauss\n");
    const unsigned n_max = 1<<6;
    real vec[n_max] = { 0 };
    for ( int i = 0; i < 10000000; ++i )
    {
        unsigned n = RNG.pint(n_max);
        RNG.gauss_set(vec, n);
    }
}


void test_exponential()
{
    unsigned cnt = 1<<29;
    real avg = 0;
    real var = 0;
    for ( unsigned i = 0; i < cnt; ++i )
    {
        real x = RNG.exponential();
        real y = RNG.exponential();
        real z = RNG.exponential();
        real t = RNG.exponential();
        avg += x + y + z + t;
        var += x*x + y*y + z*z + t*t;
    }
    cnt *= 4;
    avg /= cnt;
    var = var/cnt - avg * avg;
    printf("EXPONENTIAL  avg = %.12e   var = %.12e\n", avg, var);
}


void test_poisson(unsigned int sup)
{
    for ( unsigned n = 0; n < sup; ++n )
    {
        int x = (int)(RNG.gauss() * sqrt(n) + n);
        printf("%10i %9i %9i %9i\n", n, RNG.poisson_knuth(n), RNG.poisson(n), x);
    }
}


//==========================================================================
//test 3 methods to generate a random event time, when the rate varies in time
// F. Nedelec, Oct 2005

//this is our standard method: 64s CPU
int method1(const int maxTime, const real rate[])
{
    for ( int ii=0; ii<maxTime; ++ii )
    {
        if (RNG.test(rate[ii])) return ii;
    }
    return maxTime;
}

//this is 'exact' and very slow: 370s CPU (an exponential at each step!)
int method2(const int maxTime, const real rate[])
{
    for ( int ii=0; ii<maxTime; ++ii )
    {
        if ( RNG.preal() < (1.-exp(-rate[ii])) )
            return ii;
    }
    return maxTime;
}

//this is exact, and the fastest method: 10s CPU!
int method3(const int maxTime, const real rate[])
{
    real T = -log( RNG.preal() );
    for ( int ii=0; ii<maxTime; ++ii )
    {
        T -= rate[ii];
        if ( T < 0 ) return ii;
    }
    return maxTime;
}


int testGillespie(const int method)
{
    //test new idea for gillespie with changing rate (Oct 2005)
    const int maxTime = 200;
    real rate[maxTime];
    for ( int ii=0; ii<maxTime; ++ii )
        rate[ii] = ( ii % 10 ) / 30.0;
    
    int bins[3][maxTime+1];
    for ( int ii=0; ii<=maxTime; ++ii )
    {
        bins[0][ii] = 0;
        bins[1][ii] = 0;
        bins[2][ii] = 0;
    }
    
    const int nbSamples = 1000000;
    const int subSamples = 10;
    int result;
    switch( method )
    {
        case 0:
            for ( int ii=0; ii<nbSamples; ++ii )
            {
                bins[0][ method1(maxTime, rate) ]++;
                bins[1][ method2(maxTime, rate) ]++;
                bins[2][ method3(maxTime, rate) ]++;
            }
            break;
            
        case 1:
            printf("method 1:");
            for ( int ii=0; ii<nbSamples; ++ii )
                for ( int jj=0; jj<subSamples; ++jj )
                    result = method1(maxTime, rate);
            return result;
            
        case 2:
            printf("method 2:");
            for ( int ii=0; ii<nbSamples; ++ii )
                for ( int jj=0; jj<subSamples; ++jj )
                    result = method2(maxTime, rate);
            return result;
            
        case 3:
            printf("method 3:");
            for ( int ii=0; ii<nbSamples; ++ii )
                for ( int jj=0; jj<subSamples; ++jj )
                    result = method3(maxTime, rate);
            return result;
    }
    
    
    FILE* file = fopen("test.out", "w");
    for ( int ii=0; ii<=maxTime; ++ii )
        fprintf(file, "%4i   %6i %6i %6i\n", ii, bins[0][ii], bins[1][ii], bins[2][ii]);
    fclose(file);
    return 0;
}


//==========================================================================
int main(int argc, char* argv[])
{
    if ( argc > 1 )
        RNG.seed(strtol(argv[1], 0, 10));
    else
        RNG.seedTimer();
    
    real rate = 0;
    if ( argc > 1 )
        rate = strtod(argv[1], 0);

    switch ( 6 )
    {
        case 0:
            test_poisson(1024);
            break;
            
        case 1:
            test_exponential();
            test_uniform();
            test_gauss();
            break;
    
        case 2:
            testGillespie(rate);
            break;

        case 3:
            for ( int kk=0; kk < 11; ++kk )
                test_test(rate*kk, 5000000);
            break;
            
        case 4:
            printf("sizeof(uint32_t) = %lu\n", sizeof(uint32_t));
            test_int();
            test_float();
            break;
            
        case 5:
            speed_test();
            break;
            
        case 6:
            silly_test();
            break;
    }
    
    printf("done\n");
    return EXIT_SUCCESS;
}


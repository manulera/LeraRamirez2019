// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/*
 A test for the Floating Point Exceptions (Signal)
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <csignal>
#include <cmath>

/*
 icpc --help
 
 -fp-trap=<arg>[,<arg>,...]
 control floating point traps at program start.  <arg> can be of the
 following values
 [no]divzero   - [Do not] trap on division by zero
 [no]inexact   - [Do not] trap on inexact result
 [no]invalid   - [Do not] trap on invalid operation
 [no]overflow  - [Do not] trap on overflow
 [no]underflow - [Do not] trap on underflow
 [no]denormal  - [Do not] trap on denormal
 all           - enable trap on all of the above
 none          - trap on none of the above
 common        - trap on most commonly used IEEE traps
 (invalid, division by zero, overflow)
 -fp-trap-all=<arg>[,<arg>,...]
 control floating point traps in every routine.  <arg> can be of the
 values specified in -fp-trap
 */

typedef double real;

void fpe_handler(int sig)
{
    std::cout << std::endl;
    psignal(sig, "Cytosim");
    exit(sig);
}

void modulo()
{
    std::cout << "   x    fmod remainder";
    for ( real x = -4; x <= 4; x += 0.5 )
    {
        std::cout << "\n" << std::setw(5) << x;
        std::cout << "  " << std::setw(5) << fmod(x, 2.0);
        std::cout << "  " << std::setw(5) << remainder(x, 2.0);
    }
    std::cout << std::endl;
}

void infinities()
{
    std::cout << "0   < inf = " << ( 0 < INFINITY ) << std::endl;
    std::cout << "inf < inf = " << ( INFINITY < INFINITY ) << std::endl;
    real z = 0;
    real y = 0.0 / z;
    real x = 1.0 / z;
    std::cout << " 1.0/0.0 = " << x << std::endl;
    std::cout << " 0.0/0.0 = " << y << std::endl;
}

void print_numbers()
{
    std::cout << " 1.0 / 0 = " <<  1.0 / 0 << std::endl;
    std::cout << "-1.0 / 0 = " << -1.0 / 0 << std::endl;
    std::cout << " 0.0 / 0 = " <<  0.0 / 0 << std::endl;
    std::cout << "-log(0)  = " << -log(0.0) << std::endl;
#if ( 0 )
    std::cout << "absf(-2) = " << sMath::absf(-2.0) << std::endl;
    std::cout << "absf(-1) = " << sMath::absf(-1.) << std::endl;
    std::cout << "absf(+1) = " << sMath::absf(+1.) << std::endl;
    std::cout << "absf(+2) = " << sMath::absf(+2.) << std::endl;
#endif
}


void read_numbers(std::string const& str)
{
    std::stringstream is(str);
    double x, y, z;
    if ( is >> x )
        std::cout << "read x = " << x << std::endl;
    if ( is >> y )
        std::cout << "read y = " << y << std::endl;
    if ( is >> z )
        std::cout << "read z = " << z << std::endl;
    
    char tmp[128] = { 0 };
    is.clear();
    is.readsome(tmp, sizeof(tmp));
    std::cout << "remaining >" << tmp << std::endl;
}


int main()
{
    //read_numbers("1 2 nnnn .0 aaaaa");
    
    if ( signal(SIGFPE, fpe_handler) == SIG_ERR )
    {
        std::cout << "Could not register SIGFPE handler\n";
        return EXIT_FAILURE;
    }
    infinities();
    print_numbers();
    std::cout << "test completed" << std::endl;
    return 0;
}

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 @file real.h
 SYNOPSIS: we define and use "real" to be able to easily change
 the floating point precision, depending on the application.
 REAL_EPSILON is a lower limit of the precision achieved.
 
 It is possible to select double or single precision by editing this file:
 - To use float, define REAL_IS_FLOAT
 - if REAL_IS_FLOAT is not defined, double precision will be used
 .
 
 Cytosim is faster in single precision, but the iterative solver used
 in Meca::solve() (conjugate-gradient) may fail in adverse conditions.
 
 It is safer, and STRONGLY ADVISED therefore, to use double precision,
 and to not define 'REAL_IS_FLOAT'
 */

#ifndef REAL
#define REAL


#include <cmath>
#include <cfloat>


//#define REAL_IS_FLOAT


#ifdef REAL_IS_FLOAT
    /// real is an alias to float
    typedef float real;
    const real REAL_EPSILON = 128 * FLT_EPSILON;
#else
    /// real is an alias to double
    typedef double real;
    const real REAL_EPSILON = 128 * DBL_EPSILON;
#endif


/// square of the argument: `x * x`
inline real sqr(const real x) { return x * x; }

/// cube of the argument: `x * x * x`
inline real cub(const real x) { return x * x * x; }

/// sign of a float
inline real sign(const float x) { return copysignf(1.0f, x); }

/// sign of a float
inline real sign(const double x) { return copysign(1.0, x); }

/// absolute value of `x`
//inline real fabs(const real x) { return std::fabs(x); }

#if ( 0 )

#include <iostream>
#include <sstream>
#include "tokenizer.h"
/**
 Redefining the standard extraction operator is normally not necessary,
 but this can be useful to track certain bugs
 */
inline std::istringstream& operator >> (std::istringstream& iss, real& x)
{
    std::string str = Tokenizer::get_real(iss);
    if ( str.empty() )
        return iss;
    try {
        x = std::stod(str);
        //std::clog << " using custom >> |" << str << "| -> " << x << std::endl;
    }
    catch ( std::invalid_argument & e ) {
        std::cerr << " error in custom >> |" << str << "| " << std::endl;
    }
    return iss;
}
#endif

#endif

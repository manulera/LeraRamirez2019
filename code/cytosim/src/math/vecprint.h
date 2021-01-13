// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef VECPRINT_H
#define VECPRINT_H

#include <iostream>
#include <iomanip>
#include <cmath>


/// Templated functions to print Vectors and Matrices with minimal formatting
namespace VecPrint
{
    /// print a vector on a line
    template< typename T >
    void vecPrint(std::ostream & os, unsigned m, const T* X, int digits = 3)
    {
        if ( X == 0 ) {
            os << " void\n";
            return;
        }
        
        char str[32], fmt[32];
        snprintf(fmt, sizeof(fmt), " %%%i.%if", digits+5, digits);
        for ( unsigned ii = 0; ii < m; ++ii )
        {
            snprintf(str, sizeof(str), fmt, X[ii]);
            os << str;
        }
        os << std::endl;
    }
    
    
    /// print a vector in column format
    template< typename T >
    void vecDump(std::ostream & os, unsigned m, const T* X, int digits = 8)
    {
        if ( X == 0 ) {
            os << " void\n";
            return;
        }
        
        char str[32], fmt[32];
        snprintf(fmt, sizeof(fmt), " %%%i.%ie", digits+6, digits);
        for ( unsigned ii = 0; ii < m; ++ii )
        {
            snprintf(str, sizeof(str), fmt, X[ii]);
            os << str << std::endl;
        }
    }
    
    
    /// print matrix `mat` of size m*n, and leading dimension `ldd` with precision 'digits'
    template< typename T >
    void matPrint(std::ostream & os, unsigned m, unsigned n, const T* mat, unsigned ldd, int digits = 3)
    {
        const real threshold = 0.01;
        if ( mat == 0 ) {
            os << " void\n";
            return;
        }
        
        char str[32], fmt[32], zer[32];
        snprintf(fmt, sizeof(fmt), " %%%i.%if", digits+6, digits);
        snprintf(zer, sizeof(zer), " %%%is", digits+6);
       
        for ( unsigned ii = 0; ii < m; ++ii )
        {
            for ( unsigned jj = 0; jj < n; ++jj )
            {
                T val = mat[ii+ldd*jj];
                if ( std::abs(val) < threshold )
                    snprintf(str, sizeof(str), zer, ".");
                else
                    snprintf(str, sizeof(str), fmt, mat[ii+ldd*jj]);
                os << str;
            }
            os << std::endl;
        }
    }
    
    
    /// print a matrix in sparse format: line_index, column_index, value 
    template< typename T >
    void matSparsePrint(std::ostream & os, unsigned m, unsigned n, const T* mat, unsigned ldd, int digits = 8)
    {
        if ( mat == 0 ) {
            os << " void\n";
            return;
        }
        
        char str[32], fmt[32];
        snprintf(fmt, sizeof(fmt), " %%%i.%if\n", digits+6, digits);
        for (unsigned ii = 0; ii < m; ++ii )
            for (unsigned jj = 0; jj < n; ++jj )
            {
                snprintf(str, sizeof(str), fmt, mat[ii+ldd*jj]);
                os << ii << " " << jj << " " << str;
            }
        os << std::endl;
    }
    
    
    /// print a matrix in sparse format, but adding `off` to all line and column indices
    template< typename T >
    void matSparsePrintOffset(std::ostream & os, unsigned m, unsigned n, const T* mat, unsigned ldd, int off, int digits = 8)
    {
        if ( mat == 0 ) {
            os << " void\n";
            return;
        }
        
        
        char str[32], fmt[32];
        snprintf(fmt, sizeof(fmt), " %%%i.%if\n", digits+6, digits);
        for (unsigned ii = 0; ii < m; ++ii )
            for (unsigned jj = 0; jj < n; ++jj )
            {
                snprintf(str, sizeof(str), fmt, mat[ii+ldd*jj]);
                os << ii+off << " " << jj+off << str;
            }
        os << std::endl;
    }
}

#endif

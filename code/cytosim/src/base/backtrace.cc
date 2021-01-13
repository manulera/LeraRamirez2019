// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 24/04/2010.


#include "backtrace.h"

// You can disable backtrace by changing the line below into '#if ( 0 )'
#ifdef __GNUC__


#include <execinfo.h>
#include <stdlib.h>

/**
 * print the current call-stack using functions from the GNU C Library:
 * - backtrace()
 * - backtrace_symbols()
 * .
 * provided by <execinfo.h>
 */
void print_backtrace(FILE * out)
{
    void* callstack[128];
    int size = backtrace(callstack, 128);
#if ( 1 )
    char** strs = backtrace_symbols(callstack, size);

    fprintf(out, "Cytosim execution stack:\n");
    for ( size_t ii = 1; ii < size; ++ii )
        fprintf(out, "      %s\n", strs[ii]);
    
    free(strs);
#else
    backtrace_symbols_fd(callstack, size, out);
#endif
}

#else

void print_backtrace(FILE * out)
{
    fprintf(out, "Execution stack information unavailable\n");
}

#endif


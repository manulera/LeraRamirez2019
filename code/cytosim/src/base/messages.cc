// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "messages.h"

//the global instantiation Messages used for input/output
Messages MSG;

#include "assert_macro.h"
#include <cstdlib>
#include <cctype>
#include <iostream>

/**
 default verbose level = 4
 */
Messages::Messages() : FileWrapper(stdout)
{
    verbose_ = 4;
    num_warnings_ = 0;
}


void Messages::operator()(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    vfprintf(mFile, fmt, args);
    va_end(args);
}



void Messages::operator()(int level, const char* fmt, ...)
{
    if ( verbose_ >= level )
    {
        va_list args;
        va_start(args, fmt);
        vfprintf(mFile, fmt, args);
        va_end(args);
    }
}


void Messages::warning(const char* fmt, ...)
{
    if ( verbose_ >= 0  &&  num_warnings_ < max_warnings )
    {
        fprintf(mFile, "warning: ");
        va_list args;
        va_start(args, fmt);
        vfprintf(mFile, fmt, args);
        va_end(args);
        
        if (++num_warnings_ >= max_warnings)
            fprintf(mFile, "warning messages are now silent\n");
    }
}


void Messages::warning(std::string const& str)
{
    if ( verbose_ >= 0  &&  num_warnings_ < max_warnings )
    {
        fprintf(mFile, "warning: %s", str.c_str());
        if (++num_warnings_ >= max_warnings)
            fprintf(mFile, "warning messages are now silent\n");
    }
}


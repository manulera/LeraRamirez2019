// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
//  Created by François Nédélec on 27/9/08.


#include <ctime>

/// A set of functions related to time
/**
 Functions to get wall-time, and processor-time,
 calling the C-standard library.
 */
namespace TicToc
{    
    
    /// set current date in short format, `buf` should be 26 character long or more
    void    get_date(char * buf, size_t buf_size);
    
    /// set current date in short format, `buf` should be 26 character long or more
    void    get_date(char * buf, size_t buf_size, bool no_year);

    
    /// apprimately the number of days since Jan 1 2000
    int     days_since_2000();
    
    /// year
    int     year();

    /// day of the year (0-365)
    int     day_of_the_year();
    
    /// hour of the day (0-23)
    int     hours_today();

    /// number of seconds since midnight
    long    seconds_today();
    
    /// number of milli-seconds since midnight
    long    milli_seconds_today();
 
    
    /// CPU time in short format, `buf` must be   
    clock_t processor_time(char * buf, size_t buf_size, clock_t, double& sum);
    
    /// call to start timer
    void    tic();
    
    /// return number of milli-seconds elapsed since last call to `tic()`
    double  toc();

    /// call to report time elapsed since 'tic'. Will print a message if msg!=0
    double  toc(const char * msg, const char * nl = 0);
    
};


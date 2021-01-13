// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by François Nédélec on 27/9/08.


#include "tictoc.h"
#include <cstdio>
#include <ctime>
#include <cstring>
#include <sys/time.h>

#pragma mark Wall time

/**
 This calls the C-library functions time() and ctime_r()
 */
void TicToc::get_date(char * buf, size_t buf_size)
{
    if ( buf_size > 25 )
    {
        time_t now = time(0);
        //asctime_r(localtime(&now), buf);
#if ( 0 )
        ctime_r(&now, buf);
#else
        strncpy(buf, ctime(&now), buf_size);
#endif
        // remove new line:
        buf[24]='\0';
    }
    // terminate string:
    buf[buf_size-1] = '\0';
}


void TicToc::get_date(char * buf, size_t buf_size, bool no_year)
{
    if ( buf_size > 25 )
    {
        TicToc::get_date(buf, buf_size);
        // remove year:
        if ( no_year )
            buf[19]='\0';
    }
}


int  TicToc::days_since_2000()
{
    time_t now = time(0);
    tm * loc = localtime(&now);
    return loc->tm_year - 2000 + loc->tm_yday;
}


int  TicToc::year()
{
    time_t now = time(0);
    tm * loc = localtime(&now);
    return loc->tm_year;
}


int  TicToc::day_of_the_year()
{
    time_t now = time(0);
    tm * loc = localtime(&now);
    return loc->tm_yday;
}


int  TicToc::hours_today()
{
    time_t now = time(0);
    tm * loc = localtime(&now);
    return loc->tm_hour;
}


long TicToc::seconds_today()
{
    time_t now = time(0);
    tm * loc = localtime(&now);
    return loc->tm_sec + 60 * ( loc->tm_min + 60 * loc->tm_hour );
}


long TicToc::milli_seconds_today()
{
    long msec = seconds_today() * 1e3;
    struct timeval now;
    gettimeofday(&now, 0);
    return msec + now.tv_usec * 1e-3;
}


#pragma mark CPU time


/**
 This calls the C-library function clock()
 */
clock_t TicToc::processor_time(char * buf, size_t buf_size, clock_t old, double& sum)
{
    clock_t now = clock();
    double sec = double( now - old ) / CLOCKS_PER_SEC;
    sum += sec;
    snprintf(buf, buf_size, "%10.0fms  %10.0fs", 1000 * sec, sum);
    return now;
}


#if ( 1 )

struct timeval tic_t;

void TicToc::tic()
{
    gettimeofday(&tic_t, 0);
}

double TicToc::toc()
{
    timeval toc_t;
    gettimeofday(&toc_t, 0);
    return (toc_t.tv_sec-tic_t.tv_sec)*1e3 + (toc_t.tv_usec-tic_t.tv_usec)*1e-3;
}

#else

#include <ctime>
clock_t tic_t;

void TicToc::tic()
{
    tic_t = clock();
}

double TicToc::toc()
{
    return (1000 * ( clock() - tic_t ))/CLOCKS_PER_SEC;
}

#endif


double TicToc::toc(const char * msg, const char * nl)
{
    double elapsed = toc();
    if ( nl )
        printf("%s : %4.0f ms %s", msg, elapsed, nl);
    else
        printf("%s : %4.0f ms ", msg, elapsed);
    return elapsed;
}



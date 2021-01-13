// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 This is a program to analyse simulation results:
 it reads a trajectory-file, and print some data from it.
*/

#include <fstream>
#include <sstream>

#include "stream_func.h"
#include "frame_reader.h"
#include "iowrapper.h"
#include "glossary.h"
#include "messages.h"
#include "splash.h"
#include "parser.h"
#include "simul.h"

Simul simul;
int verbose = 1;

void help(std::ostream& os)
{
    os << "Synopsis: generate reports/statistics about cytosim's objects\n";
    os << "\n";
    os << "Syntax:\n";
    os << "       report WHAT [OPTIONS]\n";
    os << "Supported options:\n";
    os << "       prefix = time]\n";
    os << "       precision = INTEGER\n";
    os << "       verbose = 0\n";
    os << "       frame = POSITIVE_INTEGER\n";
    os << "       output = FILE_NAME\n";
    os << "\n";
    os << "This will generate reports using Simul::report()\n";
    os << "See the documentation of Simul::report() for a list of possible values for WHAT\n";
    os << "\n";
}

//------------------------------------------------------------------------------

void report_raw(std::ostream& os, std::string const& what, int frm, Glossary& opt)
{
    if ( verbose > 0 )
    {
        os << "\n% frame   " << frm;
        simul.report(os, what, opt);
    }
    else
    {
        std::stringstream ss;
        simul.report(ss, what, opt);
        StreamFunc::skip_lines(os, ss, '%');
    }
}


void report_prefix(std::ostream& os, std::string const& what, int frm, Glossary& opt)
{
    char prefix[256] = { 0 };
    snprintf(prefix, sizeof(prefix), "%9.3f ", simul.time());
    
    std::stringstream ss;
 
    if ( verbose )
    {
        os << "% frame   " << frm << '\n';
        simul.report(ss, what, opt);
        StreamFunc::prefix_lines(os, ss, prefix, '%', 0);
    }
    else
    {
        simul.report(ss, what, opt);
        StreamFunc::prefix_lines(os, ss, prefix, 0, '%');
    }
}


void report(std::ostream& os, std::string const& what, int frm, Glossary& opt)
{
    std::string prefix;
    opt.set(prefix, "prefix");

    try
    {
        if ( prefix == "time" )
            report_prefix(os, what, frm, opt);
        else
            report_raw(os, what, frm, opt);
    }
    catch( Exception & e )
    {
        std::cerr << "Aborted: " << e.what() << '\n';
        exit(EXIT_FAILURE);
    }
}


//------------------------------------------------------------------------------


int main(int argc, char* argv[])
{
    if ( argc < 2 || strstr(argv[1], "help") )
    {
        help(std::cout);
        return EXIT_SUCCESS;
    }
    
    if ( strstr(argv[1], "info") || strstr(argv[1], "--version")  )
    {
        splash(std::cout);
        print_version(std::cout);
        std::cout << " DIM = " << DIM << '\n';
        return EXIT_SUCCESS;
    }

    std::string input = simul.prop->trajectory_file;
    std::string str, what = argv[1];
    
    std::ostream * osp = &std::cout;
    std::ofstream ofs;
    
    Glossary opt;
    opt.read_strings(argc-2, argv+2);
    
    unsigned frame = 0;
    unsigned period = 1;

    opt.set(input,   ".cmo") || opt.set(input, "input");;
    opt.set(verbose, "verbose");
    opt.set(period, "period");

    FrameReader reader;
    RNG.seedTimer();

    try
    {
        simul.loadProperties();
        reader.openFile(input);
    }
    catch( Exception & e )
    {
        std::clog << "Aborted: " << e.what() << '\n';
        return EXIT_FAILURE;
    }

    if ( opt.set(str, "output") )
    {
        try {
            ofs.open(str.c_str());
        }
        catch( ... )
        {
            std::clog << "Cannot open output file\n";
            return EXIT_FAILURE;
        }
        osp = &ofs;
    }
    
    MSG.silent();
    
    
    if ( opt.has_key("frame") )
    {
        // multiple frame indices can be specified:
        unsigned s = 0;
        while ( opt.set(frame, "frame", s) )
        {
            // try to load the specified frame:
            if ( 0 == reader.loadFrame(simul, frame) )
                report(*osp, what, frame, opt);
            else
            {
                std::cerr << "Error: missing frame " << frame << '\n';
                return EXIT_FAILURE;
            }
            ++s;
        };
    }
    else
    {
        // process every 'period' frame in the file:
        while ( 0 == reader.loadNextFrame(simul) )
        {
            if ( 0 == frame % period )
                report(*osp, what, frame, opt);
            ++frame;
        }
    }
    
    /// check that all specified parameters have been used:
    std::stringstream ss;
    if ( opt.warnings(ss) > 1 )
        std::cerr << ss.str() << '\n';
    
    return EXIT_SUCCESS;
}

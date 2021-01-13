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
    os << "       reportF WHAT [verbose=0]\n";
    os << "\n";
    os << "This will generate the same reports as Simul::report()\n";
    os << "See the documentation of Simul::report() for a list of possible values for WHAT\n";
    os << "\n";
    os << "The report of each frame of the trajectory file is send to a different file\n";
}

//------------------------------------------------------------------------------

void report(std::ostream& os, std::string const& what, Glossary& opt)
{
    if ( verbose )
    {
        simul.report(os, what, opt);
    }
    else
    {
        std::stringstream ss;
        simul.report(ss, what, opt);
        StreamFunc::skip_lines(os, ss, '%');
    }
}

//------------------------------------------------------------------------------


int main(int argc, char* argv[])
{
    MSG.silent();
    
    if ( argc < 2 || strstr(argv[1], "help") )
    {
        help(std::cout);
        return EXIT_SUCCESS;
    }
    
    if ( strstr(argv[1], "info") || strstr(argv[1], "--version")  )
    {
        splash(std::cout);
        print_version(std::cout);
        std::cout << "   DIM = " << DIM << std::endl;
        return EXIT_SUCCESS;
    }

    std::string input = simul.prop->trajectory_file;
    std::string str, what = argv[1];
    
    Glossary opt;
    opt.read_strings(argc-2, argv+2);
    opt.set(input, ".cmo") || opt.set(input, "input");;
    opt.set(verbose, "verbose");
    FrameReader reader;
    RNG.seedTimer();

    try
    {
        simul.loadProperties();
        reader.openFile(input);
        
        unsigned frame = 0;
        char filename[256];
        
        // load all frames in the file:
        while ( 0 == reader.loadNextFrame(simul) )
        {
            snprintf(filename, sizeof(filename), "report%04i.txt", frame);
            std::ofstream out(filename);
            report(out, what, opt);
            ++frame;
        }
    }
    catch( Exception & e )
    {
        std::cerr << "Aborted: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}

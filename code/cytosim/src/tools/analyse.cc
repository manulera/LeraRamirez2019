// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 This is a program to analyse simulation results:
 it reads a trajectory-file, and print some data from it.
*/

#include "splash.h"
#include "frame_reader.h"
#include "iowrapper.h"
#include "glossary.h"
#include "messages.h"
#include "parser.h"
#include "simul.h"

Simul simul;



void help(std::ostream& os)
{
    os << "Synopsis: generate reports/statistics about cytosim's objects\n";
    os << "\n";
    os << "Syntax:\n";
    os << "       analyse WHAT [frame=INTEGER]\n";
    os << "\n";
    os << "Analyse will generate the same reports as Simul::report()\n";
    os << "See the HTML documentation of Simul::report() for a list of possible values for WHAT\n";
    os << "\n";
    os << "The report is send to the standard output";
    os << "\n";
}

//------------------------------------------------------------------------------


int main(int argc, char* argv[])
{
    if ( argc < 2 || strstr(argv[1], "help") )
    {
        help(std::cout);
        return EXIT_SUCCESS;
    }
    
    if ( strstr(argv[1], "info") || strstr(argv[1], "version")  )
    {
        splash(std::cout);
        print_version(std::cout);
        std::cout << " DIM = " << DIM << std::endl;
        return EXIT_SUCCESS;
    }

    std::string what = argv[1];
    
    Glossary opt;
    opt.read_strings(argc-2, argv+2);
    FrameReader reader;
    
    try
    {
        simul.loadProperties();
        reader.openFile(simul.prop->trajectory_file);
    }
    catch( Exception & e )
    {
        std::clog << "Aborted: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    MSG.silent();
    unsigned frm = 0;
    
    // multiple frame indices can be specified:
    if ( opt.set(frm, "frame") )
    {
        unsigned s = 0;
        do {
            // try to load the specified frame:
            if ( 0 == reader.loadFrame(simul, frm) )
            {
                std::cout << "% frame   " << frm << '\n';
                simul.report(std::cout, what, opt);
            }
            else
            {
                std::cerr << "Error: missing frame " << frm << std::endl;
                return EXIT_FAILURE;
            }
            ++s;
        } while ( opt.set(frm, "frame", s) );
    }
    else
    {
        // process every frame in the file:
        while ( 0 == reader.loadNextFrame(simul) )
        {
            std::cout << "% frame   " << frm << '\n';
            simul.report(std::cout, what, opt);
        }
    }
    
    return EXIT_SUCCESS;
}

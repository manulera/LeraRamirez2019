// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "simul.h"
#include "parser.h"
#include "messages.h"
#include "glossary.h"
#include "exceptions.h"
#include "backtrace.h"
#include "splash.h"
#include "tictoc.h"
#include <csignal>


using std::endl;


void help(std::ostream & os = std::cout)
{
    os << " Command line options:" << endl;
    os << "    FILENAME   set config file if FILENAME ends by `.cym'" << endl;
    os << "    terminal   send messages to terminal instead of `messages.cmo'" << endl;
    os << "    info       print build options" << endl;
    os << "    help       print this message" << endl;
    os << "    -          do not splash standard output" << endl;
}

void handle_abort(int sig)
{
    fprintf(stderr, "Cytosim: abort\n");
    print_backtrace(stderr);
    exit(sig);
}


void handle_interrupt(int sig)
{
    fwrite("killed\n\n", 1, 8, MSG.file());
    exit(sig);
}

//------------------------------------------------------------------------------
//=================================  MAIN  =====================================
//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    // Catch interrupting signals:
    if ( signal(SIGINT, handle_interrupt) )
        std::cerr << "Could not register SIGINT handler\n";
    
    if ( signal(SIGTERM, handle_interrupt) )
        std::cerr << "Could not register SIGTERM handler\n";
    
    if ( signal(SIGSEGV, handle_abort) )
        std::cerr << "Could not register SIGSEGV handler\n";
    
    if ( signal(SIGABRT, handle_abort) )
        std::cerr << "Could not register SIGABRT handler\n";

    Simul simul;

    //parse the command line:
    Glossary glos;
    glos.read_strings(argc-1, argv+1);
    
    if ( glos.use_key("help") )
    {
        help();
        return EXIT_SUCCESS;
    }
    
    if ( glos.use_key("info") || glos.use_key("--version")  )
    {
        splash(std::cout);
        print_version(std::cout);
        std::cout << "   DIM = " << DIM << endl;
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
        std::cout << "   ANISOTROPIC_FIBER_DRAG" << endl;
#endif
        return EXIT_SUCCESS;
    }
    
    if ( ! glos.use_key("stdout") )
    {
        MSG.open("messages.cmo", "w");
    }        
    
    if ( ! glos.use_key("-") )
        splash(std::cout);

#ifdef CODE_VERSION
    MSG("CYTOSIM PI version %s\n", CODE_VERSION);
#else
    MSG("CYTOSIM PI\n");
#endif
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
    MSG("    with ANISOTROPIC_FIBER_DRAG\n");
#endif

    try {
        simul.initialize(glos);
    }
    catch( Exception & e ) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch(...) {
        std::cerr << "Error: an unknown exception occured during initialization" << std::endl;
        return EXIT_FAILURE;
    }
    
    glos.warnings(std::cerr);
    MSG.flush();
    
    try {
        if ( Parser(simul, 1, 1, 1, 1, 1).readConfig() )
            std::cerr << "You must specify a config file\n";
    }
    catch( Exception & e ) {
        std::cerr << "\nError: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
    catch(...) {
        std::cerr << "\nAn unknown exception occured\n";
        return EXIT_FAILURE;
    }
    
    char date[26];
    TicToc::get_date(date, sizeof(date));
    MSG("%s\n", date);
    MSG("end\n");
    MSG.close();
    return EXIT_SUCCESS;
}

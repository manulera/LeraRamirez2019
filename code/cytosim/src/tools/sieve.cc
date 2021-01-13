// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "simul.h"
#include "parser.h"
#include "glossary.h"
#include "iowrapper.h"
#include "exceptions.h"

void help()
{
    printf("Synopsis:\n");
    printf("   `sieve` is program to manipulate cytosim trajectory file.\n");
    printf("   It reads a trajectory files, and produces another trajectory file\n");
    printf("\n");
    printf("   The file is written in the latest format, either binary or text-based.\n");
    printf("   A category of objects can be removed with option skip=WHAT.\n");
    printf("   If the specified output already exists, data is appended to it.\n");
    printf("\n");
    printf("Usage:\n");
    printf("    sieve input_file output_file [options]\n\n");
    printf("Possible options:\n");
    printf("    binary=0   generate output in text format\n");
    printf("    binary=1   generate output in binary format\n");
    printf("    verbose=?  set the verbose level\n");
    printf("    skip=WHAT  remove all objects of class WHAT\n");
    printf("\n");
    printf("Example:\n");
    printf("    sieve objects.cmo objects.txt binary=0\n");
    printf("(DIM = %i)\n", DIM);
}


int main(int argc, char* argv[])
{
    if ( argc < 3 )
    {
        help();
        return EXIT_SUCCESS;
    }

    Simul simul;
    Glossary glos;
    
    std::string input  = argv[1];
    std::string output = argv[2];
    glos.read_strings(argc-3, argv+3);
    
    int verbose = 0;
    glos.set(verbose, "verbose");
    
    ObjectSet * skip_set = 0;
    std::string skip;
    if ( glos.set(skip, "skip") )
       skip_set = simul.findSet(skip);
    
    bool binary = true;
    glos.set(binary, "binary");

    
    Inputter in;
    in.vectorSize(DIM);
    try {
        simul.loadProperties();
        in.open(input.c_str(), "rb");
    }
    catch( Exception & e ) {
        std::cerr << "Error opening input file `" << input << "' :" << std::endl;
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }
    
    //std::clog << ">>>>>> Copying `" << input << "' -> `" << output << "'" << std::endl;

    int frm = 0;

    while ( in.good() )
    {
        try {
            if ( simul.reloadObjects(in) )
                return EXIT_SUCCESS;
        }
        catch( Exception & e ) {
            std::clog << "Error in frame " << frm << ":\n";
            std::clog << "    " << e.what() << std::endl;
        }

        ++frm;

        if ( skip_set )
            skip_set->erase();
            
        //simul.reportInventory(std::cout);
        //std::clog << "\b\b\b\b\b" << std::setw(5) << cnt;
            
        try {
            simul.writeObjects(output.c_str(), binary, true);
        }
        catch( Exception & e ) {
            std::clog << "Error writing `" << output << "' :\n";
            std::clog << "    " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
    }
    return 0;
}

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 This is a program to export simulation objects, 
 in a format that is easy to read in Blender
 F. Nedelec, 20 Sept. 2017
*/

#include <fstream>
#include "frame_reader.h"
#include "iowrapper.h"
#include "glossary.h"
#include "parser.h"
#include "simul.h"
#include "fiber_disp.h"

/*
 \todo
 
 Store the points of the filaments as they are generated.
 For the Couple, relocate Hands to the closest filament-point.
 
 */

bool binary = 0;
FILE * file = 0;
int style = 2;


void help(std::ostream& os)
{
    os << "Synopsis: generate simple binary files with cytosim's objects\n";
    os << "\n";
    os << "Syntax:\n";
    os << "       cymart [INPUTFILE] [binary=0]\n";
    os << "\n";
    os << "Each frame of the trajectory file is sent to a separate file\n";
    std::cout << "   DIM = " << DIM << '\n';
}

//------------------------------------------------------------------------------


void savePoint(int dat[3], Vector2 const& pos)
{
    float vec[3] = { 0 };
    pos.put(vec);
    
    if ( binary )
    {
        fwrite(dat, 3, sizeof(int), file);
        fwrite(vec, 3, sizeof(float), file);
    }
    else
    {
        fprintf(file, "%i %i %i ", dat[0], dat[1], dat[2]);
        fprintf(file, "%.6f %.6f %.6f\n", vec[0], vec[1], vec[2]);
    }
}


void savePoint(int dat[3], Vector3 const& pos)
{
    float vec[3] = { 0 };
    pos.put(vec);
    
    if ( binary )
    {
        fwrite(dat, 3, sizeof(int), file);
        fwrite(vec, 3, sizeof(float), file);
    }
    else
    {
        fprintf(file, "%i %i %i ", dat[0], dat[1], dat[2]);
        fprintf(file, "%.6f %.6f %.6f\n", vec[0], vec[1], vec[2]);
    }
}


/// save filament backbone
void saveBackbone(Fiber const& fib)
{
    int dat[3] = { 2, fib.identity(), fib.prop->index() };
    
    for ( unsigned p = 0; p < fib.nbPoints(); ++p )
        savePoint(dat, fib.posP(p));
}


/// save double helical filament
void saveActin(Fiber const& fib)
{
    int dat[3] = { 1, fib.identity(), 0 };

    // axial translation between two sucessive monomers:
    const real dab = 0.00275;
    // distance from central axis to center of monomers
    real off = 0.0045 - dab;
    
    // rotation angle between consecutive monomers
    const real dan = 166 * M_PI / 180;
    const real cs = cos(dan);
    const real sn = sin(dan);
    
    real ab = 0;
    Vector3 p; //position of monomer;
    Vector3 d(fib.dirEndM(), DIM); // unit tangent to centerline
    Vector3 n(fib.normal, DIM);    // normal direction
    if ( n.normSqr() < 0.5 )
        n = d.orthogonal(1.0);
    n = d.orthogonal(n, 1.0);
    fib.normal.get(n, DIM);
    //std::clog << fib.reference() << " " << n << "    " << n.normSqr() << '\n';
    
    int cnt = 0;
    // rotate normal until we reach the MINUS_END
    while ( ab < fib.abscissaM() )
    {
        ++cnt;
        n = d.rotate(n, cs, sn);
        ab += dab;
    }
    
    // draw the monomers until the PLUS_END
    while ( ab < fib.abscissaP() )
    {
        d = Vector3(fib.dir(ab), DIM);
        n = d.rotate(n, cs, sn);
        p = Vector3(fib.pos(ab), DIM) + off * n;

        // use different tones to individualize the two strands:
        dat[2] = 1 + ( cnt & 1 );
        
        // change color for the last 2 monomers:
        if ( ab + dab > fib.abscissaP() )
            dat[2] += 2;

        savePoint(dat, p);
        ab += dab;
        ++cnt;
    }
}



/// save tubular structure made of 13-protofilaments
void saveMicrotubule(Fiber const& fib)
{
    real da[] = {0,0.000923,0.001846,0.002769,0.003692,0.004615,0.005538,0.006461,0.007384,0.008308,0.009231,0.010154,0.011077};
    real dx[] = {0.8855,0.5681,0.1205,-0.3546,-0.7485,-0.9709,-0.9709,-0.7485,-0.3546,0.1205,0.5681,0.8855,1.0000};
    real dy[] = {-0.4647,-0.8230,-0.9927,-0.9350,-0.6631,-0.2393,0.2393,0.6631,0.9350,0.9927,0.8230,0.4647,0};

    int dat[3] = { 1, fib.identity(), 0 };

    // axial translation between two sucessive monomers:
    const real dab = 0.004;
    // enlarged radius of monomers makes them overlap slighlty
    const real rad = 0.7 * dab;
    // distance from central axis to center of monomers
    real off = 0.025 / 2 - rad;
    
    real ab = dab * ceil( fib.abscissaM() / dab );
    Vector3 d(fib.dir(ab), DIM);   // unit tangent vector
    Vector3 n(fib.normal, DIM);    // normal direction
    
    // adjust normal direction:
    if ( n.normSqr() < 0.5 )
        n = d.orthogonal(1.0);
    n = d.orthogonal(n, 1.0);
    fib.normal.get(n, DIM);
    
    const real abmax = fib.abscissaP();

    int cnt = 0;
    while ( ab < abmax )
    {
        d = Vector3(fib.dir(ab), DIM);
        Vector3 p(fib.pos(ab), DIM);
        
        // adjust 'n' to keep it orthogonal to 'd':
        n = d.orthogonal(n, 1.0);
        
        // set two vectors orthogonal to 'd' of norm 'off':
        Vector3 e = n * off;
        Vector3 f = cross(d, e);
        
        // color of alpha subunit:
        dat[2] = 1;
        
        real up = abmax - ( cnt & 1 ? 0 : dab );
  
        for ( int i = 0; i < 13; ++i )
        {
            if ( ab + da[i] < up )
            {
                if ( cnt & 1 )
                {
                    // change color for beta-tubulin near the end:
                    if ( ab + da[i] + 4 * dab > abmax )
                        dat[2] = 3;
                    else
                        dat[2] = 2;
                }
                savePoint(dat, p+dx[i]*e+dy[i]*f+da[i]*d);
            }
        }
        
        ab += dab;
        ++cnt;
    }
}

//------------------------------------------------------------------------------

// radius of microtubule
const real rad = 0.0125;


void saveLink(Couple const* cop)
{
    int dat[3] = { 7, cop->identity(), 0 };
    float vec[6] = { 0 };
    
    Vector p1 = cop->posHand1();
    Vector p2 = cop->posHand2();

    Vector dir = p2 - p1;
    real dn = dir.norm();
    
    if ( dn > 1e-5 )
    {
        // position the heads at the surface of the filaments:
        p1 += dir * ( rad / dn );
        p2 -= dir * ( rad / dn );
    }

    p1.put(vec);
    p2.put(vec+3);

    if ( binary )
    {
        fwrite(dat, 3, sizeof(int), file);
        fwrite(vec, 6, sizeof(float), file);
    }
    else
    {
        fprintf(file, "%i %i %i ", dat[0], dat[1], dat[2]);
        fprintf(file, "%.6f %.6f %.6f %.6f %.6f %.6f\n", vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]);
    }
}

//------------------------------------------------------------------------------


int main(int argc, char* argv[])
{
    Simul simul;
    
    if ( argc > 1 && strstr(argv[1], "help") )
    {
        help(std::cout);
        return EXIT_SUCCESS;
    }

    Glossary opt;
    opt.read_strings(argc, argv);

    opt.set(binary, "binary");
    opt.set(style, "style");

    std::string input = simul.prop->trajectory_file;
    opt.set(input, ".cmo") || opt.set(input, "input");;

    FrameReader reader;
    RNG.seedTimer();

    try
    {
        simul.loadProperties();
        reader.openFile(input);
        
        unsigned frame = 0;
        char filename[256];
        
        // process all frames in the file:
        while ( 0 == reader.loadNextFrame(simul) )
        {
            // create file:
            if ( binary )
            {
                snprintf(filename, sizeof(filename), "frame%04i.bin", frame);
                file = fopen(filename, "wb");
            }
            else
            {
                snprintf(filename, sizeof(filename), "frame%04i.txt", frame);
                file = fopen(filename, "w");
            }
            if ( file == 0 || ferror(file) )
            {
                std::cerr << "Aborted: could not create file " << filename << std::endl;
                return EXIT_FAILURE;
            }
            
            // save fibers in the natural order:
            Fiber const* fib = simul.fibers.firstID();
            while ( fib )
            {
                //if ( ! binary ) fprintf(file, "%% fiber %li\n", fib->identity());
             
                int s = style;
                if ( fib->prop->disp )
                    s = fib->prop->disp->style;
                
                switch ( s )
                {
                    case 0: saveBackbone(*fib);    break;
                    case 1: saveActin(*fib);       break;
                    case 2: saveMicrotubule(*fib); break;
                }
                fib = simul.fibers.nextID(fib);
            }
            
            // save links:
            Couple const* cop = simul.couples.firstID();
            while ( cop )
            {
                if ( cop->linking() )
                    saveLink(cop);
                cop = simul.couples.nextID(cop);
            }

            fclose(file);
            file = 0;
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

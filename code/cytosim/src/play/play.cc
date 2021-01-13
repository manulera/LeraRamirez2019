// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "glossary.h"
#include "offscreen.h"
#include "saveimage.h"
#include "filepath.h"
#include "splash.h"

#include "opengl.h"
#include "player.cc"

using Player::simul;
using Player::simThread;
using Player::DP;
using Player::PP;
using std::endl;

extern bool functionKey[];

/// different modes:
enum { ONSCREEN, OFFSCREEN_IMAGE, OFFSCREEN_MOVIE, ONSCREEN_MOVIE };

//------------------------------------------------------------------------------

void help(std::ostream & os = std::cout)
{
    os << " Options can be specified on the command line while invoking cytosim.\n"
          "\n Global options:\n"
          "         dir=PATH                   change working directory to specified PATH\n"
          "         live                       enter live simulation mode directly\n"
          "         FILE.cym                   specify input config file\n"
          "         FILE.cmo                   specify trajectory file\n"
          "         FILE.cmi                   specify display configuration file\n"
          "         PARAMETER=value            set parameter value (example size=512)\n"
          "\n Rendering:\n"
          "         image frame=INT            render specified frame offscreen\n"
          "         image frame=INT,INT,...    render several frames offscreen\n"
          "         image magnify=INT          render frames at higher resolution\n"
          "         movie                      render all frames offscreen\n"
          "         movie=on                   render all frames onscreen\n"
          "         movie period=INT           render one frame every INT frames\n"
          "\n Help:\n"
          "         help                       print this help\n"
          "         parameters                 print a list of parameters\n"
          "                                    (there should be no space around the equal sign)\n";
}


/// kill the slave thread and free memory
void exit_handler()
{
    Player::simThread.stop();
    Player::dproperties.erase();
    gle::release();
    //fprintf(stderr, "Goodbye...\n");
}


/// swap multisample and normal buffers
void blitBuffers(GLuint normal, GLuint multi, GLint W, GLint H)
{
    //std::clog << "blitting multisample buffer\n";
    glBindFramebuffer(GL_READ_FRAMEBUFFER, multi);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, normal);
    glBlitFramebuffer(0, 0, W, H, 0, 0, W, H, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, multi);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, normal);
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    bool goLive = false;
    int  mode = ONSCREEN;
    int  magnify = 1;
    
    MSG.quiet();
    Glossary glos;
    
    try {
        glos.read_strings(argc-1, argv+1);
    }
    catch( Exception & e ) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    // check for major options:
    
    if ( glos.use_key("help") )
    {
        splash(std::cout);
        help();
        return EXIT_SUCCESS;
    }

    if ( glos.use_key("info") || glos.use_key("--version") )
    {
        splash(std::cout);
        print_version(std::cout);
        if ( SaveImage::supported("png") )
            std::cout << "   PNG enabled\n";
        std::cout << "   DIM = " << DIM << endl;
        return EXIT_SUCCESS;
    }
    
    if ( glos.use_key("parameters") )
    {
        splash(std::cout);
        Player::writePlayParameters(std::cout, false);
        return EXIT_SUCCESS;
    }

    if ( glos.use_key("live") )
        goLive = true;
    
    if ( glos.use_key("image") )
        mode = OFFSCREEN_IMAGE;

    if ( glos.use_key("poster") )
    {
        mode = OFFSCREEN_IMAGE;
        magnify = 3;
    }
    
    if ( glos.value_is("movie", 0, "on") )
        mode = ONSCREEN_MOVIE;
    else if ( glos.use_key("movie") )
        mode = OFFSCREEN_MOVIE;
    
    // get image over-sampling:
    glos.set(magnify, "magnify") || glos.set(magnify, "magnification");

    // initialize the first frame index
    int frm = 0;
    glos.set(frm, "frame");

    std::string dir;
    // change working directory if specified:
    if ( glos.set(dir, "dir") )
    {
        FilePath::change_dir(dir);
        //std::cerr << "Cytosim working directory is " << FilePath::get_dir() << std::endl;
    }
    
    glApp::setDimensionality(DIM);

    try
    {
        Player::simul.initialize(glos);
    }
    catch( Exception & e )
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    // secondary configuration file to adjust display parameters:
    std::string setup = goLive ? simul.prop->config : simul.prop->property_file;
    View& view = glApp::views[0];

    try
    {
        // obtain setup from the command line:
        glos.set(setup, "setup") || glos.set(setup, ".cmi") || glos.set(setup, ".cms");
        
        //std::cerr << "cytosim graphical setup file: " << setup << std::endl;
        
        // extract first specification of "simul:display" string from the setup file
        if ( FilePath::is_file(setup) )
            Parser(simul, 0, 0, 0, 0, 0).readConfig(setup);
        
        // settings in the file are read, but do not overwrite the command-line options:
        glos.read(simul.prop->display, 1);
        simul.prop->display_fresh = false;
        
        // read parameters:
        glApp::GP.read(glos);
        view.read(glos);
        DP.read(glos);
        PP.read(glos);
    }
    catch( Exception & e )
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    /*
     At exit, exit_handler() will kill the slave thread, 
     before any object is destroyed to prevent the slave thread
     from accessing corrupted data.
     */
    atexit(exit_handler);
    
    //---------Open trajectory file in playback mode

    if ( ! goLive || frm > 0 )
    {
        try
        {
            std::string file = simul.prop->property_file;
            
            if ( !FilePath::is_file(file) )
            {
                std::cerr << "Aborted: could not find Cytosim file `" << file << "'\n";
                return EXIT_FAILURE;
            }

            Parser(simul, 1, 1, 0, 0, 0).readConfig(file);
            
            // read 'setup' file again to overwrite values set in 'properties_file'
            if ( file != setup )
                Parser(simul, 0, 0, 0, 0, 0).readConfig(setup);
            
            simThread.openFile();
            
            // load requested frame from trajectory file:
            if ( simThread.loadFrame(frm) )
            {
                // if eof happened, reload last frame in file:
                simThread.loadFrame(-1);
                if ( simThread.currFrame() > 0 )
                {
                    frm = simThread.currFrame();
                    std::cerr << "Warning: could only load frame " << frm << '\n';
                }
            }
        }
        catch( Exception & e ) {
            std::cerr << "\nError: " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
    }
    
#ifdef __linux__
    // it is necessary under Linux to initialize GLUT to display fonts
    glutInit(&argc, argv);
#endif
    
    //-------- off-screen (non interactive) rendering -------
    
    if ( mode == OFFSCREEN_IMAGE || mode == OFFSCREEN_MOVIE )
    {
        const int W = view.width() * magnify;
        const int H = view.height() * magnify;
        
        //std::cerr << W << "x" << H << << " downsample " << PP.downsample << '\n';
        
        if ( !OffScreen::openContext() )
        {
            std::cerr << "Failed to create off-screen context" << std::endl;
            return EXIT_FAILURE;
        }
        GLuint fbo = OffScreen::createBuffer(W, H, 0);
        if ( !fbo )
        {
            std::cerr << "Failed to create off-screen pixels" << std::endl;
            return EXIT_FAILURE;
        }
        GLuint multi = 0;        
        if ( view.multisample > 1 )
        {
            multi = OffScreen::createBuffer(W, H, view.multisample);
        }
        
        gle::initialize();
        Player::initStyle(DP.style);

        // prepare View:
        //view.reshape(W, H);
        view.initGL();

        if ( mode == OFFSCREEN_IMAGE )
        {
            unsigned inx = 0;
            // it is possible to specify multiple frame indices:
            do {
                Player::readFrame(frm);
                // only save requested frames:
                if ( simThread.currFrame() == frm )
                {
                    Player::displayScene(magnify);
                    if ( multi )
                        blitBuffers(fbo, multi, W, H);
                    if ( magnify > 1 )
                        Player::saveImage("poster", frm);
                    else
                        Player::saveImage("image", frm);
                }
            } while ( glos.set(frm, "frame", ++inx) );
        }
        else if ( mode == OFFSCREEN_MOVIE )
        {
            // save every PP.period
            unsigned s = PP.period;
            do {
                if ( ++s >= PP.period )
                {
                    Player::displayScene(magnify);
                    if ( multi )
                        blitBuffers(fbo, multi, W, H);
                    Player::saveImage("movie", frm++);
                    s = 0;
                }
            } while ( 0 == simThread.loadNextFrame() );
        }
        
        if ( multi )
            OffScreen::releaseBuffer();
        OffScreen::releaseBuffer();
        OffScreen::closeContext();
        glos.warnings(std::cerr);
        return EXIT_SUCCESS;
    }
    
    glos.warnings(std::cerr);

    //--------- initialize and open window

#ifndef __linux__
    glutInit(&argc, argv);
#endif
    
    //link the Function-Keys to allow user controls:
    glApp::bindFunctionKeys(functionKey);
    //register all the GLUT callback functions:
    glApp::actionFunc(Player::processMouseClick);
    glApp::actionFunc(Player::processMouseDrag);
    glApp::normalKeyFunc(Player::processNormalKey);
    glApp::createWindow(Player::displayLive);
    
    if ( mode == ONSCREEN_MOVIE )
    {
        PP.exit_at_eof = true;
        PP.save_images = true;
        PP.play = 1;
    }
    
    //-------- initialize graphics -------

    try
    {
        gle::initialize();
        Player::initStyle(DP.style);
        Player::buildMenus();
        
        //register a first timer callback:
        glutTimerFunc(100, Player::timerCallback, 0);
    }
    catch ( Exception & e )
    {
        printf("Initialization error: %s\n", e.what());
        return EXIT_FAILURE;
    }
    
    
    if ( goLive )
    {
        try
        {
            // read 'setup' file first, but values will be overwritten by 'config.cym'
            if ( FilePath::is_file(setup) )
                Parser(simul, 0, 0, 0, 0, 0).readConfig(setup);
        
            //set PP.live to drive the timer in play:
            PP.live = 1;
            simThread.period(PP.period);
            
            int res = 0;
            if ( frm > 0 )
                res = simThread.extend();
            else
                res = simThread.start();
            
            if ( res )
            {
                std::cerr << "Error: could not start live thread" << std::endl;
                return EXIT_FAILURE;
            }
        }
        catch( Exception & e ) {
            std::cerr << "\nError: " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
    }
    
    //start the GLUT event handler:
    glutMainLoop();

    return EXIT_SUCCESS;
}

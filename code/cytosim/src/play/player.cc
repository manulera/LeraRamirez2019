// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "player.h"
#include "opengl.h"
#include "iowrapper.h"
#include "filepath.h"
#include "gle.h"
#include "glapp.h"

using namespace gle;
using glApp::flashText;

#include "fiber_disp.h"
#include "line_disp.h"
#include "point_disp.h"

#include "player_disp.cc"
#include "player_keys.cc"
#include "player_mouse.cc"
#include "player_menus.cc"

//------------------------------------------------------------------------------
#pragma mark I/O

void Player::readFrame(int f)
{
    simThread.loadFrame(f);
}


void Player::previousFrame()
{
    if ( simThread.currFrame() > 0 )
        simThread.loadFrame(simThread.currFrame()-1);
    else {
        if ( PP.loop )
            simThread.loadFrame(-1);
        else
            stop();
    }
}

/**
 Reads the next frame from the current file position.
 */
void Player::nextFrame()
{    
    try
    {
        if ( simThread.loadNextFrame() )
        {
            if ( PP.exit_at_eof )
                exit(EXIT_SUCCESS);
            if ( PP.loop )
                simThread.loadFrame(0);
            else
            {
                flashText("end-of-file\n");
                stop();
            }
        }
    }
    catch( Exception & e )
    {
        flashText("Error:\n %s", e.what());
        if ( simThread.eof() )
            stop();
    }
}

//------------------------------------------------------------------------------
#pragma mark Commands


void Player::rewind()
{
    if ( simThread.goodFile() )
    {
        simThread.loadFrame(0);
        stop();
    }
}


bool Player::startPlayback()
{
    if ( PP.play != 1  && ! PP.live )
    {
        //rewind if the end of the file was reached:
        if ( simThread.eof() )
            simThread.loadFrame(0);
        PP.play = 1;
        return true;
    }
    return false;
}


bool Player::startBackward()
{
    if ( PP.play != -1 )
    {
        if ( simThread.currFrame() == 0 )
            readFrame(-1);
        else
            flashText("Play reverse");
        PP.play = -1;
        return true;
    }
    return false;
}


void Player::accelerate()
{
    PP.delay /= 2;
    //the delay should be compatible with graphic refresh rates:
    const unsigned int min_delay = 1;
    if ( PP.delay < min_delay )
    {
        PP.delay = min_delay;
        if ( PP.live )
            flashText("Delay is %i ms! use 'A' to jump frames", PP.delay);
        else
            flashText("Delay is %i ms!", PP.delay);
    }
    else {
        flashText("Delay %i ms", PP.delay);
    }
}


void Player::stop()
{
    PP.play = 0;
    PP.live = 0;
    PP.save_images = 0;
}


void Player::startstop()
{
    if ( PP.live )
        PP.live = 0;
    else if ( simThread.goodFile() )
    {
        if ( !PP.play )
            startPlayback();
        else
            stop();
    }
    else
        PP.live = 1;
}


void Player::startLive()
{
    glApp::displayFunc(displayLive);
    if ( 0 == simThread.extend() )
        flashText("Extend simulation...");
    PP.live = 1;
}


void Player::stepLive()
{
    simThread.signal();
}


void Player::restartLive()
{
    try
    {
        simThread.stop();
        simThread.clear();
        dproperties.erase();
        FDisp = 0;
        simThread.start();
    }
    catch( Exception & e ) {
        flashText("Error: %s", e.what());
    }
}


//------------------------------------------------------------------------------
#pragma mark -


/**
 set current FiberDisp pointer FDisp
 */
void Player::setPointers(bool next)
{
    Property * val = 0;
    if ( FDisp == 0 )
        val = dproperties.find_next("fiber:display", 0);
    
    if ( next )
    {
        // change FD, allowing access to different FiberDisp
        val = dproperties.find_next("fiber:display", FDisp);
    }
    
    if ( val  &&  val != FDisp )
    {
        if ( next )
            flashText("Addressing `%s'", val->name().c_str());
        FDisp = static_cast<FiberDisp*>(val);
    }
}


/**
 Write global parameters that control the display:
 - GlappProp
 - DisplayProp
 - PlayProp
 .
 */
void Player::writePlayParameters(std::ostream & os, bool prune)
{
    os << "set simul:display *" << std::endl;
    os << "{" << std::endl;
    if ( glApp::views.size() > 0 )
    {
        View& view = glApp::currentView();
        view.write_diff(os, prune);
    }
    DP.write_diff(os, prune);
    GP.write_diff(os, prune);
    //output parameters for the main view:
    PP.write_diff(os, prune);
    os << "}" << std::endl;
}

/**
 Write all the parameters that control the display:
 - GlappProp
 - DisplayProp
 - PlayProp
 - ObjectDisp
 .
 */
void Player::writeDisplayParameters(std::ostream & os, bool prune)
{
    dproperties.write(os, prune);
}


//------------------------------------------------------------------------------
#pragma mark -

void Player::timerCallback(const int value)
{
    //std::clog << "Player::timer " << value << std::endl;
    //glApp::postRedisplay();

    if ( PP.live && simThread.alive() )
    {
        if ( PP.save_images && 0 == simThread.trylock() )
        {
            displayScene(1);
            saveImage("image", PP.image_index++);
            simThread.unlock();
        }
        
        simThread.signal();
    }
    else if ( PP.play )
    {
        if ( PP.save_images )
        {
            displayScene(1);
            saveImage("movie", simThread.currFrame());
        }

        if ( PP.play == 1 )
        {
            // skip PP.period frames, and at least one
            for ( unsigned s = 1; s < PP.period; ++s )
                nextFrame();
            nextFrame();
        }
        else if ( PP.play == -1 )
            previousFrame();
        
        glApp::postRedisplay();
    }
    
    /*
     Register the next timer callback
     in idle mode, we use a long time-interval
     */
    if ( PP.play || PP.live )
        glutTimerFunc(PP.delay, timerCallback, 2);
    else
        glutTimerFunc(30, timerCallback, 1);
}



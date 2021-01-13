// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef PLAYER_H
#define PLAYER_H

///this turns on display code in simul.h
#define  DISPLAY

#include "glapp.h"
#include "simul.h"
#include "parser.h"
#include "display.h"
#include "sim_thread.h"
#include "display_prop.h"
#include "play_prop.h"

class FiberDisp;


///GUI for Cytosim. Display is done by Display, most mouse handling by glApp
namespace Player
{
    /// container for display properties of objects
    PropertyList dproperties;
    
    /// the display parameters
    DisplayProp DP("*");
    
    /// the parameters for play
    PlayProp  PP("*");
    
    /// the 'current' FiberDisp on which any change is applied
    FiberDisp * FDisp = 0;
    
    /// SimThread to control the live simulation
    SimThread simThread(glApp::postRedisplay);
    
    /// an alias to current Simul
    Simul& simul = simThread.sim();
    
    /// the Display object
    Display * mDisplay = 0;
    
    /// set FDisp;
    void setPointers(bool);

    //---------------------------------COMMANDS---------------------------------
    
    /// reset view, without changing the current frame
    void rewind();
   
    /// start animation
    bool startPlayback();
    
    /// accelerate animation if possible
    void accelerate();
    
    /// start reverse animation
    bool startBackward();
    
    /// start live simulation
    void startLive();

    /// stop animation
    void stop();

    /// start or stop animation
    void startstop();
    
    /// allow one step of the simulation engine
    void stepLive();
    
    /// reset the sim-state and timer
    void restartLive();

    /// load specified frame
    void readFrame(int);

    /// load previous frame
    void previousFrame();
    
    /// go to the next frame, returns 1 if EOF is encountered
    void nextFrame();
    
    /// write global display parameters
    void writePlayParameters(std::ostream & out, bool prune);

    /// write Object display parameters
    void writeDisplayParameters(std::ostream & out, bool prune);

    //-----------------------------------GUI----------------------------------
    
    /// returns a string with some help on what pressing keys do
    void helpKeys(std::ostream &);
    
    /// GLUT callback function for timed events
    void timerCallback(int);

    /// build all the menus from scratch
    void buildMenus();
    
    /// GLUT callback function for most keys
    void processNormalKey(unsigned char, int = 0, int = 0);

    /// GLUT callback function for mouse motion, when one button is pressed
    void processMouseClick(int, int, const Vector3&, int);
    
    /// GLUT callback function for mouse motion, when no button is pressed
    void processMouseDrag(int, int, Vector3&, const Vector3&, int);
    
    //-----------------------------DISPLAY------------------------------------
    
    /// initialize display with given style
    void initStyle(int);
    
    /// build message that appears on top
    std::string buildMessageT(Simul const&);
    
    /// build message that appears on bottom
    std::string buildMessageB(Simul const&);
    
    /// build central message
    std::string buildMessage(int);

    
    /// set View::focus and quat to match the center of gravity of the Fibers
    void autoTrack(FiberSet const&, View&);
    
    /// adjust the viewing area
    int  autoScale(SpaceSet const&, View&);

    /// adjust the model view and load frame if asked
    void prepareDisplay(View&, Simul const&, real mag=1);
    
    /// display cytosim state and message
    void displayCytosim(Simul const&);
    
    /// display function for on-screen rendering
    void displayLive();

    /// display function for off-screen rendering
    void displayScene(int mag);
    
    
    /// save the displayed image in a graphic file
    int  saveImage(const char * root, unsigned indx);
    
    /// save high-resolution image using composite method
    int  saveMagnifiedImage(int mag, const char* filename, const char* format, int downsample=1);
    
    /// save a high-resolution image using composite method
    int  saveMagnifiedImage(int mag, const char* root, unsigned indx, int downsample=1);

}


#endif

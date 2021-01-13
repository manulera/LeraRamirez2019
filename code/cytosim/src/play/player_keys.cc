// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "fiber_prop.h"
#include "fiber_disp.h"
#include "glut.h"

/**
 It is possible to save images or many images by pressing a key..
 disabling this capacity in public distribution might be safer:
 It can be done by disabling ENABLE_WRITE below:
 */
#define ENABLE_WRITE


/// amount added/subtracted to size/width when a key is pressed
const real GRAIN = 0.5;

//------------------------------------------------------------------------------

/**
 Change size for all visible PointDisp.
 This sets PointDisp:width to `s` and PointDisp::size to `2*s`.
 */
void setPointDispSizeWidth(PropertyList const& plist, DisplayProp& DP, real s)
{
    for ( PropertyList::const_iterator pi = plist.begin(); pi < plist.end(); ++pi )
    {
        PointDisp * dsp = static_cast<PointDisp*>(*pi);
        if ( dsp->visible )
        {
            dsp->size = 2*s;
            dsp->width = s;
        }
    }
}


/**
 Change size (add inc) for all visible PointDisp.
*/
void changePointDispSize(PropertyList const& plist, DisplayProp& DP, real inc)
{
    for ( PropertyList::const_iterator pi = plist.begin(); pi < plist.end(); ++pi )
    {
        PointDisp * dsp = static_cast<PointDisp*>(*pi);
        if ( dsp->visible )
        {
            real s = inc * ( 1 + round(dsp->size/inc) );
            if ( s > 0 )
                dsp->size = s;
        }
    }
    
    // also change the global value:
    real s = inc * ( 1 + round(DP.point_size/inc) );
    if ( s > 0 )
        DP.point_size = s;
}

/**
 Change width (add inc) for all visible PointDisp.
 */
void changePointDispWidth(PropertyList const& plist, DisplayProp& DP, real inc)
{
    for ( PropertyList::const_iterator pi = plist.begin(); pi < plist.end(); ++pi )
    {
        PointDisp * dsp = static_cast<PointDisp*>(*pi);
        if ( dsp->visible )
        {
            real s = dsp->width + inc;
            if ( s > 0 )
                dsp->width = s;
        }
    }
    
    // also change the global value:
    real s = DP.line_width + inc;
    if ( s > 0 )
       DP.line_width = s;
}


void changePointDispStyle(PointDisp* pod)
{
    if ( pod == 0 )
        return;

    int & style = pod->style;
    style = ( style + 1 ) % 8;
    flashText("%s: style = %i", pod->name().c_str(), style);
}


void changePointDispStyle(const PropertyList& plist)
{
    for ( PropertyList::const_iterator pi = plist.begin(); pi < plist.end(); ++pi )
        changePointDispStyle(static_cast<PointDisp*>(*pi));
}


void setVisible(const PropertyList& plist, int val)
{
    for ( PropertyList::const_iterator pi = plist.begin(); pi < plist.end(); ++pi )
    {
        PointDisp * dsp = static_cast<PointDisp*>(*pi);
        dsp->visible = val;
    }
    
    if ( val )
        flashText("All hands are visible");
    else
        flashText("No hand is visible");
}


void shuffleVisible(const PropertyList& plist)
{
    bool all_on = true;
    PointDisp * on = 0;
    
    // find first one which is visible:
    for ( PropertyList::const_iterator pi = plist.begin(); pi < plist.end(); ++pi )
    {
        PointDisp * dsp = static_cast<PointDisp*>(*pi);
        if ( dsp->visible )
        {
            // choose follower:
            if ( !on && pi+1 < plist.end() )
                on = static_cast<PointDisp*>(*(pi+1));
        }
        else
            all_on = false;
    }
    
    if ( all_on )
    {
        // choose first:
        if ( plist.size() )
            on = static_cast<PointDisp*>(*plist.begin());
    }

    bool set = !on;
    
    for ( PropertyList::const_iterator pi = plist.begin(); pi < plist.end(); ++pi )
        static_cast<PointDisp*>(*pi)->visible = set;

    if ( on )
    {
        on->visible = true;
        flashText("Only `%s' is visible", on->name().c_str());
    }
    else
    {
        flashText("All hands are visible");
    }
}


void changeSingleSelect(DisplayProp& DP)
{
    unsigned int & select = DP.single_select;
    switch( select )
    {
        case 3:   select = 0;    flashText("single:select=0: hidden");      break;
        case 0:   select = 2;    flashText("single:select=2: bound only");  break;
        case 2:   select = 1;    flashText("single:select=1: free only");   break;
        default:  select = 3;    flashText("single:select=3: all");         break;
    }
}


void changeCoupleSelect(DisplayProp& DP)
{
    unsigned int & select = DP.couple_select;
    switch( select )
    {
        case 7:   select = 0;    flashText("couple:select=0: hidden");        break;
        case 0:   select = 2;    flashText("couple:select=2: bound only");    break;
        case 2:   select = 4;    flashText("couple:select=4: bridging only"); break;
        case 4:   select = 1;    flashText("couple:select=1: free only");     break;
        default:  select = 7;    flashText("couple:select=7: all");           break;
    }
}

void changeCoupleSelect2(DisplayProp& DP)
{
    unsigned int & select = DP.couple_select;
    if ( select & 8 )
    {
        select = 16+4;
        flashText("couple: parallel bridging only");
    }
    else if ( select & 16 )
    {
        select = 7;
        flashText("couple: all");
    }
    else
    {
        select = 8+4;
        flashText("couple: antiparallel bridging only");
    }
}

//---------------------------------------------------------------------
#pragma mark - Fibers

void changeExclude(FiberDisp* fd, bool mod)
{
    if ( mod )
        fd->exclude >>= 2;
    fd->exclude = ( fd->exclude + 1 ) % 4;
    if ( mod )
        fd->exclude <<= 2;
    
    switch( fd->exclude )
    {
        case 0: flashText("All fibers");                break;
        case 1: flashText("Right-pointing fibers");     break;
        case 2: flashText("Left-pointing fibers");      break;
        case 3: flashText("No fibers");                 break;
        case 4: flashText("Counter-clockwise fibers");  break;
        case 8: flashText("Clockwise fibers");          break;
        case 12: flashText("No fibers");                break;
    }
}


void changeExplode(FiberDisp* fd)
{
    fd->explode = ! fd->explode;
    flashText("Fiber:explode = %i", fd->explode);
}


void changeColoring(FiberDisp* fd)
{
    fd->coloring = ( fd->coloring + 1 ) % 6;
    switch( fd->coloring )
    {
        case FiberDisp::COLORING_OFF:       flashText("Fibers: no coloring");           break;
        case FiberDisp::COLORING_RANDOM:    flashText("Fibers: coloring");              break;
        case FiberDisp::COLORING_DIRECTION: flashText("Fibers: coloring by direction"); break;
        case FiberDisp::COLORING_MARK:      flashText("Fibers: coloring by mark");      break;
        case FiberDisp::COLORING_FLECK:     flashText("Fibers: coloring by cluster");   break;
        case FiberDisp::COLORING_AGE:       flashText("Fibers: coloring by age");   break;
    }
}


void changeMask(FiberDisp* fd)
{
    fd->mask_bitfield = RNG.pint_bits(fd->mask);
    flashText("fiber:mask_bitfield=0x%X (%i bits)", fd->mask_bitfield, fd->mask);
}


void increaseMask(FiberDisp* fd)
{
    // increase mask, but only allow up to 10 bits to be on:
    fd->mask = ( fd->mask + 1 ) % 11;
    changeMask(fd);
}


void changePointStyle(FiberDisp* fd)
{
    fd->point_style = ( fd->point_style + 1 ) % 3;
    switch ( fd->point_style )
    {
        case 0:  flashText("Fibers: no points");     break;
        case 1:  flashText("Fibers: model points");  break;
        case 2:  flashText("Fibers: arrowheads");    break;
        case 3:  flashText("Fibers: center point");  break;
    }
}


void changeLineStyle(FiberDisp* fd)
{
    fd->line_style = ( fd->line_style + 1 ) % 5;
    switch ( fd->line_style )
    {
        case 0:  flashText("Fibers: no lines");       break;
        case 1:  flashText("Fibers: lines");          break;
        case 2:  flashText("Fibers: axial tensions"); break;
        case 3:  flashText("Fibers: curvature");      break;
        case 4:  flashText("Fibers: orientation");    break;
    }
}


void changeSpeckleStyle(FiberDisp* fd)
{
    fd->speckle_style = ( fd->speckle_style + 1 ) % 3;
    switch ( fd->speckle_style )
    {
        case 0:  flashText("Fibers: no speckles");       break;
        case 1:  flashText("Fibers: random speckles");   break;
        case 2:  flashText("Fibers: regular speckles");  break;
    }
}


void changeLatticeStyle(FiberDisp* fd)
{
    fd->lattice_style = ( 1 + fd->lattice_style ) % 3;
    flashText("Fibers: lattice_style=%i", fd->lattice_style);
}


void changePointSize(FiberDisp* fd, int dir, real grain)
{
    int s = dir + round( fd->point_size / grain );
    
    if ( s > 0 )
        fd->point_size = grain * s;
    
    flashText("%s: point_size=%0.2f", fd->name().c_str(), fd->point_size);
}


void changeLineWidth(FiberDisp* fd, int dir, real grain, bool scale_pointsize)
{
    real & width = fd->line_width;
    int s = dir + round( width / grain );
    
    if ( s > 0 )
    {
        if ( scale_pointsize  &&  width > 0 )
        {
            real sc = s * grain / width;
            fd->point_size  *= sc;
            fd->end_size[0] *= sc;
            fd->end_size[1] *= sc;
        }
        width = grain * s;
        flashText("Fibers: line_width %0.2f", width);
    }
}

//------------------------- Fibers ends ------------------------------

void changeTipStyle(FiberDisp* fd)
{
    int * style = fd->end_style;
    // showing the plus ends -> the minus ends -> both -> none
    switch( (style[1]?1:0) + (style[0]?2:0) )
    {
        case 0:
            style[0] = 2;
            style[1] = 0;
            break;
        case 1:
            style[0] = 0;
            style[1] = 0;
            break;
        case 2:
            style[0] = 2;
            style[1] = 4;
            break;
        case 3:
        default:
            style[0] = 0;
            style[1] = 4;
            break;
    }

    switch( (style[0]?1:0) + (style[1]?2:0) )
    {
        case 0: flashText("Fibers: no ends");    break;
        case 1: flashText("Fibers: plus-ends");  break;
        case 2: flashText("Fibers: minus-ends"); break;
        case 3: flashText("Fibers: both ends");  break;
    }
}


void changeTipSize(FiberDisp* fd, real inc)
{
    real * size = fd->end_size;
    if ( size[0] + 2*inc > 0 ) size[0] += 2*inc;
    if ( size[1] + inc   > 0 ) size[1] += inc;
    flashText("Fibers: end_size %.1f %.1f", size[0], size[1]);
}


//------------------------------------------------------------------------------
//---------------------------- keyboard commands -------------------------------
//------------------------------------------------------------------------------
#pragma mark - Keyboard Commands

/// provide minimal on-screen summary of the most important key combinations
void Player::helpKeys(std::ostream & os)
{
    os << "                          Keyboard Commands\n";
    os << "\n";
    os << "    SPACE      Start-stop animation or replay\n";
    os << "    < >        Show previous; show next frame ( , . also works)\n";
    os << "    u i o p    Play reverse; stop; play slower; play faster\n";
    os << "    z          Rewind to first frame / Restart live simulation\n";
    os << "    ALT-SPACE  Reset view (i.e. zoom, translation, rotation)\n";
    os << "    f F        Toggle full-screen mode; maximize window size\n";
    os << "    CMD-Q      Quit\n";
    os << "\nFibers\n";
    os << "    `          Address another type of fibers for modifications\n";
    os << "    1          Change display: line / color-coded tension / hide\n";
    os << "    3 4        Decrease; Increase line width (ALT: point size)\n";
    os << "    !          Change tip display: none / plus / both / minus\n";
    os << "    # $        Decrease; increase fiber_end display size\n";
    os << "    2 @        Change speckle display; change lattice style\n";
    os << "    c d w      Coloring, Right/left-pointing, Fractional masking\n";
    os << "    t T        Toggle auto-tracking: 't':nematic; 'T':polar mode\n";
    os << "\nBeads - Solids - Spheres\n";
    os << "    5          Switch between different bead/sphere display style\n";
    os << "    %          Change first bead/solid display style\n";
    os << "    * (        Decrease; increase point size\n";
    os << "\nSingles - Couples\n";
    os << "    u          Toggle Hands visibility flags\n";
    os << "    6          Change Single selection mode\n";
    os << "    7 8        Both keys change the Couple selection mode;\n";
    os << "    9 0        Decrease; Increase point size\n";
    os << "    ( )        Decrease; Increase line width\n";
    os << "\nSimulation\n";
    os << "    a s        Start live mode; Perform one simulation step\n";
    os << "    A S        Double nb-steps/display; set nb-steps/display = 1\n";
    os << "    g G        Delete all mouse-controlled hands; Release current hand\n";
    os << "\nInput/Output\n";
    os << "    r          Read parameter file and update simulation\n";
    os << "    R          Write display parameters to terminal\n";
#ifdef ENABLE_WRITE
    os << "    y Y        Save current image; Play and save all images in file\n";
#endif
}


//------------------------------------------------------------------------------
void Player::processNormalKey(const unsigned char key, const int x, const int y)
{
    //std::cerr << "processing key `" << key << "'\n";
    
    setPointers(0);
    // the view associated with the current window
    View & view = glApp::currentView();

    // execute the custom piece of code (magic_key / magic_code)
    for ( int k = 0; k < PlayProp::NB_MAGIC_KEYS; ++k )
    {
        if ( key == PP.magic_key[k] )
        {
            std::istringstream iss(PP.magic_code[k]);
            simThread.execute(iss);
            glApp::postRedisplay();
            return;
        }
    }
    
    const bool altKeyDown = glutGetModifiers() & GLUT_ACTIVE_ALT;
    /*
     In the switch below:
     - use break if the display need to be refreshed,
     - otherwise, use return.
    */
    switch (key)
    {
        
        case 'h':
        {
            GP.show_message = ( GP.show_message + 1 ) % 6;
            GP.message = buildMessage(GP.show_message);
        } break;
            
#ifdef ENABLE_WRITE
            
        case 'y':
            // save current image
            displayScene(1);
            saveImage("image", PP.image_index++);
            //saveMagnifiedImage(3, "image", PP.image_index++, 3);
            return;
            
        case 'Y':
            // start to save all images
            if ( PP.save_images == 0 )
            {
                PP.save_images = 1;
                startPlayback();
            }
            else
            {
                PP.save_images = 0;
            }
            break;
            
#endif

        //------------------------- Global controls -----------------------------------
        
        case 'r': {
            try {
                simThread.reloadParameters();
                flashText("Parameters reloaded");
            }
            catch( Exception & e ) {
                flashText("Error in config: %s", e.what());
                PP.live = 0;
            }
        } break;

        case 'R':
        {
            if ( altKeyDown )
                simThread.writeProperties(std::cout, true);
            else
            {
                writePlayParameters(std::cout, true);
                std::cout << std::endl;
                writeDisplayParameters(std::cout, true);
            }
        } break;
            
        case 'z':
        case 'Z':
            if ( simThread.goodFile() )
                rewind();
            else
                restartLive();
            break;
            
        case 'a':
            if ( altKeyDown )
            {
                initStyle(1);
                flashText("Style 1");
            }
            else
            {
                startLive();
            }
            break;
            
        case 'A':
            PP.period = 2 * PP.period;
            simThread.period(PP.period);
            break;
            
        case 's':
            if ( altKeyDown )
            {
                initStyle(2);
                flashText("Style 2");
            }
            else
            {
                if ( simThread.alive() )
                    simThread.signal();
                stop();
            }
            break;
            
        case 'S':
            PP.period = 1;
            simThread.period(PP.period);
            break;
            
        case 'G':
            simThread.releaseHandle();
            break;
            
        case 'g':
            simThread.deleteHandles();
            flashText("Deleted mouse-controled handles");
            break;
            
        //------------------------- play / stop / reverse -----------------------------
            
        case '<':
        case ',':
            if ( PP.play == 1 )
                stop();
            else
                previousFrame();
            break;
            
        case '>':
        case '.':
            if ( PP.play == -1 )
                stop();
            else
                nextFrame();
            break;
            
        case 'i':
            PP.toggleReport(0);
            break;
            
        case 'I':
            PP.toggleReport(1);
            break;
            
        case 'o':
            if ( PP.delay < 1 << 13 )
                PP.delay *= 2;
            flashText("Delay %i ms", PP.delay);
            return;
            
        case 'O':
            if ( !startBackward() )
                accelerate();
            return;
            
        case 'p':
            if ( !startPlayback() )
                accelerate();
            return;
        
        case ' ':
            if ( altKeyDown )
            {
                glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
                view.reset();
                flashText("");
                break;
            }
            startstop();
            return;
            
        //------------------------------ Fibers ---------------------------------------
           
        case '`': {
            setPointers(1);
        } break;
            
        case 't':
            view.track_fibers &= 1;
            flashText("view.track_fibers = %i (translation)", view.track_fibers);
            break;
            
        case 'T':
            view.track_fibers &= 3;
            flashText("view.track_fibers = %i (rotation)", view.track_fibers);
            break;
            
        case 'd':
            if ( altKeyDown )
            {
                initStyle(3);
                flashText("Style 3");
            }
            else
            {
                if ( FDisp )
                    changeExclude(FDisp, 0);
            }
            break;
            
        case 'D':
            if ( FDisp )
                changeExclude(FDisp, 1);
            break;
            
        case 'e':
            if ( FDisp )
                changeExplode(FDisp);
            break;
            
        case 'W':
            if ( FDisp )
                changeMask(FDisp);
            break;
            
        case 'w':
            if ( FDisp )
                increaseMask(FDisp);
            break;

        case 'c':
            if ( FDisp )
                changeColoring(FDisp);
            break;
            
        case '1':
            if ( FDisp )
            {
                if ( altKeyDown )
                    changePointStyle(FDisp);
                else 
                    changeLineStyle(FDisp);
            }
            break;
            
        case '2':
            if ( FDisp )
                changeSpeckleStyle(FDisp);
            break;
            
        case '@':
            if ( FDisp )
                changeLatticeStyle(FDisp);
            break;
            
        case '3':
            if ( FDisp )
            {
                if ( altKeyDown)
                    changePointSize(FDisp, -1, GRAIN);
                else
                    changeLineWidth(FDisp, -1, GRAIN, true);
            }
            break;
            
        case '4':
            if ( FDisp )
            {
                if ( altKeyDown )
                    changePointSize(FDisp, +1, GRAIN);
                else
                    changeLineWidth(FDisp, +1, GRAIN, true);
            }
            break;
            
        case '!':
            if ( FDisp )
                changeTipStyle(FDisp);
            break;
            
        case '#':
            if ( FDisp )
                changeTipSize(FDisp, -GRAIN);
            break;
                
        case '$':
            if ( FDisp )
                changeTipSize(FDisp, +GRAIN);
            break;
            
        //------------------------ Solid / Sphere -------------------------------------
  
        case '5':
            changePointDispStyle(dproperties.find_all("bead:display", "solid:display", "sphere:display"));
            break;
        
        case '%':
        {
            DP.point_size *= 2;
            if ( DP.point_size > 8 )
                DP.point_size = 0.5;
            setPointDispSizeWidth(dproperties.find_all("bead:display", "solid:display", "sphere:display"), DP, DP.point_size);
            flashText("Point size %.1f", DP.point_size);
        }  break;
        
        //------------------------ Single/Couple + Hands ------------------------------
           
        case '6':
            changeSingleSelect(DP);
            break;
            
        case '7':
            changeCoupleSelect2(DP);
            break;

        case '8':
            changeCoupleSelect(DP);
            break;
            
        case 'u':
            shuffleVisible(dproperties.find_all("hand:display"));
            break;
            
        case 'U':
            setVisible(dproperties.find_all("hand:display"), 1);
            break;

        case '9':
            changePointDispSize(dproperties.find_all("hand:display"), DP, -GRAIN);
            flashText("Point size %.1f", DP.point_size);
            break;
            
        case '0':
            changePointDispSize(dproperties.find_all("hand:display"), DP, +GRAIN);
            flashText("Point size %.1f", DP.point_size);
            break;
            
        case '(':
            changePointDispWidth(dproperties.find_all("hand:display"), DP, -GRAIN);
            flashText("Line width %.1f", DP.line_width);
            break;
            
        case ')':
            changePointDispWidth(dproperties.find_all("hand:display"), DP, +GRAIN);
            flashText("Line width %.1f", DP.line_width);
            break;
            
        default:
            // other keys are handles by glApp
            glApp::processNormalKey(key, 0, 0);
            return;
    }
    
    // if break was called, redraw the scene:
    glApp::postRedisplay();
    // rebuild the menus, that might have changed:
    buildMenus();
}


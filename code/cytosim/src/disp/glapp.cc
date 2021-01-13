// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// F. Nedelec started glApp in Dec 2005

#include "glut.h"
#include "glapp.h"
#include "tictoc.h"
#include "filepath.h"
#include <cstdarg>
#include "exceptions.h"
#include "saveimage.h"
#include "glossary.h"
#include "gle.h"


using namespace gle;


namespace glApp
{    
    /*
     Most of the variables below are internal to glApp, and do not need to
     be accessed by the user. This is why they are only declared here and not
     in the header file.
    */
    
    glAppProp  GP("*");
    
    std::vector<View> views;
    
    int        mDIM = 3;            ///< current dimensionality
    int        mZoomWindow = -1;    ///< Window-ID for detailed view
    int        specialKeys = 0;     ///< state of special keys given by GLUT
    
    bool       fKeyArray[17];       ///< default array of function-key states
    bool     * fKey = fKeyArray;    ///< pointer to function-key states changed by bindFunctionKeys()
    char       fKeyString[4*17];    ///< String to display function-key states
    void       buildfKeyString();   ///< build fKeyString[]
    
    // --------------- MOUSE
    
    /// actions that can be performed with the mouse
    enum UserMode
    {
        MOUSE_ROTATE,
        MOUSE_TRANSLATE,
        MOUSE_ACTIVE,
        MOUSE_TRANSLATEZ,
        MOUSE_SPIN,
        MOUSE_ZOOM,
        MOUSE_MAGNIFIER,
        MOUSE_SET_ROI,
        MOUSE_MOVE_ROI,
        MOUSE_SELECT,
        MOUSE_PASSIVE
    };
    
    /// Specifies in which dimensionality each action is valid
    int actionValidity[] = { 3, 1, 1, 3, 2, 1, 1, 1, 4, 4, 4, 0 };
    
    void         nextUserMode(int dir);
    
    /// the current mode (this decide which action is done with the mouse)
    UserMode     userMode = MOUSE_ROTATE;

    
    View         savedView("backup");
    UserMode     mouseAction = MOUSE_TRANSLATE;  ///< the action being performed by the mouse
    Vector3      mouseDown;                 ///< position where mouse button was pressed down
    GLint        mouseX, mouseY;            ///< current position of mouse in pixels
    Vector3      depthAxis;                 ///< vector normal to desired rotation
    Vector3      mouseAxis;                 ///< axis of rotation for MOUSE_SPIN
    real         mouseZoomScalar;           ///< accessory scalar for zooming
    Vector3      ROI[2], ROIdown[2];        ///< Regions of interest selected with the mouse
    int          savedWindowPos[4] = { 512, 512, 10, 10 };
    
    real         nearZ  = 0;       ///< normalized device Z-coordinate-Z of front-plane
    real         midZ   = 0.5;     ///< normalized device Z-coordinate-Z of middle
    real         farZ   = 1.0;     ///< normalized device Z-coordinate-Z of back-plane

    unsigned int imageIndex = 0;   ///< index for image name
    
    /// function pointer for shift-click actions
    void (*mouseClickCallback)(int, int, const Vector3 &, int) = 0;
    
    /// function pointer for shift-motion actions
    void (*mouseDragCallback)(int, int, Vector3 &, const Vector3 &, int) = 0;
    
    /// function pointer for shift-click actions
    void (*normalKeyCallback)(unsigned char, int, int) = processNormalKey;
    
    /// function pointer for shift-motion actions
    void (*specialKeyCallback)(int, int, int) = processSpecialKey;
    
    // --------------- GRAPHICS

    unsigned int flash_count = 0;
    std::string  flash_text;
    
    void         dummyDisplay();
    void         (*displayCallback)() = dummyDisplay;
    
    void         setROI(Vector3);
    void         setROI(Vector3, Vector3);
    bool         insideROI(Vector3);    
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 initialize views[0]
 */
void glApp::initialize()
{
    for ( int ii = 0; ii < 16; ++ii )
        fKey[ii] = false;
    
    /// start with one single view
    views.clear();
    View tmp("view_base");
    views.push_back(tmp);
}


/**
 This will disable OpenGL depth-test for DIM<3
 */
void glApp::setDimensionality(const int d)
{
    if ( mDIM != d )
    {
        //flashText("dimensionality changed to %i", d);
        mDIM = d;
        userMode = ( d == 3 ) ? MOUSE_ROTATE : MOUSE_TRANSLATE;
    }
    
    if ( 0 == views.size() )
        initialize();
}


void glApp::displayFunc(void (*func)())
{
    displayCallback = func;
}


/**
 The first argument is a display callback.
 The second argument is the dimensionality (2 or 3) that appears to the user.
 
 glutInit() should be called before this function.
 */
void glApp::createWindow(void (*func)(), Glossary* glos)
{    
    displayCallback = func;

    if ( glos )
    {
        GP.read(*glos);
        views[0].read(*glos);
    }
    
    std::string name("cytosim - "+FilePath::get_dir());
    createWindow(name.c_str());
    
    //enter full-screen mode directly if asked:
    if ( GP.full_screen )
        enterFullScreen();
}

//------------------------------------------------------------------------------
void glApp::enterFullScreen()
{
    GP.full_screen = true;
    //save the current window size and position:
    savedWindowPos[0] = glutGet(GLUT_WINDOW_WIDTH);
    savedWindowPos[1] = glutGet(GLUT_WINDOW_HEIGHT);
    savedWindowPos[2] = glutGet(GLUT_WINDOW_X);
    savedWindowPos[3] = glutGet(GLUT_WINDOW_Y);
    //invoke full screen from GLUT
    glutFullScreen();
    //std::clog << "Fullscreen window " << glutGetWindow() << std::endl;
}


void glApp::exitFullScreen()
{
    GP.full_screen = false;
    //restore saved window dimensions:
    glutReshapeWindow(savedWindowPos[0], savedWindowPos[1]);
    glutPositionWindow(savedWindowPos[2], savedWindowPos[3]);
}


void glApp::switchFullScreen()
{
    if ( GP.full_screen )
        exitFullScreen();
    else
        enterFullScreen();
}

#if !defined(GLUT_WINDOW_SCALE)
#    define GLUT_WINDOW_SCALE 199
#endif

/**
 Adjust the size of window to maximize the vertical or horizontal dimension,
 without changing the aspect ratio of the window.
 */
void glApp::maximizeDisplay()
{
    int s = glutGet(GLUT_WINDOW_SCALE);
    if ( s < 0 ) s = 1;
 
    const int maxWidth  = s * glutGet(GLUT_SCREEN_WIDTH);
    const int maxHeight = s * ( glutGet(GLUT_SCREEN_HEIGHT) - 49 );

    const real zx = (real)maxWidth  / glutGet(GLUT_WINDOW_WIDTH);
    const real zy = (real)maxHeight / glutGet(GLUT_WINDOW_HEIGHT);
    
    glutPositionWindow(0, 45);
    if ( zx < zy )
        glutReshapeWindow(maxWidth, zx*glutGet(GLUT_WINDOW_HEIGHT));
    else
        glutReshapeWindow(zy*glutGet(GLUT_WINDOW_WIDTH), maxHeight);
}

//------------------------------------------------------------------------------
#pragma mark -

int glApp::createWindow(const char * window_name)
{
    View & view = glApp::currentView();
    
    std::ostringstream oss;
    oss << "rgba";
    if ( GP.buffered )      oss << " double";
    if ( view.depth_test )  oss << " depth";
    if ( view.stencil )     oss << " stencil~" << view.stencil;
    if ( view.multisample ) oss << " samples~" << view.multisample;
    if ( view.retina )      oss << " hidpi";

    //std::clog << "GLUT string mode " << oss.str() << std::endl;
    
    // set GLUT display mode:
    std::string str = oss.str();
    glutInitDisplayString(str.c_str());

    // set window size:
    glutInitWindowSize(view.width(), view.height());
    
    // set window position:
    view.window_position[0] += 20;
    view.window_position[1] += 20;
    glutInitWindowPosition(view.window_position[0], view.window_position[1]);
    
    int win = 0;
    if ( window_name )
        win = glutCreateWindow(window_name);
    else
        win = glutCreateWindow("cytosim");
    
    assert_true( win > 0 );
    //std::clog << "new window " << win << std::endl;

    //duplicate current View:
    if ( win >= (int)views.size() )
        views.resize(win+1, view);
    
    views[win].window(win);
    views[win].initGL();

    glutReshapeFunc(resizeWindow);
    glutKeyboardFunc(normalKeyCallback);
    glutSpecialFunc(specialKeyCallback);
    glutMouseFunc(processMouseClick);
    glutMotionFunc(processMouseDrag);
    glutPassiveMotionFunc(processPassiveMouseMotion);
    attachMenu(GLUT_RIGHT_BUTTON);

    if ( win <= 1 )
        glutDisplayFunc(displayMain);
    else
        glutDisplayFunc(displayPlain);
    
    //register call back for traveling:
    travelingTimer(win);
    return win;
}


/**
 This will not destroy the main window
 */
void glApp::destroyWindow(int win)
{
    if ( win == mZoomWindow )
        mZoomWindow = -1;
    
    if ( 1 < win  &&  win < (int)views.size()  &&  views[win].window() > 0 )
    {
        //std::clog << "Destroy window " << win << std::endl;
        assert_true( views[win].window() == win );
        glutDestroyWindow(views[win].window());
        views[win].window(0);
    }
}


void glApp::resizeWindow(int w, int h)
{
    unsigned win = glutGetWindow();
    views[win].reshape(w, h);
    flashText("Resize window %i %i", w, h);
    glClear(GL_COLOR_BUFFER_BIT);
    glutPostRedisplay();
}


void glApp::setScale(real s)
{
    if ( 0 == views.size() )
        initialize();

    views[0].setScale(s);
    
    if ( views.size() > 1 )
    {
        int win = glutGetWindow();
        // update all window-associated views:
        for ( unsigned n = 1; n < views.size(); ++n )
        {
            View & view = views[n];
            if ( view.window() > 0 )
            {
                glutSetWindow(view.window());
                view.setScale(s);
                view.set();
            }
        }
        glutSetWindow(win);
    }
}



/**
 Export pixels currently stored in graphical memory, into a file called `name`,
 in the particular image format specified in `format`.
 This does not require access to the simulation world.
 */
int glApp::savePNG(char const* name)
{
    return SaveImage::saveImage(name, "png");
}


#include "offscreen.h"

int glApp::savePNG(char const* name, unsigned mag, unsigned downsample)
{
    View view = currentView();
    view.reshape(view.window_size[0]*mag, view.window_size[1]*mag);
    
    if ( OffScreen::createBuffer(view.width(), view.height(), 0) )
    {
        view.initGL();
        displayPlain();
        int res = SaveImage::saveImage(name, "png", downsample);
        OffScreen::releaseBuffer();
        return res;
    }
    else
    {
        std::cerr << "Failed to create off-screen pixels" << std::endl;
        return 1;
    }
}


//------------------------------------------------------------------------------
#pragma mark -


/// this works even if no window is open
View& glApp::currentView()
{
    assert_true( views.size() > 0 );
    
    if ( views.size() <= 1 )
        return views[0];
    else
        return views[glutGetWindow()];
}


void glApp::resetView()
{
    assert_true( glutGetWindow() < (int)views.size() );
    glApp::currentView().reset();
}


void glApp::resetAllViews()
{
    for ( unsigned int n = 0; n < views.size(); ++n )
        views[n].reset();
}

//------------------------------------------------------------------------------
void glApp::travelingTimer(const int win)
{
    View & view = views[win];

    if ( view.traveling )
    {
        assert_true( view.window() > 0 );
        glutTimerFunc(view.traveling, travelingTimer, win);
        glutSetWindow(view.window());
        view.travelingMotion(0.001*view.traveling);
        glutPostRedisplay();
    }
    else
        glutTimerFunc(100, travelingTimer, win);
}

//------------------------------------------------------------------------------
#pragma mark -

/** Only 2D */
bool glApp::insideROI(Vector3 pos)
{
    bool inX = ( ROI[0].XX < pos.XX  &&  pos.XX < ROI[1].XX );
    bool inY = ( ROI[0].YY < pos.YY  &&  pos.YY < ROI[1].YY );
    return inX && inY;
}

/** Only 2D */
void glApp::setROI(Vector3 a)
{
    ROI[0] = a;
    ROI[1] = a;
}

/** Only 2D */
void glApp::setROI(Vector3 a, Vector3 b)
{
    ROI[0].XX = ( a.XX < b.XX ? a.XX : b.XX );
    ROI[1].XX = ( a.XX < b.XX ? b.XX : a.XX );
    ROI[0].YY = ( a.YY < b.YY ? a.YY : b.YY );
    ROI[1].YY = ( a.YY < b.YY ? b.YY : a.YY );
    ROI[0].ZZ = ( a.ZZ < b.ZZ ? a.ZZ : b.ZZ );
    ROI[1].ZZ = ( a.ZZ < b.ZZ ? b.ZZ : a.ZZ );
}


void glApp::matchROI(int win)
{
    //std::clog << " zoom on ROI " << win << " wind = " << glutGetWindow() << std::endl;
    if ( 0 < win && win < (int)views.size() && views[win].window() > 0 )
    {
        int w = glutGetWindow();
        glutSetWindow(views[win].window());
        views[win].matchROI(ROI[0], ROI[1]);
        glutPostRedisplay();
        glutSetWindow(w);
    }
}

//------------------------------------------------------------------------------
//------------------------------ keyboard commands -----------------------------
//------------------------------------------------------------------------------
#pragma mark -

void glApp::help(std::ostream & os)
{
    os << "                           Mouse Controls\n";
    os << "\n";
    os << " The display can be changed with mouse click-and-drag, and the outcome,\n";
    os << " depends on the current 'mode' which is changed by pressing the TAB key:\n";
    os << "         Rotate                                          (3D only)\n";
    os << "         Translate in XY\n";
    os << "         Action                             (click to grab fibers)\n";
    os << "         Translate in XZ                                 (3D only)\n";
    os << "         Spin                               (Rotation in XY plane)\n";
    os << "         Zoom  (Click in central part of window, and drag outward)\n";
    os << "         Select/move region-of-interest\n";
    os << "\n";
    os << "  In the default mode, you can SHIFT-CLICK to grab the filaments!\n";
    os << "  A Menu is accessed by a right click\n";
    os << "  Optionally, you might be able to zoom in/out with the mouse-wheel\n";
    os << "\n";
    os << "                           Keyboard Controls\n\n";
    os << " + -         Zoom in and out; hold SHIFT for fine motion\n";
    os << " arrow keys  Translate in XY; hold SHIFT for fine motion\n";
    os << " z           Reset view and refresh display\n";
    os << " h           Hide/show help\n";
    os << " b           Show/hide a 10 um scale bar\n";
    os << " f           toggle fullscreen mode\n";
    os << " y           save PNG image\n";
}


//------------------------------------------------------------------------------
void glApp::nextUserMode(int dir)
{
    unsigned u = userMode;
    do {
        u = ( u + dir + MOUSE_PASSIVE ) % MOUSE_PASSIVE;
    } while ( actionValidity[u] > mDIM );

    userMode = (UserMode)u;
    switch ( u )
    {
        case MOUSE_ROTATE:    flashText("Mouse: Rotate");       break;
        case MOUSE_TRANSLATE: flashText("Mouse: Translate");    break;
        case MOUSE_ACTIVE:    flashText("Mouse: Active");       break;
        case MOUSE_TRANSLATEZ:flashText("Mouse: Translate-Z");  break;
        case MOUSE_SPIN:      flashText("Mouse: Spin");         break;
        case MOUSE_ZOOM:      flashText("Mouse: Zoom");         break;
        case MOUSE_MAGNIFIER: flashText("Mouse: Magnifier");    break;
        case MOUSE_SET_ROI:   flashText("Mouse: Select ROI");   break;
        case MOUSE_MOVE_ROI:  flashText("Mouse: Move ROI");     break;
        case MOUSE_SELECT:    flashText("Mouse: Select");       break;
        case MOUSE_PASSIVE:   flashText("Mouse: Passive");      break;
    }
}


///\todo flexible key assignment map to accomodate different keyboard layouts
void glApp::processNormalKey(unsigned char c, int, int)
{
    View & view = glApp::currentView();
    
    /* In the switch below:
     - use 'break' if the display needs a refresh
     - use 'return' if redrawing is not necessary.
    */
    switch (c)
    {
        case 17:
            if ( glutGetModifiers() &  GLUT_ACTIVE_CTRL )
                exit(EXIT_SUCCESS);
            break;
        
            
        case 9:          // ascii 9 is the horizontal tab
        case 25:         // ascii 25 is shift-tab
            nextUserMode(c==9 ? 1 : -1);
            break;
        
        
        case 27:             // ascii 27 is the escape key
            if ( GP.full_screen )
                exitFullScreen();
            else
                destroyWindow(glutGetWindow());
            break;
        
            
        case 'f':
            switchFullScreen();
            break;
        
        case 'F':
            maximizeDisplay();
            break;

        case 'z':
            view.reset();
            postRedisplay();
            break;
        
        
        case 'n':
            view.slice = ( view.slice + 1 ) % 4;
            flashText("view:slice = %i", view.slice);
            break;
        
        
        case 'N':
            if ( userMode == MOUSE_SET_ROI )
            {
                mZoomWindow = createWindow("view");
                matchROI(mZoomWindow);
            }
            else
            {
                createWindow();
            }
            break;
        
        
        case 'b':
            view.scale_bar_mode = ( view.scale_bar_mode + 1 ) % 4;
            break;
        
            
        case 'h':
            GP.show_message = ( GP.show_message + 1 ) % 3;
            if ( GP.show_message == 2 ) {
                std::ostringstream oss;
                help(oss);
                GP.message = oss.str();
            }
            else
                GP.message = "Please, visit www.cytosim.org";
            break;
        
        
        case 'x':
            view.show_axes = ( view.show_axes ? 0 : mDIM );
            break;
            
            
        case 'y': {
            char name[1024] = { 0 };
            snprintf(name, sizeof(name), "image%04i.png", imageIndex++);
            savePNG(name);
        } break;
        
        case 'Y': {
            char name[1024] = { 0 };
            snprintf(name, sizeof(name), "image%04i.png", imageIndex++);
            savePNG(name, 4, 2);
        } break;
            
        //------------------------- zooming & window resizing:
            
        case '_':
            view.zoom_out(1.0352649238);
            break;
        case '-':
            view.zoom_out(1.4142135623);
            break;
        
        
        case '+':
            view.zoom_in(1.0352649238);
            break;
        case '=':
            view.zoom_in(1.4142135623);
            break;
        
        
        default:
            flashText("glapp ignored key %i [%c]", c, c);
            return;
    }
    
    //if break was used, redisplay is needed:
    postRedisplay();
    buildMenu();
}


void glApp::normalKeyFunc(void (*func)(unsigned char, int, int))
{
    normalKeyCallback = func;
}

//------------------------------------------------------------------------------
void glApp::buildfKeyString()
{
    //rebuild the function key string:
    fKey[0] = false;
    assert_true( sizeof(fKeyString) >= 4*12 );
    for ( int ii = 0; ii < 12; ++ii ) {
        if ( fKey[ii+1] ) {
            fKeyString[3*ii  ] = 'F';
            fKeyString[3*ii+1] = (char)('1'+ii);
            fKey[0] = true;
        }
        else {
            fKeyString[3*ii  ] = ' ';
            fKeyString[3*ii+1] = ' ';
        }
        fKeyString[3*ii+2] = ' ';
    }
    
    if ( fKey[0] )
        fKeyString[3*9+2] = '\0';
    else
        fKeyString[0] = '\0';
}

void glApp::bindFunctionKeys(bool * f)
{
    fKey = f;
    buildfKeyString();
}

bool glApp::functionKey(int k)
{
    if ( 0 < k && k < 17 )
        return fKey[k];
    else
        return false;
}

void glApp::toggleFunctionKey(int key)
{
    switch ( key )
    {
        case GLUT_KEY_F1:    fKey[1]  = !fKey[1];   break;
        case GLUT_KEY_F2:    fKey[2]  = !fKey[2];   break;
        case GLUT_KEY_F3:    fKey[3]  = !fKey[3];   break;
        case GLUT_KEY_F4:    fKey[4]  = !fKey[4];   break;
        case GLUT_KEY_F5:    fKey[5]  = !fKey[5];   break;
        case GLUT_KEY_F6:    fKey[6]  = !fKey[6];   break;
        case GLUT_KEY_F7:    fKey[7]  = !fKey[7];   break;
        case GLUT_KEY_F8:    fKey[8]  = !fKey[8];   break;
        case GLUT_KEY_F9:    fKey[9]  = !fKey[9];   break;
        case GLUT_KEY_F10:   fKey[10] = !fKey[10];  break;
        case GLUT_KEY_F11:   fKey[11] = !fKey[11];  break;
        case GLUT_KEY_F12:   fKey[12] = !fKey[12];  break;
        
        default:
            //MSG("ignored input key %i mouse=( %i %i )\n", key, x, y);
            return;
    }
    
    buildfKeyString();
    postRedisplay();
}

/**
 arrow-keys controls translation, and
 arrow-keys with 'ALT' pressed controls rotation.
 
 motion is reduced by holding down SHIFT.
 */
void glApp::processSpecialKey(int key, int, int)
{
    View & view = glApp::currentView();

    Vector3 vec, dxy(0, 0, 0);
    
    real F = ( glutGetModifiers() & GLUT_ACTIVE_SHIFT ) ? 0.0625 : 1;

    switch ( key )
    {
        case GLUT_KEY_HOME:      view.reset();            glutPostRedisplay(); return;
        case GLUT_KEY_PAGE_UP:   view.zoom_in(1.4142);    glutPostRedisplay(); return;
        case GLUT_KEY_PAGE_DOWN: view.zoom_out(1.4142);   glutPostRedisplay(); return;
        case GLUT_KEY_LEFT:      dxy.set(-F,0,0);         break;
        case GLUT_KEY_RIGHT:     dxy.set(+F,0,0);         break;
        case GLUT_KEY_DOWN:      dxy.set(0,-F,0);         break;
        case GLUT_KEY_UP:        dxy.set(0,+F,0);         break;
        default:                 toggleFunctionKey(key);  return;
    }

    // inverse the rotation of the current view:
    Quaternion<real> rot = view.rotation.conjugated();
    
    if ( glutGetModifiers() & GLUT_ACTIVE_ALT )
    {
        // Rotate view
        rot.rotateVector(vec, cross(Vector3(0, 0, 1), dxy));
        rot.setFromAxis(vec, F * (M_PI/8));
        view.rotate_by(rot);
    }
    else
    {
        // Translate view
        if ( userMode == MOUSE_TRANSLATEZ )
            dxy.set(dxy.XX, 0, dxy.YY);
        rot.rotateVector(vec, dxy);
        //std::clog << "vec " << vec << "\n";
        view.move_by((128*view.pixelSize())*vec);
    }
    glutPostRedisplay();
}


void glApp::specialKeyFunc(void (*func)(int, int, int))
{
    specialKeyCallback = func;
}

//------------------------------------------------------------------------------
#pragma mark -

int buildFogMenu()
{
    static int menu = 0;
    if ( menu == 0 )
    {
        menu = glutCreateMenu(glApp::processMenuEvent);
        glutAddMenuEntry("Disable",          100);
        glutAddMenuEntry("Linear ",          101);
        glutAddMenuEntry("Exponential 1/16", 102);
        glutAddMenuEntry("Exponential 1/8",  103);
        glutAddMenuEntry("Exponential 1/4",  104);
        glutAddMenuEntry("Exponential 1/2",  105);
        glutAddMenuEntry("Exponential 1",    106);
        glutAddMenuEntry("Exponential 2",    107);
        glutAddMenuEntry("Exponential 4",    108);
        glutAddMenuEntry("Exponential 8",    109);
        glutAddMenuEntry("Exponential 16",   110);
    }
    return menu;
}

int buildWindowSizeMenu()
{
    static int menu = 0;
    if ( menu == 0 )
    {
        menu = glutCreateMenu(glApp::processMenuEvent);
        glutAddMenuEntry("256x256",   200);
        glutAddMenuEntry("384x384",   201);
        glutAddMenuEntry("512x256",   202);
        glutAddMenuEntry("512x384",   203);
        glutAddMenuEntry("512x512",   204);
        glutAddMenuEntry("768x768",   205);
        glutAddMenuEntry("1024x512",  206);
        glutAddMenuEntry("1024x768",  207);
        glutAddMenuEntry("1024x1024", 208);
        glutAddMenuEntry("1536x1536", 209);
        glutAddMenuEntry("426x240 (240p)",    210);
        glutAddMenuEntry("640x360 (360p)",    211);
        glutAddMenuEntry("854x480 (480p)",    212);
        glutAddMenuEntry("1280x720 (720p)",   213);
        glutAddMenuEntry("1920x1080 (1080p)", 214);
        glutAddMenuEntry("2560x1440 (1440p)", 215);
    }
    return menu;
}


int buildClipMenu()
{
    static int menu = 0;
    if ( menu == 0 )
    {
        menu = glutCreateMenu(glApp::processMenuEvent);
        glutAddMenuEntry("Disable",    300);
        
        glutAddMenuEntry(" X > 0",     301);
        glutAddMenuEntry(" X < 0",     302);
        glutAddMenuEntry("-1 < X < 1", 303);
        
        glutAddMenuEntry(" Y > 0",     311);
        glutAddMenuEntry(" Y < 0",     312);
        glutAddMenuEntry("-1 < Y < 1", 313);
        
        glutAddMenuEntry(" 0 < Z",     321);
        glutAddMenuEntry(" Z < 0",     322);
        glutAddMenuEntry("-1 < Z < 1", 323);
        glutAddMenuEntry("-0.5 < Z < 0.5", 324);
    }
    return menu;
}


int glApp::buildMenu()
{
    static int menu = 0;
    static int menu1, menu2, menu3;
    
    //std::clog << "buildMenu" << std::endl;
    if ( menu )
        clearMenu(menu);
    else {
        menu1 = buildFogMenu();
        menu2 = buildWindowSizeMenu();
        menu3 = buildClipMenu();
        menu  = glutCreateMenu(processMenuEvent);
    }
    
    glutAddSubMenu("Fog",            menu1);
    glutAddSubMenu("Window Size",    menu2);
    glutAddSubMenu("Slice",          menu3);
    glutAddMenuEntry("Reset View",         1);
    glutAddMenuEntry("Show/hide Scalebar", 2);
    glutAddMenuEntry("Show/hide XYZ-axes", 3);
    glutAddMenuEntry(GP.full_screen?"Exit Fullscreen":"Enter Fullscreen", 4);
    glutAddMenuEntry(mDIM==2?"Use 3D Controls":"Use 2D Controls", 7);
    glutAddMenuEntry("New Window",   10);
    glutAddMenuEntry("Close Window", 11);
    glutAddMenuEntry("Quit",         20);
    
    return menu;
}

//------------------------------------------------------------------------------

void glApp::clearMenu(int menu)
{
    glutSetMenu(menu);
    const int mx = glutGet(GLUT_MENU_NUM_ITEMS);
    for ( int m = mx; m > 0; --m )
        glutRemoveMenuItem(m);
    assert_true( glutGet(GLUT_MENU_NUM_ITEMS) == 0 );
}

void glApp::attachMenu(int b)
{
    buildMenu();
    assert_true( b==GLUT_LEFT_BUTTON || b==GLUT_MIDDLE_BUTTON || b==GLUT_RIGHT_BUTTON );
    glutAttachMenu(b);
}

void glApp::processMenuEvent(int item)
{
    View & view = glApp::currentView();
    switch( item )
    {
        case 0:   return;
        case 1:   view.reset();                      break;
        case 2:   view.scale_bar_mode = ! view.scale_bar_mode;    break;
        case 3:   view.show_axes = ( view.show_axes ? 0 : mDIM ); break;
        case 4:   switchFullScreen();                break;
        case 7:   setDimensionality(mDIM==2?3:2);    break;

        case 10:  createWindow();                    break;
        case 11:  destroyWindow(glutGetWindow());    break;
        
        case 20:  exit(EXIT_SUCCESS);                break;
            
        case 100: view.enableFog(0, 0);              break;
        case 101: view.enableFog(1, 0);              break;
        case 102: view.enableFog(2, 0.0625);         break;
        case 103: view.enableFog(2, 0.125);          break;
        case 104: view.enableFog(2, 0.25);           break;
        case 105: view.enableFog(2, 0.5);            break;
        case 106: view.enableFog(2, 1);              break;
        case 107: view.enableFog(2, 2);              break;
        case 108: view.enableFog(2, 4);              break;
        case 109: view.enableFog(2, 8);              break;
        case 110: view.enableFog(2, 16);             break;
            
        case 200: glutReshapeWindow(256, 256);       break;
        case 201: glutReshapeWindow(384, 384);       break;
        case 202: glutReshapeWindow(512, 256);       break;
        case 203: glutReshapeWindow(512, 384);       break;
        case 204: glutReshapeWindow(512, 512);       break;
        case 205: glutReshapeWindow(768, 768);       break;
        case 206: glutReshapeWindow(1024, 512);      break;
        case 207: glutReshapeWindow(1024, 768);      break;
        case 208: glutReshapeWindow(1024, 1024);     break;
        case 209: glutReshapeWindow(1536, 1536);     break;
        case 210: glutReshapeWindow(426, 240);       break;
        case 211: glutReshapeWindow(640, 360);       break;
        case 212: glutReshapeWindow(854, 480);       break;
        case 213: glutReshapeWindow(1280, 720);      break;
        case 214: glutReshapeWindow(1920, 1080);     break;
        case 215: glutReshapeWindow(2560, 1440);     break;
        
        case 300:
            view.disableClipPlane(0);
            view.disableClipPlane(1);
            break;
            
        case 301:
            view.enableClipPlane(0, Vector3(+1,0,0), 0);
            view.disableClipPlane(1);
            break;
            
        case 302:
            view.enableClipPlane(0, Vector3(-1,0,0), 0);
            view.disableClipPlane(1);
            break;
            
        case 303:
            view.enableClipPlane(0, Vector3(+1,0,0), 1);
            view.enableClipPlane(1, Vector3(-1,0,0), 1);
            break;
 
        case 311:
            view.enableClipPlane(0, Vector3(0,+1,0), 0);
            view.disableClipPlane(1);
            break;
            
        case 312:
            view.enableClipPlane(0, Vector3(0,-1,0), 0);
            view.disableClipPlane(1);
            break;
            
        case 313:
            view.enableClipPlane(0, Vector3(0,+1,0), 1);
            view.enableClipPlane(1, Vector3(0,-1,0), 1);
            break;
 
        case 321:
            view.enableClipPlane(0, Vector3(0,0,+1), 0);
            view.disableClipPlane(1);
            break;
            
        case 322:
            view.enableClipPlane(0, Vector3(0,0,-1), 0);
            view.disableClipPlane(1);
            break;
            
        case 323:
            view.enableClipPlane(0, Vector3(0,0,+1), 1);
            view.enableClipPlane(1, Vector3(0,0,-1), 1);
            break;

        case 324:
            view.enableClipPlane(0, Vector3(0,0,+1), 0.5);
            view.enableClipPlane(1, Vector3(0,0,-1), 0.5);
            break;

        default: ABORT_NOW("unknown menu item");
    }
    glutPostRedisplay();
    buildMenu();
}

//------------------------------------------------------------------------------
//--------------------------------  MOUSE  -------------------------------------
//------------------------------------------------------------------------------
#pragma mark -

void  glApp::actionFunc(void (*func)(int, int, const Vector3 &, int))
{
    mouseClickCallback = func;
}

void  glApp::actionFunc(void (*func)(int, int, Vector3 &, const Vector3 &, int))
{
    mouseDragCallback = func;
}

#if !defined(GLUT_WHEEL_UP)
#  define GLUT_WHEEL_UP    3
#  define GLUT_WHEEL_DOWN  4
#endif

//------------------------------------------------------------------------------
void glApp::processMouseClick(int button, int state, int mX, int mY)
{
    View & view = glApp::currentView();
    int winW = view.width(), winH = view.height();    

    //printf("mouse button %i (%4i %4i) state %i key %i\n", button, mx, my, state, specialKeys);

    mouseX = mX;
    mouseY = winH-mY;
    
    savedView = view;
    savedView.getMatrices();
    mouseDown = savedView.unproject(mouseX, mouseY, nearZ);

    if ( state == GLUT_UP )
    {
        if ( mouseAction == MOUSE_SET_ROI )
            matchROI(mZoomWindow);
        
        /*
         Zooming with the mouse-wheel requires an extended GLUT.
         http://iihm.imag.fr/blanch/howtos/MacOSXGLUTMouseWheel.html
         
         The zoom preserves the position pointed at by the mouse.
         */
        real wz = 1.0;

        if ( button == GLUT_WHEEL_UP )
            wz = 0.96969696;
        if ( button == GLUT_WHEEL_DOWN )
            wz = 1.031250;

        if ( wz != 1 )
        {
            /* 
             in 2D, we do not allow any shift in Z,
             and in 3d, we zoom in on the middle Z-plane
             */
            if ( mDIM == 3 )
                mouseDown = savedView.unproject(mouseX, mouseY, midZ);
            else
                mouseDown.ZZ = 0;
            view.zoom_out(wz);
            view.move_to((1.0-wz)*mouseDown+wz*view.focus);
        }

        glutSetCursor(GLUT_CURSOR_INHERIT);
        mouseAction = MOUSE_PASSIVE;
        glutPostRedisplay();
        return;
    }

    glutSetCursor(GLUT_CURSOR_CROSSHAIR);
    
    // action is primarily decided by current mode
    mouseAction = userMode;
    specialKeys = glutGetModifiers();

    // change the mouse action if the CONTROL is pressed:
    if ( specialKeys & GLUT_ACTIVE_CTRL )
    {
        switch ( mouseAction )
        {
            case MOUSE_TRANSLATE: mouseAction = (mDIM==2)?MOUSE_SPIN:MOUSE_ROTATE; break;
            case MOUSE_SPIN:      mouseAction = (mDIM==2)?MOUSE_TRANSLATE:MOUSE_TRANSLATEZ; break;
            case MOUSE_ZOOM:      mouseAction = MOUSE_TRANSLATE; break;
            case MOUSE_SET_ROI:   mouseAction = MOUSE_TRANSLATE; break;
            case MOUSE_ROTATE:    mouseAction = MOUSE_TRANSLATE; break;
            case MOUSE_TRANSLATEZ:mouseAction = (mDIM==2)?MOUSE_TRANSLATE:MOUSE_ROTATE;  break;
            default: break;
        }
    }
    
    // change the mouse action because the shift key is down:
    if ( specialKeys & GLUT_ACTIVE_SHIFT )
    {
        mouseAction = MOUSE_ACTIVE;
        specialKeys -= GLUT_ACTIVE_SHIFT;
    }
    
    switch( mouseAction )
    {
        case MOUSE_TRANSLATE:
        {
        } return;
            
            
        case MOUSE_TRANSLATEZ:
        {
            Vector3 uc = savedView.unproject(winW/2.0, winH/2.0, nearZ);
            depthAxis  = ( uc - savedView.focus ).normalized();
            Vector3 ud = savedView.unproject(winW/2.0, winH, nearZ);
            mouseAxis  = ( ud - uc ).normalized();
        } break;
            
            
        case MOUSE_ROTATE:
        {
            /* 
            Choose the amplification factor for mouse controlled rotation:
            for a value of one, the rotation exactly follows the mouse pointer 
            */
            const real amplification = 3.0;
            depthAxis  = mouseDown - savedView.focus;
            depthAxis *= amplification / depthAxis.normSqr();
        } break;
            
            
        case MOUSE_SPIN:
        {
            Vector3 up = mouseDown;
            mouseDown  = savedView.unproject(winW/2.0, winH/2.0, nearZ);
            mouseAxis  = ( mouseDown - savedView.focus ).normalized();
            depthAxis  = up - mouseDown;
        } break;
            
            
        case MOUSE_ZOOM:
        {
            real xx = mX - 0.5 * savedView.width();
            real yy = mY - 0.5 * savedView.height();
            mouseZoomScalar = sqrt( xx*xx + yy*yy );
            if ( mouseZoomScalar > 5 )
                mouseZoomScalar = 1.0 / mouseZoomScalar;
            else
                mouseZoomScalar = 0.2;
        } break;
        
        case MOUSE_MAGNIFIER:
        {
            mouseDown = savedView.unproject(mouseX, mouseY, midZ);
            if ( mDIM == 2 ) mouseDown.ZZ = 0;
            glutSetCursor(GLUT_CURSOR_NONE);
        } break;
            
        case MOUSE_SET_ROI:
        case MOUSE_MOVE_ROI:
        {
            mouseDown = savedView.unproject(mouseX, mouseY, midZ);
            if ( insideROI(mouseDown) )
            {
                ROIdown[0] = ROI[0];
                ROIdown[1] = ROI[1];
                mouseAction = MOUSE_MOVE_ROI;
            }
            if ( mouseAction == MOUSE_SET_ROI )
            {
                setROI(mouseDown);
#if ( DIM == 2 )
                flashText("Position = %.4f %.4f", ROI[0].XX, ROI[0].YY);
#else
                flashText("Position = %.4f %.4f %.4f", ROI[0].XX, ROI[0].YY, ROI[0].ZZ);
#endif
            }
        } break;
        
            
        case MOUSE_ACTIVE:
        {
            if ( mouseClickCallback )
            {
                mouseDown = savedView.unproject(mouseX, mouseY, midZ);
                //std::clog << "Action down at "<<mouseDown<<std::endl;
                mouseClickCallback(mouseX, mouseY, mouseDown, specialKeys);
            }
        }
        
        case MOUSE_SELECT:
        case MOUSE_PASSIVE:
            return;
    }
    glutPostRedisplay();
}



//------------------------------------------------------------------------------
void glApp::processMouseDrag(int mX, int mY)
{
    //printf("mouse motion (%i %i) %i\n", mx, my, mouseAction);
    View & view = glApp::currentView();
    int winH = view.height();    

    mouseX = mX;
    mouseY = winH-mY;

    Vector3 mouse = savedView.unproject(mouseX, mouseY, nearZ);

    switch( mouseAction )
    {
        case MOUSE_ROTATE:
        {
            /* we should rotate after: Q <- dQ * sQ, however dQ is defined in the 
            reference frame rotated by sQ already, so dQ = sQ * dM * inv(sQ).
            This leads to the multiplication on the right: Q <- sQ * dM. */
            Quaternion<real> q;
            q.setFromAxis( cross(depthAxis, mouse-mouseDown) );
            view.rotate_to( savedView.rotation * q );
        } break;
        
        
        case MOUSE_SPIN:
        {
            real cos = depthAxis * ( mouse - mouseDown );
            real sin = cross( depthAxis, ( mouse - mouseDown ) ) * mouseAxis;
            Quaternion<real> q;
            q.setFromAxis( mouseAxis, atan2( sin, cos ) );
            view.rotate_to( savedView.rotation * q );
        } break;
        
        
        case MOUSE_TRANSLATE:
        {
            view.move_to( savedView.focus - ( mouse - mouseDown ) );
        } break;
        
        
        case MOUSE_TRANSLATEZ:
        {
            real S = ( mouse - mouseDown ) * mouseAxis;
            Vector3 move = mouse - mouseDown - S * ( depthAxis + mouseAxis );
            view.move_to( savedView.focus - move );
        } break;
        
        
        case MOUSE_ZOOM:
        {
            real xx = mX - 0.5 * savedView.width();
            real yy = mY - 0.5 * savedView.height();
            real Z = mouseZoomScalar * sqrt( xx*xx + yy*yy );
            if ( Z > 0.001 ) view.zoom_to( savedView.zoom * Z );
        } break;
        
            
        case MOUSE_MAGNIFIER:
        {
            mouseDown = savedView.unproject(mouseX, mouseY, midZ);
            if ( mDIM == 2 ) mouseDown.ZZ = 0;
        } break;

        case MOUSE_SET_ROI:
        {
            setROI(mouseDown, savedView.unproject(mouseX, mouseY, midZ));
            Vector3 d = ROI[1] - ROI[0];
            flashText("ROI %.3fx%.3f diag. %.3f", d.XX, d.YY, d.norm());
        } break;
        
        
        case MOUSE_MOVE_ROI:
        {
            Vector3 d = savedView.unproject(mouseX, mouseY, midZ) - mouseDown;
            if ( glutGetWindow() == mZoomWindow )
                d = -d;
            ROI[0] = ROIdown[0] + d;
            ROI[1] = ROIdown[1] + d;
            flashText("ROI moved %.3f %.3f", d.XX, d.YY);
            matchROI(mZoomWindow);
        } break;
            
        
        case MOUSE_ACTIVE:
        {
            if ( mouseDragCallback )
            {
                mouse = savedView.unproject(mouseX, mouseY, midZ);
                //std::clog << "Action move at " << mouse << std::endl;
                mouseDragCallback(mouseX, mouseY, mouseDown, mouse, specialKeys);
            }
        } break;
        
        
        case MOUSE_SELECT:
        case MOUSE_PASSIVE: 
            break;
    }
    glutPostRedisplay();
}

//------------------------------------------------------------------------------
void glApp::processPassiveMouseMotion(int mx, int my)
{
    //printf("passive mouse (%i %i)\n", mx, my);
    //int x = glutGet(GLUT_WINDOW_WIDTH)-8;
    //int y = glutGet(GLUT_WINDOW_HEIGHT)-8;
}

//------------------------------------------------------------------------------
#pragma mark -

void glApp::flashText(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    char tmp[1024];
    vsnprintf(tmp, sizeof(tmp), fmt, args);
    va_end(args);

    flash_text = tmp;
    if ( flash_count == 0 )
        glutTimerFunc(100, flashTimer, 1);
    flash_count = 50;
    
    if ( views.size() > 1  &&  views[1].window()==1 )
        glutPostWindowRedisplay(1);
    else
        ;//std::clog << tmp << std::endl;
}

void glApp::flashTimer(int win)
{
    //printf("flash %i, cnt %i at %lu\n", flash_text, cnt, clock());
    
    if ( --flash_count > 0 )
    {
        //continue count-down:
        glutTimerFunc(100, flashTimer, win);
    }
    else
    {
        //count-down completed:
        flash_text = "";
        //restore the front buffer image:
        glutPostRedisplay();
    }
}


//------------------------------------------------------------------------------
//draw the Region of interest
void glApp::drawROI(Vector3 roi[2])
{
    glPushAttrib(GL_ENABLE_BIT|GL_COLOR_BUFFER_BIT);
    glEnable(GL_LINE_STIPPLE);
    glEnable(GL_COLOR_LOGIC_OP);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_LIGHTING);
    glLogicOp(GL_XOR);
    glLineStipple(1, 0x0f0f);
    glBegin(GL_LINE_LOOP);
    glVertex3d(roi[0].XX, roi[0].YY, roi[0].ZZ);
    glVertex3d(roi[1].XX, roi[0].YY, roi[0].ZZ);
    glVertex3d(roi[1].XX, roi[1].YY, roi[0].ZZ);
    glVertex3d(roi[0].XX, roi[1].YY, roi[0].ZZ);
    glEnd();
    if ( mDIM == 3 ) {
        glBegin(GL_LINE_LOOP);
        glVertex3d(roi[0].XX, roi[0].YY, roi[1].ZZ);
        glVertex3d(roi[1].XX, roi[0].YY, roi[1].ZZ);
        glVertex3d(roi[1].XX, roi[1].YY, roi[1].ZZ);
        glVertex3d(roi[0].XX, roi[1].YY, roi[1].ZZ);
        glEnd();
    }
    glPopAttrib();
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Display additional over-the-window text for user information
 */
void glApp::displayInfo(const int W, const int H)
{
#if ( 0 )
    // display frames-per-second
    static char buf[16] = "\0";
    static int nb_frames = 0;
    static long msec = TicToc::milli_seconds_today();
    ++nb_frames;
    long now = TicToc::milli_seconds_today();
    if ( now > msec + 1000 )
    {
        double fps = nb_frames * 1000.0 / ( now - msec );
        snprintf(buf, sizeof(buf), "%3.2f fps", fps);
        msec = now;
        nb_frames = 0;
    }
    glColor3f(1,1,1);
    gleDisplayText(buf, GLUT_BITMAP_8_BY_13, 0x0, 3, W, H);
#endif
    
    if ( fKeyString[0] != '\0' )
    {
        glColor3f(1,0,1);
        gleDisplayText(fKeyString, GLUT_BITMAP_HELVETICA_18, 0x0, 1, W, H);
    }
    
    if ( flash_text.size() )
    {
        glColor3f(0.6,0.6,1);
        gleDisplayText(flash_text.c_str(), GLUT_BITMAP_9_BY_15, 0x0, 2, W, H);
    }
    
    if ( GP.show_message && GP.message.size() )
    {
        glColor3f(1,1,1);
        gleDisplayText(GP.message.c_str(), GLUT_BITMAP_8_BY_13, 0x000000CC, 4, W, H);
    }
    
    if ( mouseAction == MOUSE_SET_ROI  ||  mouseAction == MOUSE_MOVE_ROI )
    {
        //draw the Region of interest
        glLineWidth(0.5);
        glColor3f(1,1,0);
        drawROI(ROI);
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void glApp::dummyDisplay()
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(0.0, 0.0, 1.0);
    gleDrawText("Maximum number of windows exceeded", GLUT_BITMAP_8_BY_13);

    if ( GP.buffered )
        glutSwapBuffers();
    else
        glFlush();        
}


/**
 This is used for any secondary window.
 It does not show the interactive feedback to user.
 */
void glApp::displayPlain()
{
    View & view = glApp::currentView();
    view.openDisplay();
    displayCallback();
    view.closeDisplay();
    
    if ( GP.buffered )
        glutSwapBuffers();
    else
        glFlush();
}


void print_cap(const char * name, GLenum cap)
{
    GLint b[10] = { 0 };
    glGetIntegerv(cap, b);
    printf("gl %12s = %i\n", name, b[0]);
}


void glApp::displayMagnifier(real Z)
{
    View view = glApp::currentView();
    
    int W = view.width();
    int H = view.height();
    
    int M = ( W > H ? H : W ) / 2;
    
#if ( 0 )
    GLint readbuf = 0, drawbuf = 0;
    glGetIntegerv(GL_READ_BUFFER, &readbuf);
    glGetIntegerv(GL_DRAW_BUFFER, &drawbuf);
    printf("normal buffers: read %i draw %i\n", readbuf, drawbuf);
#endif
    
    if ( OffScreen::createBuffer(M, M, 0) )
    {
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        view.initGL();
        view.setScale(M*view.pixelSize()/Z);
        view.reshape(M, M);
        view.move_to(mouseDown);
        view.zoom_to(1);
        
        view.openDisplay();
        displayCallback();
        view.endClipPlanes();
#if ( 0 )
        GLubyte * tmp = new GLubyte[3*M*M];
        glReadPixels(0, 0, M, M, GL_RGB, GL_UNSIGNED_BYTE, tmp);
        SaveImage::savePixels("pixels.png", "png", tmp, W, H);
        delete[] tmp;
#endif
        //print_cap("read", GL_READ_BUFFER);
        //print_cap("draw", GL_DRAW_BUFFER);

        // restore writing destination:
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        //glDrawBuffer(drawbuf);

        glDisable(GL_ALPHA_TEST);
        glBlitFramebuffer(0, 0, M, M,
                          mouseX-M/2, mouseY-M/2,
                          mouseX+M/2, mouseY+M/2,
                          GL_COLOR_BUFFER_BIT,GL_LINEAR);
        //checkError("glBlitFramebuffer()");

        OffScreen::releaseBuffer();
        glViewport(0, 0, W, H);
        glPopMatrix();
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        //glReadBuffer(readbuf);
    }
}


/**
 This is used for the main window
 */
void glApp::displayMain()
{
    View & view = glApp::currentView();
    view.set();

    view.openDisplay();
    displayCallback();
    view.closeDisplay();
    
    if ( mouseAction == MOUSE_MAGNIFIER )
        displayMagnifier(5);

    //display over-the-window features:
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    view.displayMessages();
    displayInfo(view.width(), view.height());
    
    glPopAttrib();

    if ( GP.buffered )
        glutSwapBuffers();
    else
        glFlush();
}

/**
 call glutPostRedisplay() if only the current window needs to be updated
 */
void glApp::postRedisplay()
{
    for ( unsigned n = 1; n < views.size(); ++n )
        if ( views[n].window() > 0 )
            glutPostWindowRedisplay(n);
}


void glApp::setMessage(std::string const& msg)
{
    glApp::currentView().setMessage(msg);
}


void glApp::setMessageT(std::string const& msg)
{
    glApp::currentView().setMessageT(msg);
}



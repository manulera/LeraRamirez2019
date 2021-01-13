// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "assert_macro.h"
#include "view.h"
#include "gle.h"
#include "glu_unproject.cc"

using namespace gle;

//------------------------------------------------------------------------------

View::View(const std::string& n) : ViewProp(n)
{    
    mWindowId = 0;
    
    visRegion[0] = view_size;
    visRegion[1] = view_size;
    visRegion[2] = view_size;
    
    hasMatrices = false;
    setPixelSize();
}

View::~View()
{
}


//------------------------------------------------------------------------------
#pragma mark -

/// GL_DEPTH_CLAMP_NV should be defined in OpenGL/glext.h
#define GL_DEPTH_CLAMP 0x864F

void View::initGL()
{
    glDisable(GL_STENCIL_TEST);
    glDisable(GL_ALPHA_TEST);
    glDisable(GL_DITHER);
    
    if ( 1 )
    {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
    else
        glDisable(GL_BLEND);
    
    if ( multisample > 1 )
    {
        glEnable(GL_MULTISAMPLE);
        /*
         GLint s = 0;
         glGetIntegerv(GL_MAX_SAMPLES, &s);
         std::clog << "OpenGL samples = " << samples << "  max = " << s << std::endl;
         */
    }
    else
    {
        glDisable(GL_MULTISAMPLE);
        glEnable(GL_POINT_SMOOTH);
        glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        /*
         Do not enable POLYGON_SMOOTH, which destroys joints of triangulated surfaces
         glEnable(GL_POLYGON_SMOOTH);
         glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
         */
    }
    
    if ( depth_clamp )
        glEnable(GL_DEPTH_CLAMP);
    else
        glDisable(GL_DEPTH_CLAMP);

    if ( depth_test )
    {
        glEnable(GL_DEPTH_TEST);
        //glDepthFunc(GL_LESS);
        glDepthFunc(GL_LEQUAL);
        // enable Alpha Test to discard transparent pixels:
        glEnable(GL_ALPHA_TEST);
        glAlphaFunc(GL_GREATER, 0.05);
    }
    else
    {
        //std::clog << "no depth-test" << std::endl;
        glDisable(GL_DEPTH_TEST);
    }

    set();
}


void View::openDisplay() const
{
    back_color.load_clear();
    inner_color.load_back();
    setFog(fog_type, fog_param, fog_color);
    setLights();
    setClipPlanes();
    if ( slice )
        sliceView(slice);
    gle::gleReportErrors(stderr, "in glApp::startDisplay()");
}

/**
 Displays Axes and Scalebar
 */
void View::closeDisplay() const
{
    endClipPlanes();

    if ( show_axes )
        gleDrawAxes(axes_size, show_axes);

    if ( scale_bar_mode )
    {
        glPushAttrib(GL_ENABLE_BIT);
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        scale_bar_color.load();
        drawScaleBar(scale_bar_mode, scale_bar_size);
        glPopAttrib();
    }
}

/**
 Display top and bottom info messages
 */
void View::displayMessages() const
{    
    if ( message_top.size() )
    {
        glEnable(GL_COLOR_LOGIC_OP);
        glLogicOp(GL_INVERT);
        gleDisplayText(message_top.c_str(), 0, 0x0, 3, width(), height());
        glDisable(GL_COLOR_LOGIC_OP);
    }

    if ( message_bot.size() )
    {
        std::string msg = message_bot;
        if ( msg.compare(0, 4, "none") )
        {
            // in non-interactive mode, only print the first line:
            if ( mWindowId != 1 )
                msg = msg.substr(0, msg.find('\n'));
            front_color.load();
            gleDisplayText(msg.c_str(), 0, 0x0, 0, width(), height());
        }
    }
}


/**
 Set two light sources
 */
void View::setLights(bool local) const
{
    glMatrixMode(GL_MODELVIEW);
    if ( local )
    {
        glPushMatrix();
        glLoadIdentity();
    }
    
    glShadeModel(GL_SMOOTH);
    
    GLfloat matWhite[]  = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat matGray[]   = { 0.2, 0.2, 0.2, 1.0 };
    GLfloat matBlack[]  = { 0.0, 0.0, 0.0, 1.0 };
    //GLfloat matBlue[]   = { 0.0, 0.0, 1.0, 1.0 };
    
    glMaterialfv(GL_FRONT, GL_AMBIENT,   matBlack);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   matBlack);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  matWhite);
    glMateriali (GL_FRONT, GL_SHININESS, 64);

    // set a gray color for the back-side of everything
    glMaterialfv(GL_BACK, GL_AMBIENT,  matGray);
    glMaterialfv(GL_BACK, GL_DIFFUSE,  matBlack);
    glMaterialfv(GL_BACK, GL_SPECULAR, matBlack);
    glMateriali (GL_BACK, GL_SHININESS, 32);
    
    GLfloat lightDiffuse[]  = { 0.8, 0.8, 0.8, 1.0 };
    GLfloat lightSpecular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat lModelAmbient[] = { 0.4, 0.4, 0.4, 1.0 };
    
    GLfloat light0Pos[] = { 5.0,-3.0, 3.0, 0.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, light0Pos);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  lightDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
    glEnable(GL_LIGHT0);
    
    GLfloat light1Pos[] = {-4.0, 0.0,-3.0, 0.0 };
    glLightfv(GL_LIGHT1, GL_POSITION, light1Pos);
    glLightfv(GL_LIGHT1, GL_DIFFUSE,  lightDiffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpecular);
    glEnable(GL_LIGHT1);
    
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lModelAmbient);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
    
    // let GL normalize the normals:
    glEnable(GL_NORMALIZE);

    if ( local )
        glPopMatrix();
}


//------------------------------------------------------------------------------
#pragma mark -

void View::makeProjection()
{
    //calculate the visible region in the 3 directions:
    if ( window_size[0] > window_size[1] )
    {
        visRegion[0] = view_size;
        visRegion[1] = view_size * window_size[1] / real(window_size[0]);
    }
    else
    {
        visRegion[0] = view_size * window_size[0] / real(window_size[1]);
        visRegion[1] = view_size;
    }
    visRegion[2] = view_size;
    
    //std::clog << "View::makeProjection        " << visRegion[0] << " " << visRegion[1] << " " << visRegion[2] << "\n";
    
    if ( perspective == 2 )
    {
        // this creates a stronger perspective, if eye_shift is also changed:
        //    real eye_shift[3] = { 0, 0, -2.0*view_size };
        glFrustum(-0.5*visRegion[0], 0.5*visRegion[0],
                  -0.5*visRegion[1], 0.5*visRegion[1],
                   1.0*visRegion[2],   6*visRegion[2]);
    }
    else if ( perspective )
    {
        // this creates a perspective, if eye_shift is also changed:
        //    real eye_shift[3] = { 0, 0, -2.0*view_size };
        glFrustum(-0.5*visRegion[0], 0.5*visRegion[0],
                  -0.5*visRegion[1], 0.5*visRegion[1],
                   1.0*visRegion[2],  11*visRegion[2]);
    }
    else
    {
        // The back-plane is set far back to avoid any clipping there
        //    real eye_shift[3] = { 0, 0, -0.5*view_size };
        glOrtho(-0.5*visRegion[0], 0.5*visRegion[0],
                -0.5*visRegion[1], 0.5*visRegion[1],
                 0,                    visRegion[2]);
    }
}


void View::set()
{
    //std::clog << "setProjection win " << window() << std::endl;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    makeProjection();
    glMatrixMode(GL_MODELVIEW);
    setModelView();
    setPixelSize();
}


void View::setModelView() const
{
    //std::cerr<<"setModelView win " << window() << std::endl;

    real eye_shift[3] = { 0, 0, (perspective?-2.0:-0.5)*view_size };

    // setup the OpenGL transformation matrix
    glMatrixMode(GL_MODELVIEW);
#if ( 1 )
    real mat[16];
    rotation.setOpenGLMatrix(mat, eye_shift);
    gleLoadMatrix(mat);
#else
    glLoadIdentity();
    real rot[3];
    real ang = 180 / M_PI * rotation.getAngle(rot);
    gleTranslate(eye_shift[0], eye_shift[1], eye_shift[2]);
    gleRotate(ang, rot[0], rot[1], rot[2]);
#endif
    gleScale(zoom);
    
    // point-of-focus:
    gleTranslate(-(focus+focus_shift));
    
#if ( 0 )
    std::clog << "View::setModelView eye      " << eye_shift[2] << "\n";
    std::clog << "View::setModelView rotation " << rotation << "\n";
    std::clog << "View::setModelView zoom     " << zoom << "\n";
    std::clog << "View::setModelView focus    " << focus+focus_shift << "\n";
#endif
}


/**
 This will change what is visible in the Z direction near and far
 from the observer, using clipping planes and fog.
 A a function of `mode`:
 - 0 : disabled
 - 1 : show ( Z > 0 ) with fog
 - 2 : show ( Z < 0 ) with fog
 - 3 : show a thin slice ( -a < Z < a ) where a = 10% of visRegion[2]
 .
 
 This uses GL_CLIP_PLANE2 and GL_CLIP_PLANE3
 */
void View::sliceView(int mode) const
{
    real off = -view_size * 0.5;
    
    if ( mode == 1 )
    {
        setClipPlaneEye(GL_CLIP_PLANE3, Vector3(0,0,+1), -off);
        if ( !depth_clamp )
            setFog(1, 0, fog_color);
    }
    else if ( mode == 2 )
    {
        setClipPlaneEye(GL_CLIP_PLANE3, Vector3(0,0,-1), +off);
        setFog(1, 1, fog_color);
    }
    else if ( mode == 3 )
    {
        real thk = 0.1 * visRegion[2];
        setClipPlaneEye(GL_CLIP_PLANE2, Vector3(0,0,-1), thk+off);
        setClipPlaneEye(GL_CLIP_PLANE3, Vector3(0,0,+1), thk-off);
    }
}


void View::setPixelSize()
{
    mPixelSize = visRegion[0] / ( zoom * window_size[0] );
    //std::clog << "View::pixelSize       " << mPixelSize << "\n";
}


void View::reshape(int W, int H)
{
    //std::clog << "View::reshaped " << W << " " << H << std::endl;
    window_size[0] = W;
    window_size[1] = H;
    glViewport(0, 0, W, H);
    set();
}

//------------------------------------------------------------------------------
#pragma mark -

void View::reset()
{
    zoom = 1;
    auto_scale = 1;
    focus.zero();
    focus_shift.zero();
    rotation.set(1,0,0,0);
    setModelView();
    setPixelSize();
}


void View::zoom_to(real z)
{
    //std::clog << "zoom_to " << z << " " << this << std::endl;
    zoom = z;
    setModelView();
    setPixelSize();
}


void View::setScale(real s)
{
    //std::clog << "View::setScale("<<s<<") win=" << window() << std::endl;
    view_size = s;
}


void View::matchROI(Vector3 a, Vector3 b)
{
    focus = 0.5 * ( a + b );
    real r = 0.5 * ( a - b ).norm_inf();
    
    // require 3 pixels to zoom in:
    if ( r > 3 * mPixelSize )
        zoom = view_size / r;
    
    setModelView();
    setPixelSize();
}


void View::move_to(const Vector3& d)
{
    focus = d;
    setModelView();
}


void View::move_shift(const Vector3& d)
{
    focus_shift = d;
    setModelView();
}


void View::rotate_to(const Quaternion<real> & q)
{
    rotation = q.normalized();
    setModelView();
}


/**
 This assumes that vector `dir` is normalized
 */
void View::align_with(const Vector3 & dir)
{
    // axis is obtained by vector product: axis = (1, 0, 0) ^ a
    Vector3 axis( 0, dir.ZZ, -dir.YY );
    // cosine is scalar product, sine is norm of vector-product:
    real cs = dir.XX, si = axis.norm();
    if ( si > REAL_EPSILON )
    {
        rotation.setFromAxis(axis, atan2(si, cs));
        setModelView();
    }
}


void View::travelingMotion(real dt)
{
    focus += dt * auto_translation;
    zoom = zoom * ( 1 + dt * auto_zoom );
    Quaternion<real> Q = auto_rotation.scaledAngle(dt) * rotation;
    rotation = Q.normalized();
    setModelView();
}


//------------------------------------------------------------------------------
#pragma mark -

void View::getMatrices()
{
    //get the transformation matrices, to be used for mouse control
    glGetIntegerv(GL_VIEWPORT,         mViewport);
    glGetDoublev(GL_PROJECTION_MATRIX, mProjection);
    glGetDoublev(GL_MODELVIEW_MATRIX,  mModelview);
    hasMatrices = true;
}

/**
 This set a matrix like glOrtho()
 */
void View::setProjectionOrtho(GLdouble * mat)
{
    real L =-0.5*visRegion[0];
    real R = 0.5*visRegion[0];
    real B =-0.5*visRegion[1];
    real T = 0.5*visRegion[1];
    real N = 0;
    real F = visRegion[2];

    for ( int i = 0; i < 16; ++i )
        mat[i] = 0;
    
    mat[ 0] =  2.0 / ( R - L );
    mat[ 5] =  2.0 / ( T - B );
    mat[10] = -2.0 / ( F - N );
    
    mat[12] = -( R + L ) / ( R - L );
    mat[13] = -( T + B ) / ( T - B );
    mat[14] = -( F + N ) / ( F - N );
    mat[15] = 1.0;
}


void View::dump() const
{
    GLint vp[4] = { 0 };
    glGetIntegerv(GL_VIEWPORT, vp);
    std::clog << "viewport = " << vp[0] << " " << vp[1] << " " << vp[2] << " " << vp[3] << std::endl;
}

/**
 Transforms the given window coordinates into user coordinates.
 
 It uses the matrices obtained at the last call of getMatrices(),
 or the current matrices if get_matrices == true.
 
 For more info, try `man gluUnProject`
 */
Vector3 View::unproject(GLdouble x, GLdouble y, GLdouble z, bool get_matrices)
{
    GLdouble ux = 0, uy = 0, uz = 0;
    if ( get_matrices )
    {
        GLint      vp[4];
        GLdouble   mv[16];
        GLdouble   pj[16];
        
        glGetIntegerv(GL_VIEWPORT,         vp);
        glGetDoublev(GL_PROJECTION_MATRIX, pj);
        glGetDoublev(GL_MODELVIEW_MATRIX,  mv);

        setProjectionOrtho(pj);
        myUnproject(x, y, z, mv, pj, vp, &ux, &uy, &uz);
    }
    else if ( hasMatrices )
        myUnproject(x, y, z, mModelview, mProjection, mViewport, &ux, &uy, &uz);
    else
        std::cerr << "warning: View::unproject called without matrices\n";
    
    //printf("unproject( %.2f, %.2f, %.2f ) = ( %.2f, %.2f, %.2f )\n", x, y, z, ux, uy, uz);
    return Vector3(ux, uy, uz);
}


//------------------------------------------------------------------------------
#pragma mark -

void View::setFog(int type, real param, gle_color color) const
{
    GLint gl_type = 0;
    switch( type )
    {
        case 1: gl_type = GL_LINEAR; break;
        case 2: gl_type = GL_EXP;    break;
        case 3: gl_type = GL_EXP2;   break;
        default: glDisable(GL_FOG); return;
    }
   
    glEnable(GL_FOG);
    glFogi(GL_FOG_MODE, gl_type);
    
    if ( gl_type == GL_LINEAR )
    {
        glFogf(GL_FOG_START, param*visRegion[2]);
        glFogf(GL_FOG_END, (param*2+1)*visRegion[2]);
    }
    else
    {
        glFogf(GL_FOG_DENSITY, param/visRegion[2]);
    }
    
    glFogfv(GL_FOG_COLOR, color.data());
}

void View::enableFog(const int type, const real param, gle_color color)
{
    fog_type = type;
    fog_param = param;
    fog_color = color;
}

void View::enableFog(const int type, const real param)
{
    fog_type = type;
    fog_param = param;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 The plane equations is relative to the model
 */
void View::setClipPlane(GLenum glp, Vector3 dir, real sca) const
{
    GLdouble eq[] = {dir.XX, dir.YY, dir.ZZ, sca};
    glClipPlane(glp, eq);
    glEnable(glp);
}

/**
 The plane equation is relative to the camera
 */
void View::setClipPlaneEye(GLenum glp, Vector3 dir, real sca) const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    setClipPlane(glp, dir, sca);
    glPopMatrix();
}


void View::setClipPlanes() const
{
    for ( unsigned int ix = 0; ix < NB_CLIP_PLANES; ++ix )
    {
        if ( clip_plane_mode[ix] == 1 )
            setClipPlane(GL_CLIP_PLANE0+ix, clip_plane_vector[ix], clip_plane_scalar[ix]);
        else if ( clip_plane_mode[ix] == 2 )
            setClipPlaneEye(GL_CLIP_PLANE0+ix, clip_plane_vector[ix], clip_plane_scalar[ix]);
        else
            glDisable(GL_CLIP_PLANE0+ix);
    }
}

/*
 Disable all clip planes, including the one set by sliceView()
 */
void View::endClipPlanes() const
{
    for ( unsigned int ix = 0; ix < NB_CLIP_PLANES; ++ix )
        glDisable(GL_CLIP_PLANE0+ix);
}

void View::enableClipPlane(unsigned ix, Vector3 dir, real sc, bool model)
{
    if ( ix < NB_CLIP_PLANES )
    {
        clip_plane_mode[ix]   = ( model ? 1 : 2 );
        clip_plane_vector[ix] = dir;
        clip_plane_scalar[ix] = sc;
    }
}

void View::disableClipPlane(unsigned int ix)
{
    if ( ix < NB_CLIP_PLANES )
    {
        clip_plane_mode[ix] = 0;
        glDisable(GL_CLIP_PLANE0+ix);
    }
}

bool View::hasClipPlane(unsigned int ix) const
{
    if ( ix < NB_CLIP_PLANES )
        return clip_plane_mode[ix];
    return false;
}

//------------------------------------------------------------------------------
#pragma mark -


void View::drawScaleTicksH(GLfloat d, GLfloat a, GLfloat b) const
{
    glBegin(GL_LINES);
    for ( int ii = 1; ii < 10; ++ii )
    {
        glVertex2f( ii*d, a);
        glVertex2f( ii*d, b);
        glVertex2f(-ii*d, a);
        glVertex2f(-ii*d, b);
    }
    glEnd();
}


void View::drawScaleTicksV(GLfloat d, GLfloat a, GLfloat b) const
{
    glBegin(GL_LINES);
    for ( int ii = 1; ii < 10; ++ii )
    {
        glVertex2f( a, ii*d);
        glVertex2f( b, ii*d);
        glVertex2f( a,-ii*d);
        glVertex2f( b,-ii*d);
    }
    glEnd();
}


/**
 This will draw:
 - a horizontal box of length scale, bounded with Y=a and Y=b
 - lines every scale/10, of width (b-a)/5
 - lines every scale/100, of width (b-a)/25
 - lines every scale/1000, of width (b-a)/125
 .
 */
void View::drawScaleH(GLfloat s, GLfloat a, GLfloat b) const
{
    // draw a box of length 'scale'
    glLineWidth(3);
    glBegin(GL_LINE_LOOP);
    glVertex2f(-s/2, a);
    glVertex2f(-s/2, b);
    glVertex2f( s/2, b);
    glVertex2f( s/2, a);
    glEnd();

    // draw bars
    s /= 10;
    glBegin(GL_LINES);
    for ( int ii = -5; ii <= 5; ++ii )
    {
        glVertex2f(ii*s, a);
        glVertex2f(ii*s, b);
    }
    glEnd();
    
    // draw tick marks
    GLfloat w = 2;
    char str[8] = {0};
    do {
        s /= 10;
        a /= 10;
        b /= 10;
        if ( s > 4 * pixelSize() )
        {
            glLineWidth(w);
            drawScaleTicksH(s, a, b);
            glRasterPos2f(s-6*pixelSize(), b-12*pixelSize());
            snprintf(str, sizeof(str), "%g", s);
            gleDrawText(str);
        }
        w /= 2;
    } while ( w >= 0.5 );
}


/**
 This will draw:
 - a vertical box of length scale, bounded at X=a and X=b
 - lines every scale/10, of width (b-a)/5
 - lines every scale/100, of width (b-a)/25
 - lines every scale/1000, of width (b-a)/125
 .
 */
void View::drawScaleV(GLfloat s, GLfloat a, GLfloat b) const
{
    // draw a box of length 'scale'
    glLineWidth(3);
    glBegin(GL_LINE_LOOP);
    glVertex2f(a, -s/2);
    glVertex2f(b, -s/2);
    glVertex2f(b,  s/2);
    glVertex2f(a,  s/2);
    glEnd();
    
    // draw bars
    s /= 10;
    glBegin(GL_LINES);
    for ( int ii = -5; ii <= 5; ++ii )
    {
        glVertex2f(a, ii*s);
        glVertex2f(b, ii*s);
    }
    glEnd();
    
    // draw tick marks
    GLfloat w = 2;
    char str[8] = {0};
    do {
        s /= 10;
        a /= 10;
        b /= 10;
        if ( s > 4 * pixelSize() )
        {
            glLineWidth(w);
            drawScaleTicksV(s, a, b);
            glRasterPos2f(b+pixelSize(), s-4*pixelSize());
            snprintf(str, sizeof(str), "%g", s);
            gleDrawText(str);
        }
        w /= 2;
    } while ( w >= 0.5 );
}


/**
 This will draw a centered cross with :
 - lines every scale/10, of width 1
 - lines every scale/100, of width 0.5
 - minuscule lines every scale/1000, of width 0.25
 .
 */
void View::drawScaleX(GLfloat scale) const
{
    GLfloat w = 4;
    GLfloat s =  scale;
    GLfloat a =  scale/20;
    GLfloat b = -scale/20;
    
    glLineWidth(w);
    glBegin(GL_LINES);
    glVertex2f(-scale, a);
    glVertex2f(-scale, b);
    glVertex2f( scale, a);
    glVertex2f( scale, b);
    glVertex2f(a, -scale);
    glVertex2f(b, -scale);
    glVertex2f(a,  scale);
    glVertex2f(b,  scale);
    glEnd();
    
    do {
        s /= 10;
        if ( s > 2 * pixelSize() )
        {
            glLineWidth(w);
            drawScaleTicksV(s, a, b);
            drawScaleTicksH(s, a, b);
        }
        a /= 10;
        b /= 10;
        w /= 2;
    } while ( w >= 0.5 );

    glLineWidth(0.5);
    glBegin(GL_LINES);
    glVertex2f(-scale, 0);
    glVertex2f( scale, 0);
    glVertex2f(0, -scale);
    glVertex2f(0,  scale);
    glEnd();
}


/**
 */
void View::drawScaleBar(int mode, const GLfloat scale) const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    GLfloat shift = 32 * pixelSize() * zoom;
    
    switch( mode )
    {
        case 0:
            break;
        case 1:
            gleTranslate(0, shift-0.5*visRegion[1], 0);
            gleScale(zoom);
            drawScaleH(scale, scale/10, 0);
            break;
        case 2:
            gleTranslate(0.5*visRegion[0]-shift, 0, 0);
            gleScale(zoom);
            drawScaleV(scale, -scale/10, 0);
            break;
        case 3: {
            gleScale(zoom);
            drawScaleX(scale);
        } break;
    }
    
    glPopMatrix();
}

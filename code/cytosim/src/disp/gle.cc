// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cctype>
#include "assert_macro.h"
#include "gle.h"
#include "glut.h"
#include "platonic.h"
#include "smath.h"

//#define GLE_USE_VBO


void gle::initialize()
{
    //std::clog << "gle::initialize()" << std::endl;
    /*
     GLfloat c[4*fin+1], s[4*fin+1];
     GLfloat da = M_PI_2 / (GLfloat) finesse;
     for ( int ii = 0; ii <= 4*finesse; ++ii )
     {
       GLfloat a = ii * da;
       co[ii] = cosf(a);
       si[ii] = sinf(a);
     }
     */
#ifdef GLE_USES_DISPLAY_LISTS
    initializeDL();
#endif
    
#ifdef GLE_USE_VBO
    initializeVBO();
#endif
}

void gle::release()
{
#ifdef GLE_USES_DISPLAY_LISTS
    releaseDL();
#endif

#ifdef GLE_USE_VBO
    releaseVBO();
#endif
}

//-----------------------------------------------------------------------
void gle::gleAlignX(const Vector2 & v)
{
    const GLfloat n = v.norm();
    //warning! this matrix is displayed transposed
    GLfloat mat[16] = {
        (GLfloat)v.XX, -(GLfloat)v.YY,  0,  0,
        (GLfloat)v.YY,  (GLfloat)v.XX,  0,  0,
        0,                          0,  n,  0,
        0,                          0,  0,  1 };
    glMultMatrixf(mat);
}


/**
Graphical elements are aligned in 3D along Z and this function is used
to rotate them in the XY plane for the 2D display.
The rotation is chosen such that the Y face of the rotated object points
down the Z axis. In this way, the lower part of the object is drawn first,
such that the upper half overwrites it and become the only visible part.
The display is thus correct even if DEPTH_TEST is disabled.
*/
void gle::gleAlignZ(const Vector2 & A, const Vector2 & B)
{    
    Vector2 D = B - A;
    GLfloat n = sqrt(D.XX*D.XX+D.YY*D.YY);
    if ( n < REAL_EPSILON ) return;
    //warning! this matrix is displayed transposed
    GLfloat mat[16] = {
        (GLfloat)D.YY/n, -(GLfloat)D.XX/n,   0,   0,
        0,                              0,  -1,   0,
        (GLfloat)D.XX,      (GLfloat)D.YY,   0,   0,
        (GLfloat)A.XX,      (GLfloat)A.YY,   0,   1 };
    glMultMatrixf(mat);
}


/**
 ts is the transverse scaling done in the XY plane after rotation
 */
void gle::gleAlignZ(const Vector2 & A, const Vector2 & B, real ts)
{
    Vector2 D = B - A;
    GLfloat p = ts / sqrt(D.XX*D.XX+D.YY*D.YY);
    //warning! this matrix is displayed transposed
    GLfloat mat[16] = {
        (GLfloat)D.YY*p, -(GLfloat)D.XX*p,            0,   0,
        0,                              0, -(GLfloat)ts,   0,
        (GLfloat)D.XX,      (GLfloat)D.YY,            0,   0,
        (GLfloat)A.XX,      (GLfloat)A.YY,            0,   1 };
    glMultMatrixf(mat);
}


void gle::gleRotate(const Vector3 & v1, const Vector3 & v2, const Vector3 & v3, bool inverse)
{
    GLfloat mat[16];
    for ( int ii = 0; ii < 3; ++ii )
    {
        if ( inverse )
        {
            mat[4*ii  ] = v1[ii];
            mat[4*ii+1] = v2[ii];
            mat[4*ii+2] = v3[ii];
        }
        else
        {
            mat[ii  ]   = v1[ii];
            mat[ii+4]   = v2[ii];
            mat[ii+8]   = v3[ii];
        }
        mat[ii+12]  = 0;
        mat[ii*4+3] = 0;
    }
    mat[15] = 1;
    glMultMatrixf(mat);
}


void gle::gleTransRotate(const Vector3 & v1, const Vector3 & v2,
                         const Vector3 & v3, const Vector3 & vt)
{
    //warning! this matrix is displayed here transposed
    GLfloat mat[16] = {
        (GLfloat)v1.XX, (GLfloat)v1.YY, (GLfloat)v1.ZZ, 0,
        (GLfloat)v2.XX, (GLfloat)v2.YY, (GLfloat)v2.ZZ, 0,
        (GLfloat)v3.XX, (GLfloat)v3.YY, (GLfloat)v3.ZZ, 0,
        (GLfloat)vt.XX, (GLfloat)vt.YY, (GLfloat)vt.ZZ, 1};
    glMultMatrixf(mat);
}


void gle::setClipPlane(GLenum glp, Vector1 const& dir, Vector1 const& pos)
{
    GLdouble eq[4] = { dir.XX, 0, 0, -dir*pos };
    glClipPlane(glp, eq);
}

void gle::setClipPlane(GLenum glp, Vector2 const& dir, Vector2 const& pos)
{
    GLdouble eq[4] = { dir.XX, dir.YY, 0, -dir*pos };
    glClipPlane(glp, eq);
}

void gle::setClipPlane(GLenum glp, Vector3 const& dir, Vector3 const& pos)
{
    GLdouble eq[4] = { dir.XX, dir.YY, dir.ZZ, -dir*pos };
    glClipPlane(glp, eq);
}


//-----------------------------------------------------------------------
#pragma mark - 2D Primitives


void gle::gleTriangle0()
{
    const GLfloat H = 0.5 * sqrt(3);
    glVertex2f(0, 1);
    glVertex2f(-H, -0.5);
    glVertex2f( H, -0.5);
}

void gle::gleTriangleS()
{
    glBegin(GL_TRIANGLES);
    glNormal3f(0, 0, 1);
    gleTriangle0();
    glEnd();
}

void gle::gleTriangleL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    gleTriangle0();
    glEnd();
}

//-----------------------------------------------------------------------

void gle::gleNabla0()
{
    const GLfloat H = 0.5 * sqrt(3);
    glVertex2f(0, -1);
    glVertex2f( H, 0.5);
    glVertex2f(-H, 0.5);
}

void gle::gleNablaS()
{
    glBegin(GL_TRIANGLES);
    glNormal3f(0, 0, 1);
    gleNabla0();
    glEnd();
}

void gle::gleNablaL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    gleNabla0();
    glEnd();
}

//-----------------------------------------------------------------------
void gle::gleSquare0()
{
    glVertex2f( 1,  1);
    glVertex2f(-1,  1);
    glVertex2f(-1, -1);
    glVertex2f( 1, -1);
}

void gle::gleSquareS()
{
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    gleSquare0();
    glEnd();
}

void gle::gleSquareL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    gleSquare0();
    glEnd();
}

//-----------------------------------------------------------------------
void gle::gleRectangle0()
{
    glVertex2f( 1,  0.5);
    glVertex2f(-1,  0.5);
    glVertex2f(-1, -0.5);
    glVertex2f( 1, -0.5);
}

void gle::gleRectangleS()
{
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    gleRectangle0();
    glEnd();
}

void gle::gleRectangleL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    gleRectangle0();
    glEnd();
}


//-----------------------------------------------------------------------
/// draw pentagon that has the same surface as the disc or Radius 1.
void gle::glePentagon0()
{
    const GLfloat R  = sqrt( 4 * M_PI / sqrt( 5 * ( 5 + 2 * sqrt(5))) );
    const GLfloat C1 = R * cosf(M_PI*0.1);
    const GLfloat S1 = R * sinf(M_PI*0.1);
    const GLfloat C3 = R * cosf(M_PI*0.3);
    const GLfloat S3 = R * sinf(M_PI*0.3);
    
    glVertex2f(  0,  1);
    glVertex2f(-C1,  S1);
    glVertex2f(-C3, -S3);
    glVertex2f( C3, -S3);
    glVertex2f( C1,  S1);
}

void gle::glePentagonS()
{
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f(0, 0);
    glePentagon0();
    glVertex2f(0, 1);
    glEnd();
}

void gle::glePentagonL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    glePentagon0();
    glEnd();
}

//-----------------------------------------------------------------------
/// draw hexagon that has the same surface as the disc or Radius 1.
void gle::gleHexagon0()
{
    const GLfloat R = sqrt( 2 * M_PI / ( 3 * sqrt(3) ));
    const GLfloat H = R * 0.5 * sqrt(3);
    const GLfloat X = R * 0.5;
    glVertex2f( R,  0);
    glVertex2f( X,  H);
    glVertex2f(-X,  H);
    glVertex2f(-R,  0);
    glVertex2f(-X, -H);
    glVertex2f( X, -H);
}

void gle::gleHexagonS()
{
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f(0, 0);
    gleHexagon0();
    glVertex2f(1, 0);
    glEnd();
}

void gle::gleHexagonL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    gleHexagon0();
    glEnd();
}

//-----------------------------------------------------------------------

void gle::gleCircle0(int n_seg)
{
    const GLfloat theta = 2 * M_PI / n_seg;
    const GLfloat c = cosf(theta);
    const GLfloat s = sinf(theta);

    GLfloat t;
    GLfloat x = 1;
    GLfloat y = 0;

    for( int ii = 0; ii < n_seg; ++ii )
    {
        glVertex2f(x, y);
            
        //apply the rotation matrix
        t = x;
        x = c * x - s * y;
        y = s * t + c * y;
    }
}

void gle::gleCircleL()
{
    glNormal3f(0, 0, 1);
    glBegin(GL_LINE_LOOP);
    gleCircle0(8*finesse);
    glVertex2f(1, 0);
    glEnd();
}

void gle::gleCircleS()
{
    glNormal3f(0, 0, 1);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(0, 0);
    gleCircle0(8*finesse);
    glVertex2f(1, 0);
    glEnd();
}

//-----------------------------------------------------------------------

void gle::gleStar0()
{
    const GLfloat R  = 1.2;
    const GLfloat C1 = R * cosf(M_PI*0.1);
    const GLfloat S1 = R * sinf(M_PI*0.1);
    const GLfloat C3 = R * cosf(M_PI*0.3);
    const GLfloat S3 = R * sinf(M_PI*0.3);
    const GLfloat H = -0.6;
    
    glVertex2f(    0,     R);
    glVertex2f( H*C3, -H*S3);
    glVertex2f(  -C1,    S1);
    glVertex2f( H*C1,  H*S1);
    glVertex2f(  -C3,   -S3);
    glVertex2f(    0,   H*R);
    glVertex2f(   C3,   -S3);
    glVertex2f(-H*C1,  H*S1);
    glVertex2f(   C1,    S1);
    glVertex2f(-H*C3, -H*S3);
}

void gle::gleStarS()
{
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f(0, 0);
    gleStar0();
    glVertex2f(0, 1);
    glEnd();
}

void gle::gleStarL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    gleStar0();
    glEnd();
}

//-----------------------------------------------------------------------

void gle::glePlusS()
{
    const GLfloat R = 1.1;
    const GLfloat C = 0.4;
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f( R,  C);
    glVertex2f(-R,  C);
    glVertex2f(-R, -C);
    glVertex2f( R, -C);
    glEnd();

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f( C,  R);
    glVertex2f(-C,  R);
    glVertex2f(-C,  C);
    glVertex2f( C,  C);
    glEnd();

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f( C, -C);
    glVertex2f(-C, -C);
    glVertex2f(-C, -R);
    glVertex2f( C, -R);
    glEnd();
}

void gle::glePlusL()
{
    const GLfloat R = 1.2;
    const GLfloat C = 0.6;
    
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    glVertex2f( C,  R);
    glVertex2f(-C,  R);
    glVertex2f(-C,  C);
    glVertex2f(-R,  C);
    glVertex2f(-R, -C);
    glVertex2f(-C, -C);
    glVertex2f(-C, -R);
    glVertex2f( C, -R);
    glVertex2f( C, -C);
    glVertex2f( R, -C);
    glVertex2f( R,  C);
    glVertex2f( C,  C);
    glEnd();
}


//-----------------------------------------------------------------------
#pragma mark - 3D Primitives

void gle::gleTube0(GLfloat a, GLfloat b, int fin)
{
    const GLfloat inc = M_PI / fin;
    const GLfloat max = 2 * M_PI + 0.5 * inc;

    glBegin(GL_TRIANGLE_STRIP);
    GLfloat ang = 0;
    while ( ang < max )
    {
        glNormal3f(cosf(ang), sinf(ang),  0);
        glVertex3f(cosf(ang), sinf(ang),  b);
        glVertex3f(cosf(ang), sinf(ang),  a);
        ang += inc;
    }
    glEnd();
}


void gle::gleTubeZ(GLfloat za, GLfloat ra, gle_color ca, GLfloat zb, GLfloat rb, gle_color cb)
{
    const GLfloat inc = M_PI / finesse;
    const GLfloat max = 2 * M_PI + 0.5 * inc;
    
    glPushAttrib(GL_CURRENT_BIT);
    glBegin(GL_TRIANGLE_STRIP);
    GLfloat ang = 0;
    while ( ang < max )
    {
        cb.load_load();
        glNormal3f(cosf(ang), sinf(ang), 0);
        glVertex3f(rb*cosf(ang), rb*sinf(ang), zb);
        
        ca.load_load();
        glNormal3f(cosf(ang), sinf(ang), 0);
        glVertex3f(ra*cosf(ang), ra*sinf(ang), za);
        ang += inc;
    }
    glEnd();
    glPopAttrib();
}

//-----------------------------------------------------------------------

void gle::gleTube1()
{
    gleTube0(0, 1, finesse/2);
}

void gle::gleLongTube1()
{
    gleTube0(-0.5, 1.5, finesse/2);
}

void gle::gleTube2()
{
    gleTube0(0, 1, finesse);
}

void gle::gleLongTube2()
{
    gleTube0(-0.5, 1.5, finesse);
}


void gle::gleHexTube1()
{
    gleTube0(0, 1, 3);
}

void gle::gleCylinder1()
{
    gleTube0(0, 1, finesse/2);
    gleTranslate(0,0,1);
    gleCircleS();
    gleTranslate(0,0,-1);
    glRotated(180,0,1,0);
    gleCircleS();
}

//-----------------------------------------------------------------------

#ifdef PLATONIC_H

//we use an icosahedron to approximate the sphere:
Platonic::Solid ico1(Platonic::Solid::ICOSAHEDRON, gle::finesse/4);
Platonic::Solid ico2(Platonic::Solid::ICOSAHEDRON, gle::finesse/2);
Platonic::Solid ico4(Platonic::Solid::ICOSAHEDRON, gle::finesse);
Platonic::Solid ico8(Platonic::Solid::ICOSAHEDRON, gle::finesse*2);

void gleSpherePlatonic(Platonic::Solid & ico)
{
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, ico.vertex_data());
    glNormalPointer(GL_FLOAT, 0, ico.vertex_data());
    glDrawElements(GL_TRIANGLES, 3*ico.nb_faces(), GL_UNSIGNED_INT, ico.faces_data());
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
}

void gle::gleSphere1() { gleSpherePlatonic(ico1); }
void gle::gleSphere2() { gleSpherePlatonic(ico2); }
void gle::gleSphere4() { gleSpherePlatonic(ico4); }
void gle::gleSphere8() { gleSpherePlatonic(ico8); }

#else

/// using trigonometric functions to draw a smooth ball
void gleSphereCS(unsigned fin)
{
    GLfloat c[2*fin+1], s[2*fin+1];
    for ( unsigned ii = 0; ii <= 2*fin; ++ii )
    {
        GLfloat ang = ii * M_PI / (GLfloat) fin;
        c[ii] = cosf(ang);
        s[ii] = sinf(ang);
    }
    
    for ( unsigned ii = 0; ii < fin; ++ii )
    {
        real uZ = c[ii  ], uR = s[ii  ];
        real lZ = c[ii+1], lR = s[ii+1];
        glBegin(GL_TRIANGLE_STRIP);
        for ( unsigned jj = 0; jj <= 2*fin; ++jj )
        {
            glNormal3f(c[jj]*uR, s[jj]*uR, uZ);
            glVertex3f(c[jj]*uR, s[jj]*uR, uZ);
            glNormal3f(c[jj]*lR, s[jj]*lR, lZ);
            glVertex3f(c[jj]*lR, s[jj]*lR, lZ);
        }
        glEnd();
    }
}

void gle::gleSphere1() { gleSphereCS(finesse); }
void gle::gleSphere2() { gleSphereCS(2*finesse); }
void gle::gleSphere4() { gleSphereCS(4*finesse); }
void gle::gleSphere8() { gleSphereCS(8*finesse); }

#endif


//-----------------------------------------------------------------------
#pragma mark -

/**
 Draw a cylindrical band on the equator of a sphere of radius 1.
 The band is in the XY plane. The axis of the cylinder is Z.
 The band is made of triangles indicating the clockwise direction.
 */
void gle::gleArrowedBand(const unsigned nb_triangles, GLfloat width)
{
    GLfloat a = 2 * M_PI / nb_triangles;
    GLfloat w = width * a / sqrt(3);
    GLfloat R = 1.0 / cos(a*0.5);
    
    glBegin(GL_TRIANGLES);
    glNormal3f(1, 0, 0);
    glVertex3f(1, 0, w);
    glVertex3f(1, 0,-w);
    for ( int ii = 1; ii < nb_triangles; ++ii )
    {
        GLfloat ang = ii * a;
        GLfloat c = R * cosf(ang);
        GLfloat s = R * sinf(ang);
        
        glNormal3f(c, s, 0);
        glVertex3f(c, s, 0);
        glVertex3f(c, s, w);
        glVertex3f(c, s,-w);
    }
    glNormal3f(1, 0, 0);
    glVertex3f(1, 0, 0);
    glEnd();
}


void gle::gleThreeBands(const unsigned nb_triangles)
{
    gleArrowedBand(nb_triangles, 0.25);
    glRotated(-90,1,0,0);
    gleArrowedBand(nb_triangles, 0.25);
    glRotated(90,0,1,0);
    gleArrowedBand(nb_triangles, 0.25);
}

//-----------------------------------------------------------------------
inline void icoFace(GLfloat* a, GLfloat* b, GLfloat* c)
{
    glNormal3f((a[0]+b[0]+c[0])/3.0, (a[1]+b[1]+c[1])/3.0, (a[2]+b[2]+c[2])/3.0);
    glVertex3fv(a);
    glVertex3fv(b);
    glVertex3fv(c);
}

void gle::gleIcosahedron1()
{
    const GLfloat tau=0.8506508084;      /* t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)  */
    const GLfloat one=0.5257311121;      /* one=1/sqrt(1+t^2) , unit sphere     */
    
    /* Twelve vertices of icosahedron on unit sphere */
    GLfloat pts[] = {
        +tau,  one,    0 , // 0
        -tau, -one,    0 , // 1
        -tau,  one,    0 , // 2
        +tau, -one,    0 , // 3
        +one,   0 ,  tau , // 4
        -one,   0 , -tau , // 5
        +one,   0 , -tau , // 6
        -one,   0 ,  tau , // 7
        0   ,  tau,  one , // 8
        0   , -tau, -one , // 9
        0   , -tau,  one , // 10
        0   ,  tau, -one };// 11
    
    /* The faces are ordered with increasing Z */
    glBegin(GL_TRIANGLES);
    icoFace(pts+3*5, pts+3*6,  pts+3*9);
    icoFace(pts+3*5, pts+3*11, pts+3*6);
    
    icoFace(pts+3*6, pts+3*3,  pts+3*9);
    icoFace(pts+3*2, pts+3*11, pts+3*5);
    icoFace(pts+3*1, pts+3*5,  pts+3*9);
    icoFace(pts+3*0, pts+3*6,  pts+3*11);
    
    icoFace(pts+3*0, pts+3*3,  pts+3*6);
    icoFace(pts+3*1, pts+3*2,  pts+3*5);

    icoFace(pts+3*1, pts+3*9,  pts+3*10);
    icoFace(pts+3*0, pts+3*11, pts+3*8);
    icoFace(pts+3*8, pts+3*11, pts+3*2);
    icoFace(pts+3*9, pts+3*3,  pts+3*10);

    icoFace(pts+3*0, pts+3*4,  pts+3*3);
    icoFace(pts+3*1, pts+3*7,  pts+3*2);
    
    icoFace(pts+3*0, pts+3*8,  pts+3*4);
    icoFace(pts+3*1, pts+3*10, pts+3*7);
    icoFace(pts+3*3, pts+3*4,  pts+3*10);
    icoFace(pts+3*7, pts+3*8,  pts+3*2);
     
    icoFace(pts+3*4, pts+3*8,  pts+3*7);
    icoFace(pts+3*4, pts+3*7,  pts+3*10);
    glEnd();
}

//-----------------------------------------------------------------------
void gle::gleCylinderZ()
{
    const GLfloat inc = M_PI / finesse;
    const GLfloat max = 2 * M_PI + 0.5 * inc;
    const GLfloat top = 0.5;
    const GLfloat bot = -0.5;
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f( 0, 0, -1 );
    glVertex3f( 0, 0, -1 );
    GLfloat ang = 0;
    while ( ang < max )
    {
        glVertex3f( cosf(ang), -sinf(ang), bot );
        ang += inc;
    }
    glEnd();

    glBegin(GL_TRIANGLE_STRIP);
    ang = 0;
    while ( ang < max )
    {
        glNormal3f( cosf(ang), sinf(ang), 0 );
        glVertex3f( cosf(ang), sinf(ang), top );
        glVertex3f( cosf(ang), sinf(ang), bot );
        ang += inc;
    }
    glEnd();
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f( 0, 0, 1 );
    glVertex3f( 0, 0, 0 );
    ang = 0;
    while ( ang < max )
    {
        glVertex3f( cosf(ang), sinf(ang), top );
        ang += inc;
    }
    glEnd();
}


void gle::gleCone1()
{
    const GLfloat inc = M_PI / finesse;
    const GLfloat max = 2 * M_PI + 0.5 * inc;
    
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f( 0, 0, -1 );
    glVertex3f( 0, 0, -1 );
    GLfloat ang = 0;
    while ( ang < max )
    {
        glVertex3f( cosf(ang), -sinf(ang), -1 );
        ang += inc;
    }
    glEnd();
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f( 0, 0, 1 );
    glVertex3f( 0, 0, 2 );
    GLfloat cn = 3.f/sqrt(10), sn = 1.f/sqrt(10);
    ang = 0;
    while ( ang < max )
    {
        glNormal3f( cn*cosf(ang), cn*sinf(ang), sn );
        glVertex3f( cosf(ang), sinf(ang), -1 );
        ang += inc;
    }
    glEnd();
}


void gle::gleCone1L()
{
    const GLfloat inc = M_PI / finesse;
    const GLfloat max = 2 * M_PI + 0.5 * inc;
    
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f( 0, 0, -1 );
    glVertex3f( 0, 0, -1 );
    GLfloat ang = 0;
    while ( ang < max )
    {
        glVertex3f( cosf(ang), -sinf(ang), -1 );
        ang += inc;
    }
    glEnd();
}


void gle::gleArrowTail1()
{
    const GLfloat inc = M_PI / finesse;
    const GLfloat max = 2 * M_PI + 0.5 * inc;
    
    
    GLfloat cn = M_SQRT1_2;

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f( 0, 0, -1 );
    glVertex3f( 0, 0, -0.5 );
    GLfloat ang = 0;
    while ( ang < max )
    {
        glNormal3f( -cn*cosf(ang), cn*sinf(ang), -cn );
        glVertex3f( cosf(ang), -sinf(ang), -1.5 );
        ang += inc;
    }
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    ang = 0;
    while ( ang < max )
    {
        glNormal3f( cosf(ang), sinf(ang),  0 );
        glVertex3f( cosf(ang), sinf(ang),  0.5 );
        glVertex3f( cosf(ang), sinf(ang), -1.5 );
        ang += inc;
    }
    glEnd();

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f( 0, 0, 1 );
    glVertex3f( 0, 0, 1.5 );
    ang = 0;
    while ( ang < max )
    {
        glNormal3f( cn*cosf(ang), cn*sinf(ang), cn );
        glVertex3f( cosf(ang), sinf(ang), 0.5 );
        ang += inc;
    }
    glEnd();
}

/**
 Draw three fins similar to the tail of a V2 rocket
 */
void gle::gleArrowTail2()
{
    GLfloat r = 0.1;  //bottom inner radius
    GLfloat c = 0.5, d = -0.5;
    GLfloat s = sqrt(3)/2, t = -s;
    GLfloat rc = r * c;
    GLfloat rs = r * s;
    GLfloat rt = -rs;

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(  0, -1, 0 );
    glVertex3f( rc, rt, -0.5 );
    glVertex3f(  1,  0, -1.5 );
    glVertex3f(  1,  0,  0.5 );
    glVertex3f(  0,  0,  1.5 );
    glEnd();
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(  0, +1, 0 );
    glVertex3f( rc, rs, -0.5 );
    glVertex3f(  0,  0,  1.5 );
    glVertex3f(  1,  0,  0.5 );
    glVertex3f(  1,  0, -1.5 );
    glEnd();

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(  s,  d, 0 );
    glVertex3f( rc, rt, -0.5 );
    glVertex3f(  0,  0,  1.5 );
    glVertex3f(  d,  t,  0.5 );
    glVertex3f(  d,  t, -1.5 );
    glEnd();
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(  t, c, 0 );
    glVertex3f( -r, 0, -0.5 );
    glVertex3f(  d, t, -1.5 );
    glVertex3f(  d, t,  0.5 );
    glVertex3f(  0, 0,  1.5 );
    glEnd();

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(  s, c, 0 );
    glVertex3f( rc, rs, -0.5 );
    glVertex3f(  d,  s, -1.5 );
    glVertex3f(  d,  s,  0.5 );
    glVertex3f(  0,  0,  1.5 );
    glEnd();
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(  t, d, 0 );
    glVertex3f( -r, 0, -0.5 );
    glVertex3f(  0, 0,  1.5 );
    glVertex3f(  d, s,  0.5 );
    glVertex3f(  d, s, -1.5 );
    glEnd();
    
    // closing the bottom gaps
    glBegin(GL_TRIANGLES);
    glNormal3f(  c,  t, -1 );
    glVertex3f( rc, rs, -0.5 );
    glVertex3f( -r,  0, -0.5 );
    glVertex3f(  d,  s, -1.5 );

    glNormal3f(  c,  s, -1 );
    glVertex3f( -r,  0, -0.5 );
    glVertex3f( rc, rt, -0.5 );
    glVertex3f(  d,  t, -1.5 );
    
    glNormal3f( -1,  0, -1 );
    glVertex3f( rc, rt, -0.5 );
    glVertex3f( rc, rs, -0.5 );
    glVertex3f(  1,  0, -1.5 );
    glEnd();
}


GLfloat dumbbellRadius(GLfloat z)
{
    return sinf(M_PI*z) * ( 1.3 + cosf(2*M_PI*z) );
}


void gle::gleDumbbell1()
{
    gleRevolution(dumbbellRadius);
}

GLfloat barrelRadius(GLfloat z)
{
    return sin(M_PI*z);
}

void gle::gleBarrel1()
{
    gleRevolution(barrelRadius);
}


//-----------------------------------------------------------------------
#pragma mark - Primitives with Display Lists

#ifdef GLE_USES_DISPLAY_LISTS

/// index of OpenGL display lists
GLuint gle::dlist = 0;
GLuint gle::slist = 0;



void makeDisplayList(GLint list, void (*primitive)())
{
    if ( glIsList(list) )
    {
        glNewList(list, GL_COMPILE);
        primitive();
        glEndList();
    }
    else
    {
        fprintf(stderr, "gle::makeDisplayList failed\n");
    }
}


void callDisplayList(GLint list)
{
    assert_true( glIsList(list) );
    glCallList(list);
}

void gle::initializeDL()
{
    if ( dlist == 0 )
    {
        dlist = glGenLists(12);
        
        if ( dlist == 0 )
        {
            fprintf(stderr, "gle::glGenLists(12) failed\n");
            gleReportErrors(stderr, "");
        }
        else
        {
            makeDisplayList(dlist+0, gleCircleL);
            makeDisplayList(dlist+1, gleCircleS);
            makeDisplayList(dlist+2, gleTube1);
            makeDisplayList(dlist+3, gleTube2);
            makeDisplayList(dlist+4, gleLongTube1);
            makeDisplayList(dlist+5, gleLongTube2);
            makeDisplayList(dlist+6, gleCone1);
            makeDisplayList(dlist+7, gleCylinderZ);
            makeDisplayList(dlist+8, gleDumbbell1);
            makeDisplayList(dlist+9, gleIcosahedron1);
            makeDisplayList(dlist+10, gleArrowTail1);
            makeDisplayList(dlist+11, gleArrowTail2);
        }
    }

#ifndef GLE_USE_VBO
    if ( slist == 0 )
    {
        slist = glGenLists(4);
        if ( slist == 0 )
        {
            fprintf(stderr, "gle::glGenLists(4) failed\n");
            gleReportErrors(stderr, "");
        }
        else
        {
            makeDisplayList(slist,   gleSphere1);
            makeDisplayList(slist+1, gleSphere2);
            makeDisplayList(slist+2, gleSphere4);
            makeDisplayList(slist+3, gleSphere8);
        }
    }
#endif
}

void gle::releaseDL()
{
    if ( dlist > 0 )
        glDeleteLists(dlist, 15);
}

#endif


//-----------------------------------------------------------------------
#pragma mark - Vertex Buffer Objects

#ifdef GLE_USE_VBO

GLuint ico_buffer[8] = { 0 };
GLuint ico_nfaces[4] = { 0 };


void initializeBuffers(GLuint& nfaces, GLuint buffer1, GLuint buffer2, int fin)
{
    Platonic::Solid ico(Platonic::Solid::ICOSAHEDRON, fin);
    nfaces = ico.nb_faces();
    //std::clog << " ico " << nfaces << std::endl;
    
    // upload vertex data:
    glBindBuffer(GL_ARRAY_BUFFER, buffer1);
    glBufferData(GL_ARRAY_BUFFER, 3*ico.nb_vertices()*sizeof(float), ico.vertex_data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    // upload indices:
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer2);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*ico.nb_faces()*sizeof(unsigned), ico.faces_data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


void gle::initializeVBO()
{
    if ( !glIsBuffer(ico_buffer[0]) )
    {
        glGenBuffers(8, ico_buffer);
        initializeBuffers(ico_nfaces[0], ico_buffer[0], ico_buffer[1], finesse/4);
        initializeBuffers(ico_nfaces[1], ico_buffer[2], ico_buffer[3], finesse/2);
        initializeBuffers(ico_nfaces[2], ico_buffer[4], ico_buffer[5], finesse);
        initializeBuffers(ico_nfaces[3], ico_buffer[6], ico_buffer[7], finesse*2);
    }
}


void gle::releaseVBO()
{
    if ( ico_buffer[0] && glIsBuffer(ico_buffer[0]) )
        glDeleteBuffers(8, ico_buffer);
}


void gleSphereVBO(GLuint nfaces, GLuint buffer1, GLuint buffer2)
{
    assert_true( GL_TRUE == glIsBuffer(buffer1) );
    assert_true( GL_TRUE == glIsBuffer(buffer2) );
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    glBindBuffer(GL_ARRAY_BUFFER, buffer1);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glNormalPointer(GL_FLOAT, 0, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer2);
    glDrawElements(GL_TRIANGLES, 3*nfaces, GL_UNSIGNED_INT, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
}

    
void gle::gleSphere1B()   { gleSphereVBO(ico_nfaces[0], ico_buffer[0], ico_buffer[1]); }
void gle::gleSphere2B()   { gleSphereVBO(ico_nfaces[1], ico_buffer[2], ico_buffer[3]); }
void gle::gleSphere4B()   { gleSphereVBO(ico_nfaces[2], ico_buffer[4], ico_buffer[5]); }
void gle::gleSphere8B()   { gleSphereVBO(ico_nfaces[3], ico_buffer[6], ico_buffer[7]); }


#elif defined GLE_USES_DISPLAY_LISTS

void gle::gleSphere1B()   { callDisplayList(slist+0); }
void gle::gleSphere2B()   { callDisplayList(slist+1); }
void gle::gleSphere4B()   { callDisplayList(slist+2); }
void gle::gleSphere8B()   { callDisplayList(slist+3); }

#else

void gle::gleSphere1B()   { gleSphere1(); }
void gle::gleSphere2B()   { gleSphere2(); }
void gle::gleSphere4B()   { gleSphere4(); }
void gle::gleSphere8B()   { gleSphere8(); }

#endif


//-----------------------------------------------------------------------
#pragma mark -

/**
 Draw a surface of revolution around the Z-axis.
 The surface goes from Z=0 to Z=1, and its radius is
 given by the function `radius`(z) provided as argument.
 */
void gle::gleRevolution(GLfloat (*radius)(GLfloat))
{
    GLfloat r0, z0, z1=0, r1=radius(z1), dr, dn;
    GLfloat dz = 0.25 / (GLfloat) finesse;
    
    GLfloat s[2*finesse+1], c[2*finesse+1];
    for ( int ii = 0; ii <= 2*finesse; ++ii )
    {
        GLfloat ang = ii * M_PI / (GLfloat) finesse;
        c[ii] = cosf(ang);
        s[ii] = sinf(ang);
    }
    
    for ( int jj = 0; jj <= 4*finesse; ++jj )
    {
        z0 = z1;
        r0 = r1;
        z1 = jj * dz;
        r1 = radius(z1);
        
        dr = ( r1 - r0 ) / dz;
        dn = 1.0 / sqrt( 1 + dr * dr );
        dr = dr*dn;
        
        glBegin(GL_TRIANGLE_STRIP);
        for ( int ii = 0; ii <= 2*finesse; ++ii )
        {
            glNormal3f(dn*c[ii], dn*s[ii], -dr);
            glVertex3f(r1*c[ii], r1*s[ii], z1);
            glVertex3f(r0*c[ii], r0*s[ii], z0);
        }
        glEnd();
    }
}

//-----------------------------------------------------------------------
#pragma mark - Object Placement


/**
 draw back first, and then front of object,
 GL_CULL_FACE should be enabled
 */
void gle::gleDualPass(void primitive())
{
    assert_true(glIsEnabled(GL_CULL_FACE));
    glCullFace(GL_FRONT);
    primitive();
    glCullFace(GL_BACK);
    primitive();
}


void gle::gleObject( const real radius, void (*obj)() )
{
    glPushMatrix();
    gleScale(radius);
    obj();
    glPopMatrix();
}

void gle::gleObject( const Vector1 & x, const real radius, void (*obj)() )
{
    glPushMatrix();
    gleTranslate(x);
    gleScale(radius);
    obj();
    glPopMatrix();
}

void gle::gleObject( const Vector2 & x, const real radius, void (*obj)() )
{
    glPushMatrix();
    gleTranslate(x);
    gleScale(radius);
    obj();
    glPopMatrix();
}

void gle::gleObject( const Vector3 & x, const real radius, void (*obj)() )
{
    glPushMatrix();
    gleTranslate(x);
    gleScale(radius);
    obj();
    glPopMatrix();
}


//-----------------------------------------------------------------------
void gle::gleObject(const Vector1 & a, const Vector1 & b, void (*obj)())
{
    glPushMatrix();
    if ( a.XX < b.XX )
        glRotated(  90, 0.0, 1.0, 0.0 );
    else
        glRotated( -90, 0.0, 1.0, 0.0 );
    obj();
    glPopMatrix();
}

void gle::gleObject(const Vector2 & a, const Vector2 & b, void (*obj)())
{
    glPushMatrix();
    gleAlignZ(a, b);
    obj();
    glPopMatrix();
}

void gle::gleObject(const Vector3 & a, const Vector3 & b, void (*obj)())
{
    glPushMatrix();
    Vector3 dir = b-a;
    const real dn = dir.norm();
    Vector3 P1  = dir.orthogonal(dn);
    Vector3 P2  = cross(dir, P1) / dn;
    gleTransRotate( P1, P2, dir, a );
    obj();
    glPopMatrix();
}


//-----------------------------------------------------------------------
void gle::gleObject( const Vector1 & x, const Vector1 & d, const real r, void (*obj)() )
{
    glPushMatrix();
    gleTranslate(x);
    if ( d.XX < 0 )
        glRotated(90, 0, 0, 1);
    gleScale(r);
    obj();
    glPopMatrix();
}

void gle::gleObject( const Vector2 & x, const Vector2 & d, const real r, void (*obj)() )
{
    glPushMatrix();
    gleAlignZ(x, x+d.normalized(r), r);
    obj();
    glPopMatrix();
}

void gle::gleObject( const Vector3 & x, const Vector3 & d, const real r, void (*obj)() )
{
    glPushMatrix();
    Vector3 P1 = d.orthogonal(r);
    Vector3 P2 = cross(d.normalized(), P1);
    gleTransRotate( P1, P2, d.normalized(r), x );
    obj();
    glPopMatrix();
}

//-----------------------------------------------------------------------
void gle::gleObject( const Vector1 & x, const Vector1 & d, const real r,
                     const real l, void (*obj)() )
{
    glPushMatrix();
    gleTranslate(x);
    if ( d.XX < 0 )
        glRotated(90, 0, 0, 1);
    gleScale(l,r,r);
    obj();
    glPopMatrix();
}

void gle::gleObject( const Vector2 & x, const Vector2 & d, const real r,
                     const real l, void (*obj)() )
{
    glPushMatrix();
    gleAlignZ(x, x+d.normalized(l), r);
    obj();
    glPopMatrix();
}

void gle::gleObject( const Vector3 & x, const Vector3 & d, const real r,
                     const real l, void (*obj)() )
{
    glPushMatrix();
    Vector3 P1 = d.orthogonal(r);
    Vector3 P2 = cross(d.normalized(), P1);
    gleTransRotate(P1, P2, d.normalized(l), x);
    obj();
    glPopMatrix();
}

//-----------------------------------------------------------------------
#pragma mark - Tubes


void gle::gleTube(const Vector1 & a, const Vector1 & b, real radius, void (*obj)())
{
    glPushMatrix();
    if ( a.XX < b.XX )
        glRotated(  90, 0.0, 1.0, 0.0 );
    else
        glRotated( -90, 0.0, 1.0, 0.0 );
    gleScale(1,radius,1);
    obj();
    glPopMatrix();
}

void gle::gleTube(const Vector2 & a, const Vector2 & b, real radius, void (*obj)())
{
    glPushMatrix();
    gleAlignZ(a, b, radius);
    obj();
    glPopMatrix();
}

void gle::gleTube(const Vector3 & a, const Vector3 & b, real radius, void (*obj)())
{
    glPushMatrix();
    Vector3 dir = b-a;
    Vector3 P1  = dir.orthogonal(radius);
    Vector3 P2  = cross(dir, P1).normalized(radius);
    gleTransRotate( P1, P2, dir, a );
    obj();
    glPopMatrix();
}


//-----------------------------------------------------------------------
void gle::gleTube(const Vector1 & a, real ra, gle_color ca,
                  const Vector1 & b, real rb, gle_color cb)
{
    glPushMatrix();
    gleTranslate(-a);
    if ( a.XX < b.XX )
        glRotated(  90, 0.0, 1.0, 0.0 );
    else
        glRotated( -90, 0.0, 1.0, 0.0 );
    gleTubeZ(a.XX, ra, ca, b.XX, rb, cb);
    glPopMatrix();
}

void gle::gleTube(const Vector2 & a, real ra, gle_color ca,
                  const Vector2 & b, real rb, gle_color cb)
{
    glPushMatrix();
    gleAlignZ(a, b);
    //gleTube1();
    gleTubeZ(0, ra, ca, 1, rb, cb);
    glPopMatrix();
}

void gle::gleTube(const Vector3 & a, real ra, gle_color ca,
                  const Vector3 & b, real rb, gle_color cb)
{
    glPushMatrix();
    Vector3 dir = b-a;
    Vector3 P1  = dir.orthogonal(1);
    Vector3 P2  = cross(dir, P1).normalized();
    gleTransRotate(P1, P2, dir, a);
    gleTubeZ(0, ra, ca, 1, rb, cb);
    glPopMatrix();
}


//-----------------------------------------------------------------------

void gle::gleBand(const Vector2 & a, const Vector2 & b, real rad)
{
    Vector2 d = ( b - a ).orthogonal();
    real n = d.norm();
    if ( n > 0 )
    {
        rad /= n;
        glBegin(GL_TRIANGLE_STRIP);
        gleVertex(a+rad*d);
        gleVertex(a-rad*d);
        gleVertex(b+rad*d);
        gleVertex(b-rad*d);
        glEnd();
    }
}


void gle::gleBand(const Vector1 & a, real ra,
                  const Vector1 & b, real rb)
{
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(a.XX,+ra);
    gleVertex(a.XX,-ra);
    gleVertex(b.XX,+rb);
    gleVertex(b.XX,-rb);
    glEnd();
}

void gle::gleBand(const Vector2 & a, real ra,
                  const Vector2 & b, real rb)
{
    Vector2 d = ( b - a ).orthogonal();
    real n = d.norm();
    if ( n > 0 )
    {
        d /= n;
        glBegin(GL_TRIANGLE_STRIP);
        gleVertex(a-ra*d);
        gleVertex(a+ra*d);
        gleVertex(b-rb*d);
        gleVertex(b+rb*d);
        glEnd();
    }
}

void gle::gleBand(const Vector1 & a, real ra, gle_color ca,
                  const Vector1 & b, real rb, gle_color cb)
{
    glBegin(GL_TRIANGLE_STRIP);
    ca.load();
    gleVertex(a.XX,-ra);
    gleVertex(a.XX,+ra);
    cb.load();
    gleVertex(b.XX,-rb);
    gleVertex(b.XX,+rb);
    glEnd();
}

void gle::gleBand(const Vector2 & a, real ra, gle_color ca,
                  const Vector2 & b, real rb, gle_color cb)
{
    Vector2 d = ( b - a ).orthogonal();
    real n = d.norm();
    if ( n > 0 )
    {
        d /= n;
        glBegin(GL_TRIANGLE_STRIP);
        ca.load();
        gleVertex(a+ra*d);
        gleVertex(a-ra*d);
        cb.load();
        gleVertex(b+rb*d);
        gleVertex(b-rb*d);
        glEnd();
    }
}


/**
 This will displays a rectangle if the connection is parallel,
 and a hourglass if the connection is antiparallel
 */
void gle::gleMan(Vector2 const& a, Vector2 const& da,
                 Vector2 const& b, Vector2 const& db)
{
#if ( 1 )
    Vector2 vec[6] = { b-db, b, a-da, a+da, b, b+db };
    assert_true(sizeof(vec)==12*sizeof(double));
    glVertexPointer(2, GL_DOUBLE, 0, vec);
    glEnableClientState(GL_VERTEX_ARRAY);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);
    glDisableClientState(GL_VERTEX_ARRAY);
#else
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(b-db);
    gleVertex(b);
    gleVertex(a-da);
    gleVertex(a+da);
    gleVertex(b);
    gleVertex(b+db);
    glEnd();
#endif
}

/**
 This will displays a rectangle if the connection is parallel,
 and a hourglass if the connection is antiparallel
 */
void gle::gleMan(Vector2 const& a, Vector2 const& da, gle_color ca,
                 Vector2 const& b, Vector2 const& db, gle_color cb)
{
#if ( 1 )
    GLfloat col[6*4];
    cb.put_floats(col);
    cb.put_floats(col+4);
    ca.put_floats(col+8);
    ca.put_floats(col+12);
    cb.put_floats(col+16);
    cb.put_floats(col+20);
    
    Vector2 vec[6] = { b-db, b, a-da, a+da, b, b+db };
    assert_true(sizeof(vec)==12*sizeof(double));
    
    glVertexPointer(2, GL_DOUBLE, 0, vec);
    glEnableClientState(GL_VERTEX_ARRAY);
    glColorPointer(4, GL_FLOAT, 0, col);
    glEnableClientState(GL_COLOR_ARRAY);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
#else
    glBegin(GL_TRIANGLE_STRIP);
    cb.load();
    gleVertex(b-db);
    gleVertex(b);
    ca.load();
    gleVertex(a-da);
    gleVertex(a+da);
    cb.load();
    gleVertex(b);
    gleVertex(b+db);
    glEnd();
#endif
}


/**
 This will displays a rectangle if the connection is antiparallel,
 and a hourglass if the connection is parallel
 */
void gle::gleCross(Vector2 const& a, Vector2 const& da,
                   Vector2 const& b, Vector2 const& db, real rad)
{
    glLineWidth(0.5);
    glBegin(GL_TRIANGLE_FAN);
    gleVertex(a);
    gleVertex(a-rad*da);
    gleVertex(b);
    gleVertex(b-rad*db);
    glEnd();
    glBegin(GL_TRIANGLE_FAN);
    gleVertex(a);
    gleVertex(a+rad*da);
    gleVertex(b);
    gleVertex(b+rad*db);
    glEnd();
}

void gle::gleBar(Vector3 const& a, Vector3 const& da,
                 Vector3 const& b, Vector3 const& db, real rad)
{
    Vector3 ab = ( a - b ).normalized();
    Vector3 ea = cross(ab, da);
    Vector3 eb = cross(ab, db);
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(a-rad*(da-ea));
    gleVertex(a-rad*(da+ea));
    gleVertex(b-rad*(db-eb));
    gleVertex(b-rad*(db+eb));
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(a+rad*(da-ea));
    gleVertex(a+rad*(da+ea));
    gleVertex(b+rad*(db-eb));
    gleVertex(b+rad*(db+eb));
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(a-rad*da);
    gleVertex(a+rad*da);
    gleVertex(b-rad*db);
    gleVertex(b+rad*db);
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(a-rad*da);
    gleVertex(a+rad*da);
    gleVertex(b-rad*db);
    gleVertex(b+rad*db);
    glEnd();
}


/**
 Two hexagons linked by a rectangle
 */
void gle::gleDumbbell(const Vector2 & a, const Vector2 & b, real diameter)
{
    //side of hexagon that has the same surface as the disc or Radius 1.
    const GLfloat S = sqrt( 2 * M_PI / ( 3 * sqrt(3) ));
    const GLfloat R = diameter * S;
    const GLfloat H = R * 0.5 * sqrt(3);
    const GLfloat X = R * 0.5;

    Vector2 x = ( b - a ).normalized(H);
    Vector2 y = x.orthogonal(X);
    
    glPushMatrix();
    gleTranslate(a);
    
    // this is an hexagon centered around 'a':
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(0,0);
    gleVertex(x+y);
    gleVertex(2*y);
    gleVertex(-x+y);
    gleVertex(-x-y);
    gleVertex(-2*y);
    gleVertex(x-y);
    gleVertex(x+y);
    glEnd();
    
    // a band from 'a' to 'b'
    glBegin(GL_TRIANGLE_FAN);
    gleVertex(+y+x);
    gleVertex(-y+x);
    gleVertex(b-a-y-x);
    gleVertex(b-a+y-x);
    glEnd();
    
    // an hexagon centered around 'b'
    gleTranslate(b-a);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(0,0);
    gleVertex(x+y);
    gleVertex(2*y);
    gleVertex(-x+y);
    gleVertex(-x-y);
    gleVertex(-2*y);
    gleVertex(x-y);
    gleVertex(x+y);
    glEnd();
    
    glPopMatrix();
}

//-----------------------------------------------------------------------
#pragma mark - Arrows

void gle::gleCone(const Vector1 & center, const Vector1 & dir, const real scale)
{
    real dx = scale*dir.XX, cx = center.XX;
    glBegin(GL_TRIANGLES);
    gleVertex( cx+dx+dx, 0 );
    gleVertex( cx-dx,    dx );
    gleVertex( cx-dx,   -dx );
    glEnd();
}

void gle::gleCone(const Vector2 & center, const Vector2 & dir, const real scale)
{
    real dx = scale*dir.XX,  cx = center.XX;
    real dy = scale*dir.YY,  cy = center.YY;
    glBegin(GL_TRIANGLES);
    gleVertex( cx+dx+dx, cy+dy+dy );
    gleVertex( cx-dx-dy, cy-dy+dx );
    gleVertex( cx-dx+dy, cy-dy-dx );
    glEnd();
}

void gle::gleCone(const Vector3 & center, const Vector3 & dir, const real scale)
{
    glPushMatrix();
    //build the rotation matrix, assuming dir is normalized
    Vector3 P1 = dir.orthogonal(scale);
    Vector3 P2 = cross(dir, P1);
    gleTransRotate( P1, P2, dir*scale, center );
    gleCone1B();
    glPopMatrix();
}


void gle::gleConeL(const Vector1 & center, const Vector1 & dir, const real scale)
{
    real dx = scale*dir.XX, cx = center.XX;
    glBegin(GL_LINE_STRIP);
    gleVertex( cx-dx,    dx );
    gleVertex( cx+dx+dx, 0 );
    gleVertex( cx-dx,   -dx );
    glEnd();
}

void gle::gleConeL(const Vector2 & center, const Vector2 & dir, const real scale)
{
    real dx = scale*dir.XX,  cx = center.XX;
    real dy = scale*dir.YY,  cy = center.YY;
    glBegin(GL_LINE_STRIP);
    gleVertex( cx-dx-dy, cy-dy+dx );
    gleVertex( cx+dx+dx, cy+dy+dy );
    gleVertex( cx-dx+dy, cy-dy-dx );
    glEnd();
}

void gle::gleConeL(const Vector3 & center, const Vector3 & dir, const real scale)
{
    glPushMatrix();
    //build the rotation matrix, assuming dir is normalized
    Vector3 P1 = dir.orthogonal(scale);
    Vector3 P2 = cross(dir, P1);
    gleTransRotate( P1, P2, dir*scale, center );
    gleCone1L();
    glPopMatrix();
}

//-----------------------------------------------------------------------

void gle::gleCylinder(const Vector1 & center, const Vector1 & dir, const real scale)
{
    real cx = center.XX;
    glBegin(GL_TRIANGLE_STRIP);
    real dx = 0.5 * scale * dir.XX;
    gleVertex( cx-dx, -scale );
    gleVertex( cx-dx,  scale );
    gleVertex( cx+dx, -scale );
    gleVertex( cx+dx,  scale );
    glEnd();
}

void gle::gleCylinder(const Vector2 & center, const Vector2 & dir, const real scale)
{
    real dx = scale * dir.XX, cx = center.XX - 0.5 * dx;
    real dy = scale * dir.YY, cy = center.YY - 0.5 * dy;
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex( cx+dy, cy-dx );
    gleVertex( cx-dy, cy+dx );
    gleVertex( cx+dx+dy, cy+dy-dx );
    gleVertex( cx+dx-dy, cy+dy+dx );
    glEnd();
}

void gle::gleCylinder(const Vector3 & center, const Vector3 & dir, const real scale)
{
    glPushMatrix();
    //build the rotation matrix, assuming dir is normalized
    Vector3 P1 = dir.orthogonal(scale);
    Vector3 P2 = cross(dir, P1);
    gleTransRotate( P1, P2, dir*scale, center );
    gleCylinderB();
    glPopMatrix();
}


//-----------------------------------------------------------------------

void gle::gleArrowTail(const Vector1 & center, const Vector1 & dir, const real scale)
{
    GLfloat dx = scale * dir.XX;
    GLfloat cx = center.XX - 0.5 * dx;
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f( cx,       0  );
    glVertex2f( cx-dx,   -dx );
    glVertex2f( cx+dx,   -dx );
    glVertex2f( cx+dx+dx, 0  );
    glVertex2f( cx+dx,    dx );
    glVertex2f( cx-dx,    dx );
    glEnd();
}

void gle::gleArrowTail(const Vector2 & center, const Vector2 & dir, const real scale)
{
    GLfloat dx = scale * dir.XX;
    GLfloat dy = scale * dir.YY;
    GLfloat cx = center.XX - 1.5 * dx;
    GLfloat cy = center.YY - 1.5 * dy;
    GLfloat ex = cx + 2 * dx;
    GLfloat ey = cy + 2 * dy;
    
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f( cx+dx, cy+dy );
    glVertex2f( cx+dy, cy-dx );
    glVertex2f( ex+dy, ey-dx );
    glVertex2f( ex+dx, ey+dy );
    glVertex2f( ex-dy, ey+dx );
    glVertex2f( cx-dy, cy+dx );
    glEnd();
}

void gle::gleArrowTail(const Vector3 & center, const Vector3 & dir, const real scale)
{
    glPushMatrix();
    //build the rotation matrix, assuming dir is normalized
    Vector3 P1 = dir.orthogonal(scale);
    Vector3 P2 = cross(dir, P1);
    gleTransRotate( P1, P2, dir*scale, center );
    gleArrowTail2B();
    glPopMatrix();
}

//-----------------------------------------------------------------------
void gle::gleArrow(const Vector1 & a, const Vector1 & b, real radius)
{
    glPushMatrix();
    if ( a.XX < b.XX )
        glRotated(  90, 0.0, 1.0, 0.0 );
    else
        glRotated( -90, 0.0, 1.0, 0.0 );
    gleScale(1,radius,1);
    gleTube1B();
    glTranslatef(0, 0, 1);
    glScalef(3.0, 3.0, 3*radius);
    gleCone1B();
    glPopMatrix();
}

void gle::gleArrow(const Vector2 & a, const Vector2 & b, real radius)
{
    glPushMatrix();
    gleAlignZ(a, b, radius);
    gleTube1B();
    glTranslatef(0, 0, 1);
    glScalef(3.0, 3.0, 3*radius);
    gleCone1B();
    glPopMatrix();
}

void gle::gleArrow(const Vector3 & a, const Vector3 & b, real radius)
{
    glPushMatrix();
    Vector3 dir = b-a;
    Vector3 P1  = dir.orthogonal(radius);
    Vector3 P2  = cross(dir, P1).normalized(radius);
    gleTransRotate( P1, P2, dir, a );
    gleTube1B();
    glTranslatef(0, 0, 1);
    glScalef(3.0, 3.0, 3*radius);
    gleCone1B();
    glPopMatrix();
}


//-----------------------------------------------------------------------
#pragma mark - Text


int gle::gleLineHeight(void* font)
{
    if ( font == GLUT_BITMAP_8_BY_13 )        return 13;
    if ( font == GLUT_BITMAP_9_BY_15 )        return 15;
    if ( font == GLUT_BITMAP_TIMES_ROMAN_10 ) return 11;
    if ( font == GLUT_BITMAP_TIMES_ROMAN_24 ) return 26;
    if ( font == GLUT_BITMAP_HELVETICA_10 )   return 11;
    if ( font == GLUT_BITMAP_HELVETICA_12 )   return 15;
    if ( font == GLUT_BITMAP_HELVETICA_18 )   return 22;
    return 13;
}


/**
 Compute the max width of all the lines in the given text
 This uses GLUT, which should be initialized.
*/
int gle::gleComputeTextSize(const char text[], void* font, int& lines)
{
    int width = 0;
    lines = 0;
    int w = 0;
    for (const char* c = text; *c != '\0' ; ++c)
    {
        if ( *c == '\n' )
        {
            if ( w > width ) width = w;
            ++lines;
            w = 0;
        }
        else if ( isspace(*c))
        {
            w += glutBitmapWidth(font, ' ');
        }
        else if ( isprint(*c))
        {
            w += glutBitmapWidth(font, *c);
        }
    }
    if ( w > width )
        width = w;
    if ( width > 0 && lines == 0 )
        lines = 1;
    return width;
}

//-----------------------------------------------------------------------
/**
 draw the string character per character using:
 glutBitmapCharacter()
 */
void gle::gleDrawText(const char text[], void* font, GLfloat vshift)
{
    if ( font == 0 )
    {
        font = GLUT_BITMAP_HELVETICA_12;
        vshift = ( vshift > 0 ? 1 : -1 ) * gleLineHeight(font);
    }
    if ( vshift == 0 )
        vshift = -gleLineHeight(font);

    GLfloat ori[4], pos[4];
    glGetFloatv(GL_CURRENT_RASTER_POSITION, ori);
    
    for (const char* p = text; *p; ++p)
    {
        if ( *p == '\n' )
        {
            glGetFloatv(GL_CURRENT_RASTER_POSITION, pos);
            glBitmap(0, 0, 0, 0, ori[0]-pos[0], vshift, 0);
        }
        else if ( isspace(*p) )
        {
            glutBitmapCharacter(font, ' ');
        }
        else if ( isprint(*p) )
        {
            glutBitmapCharacter(font, *p);
        }
    }
}


/**
 set the current raster position to `w`
 */
void gle::gleDrawText(const Vector3& vec, const char text[], void* font)
{
    glPushAttrib(GL_LIGHTING_BIT|GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    int lh = gleLineHeight(font);
    gleRasterPos(vec);
    //translate to center the bitmap:
    glBitmap(0,0,0,0,1,-lh/3,0);
    gleDrawText(text, font, -lh);
    glPopAttrib();
}

void gle::gleDrawText(const Vector1& w, const char text[], void* font)
{    
    gleDrawText(Vector3(w.XX, 0, 0), text, font);
}

void gle::gleDrawText(const Vector2& w, const char text[], void* font)
{
    gleDrawText(Vector3(w.XX, w.YY, 0), text, font);
}

//-----------------------------------------------------------------------

/**
 The text is displayed in the current color.
 A background rectangle is displayed only if `bcol` is visible.
 
 @code
 glColor3f(1,1,1);
 gleDisplayText(fKeyString, GLUT_BITMAP_8_BY_13, 0x0, 1);
 @endcode
 
 Possible values for `position`:
 - 0: bottom-left, text going up
 - 1: bottom-right, text going up
 - 2: top-right, text going down
 - 3: top-left, text going down
 - 4: center, text going down
 .
 
 Note: width and height are the current size of the viewport (window)
 */
void gle::gleDisplayText(const char text[], void* font, const gle_color bcol,
                         const int position, int width, int height)
{
    assert_true( width > 0 );
    assert_true( height > 0 );
    
    if ( font == 0 )
        font = GLUT_BITMAP_9_BY_15;
    
    int lineHeight = gleLineHeight(font);
    int textWidth = 0;
    int nblines = 1;

    GLint px, py;
    switch( position )
    {
        case 0: {
            //bottom-left, text going up
            px = 0.5*lineHeight;
            py = 0.5*lineHeight;
        } break;
        case 1: {
            //bottom-right, text going up
            textWidth = gleComputeTextSize(text, font, nblines);
            px = width - textWidth - 0.5*lineHeight;
            if ( px < 0 ) px = 0;
            py = 0.5*lineHeight;
        } break;
        case 2: {
            //top-right, text going down
            textWidth = gleComputeTextSize(text, font, nblines);
            px = width - textWidth - 0.5*lineHeight;
            if ( px < 0 ) px = 0;
            py = height - lineHeight;
            lineHeight *= -1;
        } break;
        default:
        case 3: {
            //top-left, text going down
            px = 0.5*lineHeight;
            py = height - lineHeight;
            lineHeight *= -1;
        } break;
        case 4: {
            //center, text going down
            textWidth = gleComputeTextSize(text, font, nblines);
            px = ( width - textWidth ) / 2;
            if ( px < 0 ) px = 0;
            py = ( height + nblines*lineHeight ) / 2;
            lineHeight *= -1;
        } break;
    }
    
    //set pixel coordinate system:
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, width, 0, height, 0, 1);

    glRasterPos2i(0, 0);
    glBitmap(0, 0, 0, 0, px, py, 0);
    
    glPushAttrib(GL_LIGHTING_BIT|GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    
    if ( bcol.visible() )
    {
        glPushAttrib(GL_CURRENT_BIT);
        int rd = std::abs(lineHeight);
        int bt = py;
        int tp = py + nblines*lineHeight;
        if ( lineHeight < 0 )
        {
            int x = tp;
            tp = bt;
            bt = x;
        }

        int rect[4] = { px-rd, bt, px+textWidth+rd, tp+1.75*rd };

        bcol.load();
        glBegin(GL_TRIANGLE_FAN);
        gleNiceRectangle(rect, 4);
        glEnd();
        
        glPopAttrib();
        
        glLineWidth(0.5);
        glBegin(GL_LINE_STRIP);
        gleNiceRectangle(rect, 4);
        glEnd();
    }
    
    gleDrawText(text, font, lineHeight);

    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glPopAttrib();
}


//-----------------------------------------------------------------------
#pragma mark - Misc

/**
 Draw an array of pixels using GL_TRIANGLE_STRIP
 
 The array rgba[] should ( nbc * width * height ) bytes,
 containing nbc-components (eg. RGBA) per pixel and
 stored by columns:
 
 @code
 load(i,j) = rgba[ nbc*(i+height*j) ]
 0 <= i < height
 0 <= j < width
 @endcode
 
 `pos` is the position of the top-left corner
 `dx` is the direction of the width
 `dy` is the direction of the height
 The magnitudes of `dx` and `dy` indicates the dimensions of a pixel.
 They may be of different magnitudes, and not necessarily orthogonal.
 */
void gle::gleDrawPixels(int width, int height, int nbc, GLubyte rgba[], Vector2 pos, Vector2 dx, Vector2 dy)
{
    glPushAttrib(GL_ENABLE_BIT|GL_POLYGON_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glDisable(GL_LIGHTING);

    Vector2 left, right;
    GLubyte * col = rgba;
    
    for ( int jj = 0; jj < width; ++jj )
    {
        left  = pos + dx * jj;
        right = left + dx;
        for ( int ii = 0; ii < height; ++ii )
        {
            if ( nbc == 3 )
                glColor3ubv(col);
            else
                glColor4ubv(col);
            glBegin(GL_TRIANGLE_STRIP);
            col += nbc;
            gle::gleVertex(left);
            gle::gleVertex(right);
            left  += dy;
            right += dy;
            gle::gleVertex(left);
            gle::gleVertex(right);
            glEnd();
        }
    }
    
    glPopAttrib();
}

//-----------------------------------------------------------------------


/**
 rectangle should be specified as [ left, bottom, right, top ]
 The rectangle will be drawn counter-clockwise
*/
void gle::gleRectangle(const int rect[4])
{
    glVertex2i(rect[0], rect[1]);
    glVertex2i(rect[2], rect[1]);
    glVertex2i(rect[2], rect[3]);
    glVertex2i(rect[0], rect[3]);
    glVertex2i(rect[0], rect[1]);
}


void gle::gleNiceRectangle(const int rect[4], const int rad)
{
    glVertex2i(rect[0], rect[1]+rad);
    glVertex2i(rect[0]+rad, rect[1]);
    glVertex2i(rect[2]-rad, rect[1]);
    glVertex2i(rect[2], rect[1]+rad);
    glVertex2i(rect[2], rect[3]-rad);
    glVertex2i(rect[2]-rad, rect[3]);
    glVertex2i(rect[0]+rad, rect[3]);
    glVertex2i(rect[0], rect[3]-rad);
    glVertex2i(rect[0], rect[1]+rad);
}


void gle::gleDrawRectangle(const int rect[4], int width, int height)
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, width, 0, height, 0, 1);
    //disable advanced features
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    
    glBegin(GL_LINE_LOOP);
    gleRectangle(rect);
    glEnd();
    
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glPopAttrib();
}


void gle::gleDrawResizeBox(int width, int height)
{
    //set the matrices
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    
    glOrtho(width, 0, 0, height, 0, 1 );
    
    //draw lines at 45 degrees
    glBegin(GL_LINES);
    glVertex2i(16, 1);    glVertex2i(1, 16);
    glVertex2i(12, 1);    glVertex2i(1, 12);
    glVertex2i(8,  1);    glVertex2i(1,  8);
    glVertex2i(4,  1);    glVertex2i(1,  4);
    glEnd();
    
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}


//-----------------------------------------------------------------------
void gle::gleDrawAxes(const GLfloat S, int dim)
{
    const GLfloat R = S * 0.1;
    
    glMatrixMode(GL_MODELVIEW);
    
    for( int ii=0; ii<dim; ++ii)
    {
        glPushMatrix();
        switch( ii )
        {
            case 0: 
                gle_color(1.0, 0.0, 0.0, 1.0).load_load();
                glRotatef(+90, 0, 1, 0);
                break;
            case 1:
                gle_color(0.0, 1.0, 0.0, 1.0).load_load();
                glRotatef(-90, 1, 0, 0);
                break;
            case 2:
                gle_color(0.0, 0.0, 1.0, 1.0).load_load();
                break;
        }
        glScalef(R/4, R/4, S-R);
        gleTube1();
        glTranslatef(0, 0, 1);
        glScalef(4, 4, R/(S-R));
        gleCone1();
        glRasterPos3f(0, 0, 2);
        glNormal3f(0, 0, 1);
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, 'X'+ii);
        glPopMatrix();
    }
    
    // display a white ball at the origin
    gle_color(1.0, 1.0, 0.0, 1.0).load_load();
    glPushMatrix();
    gleScale(R);
    gleSphere4();
    glPopMatrix();
}

//-----------------------------------------------------------------------
char const* gle::gleErrorString(GLenum code)
{
    switch ( code )
    {
        case GL_NO_ERROR:          return "GL_NO_ERROR";
        case GL_INVALID_ENUM:      return "GL_INVALID_ENUM";
        case GL_INVALID_VALUE:     return "GL_INVALID_VALUE";
        case GL_INVALID_OPERATION: return "GL_INVALID_OPERATION";
        case GL_STACK_OVERFLOW:    return "GL_STACK_OVERFLOW";
        case GL_STACK_UNDERFLOW:   return "GL_STACK_UNDERFLOW";
        case GL_OUT_OF_MEMORY:     return "GL_OUT_OF_MEMORY";
        case GL_TABLE_TOO_LARGE:   return "GL_TABLE_TOO_LARGE";
        default:                   return "GL_UNKNOWN_ERROR";
    }
}

/**
 This is similart to glutReportError,
 but the additional argument can provide useful feedback for debugging
 */
void gle::gleReportErrors(FILE * out, const char* msg)
{
    GLenum glError = glGetError();
    while ( glError != GL_NO_ERROR )
    {
        fprintf(out, "OpenGL error `%s' %s\n", gleErrorString(glError), msg);
        glError = glGetError();
    }
}


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Visual test for the rasterizer used in attachment algorithm of Cytosim
// Created by Francois Nedelec, October 2002

#include "dim.h"
#include <ctime>
#include "vector.h"
#include "glapp.h"
#include "glut.h"
#include "real.h"
#include "smath.h"
#include "random.h"
#include "rasterizer.h"
#include "vector.h"

extern Random RNG;
extern bool rasterizer_draw_things;

const int size = 50;
const unsigned MAX = 16;
unsigned n_pts = 2;

real radius  = 10;
real shift[] = {0, 0, 0};
real delta[] = {1, 1, 1};
Vector pts[MAX];

#if ( DIM == 3 )
int hit[2*size+1][2*size+1][2*size+1];
#else
int hit[2*size+1][2*size+1];
#endif

//-------------------------------------------------------------------

void newPoints()
{
    for ( unsigned i = 0; i < MAX ; ++i )
        pts[i] = (size-1) * Vector::srand();
}


void clearHits()
{
    for ( int i = 0; i <= 2*size; ++i )
    for ( int j = 0; j <= 2*size; ++j )
#if ( DIM == 2 )
        hit[i][j] = 0;
#elif ( DIM == 3 )
    for ( int k = 0; k <= 2*size; ++k )
        hit[i][j][k] = 0;
#endif
}


void paintHit(int x_inf, int x_sup, int y, int z, void*)
{
    if ( x_sup < x_inf )
        printf("mixedup order x = %i %i at y = %i\n", x_inf, x_sup, y);
    
    if ( y < -size ||  y > size ) return;
    if ( z < -size ||  z > size ) return;
    
    for ( int x = x_inf; x <= x_sup; ++x )
    {
        if ( -size <= x && x <= size )
        {
#if ( DIM == 3 )
            ++hit[x+size][y+size][z+size];
#else
            ++hit[x+size][y+size];
#endif
        }
    }
}

void paintDraw(int x_inf, int x_sup, int y, int z, void*)
{
    glColor4f(0, 1, 0, 0.5);

    glPointSize(3);
    glBegin(GL_POINTS);
    glVertex3i( x_inf, y, z );
    //glVertex3i( x_sup, y, z );
    glEnd();
    
    glLineWidth(1);
    glBegin(GL_LINES);
    glVertex3i( x_inf, y, z );
    glVertex3i( x_sup, y, z );
    glEnd();
}

void rasterize(Vector P, Vector Q, void (paint)(int, int, int, int, void*))
{
    real len = ( P - Q ).norm();
#if ( DIM == 2 )
    Rasterizer::paintFatLine2D(paint, 0, P, Q, len, radius, shift, delta);
    //Rasterizer::paintFatLine2D(paint, 0, P, Q, len, radius);
    //Rasterizer::paintBox2D(paint, 0, P, Q, radius, shift, delta);
#elif ( DIM == 3 )
    //Rasterizer::paintHexLine3D(paint, 0, P, Q, len, radius, shift, delta);
    Rasterizer::paintFatLine3D(paint, 0, P, Q, len, radius, shift, delta);
    //Rasterizer::paintBox3D(paint, 0, P, Q, radius, shift, delta);
#endif
}

//-------------------------------------------------------------------


bool inCylinder(Vector const& p, Vector const& q, Vector const& x)
{
    const Vector pq = q - p;
    const real pqn = pq.norm();
    real abs = pq * (x-p) / pqn;
    if ( abs <   0-radius ) return false;
    if ( abs > pqn+radius ) return false;
    abs /= pqn;
    if ( abs < 0 ) abs = 0;
    if ( abs > 1 ) abs = 1;
    return ( x-p-pq*abs ).normSqr() <= radius * radius;
}


bool checkPoint(Vector const& p, Vector const& q, int i, int j, int k)
{
    int in = inCylinder(p,q,Vector(i,j,k));
    
#if ( DIM == 3 )
    int ht = hit[i+size][j+size][k+size];
#else
    int ht = hit[i+size][j+size];
#endif
    
    if ( ht != in )
    {
        glPointSize(5);
        glBegin(GL_POINTS);
        if ( ht == 1 )
            glColor4f(0, 0, 1, 0.5);
        else
            glColor3f(1, 0, 0);
#if ( DIM == 3 )
        glVertex3i(i, j, k);
#else
        glVertex2i(i, j);
#endif
        glEnd();
    }
    
    return ( ht != in );
}


bool testFatLine(Vector P, Vector Q)
{
    bool res = false;
    clearHits();
    rasterize(P, Q, paintHit);
    for ( int i = -size; i <= size; ++i )
    for ( int j = -size; j <= size; ++j )
#if ( DIM == 2 )
        res |= checkPoint(P, Q, i, j, 0);
#elif ( DIM == 3 )
    for ( int k = -size; k <= size; ++k )
        res |= checkPoint(P, Q, i, j, k);
#endif
    return res;
}


void manyTest()
{
    unsigned ouf = 0;
    do {
        pts[0] = (size-1) * Vector::srand();;
        pts[1] = (size-1) * Vector::srand();;
        if ( testFatLine(pts[0], pts[1]) )
            break;
    } while ( ++ouf < 1<<10 );
}

//===================================================================


void processNormalKey(unsigned char c, int x=0, int y=0)
{
    switch (c) {
        case 27:
        case 'q':
            exit(EXIT_SUCCESS);
        case ' ':
            newPoints();
            break;
        case '0':
            glApp::resetView();
            break;
        case 'p': if ( n_pts+1 < MAX ) ++n_pts; break;
        case 'o': if ( n_pts > 2 ) --n_pts; break;
        case 'r':
            manyTest();
            break;
        case 'h':
            printf("keyboard commands:\n"
                   " space : draw a new random distribution\n"
                   " p     : increase number of points\n"
                   " o     : decrease number of points\n"
                   " r     : perform many tests as fast as possible\n");
        default:
            glApp::processNormalKey(c, 0, 0);
    }
    glApp::postRedisplay();
}

//===================================================================


void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    //--------------draw points on the grid:
#if ( DIM == 2 )
    glPointSize(1.0);
    glColor3f(0.5, 0.5, 0.5);
    glBegin(GL_POINTS);
    for ( int i = -size; i <= size; i += 1)
    for ( int j = -size; j <= size; j += 1)
        glVertex2i(i, j);
    glEnd();
#endif

    //--------------draw a grid in gray:
#if ( DIM == 2 )
    glLineWidth(0.5);
    glColor3f(0.5, 0.5, 0.5);
    glBegin(GL_LINES);
    for ( int i = -size; i <= size; i += 5)
    {
        glVertex2i(i, -size);
        glVertex2i(i, +size);
        glVertex2i(-size, i);
        glVertex2i(+size, i);
    }
    glEnd();
#endif
    
    /// draw large points:
    glPointSize(16);
    glBegin(GL_POINTS);
    glColor3f(0, 0, 1);
    for ( unsigned i = 0; i < n_pts ; ++i )
    {
#if ( DIM == 2 )
        glVertex2d(pts[i].XX, pts[i].YY);
#else
        glVertex3d(pts[i].XX, pts[i].YY, pts[i].ZZ);
#endif
    }
    glEnd();

#if ( DIM == 2 )
    if ( n_pts > 2 )
    {
        glLineWidth(1);
        unsigned int nb = Rasterizer::convexHull2D(n_pts, pts);
        Rasterizer::paintPolygon2D(paintDraw, 0, nb, pts, 0);
        
        // draw convex-hull
        glLineWidth(2);
        glColor3f(1, 1, 0);
        glBegin(GL_LINE_LOOP);
        glVertex2d(pts[0].XX, pts[0].YY);
        glColor3f(0, 0, 1);
        for ( unsigned i = 1; i < nb ; ++i )
           glVertex2d(pts[i].XX, pts[i].YY);
        glVertex2d(pts[0].XX, pts[0].YY);
        glEnd();
        
        /// draw first points:
        glPointSize(16);
        glBegin(GL_POINTS);
        glColor3f(1, 0, 1);
#if ( DIM == 2 )
        glVertex2d(pts[0].XX, pts[0].YY);
#else
        glVertex3d(pts[0].XX, pts[0].YY, pts[0].ZZ);
#endif
        glEnd();

        return;
    }
#endif
    

    Vector P = pts[0];
    Vector Q = pts[1];
    testFatLine(P,Q);
    rasterize(P,Q,paintDraw);
}

/* 
 This only work if rasterizer does not issue openGL commands */
void speedTest(unsigned cnt)
{
    clearHits();
    
    for ( unsigned n = 1; n < n_pts; ++n )
    {
        Vector P = pts[n-1];
        Vector Q = pts[n];
        real len = ( P - Q ).norm();
        for ( unsigned ii = 0; ii < cnt; ++ii )
        {
#if ( DIM == 2 )
            Rasterizer::paintFatLine2D(paintHit, 0, P, Q, len, radius, shift, delta);
#elif ( DIM == 3 )
            Rasterizer::paintFatLine3D(paintHit, 0, P, Q, len, radius, shift, delta);
#endif
        }
    }
}

//===================================================================

int main(int argc, char* argv[])
{
    RNG.seedTimer();
    newPoints();
    if ( argc > 1 )
    {
        rasterizer_draw_things = 0;
        speedTest(strtol(argv[1], 0, 10));
    }
    else
    {
        rasterizer_draw_things = 1;
        glutInit(&argc, argv);
        glApp::setDimensionality(3);
        glApp::attachMenu(GLUT_RIGHT_BUTTON);
        glApp::normalKeyFunc(processNormalKey);
        glApp::setScale(2*(size+radius+1));
        glApp::createWindow(display);

        glutMainLoop();
    }
    return EXIT_SUCCESS;
}



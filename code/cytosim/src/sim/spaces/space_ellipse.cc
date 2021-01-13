// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_ellipse.h"
#include "exceptions.h"
#include "project_ellipse.h"
#include "smath.h"


inline real sqrE(const real x) { return x*x; }


SpaceEllipse::SpaceEllipse(const SpaceProp* p)
: Space(p)
{
#ifdef HAS_SPHEROID
    mSpheroid = -1;
#endif
}


void SpaceEllipse::resize()
{
    Space::checkLengths(DIM, true);
    
#if ( DIM == 3 ) && defined HAS_SPHEROID
    mSpheroid = -1;
    
    // if any two dimensions are similar, then the ellipsoid is a spheroid
    for ( int zz = 0; zz < DIM; ++zz )
    {
        int xx = ( zz + 1 ) % DIM;
        int yy = ( zz + 2 ) % DIM;
        if ( fabs( (length(xx)-length(yy)) / (length(xx)+length(yy)) ) < REAL_EPSILON )
            mSpheroid = zz;
    }
#endif
}


void SpaceEllipse::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length(0),-length(1),-length(2));
    sup.set( length(0), length(1), length(2));
}


Vector SpaceEllipse::normalToEdge(Vector const& pos) const
{
#if (DIM == 1)
    return Vector(sign(pos.XX), 0);
#elif (DIM == 2)
    return Vector(pos.XX/length(0), pos.YY/length(1));
#else
    return Vector(pos.XX/length(0), pos.YY/length(1), pos.ZZ/length(2));
#endif
}

#if (DIM == 1)

real SpaceEllipse::volume() const
{
    return 2 * length(0);
}

bool  SpaceEllipse::inside( const real w[] ) const
{
    return (( w [0] >= -length(0) ) && ( w [0] <=  length(0) ));
}

#elif (DIM == 2)

real SpaceEllipse::volume() const
{
    return M_PI * length(0) * length(1);
}

bool  SpaceEllipse::inside( const real w[] ) const
{
    return ( sqrE( w[0] / length(0) ) + sqrE( w[1] / length(1) ) <= 1 );
}

#else

real SpaceEllipse::volume() const
{
    return 4/3.0 * M_PI * length(0) * length(1) * length(2);
}

bool SpaceEllipse::inside( const real w[] ) const
{
    return  sqrE( w[0]/length(0) ) + sqrE( w[1]/length(1) ) + sqrE( w[2]/length(2)) <= 1;
}

#endif



void SpaceEllipse::project1D( const real w[], real p[] ) const
{
    if ( w[0] >= 0 )
        p[0] =  length(0);
    else
        p[0] = -length(0);
}



void SpaceEllipse::project2D( const real w[], real p[] ) const
{
    projectEllipse(p[0], p[1], w[0], w[1], mLength[0], mLength[1]);
#if ( 0 )
    // check that results are valid numbers:
    assert_true(p[0]==p[0]);
    assert_true(p[1]==p[1]);
#endif
}


void SpaceEllipse::project3D( const real w[], real p[] ) const
{
#if ( DIM == 3 ) && defined HAS_SPHEROID
    /*
     If the ellipsoid has two equal axes, we can reduce the problem to 2D,
     because it is symmetric by rotation around the remaining axis, which
     is here indicated by 'mSpheroid'.
     */
    if ( mSpheroid >= 0 )
    {
        const int zz = mSpheroid;
        const int xx = ( zz + 1 ) % DIM;
        const int yy = ( zz + 2 ) % DIM;
        
        if ( length(xx) != length(yy) )
            throw InvalidParameter("Inconsistent mSpheroid value");
        
        //rotate point around the xx axis to bring it into the yy-zz plane:
        real pR, rr = sqrt( w[xx]*w[xx] + w[yy]*w[yy] );
        projectEllipse(pR, p[zz], rr, w[zz], length(xx), length(zz), 8*REAL_EPSILON);
        // back-rotate to get the projection in 3D:
        if ( rr > 0 ) {
            real s = pR / rr;
            p[xx] = w[xx] * s;
            p[yy] = w[yy] * s;
        }
        else {
            p[xx] = 0;
            p[yy] = 0;
        }
        return;
    }
#endif
    
    projectEllipsoid(p, w, mLength);
    
#if ( 0 )
    // check that results are valid numbers:
    assert_true(p[0]==p[0]);
    assert_true(p[1]==p[1]);
    assert_true(p[2]==p[2]);
#endif
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

// Modified from sphere by Aastha Mathur, 18th Jan 2013

#ifdef DISPLAY

#include "glut.h"
#include "gle.h"

bool SpaceEllipse::display() const
{
#if ( DIM == 1 )
    
    GLfloat X = length(0);
    
    glBegin(GL_LINES);
    glVertex2f(-X, -1);
    glVertex2f(-X,  1);
    glVertex2f( X, -1);
    glVertex2f( X,  1);
    glEnd();

#elif ( DIM == 2 )
    
    GLfloat X = length(0);
    GLfloat Y = length(1);
    glBegin(GL_LINE_LOOP);
    for ( real aa = 0; aa < 6.28; aa += 0.01 )
        glVertex2f( X*cosf(aa), Y*sinf(aa) );
    glEnd();
    
#elif ( DIM == 3 )
    
    const unsigned fin = ((DIM==2) ? 32 : 8) * gle::finesse;
    GLfloat X = length(0);
    GLfloat Y = length(1);
    GLfloat Z = length(2);
    
    GLfloat c[2*fin+1], s[2*fin+1];
    for ( unsigned ii = 0; ii <= 2*fin; ++ii )
    {
        GLfloat ang = ii * M_PI / (GLfloat) fin;
        c[ii] = cosf(ang);
        s[ii] = sinf(ang);
    }
    
    for ( unsigned ii = 0; ii < fin; ++ii )
    {
        real uX = s[ii  ]*X, uY = s[ii  ]*Y, uZ = c[ii  ]*Z;
        real lX = s[ii+1]*X, lY = s[ii+1]*Y, lZ = c[ii+1]*Z;
        glBegin(GL_TRIANGLE_STRIP);
        for ( unsigned jj = 0; jj <= 2*fin; ++jj )
        {
            glNormal3f(c[jj]*s[ii], s[jj]*s[ii], c[ii]);
            glVertex3f(c[jj]*uX, s[jj]*uY, uZ);
            glNormal3f(c[jj]*s[ii+1], s[jj]*s[ii+1], c[ii+1]);
            glVertex3f(c[jj]*lX, s[jj]*lY, lZ);
        }
        glEnd();
    }

    if ( 0 )
    {
        // add decorations:
        glPushMatrix();
        gle::gleScale(length(0), length(1), length(2));
        gle::gleThreeBands(128);
        glPopMatrix();
    }

    if ( 1 )
    {
        glLineWidth(1);
        // add decorations:
        glPushMatrix();
        gle::gleScale(length(0), length(1), length(2));
        for ( real Z = 0.1; Z < 1; Z += 0.2 )
        {
            real R = sqrt( 1 - Z*Z );
            glPushMatrix();
            glScaled(R, R, 1);
            glTranslated(0, 0, Z);
            gle::gleCircleL();
            glTranslated(0, 0, -2*Z);
            gle::gleCircleL();
            glPopMatrix();
        }
        glPopMatrix();
    }
#endif
    
    return true;
}

#else

bool SpaceEllipse::display() const
{
    return false;
}


#endif

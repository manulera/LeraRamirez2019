// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_sphere.h"
#include "exceptions.h"
#include "random.h"
#include "smath.h"
#include "meca.h"

extern Random RNG;

SpaceSphere::SpaceSphere(const SpaceProp* p)
: Space(p), radius(mLength[0]), radiusSqr(mLengthSqr[0])
{
}


void SpaceSphere::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius,-radius,-radius);
    sup.set( radius, radius, radius);
}


#if (DIM == 1)


real SpaceSphere::volume() const
{
    return 2 * radius;
}

bool SpaceSphere::inside( const real point[] ) const
{
    return point[0] * point[0] <= radiusSqr;
}

void SpaceSphere::project( const real point[], real proj[] ) const
{
    real n = point[0] * point[0];
    
    if ( n > 0 ) {
        n = radius / sqrt( n );
        proj[0] = n * point[0];
    }
    else {
        proj[0] = RNG.sflip() * radius;
    }
}

#endif


//------------------------------------------------------------------------------


#if (DIM == 2)

real SpaceSphere::volume() const
{
    return M_PI * radius * radius;
}

bool SpaceSphere::inside( const real point[] ) const
{
    return point[0] * point[0] + point[1] * point[1] <= radiusSqr;
}

void SpaceSphere::project( const real point[], real proj[] ) const
{    
    real n = point[0] * point[0] + point[1] * point[1];
    
    if ( n > 0 ) {
        n = radius / sqrt( n );
        proj[0] = n * point[0];
        proj[1] = n * point[1];
    }
    else {
        //select a random point on the surface
        real x, y;
        do {
            x = RNG.sreal();
            y = RNG.sreal();
            n = x*x + y*y;
        } while ( n > 1.0  ||  n == 0 );
        n = radius / sqrt( n );
        proj[0] = n * x;
        proj[1] = n * y;
    }
}

#endif


//------------------------------------------------------------------------------

#if (DIM == 3)

real SpaceSphere::volume() const
{
    return 4/3.0 * M_PI * radius * radius * radius;
}

bool SpaceSphere::inside( const real point[] ) const
{
    return point[0] * point[0] + point[1] * point[1] + point[2] * point[2] <= radiusSqr;
}

void SpaceSphere::project( const real point[], real proj[] ) const
{    
    real n = point[0] * point[0] + point[1] * point[1] + point[2] * point[2];
    
    if ( n > 0 ) {
        n = radius / sqrt( n );
        proj[0] = n * point[0];
        proj[1] = n * point[1];
        proj[2] = n * point[2];
    }
    else {
        //select a random point on the surface
        real x, y, z;
        do {
            x = RNG.sreal();
            y = RNG.sreal();
            z = RNG.sreal();
            n = x*x + y*y + z*z;
        } while ( n > 1.0  ||  n == 0 );
        n = radius / sqrt( n );
        proj[0] = n * x;
        proj[1] = n * y;
        proj[2] = n * z;
    }
}

#endif

//------------------------------------------------------------------------------

void SpaceSphere::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    meca.addLongPointClamp( pos, pe, Vector(0,0,0), radius, stiff );
}


void SpaceSphere::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    if ( radius > rad )
        meca.addLongPointClamp( pos, pe, Vector(0,0,0), radius-rad, stiff );
    else {
        meca.addPointClamp( pe, Vector(0,0,0), stiff );
        std::cerr << "object is too big to fit in SpaceSphere\n";
    }
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY

#include "glut.h"
#include "gle.h"

bool SpaceSphere::display() const
{
    GLfloat R = radius;

#if ( DIM <= 2 )
   
    const GLfloat da = M_PI / 360;
    
    glBegin(GL_LINE_LOOP);
    for ( real a = -M_PI; a <= M_PI;  a += da )
        glVertex2f(R*cosf(a), R*sinf(a));
    glEnd();

#else
    
    glPushMatrix();
    glScalef(R, R, R);
    gle::gleSphere8();
    gle::gleThreeBands(128);
    glPopMatrix();
    
#endif
    
    return true;
}

#else

bool SpaceSphere::display() const
{
    return false;
}


#endif

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_cylinderP.h"
#include "point_exact.h"
#include "exceptions.h"
#include "smath.h"
#include "meca.h"


SpaceCylinderP::SpaceCylinderP(const SpaceProp* p)
: Space(p), length(mLength[0]), radius(mLength[1]), radiusSqr(mLengthSqr[1])
{
    if ( DIM != 3 )
        throw InvalidParameter("cylinderP is only valid in 3D: use strip instead");
}

void SpaceCylinderP::setModulo(Modulo * mod) const
{
    mod->enable(0, length);
}

void SpaceCylinderP::resize()
{
    Space::checkLengths(2, false);
    
    if ( length <= 0 )
        throw InvalidParameter("length of cylinderP must be > 0");
}


void SpaceCylinderP::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length,-radius,-radius);
    sup.set( length, radius, radius);
}



real  SpaceCylinderP::volume() const
{
    return 2 * M_PI * length * radius * radius;
}


bool  SpaceCylinderP::inside( const real w[] ) const
{
    return ( w[1]*w[1]+ w[2]*w[2] <= radiusSqr );
}


bool  SpaceCylinderP::allInside( const real w[], const real rad ) const
{
    assert_true( rad > 0 );
    
    return ( sqrt( w[1]*w[1]+ w[2]*w[2] ) + rad <= radius );
}


Vector SpaceCylinderP::randomPlace() const
{
#if ( DIM == 3 )
    Vector2 sec = Vector2::randB(radius);
    return Vector(length*RNG.sreal(), sec.XX, sec.YY);
#else
    return Vector(length*RNG.sreal(), radius*RNG.sreal());
#endif
}

//------------------------------------------------------------------------------
void SpaceCylinderP::project( const real w[], real p[] ) const
{
    p[0] = w[0];
    
    real n = sqrt( w[1]*w[1]+ w[2]*w[2] );
    if ( n > REAL_EPSILON )
    {
        p[1] = radius * w[1] / n;
        p[2] = radius * w[2] / n;
    }
    else
    {
        real a = M_PI * RNG.sreal();
        p[1] = radius * sin(a);
        p[2] = radius * cos(a);
    }
}

//------------------------------------------------------------------------------

/**
 This applies forces towards the cylindrical surface only
 */
void SpaceCylinderP::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    meca.addLongPointClampYZ(pe, radius, stiff);
}

/**
 This applies forces towards the cylindrical surface only
 */
void SpaceCylinderP::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    real eRadius = radius - rad;
    if ( eRadius < 0 ) eRadius = 0;
    
    meca.addLongPointClampYZ(pe, eRadius, stiff);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

bool SpaceCylinderP::display() const
{
#if ( DIM == 3 )

    const int fin = 512;
    
    GLfloat L = length;
    GLfloat R = radius;
    
    glBegin(GL_TRIANGLE_STRIP);
    for ( int ii = 0; ii <= fin; ++ii )
    {
        GLfloat ang = ii * 2 * M_PI / (GLfloat) fin;
        GLfloat ca = cosf(ang), sa = sinf(ang);
        glNormal3f( 0, ca, sa );
        glVertex3f( +L, R*ca, R*sa );
        glVertex3f( -L, R*ca, R*sa );
    }
    glEnd();
    
    if ( 1 )
    {
        //draw dotted-rings to indicate periodicity
        glLineStipple(2, 0x0303);
        glEnable(GL_LINE_STIPPLE);
        glPushMatrix();
        gle::gleTranslate(L, 0, 0);
        gle::gleScale(R, R, R);
        glRotated(90, 0, 1, 0);
        gle::gleCircleL();
        gle::gleTranslate(0, 0, -2*L/R);
        glRotated(180, 0, 1, 0);
        gle::gleCircleL();
        glPopMatrix();
        glDisable(GL_LINE_STIPPLE);
    }

#endif
    return true;
}

#else

bool SpaceCylinderP::display() const
{
    return false;
}

#endif


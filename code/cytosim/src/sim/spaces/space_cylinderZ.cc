// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_cylinderZ.h"
#include "point_exact.h"
#include "exceptions.h"
#include "smath.h"
#include "meca.h"


SpaceCylinderZ::SpaceCylinderZ(const SpaceProp* p)
: Space(p), radius(mLength[0]), radiusSqr(mLengthSqr[0]), bottom(mLength[1]), top(mLength[2])
{
    if ( DIM != 3 )
        throw InvalidParameter("cylinderZ is only valid in 3D: use sphere instead");
}


void SpaceCylinderZ::resize()
{
    Space::checkLengths(1, false);
    if ( top < bottom )
        throw InvalidParameter("bottom must be lower than top ( size[1] <= size[2] )");
}


void SpaceCylinderZ::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius,-radius, bottom);
    sup.set( radius, radius, top);
}



real  SpaceCylinderZ::volume() const
{
    return 2 * M_PI * ( top - bottom ) * radius * radius;
}


bool  SpaceCylinderZ::inside( const real w[] ) const
{
    if ( w[2] < bottom || w[2] > top )
        return false;
    return ( w[0]*w[0] + w[1]*w[1] <= radiusSqr );
}

bool  SpaceCylinderZ::allInside( const real w[], const real rad ) const
{
    assert_true( rad > 0 );
    
    if ( w[2] < bottom || w[2] > top )
        return false;
    return ( sqrt( w[0]*w[0]+ w[1]*w[1] ) + rad <= radius );
}

Vector SpaceCylinderZ::randomPlace() const
{
    Vector2 sec = Vector2::randB(radius);
    return Vector(sec.XX, sec.YY, bottom+RNG.preal()*(top-bottom));
}

//------------------------------------------------------------------------------
void SpaceCylinderZ::project( const real w[], real p[] ) const
{
    int inZ = 1;
    
    p[0] = w[0];
    p[1] = w[1];
    p[2] = w[2];
    
    if ( w[2] > top )
    {
        p[2] = top;
        inZ = 0;
    }
    else if ( w[2] < bottom )
    {
        p[2] = bottom;
        inZ = 0;
    }
    
    real n = sqrt( w[0]*w[0]+ w[1]*w[1] );
    
    if ( n > radius )
    {
        n = radius / n;
        p[0] = n * w[0];
        p[1] = n * w[1];            
    }
    else
    {
        if ( inZ )
        {
            if ( top - w[2] < radius - n )
                p[2] = top;
            else if ( w[2] - bottom < radius - n )
                p[2] = bottom;
            else
            {
                n = radius / n;
                p[0] = n * w[0];
                p[1] = n * w[1];
            }
        }
    }
}

//------------------------------------------------------------------------------

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinderZ::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff,
                                    const real rad, const real bot, const real top)
{
#if ( DIM == 3 )
    bool cap = false;
    bool cyl = false;
    real p;

    // inside cylinder radius
    if ( 2 * pos.ZZ > top + bot )
    {
        p = top;
        cap = ( pos.ZZ > top );
    }
    else
    {
        p = bot;
        cap = ( pos.ZZ < bot );
    }
    
    real dis = pos.XX*pos.XX + pos.YY*pos.YY;
    
    if ( rad*rad < dis )
    {
        // outside cylinder in XY plane
        cyl = true;
    }
    else if ( ! cap )
    {
        // inside cylinder in XY plane and also inside in Z:
        if ( std::abs( pos.ZZ - p ) > rad - sqrt(dis) )
            cyl = true;
        else
            cap = true;
    }
    
    if ( cap )
    {
        const Matrix::index_type inx = 2 + DIM * pe.matIndex();
        meca.mC(inx, inx) -= stiff;
        meca.base(inx)    += stiff * p;
    }
    
    if ( cyl )
        meca.addLongPointClampXY(pe, rad, stiff);
#endif
}


/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinderZ::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, radius, bottom, top);
}

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinderZ::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    real R = radius - rad;
    if ( R < 0 ) R = 0;
    
    real T = top - rad;
    real B = bottom + rad;
    
    if ( B > T )
    {
        B = 0.5 * ( top + bottom );
        T = B;
    }
    
    setInteraction(pos, pe, meca, stiff, R, B, T);
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"

bool SpaceCylinderZ::display() const
{
#if ( DIM == 3 )
    
    const int fin = 512;
    GLfloat T = top;
    GLfloat B = bottom;
    GLfloat R = radius;
    
    GLfloat c[fin+1], s[fin+1];
    for ( int ii = 0; ii <= fin; ++ii )
    {
        GLfloat ang = ii * 2 * M_PI / (GLfloat) fin;
        c[ii] = cosf(ang);
        s[ii] = sinf(ang);
    }
    
    glBegin(GL_TRIANGLE_STRIP);
    //display strips along the side of the volume:
    for ( int sc = 0; sc <= fin; ++sc )
    {
        GLfloat ca = c[sc], sa = s[sc];
        glNormal3f(ca, sa, 0);
        glVertex3f(R*ca, R*sa, T);
        glVertex3f(R*ca, R*sa, B);
    }
    glEnd();
    
    // draw top cap:
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, +1);
    glVertex3f(0, 0,  T);
    for ( int sc = 0; sc <= fin; ++sc )
        glVertex3f(R*c[sc], R*s[sc], T);
    glEnd();
    
    // draw bottom cap:
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, -1);
    glVertex3f(0, 0,  B);
    for ( int sc = 0; sc <= fin; ++sc )
        glVertex3f(-R*c[sc], R*s[sc], B);
    glEnd();
    
#endif
    return true;
}

#else

bool SpaceCylinderZ::display() const
{
    return false;
}

#endif

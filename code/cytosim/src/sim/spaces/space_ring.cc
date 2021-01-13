// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_ring.h"
#include "point_exact.h"
#include "exceptions.h"
#include "smath.h"
#include "meca.h"


SpaceRing::SpaceRing(const SpaceProp* p)
: Space(p), length(mLength[0]), radius(mLength[1]), radiusSqr(mLengthSqr[1])
{
    if ( DIM != 3 )
        throw InvalidParameter("ring is only valid in 3D: use rectangle instead");
}


void SpaceRing::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length,-radius,-radius);
    sup.set( length, radius, radius);
}


real  SpaceRing::volume() const
{
    return 2 * M_PI * length * radius * radius;
}

//------------------------------------------------------------------------------
bool  SpaceRing::inside( const real w[] ) const
{
    if ( fabs(w[0]) > length )
        return false;
    return ( w[1]*w[1]+ w[2]*w[2] <= radiusSqr );
}

bool  SpaceRing::allInside( const real w[], const real rad ) const
{
    assert_true( rad > 0 );
    
    if ( fabs(w[0]) > length-rad )
        return false;
    return ( sqrt( w[1]*w[1]+ w[2]*w[2] ) + rad <= radius );
}

//------------------------------------------------------------------------------
/**
 Project always on the surface of the cylinder
 */
void SpaceRing::project( const real w[], real p[] ) const
{    
    if ( w[0] >  length )
        p[0] =  length;
    else if ( w[0] < -length )
        p[0] = -length;
    else
        p[0] = w[0];
    
    real n = sqrt( w[1]*w[1]+ w[2]*w[2] );
    
    if ( n > 0 )
    {
        n = radius / n;
        p[1] = n * w[1];
        p[2] = n * w[2];            
    }
    else
    {
        p[1] = radius;
        p[2] = 0;
    }
}

//------------------------------------------------------------------------------

/**
 This applies a force directed to the surface of the cylinder
 */
void SpaceRing::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff, const real len, const real rad)
{
    const Matrix::index_type inx = DIM * pe.matIndex();

    if ( pos.XX > len )
    {
        meca.mC(inx, inx) -= stiff;
        meca.base(inx)    += stiff * len;
    }
    else if ( pos.XX < -len )
    {
        meca.mC(inx, inx) -= stiff;
        meca.base(inx)    -= stiff * len;
    }
    
    meca.addLongPointClampYZ(pe, rad, stiff);
}


/**
 This applies a force directed to the surface of the cylinder
 */
void SpaceRing::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length, radius);
}

/**
 This applies a force directed to the surface of the cylinder
 */
void SpaceRing::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length, radius);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"

bool SpaceRing::display() const
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
    
#endif
    return true;
}

#else

bool SpaceRing::display() const
{
    return false;
}

#endif


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_capsule.h"
#include "point_exact.h"
#include "exceptions.h"
#include "smath.h"
#include "meca.h"



SpaceCapsule::SpaceCapsule(const SpaceProp* p)
: Space(p), length(mLength[0]), radius(mLength[1]), radiusSqr(mLengthSqr[1])
{
    if ( DIM == 1 )
        throw InvalidParameter("capsule is only defined for DIM = 2 or 3");
}


void SpaceCapsule::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius-length,-radius,-radius);
    sup.set( radius+length, radius, radius);
}


real SpaceCapsule::volume() const
{
#if (DIM == 3)
    return ( length + (2/3.0) * radius ) * radiusSqr * ( 2 * M_PI );
#else
    return 4 * length * radius + M_PI * radiusSqr;
#endif
}


bool SpaceCapsule::inside( const real w[] ) const
{
    real nrm, x = fabs( w[0] );
    
    if ( x > length )
        nrm = sqr( x - length );
    else
        nrm = 0;
    
#if ( DIM == 2 )
    nrm += w[1] * w[1];
#else
    nrm += w[1] * w[1] + w[2] * w[2];
#endif
    
    return ( nrm <= radiusSqr );
}


bool SpaceCapsule::allInside( const real w[], const real rad ) const
{
    assert_true( rad > 0 );
    
    real nrm, x = fabs( w[0] );
    
    if ( x > length )
        nrm = sqr( x - length );
    else
        nrm = 0;
    
#if ( DIM == 2 )
    nrm += w[1] * w[1];
#else
    nrm += w[1] * w[1] + w[2] * w[2];
#endif
    
    return ( nrm <= sqr(radius-rad) );
}

//------------------------------------------------------------------------------
void SpaceCapsule::project( const real w[], real p[] ) const
{
    real nrm;
#if ( DIM == 2 )
    nrm = w[1] * w[1];
#else
    nrm = w[1] * w[1] + w[2] * w[2];
#endif
    
    //calculate the projection on the axis, within boundaries:
    if ( w[0] >  length )
    {
        nrm  += sqr( w[0] - length );
        //normalize from this point on the axis
        if ( nrm > 0 ) nrm = radius / sqrt( nrm );
        
        p[0] = length + nrm * ( w[0] - length );
    }
    else
    {
        if ( w[0] < -length )
        {
            nrm  += sqr( length + w[0] );
            //normalize from this point on the axis
            if ( nrm > 0 ) nrm = radius / sqrt( nrm );
            
            p[0]  = -length + nrm * ( w[0] + length );
        }
        else
        {
            //normalize from this point on the axis
            if ( nrm > 0 ) nrm = radius / sqrt( nrm );
            
            p[0] = w[0];
        }
    }
    
    if ( nrm > 0 )
    {
        p[1] = nrm * w[1];
#if ( DIM == 3 )
        p[2] = nrm * w[2];
#endif
    }
    else
    {
        //we project on a arbitrary point on the cylinder
        p[1] = radius;
#if ( DIM == 3 )
        p[2] = 0;
#endif
    }
}



Vector SpaceCapsule::randomPlace() const
{
    unsigned long nb_trials = 1<<13;
    unsigned long ouf = 0;
    Vector res;

    do {
        
#if ( DIM == 3 )
        Vector2 sec = Vector2::randB(radius);
        res.set((length+radius)*RNG.sreal(), sec.XX, sec.YY);
#elif ( DIM == 2 )
        res.set((length+radius)*RNG.sreal(), radius*RNG.sreal());
#else
        res.set((length+radius)*RNG.sreal());
#endif
        
        if ( ++ouf > nb_trials )
        {
            std::clog << "placement failed in SpaceCapsule::randomPlace()" << std::endl;
            return Vector(0,0,0);
        }
        
    } while ( ! inside(res) );
    
    return res;
}

//------------------------------------------------------------------------------

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCapsule::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff, const real len, const real rad)
{
    if ( pos.XX > len )
        meca.addLongPointClamp( pos, pe, Vector( len,0,0), rad, stiff );
    else if ( pos.XX < -len )
        meca.addLongPointClamp( pos, pe, Vector(-len,0,0), rad, stiff );
    else
        meca.addLongPointClampYZ( pe, rad, stiff );
}


/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCapsule::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length, radius);
}

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCapsule::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    if ( rad < radius )
        setInteraction(pos, pe, meca, stiff, length, radius-rad);
    else
        setInteraction(pos, pe, meca, stiff, length, 0);
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

bool SpaceCapsule::display() const
{
    //number of sections in the quarter-circle
    const int fin = ((DIM==2) ? 32 : 8) * gle::finesse;
    
    GLfloat c[4*fin+1], s[4*fin+1];
    for ( int ii = 0; ii <= 4*fin; ++ii )
    {
        GLfloat ang = ii * M_PI_2 / (GLfloat) fin;
        c[ii] = cosf(ang);
        s[ii] = sinf(ang);
    }
    
    GLfloat L = length;
    GLfloat R = radius;
    
#if ( DIM <= 2 )
    
    //display a loop in X/Y plane
    glBegin(GL_LINE_LOOP);
    for ( int ii = 0;     ii <= 2*fin; ++ii )
        glVertex2f( +L+R*s[ii], R*c[ii] );
    for ( int ii = 2*fin; ii <= 4*fin; ++ii )
        glVertex2f( -L+R*s[ii], R*c[ii] );
    glEnd();
    
#else
    
    //display strips along the side of the volume:
    for ( int sc = 0; sc < 4*fin; ++sc )
    {
        //compute the transverse angles:
        GLfloat ctb  = c[sc  ],   stb  = s[sc  ];
        GLfloat cta  = c[sc+1],   sta  = s[sc+1];
        GLfloat ctbR = R*ctb,     stbR = R*stb;
        GLfloat ctaR = R*cta,     staR = R*sta;
        
        //draw one strip of the oval:
        glBegin(GL_TRIANGLE_STRIP);
        for ( int ii=0; ii <= fin; ++ii )
        {
            GLfloat ca = c[ii], sa = s[ii];
            glNormal3f( ca, cta*sa, sta*sa );
            glVertex3f( +L+R*ca, ctaR*sa, staR*sa );
            glNormal3f( ca, ctb*sa, stb*sa );
            glVertex3f( +L+R*ca, ctbR*sa, stbR*sa );
        }
        for ( int ii=fin; ii >= 0; --ii)
        {
            GLfloat ca = -c[ii], sa = s[ii];
            glNormal3f( ca, cta*sa, sta*sa );
            glVertex3f( -L+R*ca, ctaR*sa, staR*sa );
            glNormal3f( ca, ctb*sa, stb*sa );
            glVertex3f( -L+R*ca, ctbR*sa, stbR*sa );
        }
        glEnd();
    }
    
    if ( 1 )
    {
        //draw 2 rings on the surface
        glPushMatrix();
        gle::gleTranslate(L, 0, 0);
        gle::gleScale(R, R, R);
        glRotated(90, 0, 1, 0);
        gle::gleArrowedBand(24, 0.25);
        gle::gleTranslate(0, 0, -2*L/R);
        glRotated(180, 0, 1, 0);
        gle::gleArrowedBand(24, 0.25);
        glPopMatrix();
    }

#endif
    return true;
}

#else

bool SpaceCapsule::display() const
{
    return false;
}

#endif

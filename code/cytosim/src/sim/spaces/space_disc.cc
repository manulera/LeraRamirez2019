// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "space_disc.h"
#include "exceptions.h"
#include "random.h"
#include "smath.h"
#include "meca.h"

extern Random RNG;

SpaceDisc::SpaceDisc(const SpaceProp* p)
: Space(p), radius(mLength[0])
{
    if ( DIM != 2 )
        throw InvalidParameter("disc is only usable in 2D");
    rForce = 0;
}


void SpaceDisc::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius,-radius,-radius);
    sup.set( radius, radius, radius);
}


#if (DIM != 2)


real SpaceDisc::volume() const
{
    return 0;
}

bool SpaceDisc::inside( const real point[] ) const
{
    return false;
}

void SpaceDisc::project( const real point[], real proj[] ) const
{
    for ( int d = 0; d < DIM; ++d)
        proj[d] = 0;
}

#else


real SpaceDisc::volume() const
{
    return M_PI * radius * radius;
}

bool SpaceDisc::inside( const real point[] ) const
{
    return point[0] * point[0] + point[1] * point[1] <= radius * radius;
}

void SpaceDisc::project( const real point[], real proj[] ) const
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

/// add interactions to a Meca
void SpaceDisc::setInteractions(Meca &, FiberSet const&) const
{
    rForce = 0;
}


void SpaceDisc::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    meca.addLongPointClamp(pos, pe, Vector(0,0,0), radius, stiff);
    rForce += stiff * ( pos.norm() - radius );
}


void SpaceDisc::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    if ( radius > rad )
    {
        meca.addLongPointClamp(pos, pe, Vector(0,0,0), radius-rad, stiff);
        rForce += stiff * ( rad + pos.norm() - radius );
    }
    else {
        meca.addPointClamp( pe, Vector(0,0,0), stiff );
        std::cerr << "object is too big to fit in SpaceDisc\n";
        rForce += 2 * stiff * ( rad - radius );
    }
}


void SpaceDisc::step()
{
    real dr = prop->mobility_dt * rForce;
    std::clog << "SpaceDisc:  radius " << std::setw(12) << radius << " force " << rForce << " delta_radius " << dr << "\n";
    radius += dr;
}


#ifdef DISPLAY

#include "opengl.h"

bool SpaceDisc::display() const
{
#if ( DIM <= 2 )
    GLfloat R = radius;

    const GLfloat da = M_PI / 360;
    
    glBegin(GL_LINE_LOOP);
    for ( real a = -M_PI; a <= M_PI;  a += da )
        glVertex2f(R*cosf(a), R*sinf(a));
    glEnd();
    
#endif
    
    return true;
}

#else

bool SpaceDisc::display() const
{
    return false;
}


#endif

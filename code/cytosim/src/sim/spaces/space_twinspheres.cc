// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_twinspheres.h"
#include "exceptions.h"
#include "random.h"
#include "smath.h"
#include "meca.h"
#include "object_set.h"
#include "simul.h"

extern Random RNG;

SpaceTwinSpheres::SpaceTwinSpheres(const SpaceProp* p)
: Space(p), radiusL(mLength[0]), radiusLSqr(mLengthSqr[0]),
radiusR(mLength[1]), radiusRSqr(mLengthSqr[1]), overlap(mLength[2])
{
    if ( DIM == 1 )
        throw InvalidParameter("twin_spheres is not usable in 1D");
}


void SpaceTwinSpheres::resize()
{
    Space::checkLengths(3, true);
    
    real X = radiusL + radiusR - overlap;
    if ( X < 0 )
        throw InvalidParameter("twin_spheres:length[2] is too large");

    cenL = -( radiusL*radiusL - radiusR*radiusR + X*X ) / ( 2 * X );
    cenR = X + cenL;

    if ( cenL < -radiusL )
        throw InvalidParameter("twin_spheres has incompatible dimensions");
    if ( cenR >  radiusR )
        throw InvalidParameter("twin_spheres has incompatible dimensions");

    // so the radius of the base circle is  a=sqrt(h(2R-h)),
    neck = sqrt( radiusL*radiusL - cenL*cenL );
}


void SpaceTwinSpheres::boundaries(Vector& inf, Vector& sup) const
{
    real R = std::max(radiusL, radiusR);
    inf.set(cenL-radiusL, -R,-R);
    sup.set(cenR+radiusR,  R, R);
}


#if (DIM == 1)

real SpaceTwinSpheres::volume() const
{
    return 0;
}

bool SpaceTwinSpheres::inside( const real point[] ) const
{
    return false;
}

real SpaceTwinSpheres::projectS(real rad, real cen, const real point[], real proj[])
{
    proj[0] = 0;
    return 0;
}


real SpaceTwinSpheres::projectC(real rad, const real point[], real proj[])
{
    proj[0] = 0;
    return 0;
}

#endif


//------------------------------------------------------------------------------


#if (DIM == 2)

real SpaceTwinSpheres::volume() const
{
    /*
     The area of the circular segment is equal to the area of the
     circular sector minus the area of the triangular portion:
     A = \frac{R^2}{2} ( \theta - sin(\theta) ).
    */
    real thetaL = 2 * asin( neck / radiusL );
    real thetaR = 2 * asin( neck / radiusR );

    real capL = 0.5 * ( thetaL - sin(thetaL) );
    real capR = 0.5 * ( thetaR - sin(thetaR) );
    
    real res = 0;
    if ( cenL < 0 )
        res += radiusL*radiusL * ( M_PI - capL );
    else
        res += radiusL*radiusL * capL;
    
    if ( cenR > 0 )
        res += radiusR*radiusR * ( M_PI - capR );
    else
        res += radiusR*radiusR * capR;
    
    return res;
}

bool SpaceTwinSpheres::inside( const real point[] ) const
{
    if ( point[0] > 0 )
    {
        real X = point[0] - cenR;
        return X * X + point[1] * point[1] <= radiusRSqr;
    }
    else
    {
        real X = point[0] - cenL;
        return X * X + point[1] * point[1] <= radiusLSqr;
    }
}

real SpaceTwinSpheres::projectS(real rad, real cen, const real point[], real proj[])
{
    real x = point[0] - cen;
    real n = x * x + point[1] * point[1];
    
    if ( n > 0 )
    {
        n = sqrt(n);
        real s = rad / n;
        proj[0] = s * x + cen;
        proj[1] = s * point[1];
        return fabs(rad-n);
    }
    else
    {
        //select a random point on the surface
        real x, y;
        do {
            x = RNG.sreal();
            y = RNG.sreal();
            n = x*x + y*y;
        } while ( n > 1.0  ||  n == 0 );
        n = rad / sqrt(n);
        proj[0] = n * x + cen;
        proj[1] = n * y;
        return rad;
    }
}


real SpaceTwinSpheres::projectC(real rad, const real point[], real proj[])
{
    proj[0] = 0;
    proj[1] = ( point[1] > 0 ) ? rad : -rad;
    return sqrt( point[0] * point[0] + (point[1]-proj[1]) * (point[1]-proj[1]));
}


#endif


//------------------------------------------------------------------------------

#if (DIM == 3)

real SpaceTwinSpheres::volume() const
{
    /*
     h = epaisseur du cap, R = rayon sphere
     V_cap=1/3 * pi * h^2 * (3*R-h).
     */

    real hL = radiusL + cenL;
    real hR = radiusR - cenR;
    return M_PI / 3.0 * ( 4 * ( radiusL*radiusL*radiusL + radiusR*radiusR*radiusR )
                         - hL * hL * ( 3 * radiusL - hL ) - hR * hR * ( 3 * radiusR - hR ) );
}


bool SpaceTwinSpheres::inside( const real point[] ) const
{
    if ( point[0] > 0 )
    {
        real X = point[0] - cenR;
        return X * X + point[1] * point[1] + point[2] * point[2] <= radiusRSqr;
    }
    else
    {
        real X = point[0] - cenL;
        return X * X + point[1] * point[1] + point[2] * point[2] <= radiusLSqr;
    }
}

real SpaceTwinSpheres::projectS(real rad, real cen, const real point[], real proj[])
{
    real x = point[0] - cen;
    real n = x * x + point[1] * point[1] + point[2] * point[2];
    
    if ( n > 0 )
    {
        n = sqrt(n);
        real s = rad / n;
        proj[0] = s * x + cen;
        proj[1] = s * point[1];
        proj[2] = s * point[2];
        return fabs(rad-n);
    }
    else
    {
        //select a random point on the surface
        real x, y, z;
        do {
            x = RNG.sreal();
            y = RNG.sreal();
            z = RNG.sreal();
            n = x*x + y*y + z*z;
        } while ( n > 1.0  ||  n == 0 );
        n = rad / sqrt(n);
        proj[0] = n * x + cen;
        proj[1] = n * y;
        proj[2] = n * z;
        return rad;
    }
}


real SpaceTwinSpheres::projectC(real rad, const real point[], real proj[])
{
    real n = point[1] * point[1] + point[2] * point[2];
    
    if ( n > 0 )
    {
        n = rad / sqrt(n);
        proj[0] = 0;
        proj[1] = n * point[1];
        proj[2] = n * point[2];
    }
    else
    {
        //select a random point on the surface
        real y, z;
        do {
            y = RNG.sreal();
            z = RNG.sreal();
            n = y*y + z*z;
        } while ( n > 1.0  ||  n == 0 );
        n = rad / sqrt(n);
        proj[0] = 0;
        proj[1] = n * y;
        proj[2] = n * z;
    }
    return sqrt( point[0] * point[0]
                + (point[1]-proj[1]) * (point[1]-proj[1])
                + (point[2]-proj[2]) * (point[2]-proj[2]) );
}


#endif


void SpaceTwinSpheres::project( const real point[], real proj[] ) const
{
    real p[DIM];
    
    real dC = projectC(neck, point, proj);
    
    real dR = projectS(radiusR, cenR, point, p);
    if ( p[0] > 0  &&  dR < dC )
    {
        for ( int d = 0; d < DIM; ++d )
            proj[d] = p[d];
        dC = dR;
    }
    
    real dL = projectS(radiusL, cenL, point, p);
    if ( p[0] < 0  &&  dL < dC )
    {
        for ( int d = 0; d < DIM; ++d )
            proj[d] = p[d];
    }
}

//------------------------------------------------------------------------------
void SpaceTwinSpheres::setInteractions(Meca& meca, FiberSet const& fibers) const
{
    ABORT_NOW("unfinished");
    Vector dir(1,0,0);
    
    for ( Fiber * fib=fibers.first(); fib; fib=fib->next() )
    {
        for ( unsigned s = 0; s < fib->nbSegments() ; ++s )
        {
            real abs = fib->planarIntersect(s, dir, 0);
            if ( 0 <= abs  &&  abs < 1 )
                ;// meca.addLongPointClampYZ(PointInterpolated(seg, abs), Vector(0,0,0), neck, 100)
        }
    }
}

void SpaceTwinSpheres::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    ABORT_NOW("unfinished");
    if ( pos.XX > 0 )
        meca.addLongPointClamp(pos, pe, Vector(cenR,0,0), radiusR, stiff);
    else
        meca.addLongPointClamp(pos, pe, Vector(cenL,0,0), radiusL, stiff);
}


void SpaceTwinSpheres::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    ABORT_NOW("unfinished");
    if ( pos.XX > 0 )
    {
        if ( radiusR > rad )
            meca.addLongPointClamp(pos, pe, Vector(cenR,0,0), radiusR-rad, stiff);
        else
        {
            meca.addPointClamp(pe, Vector(cenR,0,0), stiff);
            std::cerr << "object is too big to fit in SpaceTwinSpheres\n";
        }
    }
    else
    {
        if ( radiusL > rad )
            meca.addLongPointClamp(pos, pe, Vector(cenL,0,0), radiusL-rad, stiff);
        else
        {
            meca.addPointClamp(pe, Vector(cenL,0,0), stiff);
            std::cerr << "object is too big to fit in SpaceTwinSpheres\n";
        }
    }
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY

#include "glut.h"
#include "gle.h"

bool SpaceTwinSpheres::display() const
{
    GLfloat L = radiusL;
    GLfloat R = radiusR;

#if ( DIM <= 2 )
    
    const GLfloat da = M_PI / 360;

    glBegin(GL_LINE_LOOP);
    
    for ( GLfloat a = -M_PI; a <= M_PI; a += da )
        if ( cenR+R*cosf(a) > 0 )
            glVertex2f(cenR+R*cosf(a), R*sinf(a));
    
    for ( GLfloat a = 0; a <= 2*M_PI; a += da )
        if ( cenL+L*cosf(a) < 0 )
            glVertex2f(cenL+L*cosf(a), L*sinf(a));
    
    glEnd();

#else
    
    const GLenum glp = GL_CLIP_PLANE5;
    glEnable(glp);
    
    //right side:
    GLdouble planeX[] = { +1, 0, 0, 0 };
    glClipPlane(glp, planeX);

    glPushMatrix();
    glTranslated(cenR, 0, 0);
    glScalef(R, R, R);
    gle::gleSphere8();
    glPopMatrix();

    //left side:
    planeX[0] = -1;
    glClipPlane(glp, planeX);

    glPushMatrix();
    glTranslated(cenL, 0, 0);
    glScalef(L, L, L);
    gle::gleSphere8();
    glPopMatrix();
    glDisable(glp);

#endif
    
    return true;
}

#else

bool SpaceTwinSpheres::display() const
{
    return false;
}


#endif

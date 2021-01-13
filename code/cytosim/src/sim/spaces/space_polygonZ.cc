// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_polygonZ.h"
#include "point_exact.h"
#include "exceptions.h"
#include "glossary.h"
#include "polygon.h"
#include "random.h"
#include "meca.h"
#include <fstream>

extern Random RNG;


SpacePolygonZ::SpacePolygonZ(SpaceProp const* p, Glossary& opt)
: Space(p)
{
    mVolume = 0;
    
    if ( DIM != 3 )
        throw InvalidParameter("polygonZ is only usable in 3D");
    
    mPoly.read(p->shape_spec);
    
#if ( DIM > 2 )
    Vector vec;
    if ( opt.set(vec, "translate") )
        mPoly.translate(vec.XX, vec.YY);

    real len;
    if ( opt.set(len, "inflate") && len > 0 )
        mPoly.inflate(len);
#endif
    
    resize();
}


SpacePolygonZ::~SpacePolygonZ()
{
}


/**
 The volume is estimated with a monte-carlo approach,
 but taking into account that the volume is axisymmetric.
 The result should be more precise than Space::estimateVolume()
 */
real SpacePolygonZ::estimateVolumeZ(unsigned long cnt) const
{
    real box[4];
    mPoly.find_extremes(box);
    
    real H = box[3] - box[2];

    real Ri = box[0];
    real Ro = box[1];
    if ( Ri < 0 )
    {
        box[0] = 0;
        Ri = 0;
    }
    real W = Ro - Ri;
    
    real in = 0, out = 0;
    for ( unsigned long i = 0; i < cnt; ++i )
    {
        real x = box[0] + W * RNG.preal();
        real y = box[2] + H * RNG.preal();
        if ( mPoly.inside(x, y, 1) )
            in  += x;
        else
            out += x;
    }
    
    return M_PI * H * ( Ro*Ro - Ri*Ri ) * in / ( in + out );
}


/**
 recalculate bounding box, volume
 and points offsets that are used to project
 */
void SpacePolygonZ::resize()
{
    if ( mPoly.surface() < 0 )
    {
        //std::clog << "flipping clockwise polygon `" << prop->shape_spec << "'" << std::endl;
        mPoly.flip();
    }

    if ( mPoly.complete(REAL_EPSILON) )
        throw InvalidParameter("unfit polygon: consecutive points may overlap");

    real box[4];
    mPoly.find_extremes(box);
    mInf.set(-box[1],-box[1], box[2]);
    mSup.set( box[1], box[1], box[3]);

    mVolume = estimateVolumeZ(1<<17);
}


bool SpacePolygonZ::inside( const real w[] ) const
{
    real R = sqrt( w[0]*w[0]+ w[1]*w[1] );

    return mPoly.inside(R, w[2], 1);
}


void SpacePolygonZ::project( const real w[], real p[] ) const
{
    real P, R = sqrt( w[0]*w[0]+ w[1]*w[1] );
    int hit;

    mPoly.project(R, w[2], P, p[2], hit);
    
    p[0] = w[0] * P / R;
    p[1] = w[1] * P / R;
}


/**
 The current procedure tests the model-points of fibers against the segments of the polygon.
 This fails for non-convex polygon since the re-entrant corners can intersect the fibers.
 
 @todo Also project re-entrant polygon corners on the segments of the Fiber.
 */
void SpacePolygonZ::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    //Space::setInteraction(pos, pe, meca, stiff); return;
#if ( DIM == 3 )
    real P, R = sqrt( pos.XX * pos.XX + pos.YY * pos.YY );
    int hit;
    
    Vector prj;

    int edg = mPoly.project(R, pos.ZZ, P, prj.ZZ, hit);
    real nX = -mPoly.pts_[hit].dy;
    real nY =  mPoly.pts_[hit].dx;

    prj.XX = pos.XX * P / R;
    prj.YY = pos.YY * P / R;
    
    if ( edg )
    {
        Vector dir( nX * pos.XX / R, nX * pos.YY / R, nY );
        meca.addPlaneClamp(pe, prj, dir, stiff);
    }
    else
    {
        Vector dir( -pos.YY / R, pos.XX / R, 0 );
        meca.addLineClamp(pe, prj, dir, stiff);
    }
#endif
}
                       

void SpacePolygonZ::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    //setInteraction(pos, pe, meca, stiff);
    std::cerr << "unfinished SpacePolygonZ::setInteractions(with radius)\n";
}


void SpacePolygonZ::setInteractions(Meca & meca, FiberSet const& fibers) const
{
    /// @todo add interactions between fibers and reentrant corners!
#if ( 0 )
    real stiffness = 0;
    Vector dir(0,0,1);
    
    for (Fiber * fib=fibers.first(); fib; fib=fib->next() )
    {
        for ( unsigned s = 0; s < fib->nbSegments() ; ++s )
        {
            //project on point on segment
            if ( 0 <= abs  &&  abs < 1 )
                ;// meca.addLongPointClampYZ(PointInterpolated(seg, abs), Vector(0,0,0), neck, 100)
        }
    }
#endif
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

void SpacePolygonZ::displayZ(bool rings) const
{
    const unsigned npts = mPoly.nbPoints();
    Polygon::Point2D const* pts = mPoly.pts_;

    const size_t FIN = 8 * gle::finesse;
    GLfloat c[FIN+1], s[FIN+1];
    for ( unsigned ii = 0; ii <= FIN; ++ii )
    {
        GLfloat ang = ii * 2 * M_PI / (GLfloat) FIN;
        c[ii] = cosf(ang);
        s[ii] = sinf(ang);
    }
    
    //display surface
    for ( unsigned n=1; n <= npts; n++ )
    {
        // do not display special edges
        if ( pts[n-1].color )
            continue;
        
        real R1 = pts[n-1].xx;
        real Z1 = pts[n-1].yy;
        real R2 = pts[n].xx;
        real Z2 = pts[n].yy;
        
        real nX =  pts[n-1].dy;
        real nY = -pts[n-1].dx;
        
        if ( R1 >= 0 && R2 >= 0 )
        {
            glBegin(GL_TRIANGLE_STRIP);
            for ( int ii = 0; ii <= FIN; ++ii )
            {
                glNormal3f( nX*c[ii], nX*s[ii], nY );
                glVertex3f( R2*c[ii], R2*s[ii], Z2 );
                glVertex3f( R1*c[ii], R1*s[ii], Z1 );
            }
            glEnd();
        }
    }
    
    //display rings around:
    if ( rings )
    {
        glLineWidth(0.5);
        for ( unsigned n=0; n < npts; n++ )
        {
            real R = pts[n].xx;
            real Z = pts[n].yy;
            if ( R > 0 )
            {
                glBegin(GL_LINE_LOOP);
                for ( int ii = 0; ii <= FIN; ++ii )
                    gle::gleVertex( R*c[ii], R*s[ii], Z );
                glEnd();
            }
        }
    }
}


bool SpacePolygonZ::display() const
{
    displayZ(0);
    return true;
}

#else

bool SpacePolygonZ::display() const
{
    return false;
}

#endif

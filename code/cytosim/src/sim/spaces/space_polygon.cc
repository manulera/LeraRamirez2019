// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_polygon.h"
#include "point_exact.h"
#include "exceptions.h"
#include "glossary.h"
#include "polygon.h"
#include "meca.h"
#include <fstream>


SpacePolygon::SpacePolygon(SpaceProp const* p, Glossary& opt)
: Space(p)
{
    mVolume = 0;
    
    if ( DIM == 1 )
        throw InvalidParameter("polygon is not usable in 1D");
    
    mPoly.read(p->shape_spec);
    
    if (  mPoly.surface() < 0 )
    {
        std::clog << "flipping clockwise polygon `" << prop->shape_spec << "'" << std::endl;
        mPoly.flip();
    }

    real x;
    if ( opt.set(x, "inflate") )
        mPoly.inflate(x);
    
    resize();
}


SpacePolygon::~SpacePolygon()
{
}


/**
 recalculate bounding box, volume
 and points offsets that are used to project
 */
void SpacePolygon::resize()
{
    //in 3D one dimension should be specified: height
    checkLengths(DIM==3, true);

    mHeight = mLength[0];
    mVolume = mPoly.surface();
    assert_true( mPoly.surface() > 0 );
    
#if ( DIM == 3 )
    // total height = half_height
    mVolume *= 2 * mHeight;
#endif
    
    if ( mPoly.complete(REAL_EPSILON) )
        throw InvalidParameter("unfit polygon: consecutive points may overlap");

    real box[4];
    mPoly.find_extremes(box);
    mInf.set(box[0], box[2], 0);
    mSup.set(box[1], box[3], 0);
}


bool SpacePolygon::inside( const real w[] ) const
{
#if ( DIM == 3 )
    if ( w[2] < -mHeight  ||  w[2] > mHeight )
        return false;
#endif
    return mPoly.inside(w[0], w[1], 1);
}


void SpacePolygon::project( const real w[], real p[] ) const
{    
#if ( DIM == 1 )
    
    p[0] = w[0];
    
#elif ( DIM == 2 )
    
    int hit;
    mPoly.project(w[0], w[1], p[0], p[1], hit);
    
    
#elif ( DIM == 3 )
    
    if ( fabs(w[2]) > mHeight )
    {
        if ( mPoly.inside(w[0], w[1], 1) )
        {
            // too high or too low in the Z axis, but inside XY
            p[0] = w[0];
            p[1] = w[1];
        }
        else
        {
            // outside in Z and XY
            int hit;
            mPoly.project(w[0], w[1], p[0], p[1], hit);
        }
        p[2] = (w[2]>0) ? mHeight : -mHeight;
    }
    else
    {
        int hit;
        mPoly.project(w[0], w[1], p[0], p[1], hit);
        
        if ( mPoly.inside(w[0], w[1], 1) )
        {
            // inside in the Z axis and the XY polygon: compare distances
            
            real hdis = (w[0]-p[0])*(w[0]-p[0]) + (w[1]-p[1])*(w[1]-p[1]);
            
            //we are inside in both the Z and XY sense, we compare the distances
            //to the top/bottom planes, and to the sides of the polygon
            //calculate the distance to the top/bottom planes:
            real vdis = mHeight - fabs(w[2]);
            if ( vdis * vdis < hdis )
            {
                p[0] = w[0];
                p[1] = w[1];
                p[2] = (w[2]>0) ? mHeight : -mHeight;
                return;
            }
            else
            {
                p[2] = w[2];
            }
        }
        else 
        {
            // outsize in XY, inside in Z
            p[2] = w[2];
        }
    }
    
#endif
}


/**
 The current procedure tests the model-points of fibers against the segments of the polygon.
 This fails for non-convext polygon since the re-entrant corners can intersect the fibers.
 
 @todo Also project re-entrant polygon corners on the segments of the Fiber.
 */
void SpacePolygon::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{    
#if ( DIM > 1 )
    Matrix::index_type inx = DIM * pe.matIndex();
    
    int hit;
    real pX, pY;
    int edg = mPoly.project(pos.XX, pos.YY, pX, pY, hit);
    real nX = -mPoly.pts_[hit].dy;
    real nY =  mPoly.pts_[hit].dx;
    
#if ( DIM == 3 )
    
    if ( pos.ZZ >= mHeight )
    {
        meca.mC(inx+2, inx+2) -= stiff;
        meca.base(inx+2)      += stiff * mHeight;
        if ( mPoly.inside(pos.XX, pos.YY, 1) )
            return;
    }
    else if ( pos.ZZ <= -mHeight )
    {
        meca.mC(inx+2, inx+2) -= stiff;
        meca.base(inx+2)      -= stiff * mHeight;
        if ( mPoly.inside(pos.XX, pos.YY, 1) )
            return;
    }
    else
    {
        // Compare distance to top/bottom plate, and distance to edge in XY plane
        real vdis = mHeight - fabs(pos.ZZ);
        real hdis = (pos.XX-pX)*(pos.XX-pX) + (pos.YY-pY)*(pos.YY-pY);
        
        if ( vdis * vdis < hdis  &&  mPoly.inside(pos.XX, pos.YY, 1) )
        {
            if ( pos.ZZ >= 0 )
            {
                meca.mC(inx+2, inx+2) -= stiff;
                meca.base(inx+2)      += stiff * mHeight;
            }
            else
            {
                meca.mC(inx+2, inx+2) -= stiff;
                meca.base(inx+2)      -= stiff * mHeight;
            }
            return;
        }
    }

#endif

    if ( edg )
    {
        // projection on an edge of normal (nX, nY) already normalized
        const real pr = ( pX * nX + pY * nY ) * stiff;
        
        meca.mC(inx  , inx  ) -= nX * nX * stiff;
        meca.mC(inx  , inx+1) -= nX * nY * stiff;
        meca.mC(inx+1, inx+1) -= nY * nY * stiff;
        
        meca.base(inx  )  += nX * pr;
        meca.base(inx+1)  += nY * pr;
    }
    else
    {
        // projection on a vertex:
#if ( DIM == 2 )
        meca.mB(pe.matIndex(), pe.matIndex()) -= stiff;
#elif ( DIM == 3 )
        meca.mC(inx,   inx  ) -= stiff;
        meca.mC(inx+1, inx+1) -= stiff;
#endif
        meca.base(inx  )  += stiff * pX;
        meca.base(inx+1)  += stiff * pY;
    }
#endif
}


void SpacePolygon::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    //setInteraction(pos, pe, meca, stiff);
    std::cerr << "unfinished SpacePolygon::setInteractions(with radius)\n";
}

#include "fiber_locus.h"
#include "fiber_set.h"

void SpacePolygon::setInteractions(Meca & meca, FiberSet const& fibers) const
{
#if ( 0 )
    /// WORK IN PROGRESS
    Polygon::Point2D const* pts = mPoly.pts_;
    const int n_pik = 2;
    const int inx[n_pik] = { 0, 100 };
    Vector pik[n_pik];
    
    for ( int i = 0; i < n_pik; ++i )
        pik[i].set(pts[inx[i]].xx, pts[inx[i]].yy, 0);
    
    for ( Fiber * fib=fibers.first(); fib; fib=fib->next() )
    {
        real ls = fib->segmentation();
        for ( unsigned seg = 0; seg < fib->nbSegments() ; ++seg )
        {
            FiberLocus loc(fib, seg);
            for ( int i = 0; i < n_pik; ++i )
            {
                real dis;
                real abs = loc.projectPoint(pik[i], abs, dis);
                if ( 0 <= abs  &&  abs < ls )
                {
                    if ( !inside(loc.pos(abs)) || !inside(loc.pos1()) || !inside(loc.pos2()) )
                        meca.addPointClamp(PointInterpolated(loc, abs), pik[i], 100);
                }
            }
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

bool SpacePolygon::display() const
{
    const unsigned npts = mPoly.nbPoints();
    Polygon::Point2D const* pts = mPoly.pts_;

#if ( DIM == 2 )
    
    glEnable(GL_STENCIL_TEST);
    glClearStencil(1);
    glClear(GL_STENCIL_BUFFER_BIT);
    glStencilFunc(GL_EQUAL, 1, ~0);
    glStencilOp(GL_KEEP, GL_ZERO, GL_ZERO);
    
    //display points
    GLfloat s = 1;
    glGetFloatv(GL_LINE_WIDTH, &s);
    glPointSize(s);
    glBegin(GL_POINTS);
    for ( unsigned n=0; n < npts; ++n )
        gle::gleVertex(pts[n].xx, pts[n].yy);
    glEnd();

    //display polygon
    glBegin(GL_LINE_LOOP);
    for ( unsigned n=0; n < npts; ++n )
        gle::gleVertex(pts[n].xx, pts[n].yy);
    glEnd();
    
    glClear(GL_STENCIL_BUFFER_BIT);
    glDisable(GL_STENCIL_TEST);

#elif ( DIM == 3 )
    
    // display bottom
    glLineWidth(2);
    glBegin(GL_LINE_LOOP);
    for ( unsigned n=0; n < npts; ++n )
        gle::gleVertex(pts[n].xx, pts[n].yy, -mHeight);
    glEnd();
    
    // display top
    glBegin(GL_LINE_LOOP);
    for ( unsigned n=npts; n > 0; --n )
        gle::gleVertex(pts[n].xx, pts[n].yy,  mHeight);
    glEnd();
    
    // display sides
    real Z = mHeight;
    glBegin(GL_TRIANGLE_STRIP);
    for ( unsigned n=0; n <= npts; ++n )
    {
        gle::gleVertex(pts[n].xx, pts[n].yy, Z);
        gle::gleVertex(pts[n].xx, pts[n].yy,-Z);
    }
    glEnd();
    
#endif
    
#if ( 0 )
    // display points:
    glColor3f(1,1,1);
    glPointSize(3);
    glBegin(GL_POINTS);
    for ( unsigned n=0; n < npts; ++n )
        gle::gleVertex(pts[n].xx, pts[n].yy);
    glEnd();
#endif
#if ( 0 )
    // indicate index of each point:
    char tmp[8];
    for ( unsigned n=0; n < npts; ++n )
    {
        snprintf(tmp, sizeof(tmp), "%i", n);
        Vector p(pts[n].xx, pts[n].yy, mHeight);
        gle::gleDrawText(p, tmp, 0);
    }
#endif
    

    return true;
}

#else

bool SpacePolygon::display() const
{
    return false;
}

#endif

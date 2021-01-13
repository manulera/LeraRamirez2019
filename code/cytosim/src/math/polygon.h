// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef POLYGON_H
#define POLYGON_H

#include "assert_macro.h"
#include <fstream>

#ifndef REAL
    #include "real.h"
#endif

/// Basic data structures and function associated with closed polygons
class Polygon
{
public:

    /// holds the coordinates of a point in 2D + info about the segment
    struct Point2D 
    {
        real xx, yy;   ///< coordinates of point
        int  color;    ///< this is used to indicate the type of edge
        real dx, dy;   ///< normalized direction to next point
        real len;      ///< distance to next point
        
        /// constructor
        Point2D() {}
        
        /// set coordinates of point
        Point2D(real sx, real sy)
        {
            xx = sx;
            yy = sy;
            color = 0;
        }
        
        /// test if given point overlap with *this
        bool operator == (const Point2D& p)
        {
            return ( xx == p.xx  &&  yy == p.yy );
        }
    };
    
    /// list of points. The array is allocated to hold index = 1+npts_
    Point2D* pts_;
    
    /// number of points
    unsigned npts_;
    
public:
    
    /// constructor
    Polygon();
    
    /// destructor
    ~Polygon();
    
    /// number of points
    unsigned nbPoints() const { return npts_; }
    
    /// set number of points and allocate memory
    void     resize(unsigned s);
    
    /// return copy of point at index `inx`
    Point2D  point(unsigned inx) { assert_true( inx < npts_ ); return pts_[inx]; }

    /// set coordinates of point at index `inx`:
    void     set(unsigned inx, real x, real y);
    
    /// subfunction
    static unsigned read(std::istream&, Point2D *pts, unsigned pts_size);
    
    /// read polygon from stream
    void     read(std::istream&);
    
    /// read polygon from file
    void     read(std::string const&);
    
    /// write a polygon to stream
    void     write(std::ostream&) const;

    /// flip the order of the points
    void     flip();
    
    /// move all points by given amount
    void     translate(real dx, real dy);
    
    /// move all points by given amount
    void     scale(real sx, real sy);
    
    /// move all segments sideways, to uniformly increase the surface of polygon
    void     inflate(real eps);
    
    /// calculate the offsets necessary for the other functions. Return 0 if OK
    int      complete(real epsilon);
    
    /// tell if a point is inside a polygon
    int      inside(real x, real y, int edge, real threshold = REAL_EPSILON) const;
    
    /// calculate the projection (pX, pY) of the point (x,y) on a polygon
    int      project(real x, real y, real& pX, real& pY, int& hit) const;
    
    /// calculate the bounding box [xmin, xmax, ymin, ymax] of a polygon
    void     find_extremes(real box[4]) const;
    
    /// calculate the surface of a polygon
    real     surface() const;
    
    /// printout
    void     dump(std::ostream&) const;
};

#endif


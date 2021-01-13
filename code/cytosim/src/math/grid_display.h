// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created 09/03/2015 by Francois Nedelec

#ifndef GRID_DISPLAY_H
#define GRID_DISPLAY_H


#include "opengl.h"
#include "vector.h"

/// display the edges of a 1D grid using OpenGL
/**
 This uses the current OpenGL color and line width.
 */
template <typename INDEX>
void drawEdges(GridBase<1, INDEX> const& grid)
{
    const float u =  0.5;
    const float d = -0.5;
    glBegin(GL_LINES);
    for ( real ix = 0; ix <= grid.dim(0); ++ix )
    {
        float x = grid.position(0, ix);
        glVertex2f(x, d);
        glVertex2f(x, u);
    }
    glEnd();
}


/// display the edges of a 2D grid using OpenGL
/**
 This uses the current OpenGL color and line width.
 */
template <typename INDEX>
void drawEdges(GridBase<2, INDEX> const& grid)
{
    float i = grid.inf(0);
    float s = grid.sup(0);
    glBegin(GL_LINES);
    for ( float iy = 0; iy <= grid.dim(1); ++iy )
    {
        real y = grid.position(1, iy);
        glVertex2f(i, y);
        glVertex2f(s, y);
    }
    glEnd();
    
    i = grid.inf(1);
    s = grid.sup(1);
    glBegin(GL_LINES);
    for ( float ix = 0; ix <= grid.dim(0); ++ix )
    {
        float x = grid.position(0, ix);
        glVertex2f(x, i);
        glVertex2f(x, s);
    }
    glEnd();
}


/// display the edges of a 3D grid using OpenGL
/**
 This uses the current OpenGL color and line width.
 */
template <typename INDEX>
void drawEdges(GridBase<3, INDEX> const& grid)
{
    float i = grid.inf(0);
    float s = grid.sup(0);
    glBegin(GL_LINES);
    for ( float iy = 0; iy <= grid.dim(1); ++iy )
    for ( float iz = 0; iz <= grid.dim(2); ++iz )
    {
        float y = grid.position(1, iy);
        float z = grid.position(2, iz);
        glVertex3f(i, y, z);
        glVertex3f(s, y, z);
    }
    glEnd();
    
    i = grid.inf(1);
    s = grid.sup(1);
    glBegin(GL_LINES);
    for ( float ix = 0; ix <= grid.dim(0); ++ix )
    for ( float iz = 0; iz <= grid.dim(2); ++iz )
    {
        float x = grid.position(0, ix);
        float z = grid.position(2, iz);
        glVertex3f(x, i, z);
        glVertex3f(x, s, z);
    }
    glEnd();
    
    i = grid.inf(2);
    s = grid.sup(2);
    glBegin(GL_LINES);
    for ( float ix = 0; ix <= grid.dim(0); ++ix )
    for ( float iy = 0; iy <= grid.dim(1); ++iy )
    {
        float x = grid.position(0, ix);
        float y = grid.position(1, iy);
        glVertex3f(x, y, i);
        glVertex3f(x, y, s);
    }
    glEnd();
}


//------------------------------------------------------------------------------
#pragma mark -


/// display the values stored in the cells of a 1D grid using OpenGL
/**
 OpenGL color is to be specified by the function `set_color`:
 bool set_color(CELL const&, Vector2 const&);
 The return value of this function is used to enable the display of
 individual cells: true = display; false = skip.
 */
template <typename CELL, typename INDEX>
void drawValues(Grid<CELL, 1, INDEX> const& grid,
                bool set_color(CELL const&, Vector1 const&))
{
    float d = grid.cellWidth(0);
    float e = 2;
    
    // paint all cells one by one
    for ( INDEX c = 0; c < grid.dim(0); ++c )
    {
        float x = grid.position(0, c);
        if ( set_color(grid[c], Vector1(x,0)) )
        {
            glBegin(GL_TRIANGLE_STRIP);
            glVertex2f(x  , -e);
            glVertex2f(x+d, -e);
            glVertex2f(x  ,  e);
            glVertex2f(x+d,  e);
            glEnd();
        }
    }
}


/// display the values stored in the cells of a 2D grid using OpenGL
/**
 OpenGL color is to be specified by the function `set_color`:
 bool set_color(CELL const&, Vector2 const&);
 The return value of this function is used to enable the display of 
 individual cells: true = display; false = skip.
 */
template <typename CELL, typename INDEX>
void drawValues(Grid<CELL, 2, INDEX> const& grid,
                bool set_color(CELL const&, Vector2 const&))
{
    float d = 0.5 * grid.cellWidth(0);
    float e = 0.5 * grid.cellWidth(1);
    
    // paint all cells one by one
    for ( INDEX c = 0; c < grid.nbCells(); ++c )
    {
        Vector2 w;
        grid.setPositionFromIndex(w, c, 0.5);
        if ( set_color(grid[c], w) )
        {
            glBegin(GL_TRIANGLE_STRIP);
            glVertex2f(w.XX-d, w.YY-e);
            glVertex2f(w.XX+d, w.YY-e);
            glVertex2f(w.XX-d, w.YY+e);
            glVertex2f(w.XX+d, w.YY+e);
            glEnd();
        }
    }
}



/// display the slice of a 3D grid in a plane parallel to XY at `Z = z_pos`
/**
 OpenGL color is to be specified by the function `set_color`:
 bool set_color(CELL const&, Vector2 const&);
 The return value of this function is used to enable the display of
 individual cells: true = display; false = skip.
 */
template <typename CELL, typename INDEX>
void drawValues(Grid<CELL, 3, INDEX> const& grid,
                bool set_color(CELL const&, Vector3 const&),
                real z_pos = 0)
{
    assert_true(grid.hasCells());
    
    real d = 0.5 * grid.cellWidth(0);
    real e = 0.5 * grid.cellWidth(1);
    
    INDEX z = grid.index(2, z_pos);
    
    for ( INDEX y = 0; y < grid.dim(1); ++y )
    for ( INDEX x = 0; x < grid.dim(0); ++x )
    {
        Vector3 w(grid.position(0, x+0.5), grid.position(1, y+0.5), z_pos);
        if ( set_color(grid.icell3D(x,y,z), w) )
        {
            glBegin(GL_TRIANGLE_STRIP);
            glVertex3f(w.XX-d, w.YY-e, z_pos);
            glVertex3f(w.XX+d, w.YY-e, z_pos);
            glVertex3f(w.XX-d, w.YY+e, z_pos);
            glVertex3f(w.XX+d, w.YY+e, z_pos);
            glEnd();
        }
    }
}


/// display the slice of a 3D grid in a plane parallel to Y: `Y=pos`
/**
 OpenGL color is to be specified by the function `set_color`:
 bool set_color(CELL const&, Vector2 const&);
 The return value of this function is used to enable the display of
 individual cells: true = display; false = skip.
 */
template <typename CELL, typename INDEX>
void drawValuesXZ(Grid<CELL, 3, INDEX> const& grid,
                  bool set_color(CELL const&, Vector3 const&),
                  real pos)
{
    assert_true(grid.hasCells());
    
    real d = 0.5 * grid.cellWidth(0);
    real e = 0.5 * grid.cellWidth(2);
    
    INDEX y = grid.index(1, pos);
    
    for ( INDEX z = 0; z < grid.dim(2); ++z )
    for ( INDEX x = 0; x < grid.dim(0); ++x )
    {
        Vector3 w(grid.position(0, x+0.5), pos, grid.position(2, y+0.5));
        if ( set_color(grid.icell3D(x,y,z), w) )
        {
            glBegin(GL_TRIANGLE_STRIP);
            glVertex3f(w.XX-d, pos, w.ZZ-e);
            glVertex3f(w.XX+d, pos, w.ZZ-e);
            glVertex3f(w.XX-d, pos, w.ZZ+e);
            glVertex3f(w.XX+d, pos, w.ZZ+e);
            glEnd();
        }
    }
}



// display the slice of a 3D grid in a plane parallel to X: `X=pos`
/**
 OpenGL color is to be specified by the function `set_color`:
 bool set_color(CELL const&, Vector2 const&);
 The return value of this function is used to enable the display of
 individual cells: true = display; false = skip.
 */
template <typename CELL, typename INDEX>
void drawValuesYZ(Grid<CELL, 3, INDEX> const& grid,
                  bool set_color(CELL const&, Vector3 const&),
                  real pos)
{
    assert_true(grid.hasCells());
    
    real d = 0.5 * grid.cellWidth(0);
    real e = 0.5 * grid.cellWidth(2);
    
    INDEX x = grid.index(0, pos);
    
    for ( INDEX z = 0; z < grid.dim(2); ++z )
    for ( INDEX y = 0; y < grid.dim(1); ++y )
    {
        Vector3 w(pos, grid.position(1, y+0.5), grid.position(2, z+0.5));
        if ( set_color(grid.icell3D(x,y,z), w) )
        {
            glBegin(GL_TRIANGLE_STRIP);
            glVertex3f(pos, w.YY-d, w.ZZ-e);
            glVertex3f(pos, w.YY+d, w.ZZ-e);
            glVertex3f(pos, w.YY-d, w.ZZ+e);
            glVertex3f(pos, w.YY+d, w.ZZ+e);
            glEnd();
        }
    }
}




/// display a slice of the field in a plane perpendicular to 'dir'
/**
 OpenGL color is to be specified by the function `set_color`:
 bool set_color(CELL const&, Vector2 const&);
 The return value of this function is ignored.
 */
template <typename CELL, typename INDEX>
void drawValues(Grid<CELL, 3, INDEX> const& grid,
                bool set_color(CELL const&, Vector3 const&),
                Vector3 const& dir, real z_pos)
{
    assert_true(grid.hasCells());
    
    // this defines the finesse of the triangular mesh:
    real n = 0.2 * grid.minimumWidth(1);
    int m = grid.radius() / n;
    
    Vector3 dx, dy;
    dir.orthonormal(dx, dy);
    dx *= n;
    dy *= n;
    
    Vector3 dh = dy * cos(M_PI/6);
    Vector3 w, a;
    
    for ( int y = -m; y <= m; y+=2 )
    {
        a = y * dh + z_pos * dir;
        glBegin(GL_TRIANGLE_STRIP);
        for ( int x = -m; x <= m; ++x )
        {
            w = a + x * dx;
            set_color(grid.interpolate3D(w), w);
            gle::gleVertex(w);
            
            w = a + ( x + 0.5 ) * dx + dh;
            set_color(grid.interpolate3D(w), w);
            gle::gleVertex(w);
        }
        glEnd();
        glBegin(GL_TRIANGLE_STRIP);
        for ( int x = -m; x <= m; ++x )
        {
            w = a + x * dx + dh + dh;
            set_color(grid.interpolate3D(w), w);
            gle::gleVertex(w);
            
            w = a + ( x + 0.5 ) * dx + dh;
            set_color(grid.interpolate3D(w), w);
            gle::gleVertex(w);
        }
        glEnd();
   }
}


#endif



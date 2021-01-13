// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Francois Nedelec; Created 07/03/2015. nedelec@embl.de

#ifndef GRID_BASE_H
#define GRID_BASE_H

#include "assert_macro.h"
#include "exceptions.h"
#include <cstdio>
#include <cmath>
#include "real.h"

///\def flag to disable support for periodic boundaries conditions
/**
Periodic boundary conditions are supported using a function-pointer.
 If keyword NO_PERIODIC_SUPPORT is defined however, inlined functions are 
 used instead, which might be faster, but Periodic boundary are then not supported.
 */
//#define NO_PERIODIC_SUPPORT


///Divides a rectangle of dimensionality ORD into regular voxels
/** 
Grid<int ORD, typename INDEX> creates a regular lattice over a rectangular
region of space of dimensionality ORD, initialized by setDimensions().

Functions are provided to convert from the real space coordinates (of type real)
into an index (of type INDEX) usable to access a one-dimensional C-array.
The cells are ordered successively, the first dimension (X) varying the fastest
i.e. cell[ii+1] will in most cases be located on the right of cell[ii], although
if cell[ii] is on the right edge, then cell[ii+1] is on the symmetric edge. 

\par Access:

Cells can be accessed in three ways:
 - Position:      a set of real       operator()( real[] ), or operator(real, real, real)
 - Index:         one integer         operator[](int index)
 - Coordinates:   a set of integer    function cell(int[]), or cell(int,int,int)
.
Valid indices are [0...nbCells()-1], where nbCells() is calculated by setDimensions().
If a position lies outside the rectangular region where the grid is defined,
index(real[]) returns the index of the closest voxel.

Functions to convert between the three types are provided:
 - index()
 - pack()
 - setCoordinatesFromIndex(),
 - setCoordinatesFromPosition()
 - setPositionFromCoordinates()
 - setPositionFromIndex()
.

\par Indices:

The grid is initialized by setDimensions(inf, sup, nbCells), which calculates:
  - cWidth[d] = ( sup[d] - inf[d] ) / nbCells[d], for d in [0, ORD[

The coordinates of a cell at position pos[] are:
  - c[d] = int(  ( pos[d] - inf[d] ) / cWidth[d] )

and its index is
  - with ORD==1: index = c[0]
  - with ORD==2: index = c[0] + nbcells[0] * c[1]
  - with ORD==3: index = c[0] + nbcells[0] * ( c[1] + nbcells[1] * c[2] )
  - etc.
.
    
For a 4x4 2D grid, the index are like this:
@code
12  13  14  15
8    9  10  11
4    5   6   7
0    1   2   3
@endcode

\par Neighborhood:

The class also provides information on which cells surround each cell:
 - createSquareRegions(range) calculates square regions of size range
   ( range==1 gives nearest neighbors ).
 - createRoundRegions(range) calculates round regions of size range
 - createSideRegions(range)
.
After calling one of the above function, getRegion(offsets, index) will set 'offsets'
to point to an array of 'index offsets' for the cell referred by 'index'.
A zero offset value (0) is always first in the list and refers to self.
In the example above:
    - for index = 0 it would return { 0 1 4 5 }
    - for index = 5 it would return { 0 -1 1 -5 -4 -3 3 4 5 }
.
You obtain the cell-indices of the neighboring cells by adding offsets[n] to 'index':
Example:
@code
    CELL * cell = & myGrid.icell(indx);
    nb_neighbors = myGrid.getRegion(region, indx);
    for ( int n = 1; n < nb_neighbors; ++n ) 
    {
        Cell & neighbor = cell[region[n]];
        ...
    }
@endcode
*/

///\todo add Grid<> copy constructor and copy assignment

template <int ORD, typename INDEX>
class GridBase
{
    
    /// Disabled copy constructor
    GridBase<ORD, INDEX>(GridBase<ORD, INDEX> const&);
    
    /// Disabled copy assignment
    GridBase<ORD, INDEX>& operator=(GridBase<ORD, INDEX> const&);

protected:
    
#ifdef NO_PERIODIC_SUPPORT
    
    /// return closest integer to c in the segment [ 0, s-1 ]
    inline INDEX imageI(const int& s, const int& c) const
    {
        return c <= 0 ? 0 : ( c >= s ? s-1 : c );
    }

    /// return f modulo s in [ 0, s [
    inline INDEX imageF(const int& s, const real& f) const
    {
        if ( f > 0 )
        {
            int c = (int) f;
            return ( c >= s ? s-1 : c );
        }
        return 0;
    }

#else
    
    /// return closest integer to c in the segment [ 0, s-1 ]
    static INDEX imageIB(int s, int c)
    {
        return c <= 0 ? 0 : ( c >= s ? s-1 : c );
    }

    /// return c modulo s in [ 0, s [
    static INDEX imageIP(int s, int c)
    {
        while ( c <  0 )  c += s;
        while ( c >= s )  c -= s;
        return c;
    }

    /// return closest integer to c in the segment [ 0, s-1 ]
    static INDEX imageFB(int s, real c)
    {
        return c <= 0 ? 0 : ( c >= s ? s-1 : (int)c );
    }
    
    /// return c modulo s in [ 0, s [
    static INDEX imageFP(int s, real c)
    {
        while ( c <  0 )  c += s;
        int i = (int)c;
        while ( i >= s )  i -= s;
        return i;
    }

    /// pointer to imageIB() or imageIP()
    INDEX (*imageI)(int, int);
    /// pointer to imageFB() or imageFP()
    INDEX (*imageF)(int, real);
    
#endif

public:
    
    /// the type for indices (=INDEX)
    typedef INDEX index_type;

protected:
    
    /// allocated size of array cells[]
    INDEX   allocated;
   
    /// Total number of cells in the map; size of cells[]
    INDEX   nCells;
    
    /// The number of cells in each dimension
    INDEX   gDim[ORD];
    
    /// Offset between two consecutive cells along each dimension
    INDEX   gStride[ORD];
    
    /// The position of the inferior edge (min) in each dimension
    real    gInf[ORD];
    
    /// The position of the superior edge (max) in each dimension
    real    gSup[ORD];
    
    /// cWidth[d] = ( gSup[d] - inf[d] ) / gDim[d]
    real    cWidth[ORD];
    
    /// cDelta[d] = 1.0 / cWidth[d]
    real    cDelta[ORD];
    
    /// The volume occupied by one cell
    real    cVolume;

    
//--------------------------------------------------------------------------
#pragma mark -
public:
    
    /// constructor
    GridBase() : gDim(), gInf(), gSup(), cWidth(), cDelta()
    {
        allocated   = 0;
        nCells      = 0;
        regionsEdge = 0;
        regions     = 0;
        cVolume     = 0;
#ifndef NO_PERIODIC_SUPPORT
        imageI      = imageIB;
        imageF      = imageFB;
#endif
    }
    
    /// Free memory
    void destroy()
    {
        deleteRegions();
    }
    
    /// Destructor
    virtual ~GridBase()
    {
        destroy();
    }
    
    
#ifdef NO_PERIODIC_SUPPORT
    
    /// true if boundary conditions are periodic
    bool periodic()
    {
        return false;
    }
    
    /// change boundary conditions
    void periodic(const bool p)
    {
        if ( p )
        {
            throw InvalidParameter("grid.h was compiled with NO_PERIODIC_SUPPORT");
        }
    }

#else
    
    /// true if boundary conditions are periodic
    bool periodic()
    {
        return imageI == imageIP;
    }
    
    /// change boundary conditions
    void periodic(const bool p)
    {
        if ( p )
        {
            imageI = imageIP;
            imageF = imageFP;
        }
        else
        {
            imageI = imageIB;
            imageF = imageFB;
        }
    }
    
#endif
    
    //--------------------------------------------------------------------------
    /// specifies the area covered by the Grid
    /**
     the edges of the area are specified in dimension `d` by 'infs[d]' and 'sups[d]',
     and the number of cells by 'nbcells[d]'.
     */
    void setDimensions(const real infs[ORD], real sups[ORD], const int nbcells[ORD])
    {
        nCells = 1;
        cVolume = 1;
        
        for ( int d = 0; d < ORD; ++d )
        {
            if ( nbcells[d] <= 0 )
                throw InvalidParameter("Cannot build a grid as nbcells[] is <= 0");
            
            if ( infs[d] > sups[d] )
            {
                if ( infs[d] > sups[d] + REAL_EPSILON )
                    throw InvalidParameter("Cannot build a grid as sup[] < inf[]");
                sups[d] = infs[d];
            }
            
            if ( infs[d] == sups[d] )
                throw InvalidParameter("Cannot build a grid as sup[] == inf[]");
            
            gStride[d] = nCells;
            nCells    *= nbcells[d];
            gDim[d]    = nbcells[d];
            gInf[d]    = infs[d];
            gSup[d]    = sups[d];
            cWidth[d]  = ( gSup[d] - gInf[d] ) / real( gDim[d] );
            cDelta[d]  = real( gDim[d] ) / ( gSup[d] - gInf[d] );
            cVolume   *= cWidth[d];
        }
    }
    
    ///true if setDimensions() was called
    bool hasDimensions()
    {
        return nCells > 0;
    }
    
    //--------------------------------------------------------------------------
#pragma mark -

    /// total number of cells in the map
    INDEX           nbCells()           const { return nCells; }
    INDEX           dim()               const { return nCells; }

    /// number of cells in dimensionality d
    INDEX           nbCells(int d)      const { return gDim[d]; }
    INDEX           dim(int d)          const { return gDim[d]; }
    
    /// offset to the next cell in the direction `d`
    INDEX           stride(int d)       const { return gStride[d]; }
    
    /// position of the inferior (left/bottom/etc) edge
    const real*     inf()               const { return gInf;    }
    real            inf(int d)          const { return gInf[d]; }
    
    /// position of the superior (right/top/etc) edge
    const real*     sup()               const { return gSup;    }
    real            sup(int d)          const { return gSup[d]; }
    
    /// the widths of a cell
    const real*     delta()             const { return cDelta;    }
    real            delta(int d)        const { return cDelta[d]; }
    
    const real*     cellWidth()         const { return cWidth;    }
    real            cellWidth(int d)    const { return cWidth[d]; }

    /// position in dimension `d`, of the cell of index `c`
    real            position(int d, real c) const { return gInf[d] + c * cWidth[d]; }
    
    /// index in dimension `d` corresponding to position `w`
    int             index(int d, real w) const { return ( w - gInf[d] ) * cDelta[d]; }
    
    /// the volume of a cell
    real            cellVolume()        const { return cVolume; }

    /// the length of the diagonal of a cell = sqrt( sum(cWidth[d]^2) )
    real diagonalLength() const
    {
        real res = cWidth[0] * cWidth[0];
        for ( int d = 1; d < ORD; ++d )
            res += cWidth[d] * cWidth[d];
        return sqrt( res );
    }
    
    /// the smallest cell width, along dimensions that have more than `min_size` cells
    real minimumWidth(unsigned int min_size) const
    {
        real res = 0;
        for ( int d = 0; d < ORD; ++d )
            if ( gDim[d] > min_size )
                res = cWidth[d];
        for ( int d = 0; d < ORD; ++d )
            if ( gDim[d] > min_size  &&  cWidth[d] < res )
                res = cWidth[d];
        return res;
    }
    
    /// radius of the minimal sphere placed in (0,0,0) that entirely covers all cells
    real radius() const
    {
        real res = 0;
        for ( int d = 0; d < ORD; ++d )
        {
            real m = gSup[d] > -gInf[d] ? gSup[d] : -gInf[d];
            res += m * m;
        }
        return sqrt(res);
    }

    //--------------------------------------------------------------------------
#pragma mark - Conversion

    /// checks if coordinates are inside the box
    bool inside(const int coord[ORD]) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            if ( coord[d] < 0 || (index_type)coord[d] >= gDim[d] )
                return false;
        }
        return true;
    }
    
    /// checks if point is inside the box
    bool inside(const real w[ORD]) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            if ( w[d] < gInf[d] || w[d] >= gSup[d] )
                return false;
        }
        return true;
    }
    
    /// periodic image
    void bringInside(int coord[ORD]) const
    {
        for ( int d = 0; d < ORD; ++d )
            coord[d] = imageI(gDim[d], coord[d]);
    }
    
    /// conversion from index to coordinates
    void setCoordinatesFromIndex(int coord[ORD], INDEX indx) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            coord[d] = indx % gDim[d];
            indx    /= gDim[d];
        }
    }
    
    /// conversion from Position to coordinates (offset should be in [0,1])
    void setCoordinatesFromPosition(int coord[ORD], const real w[ORD], const real offset=0) const
    {
        for ( int d = 0; d < ORD; ++d )
            coord[d] = imageF(gDim[d], offset+(w[d]-gInf[d])*cDelta[d]);
    }

    /// conversion from Index to Position (offset should be in [0,1])
    void setPositionFromIndex(real w[ORD], INDEX indx, real offset) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            w[d] = gInf[d] + cWidth[d] * ( offset + indx % gDim[d] );
            indx /= gDim[d];
        }
    }
    
    /// conversion from Coordinates to Position (offset should be in [0,1])
    void setPositionFromCoordinates(real w[ORD], const int coord[ORD], real offset=0) const
    {
        for ( int d = 0; d < ORD; ++d )
            w[d] = gInf[d] + cWidth[d] * ( offset + coord[d] );
    }

    /// conversion from coordinates to index
    INDEX pack(const int coord[ORD]) const
    {
        INDEX inx = imageI(gDim[ORD-1], coord[ORD-1]);
        
        for ( int d = ORD-2; d >= 0; --d )
            inx = gDim[d] * inx + imageI(gDim[d], coord[d]);
        
        return inx;
    }
    
    
    /// returns the index of the cell whose center is closest to the point w[]
    INDEX index(const real w[ORD]) const
    {
        INDEX inx = imageF(gDim[ORD-1], (w[ORD-1]-gInf[ORD-1])*cDelta[ORD-1]);
        
        for ( int d = ORD-2; d >= 0; --d )
            inx = gDim[d] * inx + imageF(gDim[d], (w[d]-gInf[d])*cDelta[d]);
        
        return inx;
    }

    
    /// returns the index of the cell whose center is closest to the point w[]
    INDEX index(const real w[ORD], const real offset) const
    {
        INDEX inx = imageF(gDim[ORD-1], offset+(w[ORD-1]-gInf[ORD-1])*cDelta[ORD-1]);
        
        for ( int d = ORD-2; d >= 0; --d )
            inx = gDim[d] * inx + imageF(gDim[d], offset+(w[d]-gInf[d])*cDelta[d]);
        
        return inx;
    }

    
    /// return cell that is next to `c` in the direction `d`
    INDEX next(INDEX c, int dd) const
    {
        int coord[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            coord[d] = c % gDim[d];
            c       /= gDim[d];
        }

        coord[dd] = imageI(gDim[dd], coord[dd]+1);

        INDEX inx = coord[ORD-1];
        
        for ( int d = ORD-2; d >= 0; --d )
            inx = gDim[d] * inx + coord[d];
        
        return inx;
    }
    
    /// convert coordinate to array index, if ORD==1
    INDEX pack1D(const int x) const
    {
        return imageI(gDim[0], x);
    }
    
    /// convert coordinate to array index, if ORD==2
    INDEX pack2D(const int x, const int y) const
    {
        return imageI(gDim[0], x) + gDim[0]*imageI(gDim[1], y);
    }
    
    /// convert coordinate to array index, if ORD==3
    INDEX pack3D(const int x, const int y, const int z) const
    {
        return imageI(gDim[0], x) + gDim[0]*( imageI(gDim[1], y) + gDim[1]*imageI(gDim[2], z) );
    }

    
    /// return index of cell corresponding to position (x), if ORD==1
    INDEX index1D(const real x) const
    {
        return imageI(gDim[0], (x-gInf[0])*cDelta[0]);
    }
    
    /// return index of cell corresponding to position (x, y), if ORD==2
    INDEX index2D(const real x, const real y) const
    {
        return imageI(gDim[0], (x-gInf[0])*cDelta[0])
               + gDim[0] * imageI(gDim[1], (y-gInf[1])*cDelta[1]);
    }
    
    /// return index of cell corresponding to position (x, y, z), if ORD==3
    INDEX index3D(const real x, const real y, const real z) const
    {
        return imageI(gDim[0], (x-gInf[0])*cDelta[0])
               + gDim[0]*( imageI(gDim[1], (y-gInf[1])*cDelta[1])
               + gDim[1]*  imageI(gDim[2], (z-gInf[2])*cDelta[2]) );
    }

    //--------------------------------------------------------------------------

#pragma mark - Regions

    /** For any cell, we can find the adjacent cells by adding 'index offsets'
    However, the valid offsets depends on wether the cell is on a border or not.
    For each 'edge', a list of offsets and its gDim are stored.*/
    
private:
    
    /// array of index offset to neighbors, for each edge type
    int  * regionsEdge;
    
    /// pointers to regionsEdge[], as a function of cell index
    int ** regions;
    
private:
    
    /// calculate the edge-characteristic from the size `s`, coordinate `c` and range `r`
    static INDEX edge_value(const int s, const int r, const int c)
    {
        if ( c < r )
            return r - c;
        else if ( c + r + 1 > s )
            return 2 * r + c - s + 1;
        else
            return 0;
    }
    
    /// caculate the edge characteristic from the coordinates of a cell and the range vector
    INDEX edgeFromCoordinates(const int coord[ORD], const int range[ORD]) const
    {
        INDEX e = 0;
        for ( int d = ORD-1; d >= 0; --d )
        {
            e *= 2 * range[d] + 1;
            e += edge_value(gDim[d], range[d], coord[d]);
        }
        return e;
    }
    
    
    int * makeRectangularGrid(int& cmx, const int range[ORD])
    {
        cmx = 1;
        for ( int d = 0; d < ORD; ++d )
            cmx *= ( 2 * range[d] + 1 );
        int * ccc = new int[ORD*cmx];
        
        int nb = 1;
        for ( int d = 0; d < ORD; ++d )
        {
            int h = 0;
            for ( ; h < nb && h < cmx; ++h )
                ccc[ORD*h+d] = 0;
            for ( int s = -range[d]; s <= range[d]; ++s )
            {
                if ( s != 0 )
                {
                    for ( int n = 0; n < nb && h < cmx; ++n, ++h )
                    {
                        for ( int e = 0; e < d; ++e )
                            ccc[ORD*h+e] = ccc[ORD*n+e];
                        ccc[ORD*h+d] = s;
                    }
                }
            }
            nb = h;
        }
        assert_true(nb==cmx);
        return ccc;
    }
    
    
    /// calculate cell index offsets between 'ori' and 'ori+shift'
    int calculateOffsets(int offsets[], int shift[], int cnt, int ori[], bool positive)
    {
        int nb = 0;
        int cc[ORD];
        int ori_indx = (int)pack(ori);
        for ( int ii = 0; ii < cnt; ++ii )
        {
            for ( int d = 0; d < ORD; ++d )
                cc[d] = ori[d] + shift[ORD*ii+d];
            int off = (int)pack(cc) - ori_indx;
            
            bool add = ( positive ? off >= 0 : true );
            if ( periodic() )
            {
                //check that cell is not already included:
                for ( int n = 0; n < nb; ++n )
                    if ( offsets[n] == off )
                    {
                        add = false;
                        break;
                    }
            }
            else 
                add &= inside(cc);
            
            if ( add )
                offsets[nb++] = off;
        }
        return nb;
    }
    
   
    /// create regions in the offsets buffer
    /**
     Note: the range is taken specified in units of cells: 1 = 1 cell
     @todo: specify range in calculateRegion() as real distance!
     */
    void createRegions(int * ccc, const int regMax, const int range[ORD], bool positive)
    {
        INDEX edgeMax = 0;
        for ( int d = ORD-1; d >= 0; --d )
        {
            edgeMax *= 2 * range[d] + 1;
            edgeMax += 2 * range[d];
        }
        ++edgeMax;
        
        //allocate and reset arrays:
        deleteRegions();
        
        regions     = new int*[nCells];
        regionsEdge = new int[edgeMax*(regMax+1)];
        for ( INDEX e = 0; e < edgeMax*(regMax+1); ++e )
            regionsEdge[e] = 0;
        
        int ori[ORD];
        for ( INDEX indx = 0; indx < nCells; ++indx )
        {
            setCoordinatesFromIndex(ori, indx);
            INDEX e = edgeFromCoordinates(ori, range);
            assert_true( e < edgeMax );
            int * reg = regionsEdge + e * ( regMax + 1 );
            if ( reg[0] == 0 )
            {
                // calculate the region for this new edge-characteristic
                reg[0] = calculateOffsets(reg+1, ccc, regMax, ori, positive);
                //printf("edge %i has region of %i cells\n", e, reg[0]);
            }
            else if ( 0 )
            {
                // compare result for a different cell of the same edge-characteristic
                int * rig = new int[regMax+1];
                rig[0] = calculateOffsets(rig+1, ccc, regMax, ori, positive);
                if ( rig[0] != reg[0] )
                    ABORT_NOW("inconsistent region size");
                for ( int s = 1; s < rig[0]+1; ++s )
                    if ( rig[s] != reg[s] )
                        ABORT_NOW("inconsistent region offsets");
                delete [] rig;
            }
            regions[indx] = reg;
        }
    }
    
    /// accept within a certain diameter
    bool reject_disc(const int c[ORD], real radius)
    {
        real dsq = 0;
        for ( int d = 0; d < ORD; ++d ) 
            dsq += cWidth[d] * c[d] * cWidth[d] * c[d];
        return ( dsq > radius * radius );
    }
    
    /// accept within a certain diameter
    bool reject_square(const int c[ORD], real radius)
    {
        for ( int d = 0; d < ORD; ++d ) 
            if ( fabs( cWidth[d] * c[d] ) > radius )
                return true;
        return false;
    }
    
    /// used to remove ORD coordinates in array `ccc`
    void remove_entry(int * ccc, int& cmx, const int s)
    {
        --cmx;
        for ( int x = ORD*s; x < ORD*cmx; ++x )
            ccc[x] = ccc[x+ORD];
    }        
    
public:
    
    /// create regions which contains cells at a distance 'range' or less
    /**
     Note: the range is specified in real units.
     the region will cover an area of space that is approximately square.
     */
    void createSquareRegions(const real radius)
    {
        int range[ORD];
        for ( int d = 0; d < ORD; ++d )
            range[d] = ceil( radius / cWidth[d] );
        int cmx = 0;
        int * ccc = makeRectangularGrid(cmx, range);
        
        for ( int s = cmx-1; s >= 0 ; --s )
            if ( reject_square(ccc+ORD*s, radius) )
                remove_entry(ccc, cmx, s);
        
        createRegions(ccc, cmx, range, false);
        delete [] ccc;
    }
    
    /// create regions which contains cells at a distance 'range' or less
    /**
     Note: the range is specified in real units.
     The region covers an area of space that is approximately circular.
     */
    void createRoundRegions(const real radius)
    {
        int range[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            assert_true( cWidth[d] > 0 );
            range[d] = ceil( radius / cWidth[d] );
        }
        int cmx = 0;
        int * ccc = makeRectangularGrid(cmx, range);
       
        for ( int s = cmx-1; s >= 0 ; --s )
            if ( reject_disc(ccc+ORD*s, radius) )
                remove_entry(ccc, cmx, s);
        
        createRegions(ccc, cmx, range, false);
        delete [] ccc;
    }

    /// regions that only contain cells of greater index.
    /**
     This is suitable for pair-wise interaction of particles, since it can
     be used to go through the cells one by one such that at the end, all
     pairs of cells have been considered only once.

     Note: the radius is taken specified in units of cells: 1 = 1 cell
     */
    void createSideRegions(const int radius)
    {
        int range[ORD];
        for ( int d = 0; d < ORD; ++d )
            range[d] = radius;
        int cmx = 0;
        int * ccc = makeRectangularGrid(cmx, range);
        createRegions(ccc, cmx, range, true);
        delete [] ccc;
    }
    
    /// true is createRegions() or createRoundRegions() was called
    bool hasRegions() const
    {
        return ( regions != 0  &&  regionsEdge != 0 );
    }
    
    /// set region array 'offsets' for given cell index
    /**
     A zero offset is always first in the list.
     //\returns the size of the list.
     
     \par Example:
     @code
     CELL * cell = & myGrid.icell(indx);
     n_offset = myGrid.getRegion(offset, indx);
     for ( int n = 1; n < n_offset; ++n )
     {
         Cell & neighbor = cell[offset[n]];
         ...
     }
     @endcode
     
     Note: you must call createRegions() first
    */
    int getRegion(int*& offsets, const INDEX indx) const
    {
        assert_true( hasRegions() );
        offsets = regions[indx]+1;
        assert_true( offsets[0] == 0 );
        return regions[indx][0];
    }
    
    /// free memory occupied by the regions
    void deleteRegions()
    {
        if ( regions )
            delete [] regions;
        regions = 0;
        
        if ( regionsEdge )
            delete [] regionsEdge;
        regionsEdge = 0;
    }

#pragma mark -

    /// write some info on the grid
    void printSummary(std::ostream& os, std::string str)
    {
        os << str << " of dim " << ORD << " has " << allocated << " cells:" << std::endl;
        for ( int d = 0; d < ORD; ++d )
            os << "     [ " << gInf[d] << " " << gSup[d] << " ] / " << gDim[d] << " = " << cWidth[d] << std::endl;
    }
};


#endif

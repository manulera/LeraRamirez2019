// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "assert_macro.h"
#include "point_grid.h"
#include "exceptions.h"
#include "messages.h"
#include "modulo.h"
#include "space.h"
#include "meca.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

PointGrid::PointGrid()
: max_diameter(0)
{
}


void PointGrid::setGrid(Space const* space, Modulo const* modulo, real min_step)
{
    if ( min_step <= REAL_EPSILON )
        return;
    
    Vector inf, sup;
    space->boundaries(inf, sup);
    
    bool periodic = false;
    int n_cell[3];
    for ( int d = 0; d < DIM; ++d )
    {
        real n = ( sup[d] - inf[d] ) / min_step;
        
        if ( n < 0 )
            throw InvalidParameter("invalid space:boundaries");
        
        if ( modulo  &&  modulo->isPeriodic(d) )
        {
            //adjust the grid to match the edges
            n_cell[d] = (int)floor(n);
            if ( n_cell[d] <= 0 )
                n_cell[d] = 1;
            periodic = true;
        }
        else
        {
            //add a border in any dimension which is not periodic
            n_cell[d] = (int)ceil(n) + 2;
            n = n_cell[d] * 0.5 * min_step;
            real mid = inf[d] + sup[d];
            inf[d] = mid - n;
            sup[d] = mid + n;
        }
    }
    
    //create the grid using the calculated dimensions:
    pGrid.periodic(periodic);
    pGrid.setDimensions(inf, sup, n_cell);
    pGrid.createCells();

    //Create side regions suitable for pairwise interactions:
    pGrid.createSideRegions(1);

    //The maximum allowed diameter of particles is half the minimum cell width
    max_diameter = pGrid.minimumWidth(1);

    //report the grid size used
    //pGrid.printSummary(std::clog, "StericGrid");
}


//------------------------------------------------------------------------------
#pragma mark -

/// include verifications that the grid is appropriate for the particule radius
#define CHECK_RANGE 1


#if ( NB_STERIC_PANES != 1 )

void PointGrid::add(unsigned pan, PointExact const& pe, real rd, real rg) const
{
    if ( pan == 0 || pan > NB_STERIC_PANES )
        throw InvalidParameter("object:steric is out-of-range");

    Vector w = pe.pos();
    point_list(w, pan).new_val().set(pe, rd, rg, w);
    
#if ( CHECK_RANGE )
    //we check that the grid would correctly detect collision of two particles
    if ( max_diameter < 2 * rg )
    {
        std::ostringstream oss;
        oss << "simul:steric_max_range is too short" << std::endl;
        oss << PREF << "steric_max_range should be greater than 2 * ( particle-radius + extra-range )" << std::endl;
        oss << PREF << "= " << 2 * rg << " for some particles" << std::endl;
        throw InvalidParameter(oss.str());
    }
#endif
}


void PointGrid::add(unsigned pan, FiberLocus const& fl, real rd, real rg) const
{
    if ( pan == 0 || pan > NB_STERIC_PANES )
        throw InvalidParameter("object:steric is out-of-range");

    // link in the cell containing the middle of the segment:
    Vector w = fl.center();
    locus_list(w, pan).new_val().set(fl, rd, rg);
    
#if ( CHECK_RANGE )
    //we check that the grid would correctly detect collision of two segments
    //along the diagonal, corresponding to the worst-case scenario
    real diag = sqrt( fl.len() * fl.len() + 4 * rg * rg );
    if ( max_diameter < diag )
    {
        std::ostringstream oss;
        oss << "simul:steric_max_range is too short" << std::endl;
        oss << PREF << "steric_max_range should be greater than sqrt( sqr(segment_length) + 4*sqr(range) )" << std::endl;
        oss << PREF << "where normally segment_length ~ 4/3 segmentation" << std::endl;
        oss << PREF << "= " << diag << " for some fibers" << std::endl;
        throw InvalidParameter(oss.str());
    }
#endif
}


#endif


//------------------------------------------------------------------------------
#pragma mark - Steric functions


/**
 This is used to check two spherical objects:
 Solid/Bead/Sphere or Fiber-tip
 
 The force is applied if the objects are closer than the
 sum of their radiuses.
 */
void PointGrid::checkPP(Meca& meca, PointGridParam const& pam,
                        FatPoint const& aa, FatPoint const& bb) const
{
    //std::clog << "   PP- " << bb.pe << " " << aa.pe << std::endl;
    const real len = aa.radius + bb.radius;
    Vector vab = bb.pos - aa.pos;
    
    if ( modulo )
        modulo->fold(vab);
    
    if ( vab.normSqr() < len*len )
        meca.interLongLink(aa.pe, bb.pe, len, pam.stiff_push);
}


/**
 This is used to check a segment of a fiber against a spherical object:
 Solid/Bead/Sphere or Fiber-tip.
 
 The force is applied if the objects are closer than the sum of their radiuses.
 */
void PointGrid::checkPL(Meca& meca, PointGridParam const& pam,
                        FatPoint const& aa, FatLocus const& bb) const
{
    //std::clog << "   PL- " << bb.fl << " " << aa.pe << std::endl;
    const real len = aa.radius + bb.radius;
    
    // get position of point with respect to segment:
    real dis2;
    real abs = bb.fl.projectPoint0(aa.pos, dis2);
    
    if ( 0 <= abs )
    {
        if ( abs <= bb.fl.len() )
        {
            if ( dis2 < len*len )
            {
                PointInterpolated bi(bb.fl, abs);
                meca.interSideSlidingLink(bi, aa.pe, len, pam.stiff_push);
            }
        }
        else
        {
            if ( bb.isLast() )
                checkPP(meca, pam, aa, bb.point2());
        }
    }
    else
    {
        if ( bb.isFirst() )
            checkPP(meca, pam, aa, bb.point1());
        else
        {
            /* we check the projection to the previous segment,
             and if it falls on the right of it, then we interact with the node */
            Vector vab = aa.pos - bb.fl.pos1();
            
            if ( modulo )
                modulo->fold(vab);
            
            if ( vab * bb.fl.fiber()->diffPoints(bb.fl.point()-1) >= 0 )
            {
                if ( vab.normSqr() < len*len )
                    meca.interLongLink(aa.pe, bb.fl.exact1(), len, pam.stiff_push);
            }
        }
    }
}


/**
 This is used to check a segment of a fiber against the non-end model-point of a fiber.
 
 The interaction is applied only if the model-point projects 'inside' the segment.
 */
void PointGrid::checkLL1(Meca& meca, PointGridParam const& pam,
                         FatLocus const& aa, FatLocus const& bb) const
{
    //std::clog << "   LL1 " << aa.fl << " " << bb.point1() << std::endl;
    const real ran = aa.range + bb.radius;
    
    // get position of bb.point1() with respect to segment 'aa'
    real dis2 = INFINITY;
    real abs = aa.fl.projectPoint0(bb.fl.pos1(), dis2);
    
    if ( dis2 < ran*ran )
    {
        /*
         bb.point1() projects inside segment 'aa'
         */
        assert_true( 0 <= abs  &&  abs <= aa.fl.len() );
        const real len = aa.radius + bb.radius;
        PointInterpolated ai(aa.fl, abs);
        if ( dis2 > len*len )
            meca.interSideSlidingLink(ai, bb.fl.exact1(), len, pam.stiff_pull);
        else
            meca.interSideSlidingLink(ai, bb.fl.exact1(), len, pam.stiff_push);
    }
    else if ( abs < 0 )
    {
        if ( aa.isFirst() )
        {
            /*
             Check the projection of aa.point1(),
             on the segment represented by 'bb'
             */
            if ( &bb < &aa  &&  bb.isFirst() )
            {
                Vector vab = bb.fl.pos1() - aa.fl.pos1();
                
                if ( modulo )
                    modulo->fold(vab);
                
                const real len = aa.radius + bb.radius;
                if ( vab.normSqr() < len*len  &&  vab * bb.fl.diff() >= 0 )
                    meca.interLongLink(aa.fl.exact1(), bb.fl.exact1(), len, pam.stiff_push);
            }
        }
        else
        {
            /*
             Check the projection to the segment located before 'aa',
             and interact if 'bb.point1()' falls on the right side of it
             */
            Vector vab = bb.fl.pos1() - aa.fl.pos1();
            
            if ( modulo )
                modulo->fold(vab);
            
            if ( vab * aa.fl.fiber()->diffPoints(aa.fl.point()-1) >= 0 )
            {
                const real d = vab.normSqr();
                if ( d < ran*ran )
                {
                    const real len = aa.radius + bb.radius;
                    if ( d > len*len )
                        meca.interLongLink(aa.fl.exact1(), bb.fl.exact1(), len, pam.stiff_push);
                    else
                        meca.interLongLink(aa.fl.exact1(), bb.fl.exact1(), len, pam.stiff_push);
                }
            }
        }
    }
}


/**
 This is used to check a segment of a fiber against the non-end model-point of a fiber.
 
 The interaction is applied only if the model-point projects 'inside' the segment.
 */
void PointGrid::checkLL2(Meca& meca, PointGridParam const& pam,
                         FatLocus const& aa, FatLocus const& bb) const
{
    //std::clog << "   LL2 " << aa.fl << " " << bb.point2() << std::endl;
    const real ran = aa.range + bb.radius;
    
    // get position of bb.point2() with respect to segment 'aa'
    real dis2 = INFINITY;
    real abs = aa.fl.projectPoint0(bb.fl.pos2(), dis2);
    
    if ( dis2 < ran*ran )
    {
        /*
         bb.point2() projects inside segment 'aa'
         */
        assert_true( 0 <= abs  &&  abs <= aa.fl.len() );
        const real len = aa.radius + bb.radius;
        PointInterpolated ai(aa.fl, abs);
        if ( dis2 > len*len )
            meca.interSideSlidingLink(ai, bb.fl.exact2(), len, pam.stiff_pull);
        else
            meca.interSideSlidingLink(ai, bb.fl.exact2(), len, pam.stiff_push);
    }
    else if ( abs < 0 )
    {
        /*
         Check the projection to the segment located before 'aa',
         and interact if 'bb.point1()' falls on the right side of it
         */
        Vector vab = bb.fl.pos2() - aa.fl.pos1();
        
        if ( modulo )
            modulo->fold(vab);
        
        if ( aa.isFirst() )
        {
            assert_true(bb.isLast());
            const real len = aa.radius + bb.radius;
            if ( vab.normSqr() < len*len  &&  vab * bb.fl.diff() <= 0 )
                meca.interLongLink(aa.fl.exact1(), bb.fl.exact2(), len, pam.stiff_push);
        }
        else
        {
            if ( vab * aa.fl.fiber()->diffPoints(aa.fl.point()-1) >= 0 )
            {
                const real d = vab.normSqr();
                if ( d < ran*ran )
                {
                    const real len = aa.radius + bb.radius;
                    if ( d > len*len )
                        meca.interLongLink(aa.fl.exact1(), bb.fl.exact2(), len, pam.stiff_push);
                    else
                        meca.interLongLink(aa.fl.exact1(), bb.fl.exact2(), len, pam.stiff_push);
                }
            }
        }
    }
    else if ( &bb < &aa  &&  aa.isLast()  &&  abs > aa.fl.len() )
    {
        /*
         Check the projection of aa.point2(),
         on the segment represented by 'bb'
         */
        assert_true(bb.isLast());
        
        Vector vab = bb.fl.pos2() - aa.fl.pos2();
        
        if ( modulo )
            modulo->fold(vab);
        
        const real len = aa.radius + bb.radius;
        if ( vab.normSqr() < len*len  &&  vab * bb.fl.diff() <= 0 )
            meca.interLongLink(aa.fl.exact2(), bb.fl.exact2(), len, pam.stiff_push);
    }
}


/**
 This is used to check two FiberLocus, that each represent a segment of a Fiber.
 The segments are tested for intersection in 3D.
 */
void PointGrid::checkLL(Meca& meca, PointGridParam const& pam,
                        FatLocus const& aa, FatLocus const& bb) const
{
    //std::clog << "LL " << aa.fl << " " << bb.fl << std::endl;
    checkLL1(meca, pam, aa, bb);
    
    if ( aa.isLast() )
        checkLL2(meca, pam, bb, aa);
        
    checkLL1(meca, pam, bb, aa);
      
    if ( bb.isLast() )
        checkLL2(meca, pam, aa, bb);
  
#if ( DIM == 3 )
    
    const real ran = std::max(aa.range+bb.radius, aa.radius+bb.range);

    /* in 3D, we use shortestDistance() to calculate the closest distance
     between two segments, and use the result to build an interaction */
    real a, b, d;
    if ( aa.fl.shortestDistance(bb.fl, a, b, d)  &&  d < ran*ran )
    {
        const real len = aa.radius + bb.radius;
        
        PointInterpolated ai(aa.fl, a);
        PointInterpolated bi(bb.fl, b);

        //std::clog << "steric distance " << d << "  " << ai << " " << bi <<"\n";
     
        if ( d > len*len )
            meca.interSideSlidingLink(ai, bi, len, pam.stiff_pull);
        else
            meca.interSideSlidingLink(ai, bi, len, pam.stiff_push);
    }
    
#endif
}


//------------------------------------------------------------------------------
#pragma mark -


/**
 This will consider once all pairs of objects from the given lists
 */
void PointGrid::setInteractions(Meca& meca, PointGridParam const& pam,
                                FatPointList & fpl, FatLocusList & fll) const
{
    for ( FatPoint* ii = fpl.begin(); ii < fpl.end(); ++ii )
    {
        for ( FatPoint* jj = ii+1; jj < fpl.end(); ++jj )
            checkPP(meca, pam, *ii, *jj);
        
        for ( FatLocus* kk = fll.begin(); kk < fll.end(); ++kk )
            checkPL(meca, pam, *ii, *kk);
    }
    
    for ( FatLocus* ii = fll.begin(); ii < fll.end(); ++ii )
    {
        for ( FatLocus* jj = ii+1; jj < fll.end(); ++jj )
            checkLL(meca, pam, *ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated
 */
void PointGrid::setInteractions(Meca& meca, PointGridParam const& pam,
                                FatPointList & fpl1, FatLocusList & fll1,
                                FatPointList & fpl2, FatLocusList & fll2) const
{
    assert_true( &fpl1 != &fpl2 );
    assert_true( &fll1 != &fll2 );
    
    for ( FatPoint* ii = fpl1.begin(); ii < fpl1.end(); ++ii )
    {
        for ( FatPoint* jj = fpl2.begin(); jj < fpl2.end(); ++jj )
            checkPP(meca, pam, *ii, *jj);
        
        for ( FatLocus* kk = fll2.begin(); kk < fll2.end(); ++kk )
            checkPL(meca, pam, *ii, *kk);
    }
    
    for ( FatLocus* ii = fll1.begin(); ii < fll1.end(); ++ii )
    {
        for ( FatPoint* jj = fpl2.begin(); jj < fpl2.end(); ++jj )
            checkPL(meca, pam, *jj, *ii);
        
        for ( FatLocus* kk = fll2.begin(); kk < fll2.end(); ++kk )
            checkLL(meca, pam, *ii, *kk);
    }
}



#if ( NB_STERIC_PANES == 1 )

/**
 Check interactions between objects contained in the grid.
 */
void  PointGrid::setInteractions(Meca& meca, PointGridParam const& pam) const
{
    assert_true(pam.stiff_push >= 0);
    assert_true(pam.stiff_pull >= 0);
    //std::clog << "----" << std::endl;

    // scan all cells to examine each pair of particles:
    for ( unsigned inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        // We consider each pair of objects (ii, jj) only once:
        
        FatPointList & baseP = point_list(inx);
        FatLocusList & baseL = locus_list(inx);
        
        setInteractions(meca, pam, baseP, baseL);
        
        for ( int reg = 1; reg < nr; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg]);
            FatLocusList & sideL = locus_list(inx+region[reg]);
            
            setInteractions(meca, pam, baseP, baseL, sideP, sideL);
        }
    }
}


#else

/**
 Check interactions between the FatPoints contained in Pane `pan`.
 */
void  PointGrid::setInteractions(Meca& meca, PointGridParam const& pam,
                                 const unsigned pan) const
{
    assert_true(pam.stiff_push >= 0);
    assert_true(pam.stiff_pull >= 0);
    
    // scan all cells to examine each pair of particles:
    for ( unsigned inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        // We consider each pair of objects (ii, jj) only once:
        
        FatPointList & baseP = point_list(inx, pan);
        FatLocusList & baseL = locus_list(inx, pan);
        
        setInteractions(meca, pam, baseP, baseL);

        for ( int reg = 1; reg < nr; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg], pan);
            FatLocusList & sideL = locus_list(inx+region[reg], pan);
            
            setInteractions(meca, pam, baseP, baseL, sideP, sideL);
        }
    }
}



/**
 Check interactions between the FatPoints contained in Panes `pan1` and `pan2`,
 where ( pan1 != pan2 )
 */
void  PointGrid::setInteractions(Meca& meca, PointGridParam const& pam,
                                 const unsigned pan1, const unsigned pan2) const
{
    assert_true(pam.stiff_push >= 0);
    assert_true(pam.stiff_pull >= 0);
    assert_true(pan1 != pan2);
    
    // scan all cells to examine each pair of particles:
    for ( unsigned inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);

        // We consider each pair of objects (ii, jj) only once:
        
        FatPointList & baseP = point_list(inx, pan1);
        FatLocusList & baseL = locus_list(inx, pan1);

        for ( int reg = 0; reg < nr; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg], pan2);
            FatLocusList & sideL = locus_list(inx+region[reg], pan2);

            setInteractions(meca, pam, baseP, baseL, sideP, sideL);
        }
        
        FatPointList & baseP2 = point_list(inx, pan2);
        FatLocusList & baseL2 = locus_list(inx, pan2);
        
        for ( int reg = 1; reg < nr; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg], pan1);
            FatLocusList & sideL = locus_list(inx+region[reg], pan1);
            
            setInteractions(meca, pam, baseP2, baseL2, sideP, sideL);
        }
    }
}


#endif


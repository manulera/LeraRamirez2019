// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "assert_macro.h"
#include "rasterizer.h"
#include "fiber_grid.h"
#include "exceptions.h"
#include "fiber_locus.h"
#include "fiber_binder.h"
#include "messages.h"
#include "space.h"
#include "modulo.h"
#include "hand.h"
#include "hand_prop.h"
#include "simul.h"
#include "sim.h"

extern Random RNG;
extern Modulo const* modulo;

#if ( 0 )
// this includes a naive implementation, which is slow but helpful for debugging
#   include "fiber_grid2.cc"
#else

/**
 Creates a grid where the dimensions of the cells are `max_step` at most.
 If the numbers of cells that need to be created is greater than `max_nb_cells`,
 the function returns 1 without building the grid.
 The return value is zero in case of success.
 
 The algorithm works with any value of `max_step` (the results are always correct),
 but `max_step` affects the efficiency (speed) of the algorithm:
 -if `max_step` is too small, paintGrid() will be slow,
 -if `max_step` is too large, tryToAttach() will be slow.
 A good compromise is to set `max_step` equivalent to the attachment distance,
 or at least to the size of the segments of the Fibers.
 */
unsigned FiberGrid::setGrid(const Space * space, real max_step)
{
    if ( max_step <= 0 )
        throw InvalidParameter("simul:binding_grid_step should be > 0");
    
    //reset gridRange to trigger an error if paintGrid() is not called:
    gridRange = -1;
    
    Vector inf, sup;
    space->boundaries(inf, sup);
    
    bool periodic = false;
    int n_cell[3] = { 1, 1, 1 };
    
    for ( int d = 0; d < DIM; ++d )
    {
        n_cell[d] = (int) ceil( ( sup[d] - inf[d] ) / max_step );
        
        if ( n_cell[d] < 0 )
            throw InvalidParameter("invalid space:boundaries");
        
        if ( modulo  &&  modulo->isPeriodic(d) )
        {
            //adjust the grid to match the edges exactly
            periodic = true;
        }
        else
        {
            //extend the grid by one cell on each side
            inf[d]    -= max_step;
            sup[d]    += max_step;
            n_cell[d] += 2;
        }
        
        if ( n_cell[d] <= 0 )
            n_cell[d] = 1;
    }

    //create the grid using the calculated dimensions:
    if ( periodic )
        fGrid.periodic(true);
    
    fGrid.setDimensions(inf, sup, n_cell);
    return fGrid.nbCells();
}


void FiberGrid::createGrid()
{
    fGrid.createCells();
    
    //fGrid.printSummary(std::cerr, "FiberGrid");
}


unsigned FiberGrid::hasGrid() const
{
    return fGrid.hasCells();
}


void FiberGrid::setGridRange(PropertyList const& properties)
{
    PropertyList plist = properties.find_all("hand");
    
    gridRange = 0;
    for ( PropertyList::const_iterator n = plist.begin(); n != plist.end(); ++n )
    {
        HandProp const* hap = static_cast<HandProp const*>(*n);
        gridRange = std::max(gridRange, hap->binding_range);
    }
    
    //std::clog << "FiberGrid:gridRange = " << gridRange << "\n";
}

//------------------------------------------------------------------------------
#pragma mark - Paint

/** 
 paintCell(x,y,z) adds a Segment to the SegmentList associated with
 the grid point (x,y,z). 
 It is called by the rasterizer function paintFatLine().
 
 This version uses the fact that cells with consecutive
 X-coordinates should be consecutive also in the Grid
 */

struct PaintJob
{
    FiberLocus      const* segment;
    FiberGrid::grid_type * grid;
};


void paintCell(const int x_inf, const int x_sup, const int y, const int z, void * arg)
{
    PaintJob const* job = static_cast<PaintJob const*>(arg);
    FiberLocus const* seg = job->segment;
    FiberGrid::grid_type * grid = job->grid;
    //printf("paint %p in (%i to %i, %i, %i)\n", seg, x_inf, x_sup, y, z);

#if   ( DIM == 1 )
    FiberGrid::SegmentList * inf = & grid->icell1D( x_inf );
    FiberGrid::SegmentList * sup = & grid->icell1D( x_sup );
#elif ( DIM == 2 )
    FiberGrid::SegmentList * inf = & grid->icell2D( x_inf, y );
    FiberGrid::SegmentList * sup = & grid->icell2D( x_sup, y );
#elif ( DIM == 3 )
    FiberGrid::SegmentList * inf = & grid->icell3D( x_inf, y, z );
    FiberGrid::SegmentList * sup = & grid->icell3D( x_sup, y, z );
#endif
    
    for ( FiberGrid::SegmentList * list = inf; list <= sup; ++list )
        list->push_back(seg);
}


/** 
 paintCellPeriodic(x,y,z) adds a Segment in the SegmentList associated with
 the grid point (x,y,z). 
 It is called by the rasterizer function paintFatLine()
 */

void paintCellPeriodic(const int x_inf, const int x_sup, const int y, const int z, void * arg)
{
    PaintJob const* job = static_cast<PaintJob const*>(arg);
    FiberLocus const* seg = job->segment;
    FiberGrid::grid_type * grid = job->grid;
    //printf("paint %p in (%i to %i, %i, %i)\n", seg, x_inf, x_sup, y, z);
    
    for ( int x = x_inf; x <= x_sup; ++x )
    {
        //@todo write/call a specialized function for periodic: icellP1D
#if   ( DIM == 1 )
        grid->icell1D( x ).push_back(seg);
#elif ( DIM == 2 )
        grid->icell2D( x, y ).push_back(seg);
#elif ( DIM == 3 )
        grid->icell3D( x, y, z ).push_back(seg);
#endif
    }
}


/**
paintGrid(first_fiber, last_fiber) links all segments found in 'fiber' and its
 descendant, in the point-list GP that match distance(GP, segment) < H.
 
 'H' is calculated such that tryToAttach() finds any segment closer than 'gridRange':
 
 To determine H, we start from a relation on the sides of a triangle:
 (A) distance( GP, segment ) < distance( GP, X ) + distance( X, segment )
 where GP (grid-point) is the closest point on the grid to X.
 
 Since GP in tryToAttach() is the closest point on fGrid to X, we have:
 (B) distance( GP, X ) < 0.5 * fGrid.diagonalLength()
 
 Thus to find all rods for which:
 (B) distance( X, segment ) < gridRange
 we simply use:
 H =  max_range + 0.5 * fGrid.diagonalLength();
 
 Note: H is calculated by paintGrid() and gridRange by setGridRange().
 
 Linking all segments is done in an inverse way:
 for each segment, we cover all points of the grid inside a volume obtained
 by inflating the segment by the length H. We use for that the raterizer which
 calls the function paint() above.
 */

void FiberGrid::paintGrid(const Fiber * first, const Fiber * last)
{
    assert_true(hasGrid());
    assert_true(gridRange >= 0);
    
    fGrid.clear();
    const real* offset = fGrid.inf();
    const real* deltas = fGrid.delta();
    real width = gridRange + 0.5 * fGrid.diagonalLength();
    
    //define the painting function used:
    void (*paint)(int, int, int, int, void*) = modulo ? paintCellPeriodic : paintCell;
    
    PaintJob job;
    job.grid = &fGrid;

    for ( const Fiber * fib = first; fib != last ; fib=fib->next() )
    {
        Vector Q, P = fib->posP(0);
        real S = fib->segmentation();
        
        for ( unsigned pp = 1; pp < fib->nbPoints(); ++pp )
        {
            job.segment = &(fib->locus(pp-1));

            if ( pp & 1 )
                Q = fib->posP(pp);
            else
                P = fib->posP(pp);
            
#if   (DIM == 1)
            Rasterizer::paintFatLine1D(paint, &job, P, Q, width, offset, deltas);
#elif (DIM == 2)
            Rasterizer::paintFatLine2D(paint, &job, P, Q, S, width, offset, deltas);
#elif (DIM == 3)
            //Rasterizer::paintHexLine3D(paint, &job, P, Q, S, width, offset, deltas);
            Rasterizer::paintFatLine3D(paint, &job, P, Q, S, width, offset, deltas);
            //Rasterizer::paintBox3D(paint, &job, P, Q, width, offset, deltas);
#endif
        }
    }
}



//------------------------------------------------------------------------------
#pragma mark - Access

/**
 This will bind the given Hand to any Fiber found within `binding_range`, with a
 probability that is encoded in `prob`.
 The test is `RNG.pint() < prob`, and with 'prob = 1<<30', the chance is 1/4.
 The result is thus stochastic, and will depend on the number of Fiber
 within the range, but it will saturate if there are more than '4' possible targets.
 
 NOTE:
 The distance at which Fibers are detected is limited to the range given in paintGrid()
 */
void FiberGrid::tryToAttach(Vector const& place, Hand& ha) const
{
    assert_true( hasGrid() );
    
    //get the grid node list index closest to the position in space:
    const unsigned indx = fGrid.index(place, 0.5);
    
    //get the list of rods associated with this cell:
    SegmentList & segments = fGrid.icell(indx);
   
    if ( segments.empty() )
        return;

    //randomize the list, to make attachments more fair:
    //this might not be necessary, since the MT list is already mixed
    segments.mix(RNG);
    
    //std::clog << "tryToAttach has " << segments.size() << " segments\n";
    real overlap_affinity = 1;
    if (ha.prop->overlap_affinity && ha.otherHand() and ha.otherHand()->attached()) {
        overlap_affinity = ha.prop->overlap_affinity;
    }
    
    for ( SegmentList::iterator si = segments.begin(); si < segments.end(); ++si )
    {
        FiberLocus const* loc = *si;
        
        real dis = INFINITY;
        // Compute the distance from the hand to the rod, and abscissa of projection:
        real abs = loc->projectPoint(place, dis);      // always works
        //real abs = loc->projectPointF(place, dis);    // faster, but not compatible with periodic boundaries
        
        /* 
         Compare to the maximum attachment range of the hand,
         and compare a newly tossed random number with 'prob'
         */
#ifndef TRICKY_HAND_ATTACHMENT
        
        if ( dis < ha.prop->binding_range_sqr && RNG.test(ha.prop->binding_rate_dt*overlap_affinity) )
#else
        if ( dis < ha.prop->binding_range_sqr && RNG.flip_4th() )
#endif
        {
            Fiber * fib = const_cast<Fiber*>(loc->fiber());
            FiberBinder pos(fib, loc->abscissa1()+abs);

#ifdef MULTI_LATTICE
#ifdef TRAP_SINGLES
            if (!ha.otherHand())
            {
                if (ha.trappedHaMon())
                {
                    ha.set_lat_id(ha.trappedHaMon()->partnerLattice(&ha, fib));
                }
                else
                    ha.random_multi_lattice();
            }
            else if(!ha.otherHand()->attached())
                ha.random_multi_lattice();
#else
            // If it is a single or a FF couple, asign a random value of lattice when binding to the fiber
            if (!ha.otherHand()||!ha.otherHand()->attached())
            {
                ha.random_multi_lattice();
            }
#endif
#endif
            int allowed =ha.attachmentAllowed(pos);
            if ( allowed==1)
            {
                ha.attach(pos);
                return;
            }
            // Here its the other hand that should move to the right lattice, at least for the case of multi_lattice
            else if (allowed==2 && ha.otherHand()->attachmentSecondTry(pos))
            {
                ha.attach(pos);
                return;
            }
        }
#ifndef TRICKY_HAND_ATTACHMENT
        if (ha.prop->overlap_affinity && dis<ha.prop->binding_range_sqr) {
            overlap_affinity=ha.prop->overlap_affinity;
        }
#endif
        
    }
}



/**
 This function is limited to the range given in paintGrid();
 */
FiberGrid::SegmentList FiberGrid::nearbySegments(Vector const& place, const real D, Fiber * exclude)
{
    if ( gridRange < 0 )
        throw InvalidParameter("the Grid was not initialized");
    if ( gridRange < D )
    {
        printf("gridRange = %.4f < range = %.4f\n", gridRange, D);
        throw InvalidParameter("the Grid maximum distance was exceeded");
    }
    
    SegmentList res;
    
    //get the grid node list index closest to the position in space:
    const unsigned indx = fGrid.index(place, 0.5);
    
    //get the list of rods associated with this cell:
    SegmentList & segments = fGrid.icell(indx);
    
    const real DD = D*D;
    for ( SegmentList::iterator si = segments.begin(); si < segments.end(); ++si )
    {
        FiberLocus const* loc = *si;
        
        if ( loc->fiber() == exclude ) 
            continue;
        
        real dis = INFINITY;
        loc->projectPoint(place, dis);
        
        if ( dis < DD )
            res.push_back(loc);
    }
    
    return res;
}



FiberLocus const* FiberGrid::closestSegment(Vector const& place)
{
    //get the cell index from the position in space:
    const unsigned indx = fGrid.index(place, 0.5);
    
    //get the list of rods associated with this cell:
    SegmentList & segments =  fGrid.icell(indx);
    
    FiberLocus const* res = 0;
    real closest = 4 * gridRange * gridRange;
    
    for ( SegmentList::iterator si = segments.begin(); si < segments.end(); ++si )
    {
        FiberLocus const* loc = *si;
        
        //we compute the distance from the hand to the candidate rod,
        //and compare it to the best we have so far.
        real dis = INFINITY;
        loc->projectPoint(place, dis);
        
        if ( dis < closest )
        {
            closest = dis;
            res = loc;
        }
    }
    return res;
}


#endif


//============================================================================
//===                        TEST  ATTACHMENT                             ====
//============================================================================
#pragma mark - Test

#include <map>
#include "simul.h"

/**
Function testAttach() is given a position in space,
 it calls tryToAttach() from this position to check that:
 - attachement has equal probability to all targets,
 - no target is missed,
 - attachment are not made to targets that are beyond binding_range
 */
void FiberGrid::testAttach(FILE * out, const Vector pos, Fiber * start, HandProp const* hp)
{
    //create a test motor with a dummy HandMonitor:
    HandMonitor hm;
    Hand ha(hp, &hm);
    real dsq = hp->binding_range_sqr;
    
    typedef std::map < FiberLocus const*, int > map_type;
    map_type hits;
    
    //go through all the segments to find those close enough from pos:
    for ( Fiber * fib=start; fib; fib=fib->next() )
    {
        for ( unsigned p = 0; p < fib->nbSegments(); ++p )
        {
            FiberLocus const& loc = fib->locus(p);
            real dis = INFINITY;
            loc.projectPoint(pos, dis);
            
            if ( dis < dsq )
                hits[&loc] = 0;
        }
    }
    
    const unsigned targets = hits.size();
    
    if ( targets == 0 )
    {
        //fprintf(out, "no target here\n");
        return;
    }
    
    //call tryTyAttach NB times to check to which rods the Hand binds:
    const unsigned NB = 100 * targets;
    for ( unsigned n = 0; n < NB; ++n )
    {
        // we set 'prob' to bind immediately
        tryToAttach(pos, ha);
        if ( ha.attached() )
        {
            PointInterpolated inter = ha.fiber()->interpolate(ha.abscissa());
            FiberLocus const& loc = ha.fiber()->locus(inter.point1());
            
            if ( hits.find(&loc) != hits.end() )
                ++hits[&loc];
            else
                hits[&loc] = -2;
            
            ha.detach();
        }
    }
    
    //detect segments that have been missed or mistargeted:
    int verbose = 1;
    for ( map_type::const_iterator it = hits.begin(); it != hits.end(); ++it )
    {
        if ( it->second <= 50 )
            verbose = 1;
        if ( it->second < 0 )
            verbose = 2;
    }
    
    if ( verbose )
    {
        // print a summary of all targets:
        fprintf(out, "FiberGrid::testAttach %i target(s) within %.3f um of", targets, hp->binding_range);
        pos.print(out);

#if ( 0 )
        const unsigned indx = fGrid.index(pos, 0.5);
        SegmentList & segments = fGrid.icell(indx);
        for ( SegmentList::iterator si = segments.begin(); si < segments.end(); ++si )
        {
            FiberLocus const* loc = *si;
            fprintf(out, "\n    target f%04ld:%02i", loc->fiber()->identity(), loc->point());
        }
#endif
        
        //report for all the rods that were targeted:
        for ( map_type::const_iterator it = hits.begin(); it != hits.end(); ++it )
        {
            FiberLocus const* loc = it->first;
            Fiber const* fib = loc->fiber();
            real dis = INFINITY;
            real abs = loc->projectPoint(pos, dis);
            
            fprintf(out, "\n    rod f%04ld:%02i at %5.3f um, abs %+.2f : ", fib->identity(), loc->point(), dis, abs);
            if ( hits[loc] == 0 )
                fprintf(out, "missed");
            else if ( hits[loc] < 0 )
                fprintf(out, "found, although out of range");
            else if ( hits[loc] > 0 )
                fprintf(out, "%-3i hits", hits[loc]);
        }
        fprintf(out, "\n");
    }
}


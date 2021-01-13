// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/** 
 This file implements a 'dummy' grid using STL code, which can be used as a reference.
 For each position, it calculates the geometrical distance to all fiber segments.
 This is algorithmically the slowest method, but it is simple and most likely correct!
 It is useful to get a ground truth and evaluate more advanced methods.
 */


#include <algorithm>

typedef std::vector <FiberLocus const*> SegmentVector;

/// a list containing all segments, as a global variable
SegmentVector allSegments;


unsigned FiberGrid::setGrid(const Space *, real)
{
    PRINT_ONCE("Cytosim is not using a grid to find attachement to Fiber!\n");
    gridRange = 0;
    return 0;
}

void FiberGrid::paintGrid(const Fiber * first, const Fiber * last)
{
    allSegments.clear();
    //we go through all the segments
    for ( const Fiber * fb = first ; fb != last ; fb=fb->next() )
    {
        for ( unsigned int sg = 0; sg < fb->nbSegments(); ++sg )
            allSegments.push_back( &(fb->locus(sg)) );
    }
}

void FiberGrid::createGrid()
{
}

unsigned FiberGrid::hasGrid() const
{
    return 1;
}

void FiberGrid::setGridRange(PropertyList const&)
{
}

void FiberGrid::tryToAttach(Vector const& place, Hand& ha) const
{
    // randomize the list order
    std::random_shuffle( allSegments.begin(), allSegments.end() );

    // test all segments:
    for ( SegmentVector::iterator seg = allSegments.begin(); seg < allSegments.end(); ++seg )
    {
        FiberLocus const* loc = *seg;
                
        real dis = INFINITY;
        // Compute the distance from the hand to the rod, and abscissa of projection:
        real abs = loc->projectPoint(place, dis);
        
        /*
         Compare to the maximum attachment range of the hand,
         and compare a newly tossed random number with 'prob'
         */
#ifndef TRICKY_HAND_ATTACHMENT
        if ( dis < ha.prop->binding_range_sqr && RNG.preal() < ha.prop->binding_rate_dt )
#else
        if ( dis < ha.prop->binding_range_sqr && RNG.flip_4th() )
#endif
        {
            Fiber * fib = const_cast<Fiber*>(loc->fiber());
            FiberBinder bind(fib, loc->abscissa1()+abs);
            
            if ( ha.attachmentAllowed(bind) )
            {
                ha.attach(bind);
                return;
            }
        }
    }
}


/**
 This function is limited to the range given in paintGrid();
 */
FiberGrid::SegmentList FiberGrid::nearbySegments( Vector const& place, const real D, Fiber * exclude )
{
    SegmentList res;
    
    const real DD = D * D;
    for ( SegmentVector::iterator seg = allSegments.begin(); seg < allSegments.end(); ++seg )
    {
        FiberLocus const* loc = *seg;
             
        if ( loc->fiber() == exclude )
            continue;
        
        real dis = INFINITY;
        // Compute the distance from the hand to the rod:
        loc->projectPoint(place, dis);
        
        if ( dis < DD )
            res.push_back(loc);
    }
    
    return res;
}

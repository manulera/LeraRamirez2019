// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 return the maximum segmentation of all existing FiberProp,
 multiplied by 0.5
 */
#include "couple_set.h"
#include "single_set.h"

real Simul::estimateFiberGridStep() const
{
    real res = 0;
    
    PropertyList plist = properties.find_all("fiber");
    
    for ( PropertyList::iterator pi = plist.begin(); pi < plist.end(); ++pi )
    {
        FiberProp * fp = static_cast<FiberProp*>(*pi);
        if ( res < fp->segmentation )
            res = fp->segmentation;
    }
    
    return res * 0.5;
}



/**
 The FiberGrid is used to quickly find the fibers that are close to any point.
 Procedure:
 1. if binding_grid_step is not set, attempt to find a suitable value for it,
 2. if the number of cells is superior to 1e5, double the step size,
 2. initialize the grid with this calculated step size.
 */
void Simul::setFiberGrid(Space const* spc) const
{
    assert_true(spc);
    real step = prop->binding_grid_step;
    
    // try to find cell size from the filaments
    if ( step <= 0 )
        step = estimateFiberGridStep();

    /// otherwise, try to get it from the space
    if ( step <= 0 )
        step = spc->max_extension() / 128;
    
    assert_true( step > 0 );

    // increase the cell size until we get reasonable memory requirements:
    const unsigned max = 1<<17;
    while ( 1 )
    {
        if ( fiberGrid.setGrid(spc, step) < max )
            break;
        step *= 2;
    }
    
    // create the grid:
    fiberGrid.createGrid();
    prop->binding_grid_step = step;
    //std::clog << "simul:binding_grid_step = " << step << std::endl;
    
    fiberGrid.setGridRange(properties);
}



/**
 Will pepare the simulation engine to make it ready to make a step():
 - set FiberGrid used for attachment of Hands,
 - set StericGrid
 - call complete() for all registered Property
 .
 The simulated objects should not be changed.
 
 */
void Simul::prepare()
{
    if ( sSpace == 0 )
        throw InvalidSyntax("A space must be defined first!");

    // make sure properties are ready for simulations:
    prop->strict = 1;
    prop->complete(this);
    
    // prepare grid for attachments:
    setFiberGrid(sSpace);
    
    // this is necessary for diffusion in Field:
    fields.prepare();
    
    // this prepares for 'fast_diffusion=1':
    singles.prepare(properties);
    couples.prepare(properties);
    
    prop->strict = 0;
}

/**
 This is the master Monte-Carlo step function.
 
 Lists are mixed such that objects are considered in a different
 and random order at each step, to avoid biais in the simulation

 step() is called for every list, i.e. for every Object
 */
void Simul::step()
{
    // increment time:
    prop->time += prop->time_step;

    // mix every list
    events.mix();
    organizers.mix();
    beads.mix();
    solids.mix();
    fibers.mix();
    spheres.mix();
    couples.mix();
    singles.mix();
    spaces.mix();
    
    // Monte-Carlo step for every object
    events.step();
    organizers.step();
    fields.step();
    spaces.step();
    spheres.step();
    beads.step();
    solids.step();
    fibers.step();
   
    // distribute fibers over the attachment grid:
    fiberGrid.paintGrid(fibers.first(), 0);
    
    
#if ( 0 )
    
    // This code continuously tests the binding algorithm.
    
    if ( HandProp::binding_range_max > 0 ) 
    {
        HandProp hp("test_binding");
        hp.binding_rate  = 1;
        hp.binding_range = HandProp::binding_range_max;
        hp.bind_also_ends = true;
        hp.complete(this);
        
        for ( unsigned cnt = 0; cnt < 16; ++cnt )
        {
            Vector pos = sSpace->randomPlace();
            fiberGrid.testAttach(stdout, pos, fibers.first(), &hp);
        }
    }
    
#endif
    
    // step Hand-containing objects that may attach to fibers:
    for (int i=0; i<prop->handmonitor_pace; i++)
    {
#ifdef TRAP_SINGLES
#if (TRAP_SINGLES==2)
        trap_solution();
#endif
#endif
        singles.step(fibers, fiberGrid);
        couples.step(fibers, fiberGrid);
        singles.checkTrapped();
    }
    
    
    // ADDED BY MANU
    spindle_watch();
}


#ifdef TRAP_SINGLES
#if (TRAP_SINGLES==2)
void Simul::trap_solution()
{
    couples.populateTrapList();
    singles.populateTrapList();
    couples.trap(&singles.simulTrapReserve);
    couples.untrap();
    couples.clearTrapReserve();
    singles.clearTrapReserve();
}
#endif
#endif

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 Calculate the grid size automatically for dynamic Fiber,
 Bead, Sphere and Solid.
 
 This function can be used to set SimulProp::steric_max_range.\n
 
 We assume that Fiber::adjustSegmentation() is used, ensuring that
 ( actual segmentation ) < ( 4/3 * FiberProp::segmentation ).
 */
real Simul::estimateStericRange() const
{
    real ran = 0;
    real len = 0;
    
    ///check all FiberProp with enabled steric:
    PropertyList plist = properties.find_all("fiber");
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        FiberProp const* fp = static_cast<FiberProp*>(*ip);
        if ( fp->steric )
        {
            // The maximum length of a segment is 4/3 * segmentation
            real x = 4.1/3 * fp->segmentation;
            if ( len < x )
                len = x;
            
            // check extended range of interaction
            x = fp->steric_radius + fp->steric_range;
            if ( ran < x )
                ran = x;
        }
    }
    
    ///verify against the actual segmentations of the Fibers:
    for ( Fiber const* fib=fibers.first(); fib; fib=fib->next() )
    {
        if ( fib->prop->steric )
        {            
            // also check the actual segmentation
            real x = fib->segmentation();
            if ( len < x )
                len = x;
        }
    }

    /*
     The interaction can be aligned with the fiber, and we must add the distances:
     2 * range if two fibers of radius 'range' interact.
     + 2 * ( len / 2 ) since len/2 is the distance between the center of the segment
     and its most distal point.
     */
    ran =  len + 2*ran;
    
    
    for ( Sphere const* sp=spheres.first(); sp; sp=sp->next() )
    {
        real d = 2 * sp->radius() + sp->prop->steric_range;
        if ( sp->prop->steric && ran < d )
            ran = d;
    }
    
    for ( Bead const* bd=beads.first(); bd; bd=bd->next() )
    {
        real d = 2 * bd->radius() + bd->prop->steric_range;
        if ( bd->prop->steric && ran < d )
            ran = d;
    }
    
    for ( Solid const* so=solids.first(); so; so=so->next() )
    {
        if ( so->prop->steric )
        {
            for ( unsigned p = 0; p < so->nbPoints(); ++p )
            {
                real d = 2 * so->radius(p) + so->prop->steric_range;
                if ( ran < d )
                    ran = d;
            }
        }
    }
    
    if ( ran < REAL_EPSILON )
        PRINT_ONCE("Warning: could not estimate simul:steric_max_range automatically!\n");
    
    return ran;
}




void Simul::setStericGrid(Space const* spc) const
{
    assert_true(spc);
    real& range = prop->steric_max_range;
    
    if ( range <= 0 ) 
    {
        range = estimateStericRange();
        //MSG("auto setting simul:steric_max_range=%.3f\n", range);
    }
    
    if ( range > 0 ) 
        stericGrid.setGrid(spc, modulo, range);
    else
        throw InvalidParameter("simul:steric is enabled, but simul:steric_max_range was not set");
}


/**
 The prop->steric of each object is a bit-field that
 specify one or more 'pane' where the object is present.
 The different panes are then treated consecutively and independently, 
 and only objects in the same pane may interact.
 
 @code
 for ( int pane=1; pane<=2 && pane<=prop->steric; ++pane )
 {
     if ( obj->prop->steric & pane )
     ...
 }
 @endcode
 
 With this mechanism, the user can flexibly configure which objects
 may see each other and thus control the steric interactions.
 
 At present, we only support 1 pane (Simul property steric).
 This can be extended if necessary, but the steric_stiffness[]
 properties should be extended as well.
 */
void Simul::setStericInteractions(Meca& meca) const
{
    if ( !stericGrid.hasGrid() )
        setStericGrid(sSpace);

    // clear grid
    stericGrid.clear();
    
    // distribute Fiber-points on the grid
    for ( Fiber* fib=fibers.first(); fib; fib=fib->next() )
    {
        if ( fib->prop->steric )
        {
            const real rad = fib->prop->steric_radius;        // equilibrium radius
            const real ran = rad + fib->prop->steric_range;   // extended range of interaction
        
            // include segments, in the cell associated with their center
            for ( unsigned r = 0; r < fib->nbSegments(); ++r )
#if ( NB_STERIC_PANES == 1 )
                stericGrid.add(FiberLocus(fib, r), rad, ran);
#else
                stericGrid.add(fib->prop->steric, FiberLocus(fib, r), rad, ran);
#endif
        }
    }
    
    // include Spheres
    for ( Sphere* sp=spheres.first(); sp; sp=sp->next() )
    {
        if ( sp->prop->steric )
#if ( NB_STERIC_PANES == 1 )
            stericGrid.add(PointExact(sp, 0), sp->radius(), sp->radius()+sp->prop->steric_range);
#else
            stericGrid.add(sp->prop->steric, PointExact(sp, 0), sp->radius(), sp->radius()+sp->prop->steric_range);
#endif
    }
    
    // include Beads
    for ( Bead* bd=beads.first(); bd; bd=bd->next() )
    {
        if ( bd->prop->steric )
#if ( NB_STERIC_PANES == 1 )
            stericGrid.add(PointExact(bd, 0), bd->radius(), bd->radius()+bd->prop->steric_range);
#else
            stericGrid.add(bd->prop->steric, PointExact(bd, 0), bd->radius(), bd->radius()+bd->prop->steric_range);
#endif
    }
        
    // include Points that have a radius from Solids
    for ( Solid* so=solids.first(); so; so=so->next() )
    {
        if ( so->prop->steric )
        {
            for ( unsigned int pp = 0; pp < so->nbPoints(); ++pp )
            {
                if ( so->radius(pp) > REAL_EPSILON )
#if ( NB_STERIC_PANES == 1 )
                    stericGrid.add(PointExact(so, pp), so->radius(pp), so->radius(pp)+so->prop->steric_range);
#else
                    stericGrid.add(so->prop->steric, PointExact(so, pp), so->radius(pp), so->radius(pp)+so->prop->steric_range);
#endif
            }
        }
    }
    
    /// create parameters
    PointGridParam pam(prop->steric_stiffness_push[0], prop->steric_stiffness_pull[0]);
    
#if ( NB_STERIC_PANES == 1 )
    
    stericGrid.setInteractions(meca, pam);

#elif ( NB_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    stericGrid.setInteractions(meca, pam, 1);
    // add steric interactions between panes 1 and 2:
    stericGrid.setInteractions(meca, pam, 1, 2);
    //stericGrid.setInteractions(meca, pam, 2, 1);

#else
    
    // add steric interactions between different panes:
    for ( unsigned p = 1; p <= NB_STERIC_PANES; ++p )
        stericGrid.setInteractions(meca, pam, p);

#endif
}


//------------------------------------------------------------------------------
/**
 This will:
 - Register all Mecables in the Meca: Fiber Solid Bead and Sphere
 - call setInteractions() for all objects in the system,
 - call setStericInteractions() if prop->steric is true.
 .
 */
void Simul::setInteractions(Meca & meca) const
{
    // prepare the meca, and register Mecables
    meca.clear();
    
    for ( Fiber  * mt= fibers.first(); mt ; mt=mt->next() )
        meca.add(mt);
    for ( Solid  * so= solids.first(); so ; so=so->next() )
        meca.add(so);
    for ( Bead   * bd=  beads.first(); bd ; bd=bd->next() )
        meca.add(bd);
    for ( Sphere * se=spheres.first(); se ; se=se->next() )
        meca.add(se);
    
    meca.prepare(prop);
    
    // add interactions
    
    for ( Space * sp=spaces.first(); sp; sp=sp->next() )
        sp->setInteractions(meca, fibers);
    
    for ( Fiber * mt=fibers.first(); mt ; mt=mt->next() )
        mt->setInteractions(meca);
    
    for ( Solid * so=solids.first(); so ; so=so->next() )
        so->setInteractions(meca);
    
    for ( Bead * bd=beads.first(); bd ; bd=bd->next() )
        bd->setInteractions(meca);
    
#if ( 0 )
    PRINT_ONCE("AD-HOC STRING FORCE IS ENABLED");
    // attach beads together in a string:
    for ( Bead * bd=beads.first(); bd ; bd=bd->next() )
    {
        Bead * nx = beads.find(bd->identity()+1);
        if ( nx )
        {
            real d = bd->radius() + nx->radius();
            meca.interLongLink(PointExact(bd,0), PointExact(nx, 0), d, 100);
        }
    }
#endif

    for ( Sphere * sp=spheres.first(); sp ; sp=sp->next() )
        sp->setInteractions(meca);
    
#if ( 0 )
    PRINT_ONCE("AD-HOC REPULSIVE FORCE IS ENABLED");
    // add pairwise repulsive force between Spheress:
    for ( Sphere * sp=spheres.first(); sp ; sp=sp->next() )
    for ( Sphere * ss=sp->next(); ss ; ss=ss->next() )
        meca.interCoulomb(PointExact(ss,0), PointExact(sp, 0), 20);
#endif

    for ( Single * si=singles.firstA(); si ; si=si->next() )
        si->setInteractions(meca);


    for ( Couple * cx=couples.firstAA(); cx ; cx=cx->next() )
        cx->setInteractions(meca);

//Removed by Manu, because we need it for the trapper
//#ifdef NEW_DANGEROUS_CONFINEMENTS
    for ( Couple * cx=couples.firstAF(); cx ; cx=cx->next() )
        cx->setInteractionsAF(meca);

    for ( Couple * cx=couples.firstFA(); cx ; cx=cx->next() )
        cx->setInteractionsFA(meca);
//#endif
    
    for ( Organizer * as = organizers.first(); as; as=as->next() )
        as->setInteractions(meca);
    
    // add steric interactions
    if ( prop->steric && sSpace )
        setStericInteractions(meca);
    
#if ( 0 )
    PRINT_ONCE("AD-HOC CALIBRATED FORCE IS ENABLED");
    // add calibrated forces, for testing rotation
    for ( Fiber * fib = fibers.first(); fib; fib = fib->next() )
        sMeca.addTorqueClamp(fib->interpolateCenter(), Vector(0,1,0), 1);
#endif
#if ( 0 )
    PRINT_ONCE("AD-HOC CALIBRATED FORCE IS ENABLED");
    // add calibrated force to test rotation of spheres:
    Vector force(0,1,0);
    for ( Sphere * sph = spheres.first(); sph; sph = sph->next() )
    {
        sMeca.addForce(PointExact(sph, 1), -force);
        sMeca.addForce(PointExact(sph, 2), +force);
    }
#endif
}

//------------------------------------------------------------------------------
void Simul::solve()
{
    setInteractions(sMeca);
    
    if ( prop->precondition & 7 )
        precond = prop->precondition & 7;
    
    clock_t cpu = clock();
    
    // solve the system:
    sMeca.solve(prop, precond);
    
    // Automatic selection of preconditionning method:
    if ( prop->precondition == 8 )
    {
        const unsigned N_TEST = 6;  //this must be a multiple of 6
        const unsigned PERIOD = 256;
        
        //automatically select the preconditionning mode:
        //by trying each methods N_STEP steps, adding CPU time and use fastest.
        if ( precondCounter++ < N_TEST )
        {
            precondCPU[precond] += clock() - cpu;
            
            //std::clog << " precond "<<precond<<" iter " << iter << " CPU " << cpu << "\n";
            
            if ( precondCounter == N_TEST )
            {
                // if the differential of times is significant, use the fastest method
                // but otherwise, select the simplest method:
                if ( precondCPU[0] > precondCPU[1] + 32 )
                    precond = 1;
                else
                    precond = 0;

                if ( prop->verbose )
                {
                    std::clog << " precond 0 sum " << precondCPU[0] << "\n";
                    std::clog << "         1 sum " << precondCPU[1] << "\n";
                    std::clog << "         2 sum " << precondCPU[2] << "\n";
                    std::clog << " ----> " << precond << "\n";
                }
            }
            else
            {
                //alternate betwen methods
                precond = ( 1 + precond ) % 2;
            }
        }
        else if ( precondCounter > PERIOD )
        {
            precondCPU[0] = 0;
            precondCPU[1] = 0;
            precondCPU[2] = 0;
            precondCPU[3] = 0;
            precondCounter = 0;
        }
    }
}


void Simul::computeForces() const
{
    //Meca mec;
    try {
        prop->complete(this);
        setInteractions(sMeca);
        sMeca.computeForces();
    }
    catch ( Exception & e )
    {
        std::clog << "cytosim could not compute the forces:\n";
        std::clog << "   " << e.what() << std::endl;
    }
}

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "fiber.h"
#include "messages.h"
#include "glossary.h"
#include "iowrapper.h"
#include "fiber_locus.h"
#include "fiber_binder.h"
#include "fiber_prop.h"
#include "simul_prop.h"
#include "simul.h"
#include "sim.h"
#include <algorithm>

//#define OLD_FIBER_SQUEEZE
#define NEW_COLINEAR_FORCE
//#define NEW_SEVER_KINKED_FIBERS

//------------------------------------------------------------------------------

void Fiber::step()
{
    //add single that act like glue
    if ( prop->glue )
    {
        setGlue(frGlue, PLUS_END, prop->confine_space_ptr, prop->glue);
    }

    if ( frLattice )
    {
        /// get address of array of sites:
        FiberLattice::site_type * site = frLattice->data();

    if ( prop->lattice_binding_rate > 0 )
    {
        Field * field = prop->field_ptr;
        // we want roughly one point per cell:
        const real spread = field->cellWidth();
        // each point represents a Fiber chunk of length 'spread':
        const real rate = prop->lattice_binding_rate * spread / field->cellVolume();
        // fraction of the cell content that will bind in one time_step:
        const real frac = 1 - exp( -rate * prop->time_step );
        
        real abs = spread * RNG.exponential();
        const real len = length();
        // stochastic sampling with distance 'spread' along the Fiber:
        while ( abs < len )
        {
            Vector pos = posM(abs);
            real& cell = field->cell(pos);
            assert_true( cell >= 0 );
            
            // amount to be transfered:
            real mass = cell * frac;

            cell -= mass;
            frLattice->site(abs+abscissaM()) += mass;
            
            abs += spread * RNG.exponential();
        }
    }

    if ( prop->lattice_flux_speed )
    {
        assert_true(prop->field_ptr);
        const real fac = prop->lattice_flux_speed * prop->time_step / prop->lattice_unit;
        
        if ( fabs(fac) > 1 )
            throw InvalidParameter("lattice_flux_speed * time_step / lattice_unit is too high");
        
        const int inf = frLattice->index(abscissaM());
        const int sup = frLattice->index(abscissaP());
        
        if ( fac < 0 )
        {
            real s = site[inf];
            
            for ( int h = inf; h < sup; ++h )
                site[h] -= fac * ( site[h+1] - site[h] );
            
            prop->field_ptr->cell(posEndM()) -= fac * s;
            site[sup] += fac * site[sup];
        }
        else
        {
            real s = site[sup];
            
            for ( int h = sup; h > inf; --h )
                site[h] -= fac * ( site[h] - site[h-1] );
            
            prop->field_ptr->cell(posEndP()) += fac * s;
            site[inf] -= fac * site[inf];
        }
        
        //std::clog << "lattice sum = " << frLattice->sum() << "   ";
    }

    if ( prop->lattice_unbinding_rate > 0 )
    {
        releaseLattice(prop->field_ptr, prop->lattice_unbinding_prob);
    }

    if ( prop->lattice_cut_fiber )
    {
        const real fac = 1.0 / prop->time_step;
        
        const real uni = frLattice->unit();
        const int  inf = frLattice->index(abscissaM());
        const int  sup = frLattice->index(abscissaP());
        
        assert_true( inf >= frLattice->inf() );
        assert_true( sup <  frLattice->sup() );
        assert_true( frLattice->abscissa(inf) <= abscissaM() );
        assert_true( frLattice->abscissa(sup+1) >= abscissaP() );
        
        int h = inf;
        real val = fac * RNG.exponential();

        real ai  = abscissaM();
        real as  = frLattice->abscissa(h+1);
        if ( as > abscissaP() )
            as = abscissaP();
        
        while ( h <= sup )
        {
            assert_true( site[h] >= 0 );
            val -= site[h];
            //assert_true( ai >= abscissaM() );
            //assert_true( as <= abscissaP() );

            while ( val < 0 )
            {
                /*
                 Since val < 0 and 0 <= val+site[h],
                 then -val/site[h] <= 1
                 hence   0 < -val/site[h] < 1
                */
                real abs = ai - ( as - ai ) * val / site[h];
                
                assert_true( abs >= ai - REAL_EPSILON );
                assert_true( abs <= as + REAL_EPSILON );
                
                sever(abs, STATE_RED, STATE_GREEN);
                val += fac * RNG.exponential();
            }
            
            ai = as;
            if ( ++h == sup )
                as = abscissaP();
            else
                as += uni;
        }
    }
    }
    
#ifdef NEW_CHEW_FIBERS
    if ( frChewP > 0 )
    {
        if ( frChewP > prop->max_chewing_speed_dt )
            frChewP = prop->max_chewing_speed_dt;
        
        if ( length() - frChewM - frChewP < prop->min_length )
        {
            delete(this);
            return;
        }
        growP(-frChewP);
        frChewP = 0;
    }
    if ( frChewM > 0 )
    {
        if ( frChewM > prop->max_chewing_speed_dt )
            frChewM = prop->max_chewing_speed_dt;

        if ( length() - frChewM < prop->min_length )
        {
            delete(this);
            return;
        }
        growM(-frChewM);
        frChewM = 0;
    }
#endif
    
#ifdef NEW_SEVER_KINKED_FIBERS
    /**
     Cut the filaments at any position where two segments
     make an angle above PI/2
     */
    PRINT_ONCE("NEW_SEVER_KINKED_FIBERS\n");
    severKinks();
#endif

    // perform the cuts that were registered by sever()
    if ( futureCuts.size() )
        severNow();
    
    // delete self if shorter than 'FiberProp::min_length'
    if ( length() < prop->min_length )
    {
        delete(this);
        return;
    }
    
    if ( needUpdate )
    {
        update();
        needUpdate = false;
    }
}


//------------------------------------------------------------------------------

/**
 The rest of the initialization is done in FiberProp::newFiber(),
 and other newFiber() functions where the initial length is known.
 */
Fiber::Fiber(FiberProp const* p)
: frLattice(0), frGlue(0), frBirthTime(0), frChewM(0), frChewP(0), prop(p), disp(0)
{
#ifdef MULTI_LATTICE
    // This seems to be a most correct initialization of an empty array
    for (int i = 0; i<4; i++) {
        multi_lattice[i]=nullptr;
    }
#endif
    if ( prop )
    {
        segmentation(prop->segmentation);
        
        if ( prop->lattice )
        {
            //std::clog << reference() <<  " new Lattice\n";
            frLattice = new FiberLattice(prop->lattice_unit);
        }
    }
}


Fiber::~Fiber()
{
    detachBinders();
    
    if ( frLattice )
    {
        if ( prop->field_ptr )
            releaseLattice(prop->field_ptr);
        delete frLattice;
        frLattice = 0;
    }
    
    // if linked in SingleSet, frGlue could be deleted twice
    // when the simulation ends and all objects are deleted. 
    if ( frGlue )
    {
        if ( frGlue->linked() )
            prop->glue_set->remove(frGlue);
        delete(frGlue);
        frGlue = 0;
    }
    
    if ( disp )
    {
        /*
         Note: the destructor will not be called here, which is OK
         if LineDisp is a trivial type that does not allocate resources
         */
        free(disp);
        disp = 0;
    }
    
    prop = 0;
}


real Fiber::age() const
{
    return objset()->simul.time() - frBirthTime;
}


unsigned Fiber::allocatePoints(const unsigned nbp)
{
    unsigned ms = RigidFiber::allocatePoints(nbp);
    /*
     if RigidFiber::allocatePoints() allocated memory, it will return the 
     size of the new array, and we allocate the same size for other arrays.
     */
    if ( ms )
    {
        //std::clog << "Fiber::allocatePoints " << ms << std::endl;
        frLocuses.resize(ms);
        for ( unsigned ii = 0; ii < ms-1; ++ii )
            frLocuses[ii] = FiberLocus(this, ii);
    }
    return ms;
}


FiberLocus & Fiber::locus(const unsigned int pos) const
{
    assert_true( pos < frLocuses.size() );
    assert_true( frLocuses[pos].fiber() == this );
    
    return frLocuses[pos];
}


real Fiber::projectPoint(Vector const& w, real & dis) const
{
    // initialize with the minus-end:
    dis = w.distanceSqr(posP(0));
    real abs = 0, len = segmentation();
    
    // try all segments
    for ( unsigned int ii = 0; ii < nbSegments(); ++ii )
    {
        //check the segment:
        FiberLocus s(this, ii);
        real d = INFINITY;
        real a = s.projectPoint0(w, d);
        if ( len < a )
        {
            // test exact point
            real e = w.distanceSqr(posP(ii+1));
            if ( e < dis ) {
                abs = abscissaPoint(ii+1);
                dis = e;
            }
        }
        else if ( 0 <= a  &&  d < dis )
        {
            //the projection is the best found so far
            abs = abscissaPoint(ii) + a;
            dis = d;
        }
    }
    
    return abs;
}

//------------------------------------------------------------------------------
#pragma mark -

void Fiber::flip()
{
    // flip all the points:
    Filament::flip();
    
    /* update abscissa of Binders to keep them in place:
     new - minus = plus - old
     new = plus + minus - old
     */
    real abs = abscissaM() + abscissaP();
    Node * hi = frBinders.front();
    while ( hi )
    {
        FiberBinder * ha = static_cast<FiberBinder*>(hi);
        hi = hi->next();
        ha->relocate(abs-ha->abscissa());
    }
}


/**
 A portion of size `len` that includes the MINUS_END is removed.
 The FiberBinders bound within the deleted portion are detached.
 */
void Fiber::cutM(real len)
{
    real abs = abscissaM() + len;
    
    Filament::cutM(len);
    
    //update FiberBinders:
    Node * hi = frBinders.front();
    while ( hi )
    {
        Hand * ha = static_cast<Hand*>(hi);
        hi = hi->next();
        if ( ha->abscissa() < abs )
            ha->detach();
    }
}


/**
 A portion of size `len` that includes the PLUS_END is removed.
 The FiberBinders bound within the deleted portion are detached.
 */
void Fiber::cutP(real len)
{
    real abs = abscissaP() - len;
    
    Filament::cutP(len);
    
    //update FiberBinders:
    Node * hi = frBinders.front();
    while ( hi )
    {
        Hand * ha = static_cast<Hand*>(hi);
        hi = hi->next();
        if ( ha->abscissa() > abs )
            ha->detach();
    }
}



/**
 The Fiber is cut at point P
 - A new Fiber is created from the section [ P , PLUS_END ],
 - all FiberBinder attached to this section are transferred to the new Fiber,
 - Lattice content is also transferred,
 - a pointer to the new Fiber is returned, which should be added to the Simul
 .
 @return zero, if `pti` is not an internal point
 */
Fiber* Fiber::severPoint(unsigned int pti)
{
    if ( pti == 0  ||  pti >= lastPoint() )
        return 0;
    
    real abs = abscissaPoint(pti);

    // create a new Fiber of the same kind:
    Fiber* fib = prop->newFiber();
    assert_true( fib->prop == prop );
    
    // copy the Filament part of the object:
    *(static_cast<Filament*>(fib)) = *this;
    
    // the signature on both pieces should be conserved:
    fib->signature(signature());
    fib->birthTime(birthTime());

    assert_true( fib->abscissaM() == abscissaM() );
    // remove MINUS_END portion on new piece:
    fib->truncateM(pti);
    assert_true(fib->abscissaM() == abs);
    
    if ( lattice() )
    {
        assert_true( fib->lattice() );
        assert_true( lattice()->unit() == fib->lattice()->unit() );
        
        // transfer Lattice values located above the cut
        fib->lattice()->takeP(lattice(), ceil( abs / lattice()->unit() ));
    }
    
    // remove PLUS_END portion on self
    truncateP(pti);
    
    // transfer FiberBinders above point P
    // their abscissa should not change in the transfer
    Node * nd = frBinders.front();
    while( nd )
    {
        FiberBinder* ha = static_cast<FiberBinder*>(nd);
        nd = nd->next();
        if ( ha->abscissa() > abs )
            ha->relocate(fib);
        else
            ha->updateBinder();
    }
    
    return fib;
}


/**
The Fiber is cut at distance `abs` from its MINUS_END:
 - current Fiber is truncated to keep only the section [ MINUS_END , abs ],
 - A new Fiber is created inheriting the other section [ abs , PLUS_END ],
 - FiberBinder within the severed section are transfered to the new Fiber,
 - lattice substances are also transfered,
 .
 A pointer to the new Fiber is returned (containing the PLUS_END), but this
 pointer may be zero, if `abs` was not within the valid range of abscissa.
 If a new Fiber was created, it should be added to the FiberSet.
 */
Fiber* Fiber::severM(real abs)
{
    if ( abs <= REAL_EPSILON || abs + REAL_EPSILON >= length() )
        return 0;
    
    // create a new Fiber of the same kind:
    Fiber* fib = prop->newFiber();
    assert_true( fib->prop == prop );

    // copy the Filament part of the object:
    *(static_cast<Filament*>(fib)) = *this;
    
    // the signature on both pieces should be conserved:
    fib->signature(signature());
    fib->birthTime(birthTime());

    assert_small(fib->abscissaM() - abscissaM());
    // remove MINUS_END portion on new piece
    fib->Filament::cutM(abs);

    if ( lattice() )
    {
        assert_true( fib->lattice() );
        assert_true( lattice()->unit() == fib->lattice()->unit() );
        
        // ensure valid range:
        fib->lattice()->setRange(abscissaM()+abs, abscissaP());
  
        // transfer Lattice values located above the cut:
        fib->lattice()->takeP(lattice(), ceil( (abscissaM()+abs) / lattice()->unit() ));
    }
    
    assert_small(fib->abscissaM()-abs-abscissaM());
    
    // remove PLUS_END portion on self
    Filament::cutP(length()-abs);
    
    assert_small(fib->abscissaM()-abscissaP());

    // transfer all FiberBinders above cut to new piece
    // their abscissa should not change in this transfer
    const real edge = abs + abscissaM();
    Node * nd = frBinders.front();
    while( nd )
    {
        FiberBinder* ha = static_cast<FiberBinder*>(nd);
        nd = nd->next();
        if ( ha->abscissa() > edge )
            ha->relocate(fib);
        else
            ha->updateBinder();
    }

    return fib;
}


/**
 Perform the cuts that were registered in `futureCuts` by Fiber::sever(),
 and clear `futureCuts`.
 This deletes Fibers that are shorter than FiberProp::min_length
*/
void Fiber::severNow()
{
    /**
     The std::set keeps its objects always in order of descending abscissa,
     which is essential here to avoid data loss.
     */
    
    // cut starting from highest abscissa
    for ( std::set<SeverPos>::const_iterator cut = futureCuts.begin(); cut != futureCuts.end(); ++cut )
    {
        if ( cut->abs - abscissaM() <= prop->min_length )
        {
            if ( cut->abs <= abscissaM() )
                std::cerr << " severed abscissa " << cut->abs << " is below MINUS_END " << abscissaM() << std::endl;
            else
                cutM(cut->abs-abscissaM());
            /*
             because we have deleted the MINUS_END section, 
             and following cuts in the list will be of lower abscissa,
             they should not be processed.
             */
            break;
        }
        else if ( abscissaP() - cut->abs <= prop->min_length )
        {
            if ( cut->abs >= abscissaP() )
                std::cerr << " severed abscissa " << cut->abs << " is above PLUS_END " << abscissaP() << std::endl;
            else
                cutP(abscissaP()-cut->abs);
        }
        else
        {
            Fiber * frag = severM(cut->abs-abscissaM());
            
            // special case where the PLUS_END section is simply deleted
            if ( cut->staM == STATE_BLACK )
            {
                delete(frag);
                continue;
            }

            //add new fragment to simulation:
            objset()->add(frag);

            if ( frag )
            {
                // check that ends spatially match:
                assert_small((frag->posEndM() - posEndP()).norm());
                
                try {
                    // old PLUS_END converves its state:
                    frag->setDynamicStateP(dynamicStateP());
                    
                    // new ends are set as wished:
                    this->setDynamicStateP(cut->staP);
                    frag->setDynamicStateM(cut->staM);
                }
                catch ( Exception & e )
                {
                    e << ", when cutting fiber " << reference();
                    throw;
                }
            
#ifdef LOGGING
                std::clog << "severed " << reference() << " at abscissa " << s->abs;
                std::clog << "   creating " << frag->reference();
                std::clog << "   position " << frag->posEndM() << std::endl;
#endif
                //std::clog << " severed at X = " << frag->posEndM().XX << std::endl;
            }
            else
            {
                std::clog << " sever abscissa " << cut->abs << " is out of range";
                std::clog << " [ " << abscissaM() << "   " << abscissaP() << " ]" << std::endl;
            }
        }
    }
    futureCuts.clear();
}


/**
 Cut the fiber if the angle made by consecutive segments is acute ( cosine < 0 ).
 The Fiber may be severed multiple times, at every model-points where there
 is a kink. Fiber parts are added to the FiberSet to which the current Fiber belongs.
 */
void Fiber::severKinks()
{
    /*
     We must consider points in reverse order,
     because severPoint() removes the distal part
    */
    for ( int p = nbPoints()-2; p > 0 ; --p )
    {
        if ( diffPoints(p-1) * diffPoints(p) < 0 )
            objset()->add(severPoint(p));
    }
}


void Fiber::planarCut(Vector const& n, const real a, int stateP, int stateM)
{
    Array<real> cuts;
    
    /*
     The cuts should be processed in order of decreasing abscissa,
     hence we check intersections from PLUS_END to MINUS_END
    */
    for ( int s = lastSegment(); s >=0 ; --s )
    {
        real abs = planarIntersect(s, n, a);
        if ( 0 <= abs  &&  abs < 1 )
            cuts.push_back(abscissaPoint(s+abs));
    }
    
    for ( real * s = cuts.begin(); s < cuts.end(); ++s )
    {
        Fiber * fib = severNow(*s);
        if ( fib )
        {
            // old PLUS_END converves its state:
            fib->setDynamicStateP(dynamicStateP());
            // dynamic of new ends are set as usual:
            setDynamicStateP(stateP);
            fib->setDynamicStateM(stateM);
            objset()->add(fib);
        }
    }
}


/**
 The `fib` is added past the PLUS_END of `*this`,
 The dynamic state at PLUS_END is set to fib->dynamicStateP(),
 FiberBinders are transfered.
 */
void Fiber::join(Fiber * fib)
{
    assert_true( fib );
    // the two fibers should be of the same kind:
    assert_true( prop == fib->prop );
    
    // shift in abscissa must be calculated before joining
    real shift = abscissaP() - fib->abscissaM();

    // join backbones
    Filament::join(fib);
    
    //transfer dynamic state of PLUS_END:
    setDynamicStateP(fib->dynamicStateP());

    // transfer all FiberBinder
    Node * nd = fib->frBinders.front();
    while ( nd )
    {
        FiberBinder* ha = static_cast<FiberBinder*>(nd);
        nd = nd->next();
        ha->relocate(this, ha->abscissa()+shift);
    }
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 From "Random Walks in Biology" by HC. Berg, Princeton University Press,
 drag coefficients for an ellipsoid are,

 @code
 drag_transverse = 2*drag_parallel = 4*PI*L*visc / log(length/radius)
 @endcode

 We should average the mobility coefficients:  speed = mu * f
 mu_X = mu_parallel   = 2 * mu
 mu_Y = mu_transverse = mu
 mu_Z = mu_transverse = mu
 Hence:
 mu_averaged = ( mu + mu + 2*mu ) / 3 = 4/3 * mu.
 drag_averaged = 3*PI*length*viscosity / log(length/radius)

APPROXIMATE FORMULA FOR ELLIPSOIDAL PARTICLE
Clift R, Grace JR, Weber ME. Bubbles, drops, and particles: Courier Corporation; 2005.

 @code
 aspect = length / diameter;
 drag = 3.0 * M_PI * viscosity * diameter * ( 3 + 2 * length/diameter ) / 5.0;
 @endcode
 */



/** 
 Fiber::setDragCoefficientVolume() calculates the mobility for the entire fiber, 
 considering that the cylinder is straight and moving in a infinite fluid.
 fiber:hydrodynamic_radius[1] is a hydrodynamic cutoff that makes the
 drag coefficient proportional to length beyond the cutoff.
 
 The formula for a cylinder is taken from:\n
 <em>
 Tirado and de la Torre. J. Chem. Phys 71(6) 1979 \n
 http://link.aip.org/link/doi/10.1063/1.438613 \n
 </em>

 We calculate the translational drag coefficient averaged over all possible configurations:
 @code
   aspect = length / diameter;
   Ct =  0.312 + 0.565/aspect - 0.100/(aspect*aspect);
   drag_cylinder = 3*M_PI*viscosity*length / ( log(aspect) + Ct );
 @endcode
 
 The rotational diffusion coefficient is given by:
 @code
   Cr = -0.662 + 0.917/aspect - 0.050/(aspect*aspect);
   drag_rotation = M_PI*viscosity*length / ( log(aspect) + Cr )
 @endcode
 
 If the length is shorter than the diameter, the formula above fails and may even give negative result.
 Hence we also calculate the drag of a sphere with the same radius as the cylinder:
 @code
   drag_sphere = 6*PI*visc*R
 @endcode
 We use the maximum value between 'drag_sphere' and 'drag_cylinder'.
 */
real Fiber::dragCoefficientVolume()
{
    real len = length();
    assert_true( len > 0 );
    
    // hydrodynamic cut-off on length:
    real lenc = len;
    assert_true( prop->hydrodynamic_radius[1] > 0 );
    
    if ( lenc > prop->hydrodynamic_radius[1] )
        lenc = prop->hydrodynamic_radius[1];
    
    if ( lenc < prop->hydrodynamic_radius[0] )
        lenc = prop->hydrodynamic_radius[0];
    
    //Stoke's for a sphere:
    assert_true( prop->hydrodynamic_radius[0] > 0 );
    real drag_sphere = 6 * prop->hydrodynamic_radius[0];
    
#ifdef NEW_ANISOTROPIC_FIBER_DRAG
    const real pref = 4;
#else
    const real pref = 3;
#endif

#if ( 0 )
    /*
     For an ellipsoid,  
     drag_transverse = 2*drag_parallel = 4*PI*L*visc / log(length/radius)
     We should average the mobility coefficients:  speed = mu * f
       mu_X = mu_parallel   = 2 * mu
       mu_Y = mu_transverse = mu
       mu_Z = mu_transverse = mu
     Hence:
       mu_averaged = ( mu + mu + 2*mu ) / 3 = 4/3 * mu.
     drag_averaged = 3*PI*length*viscosity / log(length/radius)
     See for example "Random Walks in Biology" by HC. Berg, Princeton University Press.
     */
    
    // length below which the formula is not valid:
    real min_len = exp( 1 + log(prop->hydrodynamic_radius[0]) );

    real drag_cylinder = pref * len / log( lenc / prop->hydrodynamic_radius[0] );
#else
    /*
     Tirado and de la Torre. J. Chem. Phys 71(6) 1979
     give the averaged translational friction coefficient for a cylinder:
     3*PI*length*viscosity / ( log( length/diameter ) + 0.312 )
     */
    
    // length below which the formula is not valid:
    real min_len = exp( 1 - 0.312 + log(2*prop->hydrodynamic_radius[0]) );

    real drag_cylinder = pref * len / ( log( 0.5 * lenc / prop->hydrodynamic_radius[0] ) + 0.312 );
#endif

    real drag;
    
    if ( len < min_len )
        drag = drag_sphere;
    else
    {
        // use largest drag coefficient
        drag = ( drag_cylinder > drag_sphere ? drag_cylinder : drag_sphere );
    }
    
    //MSG("Drag coefficient of Fiber in infinite fluid = %.1e\n", drag);
    //std::clog << "Fiber L = " << len << "  bulk_drag = " << drag << std::endl;;

    return M_PI * prop->viscosity * drag;
}



/**
 Fiber::setDragCoefficientSurface() uses a formula calculated by F. Gittes in:\n
 <em>
 Hunt et al. Biophysical Journal (1994) v 67 pp 766-781 \n
 http://dx.doi.org/10.1016/S0006-3495(94)80537-5 \n
 </em>
 
 It applies to a cylinder moving parallel to its axis and near an immobile surface:
 @code
   drag_per_unit_length = 2 &pi &eta / acosh(h/r)
 @endcode
 
 With:
 - r = cylinder radius,
 - h = distance between cylinder bottom and surface,
 - &eta = viscosity of the fluid.
 
 If the cylinder is exactly touching the surface, `h=0` and the drag coefficient is infinite.
 
 The drag coefficient for motion perpendicular to the cylinder axis would be twice higher,
 but for gliding assays, the parallel drag coefficient is the appropriate choice.  
 
 Note that this is usually equivalent to the approximate formula:
 @code
   drag_per_unit_length = 2 &pi &eta / log(2*h/r)
 @endcode
 because 
 @code
   acosh(x) = ln[ x + sqrt(x^2-1)) ] ~ ln[2x] if x >> 1
 @endcode

 Hunt et al. also credit this reference for the formula:\n
 <em>
 The slow motion of a cylinder next to a plane wall.
 Jeffrey, D.J. & Onishi, Y. (1981) Quant. J. Mech. Appl. Math. 34, 129-137.
 </em>
*/
real Fiber::dragCoefficientSurface()
{
    real len = length();    
    
    if ( prop->cylinder_height <= 0 )
        throw InvalidParameter("fiber:surface_effect[1] (height above surface) must set and > 0!");
    
    // use the higher drag: perpendicular to the cylinder (factor 2)
    real drag = 2 * M_PI * prop->viscosity * len / acosh( 1 + prop->cylinder_height/prop->hydrodynamic_radius[0] );
    
    //MSG("Drag coefficient of Fiber near a planar surface = %.1e\n", drag);
    //std::clog << "Drag coefficient of Fiber near a planar surface = " << drag << std::endl;

    return drag;
}



/**
 Calculate drag coefficient from two possible formulas
@code
 if ( fiber:surface_effect )
    drag = dragCoefficientSurface();
 else
    drag = dragCoefficientVolume();
@endcode
 */
void Fiber::setDragCoefficient()
{
    real drag;
    
    if ( prop->surface_effect )
    {
        drag = dragCoefficientSurface();
        if ( 0 )
        {
            real d = dragCoefficientVolume();
            std::clog << "Drag coefficient of Fiber near a planar surface amplified by " << drag/d << std::endl;
        }
    }
    else
        drag = dragCoefficientVolume();

    //the forces are distributed equally on all points, hence we multiply by nbPoints()
    assert_true( nbPoints() > 0 );
    rfDragPoint = drag / nbPoints();
    
#if ( 0 )
    std::clog << "Fiber L = " << std::setw(7) << length();
    std::clog << " drag = " << drag << " drag_point " << rfDragPoint << std::endl;
#endif
}


void Fiber::prepareMecable()
{
    setDragCoefficient();
    storeDirections();
    makeProjection();

    assert_true( rfDragPoint > REAL_EPSILON );

    // the scaling of the bending elasticity depends on the length of the segments
    rfRigidity = prop->rigidity / segmentationCub();
    
#if ( 0 )
    real energy = bendingEnergy();
    real euler = M_PI * M_PI * prop->rigidity / ( length() * length() );
    std::clog << "Euler buckling = " << euler << "    ";
    std::clog << "Bending energy = " << energy << std::endl;
#endif
}


//------------------------------------------------------------------------------

void Fiber::setInteractions(Meca & meca) const
{
#ifdef OLD_FIBER_SQUEEZE
    if ( prop->squeeze == 1 )
    {
        // squeezing force in the YZ-plane:
        const real f = prop->squeeze_force;
        const real r = prop->squeeze_range;
        for ( unsigned pp = 0; pp < nbPoints(); ++pp )
        {
#if ( DIM == 2 )
            unsigned ii = DIM * ( matIndex() + pp ) + 1;
            real p = posP(pp).YY;
            if ( p > r )
                meca.base(ii) -= f;
            else if ( p < -r )
                meca.base(ii) += f;
            else
                meca.mC(ii, ii) -= f / r;
#elif ( DIM == 3 )
            unsigned jj = DIM * ( matIndex() + pp );
            Vector p = posP(pp);
            if ( p.norm() < r )
            {
                meca.mC(jj+1, jj+1) -= f / r;
                meca.mC(jj+2, jj+2) -= f / r;
            }
            else
            {
                Vector n = p.normalized(f);
                meca.base(jj+1) -= n.YY;
                meca.base(jj+2) -= n.ZZ;
            }
#endif
        }
    }
#endif
    
    
#ifdef NEW_COLINEAR_FORCE
    /*
     add a length-dependent force acting parallel to the filament.
     A force proportional to the length of the segments is applied 
     on the model-points.
     */
    if ( prop->colinear_force )
    {
        Matrix::index_type inx = DIM * matIndex();
        real s = 0.5 * prop->colinear_force * segmentation();
        for ( unsigned pp = 0; pp < nbSegments(); ++pp )
        {
            Vector f = s * dirPoint(pp);
            meca.addForce(inx+DIM*pp    , f);
            meca.addForce(inx+DIM*pp+DIM, f);
        }
    }
#endif

#if ( 0 )
    /// add confinement forces with a second Space
    Space const* spc = prop->confine_outside_ptr;
    if ( spc )
    {
        for ( unsigned n = 0; n < nbPoints(); ++n )
        {
            Vector pos = posP(n);
            if ( spc->inside(pos) )
                spc->setInteraction(pos, PointExact(this, n), meca, prop->confine_stiffness);
        }
    }
#endif
    
    
    switch ( prop->confine )
    {
        case CONFINE_OFF:
            break;
        
        case CONFINE_INSIDE:
        {
            Space const* spc = prop->confine_space_ptr;
            
            for ( unsigned pp = 0; pp < nbPoints(); ++pp )
            {
                Vector pos = posP(pp);
                if ( spc->outside(pos) )
                    spc->setInteraction(pos, PointExact(this, pp), meca, prop->confine_stiffness);
            }
        } break;
        
        case CONFINE_OUTSIDE:
        {
            Space const* spc = prop->confine_space_ptr;
            
            for ( unsigned pp = 0; pp < nbPoints(); ++pp )
            {
                Vector pos = posP(pp);
                if ( spc->inside(pos) )
                    spc->setInteraction(pos, PointExact(this, pp), meca, prop->confine_stiffness);
            }
        } break;
        
        case CONFINE_ON:
        {
            Space const* spc = prop->confine_space_ptr;
            
            for ( unsigned pp = 0; pp < nbPoints(); ++pp )
                spc->setInteraction(posP(pp), PointExact(this, pp), meca, prop->confine_stiffness);
        } break;
        
        case CONFINE_MINUS_END:
        {
            Space const* spc = prop->confine_space_ptr;
            
            unsigned pp = 0;
            spc->setInteraction(posP(pp), PointExact(this, pp), meca, prop->confine_stiffness);
        } break;
        
        case CONFINE_PLUS_END:
        {
            Space const* spc = prop->confine_space_ptr;
            
            unsigned pp = lastPoint();
            spc->setInteraction(posP(pp), PointExact(this, pp), meca, prop->confine_stiffness);
        } break;

        case CONFINE_BOTH_ENDS:
        {
            Space const* spc = prop->confine_space_ptr;
            
            spc->setInteraction(posP(0), PointExact(this, 0), meca, prop->confine_stiffness);
            const unsigned L = lastPoint();
            spc->setInteraction(posP(L), PointExact(this, L), meca, prop->confine_stiffness);
        } break;

        case CONFINE_PLUS_OUT:
        {
            Space const* spc = prop->confine_space_ptr;
            
            unsigned pp = lastPoint();
            Vector pos = posP(pp);
            if ( spc->inside(pos) )
                spc->setInteraction(pos, PointExact(this, pp), meca, prop->confine_stiffness);
        } break;
        
        default:
            throw InvalidParameter("Invalid bead::confine");
    }
}


//------------------------------------------------------------------------------
#pragma mark -

void Fiber::addBinder(FiberBinder * fb)
{
    assert_true(fb->fiber() == this);
    frBinders.push_back(fb);
}


void Fiber::removeBinder(FiberBinder * fb)
{
    assert_true(fb->fiber() == this);
    frBinders.pop(fb);
}


FiberBinder * Fiber::firstBinder() const
{
    return static_cast<FiberBinder*>(frBinders.front());
}


void Fiber::detachBinders()
{
    //we iterate one step forward, because updating might lead to detachment:
    Node * hi = frBinders.front();
    while ( hi )
    {
        Hand * ha = static_cast<Hand*>(hi);
        hi = hi->next();
        ha->detach();
    }
    
}


int comp_abscissa(Node const* a, Node const* b)
{
    real aa = static_cast<FiberBinder const*>(a)->abscissa();
    real bb = static_cast<FiberBinder const*>(b)->abscissa();
    if ( aa < bb )
        return -1;
    if ( bb < aa )
        return 1;
    return 0;
}

/**
 This sorts the Binders in order of increasing abscissa
 Sorting is done (unefficiently) by insertion
 */
void Fiber::sortBinders()
{
    frBinders.sort(comp_abscissa);
}


unsigned Fiber::nbBinders(unsigned (*count)(const FiberBinder &)) const
{
    unsigned res = 0;
    
    Node * hi = frBinders.front();
    while ( hi ) {
        FiberBinder * ha = static_cast<FiberBinder*>(hi);
        hi = hi->next();
        res += count(*ha);
    }
    
    //printf("countBinders(%p) = %u\n", count, res);
    return res;
}


unsigned Fiber::nbBindersInRange(real a, real b, const FiberEnd ref) const
{
    if ( b > a )
        return 0;
    
    // Convert to absolute abscissa:
    a = abscissaFrom(a, ref);
    b = abscissaFrom(b, ref);
    if ( b < a )
    {
        real c = a;
        a = b;
        b = c;
    }
    
    unsigned res = 0;
    for ( Node * hi = frBinders.front(); hi; hi = hi->next() )
    {
        FiberBinder const* ha = static_cast<FiberBinder*>(hi);
        res += ( a <= ha->abscissa()  &&  ha->abscissa() <= b );
    }
    
    //printf("nbBinderssInRange(%8.2f, %8.2f)=%u\n", a, b, res);
    return res;
}


unsigned Fiber::nbBindersNearEnd(const real len, const FiberEnd from) const
{
    unsigned res = 0;
    
    for ( Node * hi = frBinders.front(); hi; hi = hi->next() )
    {
        FiberBinder const* ha = static_cast<FiberBinder*>(hi);
        res += ( ha->abscissaFrom(from) < len );
    }
    
    //printf("nbBindersNearEnd(%8.2f)=%i\n", len, res);
    return res;
}

//------------------------------------------------------------------------------
#pragma mark -

unsigned Fiber::dynamicState(FiberEnd which) const
{
    if ( which == PLUS_END )
        return dynamicStateP();
    if ( which == MINUS_END )
        return dynamicStateM();
    ABORT_NOW("invalid argument value");
    return 0;
}


void Fiber::setDynamicState(const FiberEnd which, const unsigned s)
{
    if ( which == PLUS_END )
        setDynamicStateP(s);
    else if ( which == MINUS_END )
        setDynamicStateM(s);
}


real Fiber::freshAssembly(const FiberEnd which) const
{
    if ( which == PLUS_END )
        return freshAssemblyP();
    if ( which == MINUS_END )
        return freshAssemblyM();
    ABORT_NOW("invalid argument value");
    return 0;
}


/**
 Update all known Binders.
 The loop could be unrolled, or parallelized
 */
void Fiber::updateBinders() const
{
    Node * hi = frBinders.front();
    while ( hi )
    {
        FiberBinder * ha = static_cast<FiberBinder*>(hi);
        hi = hi->next();
        ha->updateBinder();
    }
}


/**
 Assuming that the length has changed, or that the abscissa of the ends have changed,
 this updates the segmentation of the fiber if needed, the position of the binders,
 and the boundaries of the Lattice if present.
 */
void Fiber::update()
{
#if ( 0 )
    std::clog << reference() << " update [ "  << std::setw(9) << std::left << abscissaM();
    std::clog << " "  << std::setw(9) << std::left << abscissaP() << " ]\n";
#endif
    
    adjustSegmentation();
    
    /*
     Update all known Binders.
     Some of the Binders may be updated more than once in a time-step,
     but that is only a small performance penalty.

     We iterate one step ahead, because `checkFiberRange()` may lead to detachment.
     The loop could be unrolled, or parallelized
     */
    //ADDED BY MANU
#if (0)
    if (prop->sweep_digits)
        sweep();
#endif
    Node * hi = frBinders.front();
    while ( hi )
    {
        Hand * ha = static_cast<Hand*>(hi);
        hi = hi->next();
        assert_true(ha->fiber()==this);
        ha->updateBinder();
        ha->checkFiberRange();
    }
    
    if ( frLattice )
    {
        frLattice->setRange(abscissaM(), abscissaP());
        
        if ( prop->field_ptr )
        {
            real sum;
            // release Lattice substance located outside the valid abscissa range
            frLattice->collectM(sum, frLattice->index(abscissaM()));
            prop->field_ptr->cell(posEndM()) += sum;
            //std::clog << " Fiber::MINUS_END releases " << sumM << std::endl;
            
            frLattice->collectP(sum, frLattice->index(abscissaP()));
            prop->field_ptr->cell(posEndP()) += sum;
            //std::clog << " Fiber::PLUS_END releases " << sumP << std::endl;
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 */
void Fiber::setLattice(real conc) const
{
    const real uni = frLattice->unit();
    const int  inf = frLattice->index(abscissaM());
    const int  sup = frLattice->index(abscissaP());
    assert_true( inf <= sup );
    
    FiberLattice::site_type * site = frLattice->data();
    
    if ( inf == sup )
    {
        //the Fiber is entirely covered by one site!
        assert_true( frLattice->abscissa(inf+1) >= abscissaP() );
        site[inf] = conc * ( abscissaP() - abscissaM() );
    }
    else
    {
        // the terminal site may be truncated
        site[inf] = conc * ( frLattice->abscissa(inf+1) - abscissaM() );
        
        for ( int h = inf+1; h < sup; ++h )
        site[h] = conc * uni;
        
        // the terminal site may be truncated
        site[sup] = conc * ( abscissaP() - frLattice->abscissa(sup) );
    }
}


/**
 Release a fraction 'frac' of the Lattice substance into the Field.
 The subtance in each Lattice site is released in a cell
 corresponding to a random position within this site.
 The factor `frac` must be between 0 and 1.
 */
void Fiber::releaseLattice(Field * fld, real frac) const
{
    assert_true( 0 <= frac && frac <= 1 );
    const real uni = frLattice->unit();
    const int  inf = frLattice->index(abscissaM());
    const int  sup = frLattice->index(abscissaP());
    assert_true( inf <= sup );
    
    FiberLattice::site_type * site = frLattice->data();
    
    if ( inf == sup )
    {
        //the Fiber is entirely covered by one site!
        assert_true( frLattice->abscissa(inf+1) >= abscissaP() );
        fld->cell(posM(RNG.preal()*length())) += frac * site[inf];
        site[inf] -= frac * site[inf];
    }
    else
    {
        // the terminal site may be truncated
        real s = frLattice->abscissa(inf+1) - abscissaM();
        fld->cell(pos(uni*(inf+1)-RNG.preal()*s)) += frac * site[inf];
        site[inf] -= frac * site[inf];
        
        for ( int h = inf+1; h < sup; ++h )
        {
            real a = frac * site[h];
            site[h] -= a;
            // we select a random position along each site to release:
            fld->cell(pos(uni*(RNG.preal()+h))) += a;
        }
        
        // the terminal site may be truncated
        s = abscissaP() - frLattice->abscissa(sup);
        fld->cell(pos(uni*sup+RNG.preal()*s)) += frac * site[sup];
        site[sup] -= frac * site[sup];
    }
}


/**
 Release all Lattice substance into the Field.
 The subtance in each Lattice site is released in a cell
 corresponding to a random position within this site.
 */
void Fiber::releaseLattice(Field * fld) const
{
    FiberLattice::site_type * site = frLattice->data();
    const real uni = frLattice->unit();
    
    //@todo Handle the terminal site differently since they are truncated
    const int sup = frLattice->sup();
    for ( int h = frLattice->inf(); h < sup; ++h )
    {
        fld->cell(pos(uni*(RNG.preal()+h))) += site[h];
        site[h] = 0;
    }
    
    //std::clog << "remaining lattice " << frLattice->sum() << std::endl;
}


void Fiber::resetLattice() const
{
    if ( frLattice )
    {
        frLattice->clear();
        
        ABORT_NOW("unfinished");
        
#if ( 0 )
        /// unfinished code
        Node * hi = frBinders.first();
        while ( hi )
        {
            FiberBinder * ha = static_cast<FiberBinder*>(hi);
            hi = hi->next();
            if ( ha->useLattice() )
                frLattice.inc(ha->site());
        }
#endif
    }
}


void Fiber::readLattice(Inputter& in)
{
    if ( frLattice )
    {
        frLattice->read(in);
        real diff = fabs( frLattice->unit() - prop->lattice_unit );
        if ( diff > REAL_EPSILON )
        {
            //throw InvalidIO("mismatch between lattice_unit and Lattice::unit()");
            //std::cerr << frLattice->unit() << "  " << prop->lattice_unit << "\n";
            frLattice->unit(prop->lattice_unit);
            //resetLattice();
        }
    }
    else
    {
        FiberLattice dummy;
        dummy.read(in);
    }
}



//------------------------------------------------------------------------------
#pragma mark -

/**
 fiber:glue=1 creates a Single if the tip of the Fiber goes outside the Space.
 The Single's hand is managed to always be attached at the tip of the Fiber.
 The Single detaches if the Fiber tip is pulled inside.
 This generates mostly a pushing force from the cortex
 */
void Fiber::setGlue1(Single* glue, const FiberEnd end, const Space * spc)
{
    assert_true(spc);
    if ( spc->inside(posEnd(end)) )
    {
        //detach immediately if the tip is inside the Space
        if ( glue->attached() )
            glue->detach();
    }
    else
    {
        if ( glue->attached() )
        {
            //always keep tracking the tip:
            glue->moveToEnd(end);
        }
        else {
            //reposition the Single base:
            glue->setPosition(spc->project(posEnd(end)));
            //attach to the MT-tip:
            glue->attachToEnd(this, end);
        }
    }
}


/**
 fiber:glue=2
 The Single's hand is managed to always be attached at the tip of the Fiber.
 The Single's hand detaches only spontaneously.
 This creates both pulling and pushing force from the cortex
 */
void Fiber::setGlue2(Single* glue, const FiberEnd end, const Space * spc)
{
    assert_true(spc);
    if ( glue->attached() )
    {
        //keep tracking the tip of the fiber while attached
        glue->moveToEnd(end);
    }
    else
    {
        // Attach a new grafted if MT-tip is outside and growing:
        if ( isGrowing(end) && spc->outside(posEnd(end)) )
        {
            //reposition the Single base:
            glue->setPosition(spc->project(posEnd(end)));
            //attach to the MT-tip:
            glue->attachToEnd(this, end);
        }
    }
}


/**
 fiber:glue=3 creates a Single at the position where the Fiber crosses the Space's edge.
 This makes an anchor point exactly at the cortex.
 The Single's Hand behaves and detaches normally.
 */
void Fiber::setGlue3(Single* glue, const Space * spc)
{    
    assert_true(spc);
    /*
     If the glue is not already attached, we first check if the fiber intersects
     the edge of the Space:
     */
    if ( ! glue->attached() )
    {
        bool in = spc->inside(posEndM());
        
        if ( in == spc->inside(posEndP()) )
            return;
        
        // find a model point that is on the other side of the Space edge:
        for ( int pp = 1; pp < nbPoints(); ++pp )
        {
            if ( spc->inside(posP(pp)) != in )
            {
                // the abscissa is interpolated using the distances of P1 and P2 to the edge
                real d1 = spc->distanceToEdge(posP(pp-1));
                real d2 = spc->distanceToEdge(posP(pp));
                if ( d1 + d2 > REAL_EPSILON )
                {
                    /* we find the abscissa corresponding to the intersection,
                     assuming that the edge is locally straight */
                    FiberBinder fs(this, abscissaPoint(pp-1+d1/(d2+d1)));
                    glue->attach(fs);
                    glue->setPosition(fs.pos());
                    break;
                }
            }
        }
    }
}


/**
 Search for a glue in the list of bound HandSingle
 this is useful when a simulation is restarted from file
 */
void Fiber::makeGlue(Single*& glue)
{
    for ( Single * gh = prop->glue_set->firstA(); gh; gh=gh->next() )
    {
        if ( gh->hand()->fiber() == this  &&  gh->mark() == identity() )
        {
            glue = gh;
            //std::clog << "found Fiber:glue for " << reference() << "\n";
            return;
        }
    }
    
    // create the Single if needed
    if ( glue == 0 )
    {
        glue = prop->glue_prop->newSingle();
        glue->mark(identity());
        prop->glue_set->add(glue);
    }
}


void Fiber::setGlue(Single*& glue, const FiberEnd end, const Space * space, int mode)
{
    assert_true(space);
    
    if ( glue == 0 )
        makeGlue(glue);
    
    // creates Single when MT interact with the cortex:
    switch( mode )
    {
        case 1:  setGlue1(glue, end, space);  break;
        case 2:  setGlue2(glue, end, space);  break;
        case 3:  setGlue3(glue, space);       break;
        default: throw InvalidParameter("invalid value of fiber:glue");
    }
    
#if ( 1 )
    // we keep the Single linked in the simulation only if it is attached:
    if ( glue->attached() )
    {
        if ( !glue->linked() )
            prop->glue_set->link(glue);
    }
    else
    {
        if ( glue->linked() )
            prop->glue_set->unlink(glue);
    }
#endif
}


//------------------------------------------------------------------------------
#pragma mark -

void Fiber::write(Outputter& out) const
{
    Filament::write(out);
    
    if ( frLattice )
    {
        out.put_char('\n');
        writeReference(out, TAG_LATTICE);
        frLattice->write(out);
    }    
}


void Fiber::read(Inputter & in, Simul& sim, Tag tag)
{
#ifdef BACKWARD_COMPATIBILITY
        
    if ( in.formatID() == 33 )
        mark(in.readUInt32());
    
    if ( tag == 'm' )
    {
        if ( in.formatID()==31 )
        {
            unsigned p = in.readUInt16();
            prop = sim.findProperty<FiberProp*>("fiber", p);
        }
        //tag = TAG;
    }
    
    if ( in.formatID() < 31 )
    {
        setDynamicStateM(in.readUInt8());
        setDynamicStateP(in.readUInt8());
    }
    
#endif
    
    if ( tag == TAG )
    {
        Filament::read(in, sim, tag);
        
        if ( length() + REAL_EPSILON < prop->min_length )
        {
            std::cerr << "Warning: fiber length < fiber:min_length";
            std::cerr << " ( " << length() << " < " << prop->min_length << " )\n";
        }

        if ( frLattice )
            frLattice->setRange(abscissaM(), abscissaP());

        frGlue = 0;
    }
    else if ( tag == TAG_LATTICE )
    {
        readLattice(in);
    }
#ifdef BACKWARD_COMPATIBILITY
    else if ( tag == 'm' )
    {
        Filament::read(in, sim, tag);
    }
#endif
    else
    {
        std::clog << "unknown Fiber TAG `" << (char)tag << "'\n";
    }
}

//ADDED by manu
/**

fiber f: pointer to the fiber where you want to attach the free hand (new)
pos_old: position of the hand that was attached
pos_new: position of the new site where the hand will try to attach

 */
bool Fiber::causesEntanglement( Vector const & pos_old, Vector const & pos_new, const Fiber * f) const
{
    
    for ( Node * hi = frBinders.front(); hi; hi = hi->next() )
    {
        // Is this static cast correct to be used here?
        // hand of the "old fiber"-> fiber where the hand was attached
        Hand const* ha_old = static_cast<Hand*>(hi);
        // hand of the "new fiber"-> fiber where the hand is trying to attach
        Hand const* ha_new = ha_old->otherHand();

        // This also tests wether the fbinder is a couple or a single by testing if ha_new is true
        // TODO: This might be a bit sloppy
        if (ha_new && ha_new->fiber()==f && ha_new->prop->lat_val && ha_old->prop->lat_val)
        {
            if ((ha_new->pos()-pos_new)*(ha_old->pos()-pos_old)<0)
                return true;
        }
    }
    return false;
}


void Fiber::sweep()
{
    if (lattice())
    {
        // Sort by abscissa and start by the closest to the plus end one
        frBinders.sort(comp_abscissa);
        Node * hi = frBinders.front();
        Hand * ha, * prv;
        
        bool sweep_any = false;
        while ( hi )
        {
            
            ha = static_cast<Hand*>(hi);
            if (ha->prop->name()=="ase_hand")
            {
                if (ha->needs_sweep())
                {
                    // The first couple of ase_1 we encounter,
                    // We try to rescue. If there is a rescue, no need to sweep
                    if (!sweep_any) {
                        if (ha->first_sweep()) {
                            break;
                        }
                    }
                    sweep_any = true;
                    prv = ha;
                }
                else if (sweep_any)
                {
                    prv->sweep();
                    break;
                }
                // If the first couple that is ase_walker is not passed
                // the abscissa, there is nothing to do
                else
                    break;
            }
            hi = hi->next();
        }
        frBinders.mix(RNG);
    }
    
}
#ifdef MULTI_LATTICE
unsigned int Fiber::get_lattice_val(Fiber * f)
{
    for (unsigned int i = 0; i<4; i++) {
        if (multi_lattice[i] == f)
        {
            return i+1;
        }
    }
    return 0;
}
unsigned int Fiber::available_lattice()
{

    for (int i = 0; i<4; i++)
    {
        // Return the first available lattice
        
        if (!multi_lattice[i])
        {
            return i +1;}
        
    }

    return 0;
}
#endif
real Fiber::measure_cap(unsigned int lat_val)
{
    // This is junk
    FiberLattice * lat = lattice();
    //Since we allow two holes
    int count = -3;
    if (lat) {
        //Count from the plus end
        int lat_site = abscissaP()/lat->unit();
        int two_holes = 0;
        while (two_holes<2) {
            if (lat->vacant(lat_site, lat_val))
                two_holes+=1;
            else
                two_holes=0;
            lat_site-=1;
            count++;
            
        }
        return count*lat->unit();
    }
    else
        return 0;
    
}

FiberBinder * Fiber::findHand(unsigned int lat_val, int site) const
{
    if (!frBinders.size()) {
        return nullptr;
    }
    for ( Node * hi = frBinders.front(); hi; hi = hi->next() )
    {
        Hand * hand = static_cast<Hand*>(hi);
        if (hand->lattice_site()==site && hand->prop->lat_val==lat_val) {
            return hand;
        }
    }
    return nullptr;
}

FiberBinder * Fiber::findHand(unsigned int lat_val, int site, int fiber_id)  const
{
    if (!frBinders.size()) {
        return nullptr;
    }
    for ( Node * hi = frBinders.front(); hi; hi = hi->next() )
    {
        Hand * hand = static_cast<Hand*>(hi);
#ifdef MULTI_LATTICE
        if (hand->lattice_site()==site && hand->prop->lat_val==lat_val && hand->get_fiberlat_id() ==fiber_id) {
#else
        if (hand->lattice_site()==site && hand->prop->lat_val==lat_val) {
#endif
            return hand;
        }
    }
    return nullptr;
}

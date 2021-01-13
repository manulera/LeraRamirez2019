// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "hand.h"
#include "hand_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "fiber_prop.h"
#include "simul.h"
#include "sim.h"
extern Random RNG;

//------------------------------------------------------------------------------

Hand::Hand(HandProp const* p, HandMonitor* m) : haMonitor(m), prop(p)
{
    // initialize in unattached state:
    nextAttach = RNG.exponential();
    nextDetach = 0;
#ifdef TRAP_SINGLES
#if (TRAP_SINGLES==1)
    trapped_haMon = 0;
#endif
#endif
}

Hand::~Hand()
{
    if ( attached() )
        detach();
    prop = 0;
}

//------------------------------------------------------------------------------
#pragma mark -

Hand * Hand::otherHand() const
{
    return haMonitor->otherHand(this);
}

/**
Checks that all the conditions required for attachment are met
 */
int Hand::attachmentAllowed(FiberBinder & fb) const
{
	assert_true( fb.attached() );
	
#ifdef MULTI_LATTICE
    // A bit sloppy, but prevents hands that do not use a lattice to not do the binding either
    

    if (!prop->lat_val)
    {
        if (otherHand() && otherHand()->attached())
        {
            Fiber * attached_fb = otherHand()->fiber();
            Fiber * new_fb = fb.fiber();
            if (!(attached_fb->get_lattice_val(new_fb) && new_fb->get_lattice_val(attached_fb))) {
                return false;
            }
        }
    }
#endif
    
    
	/*
	 Check that the two binding keys match:
	 Allow binding if the BITWISE-AND of the two keys is true
	 */
	if ( ! ( prop->binding_key & fb.fiber()->prop->binding_key ) )
        return false;
	
	// check end-on binding:
	if ( fb.abscissaFromM() < 0 )
	{
		if ( prop->bind_also_ends )
            fb.moveToEndM();
		else
			return false;
	}
	else if ( fb.abscissaFromP() < 0 )
	{
		if ( prop->bind_also_ends )
			fb.moveToEndP();
		else
            return false;
	}
    
    FiberEnd end = NO_END;

    switch ( prop->bind_only_end )
    {
        case NO_END:
            break;
        case MINUS_END:
            if ( fb.abscissaFromM() > prop->bind_end_range )
                return false;       // too far from fiber end
            end = MINUS_END;
            break;
        case PLUS_END:
            if ( fb.abscissaFromP() > prop->bind_end_range )
                return false;       // too far from fiber end
            end = PLUS_END;
            break;
        case BOTH_ENDS:
        {
            FiberEnd e = nearestEnd();
        
            if ( fb.abscissaFrom(e) > prop->bind_end_range )
                return false;
        
            // also check the other fiber end:
            FiberEnd o = ( e == PLUS_END ) ? MINUS_END : PLUS_END;
        
            // give equal chance to all ends within range:
            if ( fb.abscissaFrom(o) < prop->bind_end_range  )
                end = RNG.choice(MINUS_END, PLUS_END);
            else
                end = e;
        }
        default:
            throw Exception("Illegal value of hand:bind_only_end");
    }

    // check occupancy near the end:
	if ( end != NO_END && prop->bind_only_free_end )
	{
		if ( 0 < fb.fiber()->nbBindersNearEnd(prop->bind_end_range, end) )
            return false;
	}

	// also check the Monitor's permissions:
	return haMonitor->allowAttachment(fb);
}



void Hand::attach(FiberBinder const& fb)
{
    assert_true( unattached() );
    assert_true( fb.attached() );

    nextDetach = RNG.exponential();
    FiberBinder::locate(fb);
    haMonitor->afterAttachment(this);
}


void Hand::attachTo(Fiber * f, const real ab)
{
    assert_true(f);
    attach(FiberBinder(f, ab));
}


void Hand::attachTo(Fiber * f, const real ab, const FiberEnd ref)
{
    assert_true(f);
    attach(FiberBinder(f, f->abscissaFrom(ab, ref)));
}


void Hand::attachToEnd(Fiber * f, const FiberEnd end)
{
    assert_true(f);
    attach(FiberBinder(f, f->abscissaEnd(end)));
}



//------------------------------------------------------------------------------
#pragma mark -


/**
 This attempts to move to a different location on the same Fiber,
 but it can also lead to detachment, in particular if the location is 
 not between the fiber Ends
 */
void Hand::moveTo(const real abs)
{
    assert_true( attached() );
    
    // check if movement would lead beyond the Minus End:
    if ( abs < fbFiber->abscissaM() )
    {
        moveToEndM();
        if ( !prop->hold_growing_end )
            detach();
        return;
    }
    
    // check if movement would lead beyond the Plus End:
    if ( abs > fbFiber->abscissaP() )
    {
        moveToEndP();
        if ( !prop->hold_growing_end )
            detach();
        return;
    }
    
    fbAbs = abs;
    updateBinder();
}


void Hand::checkFiberRange()
{
    assert_true( attached() );
    if ( fbAbs < fbFiber->abscissaM() )
        handleDisassemblyM();
    
    else if ( fbAbs > fbFiber->abscissaP() )
        handleDisassemblyP();
}


void Hand::handleDisassemblyM()
{
    if ( prop->hold_shrinking_end )
        moveToEndM();
    else
        detach();
}

void Hand::handleDisassemblyP()
{
    if ( prop->hold_shrinking_end)
        moveToEndP();
    else
        detach();
}


//------------------------------------------------------------------------------
#pragma mark -

void Hand::detach()
{
    assert_true( attached() );
    haMonitor->beforeDetachment(this);
    FiberBinder::delocate();
    nextAttach = RNG.exponential();
}


/**
 Test for spontaneous detachment using Gillespie approach.
 
 @return true if the test has passed, and detach() was called.
 
 see @ref Stochastic
 */
bool Hand::testDetachment()
{
    assert_true( nextDetach >= 0 );

    nextDetach -= prop->unbinding_rate_dt;
    
    if ( nextDetach <= 0 )
    {
        detach();
        return true;
    }
    
    return false;
}


/**
 Test for spontaneous detachment using Gillespie approach.
 
 @return true if the test has passed, and detach() was called.
 
 see @ref Stochastic
 */
bool Hand::testKramersDetachment(const real force)
{
    assert_true( nextDetach >= 0 );
    
    /*
     Attention: the exponential term can easily become numerically "infinite",
     but `prop->unbinding_force_inv` is set to zero if 'unbinding_rate==0`
     */
    nextDetach -= prop->unbinding_rate_dt * exp(force*prop->unbinding_force_inv);
    if ( nextDetach <= 0 )
    {
        detach();
        return true;
    }
    return false;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Test for attachment to nearby Fibers, using the Gillespie time nextAttach
 */
void Hand::stepUnattached(const FiberGrid& grid, Vector const & pos)
{
    assert_true( unattached() );
    assert_true( nextAttach >= 0 );

#ifndef TRICKY_HAND_ATTACHMENT
    grid.tryToAttach(pos, *this);
#else
    // we test attachement with a rate that is 4x higher than the binding rate:
    nextAttach -= prop->binding_rate_dt_4;
    while ( nextAttach < 0 )
    {
        /*
         For each test, the probability to attach is 1/4, since we provide
         here a test integer '1<<30' which is 4 times smaller than '1<<32',
         which pass always. 
         The configuration of the hands & filaments will not change,
         between sucessive calls tryToAttach() here, but this should happen rarely.
         */
        grid.tryToAttach(pos, *this);
            
        if ( attached() )
            break;
        
        nextAttach += RNG.exponential();
    }
#endif
}



/**
 Test for spontaneous detachment at rate HandProp::unbinding_rate, 
 */
void Hand::stepUnloaded()
{
    assert_true( attached() );
    
    testDetachment();
}


/**
 Test for force-induced detachment following Kramers' law,
 vith basal rate HandProp::unbinding_rate, 
 and characteristic force HandProp::unbinding_force
 */
void Hand::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    if ( prop->unbinding_force_inv > 0 )
        testKramersDetachment(force.norm());
    else
        testDetachment();
}



//------------------------------------------------------------------------------
#pragma mark -


void Hand::resetTimers()
{
    // initialize the Gillespie counters:
    if ( attached() )
        nextDetach = RNG.exponential();
    else
        nextAttach = RNG.exponential();
}


void Hand::write(Outputter& out) const
{
    /*
     it is not necessary to write the property index here,
     since it is set only when Hand is created in class Single or Couple.
     */
    FiberBinder::write(out);
}


void Hand::read(Inputter & in, Simul& sim)
{
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 32 )
        prop = sim.findProperty<HandProp*>("hand",in.readUInt16());
#endif
    
    FiberBinder::read(in, sim);
    resetTimers();
}
#ifdef TRAP_SINGLES
#if (TRAP_SINGLES == 1)
Hand * Hand::trapped_hand() const {return trapped_haMon->trappedHand();};

bool Hand::trapped(){return trapped_haMon!=0;};

void Hand::untrap(){trapped_haMon = 0;};

// Set trapped_haMon to h->haMonitor and called afterTrapping(), this one is for singles to be moved to the SingleSet::trappedList
void Hand::trap(Hand * h){trapped_haMon = h->haMonitor; haMonitor->afterTrapping();};

void Hand::stepTrapped(Vector force){trapped_haMon->stepTrapped(force);};

Vector Hand::trappedOtherHandPos(){return trapped_haMon->trappedHandPos();};
#endif

#if (TRAP_SINGLES==2)

HandMonitor * Hand::trappedHaMon(){return haMonitor->trapped_haMon;};
bool Hand::trapped(){return trappedHaMon()!=0;};
void Hand::stepUnattachedTrappedAA()
{
}


#endif

#endif















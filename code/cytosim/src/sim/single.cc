// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "sim.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "glossary.h"
#include "iowrapper.h"
#include "single.h"
#include "simul.h"
#include "space.h"
#include "modulo.h"
#include "meca.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------
Single::Single(SingleProp const* p, Vector const& w)
: sPos(w), sHand(0), prop(p)
{
    if ( p == 0 )
        throw Exception("Null Single::prop");
    
    assert_true(prop->hand_prop);
    sHand = prop->hand_prop->newHand(this);
    assert_true(sHand);
    
#ifdef TRAP_SINGLES
    trapped_haMon = 0;
#endif
}


Single::~Single()
{
    if ( linked() )
        objset()->remove(this);

    if ( sHand  &&  sHand->attached() )
        sHand->detach();

    if ( sHand )
    {
        delete(sHand);
        sHand = 0;
    }
    
    prop = 0;
}

//------------------------------------------------------------------------------
#pragma mark -

void Single::afterAttachment(Hand const*)
{
    assert_true( attached() );
    // If it binds from the trapper, it should remain in the trapped list.
    if ( linked() && !trapped() )
    {
        SingleSet * set = static_cast<SingleSet*>(objset());
        set->relinkA(this);
    }
}


void Single::beforeDetachment(Hand const* h)
{
    assert_true( h == sHand );

#if ( 0 )
    // relocate Single to position where it was attached
    sPos = h->posHand();
#else
    /*
     Set position near the attachment point, but offset in the perpendicular
     direction at a random distance within the range of attachment of the Hand.
     
     This is necessary to achieve detailled balance, which in particular implies
     that rounds of binding/unbinding should not get the Singles closer to
     the Filaments to which they bind.
     */
    sPos = h->posHand() + h->dirFiber().randOrthoB(h->prop->binding_range);
#endif
    
    // If it still trapped, it should remain in the list of trapped singles even if it detached

    if ( linked() && !trapped() )
    {
        SingleSet * set = static_cast<SingleSet*>(objset());
#ifdef TRAP_SINGLES
#if (TRAP_SINGLES==1)
        if (sHand->trapped())
            untrap(true);
#endif
#endif
        set->relinkD(this);
    }
}


real Single::interactionLength() const
{
    return prop->length;
}


//------------------------------------------------------------------------------
#pragma mark -


Vector Single::position() const
{
    if ( sHand->attached() )
        return sHand->pos();
    return sPos;
}

void Single::foldPosition(const Modulo * s)
{
    s->fold(sPos);
}

void Single::randomizePosition()
{
    sPos = prop->confine_space_ptr->randomPlace();
}


void Single::stepF(const FiberGrid& grid)
{
    assert_false( sHand->attached() );

    // diffusion:
    sPos.addRand(prop->diffusion_dt);
    
    // confinement
    if ( prop->confine == CONFINE_INSIDE )
    {
        if ( !prop->confine_space_ptr->inside(sPos) )
            prop->confine_space_ptr->bounce(sPos);
        if ( modulo )
            modulo->fold(sPos);
    }
    else if ( prop->confine == CONFINE_ON )
    {
        Vector pos = sPos;
        prop->confine_space_ptr->project(pos, sPos);
    }
    
    sHand->stepUnattached(grid, sPos);

}


void Single::stepA()
{
    assert_true( sHand->attached() );
    assert_true( !hasForce() );

    sHand->stepUnloaded();
}

/**
 Add confinement force to the bound fiber
 */
void Single::setInteractions(Meca & meca) const
{
    assert_true( sHand->attached() );
    
#ifdef NEW_DANGEROUS_CONFINEMENTS
    if ( prop->confine )
    {
        const Space* spc = prop->confine_space_ptr;
        spc->setInteraction(sHand->interpolation(), meca, prop->stiffness, prop->confine);
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark -

void Single::write(Outputter& out) const
{
    sHand->write(out);
    out.writeFloatVector(sPos, DIM);
#ifdef TRAP_SINGLES
    if (!trapped())
        Object::writeNullReference(out);
    else
        trapped_haMon->writeReference(out);
#endif
}


void Single::read(Inputter & in, Simul& sim, Tag tag)
{
    try
    {
        /*
         Because the SingleSet has 2 lists where Single are stored depending
         on their bound/unbound state, we need to unlink and relink a Single here,
         since the state stored on file could be different from the current state.
         */

        ObjectSet * set = objset();
        
        // There should not be anything trapped in the beginning, so this should be fine
        if ( set )
            set->unlink(this);
        
        
        sHand->read(in, sim);
        in.readFloatVector(sPos, DIM);
        // See Single::write
#ifdef TRAP_SINGLES
        Tag tag = 0;
        Object * w = sim.readReference(in, tag);
        if (w)
            trapped_haMon = static_cast<Couple*>(w);
        else
            trapped_haMon = 0;
#endif
        if ( set )
            set->link(this);
    }
    catch( Exception & e ) {
        e << ", in Single::read()";
        throw;
    }
}

#ifdef TRAP_SINGLES

#if (TRAP_SINGLES ==1)
void Single::afterTrapping()
{
    if ( linked() )
    {
        // This should prevent this from doing anything (neither in the aList or dList)
        SingleSet * set = static_cast<SingleSet*>(objset());
        set->trap_single(this);
    }
    
}
void Single::untrap(bool first_call)
{

    // If it is the first call, we have to also untrap the partner
    if (first_call) {
        sHand->trapped_haMon->untrap(false);
    }
    sHand->untrap();
    SingleSet * set = static_cast<SingleSet*>(objset());
    set->untrap_single(this);
}


#else
void Single::untrap()
{
    trapped_haMon = 0;
    SingleSet * set = static_cast<SingleSet*>(objset());
    set->untrap_single(this);
}
void Single::trap(HandMonitor * haMon)
{
    trapped_haMon = haMon;
    // This should prevent this from doing anything (neither in the aList or dList)
    SingleSet * set = static_cast<SingleSet*>(objset());
    set->trap_single(this);
}


#endif

void Single::stepTrappedF(const FiberGrid& grid, Vector const & pos)
{
    sHand->stepUnattached(grid, pos);
}

void Single::stepTrappedA(Vector force)
{
    sHand->stepLoaded(force);
}

void Single::stepTrappedF_AA()
{
    sHand->stepUnattachedTrappedAA();
}

void Single::writeReference(Outputter& out) const
{
    Object::writeReference(out,tag());
}

#if (TRAP_SINGLES==1)
Vector Single::trappedHandPos() const
{
    assert_true(trapped());
    return sHand->pos();
}
Hand * Single::trappedOtherHand() const
{
    return trappedHand()->trapped_haMon->trappedHand();
}
Vector Single::trappedOtherHandPos() const
{
    return trappedOtherHand()->pos();
}
Hand * Single::trappedHand() const
{
    return sHand;
}

#endif

#endif




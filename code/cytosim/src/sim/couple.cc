// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "sim.h"
#include "couple.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "hand_prop.h"
#include "meca.h"
#include "simul.h"
#include "space.h"
#include "modulo.h"
#include "aster.h"
#include "aster_prop.h"

extern Random RNG;
extern Modulo const* modulo;

//------------------------------------------------------------------------------

Couple::Couple(CoupleProp const* p, Vector const& w)
: prop(p), cPos(w), cHand1(0), cHand2(0)
{
    if ( p == 0 )
        throw Exception("Null Couple::prop");

    cHand1 = prop->hand1_prop->newHand(this);
    cHand2 = prop->hand2_prop->newHand(this);
    
    assert_true( cHand1 );
    assert_true( cHand2 );
    
}


Couple::~Couple()
{
    if ( linked() )
        objset()->remove(this);

    if ( cHand1  &&  cHand1->attached() )
        cHand1->detach();
    
    if ( cHand2  &&  cHand2->attached() )
        cHand2->detach();
    
    if ( cHand1 )
    {
        delete(cHand1);
        cHand1 = 0;
    }
    if ( cHand2 )
    {
        delete(cHand2);
        cHand2 = 0;
    }
    
    prop = 0;
}


//------------------------------------------------------------------------------

void Couple::setProperty(CoupleProp * p)
{
    if ( p == 0 )
        throw Exception("Null Couple::prop");
    prop = p;
    
    if ( cHand1 )
        delete(cHand1);
    cHand1 = prop->hand1_prop->newHand(this);
    
    if ( cHand2 )
        delete(cHand2);
    cHand2 = prop->hand2_prop->newHand(this);
}


//------------------------------------------------------------------------------
#pragma mark -

real Couple::stiffness() const
{
    return prop->stiffness;
}


void Couple::setInteractions(Meca & meca) const
{
    assert_true( cHand1->attached() && cHand2->attached() );
    
    meca.interLink(cHand1->interpolation(), cHand2->interpolation(), prop->stiffness);
    
#ifdef NEW_DANGEROUS_CONFINEMENTS
    if ( prop->confine )
    {
        const Space* spc = prop->confine_space_ptr;
        spc->setInteraction(cHand1->interpolation(), meca, prop->stiffness, prop->confine);
        spc->setInteraction(cHand2->interpolation(), meca, prop->stiffness, prop->confine);
    }
#endif
}


void Couple::setInteractionsAF(Meca & meca) const
{
    assert_true( cHand1->attached() && !cHand2->attached() );
    
#ifdef NEW_DANGEROUS_CONFINEMENTS
    if ( prop->confine )
    {
        const Space* spc = prop->confine_space_ptr;
        spc->setInteraction(cHand1->interpolation(), meca, prop->stiffness, prop->confine);
    }
#endif
}


void Couple::setInteractionsFA(Meca & meca) const
{
    assert_true( !cHand1->attached() && cHand2->attached() );
    
#ifdef NEW_DANGEROUS_CONFINEMENTS
    if ( prop->confine )
    {
        const Space* spc = prop->confine_space_ptr;
        spc->setInteraction(cHand2->interpolation(), meca, prop->stiffness, prop->confine);
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Simulates:
 - diffusive motion
 - attachment
 .
 */
void Couple::stepFF(const FiberGrid& grid)
{
    // diffusive motion:
    cPos.addRand(prop->diffusion_dt);
    // confinement:
    if ( prop->confine == CONFINE_INSIDE )
    {
        /**
         @todo dirichlet boundary conditions
         Set concentration of molecules at edges of Space by letting molecules
         out, and put some back at a constant rate
         */
        if ( !prop->confine_space_ptr->inside(cPos) )
            prop->confine_space_ptr->bounce(cPos);
        if ( modulo )
            modulo->fold(cPos);
    }
    else if ( prop->confine == CONFINE_ON )
    {
        const Space* spc = prop->confine_space_ptr;
        Vector pos = cPos;
        spc->project(pos, cPos);
    }

    /*
     To attachment a Couple, we flip a coin to give equal chance to each Hand,
     as if they were sharing the two half of a spherical cap.
     Note that this divides by two the effective binding rate of the Hands.
     */
    if ( RNG.flip() )
    {
        cHand1->stepUnattached(grid, cPos);
    }
    else
    {
        if ( !prop->trans_activated )
            cHand2->stepUnattached(grid, cPos);
    }
}


/**
 Simulates:
 - attachment of cHand2
 - attached activity of cHand1
 .
 */
void Couple::stepAF(const FiberGrid& grid)
{
    //we use cHand1->pos() first, because stepUnloaded() may detach cHand1
    cHand2->stepUnattached(grid, cHand1->pos());
    cHand1->stepUnloaded();
}


/**
 Simulates:
 - attachment of cHand1
 - attached activity of cHand2
 .
 */
void Couple::stepFA(const FiberGrid& grid)
{
    //we use cHand2->pos() first, because stepUnloaded() may detach cHand2
    cHand1->stepUnattached(grid, cHand2->pos());
    cHand2->stepUnloaded();
}


/**
 Simulates:
 - attached activity of cHand1
 - attached activity of cHand2
 .
 */
void Couple::stepAA()
{
    Vector f = force();
    cHand1->stepLoaded( f);
    cHand2->stepLoaded(-f);
}


real Couple::new_gilles_t() const
{
    if ( attached1() )
    {
        if ( attached2() )
        {
            Vector f = force();
            return -log(RNG.preal())/(cHand1->propensLoaded(f) + cHand2->propensLoaded(-f));
        }
        else
            return -log(RNG.preal()) / cHand1->propensUnloaded();
    }
    else
    {
        if ( attached2() )
            return -log(RNG.preal()) / cHand2->propensUnloaded();
        else
            return 0;
    }
}


int Couple::gillestep()
{
    if ( attached1() )
    {
        if ( attached2() )
        {
            Vector f = force();
            real p_1 = cHand1->propensLoaded( f);
            real p_2 = cHand2->propensLoaded(-f);
            if ( RNG.preal()*(p_1+p_2) < p_1 )
                cHand1->gillestepLoaded(f);
            else
                cHand2->gillestepLoaded(-f);
            return 0;
        }
        else
        {
            return cHand1->gillestepUnloaded();
        }
    }
    else
    {
        if ( attached2() )
            return cHand2->gillestepUnloaded();
        else
            return 0;
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 @return:
 - True if attachment is possible
 - False if attachment is forbiden
 .
 
 If ( couple:stiff == true ), the two Hands of the Couple will refuse to be attached
 to the same segment, or to two neighboring segments on the same fiber.
 
 We cannot calculate the force of such 'degenerate' links, and they are undesired in
 most cases.
 
 */

bool Couple::allowAttachment(const FiberBinder & fb)
{
    Hand const* that = attachedHand();
    
    if ( that == 0 )
        return true;
    
    if ( prop->stiff )
    {
        if ( that->fiber() == fb.fiber()
            && fabs(fb.abscissa()-that->abscissa()) < 2*fb.fiber()->segmentation() )
        return false;
    
#if ( 0 )
        // Outdated code, 2013
        /*
         Test here if binding would create a link inside an aster, near the center:
         i.e. a link between two Fibers from the same aster, near the center
         of this aster. Such links would be improductive, and would trap the Couples.
         */
        const Organizer * org = fb.fiber()->organizer();
        
        if ( org  &&  org->tag() == Aster::TAG )
        {
            Fiber * fib = that->fiber();
            if ( fib->organizer() == org )
            {
                real rad = static_cast<const Aster*>(org)->prop->radius[0];
                if (  fb.abscissa() < rad  &&  that->abscissa() < rad )
                    return false;
            }
        }
#endif
    }
    
    switch( prop->specificity )
    {
        case CoupleProp::BIND_ALWAYS:
            return true;
            
        case CoupleProp::BIND_PARALLEL:
            if ( fb.dirFiber() * that->dirFiber() < 0.5 )
                return false;
            break;
            
        case CoupleProp::BIND_NOT_PARALLEL:
            if ( fb.dirFiber() * that->dirFiber() > 0.5 )
                return false;
            break;
  
        case CoupleProp::BIND_ANTIPARALLEL:
            if ( fb.dirFiber() * that->dirFiber() > -0.5 )
                return false;
            break;
            
        case CoupleProp::BIND_NOT_ANTIPARALLEL:
            if ( fb.dirFiber() * that->dirFiber() < -0.5 )
                return false;
            break;
            
        case CoupleProp::BIND_ORTHOGONAL:
            if ( fabs( fb.dirFiber() * that->dirFiber() ) > 0.866025 )
                return false;
            break;
            
        default:
            throw InvalidParameter("unknown couple:specificity");
    }

    //attachment is allowed by default:
    return true;
}


void Couple::afterAttachment(Hand const* h)
{
    if ( linked() )
    {
        CoupleSet * set = static_cast<CoupleSet*>(objset());
        if ( h == cHand1 )
            set->relinkA1(this);
        else
            set->relinkA2(this);
    }
    
#if ( 0 )
    if ( cHand1->attached() && cHand2->attached() )
        printf("%i %.9f %.9f\n", 3+prop->fast_diffusion, cHand1->nextAttach, cHand2->nextAttach);
    else if ( cHand1->unattached() )
        printf("%i %.9f  %.9f\n", 1+prop->fast_diffusion, cHand1->nextAttach, cHand2->nextAttach);
    else if ( cHand2->unattached() )
        printf("%i %.9f  %.9f\n", 1+prop->fast_diffusion, cHand2->nextAttach, cHand1->nextAttach);
#endif
}


void Couple::beforeDetachment(Hand const* h)
{
    assert_true(h->attached());
#if ( 0 )
    // relocate Couple to position where it was attached
    cPos = h->posHand();
#else
    /*
     Set position near the attachment point, but offset in the perpendicular
     direction at a random distance within the range of attachment of the Hand
 
     This is necessary to achieve detailled balance, which in particular implies
     that rounds of binding/unbinding should not get the Couples closer to
     the Filaments to which they bind.
     */
    cPos = h->posHand() + h->dirFiber().randOrthoB(h->prop->binding_range);
#endif
    
    if ( linked() )
    {
        CoupleSet * set = static_cast<CoupleSet*>(objset());
        if ( h == cHand1 )
            set->relinkD1(this);
        else
            set->relinkD2(this);
    }
    
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 The position is:
 - cPos if the Couple is free,
 - the position of the attached Hand if only one is attached
 - the average position of the two hands if they are both attached
.
 */
Vector Couple::position() const
{
    if ( cHand2->attached() )
    {
        if ( cHand1->attached() )
            return 0.5 * ( cHand2->pos() + cHand1->pos() );
        return cHand2->pos();
    }
    if ( cHand1->attached() )
    {
        return cHand1->pos();
    }
    return cPos;
}


void Couple::foldPosition(const Modulo * s)
{
    modulo->fold(cPos);
}

void Couple::randomizePosition()
{
    cPos = prop->confine_space_ptr->randomPlace();
}

//------------------------------------------------------------------------------
#pragma mark -

Vector Couple::force() const
{
    Vector d = cHand2->pos() - cHand1->pos();
    
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


Hand * Couple::attachedHand() const
{
    if ( attached1() )
        return cHand1;
    else if ( attached2() )
        return cHand2;
    else
        return 0;
}


Hand* Couple::otherHand(Hand const* h) const
{
    if ( h == cHand1 )
        return cHand2;
    else
        return cHand1;
}


Vector Couple::otherPosition(Hand const* h) const
{
    if ( h == cHand1 )
    {
        if ( cHand2->attached() )
            return cHand2->pos();
        throw Exception("otherPosition() called for unattached Hand2");
    }
    else
    {
        if ( cHand1->attached() )
            return cHand1->pos();
        throw Exception("otherPosition() called for unattached Hand1");
    }
}


Vector Couple::otherDirection(Hand const* h) const
{
    if ( h == cHand1 )
    {
        if ( cHand2->attached() )
            return cHand2->dirFiber();
        else
            return Vector::randU();
    }
    else
    {
        if ( cHand1->attached() )
            return cHand1->dirFiber();
        else
            return Vector::randU();
    }
}


real Couple::interactionLength() const
{
    return prop->length;
}

//------------------------------------------------------------------------------
#pragma mark -


void Couple::write(Outputter& out) const
{
    cHand1->write(out);
    cHand2->write(out);
    if ( !cHand1->attached() && !cHand2->attached() )
        out.writeFloatVector(cPos, DIM);
}


void Couple::read(Inputter & in, Simul& sim, Tag tag)
{
    try {
        
        /* 
         Because the CoupleSet has 4 lists where Couple are stored depending
         on their bound/unbound state, we need to unlink and relink a Couple here,
         since the state stored on file could be different from the current state.
         */
        ObjectSet * set = objset();
        if ( set )
            set->unlink(this);

        cHand1->read(in, sim);
        cHand2->read(in, sim);
        
        if ( cHand1->attached() || cHand2->attached() )
            cPos = position();
        else
            in.readFloatVector(cPos, DIM);

        if ( set )
            set->link(this);
    }
    catch( Exception & e ) {
        e << ", in Couple::read()";
        throw;
    }
}


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "smath.h"
#include "assert_macro.h"
#include "treadmilling_fiber.h"
#include "treadmilling_fiber_prop.h"
#include "fiber_locus.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"
#include "space.h"

extern Random RNG;

//------------------------------------------------------------------------------

TreadmillingFiber::TreadmillingFiber(TreadmillingFiberProp const* p) : Fiber(p), prop(p)
{
    mStateM  = STATE_WHITE;
    mGrowthM = 0;
    
    mStateP  = STATE_WHITE;
    mGrowthP = 0;
}


TreadmillingFiber::~TreadmillingFiber()
{
    prop = 0;
}


//------------------------------------------------------------------------------
#pragma mark -

unsigned TreadmillingFiber::dynamicStateM() const
{
    return mStateM;
}


void TreadmillingFiber::setDynamicStateM(unsigned s)
{
    if ( s!=STATE_WHITE && s!=STATE_GREEN && s!=STATE_RED )
        throw InvalidParameter("Invalid AssemblyState for TreadmillingFiber MINUS_END");
    
    mStateM = s;
}


real TreadmillingFiber::freshAssemblyM() const
{
    return mGrowthM;
}




unsigned TreadmillingFiber::dynamicStateP() const
{
    return mStateP;
}


void TreadmillingFiber::setDynamicStateP(unsigned s)
{
    if ( s!=STATE_WHITE && s!=STATE_GREEN && s!=STATE_RED )
        throw InvalidParameter("Invalid AssemblyState for TreadmillingFiber PLUS_END");
    
    mStateP = s;
}


real TreadmillingFiber::freshAssemblyP() const
{
    return mGrowthP;
}


//------------------------------------------------------------------------------

void TreadmillingFiber::step()
{    
    if ( mStateP == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceP = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        mGrowthP = prop->growing_speed_dt[0] * prop->free_polymer;
        
        assert_true(mGrowthP>=0);
        // antagonistic force (< 0) decreases assembly rate exponentially
        if ( forceP < 0  &&  prop->growing_force[0] < INFINITY )
            mGrowthP *= exp(forceP/prop->growing_force[0]);
    }
    else if ( mStateP == STATE_RED )
    {
        mGrowthP = prop->shrinking_speed_dt[0];
    }
    else
    {
        mGrowthP = 0;
    }
    
    
    if ( mStateM == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceM = projectedForceEndM();
        
        // growth is reduced if free monomers are scarce:
        mGrowthM = prop->growing_speed_dt[1] * prop->free_polymer;

        assert_true(mGrowthM>=0);
        // antagonistic force (< 0) decreases assembly rate exponentially
        if ( forceM < 0  &&  prop->growing_force[1] < INFINITY  )
            mGrowthM *= exp(forceM/prop->growing_force[1]);
    }
    else if ( mStateM == STATE_RED )
    {
        mGrowthM = prop->shrinking_speed_dt[1];
    }
    else
    {
        mGrowthM = 0;
    }
    
    
    const real len = length();
    if ( len + mGrowthP + mGrowthM < prop->min_length )
    {
        // the fiber is too short, we delete it:
        delete(this);
        return;
    }
    else if ( len + mGrowthP + mGrowthM < prop->max_length )
    {
        growM(mGrowthM);
        growP(mGrowthP);
    }
    else if ( len < prop->max_length )
    {
        // the remaining possible growth is distributed to the two ends:
        real c = ( prop->max_length - len ) / ( mGrowthM + mGrowthP );
        growM(c*mGrowthM);
        growP(c*mGrowthP);
    }

    Fiber::step();
}

                  
//------------------------------------------------------------------------------
#pragma mark -


void TreadmillingFiber::write(Outputter& out) const
{
    Fiber::write(out);

    /// write variables describing the dynamic state of the ends:
    out.put_char('\n');
    writeReference(out, TAG_DYNAMIC);
    out.writeUInt16(mStateM);
    out.writeUInt16(mStateP);
}


void TreadmillingFiber::read(Inputter & in, Simul& sim, Tag tag)
{
#ifdef BACKWARD_COMPATIBILITY
    if ( tag == TAG_DYNAMIC || ( tag == TAG && in.formatID() < 44 ) )
#else
    if ( tag == TAG_DYNAMIC )
#endif
    {
        mStateM = in.readUInt16();
        mStateP = in.readUInt16();
    }
#ifdef BACKWARD_COMPATIBILITY
    if ( tag != TAG_DYNAMIC || in.formatID() < 44 )
#else
    else
#endif
        Fiber::read(in, sim, tag);
}


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "smath.h"
#include "assert_macro.h"
#include "growing_fiber.h"
#include "growing_fiber_prop.h"
#include "fiber_locus.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"
#include "space.h"

extern Random RNG;

//------------------------------------------------------------------------------

GrowingFiber::GrowingFiber(GrowingFiberProp const* p) : Fiber(p), prop(p)
{
    mStateM = STATE_GREEN;
    mStateP = STATE_GREEN;
    mGrowthM = 0;
    mGrowthP = 0;
}


GrowingFiber::~GrowingFiber()
{
    prop = 0;
}


//------------------------------------------------------------------------------
#pragma mark -

unsigned GrowingFiber::dynamicStateM() const
{
    return mStateM;
}


unsigned GrowingFiber::dynamicStateP() const
{
    return mStateP;
}


void GrowingFiber::setDynamicStateM(unsigned s)
{
    if ( s == STATE_WHITE || s == STATE_GREEN )
        mStateM = s;
    else
        throw InvalidParameter("invalid AssemblyState for a GrowingFiber");
}


void GrowingFiber::setDynamicStateP(unsigned s)
{
    if ( s == STATE_WHITE || s == STATE_GREEN )
        mStateP = s;
    else
        throw InvalidParameter("invalid AssemblyState for a GrowingFiber");
}


real GrowingFiber::freshAssemblyM() const
{
    return mGrowthM;
}


real GrowingFiber::freshAssemblyP() const
{
    return mGrowthP;
}

//------------------------------------------------------------------------------

void GrowingFiber::step()
{

    // PLUS_END
    if ( prop->shrink_outside[0] && prop->confine_space_ptr->outside(posEndP()) )
    {
        mGrowthP = -prop->growing_speed_dt[0];
    }
    
    // ADDED BY MANU: New confining method added to keep a constant overlap length
    else if (prop->growth_space_ptr) {
        if ( prop->growth_space_ptr->outside(posEndP()) )
        {
            mGrowthP = -prop->growing_speed_dt[0];
        }
        else if ( prop->growth_space_ptr->inside(posEndP()))
        {
            mGrowthP = prop->growing_speed_dt[0];
        }
    }
    
    else if ( mStateP == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceP = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        mGrowthP = prop->growing_speed_dt[0] * prop->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if ( forceP < 0  &&  prop->growing_force[0] < INFINITY )
            mGrowthP *= exp(forceP/prop->growing_force[0]);
        
        mGrowthP += prop->growing_off_speed_dt[0];
    }
    else
    {
        mGrowthP = 0;
    }
    
    // MINUS_END
    if ( prop->shrink_outside[1] && prop->confine_space_ptr->outside(posEndM()) )
    {
        mGrowthM = -prop->growing_speed_dt[1];
    }
    else if ( mStateM == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceM = projectedForceEndM();
        
        // growth is reduced if free monomers are scarce:
        mGrowthM = prop->growing_speed_dt[1] * prop->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if ( forceM < 0  &&  prop->growing_force[1] < INFINITY )
            mGrowthM *= exp(forceM/prop->growing_force[1]);

        mGrowthM += prop->growing_off_speed_dt[1];
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
    else
    {
        mGrowthM = 0;
        mGrowthP = 0;
    }

    Fiber::step();
}



//------------------------------------------------------------------------------
#pragma mark -


void GrowingFiber::write(Outputter& out) const
{
    Fiber::write(out);

    /// write variables describing the dynamic state of the ends:
    out.put_char('\n');
    writeReference(out, TAG_DYNAMIC);
    out.writeFloat(mGrowthM);
    out.writeFloat(mGrowthP);
}


void GrowingFiber::read(Inputter & in, Simul& sim, Tag tag)
{
#ifdef BACKWARD_COMPATIBILITY
    if ( tag == TAG_DYNAMIC || ( tag == TAG && in.formatID() < 44 ) )
#else
    if ( tag == TAG_DYNAMIC )
#endif
    {
        mGrowthM = in.readFloat();
#ifdef BACKWARD_COMPATIBILITY
        if ( in.formatID() > 45 )
#endif
        mGrowthP = in.readFloat();
    }
#ifdef BACKWARD_COMPATIBILITY
    if ( tag != TAG_DYNAMIC || in.formatID() < 44 )
#else
    else
#endif
    {
#ifdef BACKWARD_COMPATIBILITY
        const real len = length();
#endif
        
        Fiber::read(in, sim, tag);
        
#ifdef BACKWARD_COMPATIBILITY
        if ( tag == TAG && in.formatID() < 46 )
        {
            // adjust growing variable
            mGrowthP = length() - len;
            mGrowthM = 0;
        }
#endif
    }
}


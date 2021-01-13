// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "smath.h"
#include "assert_macro.h"
#include "classic_fiber.h"
#include "classic_fiber_prop.h"
#include "fiber_locus.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"
#include "space.h"

extern Random RNG;

//------------------------------------------------------------------------------

ClassicFiber::ClassicFiber(ClassicFiberProp const* p) : Fiber(p), prop(p)
{
    mStateM  = STATE_WHITE;
    mStateP  = STATE_WHITE;
    mGrowthM = 0;
    mGrowthP = 0;
}


ClassicFiber::~ClassicFiber()
{
    prop = 0;
}


bool valid_state(int s)
{
    return  s==STATE_WHITE || s==STATE_GREEN || s==STATE_RED;
}


void ClassicFiber::setDynamicStateM(unsigned s)
{
    if ( !valid_state(s) )
        throw InvalidParameter("Invalid AssemblyState for ClassicFiber MINUS_END");
    
    if ( s != mStateM )
    {
        mStateM = (AssemblyState)s;
    }
}


void ClassicFiber::setDynamicStateP(unsigned s)
{
    if ( !valid_state(s) )
        throw InvalidParameter("Invalid AssemblyState for ClassicFiber PLUS_END");
    
    if ( s != mStateP )
    {
        mStateP = (AssemblyState)s;
    }
}


//------------------------------------------------------------------------------
#pragma mark -

/** 
 The catastrophe rate depends on the growth rate of the corresponding tip,
 which is itself reduced by antagonistic force. 
 The correspondance is : 1/rate = a + b * growthSpeed.
 For no force on the growing tip: rate = catastrophe_rate[0]*time_step
 For very large forces          : rate = catastrophe_rate_stalled[0]*time_step
 cf. `Dynamic instability of MTs is regulated by force`
 M.Janson, M. de Dood, M. Dogterom. JCB 2003, Figure 2 C.
 */
void ClassicFiber::step()
{
    const real len = length();

    if ( mStateM == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real force = projectedForceEndM();
        
        // growth is reduced if free monomers are scarce:
        real spd = prop->growing_speed_dt[1] * prop->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if ( force < 0  &&  prop->growing_force[1] < INFINITY )
            mGrowthM = spd * exp(force/prop->growing_force[1]) + prop->growing_off_speed_dt[1];
        else
            mGrowthM = spd + prop->growing_off_speed_dt[1];
        
        
        // catastrophe may be constant, or it may depend on the growth rate
        real cata;
        if ( prop->catastrophe_coef[1] > 0 )
            cata = prop->catastrophe_rate_stalled_dt[1] / ( 1.0 + prop->catastrophe_coef[1] * mGrowthM );
        else
            cata = prop->catastrophe_rate_dt[1];

        if ( RNG.test(cata) )
            mStateM = STATE_RED;
    }
    else if ( mStateM == STATE_RED )
    {
        mGrowthM = prop->shrinking_speed_dt[1];
        
        if ( RNG.test(prop->rescue_prob[0]) )
            mStateM = STATE_GREEN;
    }
    
    
    if ( mStateP == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real force = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        real spd = prop->growing_speed_dt[0] * prop->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if ( force < 0  &&  prop->growing_force[0] < INFINITY )
            mGrowthP = spd * exp(force/prop->growing_force[0]) + prop->growing_off_speed_dt[0];
        else
            mGrowthP = spd + prop->growing_off_speed_dt[0];
        
        
        // catastrophe may be constant, or it may depend on the growth rate
        real cata;
        if ( prop->catastrophe_coef[0] > 0 )
            cata = prop->catastrophe_rate_stalled_dt[0] / ( 1.0 + prop->catastrophe_coef[0] * mGrowthP );
        else
            cata = prop->catastrophe_rate_dt[0];
        
#ifdef NEW_LENGTH_DEPENDENT_CATASTROPHE
        /*
         Ad-hoc length dependence, used to simulate S. pombe with catastrophe_length=5
         Foethke et al. MSB 5:241 - 2009
         */
        if ( prop->catastrophe_length > 0 )
        {
            PRINT_ONCE("Using ad-hoc length-dependent catastrophe rate\n");
            cata *= length() / prop->catastrophe_length;
        }
#endif
        
#ifdef NEW_CATASTROPHE_OUTSIDE
        /*
         Catastrophe will be triggered immediately if the PLUS_END is outside
         */
        if ( prop->catastrophe_outside && prop->confine_space_ptr->outside(posEndP()) )
        {
            mStateP = STATE_RED;
        }
#endif
        
        if ( RNG.test(cata) )
            mStateP = STATE_RED;
    }
    else if ( mStateP == STATE_RED )
    {
        mGrowthP = prop->shrinking_speed_dt[0];
        
        if ( len + mGrowthP < prop->min_length )
        {
            if ( RNG.test(prop->rebirth_prob[0]) )
                mStateP = STATE_GREEN;
        }
        else
        {
            if ( RNG.test(prop->rescue_prob[0]) )
                mStateP = STATE_GREEN;
        }
    }
    
    
    if ( len + mGrowthP + mGrowthM < prop->min_length )
    {
        // the fiber is too short, we may delete it:
        if ( !prop->persistent )
        {
            delete(this);
            return;
        }
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


void ClassicFiber::write(Outputter& out) const
{
    Fiber::write(out);
    
    /// write variables describing the dynamic state of the ends:
    out.put_char('\n');
    writeReference(out, TAG_DYNAMIC);
    out.writeUInt16(mStateM);
    out.writeUInt16(mStateP);
}


void ClassicFiber::read(Inputter & in, Simul& sim, Tag tag)
{
#ifdef BACKWARD_COMPATIBILITY
    if ( tag == TAG_DYNAMIC || ( tag == TAG && in.formatID() < 44 ) )
#else
    if ( tag == TAG_DYNAMIC )
#endif
    {
        unsigned m, p;
#ifdef BACKWARD_COMPATIBILITY
        if ( in.formatID() < 42 )
            p = in.readUInt8();
        else
#endif
        {
            m = in.readUInt16();
            p = in.readUInt16();
        }

#ifdef BACKWARD_COMPATIBILITY
        if ( in.formatID() < 46 )
            setDynamicStateP(m);
        else
#endif
        {
            setDynamicStateM(m);
            setDynamicStateP(p);
        }
    }
#ifdef BACKWARD_COMPATIBILITY
    if ( tag != TAG_DYNAMIC || in.formatID() < 44 )
#else
    else
#endif
        Fiber::read(in, sim, tag);
}


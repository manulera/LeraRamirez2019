// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dynein.h"
#include "dynein_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "lattice.h"
#include "simul.h"

extern Random RNG;



Dynein::Dynein(DyneinProp const* p, HandMonitor* h)
: Digit(p,h), prop(p)
{
    ABORT_NOW("unfinished class");
    nextStep = RNG.exponential();
}


/**
 \todo simulate occurence of backward steps
 */
void Dynein::stepUnloaded()
{
    assert_true( attached() );
    
    nextStep   -= prop->stepping_rate_dt;
    
    while ( nextStep <= 0 )
    {
        assert_true( attached() );
        if ( stepM() == 2 )
        {
            nextStep = RNG.exponential();
            //immediately detach at the end of the Fiber:
            detach();
            return;
        }
        nextStep += RNG.exponential();
    }
    
    testDetachment();
}



/**
 Currently, antagonistic force only reduced the rate of forward stepping.
 However, force is also known to increase the rate of backward steps.
 \todo simulate occurence of backward steps
 */
void Dynein::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    // calculate displacement, dependent on the load along the desired direction of displacement
    real rate_step = prop->stepping_rate_dt + ( force * dirFiber() ) * prop->var_rate_dt;

    nextStep   -= rate_step;
    
    while ( nextStep <= 0 )
    {
        assert_true( attached() );
        if ( stepM() == 2 )
        {
            nextStep = RNG.exponential();
            //immediately detach at the end of the Fiber:
            detach();
            return;
        }
        nextStep += RNG.exponential();
    }

    testKramersDetachment(force.norm());
}




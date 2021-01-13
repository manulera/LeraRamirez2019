#include "dim.h"
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"
#include "common.h"
#include "simul_prop.h"
#include "swapper_prop.h"
#include "swapper.h"


void SwapperProp::clear()
{
    CoupleProp::clear();
    swap_rate = 0;
    swap_rate_dt = 0;
}

void SwapperProp::read(Glossary& glos)
{
    CoupleProp::read(glos);
    glos.set(swap_rate,           "swap_rate");
    
}

void SwapperProp::complete(Simul const* sim)
{
    CoupleProp::complete(sim);
    swap_rate_dt = swap_rate*sim->prop->time_step/sim->prop->handmonitor_pace;
}


Couple * SwapperProp::newCouple(Glossary*) const
{
    //std::clog << "CrosslinkProp::newCouple" << std::endl;
    if ( length > 0 )
        throw InvalidParameter("Swapper is not implemented for couples of length different from zero, if you want to make that, create a Swapper_long class\n");
    else
        return new Swapper(this);
}


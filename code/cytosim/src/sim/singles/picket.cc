// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "picket.h"
#include "simul.h"
#include "meca.h"
#include "modulo.h"


extern Modulo const* modulo;


Picket::Picket(SingleProp const* p, Vector const& w)
: Single(p, w)
{
#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter("single:diffusion cannot be > 0 if activity=fixed");
#endif
}


Picket::~Picket()
{
    //std::clog<<"~Picket("<<this<<")"<<std::endl;
}


void Picket::beforeDetachment(Hand const*)
{
    assert_true( attached() );

    if ( linked() )
    {
        SingleSet * set = static_cast<SingleSet*>(objset());
        set->relinkD(this);
    }
}


void Picket::stepF(const FiberGrid& grid)
{
    assert_false( sHand->attached() );

    sHand->stepUnattached(grid, sPos);
}


void Picket::stepA()
{
    assert_true( sHand->attached() );
    
    sHand->stepLoaded(force());
}


/**
 This calculates the force corresponding to addPointClamp()
 */
Vector Picket::force() const
{
    assert_true( sHand->attached() );
    Vector d = sPos - posHand();
    
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


void Picket::setInteractions(Meca & meca) const
{
    assert_true( prop->length == 0 );
    meca.addPointClamp(sHand->interpolation(), sPos, prop->stiffness);
    //meca.addLineClamp(sHand->interpolation(), sPos, sHand->dirFiber(), prop->stiffness);
}



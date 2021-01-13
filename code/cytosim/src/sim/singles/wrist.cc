// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "wrist.h"
#include "simul.h"
#include "meca.h"
#include "modulo.h"


extern Modulo const* modulo;


Wrist::Wrist(SingleProp const* sp, Mecable const* mec, const unsigned pti)
: Single(sp), mBase(mec), mOrder(1)
{
    if ( mec && pti > mec->nbPoints() )
        throw InvalidParameter("Could not anchor Single (invalid point index)");

    mPoint = pti;
    mCoef[0] = 1.0;
    
    for ( int i = 1; i < 4; ++i )
        mCoef[i] = 0;

#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter("single:diffusion cannot be > 0 if activity=anchored");
#endif
}


Wrist::Wrist(SingleProp const* sp, Mecable const* mec, unsigned ref, Vector pos)
: Single(sp), mBase(mec), mOrder(DIM+1)
{
    assert_true(mec);

    if ( ref+DIM > mec->nbPoints() )
        throw InvalidParameter("Could not anchor Single (invalid point index)");
    
    mPoint = ref;
    mCoef[0] = 1.0;
    for ( int i = 0; i < DIM; ++i )
    {
        mCoef[i+1] = pos[i];
        mCoef[0]  -= pos[i];
    }
    
    for ( int i = DIM+1; i < 4; ++i )
        mCoef[i] = 0;
    
#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter("single:diffusion cannot be > 0 if activity=anchored");
#endif
}


Wrist::~Wrist()
{
}


Vector Wrist::posFoot() const
{
    if ( mOrder == 1 )
        return mBase->posPoint(mPoint);
    else
    {
        Vector res = mCoef[0] * mBase->posPoint(mPoint);
        for ( int i = 1; i < mOrder; ++i )
            res += mCoef[i] * mBase->posPoint(mPoint+i);
        return res;
    }
}

Vector Wrist::force() const
{
    assert_true( sHand->attached() );
    Vector d = posFoot() - sHand->pos();
    
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


void Wrist::stepF(const FiberGrid& grid)
{
    assert_false( sHand->attached() );

    sHand->stepUnattached(grid, posFoot());
}


void Wrist::stepA()
{
    assert_true( sHand->attached() );
    
    sHand->stepLoaded(force());
}


void Wrist::setInteractions(Meca & meca) const
{
    if ( mOrder == 1 )
        meca.interLink(sHand->interpolation(), PointExact(mBase, mPoint), prop->stiffness);
    else
    {
        unsigned off = mBase->matIndex() + mPoint;
        unsigned pts[] = { 0, 1, 2, 3 };
#if ( DIM == 1 )
        meca.interLink2(sHand->interpolation(), off, pts, mCoef, prop->stiffness);
#elif ( DIM == 2 )
        meca.interLink3(sHand->interpolation(), off, pts, mCoef, prop->stiffness);
#elif ( DIM == 3 )
        meca.interLink4(sHand->interpolation(), off, pts, mCoef, prop->stiffness);
#endif
    }
}


void Wrist::write(Outputter& out) const
{
    sHand->write(out);
    out.writeSoftSpace();
    mBase->writeReference(out);
    out.writeUInt16(mPoint);
    for ( int i = 1; i < 4; ++i )
        out.writeFloat(mCoef[i]);
}


void Wrist::read(Inputter & in, Simul& sim, Tag tag)
{
    try
    {
        sHand->read(in, sim);
#ifdef BACKWARD_COMPATIBILITY
        if ( in.formatID() < 47 )
        {
            PointExact base;
            base.read(in, sim);
            mBase = base.mecable();
            mPoint = base.point();
            mCoef[0] = 1.0;
            for ( int i = 1; i < 4; ++i )
                mCoef[i] = 0;
            return;
        }
#endif
        Tag tag = 0;
        mBase = Simul::toMecable(sim.readReference(in, tag));
        mPoint = in.readUInt16();
        for ( int i = 1; i < 4; ++i )
            mCoef[i] = in.readFloat();
        mCoef[0] = 1.0 - mCoef[1] - mCoef[2] - mCoef[3];
    }
    catch( Exception & e ) {
        e << ", in Single::read()";
        throw;
    }
}


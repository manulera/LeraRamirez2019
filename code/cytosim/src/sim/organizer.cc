// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "organizer.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"



Organizer::~Organizer()
{
    //MSG(31, "destroying Organizer %p\n", this);
}


void Organizer::grasp(Mecable * m)
{
    Buddy::connect(m);
    mObjects.push_back(m);
}


void Organizer::grasp(Mecable * m, unsigned ix)
{
    if ( ix >= mObjects.size() )
        mObjects.resize(ix+1, 0);

    if ( m != mObjects[ix] )
    {
        Buddy::disconnect(mObjects[ix]);
        Buddy::connect(m);
    }
    
    mObjects[ix] = m;
}


void Organizer::goodbye(Buddy * b)
{
    if ( b )
    {
        //std::clog << this << " organizer loses " << b << std::endl;
        MecableList::iterator oi = std::find(mObjects.begin(), mObjects.end(), b);
        if ( oi != mObjects.end() )
            *oi = 0;
    }
}


void Organizer::addOrganized(Simul & simul)
{
    for ( MecableList::iterator oi = mObjects.begin(); oi < mObjects.end(); ++oi )
    {
        if ( ! (*oi)->linked() )
        {
            std::clog << " Registering " << (*oi)->reference() << "\n";
            simul.add(*oi);
        }
        //(*oi)->mark(identity());
    }
}

//------------------------------------------------------------------------------
/**
 \return The centroid from all the object positions
 */
Vector Organizer::position() const
{
    Vector res(0,0,0);
    for ( MecableList::const_iterator oi = mObjects.begin(); oi < mObjects.end(); ++oi )
        res += (*oi)->position();
    return res / mObjects.size();
}


Vector Organizer::positionP(unsigned ix) const
{
    Vector res(0,0,0);
    for ( MecableList::const_iterator oi = mObjects.begin(); oi < mObjects.end(); ++oi )
        res += (*oi)->posPoint(ix);
    return res / mObjects.size();
}


void Organizer::translate(Vector const& T)
{
    for ( MecableList::iterator oi = mObjects.begin(); oi < mObjects.end(); ++oi )
    {
        Mecable * mv = *oi;
        if ( mv && mv->mobile() )
        {
            mv->translate(T);
            mv->flag(0);
        }
    }
}


void Organizer::rotate(Rotation const& T)
{
    for ( MecableList::iterator oi = mObjects.begin(); oi < mObjects.end(); ++oi )
    {
        Mecable * mv = *oi;
        if ( mv && mv->mobile() )
        {
            mv->rotate(T);
            mv->flag(0);
        }
    }
}


real Organizer::dragCoefficient() const
{
    real res = 0;
    for ( MecableList::const_iterator oi = mObjects.begin(); oi < mObjects.end(); ++oi )
        res += (*oi)->dragCoefficient();
    return res;
}


//------------------------------------------------------------------------------

void Organizer::write(Outputter& out) const
{
    out.writeUInt16(mObjects.size());
    out.writeSoftNewline();
    for ( MecableList::const_iterator oi = mObjects.begin(); oi < mObjects.end(); ++oi )
    {
        out.writeSoftSpace();
        if ( *oi )
            (*oi)->writeReference(out);
        else
            Object::writeNullReference(out);
    }
}


void Organizer::read(Inputter & in, Simul& sim, Tag tag)
{
    try
    {
        unsigned nb = in.readUInt16();
        
        //std::clog << " Organizer::read with " << nb << " objects" << std::endl;
        for ( unsigned mi = 0; mi < nb; ++mi )
        {
            Tag tag = 0;
            Object * w = sim.readReference(in, tag);
            if ( w )
                grasp(Simul::toMecable(w), mi);
       }
    }
    catch( Exception & e )
    {
        e << ", in Organizer::read()";
        throw;
    }
}

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_set.h"
#include "space_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"
#include "space.h"
#include "modulo.h"

//---------------------------- GLOBAL VARIABLE ---------------------------------

/**
 This is a global variable that is initialized in Simul
 */
Modulo const* modulo = 0;

//------------------------------------------------------------------------------

Property * SpaceSet::newProperty(const std::string& kd,const std::string& nm, Glossary&) const
{
    if ( kd == "space" )
        return new SpaceProp(nm);
    return 0;
}

//------------------------------------------------------------------------------
void SpaceSet::step()
{
    for ( Space * sp = first(); sp; sp=sp->next() )
        sp->step();
}


//------------------------------------------------------------------------------
void SpaceSet::erase()
{
    ObjectSet::erase();
    
    // simul has lost its current Space:
    simul.changeSpace(0);
}

/**
 This will change the Simul current Space if it was not set
*/
void SpaceSet::add(Object * obj)
{
    assert_true(obj->tag() == Space::TAG);
    //std::clog << "SpaceSet::add " << obj << std::endl;
    ObjectSet::add(obj);
    
    Space const* spc = simul.space();
    if ( spc == 0 || obj->identity() < spc->identity() )
        simul.changeSpace(static_cast<Space*>(obj));
}

/**
 If the Simulation current Space is deleted,
 the 'oldest' remaining Space is chosen to replace it.
 */
void SpaceSet::remove(Object * obj)
{
    //std::clog << "SpaceSet::remove " << obj << std::endl;
    ObjectSet::remove(obj);

    if ( obj == simul.space() )
    {
        /*
         if the current space was deleted, use the oldest Space available
         */
        Space * spc = first();
        
        for ( Space * s=spc; s; s=s->next() )
            if ( s->identity() < spc->identity() )
                spc = s;
        
        simul.changeSpace(spc);
    }
}

//------------------------------------------------------------------------------
Object * SpaceSet::newObjectT(const Tag tag, unsigned idx)
{
    Space * obj = 0;
    if ( tag == Space::TAG )
    {
        SpaceProp * p = simul.findProperty<SpaceProp*>("space", idx);
        if ( p == 0 )
           throw InvalidIO("no space class defined with id "+sMath::repr(idx));
        Glossary tmp;
        obj = p->newSpace(tmp);
    }
    return obj;
}

/**
 The dimensions of a Space can be specified when it is created
 @code
 new space cell
 {
    dimension = 3 4
 }
 @endcode
 */
ObjectList SpaceSet::newObjects(const std::string& name, Glossary& opt)
{
    Property * p = simul.properties.find_or_die("space", name);
    SpaceProp * sp = static_cast<SpaceProp*>(p);
    Space * obj = sp->newSpace(opt);

    ObjectList res;
    res.push_back(obj);
    return res;
}

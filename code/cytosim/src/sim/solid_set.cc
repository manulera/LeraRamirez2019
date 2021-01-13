// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "solid_set.h"
#include "solid_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"
#include "wrist.h"


#if ( 0 )
void SolidSet::step()
{
    for ( Solid * o = first(); o; o=o->next() )
        o->step();
}
#endif


//------------------------------------------------------------------------------

Property* SolidSet::newProperty(const std::string& kd, const std::string& nm, Glossary&) const
{
    if ( kd == "solid" )
        return new SolidProp(kd, nm);
    return 0;
}


Object * SolidSet::newObjectT(const Tag tag, unsigned idx)
{
    if ( tag == Solid::TAG )
    {
        SolidProp * p = simul.findProperty<SolidProp*>("solid", idx);
#ifdef BACKWARD_COMPATIBILITY
        // prior to 04.2016, "bead" and "solid" were used interchangeably
        if ( p == 0 )
             p = simul.findProperty<SolidProp*>("bead", idx);
#endif
        if ( p )
            return new Solid(p);
        return 0;
    }
    return 0;
}


/**
@ref Solid::build
 */
ObjectList SolidSet::newObjects(const std::string& name, Glossary& opt)
{
    Property * p = simul.properties.find_or_die("solid", name);
    Solid * obj = new Solid(static_cast<SolidProp*>(p));
    
    ObjectList res;
    res.push_back(obj);
    res.append(obj->build(opt, simul));
    obj->fixShape();
    return res;
}

//------------------------------------------------------------------------------

void SolidSet::add(Object * obj)
{
    assert_true(obj->tag() == Solid::TAG);
    ObjectSet::add(obj);
}


void SolidSet::remove(Object * obj)
{
    ObjectSet::remove(obj);
    simul.singles.removeWrists(obj);
}


void SolidSet::foldPosition(const Modulo * s) const
{
    for ( Solid * o=first(); o; o=o->next() )
        o->foldPosition(s);
}


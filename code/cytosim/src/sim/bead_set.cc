// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "bead_set.h"
#include "bead_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "single_prop.h"
#include "wrist.h"
#include "simul.h"


Property* BeadSet::newProperty(const std::string& kd, const std::string& nm, Glossary&) const
{
    if ( kd == "bead" )
        return new BeadProp(kd, nm);
    return 0;
}


Object * BeadSet::newObjectT(const Tag tag, unsigned idx)
{
    if ( tag == Bead::TAG )
    {
        BeadProp * p = simul.findProperty<BeadProp*>("bead", idx);
        if ( p == 0 )
            throw InvalidIO("no bead class defined with id "+sMath::repr(idx));
        return new Bead(p, Vector(0,0,0), 0);
    }
    return 0;
}

/**
 @ingroup NewObject

 By definition, a Bead has one point, and one can only set the radius of the Bead:

 @code
 new bead NAME
 {
   radius = REAL
 }
 @endcode

 <h3> How to add Single </h3>

 Singles can only be attached at the center of the Bead:

 @code
 new bead NAME
 {
   radius = REAL
   attach = SINGLE_SPEC [, SINGLE_SPEC] ...
 }
 @endcode
 
 Where `SINGLE_SPEC` is string containing at most 3 words: `[INTEGER] NAME`,
 where the `INTEGER` specifies the number of Singles, `NAME` specifies their name.
 
 For example if `grafted` is the name of a Single, one can use:
 
 @code
 new bead NAME
 {
   attach = 10 grafted
 }
 @endcode

 */

ObjectList BeadSet::newObjects(const std::string& name, Glossary& opt)
{
    real rad = -1;
    if ( !opt.set(rad, "radius") || rad <= 0 )
        throw InvalidParameter("bead:radius must be specified and > 0");
        
    // possibly add some variability
    real var = 0;
    if ( opt.set(var, "radius", 1) )
    {
        real drad;
        do
            drad = var * RNG.gauss();
        while ( rad + drad < REAL_EPSILON );
        rad += drad;
    }
    
    Property * p = simul.properties.find_or_die("bead", name);
    Bead * obj = new Bead(static_cast<BeadProp*>(p), Vector(0,0,0), rad);
    
    ObjectList res;
    res.push_back(obj);
    
    std::string str;
    unsigned inx = 0;
    // attach different kinds of SINGLE
    while ( opt.set(str, "attach", inx++) )
        res.append(simul.singles.makeWrists(obj, 0, 1, str));

    return res;
}


void BeadSet::remove(Object * obj)
{
    ObjectSet::remove(obj);
    simul.singles.removeWrists(obj);
}


void BeadSet::foldPosition(const Modulo * s) const
{
    for ( Bead * o=first(); o; o=o->next() )
        o->foldPosition(s);
}

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "glossary.h"
#include "point_exact.h"
#include "organizer_set.h"
#include "nucleus.h"
#include "bundle.h"
#include "aster.h"
#include "fake.h"
#include "solid.h"
#include "simul.h"

//------------------------------------------------------------------------------

void OrganizerSet::step()
{
    for ( Organizer * as=first(); as ; as=as->next() )
        as->step();
}

//------------------------------------------------------------------------------

Property* OrganizerSet::newProperty(const std::string& kd, const std::string& nm, Glossary&) const
{
    if ( kd == "aster" )   return new AsterProp(nm);
    if ( kd == "bundle" )  return new BundleProp(nm);
    if ( kd == "nucleus" ) return new NucleusProp(nm);
    if ( kd == "fake" )    return new FakeProp(nm);
    return 0;
}


Object * OrganizerSet::newObjectT(const Tag tag, unsigned idx)
{
    if ( tag == Aster::TAG )
    {
        AsterProp * p = simul.findProperty<AsterProp*>("aster", idx);
        if ( p == 0 )
            throw InvalidIO("no aster class defined with id "+sMath::repr(idx));
        return new Aster(p);
    }
    
    if ( tag == Bundle::TAG )
    {
        BundleProp * p = simul.findProperty<BundleProp*>("bundle", idx);
        if ( p == 0 )
            throw InvalidIO("no bundle class defined with id "+sMath::repr(idx));
        return new Bundle(p);
    }
    
    if ( tag == Nucleus::TAG )
    {
        NucleusProp * p = simul.findProperty<NucleusProp*>("nucleus", idx);
        if ( p == 0 )
            throw InvalidIO("no nucleus class defined with id "+sMath::repr(idx));
        return new Nucleus(p);
    }
    
    if ( tag == Fake::TAG )
    {
        FakeProp * p = simul.findProperty<FakeProp*>("fake", idx);
        if ( p == 0 )
            throw InvalidIO("no fake class defined with id "+sMath::repr(idx));
        return new Fake(p);
    }
    
    throw InvalidIO("unknown organizer TAG `"+std::string(1,tag)+"'");
    return 0;
}


ObjectList OrganizerSet::newObjects(const std::string& name, Glossary& opt)
{
    Organizer * obj = 0;
    Property * p = simul.properties.find_or_die(name);
    
    if ( p->category() == "aster" )
        obj = new Aster(static_cast<AsterProp*>(p));
    else if ( p->category() == "bundle" )
        obj = new Bundle(static_cast<BundleProp*>(p));
    else if ( p->category() == "nucleus" )
        obj = new Nucleus(static_cast<NucleusProp*>(p));
    else if ( p->category() == "fake" )
        obj = new Fake(static_cast<FakeProp*>(p));

    ObjectList res;
    if ( obj )
    {
        res = obj->build(opt, simul);
        res.push_back(obj);
    }
    
    return res;
}

//------------------------------------------------------------------------------

void OrganizerSet::add(Object * obj)
{
    ObjectSet::add(obj);
    // we also link all dependent objects:
    static_cast<Organizer*>(obj)->addOrganized(simul);
}


Aster * OrganizerSet::findAster(const ObjectID n) const
{
    Object * obj = findID(n);
    if ( obj  &&  obj->tag() == Aster::TAG )
        return static_cast<Aster*>(obj);
    return 0;
}


Organizer * OrganizerSet::findOrganizer(const Mecable * m) const
{
    for ( Organizer * o=first(); o; o=o->next() )
        for ( unsigned i = 0; i < o->nbOrganized(); ++i )
            if ( m == o->organized(i) )
                return o;

    return 0;
}


//------------------------------------------------------------------------------
void OrganizerSet::foldPosition(const Modulo * s) const
{
    for ( Organizer * o=first(); o; o=o->next() )
        o->foldPosition(s);
}


void OrganizerSet::report(std::ostream& os) const
{
    if ( size() > 0 )
    {
        os << title() << "\n";
        PropertyList plist = simul.properties.find_all("aster", "bundle");
        for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
        {
            ObjectList olist = collect(match_property, *ip);
            os << std::setw(10) << olist.size() << " " << (*ip)->name() << '\n';
        }
        plist = simul.properties.find_all("nucleus", "fake");
        for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
        {
            ObjectList olist = collect(match_property, *ip);
            os << std::setw(10) << olist.size() << " " << (*ip)->name() << '\n';
        }
        if ( plist.size() > 1 )
            os << std::setw(10) << size() << " total\n";
    }
}



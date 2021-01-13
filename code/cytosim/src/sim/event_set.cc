// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "event_set.h"
#include "event_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"

//------------------------------------------------------------------------------

void EventSet::prepare()
{
    for ( Event * f=first(); f; f=f->next() )
    {
        f->prepare();
    }
}


void EventSet::step()
{
    for ( Event * f=first(); f; f=f->next() )
    {
        f->step();
    }
}

//------------------------------------------------------------------------------
#pragma mark -

Property* EventSet::newProperty(const std::string& kd, const std::string& nm, Glossary&) const
{
    if ( kd == "event" )
        return new EventProp(nm);
    return 0;
}


Object * EventSet::newObjectT(const Tag tag, unsigned idx)
{
    Event * obj = 0;
    if ( tag == Event::TAG )
    {
        EventProp * p = simul.findProperty<EventProp*>("event", idx);
        if ( p == 0 )
            throw InvalidIO("no event class defined with id "+sMath::repr(idx));
        obj = p->newEvent();
    }
    return obj;
}


/**
 @defgroup NewEvent How to create an Event
 @ingroup NewObject

 Specify a new Event:
 
 @code
 new event NAME
 {

 }
 @endcode
 */
ObjectList EventSet::newObjects(const std::string& name, Glossary& opt)
{
    Property * p = simul.properties.find_or_die(name);
    Event * obj = static_cast<EventProp*>(p)->newEvent();

    ObjectList res;
    res.push_back(obj);
    return res;
}


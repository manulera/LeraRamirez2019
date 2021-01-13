// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_SET_H
#define SPACE_SET_H

#include "object_set.h"
#include "space.h"
class Simul;

///a list of Space
class SpaceSet : public ObjectSet
{
public:
    
    /// creator
    SpaceSet(Simul& s) : ObjectSet(s) {}
    
    /// destructor
    virtual ~SpaceSet() {}
    
    //--------------------------
    
    /// identifies the property
    std::string title() const { return "space"; }
    
    /// create a new property for class `kind` with given name
    Property *  newProperty(const std::string& kind, const std::string& name, Glossary&) const;
    
    /// create new objects, of class `kind` and type `name`, give the options provided in `opt`
    ObjectList  newObjects(const std::string& name, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *    newObjectT(Tag, unsigned);

    //--------------------------
    
    /// add Object
    void add(Object *);
    
    /// remove Object
    void remove(Object *);

    /// erase all Object and all Property
    void erase();
    
    /// Monte-Carlo step for every Space
    void step();
    
    /// first Space
    Space * first() const
    {
        return static_cast<Space*>(nodes.front());
    }

    /// first Space with this Property
    Space * first(const Property * prop) const
    {
        return static_cast<Space*>(ObjectSet::first(prop));
    }
    
    /// last Space
    Space * last() const
    {
        return static_cast<Space*>(nodes.back());
    }

    /// return pointer to the Object of given ID, or zero if not found
    Space * findID(ObjectID n) const
    {
        return static_cast<Space*>(inventory.get(n));
    }
        
};


#endif


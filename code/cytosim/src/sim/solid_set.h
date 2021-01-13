// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SOLID_SET_H
#define SOLID_SET_H

#include "object_set.h"
#include "solid.h"

class Simul;

/// a list of Solid
class SolidSet : public ObjectSet
{
public:
    
    /// creator
    SolidSet(Simul& s) : ObjectSet(s) {}
    
    /// destructor
    virtual ~SolidSet() {}
    
    //--------------------------
    
    /// identifies the class
    std::string title() const { return "solid"; }
    
    /// create a new property for class `kind` with given name
    Property *  newProperty(const std::string& kind, const std::string& name, Glossary&) const;
    
    /// create new objects, of class `kind` and type `name`, give the options provided in `opt`
    ObjectList  newObjects(const std::string& name, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *    newObjectT(Tag, unsigned);
    
    //--------------------------
    
    /// register a Solid into the list
    void        add(Object *);
    
    /// remove from the list
    void        remove(Object *);
    
    /// first Solid
    Solid *     first() const
    {
        return static_cast<Solid*>(nodes.front());
    }
        
    /// return pointer to the Object of given ID, or zero if not found
    Solid *     findID(ObjectID n) const
    {
        return static_cast<Solid*>(inventory.get(n));
    }
    
    /// modulo the position (periodic boundary conditions)
    void        foldPosition(const Modulo *) const;
    
    /// Monte-Carlo simulation step for every Object
    void        step() {}
};


#endif


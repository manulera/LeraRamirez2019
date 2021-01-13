// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef PARTICLE_SET_H
#define PARTICLE_SET_H

#include "object_set.h"
#include "bead.h"

class Simul;

/// a list of Bead
class BeadSet : public ObjectSet
{
public:
    
    /// creator
    BeadSet(Simul& s) : ObjectSet(s) {}
    
    /// destructor
    virtual ~BeadSet() {}
    
    //--------------------------
    
    /// identifies the class
    std::string title() const { return "bead"; }
    
    /// create a new property for class `kind` with given name
    Property *  newProperty(const std::string& kind, const std::string& name, Glossary&) const;
    
    /// create new objects, of class `kind` and type `name`, give the options provided in `opt`
    ObjectList  newObjects(const std::string& name, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *    newObjectT(Tag, unsigned);
    
    //--------------------------
    
    /// remove from the list
    void        remove(Object *);
    
    /// first Object
    Bead *      first() const
    {
        return static_cast<Bead*>(nodes.front());
    }
        
    /// find object from its Number
    Bead *      findID(ObjectID n) const
    {
        return static_cast<Bead*>(inventory.get(n));
    }
    
    /// modulo the position (periodic boundary conditions)
    void        foldPosition(const Modulo *) const;
    
    /// Monte-Carlo simulation step for every Object
    void        step() {}
};


#endif


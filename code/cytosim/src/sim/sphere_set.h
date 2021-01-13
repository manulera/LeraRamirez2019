// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPHERE_SET_H
#define SPHERE_SET_H

#include "object_set.h"
#include "sphere.h"
class Simul;

///a list of Sphere
class SphereSet : public ObjectSet
{
public:
    
    ///creator
    SphereSet(Simul& s) : ObjectSet(s) {}
    
    ///destructor
    virtual ~SphereSet() {}
    
    //--------------------------
    
    /// identifies the class
    std::string title() const { return "sphere"; }
    
    /// create a new property for class `kind` with given name
    Property *  newProperty(const std::string& kind, const std::string& name, Glossary&) const;
    
    /// create new objects, of class `kind` and type `name`, give the options provided in `opt`
    ObjectList  newObjects(const std::string& name, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *    newObjectT(Tag, unsigned);

    //--------------------------
   
    /// remove object
    void        remove(Object *);

    /// first Object
    Sphere *    first() const
    {
        return static_cast<Sphere*>(nodes.front());
    }
    
    /// return pointer to the Object of given ID, or zero if not found
    Sphere *    findID(ObjectID n) const
    {
        return static_cast<Sphere*>(inventory.get(n));
    }
    
    /// modulo the position (periodic boundary conditions)
    void        foldPosition(const Modulo * s) const;
    
    /// Monte-Carlo simulation step for every Object
    void        step() {}
 };

#endif

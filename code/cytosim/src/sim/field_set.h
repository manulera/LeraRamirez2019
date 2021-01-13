// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FIELD_SET_H
#define FIELD_SET_H

#include "object_set.h"
#include "field.h"
class Simul;



/// a list of Field
/**
 FieldSet is not an usual ObjectSet, because Field is a templated class,
 and its instantiations do not derive from a unique parent 'Field'.

 
 Field is defined from templated class FieldBase
 */
class FieldSet : public ObjectSet
{
public:
    
    /// creator
    FieldSet(Simul& s) : ObjectSet(s) {}
    
    /// destructor
    virtual ~FieldSet() {}
    
    //--------------------------
    
    /// identifies the class
    std::string title() const { return "field"; }
    
    /// create a new property for class `kind` with given name
    Property *  newProperty(const std::string& kind, const std::string& name, Glossary&) const;
    
    /// create new objects, of class `kind` and type `name`, give the options provided in `opt`
    ObjectList  newObjects(const std::string& name, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *    newObjectT(Tag, unsigned);
    
    //--------------------------
    
    /// first object
    Field *     first() const
    {
        return static_cast<Field*>(nodes.front());
    }
    
    /// first object
    Field *     first(Property const* p) const
    {
        return static_cast<Field*>(ObjectSet::first(p));
    }
    
    ///  return pointer to the Object of given ID, or zero if not found
    Field *     findID(ObjectID n) const
    {
        return static_cast<Field*>(inventory.get(n));
    }
    
    /// get ready to do a step()
    void        prepare();
    
    /// Monte-Carlo simulation step for every Object
    void        step();

};


#endif

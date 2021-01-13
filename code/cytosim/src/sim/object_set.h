// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef OBJECT_SET_H
#define OBJECT_SET_H

#include "node.h"
#include "object.h"
#include "node_list.h"
#include "inventory.h"
#include <vector>

class Outputter;
class Property;
class PropertyList;
class Glossary;
class Simul;
extern Random RNG;

/// A set of Object
/**
 Encapsulates the different functions used to manage Objects.
 Pointers to the Objects are stored in two lists:
 - a doubly linked list: nodes
 - an array: inventory
 .
 The NodeList nodes is mixed at every time step,
 and thus it can be used to access the objects in a random order,
 as necessary for Monte-Carlo. 
 
 The Inventory can be used to access objects directly.
 
 Functions are used to manage:
 - object creation: newProperty(), newObjects().
 - object lists: size(), add(), remove(), erase().
 - object access: first(), find().
 - simulation: step(), mix().
 - I/O: readObject(), read(), write(), freeze(), thaw().
 .
 */
class ObjectSet
{
private:
    
    ObjectSet();

public:

    /// holds pointers to the Objects organized by ObjectID
    Inventory         inventory;
    
    /// holds pointers to the Objects in a doubly linked list
    NodeList          nodes;
    
    /// the Simul containing this ObjectSet
    Simul&            simul;
        
protected:
    
    /// mark objects from given list
    static void      flag(NodeList const&, unsigned);
    
    /// delete marked objects from given list
    static void      prune(NodeList const&, unsigned);
    
public:
    
    /// mark objects before import
    virtual void      freeze() { flag(nodes, 1); }
    
    /// delete marked objects
    virtual void      prune()  { prune(nodes, 1); }
    
    /// unmark objects after import
    virtual void      thaw()   { flag(nodes, 0); }
    
    
    /// apply translation to all Objects in ObjectList
    static void       translateObjects(ObjectList const&, Vector const&);
    
    /// apply rotation to all Objects in ObjectList
    static void       rotateObjects(ObjectList const&, Rotation const&);
    
    /// rotate all Objects around their position
    static void       revolveObjects(ObjectList const&, Rotation const&);

    /// apply Transformation to all Objects in ObjectList
    static void       moveObjects(ObjectList const&, Isometry const&);

protected:
    
    /// collect all objects
    static ObjectList collect(NodeList const&);

    /// collect objects from NodeList for which func(obj, val) == true
    static ObjectList collect(NodeList const&, bool (*func)(Object const*, void const*), void const*);

    /// write Object in NodeList to file
    static void       write(NodeList const&, Outputter&);
    
public:
    
    /// creator
    ObjectSet(Simul& s) : simul(s) { }
    
    /// destructor
    virtual ~ObjectSet() { erase(); }    
    
    //--------------------------
    
    /// identifies the category of objects stored in this set
    virtual std::string title() const = 0;

    /// create a new property for class `kind` with given name
    virtual Property *  newProperty(const std::string& kind, const std::string& name, Glossary&) const = 0;
    
    /// create new objects, of class `kind` and type `name`, give the options provided in `opt`
    virtual ObjectList  newObjects(const std::string& name, Glossary& opt) = 0;
   
    /// create a non-initialized Object with the corresponding Tag (used for reading trajectory file)
    virtual Object *    newObjectT(Tag, unsigned) = 0;
    
    //--------------------------
    
    /// register Object, and add it at the end of the list
    virtual void       add(Object *);
    
    /// add multiple Objects
    void               add(ObjectList&);
    
    /// remove Object
    virtual void       remove(Object *);

    /// remove all Objects in list
    void               remove(ObjectList&);
    
    /// link the object last in the list
    virtual void       link(Object *);
    
    /// link the object last in the list
    virtual void       unlink(Object *);

    /// remove Object, and delete it
    void               erase(Object *);
    
    /// delete all Objects in list
    void               erase(NodeList&);

    /// delete all Objects in list and forget all serial numbers
    virtual void       erase();
    
    /// number of elements
    virtual unsigned   size()                 const { return nodes.size(); }

    /// mix the order of elements in the doubly linked list nodes
    virtual void       mix()                        { nodes.mix(RNG); }
    
    /// first Object in the list
    Object *           first()                const { return static_cast<Object*>(nodes.front()); }
    
    /// return an Object which has this property
    Object *           first(Property const*) const;
    
    /// last Object
    Object *           last()                 const { return static_cast<Object*>(nodes.back()); }
    
    /// find Object of given serial-number (see Inventoried)
    Object *           findID(const ObjectID n) const { return static_cast<Object*>(inventory.get(n)); }
    
    /// return Object number `n` where negative 'n' designate the end of the list
    Object *           findObject(long n) const;
    
    /// return a list of Objects corresponding to a certain criteria (eg. 'first' or 'last')
    ObjectList         findObjects(std::string spec) const;
    
    /// collect all objects
    virtual ObjectList collect() const;
    
    /// collect objects for which ( func(obj, val) == true )
    virtual ObjectList collect(bool (*func)(Object const*, void const*), void const*) const;

    /// collect objects for which ( obj->property() == prop )
    ObjectList         collect(Property* prop) const;

    //--------------------------
    
    /// print a summary of the content (nb of objects, class)
    virtual void       report(std::ostream&) const;

    /// read one Object from file
    Object *           readObject(Inputter&, Tag, char pretag);

    /// write all Objects to file
    virtual void       write(Outputter&) const;
   
};

#endif

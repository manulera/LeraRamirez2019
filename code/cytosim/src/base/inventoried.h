// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef INVENTORIED_H
#define INVENTORIED_H


/// type for the serial-number of an Object
typedef unsigned long ObjectID;


/// Object that can be entered in a Inventory
/**
 Inventoried provides a serial-number of type ObjectID, used to identify objects in the simulation.
 A serial-number is strictly positive, and it is given only once in each class.
 
 Inventoried [and any derived class] can be registered in a Inventory.
 The Inventory keeps track of all assigned serial-numbers, 
 and can be used to retrieve back an object given its serial-number.
*/
class Inventoried
{
    
    friend class Inventory;
    
protected:
    
    /// object identifier, unique within the class defined by tag()
    ObjectID   ID_;
    
public:
    
    /// set serial number to 0
    Inventoried() : ID_(0) {}
    
    /// passive destructor
    ~Inventoried() {}
    
    
    /// change the serial number
    void     identity(ObjectID n)  { ID_ = n; }
    
    /// returns serial number (strictly positive integer, unique within each class)
    ObjectID identity()      const { return ID_; }
    
};


#endif

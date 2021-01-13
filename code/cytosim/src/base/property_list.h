// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef PROPERTY_LIST_H
#define PROPERTY_LIST_H

#include "assert_macro.h"
#include "property.h"
#include <iostream>
#include <vector>

class Simul;


/// a list of Property
class PropertyList
{
public:
    
    /// type of list used to store pointers to PROP
    typedef std::vector<Property*> vec_type;
    
    /// iterator class type
    typedef vec_type::iterator iterator;

    /// iterator class type
    typedef vec_type::const_iterator const_iterator;
    
private:
    
    /// list of non-null pointers to Properties
    vec_type vec_;
    
public:
    
    /// constructor
    PropertyList()   { }
    
    /// destructor forget things without deleting objects
    ~PropertyList()  { }
    
    //-------------------------------------------------------------------------------
    
    /// add a new Property in the list, and set its index()
    void         deposit(Property * p);
    
    /// push a new Property in the list
    void         push_back(Property * p) { vec_.push_back(p); }

    /// forget pointer to p
    void         remove(Property * p);
    
    /// delete all Property
    void         erase();
    
    /// iterator pointing to first element
    iterator     begin()   { return vec_.begin(); }
    
    /// iterator that points to a position just past the last element
    iterator     end()     { return vec_.end(); }

    /// iterator pointing to first element
    const_iterator begin() const { return vec_.begin(); }
    
    /// iterator pointing to a position past the last element
    const_iterator end()   const { return vec_.end(); }
    
    /// iterator pointing to first element
    Property*    front()   const { return vec_.front(); }
    
    /// iterator pointing to last element
    Property*    back()    const { return vec_.back(); }
    
    //-------------------------------------------------------------------------------
    
    /// true if no property are known
    bool         empty() const    { return vec_.empty(); }
    
    /// number of known Property
    unsigned int size()  const    { return vec_.size(); }

    /// number of Property of given kind
    unsigned int size(std::string const& kind) const;
    
    /// return property stored at index `n`
    Property *   operator[] (size_t n) const;

    /// apply function to all objects
    void         for_each(void func(Property *)) const;
    
    /// complete all objects
    void         complete(Simul const *) const;
    
    //-------------------------------------------------------------------------------
    
    /// return property which has the provided name, or zero if it cannot be found
    Property *   find(std::string const& name) const;
    
    /// return property which has the provided name, or zero if it cannot be found
    Property *   find_or_die(std::string const& name) const;

    /// return property which has the provided name, or zero if it cannot be found
    Property *   find(std::string const& kind, std::string const& name) const;
    
    /// return property which has the provided index, or zero if it cannot be found
    Property *   find(std::string const& kind, const unsigned index) const;

    /// return property which has the provided name, and throw exception if it cannot be found
    Property *   find_or_die(std::string const& kind, std::string const& name) const;

    /// return property which has the provided index, and throw exception if it cannot be found
    Property *   find_or_die(std::string const& kind, const unsigned index) const;
    
    /// return list of properties of the given kind
    PropertyList find_all(std::string const& kind) const;
    
    /// return list of properties matching any of the given kinds
    PropertyList find_all(std::string const& kind1, std::string const& kind2) const;

    /// return list of properties matching any of the given kinds
    PropertyList find_all(std::string const&, std::string const&, std::string const&) const;
    
    /// return a properties of the given kind, that follows the given one
    Property *   find_next(std::string const& kind, Property *) const;
    
    /// return list of properties which are not of the given kind
    PropertyList find_all_except(std::string const& kind) const;

    /// return position where `p` is stored, end() otherwise
    Property const* contains(Property const* p) const;

    //-------------------------------------------------------------------------------

    /// print names of known Property
    void         write_names(std::ostream &, std::string const&) const;
    
    /// write all Property
    void         write(std::ostream &, bool prune) const;
};

#endif


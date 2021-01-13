// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef KEY_LIST_H
#define KEY_LIST_H

#include "assert_macro.h"
#include <string>
#include <vector>
#include <iostream>
#include <sstream>


/// stores a set of pairs ( string, values ). Used by Glossary::set()
template <typename val_type>
class KeyList
{
public:
    
    /// type for a key in KeyList
    typedef std::string key_type;
    
    /// type for a key-value pair
    struct key_value 
    {
        /// key
        key_type   key;
        /// value
        val_type   val;
        /// constructor with initialization
        key_value(key_type const& k, val_type v)
        {
            key=k;
            val=v;
        }
    };
    
private:
    
    /// list of keys-value pairs
    std::vector<key_value> map;
    
public:
    
    /// constructor
    KeyList() { }
    
    /// constructor
    KeyList(key_type k0, val_type v0)
    {
        add(k0, v0);
    }
    
    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1)
    {
        add(k0, v0); add(k1, v1);
    }
    
    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1, key_type k2, val_type v2)
    {
        add(k0, v0); add(k1, v1); add(k2, v2);
    }

    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1, key_type k2, val_type v2,
            key_type k3, val_type v3)
    {
        add(k0, v0); add(k1, v1); add(k2, v2);
        add(k3, v3);
    }
    
    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1, key_type k2, val_type v2,
            key_type k3, val_type v3, key_type k4, val_type v4)
    {
        add(k0, v0); add(k1, v1); add(k2, v2);
        add(k3, v3); add(k4, v4);
    }
    
    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1, key_type k2, val_type v2,
            key_type k3, val_type v3, key_type k4, val_type v4, key_type k5, val_type v5)
    {
        add(k0, v0); add(k1, v1); add(k2, v2);
        add(k3, v3); add(k4, v4); add(k5, v5);
    }
    
    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1, key_type k2, val_type v2,
            key_type k3, val_type v3, key_type k4, val_type v4, key_type k5, val_type v5,
            key_type k6, val_type v6)
    {
        add(k0, v0); add(k1, v1); add(k2, v2);
        add(k3, v3); add(k4, v4); add(k5, v5);
        add(k6, v6);
    }

    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1, key_type k2, val_type v2,
            key_type k3, val_type v3, key_type k4, val_type v4, key_type k5, val_type v5,
            key_type k6, val_type v6, key_type k7, val_type v7)
    {
        add(k0, v0); add(k1, v1); add(k2, v2);
        add(k3, v3); add(k4, v4); add(k5, v5);
        add(k6, v6); add(k7, v7);
    }

    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1, key_type k2, val_type v2,
            key_type k3, val_type v3, key_type k4, val_type v4, key_type k5, val_type v5,
            key_type k6, val_type v6, key_type k7, val_type v7, key_type k8, val_type v8)
    {
        add(k0, v0); add(k1, v1); add(k2, v2);
        add(k3, v3); add(k4, v4); add(k5, v5);
        add(k6, v6); add(k7, v7); add(k8, v8);
    }

    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1, key_type k2, val_type v2,
            key_type k3, val_type v3, key_type k4, val_type v4, key_type k5, val_type v5,
            key_type k6, val_type v6, key_type k7, val_type v7, key_type k8, val_type v8,
            key_type k9, val_type v9)
    {
        add(k0, v0); add(k1, v1); add(k2, v2);
        add(k3, v3); add(k4, v4); add(k5, v5);
        add(k6, v6); add(k7, v7); add(k8, v8);
        add(k9, v9);
    }

    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1, key_type k2, val_type v2,
            key_type k3, val_type v3, key_type k4, val_type v4, key_type k5, val_type v5,
            key_type k6, val_type v6, key_type k7, val_type v7, key_type k8, val_type v8,
            key_type k9, val_type v9, key_type kA, val_type vA)
    {
        add(k0, v0); add(k1, v1); add(k2, v2);
        add(k3, v3); add(k4, v4); add(k5, v5);
        add(k6, v6); add(k7, v7); add(k8, v8);
        add(k9, v9); add(kA, vA);
    }
    
    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1, key_type k2, val_type v2,
            key_type k3, val_type v3, key_type k4, val_type v4, key_type k5, val_type v5,
            key_type k6, val_type v6, key_type k7, val_type v7, key_type k8, val_type v8,
            key_type k9, val_type v9, key_type kA, val_type vA, key_type kB, val_type vB)
    {
        add(k0, v0); add(k1, v1); add(k2, v2);
        add(k3, v3); add(k4, v4); add(k5, v5);
        add(k6, v6); add(k7, v7); add(k8, v8);
        add(k9, v9); add(kA, vA); add(kB, vB);
    }
    
    /// constructor
    KeyList(key_type k0, val_type v0, key_type k1, val_type v1, key_type k2, val_type v2,
            key_type k3, val_type v3, key_type k4, val_type v4, key_type k5, val_type v5,
            key_type k6, val_type v6, key_type k7, val_type v7, key_type k8, val_type v8,
            key_type k9, val_type v9, key_type kA, val_type vA, key_type kB, val_type vB,
            key_type kC, val_type vC)
    {
        add(k0, v0); add(k1, v1); add(k2, v2);
        add(k3, v3); add(k4, v4); add(k5, v5);
        add(k6, v6); add(k7, v7); add(k8, v8);
        add(k9, v9); add(kA, vA); add(kB, vB);
        add(kC, vC);
    }

    
    /// number of entries
    unsigned int size()                     const { return map.size(); }
    
    /// return n-th entry
    const key_value & operator[](size_t n)  const { return map[n]; }
    
    /// add new entry (k, v)
    void add(key_type const& k, val_type v) { map.push_back(key_value(k, v)); }

    /// this will set `val` to be the numerical value corresponding to 'key'
    /**
     Both the key or the ascii representation of its value are accepted
     */
    bool convert(val_type& val, key_type const& key) const
    {        
        std::ostringstream oss;
        for ( unsigned n = 0; n < map.size(); ++n )
        {
            // write in 'oss' the ascii representation of the value
            oss.str("");
            oss << map[n].val;
            // accept both the string key or the ascii-value
            if ( key == map[n].key  ||  key == oss.str() )
            {
                val = map[n].val;
                //std::clog << "KeyList::set  " << key << " -> " << val << std::endl;
                return true;
            }
        }
        return false;
    }
};



extern const char PREF[];


/// output operator
template <typename T>
std::ostream& operator << (std::ostream& os, const KeyList<T> & list)
{
    os << "Known values are:" << std::endl;
    for ( unsigned int n = 0; n < list.size(); ++n )
        os << PREF << list[n].key << " = " << list[n].val << std::endl;
    return os;
}

#endif



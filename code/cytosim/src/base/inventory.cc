// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "inventory.h"
#include "assert_macro.h"
#include "exceptions.h"


Inventory::Inventory()
{
    lowest    = 1;
    highest   = 0;
    allocated = 8;
    byNames   = new Inventoried*[allocated];
    
    for ( ObjectID n = 0; n < allocated; ++n )
        byNames[n] = 0;
}


Inventory::~Inventory()
{
    delete[] byNames;
}


void Inventory::allocate(ObjectID sz)
{
    const unsigned chunk = 32;
    sz = ( sz + chunk - 1 ) & ~( chunk -1 );
    
    Inventoried ** byNames_new = new Inventoried*[sz];
    
    ObjectID n = 0;
    for ( ; n < allocated; ++n )
        byNames_new[n] = byNames[n];
    while ( n < sz )
        byNames_new[n++] = 0;
    
    delete[] byNames;
    byNames   = byNames_new;
    allocated = sz;
}


//------------------------------------------------------------------------------

ObjectID Inventory::first_assigned() const
{
    ObjectID n = 1;
    while ( n < allocated )
    {
        if ( byNames[n] )
            return n;
        ++n;
    }
    return 0;
}


ObjectID Inventory::last_assigned() const
{
    ObjectID n = allocated-1;
    while ( n > 0 )
    {
        if ( byNames[n] )
            return n;
        --n;
    }
    return 0;
}


ObjectID Inventory::next_assigned(ObjectID n) const
{
    ++n;
    while ( n < allocated )
    {
        if ( byNames[n] )
            return n;
        ++n;
    }
    return 0;
}


ObjectID Inventory::first_unassigned()
{
    ObjectID n = lowest;
    
    if ( n < allocated )
    {
        if ( byNames[n] == 0 )
            return n;
    
        while ( n < allocated  &&  byNames[n] )
            ++n;
    
        lowest = n;
    }
    
    return n;
}

//------------------------------------------------------------------------------

/**
 This will assign a new serial-number for `obj`, if it does not have one.
 */
void Inventory::assign(Inventoried * obj)
{
    ObjectID & n = obj->ID_;
    
    if ( n == 0 )
        n = ++highest;
    else if ( highest < n )
        highest = n;
    
    if ( n >= allocated )
        allocate(n+1);
    
    assert_true( byNames[n] == 0 );
    
    byNames[n] = obj;
    //std::err << "Inventory::store() assigned " << n << " to " << obj << "\n";
}


void Inventory::unassign(const Inventoried * obj)
{
    ObjectID n = obj->ID_;
    assert_true( n < allocated );
    byNames[n] = 0;
    
    if ( lowest >= n )
        lowest = n;
    
    while ( byNames[highest] == 0  &&  highest > 0 )
        --highest;
}


Inventoried * Inventory::get(const ObjectID n) const
{
    if ( n < allocated )
    {
        assert_true( byNames[n]==0  ||  byNames[n]->identity()==n );
        return byNames[n];
    }
    return 0;
}


Inventoried* Inventory::first() const
{
    ObjectID n = 1;
    while ( n < allocated )
    {
        if ( byNames[n] )
            return byNames[n];
        ++n;
    }
    return 0;
}


Inventoried* Inventory::last() const
{
    ObjectID n = highest;
    while ( n > 0 )
    {
        if ( byNames[n] )
            return byNames[n];
        --n;
    }
    return 0;
}


Inventoried* Inventory::previous(Inventoried const* i) const
{
    ObjectID n = i->ID_ - 1;
    while ( n > 0 )
    {
        if ( byNames[n] )
            return byNames[n];
        --n;
    }
    return 0;
}

#include <iostream>
Inventoried* Inventory::next(Inventoried const* i) const
{
    ObjectID n = i->ID_ + 1;
    while ( n < allocated )
    {
        if ( byNames[n] )
            return byNames[n];
        ++n;
    }
    return 0;
}

//------------------------------------------------------------------------------
unsigned Inventory::count() const
{
    unsigned cnt = 0;
    for ( ObjectID n = 0; n < allocated; ++n )
        if ( byNames[n] ) ++cnt;
    return cnt;
}


void Inventory::reassign()
{
    ObjectID max = last_assigned();
    ObjectID next = 1;
    ObjectID nn   = 1;
    
    while ( nn <= max )
    {
        while ( nn <= max  &&  byNames[nn] == 0 )
            ++nn;
        if ( nn > max )
            break;
        if ( next < nn )
        {
            byNames[next] = byNames[nn];
            byNames[nn]   = 0;
            byNames[next]->identity(next);
        }
        ++next;
        ++nn;
    }
    
    lowest = next;
    highest = next-1;
}


void Inventory::clear()
{
    for ( ObjectID n = 0; n < allocated; ++n )
        byNames[n] = 0;
    //std::clog << "Inventory::forgetAll() removed " << cnt << "numbers\n";
    lowest = 1;
    highest = 0;
}


//------------------------------------------------------------------------------
std::ostream& operator << (std::ostream& os, Inventory const& inv)
{
    os << "Inventory " << &inv << std::endl;
    for ( ObjectID n = 0; n < inv.capacity(); ++n )
        os << n << " -> " << inv[n] << std::endl;
    return os;
}


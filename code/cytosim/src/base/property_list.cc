// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "property_list.h"
#include "exceptions.h"
#include <iomanip>

void PropertyList::erase()
{
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
        delete( *i );
    vec_.clear();
}

/**
 Add a property to the list. If ( p == 0 ) nothing is done.

 This function sets the index of `p` to follow the Properties of the same kind,
 that are already present in the list.
 */
void PropertyList::deposit(Property * p)
{
    if ( p )
    {
        int cnt = 0;
        for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
        {
            if ( (*i)->category() == p->category() )
                ++cnt;
            if ( (*i)->name() == p->name() )
                throw InvalidParameter("Property '"+p->name()+"' is already defined");
        }
        
        //std::clog << "Property `" << p->name() << "' is " << p->category() << " # " << cnt+1 << std::endl;
        
        vec_.push_back(p);
        p->reindex(cnt+1);
    }
}


/**
 The size of the array will be reduced by one
 */
void PropertyList::remove(Property * p)
{
    for ( iterator i = vec_.begin(); i != vec_.end(); ++i )
    {
        if ( *i == p )
        {
            vec_.erase(i);
            return;
        }
    }
}


unsigned int PropertyList::size(std::string const& cat) const
{
    unsigned res = 0;
    
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
        if ( (*i)->category() == cat )
            ++res;
    
    return res;
}


Property * PropertyList::operator[] (const size_t n) const
{
    if ( n >= vec_.size() )
    {
        std::ostringstream oss;
        oss << "out of range index " << n << " ( list-size = " << vec_.size() << " )";
        throw InvalidSyntax(oss.str());
    }
    return vec_[n];
}

//-------------------------------------------------------------------------------

void PropertyList::for_each(void func(Property *)) const
{
    //std::clog << "Running function for "<<vec_.size()<<" properties"<<std::endl;
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
        func(*i);
}

void PropertyList::complete(Simul const* sim) const
{
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
        (*i)->complete(sim);
}

Property const* PropertyList::contains(Property const* p) const
{
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
        if ( *i == p ) return p;
    return 0;
}


//-------------------------------------------------------------------------------
#pragma mark -

/** 
 returns the first match
 */
Property * PropertyList::find(std::string const& nm) const
{
    //std::clog << this << "->find(" << nm << ")" << std::endl;
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
    {
        if ( (*i)->name() == nm )
            return *i;
    }
    
    return 0;
}


/**
 returns the first match
 */
Property * PropertyList::find_or_die(std::string const& nm) const
{
    //std::clog << this << "->find_or_die(" << nm << ")" << std::endl;
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
    {
        if ( (*i)->name() == nm )
            return *i;
    }

    std::ostringstream oss;
    oss << "Unknown Property `" << nm << "'\n";
    write_names(oss, PREF);
    throw InvalidSyntax(oss.str());
    return 0;
}


/**
 returns the first match
 */
Property * PropertyList::find(std::string const& cat, std::string const& nm) const
{
    //std::clog << this << "->find(" << cat << ", " << nm << ")" << std::endl;

    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
    {
        if ( (*i)->category()==cat  &&  (*i)->name()==nm )
            return *i;
    }
    
    return 0;
}


Property * PropertyList::find(std::string const& cat, const unsigned idx) const
{
    //std::clog << this << "->find(" << cat << ", " << idx << ")" << std::endl;
    if ( idx <= 0 )
        return 0;
    
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
        if ( (*i)->category()==cat  &&  (*i)->index()==idx )
            return *i;
    
    return 0;
}


Property * PropertyList::find_or_die(std::string const& cat, std::string const& nm) const
{
    Property * res = find(cat, nm);
    
    if ( !res )
    {
        std::ostringstream oss;
        oss << "Unknown " << cat << " `" << nm << "'\n";
        write_names(oss, PREF);
        throw InvalidSyntax(oss.str());
    }
    
    return res;
}


Property * PropertyList::find_or_die(std::string const& cat, const unsigned idx) const
{
    Property * res = find(cat, idx);
    
    if ( !res )
    {
        std::ostringstream oss;
        oss << "Unknown " << cat << "(" << idx << ")\n";
        write_names(oss, PREF);
        throw InvalidSyntax(oss.str());
    }
    
    return res;
}


PropertyList PropertyList::find_all(std::string const& cat) const
{
    //std::clog << this << "->find_all(" << cat << ") " << std::endl;

    PropertyList res;
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
    {
        if ( (*i)->category() == cat )
            res.vec_.push_back(*i);
    }
    
    return res;
}


PropertyList PropertyList::find_all(std::string const& c1, std::string const& c2) const
{
    //std::clog << this << "->find_all(" << kd1 << "," << kd2 << ") " << std::endl;
    
    PropertyList res;
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
    {
        if ( (*i)->category() == c1 || (*i)->category() == c2 )
            res.vec_.push_back(*i);
    }
    
    return res;
}


PropertyList PropertyList::find_all(std::string const& c1, std::string const& c2, std::string const& c3) const
{
    //std::clog << this << "->find_all(" << kd1 << "," << kd2 << ") " << std::endl;

    PropertyList res;
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
    {
        if ( (*i)->category() == c1 || (*i)->category() == c2 || (*i)->category() == c3 )
            res.vec_.push_back(*i);
    }
    
    return res;
}


Property* PropertyList::find_next(std::string const& cat, Property * p) const
{
    //std::clog << this << "->find_next(" << cat << ") " << std::endl;
    bool found = ( p == 0 );
    
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
    {
        if ( (*i)->category() == cat )
        {
            if ( found )
                return *i;
            found = ( *i == p );
        }
    }
    
    if ( ! found ) 
        return 0;
    
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
    {
        if ( (*i)->category() == cat )
            return *i;
    }
    
    return 0;
}


PropertyList PropertyList::find_all_except(std::string const& cat) const
{
    PropertyList list;
    //std::clog << this << "->find_all_expect(" << cat << ") " << std::endl;
    
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
    {
        if ( (*i)->category() != cat )
            list.vec_.push_back(*i);
    }
    
    return list;
}


//-------------------------------------------------------------------------------
#pragma mark -

void PropertyList::write_names(std::ostream & os, std::string const& pf) const
{
    os << pf << "Known properties:\n";
    for ( const_iterator i = vec_.begin(); i != vec_.end(); ++i )
    {
        os << pf << std::setw(10);
        if ( *i )
            os << (*i)->category() << "_" << (*i)->index() << " `"<< (*i)->name() << "'";
        else
            os << "void";
        os << std::endl;
    }
}

/**
 The values identical to the default settings are skipped if prune==1
 */
void PropertyList::write(std::ostream & os, const bool prune) const
{
    const_iterator i = vec_.begin();
    
    if ( i != vec_.end() )
    {
        (*i)->write(os, prune);
    
        while ( ++i != vec_.end() )
            (*i)->write(os, prune);
    }
}


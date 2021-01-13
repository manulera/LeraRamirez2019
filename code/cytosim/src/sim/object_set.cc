// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


#include "object_set.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "modulo.h"
#include "space.h"
#include "simul.h"
#include <errno.h>

extern Modulo const* modulo;

//------------------------------------------------------------------------------

/**
 The object is added at the front of the list
 */
void ObjectSet::link(Object * obj)
{
    obj->objset(this);
    nodes.push_front(obj);
}


void ObjectSet::unlink(Object * obj)
{
    assert_true( obj->objset() == this );
    obj->objset(0);
    nodes.pop(obj);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Translate all listed movable objects ( Object::mobile()==true ) by `vec`
 */

void ObjectSet::translateObjects(ObjectList const& objs, Vector const& vec)
{
    for ( ObjectList::iterator oi = objs.begin(); oi < objs.end(); ++oi )
        (*oi)->flag(1);
    
    for ( ObjectList::iterator oi = objs.begin(); oi < objs.end(); ++oi )
    {
        Object * mv = *oi;
        if ( mv->mobile() && mv->flag() )
        {
            mv->translate(vec);
            mv->flag(0);
        }
    }
}

/**
 Apply Rotation around the origin to all movable objects in list
 */

void ObjectSet::rotateObjects(ObjectList const& objs, Rotation const& rot)
{
    for ( ObjectList::iterator oi = objs.begin(); oi < objs.end(); ++oi )
        (*oi)->flag(1);
    
    for ( ObjectList::iterator oi = objs.begin(); oi < objs.end(); ++oi )
    {
        Object * mv = *oi;
        if ( mv->mobile() && mv->flag() )
        {
            mv->rotate(rot);
            mv->flag(0);
        }
    }
}

/**
 Rotate all objects in `objs` around their position
 */

void ObjectSet::revolveObjects(ObjectList const& objs, Rotation const& rot)
{
    for ( ObjectList::iterator oi = objs.begin(); oi < objs.end(); ++oi )
    {
        Movable * mv = *oi;
        if ( mv->mobile() )
            mv->revolve(rot);
    }
}

/** 
Apply isometry to all objects
 */

void ObjectSet::moveObjects(ObjectList const& objs, Isometry const& iso)
{
    //std::clog << "moving " << objs.size() << " objects" << std::endl;
    for ( ObjectList::iterator oi = objs.begin(); oi < objs.end(); ++oi )
    {
        Movable * mv = *oi;
        if ( mv->mobile() )
        {
            //std::clog << "    moving " << (*oi)->reference() << std::endl;
            mv->rotate(iso.rotation());
            mv->translate(iso.translation());
        }
        //else std::clog << "    cannot move " << (*oi)->reference() << std::endl;
    }
}


//------------------------------------------------------------------------------
#pragma mark -

void ObjectSet::add(Object * obj)
{
    if ( !obj->linked() )
    {
        inventory.assign(obj);
        link(obj);
        //std::clog << "ObjectSet::add(" << obj->reference() << ")\n";
    }
    else
    {
        std::clog << "Warning: attempted to re-link "+obj->reference()+" \n";
    }
}


void ObjectSet::add(ObjectList & list)
{
    for ( ObjectList::iterator oi = list.begin(); oi < list.end(); ++oi )
        add(*oi);
}


void ObjectSet::remove(Object * obj)
{
    //std::clog << "ObjectSet::remove(" <<  obj->reference() << ")\n";
    inventory.unassign(obj);
    if ( obj->linked() )
        unlink(obj);
}


void ObjectSet::remove(ObjectList & list)
{
    for ( ObjectList::iterator oi = list.begin(); oi < list.end(); ++oi )
        remove(*oi);
}


void ObjectSet::erase(NodeList & list)
{
    Node * n = list.front();
    while ( n )
    {
        Node * p = n->next();
        list.pop(n);
        static_cast<Object*>(n)->objset(0);
        delete(n);
        n = p;
    }
}


void ObjectSet::erase(Object * obj)
{
    remove(obj);
    delete(obj);
}


void ObjectSet::erase()
{
    erase(nodes);
    inventory.clear();
}



/**
 Return object as follows:
 - `if ( n >  0 )`, return object with serial-number `n`.
 - `if ( n <= 0 )`, return object from the end of the list:
     -  0 is the last object
     - -1 is the penultimate
     - etc.
 .
 zero is return if object is not found
 */

Object * ObjectSet::findObject(long num) const
{
    Inventoried * res = 0;
    if ( num > 0 )
        res = inventory.get(num);
    else
    {
        res = inventory.last();
        while ( res  &&  ++num <= 0 )
            res = inventory.previous(res);
    }
    return static_cast<Object*>(res);
}

bool match_all(Object const*, void const*)
{
    return true;
}

/*
 There are several ways to designate an object. 
 For example, if the class name is 'fiber', one may use:
 - `fiber1`  for the first fiber ever created
 - `fiber2`  for the second ever created, etc.
 - `first`   indicates the oldest fiber remaining
 - `first+1` indicates the second oldest fiber remaining
 - `last`    indicates the last fiber created
 - `last-1`  indicates the last fiber created
 - `fiber-1` the penultimate fiber, etc.
 .
 */
ObjectList ObjectSet::findObjects(std::string spec) const
{
    ObjectList res;
    Object * obj = 0;
    //std::clog << "ObjectSet::findObjects " << str << std::endl;
    
#ifdef BACKWARD_COMPATIBILITY
    // an integer can indicate a serial number
    if ( isdigit(spec[0]) )
    {
        long num = strtol(spec.c_str(), 0, 10);
        obj = findObject(num);
    }
#endif
    
    // split into a word and a number:
    long num = 0;
    size_t pos = spec.find_first_of("0123456789+-");
    if ( pos != std::string::npos )
    {
        errno=0;
        std::string str = spec.substr(pos);
        num = strtol(str.c_str(), 0, 10);
        if ( errno )
            throw InvalidParameter("expected a number in `"+str+"'");
        spec.resize(pos);
        //std::cerr << "findObjects() SPLIT " << spec << "|" << num << std::endl;
    }
    
    // check for a string starting with the class name (eg. 'fiber'):
    if ( spec == "all" || spec == "any" )
    {
        res = collect(match_all, 0);
    }
    else if ( spec == title() )
    {
        if ( num )
            obj = findObject(num);
        else
            res = collect(match_all, 0);
    }
    // check if string starts with 'first'
    else if ( spec == "first" )
    {
        Inventoried* inv = inventory.first();
        while ( inv  &&  --num >= 0 )
            inv = inventory.next(inv);
        obj = static_cast<Object*>(inv);
    }
    // check if string starts with 'last'
    else if ( spec == "last" )
    {
        Inventoried* inv = inventory.last();
        while ( inv  &&  ++num <= 0 )
            inv = inventory.previous(inv);
        obj = static_cast<Object*>(inv);
    }
    // check finally all the objects:
    else
    {
        // brute force search
        for ( Object* i = first(); i; i = i->next() )
        {
            if ( spec == i->property()->name() || spec == i->property()->category() )
            {
                if ( num == 0 || num == i->identity() )
                    res.push_back(i);
            }
        }
    }
    
    if ( obj )
        res.push_back(obj);
    
    //std::cerr << title() << "::findObject() returns " << res.size() << " object(s)" << std::endl;

    return res;
}


/**
 return the first object encountered with the given property,
 but it can be any one of them, since the lists are regularly
 shuffled to randomize the order in the list.
 */
Object * ObjectSet::first(Property const* prop) const
{
    for ( Object* obj=first(); obj; obj=obj->next() )
        if ( obj->property() == prop )
            return obj;
    return 0;
}



ObjectList ObjectSet::collect(const NodeList & list)
{
    ObjectList res;
    for ( Node* n = list.front(); n; n=n->next() )
        res.push_back(static_cast<Object*>(n));
    return res;
}


ObjectList ObjectSet::collect(const NodeList & list,
                              bool (*func)(Object const*, void const*), void const* arg)
{
    ObjectList res;
    Node * n = list.front();
    while ( n )
    {
        Object * obj = static_cast<Object*>(n);
        n = n->next();
        if ( func(obj, arg) )
            res.push_back(obj);
    }
    return res;
}


ObjectList ObjectSet::collect() const
{
    return collect(nodes);
}


ObjectList ObjectSet::collect(bool (*func)(Object const*, void const*), void const* arg) const
{
    return collect(nodes, func, arg);
}


ObjectList ObjectSet::collect(Property * prop) const
{
    return collect(match_property, prop);
}


//------------------------------------------------------------------------------
#pragma mark - I/O


void ObjectSet::flag(NodeList const& list, unsigned f)
{
    for ( Node const* n=list.front(); n; n=n->next() )
    {
        Object const* o = static_cast<Object const*>(n);
        o->flag(f);
    }
}


void ObjectSet::prune(NodeList const& list, unsigned f)
{
    Node const* n = list.front();
    
    while ( n )
    {
        Node const* p = n->next();
        Object const* o = static_cast<Object const*>(n);
        if ( o->flag() == f )
            delete(o);
        n = p;
    }
}


/**
 Export all objects to file
 */
void ObjectSet::write(Outputter & out) const
{
    if ( size() > 0 )
    {
        out.put_line("#section "+title());
        write(nodes, out);
    }
}


/**
 Write Reference and object-data, for all Objects in 'nodes'
 */
void ObjectSet::write(NodeList const& list, Outputter & out)
{
    for ( Node const* n=list.front(); n; n=n->next() )
    {
        const Object * o = static_cast<const Object*>(n);
        out.put_char('\n');
        o->writeReference(out);
        o->write(out);
    }
}


/**
 Load an object from file, overwritting the current object if it is found
 */
Object * ObjectSet::readObject(Inputter& in, const Tag tag, const char pretag)
{
    unsigned ix = 0, mk = 0;
    ObjectID nb = 0;
    
    Object::readReference(in, ix, nb, mk, pretag);
    
    if ( nb == 0 )
        throw InvalidIO("Invalid (null) object reference");
    
    Object * w = findID(nb);
    
    if ( w )
    {
        //std::clog << "modifying " << w->reference() << std::endl;
        assert_true( w->identity() == nb );
        assert_true( w->property() );
        assert_true( w->linked() );
        w->mark(mk);
        try {
            w->read(in, simul, tag);
        }
        catch( Exception & e ) {
            e << ", while importing " << w->reference();
            throw;
        }
    }
    else
    {
#if ( 0 )
        fpos_t pos;
        in.get_pos(pos);
        std::clog << "- reading `" << (char)tag << sMath::repr(ix) << ":" << nb << "' at " << pos << std::endl;
#endif
        
        w = newObjectT(tag, ix);
        if ( w == 0 )
        {
            InvalidIO e("unknown object `");
            e << (char)tag << sMath::repr(ix) << ":" << nb << "' referenced in file";
            throw e;
        }
        
        w->identity(nb);
        w->mark(mk);
        //std::clog << "created " << w->reference() << std::endl;

        try {
            w->read(in, simul, tag);
        }
        catch( Exception & e ) {
            e << ", while importing " << w->reference();
            throw;
        }
    }
    return w;
}


//------------------------------------------------------------------------------


void ObjectSet::report(std::ostream& os) const
{
    if ( size() > 0 )
    {
        os << title() << '\n';
        PropertyList plist = simul.properties.find_all(title());
        for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
        {
            ObjectList olist = collect(match_property, *ip);
            os << std::setw(10) << olist.size() << " " << (*ip)->name() << '\n';
        }
        if ( plist.size() > 1 )
            os << std::setw(10) << size() << " total\n";
    }
}


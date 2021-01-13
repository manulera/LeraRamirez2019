// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef OBJECT_H
#define OBJECT_H

#include "inventoried.h"
#include "movable.h"
#include "random.h"
#include "array.h"
#include "node.h"

class Simul;
class Property;
class Inputter;
class Outputter;
class Display;
class Meca;

extern Random RNG;


/// Type for unique class identifier used to read/write from objects from file
typedef int Tag;

/// Parent class for all simulated objects
/**
 This is the interface used for writing / reading from a file.
 
 Three functions identify an Object:
 - tag() [ASCII character] identifies the class of Object.
 - property->index() [integer] identifies its Property.
 - identity() [serial-number] derived from Inventoried identifies unique instantiations.
 .
 These three qualities are concatenated in reference() and writeReference().
 
 Objects are stored in ObjectSet.
 */
class Object : public Node, public Inventoried, public Movable
{
    
private:
    
    /// integer used for custom tasks, which is recorded to file
    unsigned          mark_;
    
    /// another integer used for temporary tasks, not saved to file
    mutable unsigned  flag_;
    
    /// a random number associated with this object
    unsigned     signature_;
    
    /// upstream pointer to container class
    ObjectSet *        set_;
    
public:
    
    /// Object::TAG = 'v' represents the 'void' pointer
    static const Tag TAG = 'v';
    
    /// build a string reference by concatenating (tag, property_index, number, mark)
    static std::string strReference(Tag, int, ObjectID, int);

    /// read a reference (property_index, number, mark) from input
    static void        readReference(Inputter&, unsigned&, ObjectID&, unsigned&, Tag pretag);
    
public:
    
    /// constructor
    Object() : mark_(0), flag_(0), set_(0) { signature_ = RNG.pint(); }
    
    /// destructor
    ~Object();
    
    /// a character identifying the class of this object
    virtual Tag     tag() const { return TAG; }
    
    /// Property associated with the Object
    virtual Property const* property() const = 0;
    
    /// write Object data to file
    virtual void    write(Outputter&) const = 0;
    
    /// read Object from file, within the Simul
    virtual void    read(Inputter&, Simul&, Tag) = 0;
    
    //--------------------------

    /// returns container class
    ObjectSet *     objset() const { return set_; }
    
    /// change container class
    void            objset(ObjectSet* s) { set_ = s; }
    
    /// true if Node is registered in a container class
    bool            linked() const { return set_ != 0; }

    /// concatenation of [ tag(), property()->index(), identity() ] in plain ascii
    std::string     reference() const;

    /// write a reference, but using the provided Tag
    void            writeReference(Outputter &, Tag tag) const;
    
    /// write a reference that identifies the Object uniquely
    void            writeReference(Outputter & ow) const { writeReference(ow, tag()); }
    
    /// write a reference that does not refer to any Object
    static void     writeNullReference(Outputter &);
    
    //--------------------------

    /// get mark
    unsigned        mark()            const { return mark_; }
    
    /// set mark
    void            mark(unsigned m)        { mark_ = m; }
    
    
    /// retreive flag value
    unsigned        flag()           const  { return flag_; }
    
    /// set flag (this value is not stored in trajectory files)
    void            flag(unsigned f) const  { flag_ = f; }

    
    /// a random number that makes objects unique
    unsigned        signature()       const { return signature_; }
    
    /// set signature
    void            signature(unsigned s)   { signature_ = s; }

    //--------------------------

    /// extends Node::next(), with a cast to preserve type
    Object *        next()          const { return static_cast<Object*>(nNext); }
    
    /// extends Node::prev(), with a cast to preserve type
    Object *        prev()          const { return static_cast<Object*>(nPrev); }
};



/// return 'true' if ( obj->property() == prop )
bool match_property(Object const* obj, void const* prop);


/// a list of Object
typedef Array<Object *> ObjectList;
//typedef std::vector<Object *> ObjectList;


/// printout for debugging purpose
std::ostream& operator << (std::ostream& os, ObjectList const&);


#endif

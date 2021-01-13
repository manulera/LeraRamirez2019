// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "object.h"
#include "iowrapper.h"
#include "exceptions.h"
#include "property.h"
#include "object_set.h"
#include "sim.h"


Object::~Object()
{
    if ( set_ )
    {
        //std::cerr << "unlinking deleted Object " << this << '\n';
        set_->remove(this);
    }
}


/**
 This is an accessory function used in ObjectSet::collect()
 */
bool match_property(Object const* obj, void const* prop)
{
    return ( obj->property() == prop );
}


/**
 The ASCII reference has the format XP:N or XP:N:M, where:
 - X is one ascii character, the class Tag.
 - P is a non-null positive integer, the index of the property.
 - N is a non-null positive integer, a serial-number of the object within its class.
 - M is an integer, the optional mark of the object (it is added only if non-zero)
 .
 
 For example 'f0:21' is the fiber of property 0, number 21
 For example 'f1:10:1' is the fiber of property 1, number 10, and it is marked
*/
std::string Object::strReference(Tag tag, int pi, ObjectID nb, int mk)
{
    assert_true( pi >= 0 );
    char tmp[32];
    if ( mk == 0 )
        snprintf(tmp, sizeof(tmp), "%c%i:%lu", tag, pi, nb);
    else
        snprintf(tmp, sizeof(tmp), "%c%i:%lu:%i", tag, pi, nb, mk);
    return std::string(tmp);
}


/**
 The ASCII reference has the format XP:N or XP:N:M, where:
 - X is one ascii character, the class Tag.
 - P is a non-null positive integer, the index of the property.
 - N is a non-null positive integer, a serial-number of the object within its class.
 - M is an integer, the mark of the object (this is added only when not null)
 .
 
 For example 'f0:21' is the fiber of property 0, number 21
 For example 'f1:10:1' is the fiber of property 1, number 10, and it is marked as '1'
 */
std::string Object::reference() const
{
    return strReference(tag(), property()->index(), identity(), mark());
}


/**
 Two binary formats are used:
 - A short format:
     - 1 byte for the tag()
     - 1 byte for the index of the property
     - 2 bytes for the Number
     .
 - A long format:
     - the character '$' in binary mode
     - 1 byte for the tag()
     - 2 bytes for the index of the property
     - 4 bytes for the Number
     - 4 bytes for the mark
     .
 .
 There is only one ascii based format, as returned by reference().
 All formats are read by ObjectSet::readReference()
 */
void Object::writeReference(Outputter & out, Tag g) const
{
    assert_true( identity() > 0 );
    assert_true( property() );
    assert_true( property()->index() > 0 );
    
    if ( identity() <= 65535  &&  property()->index() <= 255  &&  mark() == 0 )
    {
        // short format
        out.put_char(g);
        out.writeUInt8(property()->index(), 0);
        out.writeUInt16(identity(), ':');
    }
    else
    {
        // long format with a pretag = '$'
        if ( out.binary() )
            out.put_char('$');
        out.put_char(g);
        out.writeUInt16(property()->index(), 0);
        out.writeUInt32(identity(), ':');
        out.writeUInt32(mark(), ':');
    }
}


void Object::writeNullReference(Outputter & out)
{
    out.put_char(TAG);
}

/**
 This must be able to read the formats written by Object::writeReference()
 */
void Object::readReference(Inputter& in, unsigned& ix, ObjectID& nb, unsigned& mk, const Tag pretag)
{
    ix = 0;
    nb = 0;
    mk = 0;

    if ( in.binary() )
    {
        if ( pretag == '$' )
        {
            // long format
            ix = in.readUInt16();
            nb = in.readUInt32();
#ifdef BACKWARD_COMPATIBILITY
            if ( in.formatID() < 34 )
                return;
            if ( in.formatID() < 39 )
                mk = in.readUInt16();
            else
#endif
            mk = in.readInt32();
        }
        else
        {
            // short format
            ix = in.readUInt8();
            nb = in.readUInt16();
        }
    }
    else
    {
        FILE * file = in.file();
        unsigned u;
        if ( 1 != fscanf(file, " %u", &u) )
            throw InvalidIO("readReference failed");
        ix = u;
        if ( in.get_char() != ':' )
            throw InvalidSyntax("missing ':'");
        if ( 1 != fscanf(file, " %u", &u) )
            throw InvalidIO("readReference failed");
        nb = u;
        int c = in.get_char();
        if ( c == ':' )
        {
            if ( 1 != fscanf(file, " %u", &u) )
                throw InvalidIO("readReference failed");
            mk = u;
        }
        else
            in.unget(c);
    }
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 45 )
        ++ix;
#endif
}



/// print a list of objects
std::ostream& operator << (std::ostream& os, ObjectList const& list)
{
    os << "ObjectList " << &list << " of size " << list.size() << "\n{\n";
    for ( ObjectList::iterator oi = list.begin(); oi < list.end(); ++oi )
        os << "   " << (*oi)->property()->name() << " " << (*oi)->reference() << '\n';
    os << "}" << '\n';
    return os;
}


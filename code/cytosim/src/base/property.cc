// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "property.h"
#include "property_list.h"
#include "glossary.h"
#include "tokenizer.h"
#include "stream_func.h"
#include <sstream>
#include <fstream>

//------------------------------------------------------------------------------

Property::Property(const std::string& n) : name_(n), index_(0)
{
    //std::clog << "new Property `" << mName << "'" << std::endl;
}


Property::~Property()
{
    //std::clog << "del Property `" << mName << "'" << std::endl;
}


//------------------------------------------------------------------------------
/**
 parse string `str` to set values of the property.
*/
void Property::read_string(std::string& str)
{
    if ( str.size() <= 0 )
        return;
    //std::clog << " Property::read_string(" << str << ") for " << name() << std::endl;
    
    Glossary glos;
    glos.read(str);
    read(glos);
}


void Property::read_file(char const* filename)
{
    Glossary glos;
    std::ifstream is(filename);
    glos.read(is);
    read(glos);
}


//------------------------------------------------------------------------------

void Property::write_diff(std::ostream & os, Property const* def) const
{
    if ( def )
    {
        std::stringstream val, ref;
        def->write_values(ref);
        write_values(val);
        StreamFunc::diff_stream(os, val, ref);
    }
    else
        write_values(os);
}


void Property::write_diff(std::ostream & os, const bool prune) const
{
    if ( prune )
    {
        Property * def = clone();
        if ( def )
        {
            def->clear();
            write_diff(os, def);
            delete(def);
            return;
        }
    }
    write_values(os);
}


bool Property::modified() const
{
    std::ostringstream ssr;
    
    Property * def = clone();
    if ( def )
    {
        def->clear();
        def->write_values(ssr);
        std::string str = ssr.str();
        delete(def);
        ssr.str("");
        write_values(ssr);
        return str.compare(ssr.str());
    }
    return true;
}


//------------------------------------------------------------------------------

/**
 This outputs the Properties in this format:
 @code
 set kind name
 {
   property_index = index
   key = values
   ...
 }
 @endcode
 */
void Property::write(std::ostream & os, const bool prune) const
{
    os << "\nset " << category();
    os << " " << name_ << '\n';
    os << "{\n";
    if ( index() > 0 )
        write_value(os, "property_index", index_);
    write_diff(os, prune);
    os << "}\n";
}

std::ostream & operator << (std::ostream& os, const Property& p)
{
    p.write(os, 0);
    return os;
}



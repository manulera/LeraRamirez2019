// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef PROPERTY_H
#define PROPERTY_H

#include "assert_macro.h"
#include <iostream>
#include <iomanip>
#include <string>

class Glossary;
class Simul;


/// A Property holds the parameters for a particular kind of objects
/**
 A Property is a list of parameters associated with a class of objects in Cytosim.
 A Property is identified by:
     - category() indicating the class (eg. `fiber`, `hand`, `single`, etc.)
     - name() which is chosen by the user (eg.`actin`, `microtubule`).
 .
 The `name` should be unique globally, to avoid ambiguity in the config file.
 
 A few important methods handle the most critical operations:
     - clear() will reset parameters to their default values,
     - read() will input parameter from a Glossary
     - complete() will compute derived parameter values, and check their consistency
     - write() will save parameter values to a file.
 .
 
 A Property defacto defines a class of Objects in Cytosim. User-accessible objects 
 have a pointer `prop` to their associated Property. If two objects A and B have the 
 same Property (A.prop == B.prop), they should then be of the the same kind.
 */
class Property
{
private:
    
    /// formatting constant
    static const int FIELD_WIDTH = 20;
    
    /// disabled default constructor:
    Property();
    
    /// the name of the property
    std::string  name_;
    
    /// numerical identifier used in output file
    unsigned     index_;

public:
    
    /// constructor must provide a name
    explicit     Property(const std::string& n);

    /// destructor
    virtual     ~Property();
    
    //-------------------------------------------------------------------------------
    
    /// the 'kind' of property (a class identifier)
    virtual std::string category()        const { return "undefined"; }
    
    //-------------------------------------------------------------------------------
    
    /// return identifier for instantiation
    std::string  name()                   const { return name_; }
    
    /// change name
    void         rename(const std::string& n)   { name_ = n; }
        
    /// true if this->name() is `n`
    bool         is_named(const std::string& n) { return ( n == name_ ); }
    
    //-------------------------------------------------------------------------------
    
    /// index, unique among all Property of similar category()
    unsigned     index()                  const { return index_; }
    
    /// set index in the array of Properties
    void         reindex(unsigned x)            { index_ = x; }
    
    //-------------------------------------------------------------------------------
    
    /// clear parameters to default values
    virtual void clear() = 0;
    
    /// return new object of the same class with identical parameter values
    /**
     The new object is created with `new` and should be destroyed with `delete`
     */
    virtual Property* clone() const = 0;
    
    /// true if at least one value is different from its default setting
    bool         modified() const;
    
    //-------------------------------------------------------------------------------
    
    /// set from a Glossary
    virtual void read(Glossary&) = 0;

    /// set from a string
    void         read_string(std::string&);

    /// read a file specified by name
    void         read_file(const char filename[]);
    
    /// read a file specified by name
    void         read_file(std::string const& str) { read_file(str.c_str()); }
   
    //-------------------------------------------------------------------------------
    
    /// set variables derived from the parameters, and check consistency of values
    /**
     The arguments provide the global SimulProp, and the list of all known Property.
     Any Property created within this function should be added to `plist`.
     complete() is usually called after read()
     */
    virtual void complete(Simul const*) {}
    
    //-------------------------------------------------------------------------------

    /// formatted output of one parameter, one value
    template<typename C>
    static  void write_value(std::ostream& os, std::string const& name, C const& c)
    {
        os << " " << std::left << std::setw(FIELD_WIDTH) << name << " = " << c << ";\n";
    }

    /// formatted output of one array parameter, `cnt` values
    template<typename C>
    static  void write_value(std::ostream& os, std::string const& name, C const* c, int cnt)
    {
        assert_true( cnt > 0 );
        os << " " << std::left << std::setw(FIELD_WIDTH) << name << " = " << c[0];
        for ( int i = 1; i < cnt; ++i )
            os << ", " << c[i];
        os << ";" << std::endl;
    }

    /// formatted output of one parameter, two values
    template<typename C, typename D>
    static  void write_value(std::ostream& os, std::string const& name, C const& c, D const& d)
    {
        os << " " << std::left << std::setw(FIELD_WIDTH) << name << " = " << c << ", " << d << ";" << std::endl;
    }

    /// formatted output of one parameter, three values
    template<typename C, typename D, typename E>
    static  void write_value(std::ostream& os, std::string const& name, C const& c, D const& d, E const& e)
    {
        os << " " << std::left << std::setw(FIELD_WIDTH) << name << " = " << c << ", " << d << ", " << e << ";" << std::endl;
    }

    /// formatted output of one parameter, four values
    template<typename C, typename D, typename E, typename F>
    static  void write_value(std::ostream& os, std::string const& name, C const& c, D const& d, E const& e, F const& f)
    {
        os << " " << std::left << std::setw(FIELD_WIDTH) << name << " = " << c << ", " << d << ", " << e << ", " << f << ";" << std::endl;
    }

    /// write values of all parameters stored in class
    virtual void write_values(std::ostream&) const = 0;
    
    /// write only values that differ from the ones specified in `ref`
    void         write_diff(std::ostream&, Property const* ref) const;

    /// if ( prune == true ), write values that differ from the default values
    void         write_diff(std::ostream&, bool prune) const;
    
    /// write header + data
    void         write(std::ostream&, bool prune = false) const;

};

/// printing operator
std::ostream & operator << (std::ostream&, const Property &);


#endif

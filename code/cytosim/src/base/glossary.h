// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef GLOSSARY_H
#define GLOSSARY_H

#include "key_list.h"
#include "exceptions.h"
#include "assert_macro.h"
#include "tokenizer.h"
#include <iostream>
#include <sstream>
#include <string>
#include <map>

/// Glossary holds a list of (key, values) where both key and values are strings
/** 
 This class is used for reading configuration files:
 - Reads a std::istream to builds a std::map of <key, record>
 - Simple syntax based on = ( ) { } “ , ; % focused on flexible value setting.
 - each `record` is a list of values
 - Provides values upon requests with function set(key, index) .
 - A counter records the usage of the values.
 .
 
 Notes:
-# There can be an arbitrary number of Keys, and an abitrary number of values for each key.
-# Values are kept as strings, and are converted at request by templated functions:
   @code
   template <typename T> int set(T & ptr, std::string key)
   @endcode
   - `key` is the name under which the value appeared,
   - The template argument `T` defines the type of the parameter,
     and the value-string is interpreted accordingly,
   - The interpreted value is stored in `ptr`,
   - Returns 1 if `ptr` was set successfully, 0 otherwise
   .
-# The method warning() can report values that have not been used, 
 or that have been used more than once.
.
 
 Class reviewed by Andre Clapson on 10.03.2011.
 
 @todo parse and instantiate values like 'random.uniform()' or 'PI*30' or '0.1/60'
*/

class Glossary
{
public:
    
    /// type for a key
    typedef std::string                   key_type;
    
    /// a string-encoded value with a counter
    struct val_type 
    {
        /// the value specified as a string
        std::string      value_;
        
        /// true if this value has been propertly set
        bool             defined_;
        
        /// counter for the number of times this value has been read
        mutable unsigned count_;
        
        /// constructor
        val_type()     { defined_=false; count_=0; }
        
        /// constructor with initialization
        val_type(std::string const& s, bool d) { value_=s; defined_=d; count_=0; }
    };
   
    /// a record is a set of values associated with a key
    typedef std::vector<val_type>         rec_type;
    
    /// type for the list of (key, record)
    typedef std::map<key_type, rec_type>  map_type;
    
    /// type of a pair (key, record)
    typedef std::pair<key_type, rec_type> entry_type;
    
    /// add right-hand-side entry to pair.second
    static void  add_value(entry_type&, std::string&, bool);

private:
    
    /// ordered list of key-values pairs
    map_type   mTerms;
    
    //-------------------------------------------------------------------------------
    
    /// true if a character (int) can constitute a value
    static bool  is_value_char(const int);
    
    /// write a value, adding enclosing parenthesis if it contains space characters
    static std::string format_value(const std::string&);
    
    /// read key and assignement operator
    static int   read_key(entry_type&, std::istream&);
    
    /// read one right-hand-side entry of an assignement
    static int   read_value(entry_type&, std::istream&);
    
    /// register a new pair into the dictionnary
    void         add_entry(entry_type&, int no_overwrite);
    
    //-------------------------------------------------------------------------------
    #pragma mark -
    
    /// returns first non-space character in null-terminated C-string
    static char const* not_space(const char s[])
    {
        while(*s)
        {
            if ( isspace(*s) )
                ++s;
            else
                return s;
        }
        return 0;
    }
    
    /// issue a warning if anything other than spaces is left in `iss`
    static void warn_remainder(const char* str, std::string const& msg)
    {
        char const* c = not_space(str);
        if ( c )
        {
            std::cerr << "Warning: ignored trailing `" << c << "' while reading:" << std::endl;
            std::cerr << PREF << msg << std::endl;
        }
    }    
    
    /// set `var` from string `val`
    template <typename T>
    void set_one(T & var, key_type const& key, std::string const& val) const
    {
        if ( val.empty() )
            throw InvalidSyntax("could not set `"+key+"' from empty string");
        
        std::istringstream iss(val);
        
        iss >> var;

        if ( iss.fail() )
            throw InvalidSyntax("could not set `"+key+"' from `"+val+"'");
        
        // check if every characters were used:
        if ( ! iss.eof() )
        {
            std::string str = iss.str();
            char const* chr = str.c_str() + iss.tellg();
            warn_remainder(chr, key + "=" + str);
        }
    }
    
    /// set enum of type T using a dictionary of correspondances
    template <typename T>
    void set_one(T & var, key_type const& key, std::string const& val, KeyList<T> const& dict) const
    {
        if ( ! dict.convert(var, val) )
        {
            std::ostringstream oss;
            oss << "could not set `"+std::string(key)+"' from `"+val+"'\n";
            oss << dict;
            throw InvalidParameter(oss.str());
        }
    }
    
    //-------------------------------------------------------------------------------
    #pragma mark -

public:
    
    /// initialize
    explicit Glossary();

    /// this constructor calls read(in)
    explicit Glossary(std::istream& in);

    /// this constructor calls read()
    explicit Glossary(const std::string&);

    //-------------------------------------------------------------------------------

    /// true if no key were set
    bool         empty()   const { return mTerms.empty(); }

    /// number of keys
    unsigned     nb_keys() const { return mTerms.size(); }
    
    /// return true if key is present, even if no value was set
    bool         has_key(key_type const&) const;
    
    /// return true if key is present, and delete key
    bool         use_key(key_type const&);
    
    /// remove given key
    void         clear(key_type const&);
    
    /// clear all entries
    void         clear() { mTerms.clear(); }
    
    /// clear usage counts for all entries
    void         clear_counts() const;
    
    /// create a new Glossary with only the given key
    Glossary     extract(key_type const&) const;

    /// return number of values associated with a key
    unsigned     nb_values(key_type const&) const;
    
    /// return true if key is present and a value was set for given index
    bool         has_value(key_type const&, unsigned) const;
    
    /// gives a pointer to the values corresponding to a key, or null if the key is not present
    rec_type *   values(key_type const&);

    /// gives a const pointer to the values corresponding to a key, or null if the key is not present
    rec_type const* values(key_type const&) const;
    
    /// return copy of value corresponding to `key[inx]`, or empty string if this value is not present
    std::string  value(key_type const&, unsigned inx) const;
    
    /// returns true if `key[inx]==val`, or false otherwise. Counter is incremented in case of match
    bool         value_is(key_type const& key, unsigned inx, std::string const& val) const;
    
    /// report unused values and values used more than `threshold` times
    static int   warnings(std::ostream&, entry_type const&, unsigned threshold=1);
    
    /// report unused values and values used multiple times
    int          warnings(std::ostream&, unsigned threshold=1) const;
    
    //-------------------------------------------------------------------------------
    #pragma mark -
    
    /// this adds a new key with value 'val': 'key=val'
    void define(key_type const& key, std::string const& val);
    
    /// define one value for the key at specified index: `key[inx]=val`.
    void define(key_type const& key, unsigned inx, std::string const& val);
    
    /// define one value from class T, for the key: `key[inx]=string_representation(val)`.
    template <typename T>
    void define(key_type const& key, unsigned inx, const T& val)
    {
        std::ostringstream oss;
        oss << val;
        define(key, inx, oss.str());
    }

    /// update the glossary to include one assignment read from stream
    void         read_entry(std::istream&, int no_overwrite = 2);

    /// update the glossary to include assignments stored in a stream
    void         read(std::istream&, int no_overwrite = 2);

    /// update the glossary to include assignments stored in a string
    void         read(const std::string&, int no_overwrite = 2);
    
    /// read file specified in path
    void         read_file(const char path[], int no_overwrite = 2);
    
    /// read a file specified by name
    void         read_file(std::string const& str, int no = 2) { read_file(str.c_str(), no); }

    /// a read for C-style command-line arguments
    void         read_strings(int argc, char* argv[], int no_overwrite = 2);
    
    //-------------------------------------------------------------------------------

    /// write [ key = value1, value2, ... ]
    static void  write(std::ostream&, std::string const& prefix, entry_type const&);

    /// write the number of time each value has been used
    static void  write_counts(std::ostream&, std::string const& prefix, entry_type const&);

    /// write all [key, values]
    void         write(std::ostream&, std::string const& prefix = "") const;
    
    //-------------------------------------------------------------------------------
    
    /// set `var` from `key[inx]`. The counter associated to the value is incremented.
    template <typename T>
    int set(T & var, key_type const& key, unsigned inx = 0) const
    {
        rec_type const* rec = values(key);
        
        if ( rec && inx < rec->size() )
        {
            val_type const& val = rec->at(inx);

            if ( val.defined_ )
            {
                set_one(var, key, val.value_);
                ++val.count_;
                return 1;
            }
        }
        
        return 0;
    }

    /// set `var` from `key[inx]`, without recording that the parameter was read.
    template <typename T>
    int peek(T & var, key_type const& key, unsigned inx = 0) const
    {
        rec_type const * rec = values(key);
        
        if ( rec && inx < rec->size() )
        {
            val_type const& val = rec->at(inx);
            
            if ( val.defined_ )
            {
                set_one(var, key, val.value_);
                return 1;
            }
        }

        return 0;
    }
    
    
    /// set `cnt` values in the array `ptr[]`, starting at `key[0]`
    template <typename T>
    int set(T * ptr, unsigned cnt, key_type const& key) const
    {
        rec_type const* rec = values(key);
        
        if ( rec == 0 )
            return 0;
        
        int set = 0;
        for ( unsigned inx = 0; inx < rec->size() && inx < cnt; ++inx )
        {
            val_type const& val = rec->at(inx);
            
            if ( val.defined_ )
            {
                set_one(ptr[inx], key, val.value_);
                ++rec->at(inx).count_;
                ++set;
            }
        }
 
        return set;
    }
   

    /// set `var` from `key[inx]`, using the dictionary `dict`
    template <typename T>
    int set(T & var, key_type const& key, KeyList<T> const& dict, unsigned inx = 0) const
    {
        rec_type const* rec = values(key);
        
        if ( rec && inx < rec->size() )
        {
            val_type const& val = rec->at(inx);
            
            if ( val.defined_ )
            {
                set_one(var, key,  val.value_, dict);
                ++val.count_;
                return 1;
            }
        }
        
        return 0;
    }
    
    
    /// check if value associated with `key` at index `inx` is made only of alpha characters
    int is_alpha(key_type const& key, unsigned inx = 0) const
    {
        rec_type const * rec = values(key);
        
        if ( rec == 0  || inx >= rec->size() )
            return 0;
        
        std::string const& str = rec->at(inx).value_;
        
        if ( str.empty() )
            return 0;
        
        for ( std::string::const_iterator s = str.begin(); s != str.end(); ++s )
        {
            if ( ! isalpha(*s) &&  *s != '_' )
                return 0;
        }
        return 1;
    }
    
    /// check if value associated with `key` at index `inx` is a number
    /**
     @returns:
     - 0 if this is not a number
     - 2 for positive integer
     - 3 for negative integer
     - 4 for positive float
     - 5 for negative float
     .
     */
    int is_number(key_type const& key, unsigned inx = 0) const
    {
        rec_type const * rec = values(key);
        
        if ( rec == 0  || inx >= rec->size() )
            return 0;
        
        std::string const& val = rec->at(inx).value_;
        
        char const* str = val.c_str();
        char * end;
        
        int i = strtol(str, &end, 10);
        if ( end > str && 0 == not_space(end) )
            return 2 + ( i < 0 );
        
        double d = strtod(str, &end);
        if ( end > str && 0 == not_space(end) )
            return 4 + ( d < 0 );
        
        return 0;
    }

};


#pragma mark -


/// special function for std::string arguments.
template <>
void Glossary::set_one(std::string& var, key_type const&, std::string const&) const;

/// special function for float
template <>
void Glossary::set_one(float& var, key_type const&, std::string const&) const;

/// special function for double
template <>
void Glossary::set_one(double& var, key_type const&, std::string const&) const;


/// input from stream
std::istream& operator >> (std::istream&, Glossary&);

/// output of one value
std::ostream & operator << (std::ostream&, const Glossary::entry_type&);

/// output of all values
std::ostream & operator << (std::ostream&, const Glossary&);

#endif



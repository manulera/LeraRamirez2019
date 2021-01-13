// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "glossary.h"
#include "stream_func.h"
#include "ansi_colors.h"
#include "filepath.h"
#include <fstream>
#include <cctype>
#include <iomanip>
#include "vector2.h"


// Use the second definition to get some verbose reports:
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG;

#define VLOG1(ARG) ((void) 0)
//#define VLOG1(ARG) std::clog << ARG;

#define VLOG2(ARG) ((void) 0)
//#define VLOG2(ARG) std::clog << ARG;


//------------------------------------------------------------------------------

Glossary::Glossary()
{
}


Glossary::Glossary(std::istream & in)
{
    read(in); 
}

Glossary::Glossary(const std::string& str)
{
    std::istringstream iss(str);
    read(iss);
}


//------------------------------------------------------------------------------
#pragma mark -

bool Glossary::has_key(key_type const& k) const
{
    return ( mTerms.end() != mTerms.find(k) );
}


bool Glossary::use_key(key_type const& k)
{
    map_type::iterator w = mTerms.find(k);
    
    if ( w != mTerms.end() )
    {
        mTerms.erase(w);
        return true;
    }
    return false;
}


void Glossary::clear(key_type const& key)
{
    map_type::iterator w = mTerms.find(key);
    
    if ( w != mTerms.end() )
        mTerms.erase(w);
}


void Glossary::clear_counts() const
{
    for ( map_type::const_iterator n = mTerms.begin(); n != mTerms.end(); ++n )
    {
        for ( unsigned v = 0; v < n->second.size(); ++v )
            n->second[v].count_ = 0;
    }
}


Glossary Glossary::extract(key_type const& key) const
{
    Glossary res;
    map_type::const_iterator w = mTerms.find(key);
    
    if ( w != mTerms.end() )
        res.mTerms[key] = w->second;
    
    return res;
}

//------------------------------------------------------------------------------
#pragma mark -

unsigned Glossary::nb_values(key_type const& k) const
{
    map_type::const_iterator w = mTerms.find(k);
    if ( w != mTerms.end() )
        return w->second.size();
    else
        return 0;
}


bool Glossary::has_value(key_type const& key, unsigned indx) const
{
    map_type::const_iterator w = mTerms.find(key);
    if ( w != mTerms.end() )
        return indx < w->second.size();
    return false;
}


Glossary::rec_type * Glossary::values(key_type const& key)
{
    map_type::iterator w = mTerms.find(key);
    return ( w == mTerms.end() ) ? 0 : &( w->second );
}


Glossary::rec_type const* Glossary::values(key_type const& key) const
{
    map_type::const_iterator w = mTerms.find(key);
    return ( w == mTerms.end() ) ? 0 : &( w->second );
}


std::string Glossary::value(key_type const& key, unsigned indx) const
{
    map_type::const_iterator w = mTerms.find(key);
    if ( w != mTerms.end() )
    {
        if ( indx < w->second.size() )
        {
            w->second[indx].count_++;
            return w->second[indx].value_;
        }
    }
    return "";
}


/**
 This is equivalement to value(key, indx) == val, except
 that the counter is incremented only if there is a match
 */
bool Glossary::value_is(key_type const& key, unsigned indx, std::string const& val) const
{
    map_type::const_iterator w = mTerms.find(key);
    if ( w != mTerms.end() )
    {
        if ( indx < w->second.size() )
        {
            if ( w->second[indx].value_ == val )
            {
                w->second[indx].count_++;
                return true;
            }
        }
    }
    return false;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 This reads a KEY followed by the assignement operator
 @returns
 - 0 if no valid key is found
 - 1 if the key is immediately followed by the '=' sign
 - 2 otherwise
*/
int Glossary::read_key(Glossary::entry_type& res, std::istream & is)
{
    std::string k = Tokenizer::get_symbol(is, false);

    if ( k.empty() )
        return 0;

    res.first = k;
    
    int op = Tokenizer::get_character(is);
  
    VLOG2("Glossary::   KEY |" << res.first << "|\n");

    if ( op == '=' )
        return 1;

    return 2;
}


/**
 push value at the end of `res.second`
 */
void Glossary::add_value(Glossary::entry_type& res, std::string& str, bool def)
{
    //remove any space at the end of the string:
    Tokenizer::trim(str);
    
    VLOG2("Glossary::" << std::setw(20) << res.first << "[" << res.second.size() << "] = |" << str << "|\n");

    res.second.push_back(val_type(str, def));
}


/**
 return true if `c` can constitute a value.
 space are allowed since vectors are read as sets of space-separated values
*/
bool Glossary::is_value_char(const int c)
{
   return isalnum(c) || c==' ' || c=='/' || c == '#' || c==':' || c=='\t' || c=='_' || c=='.' || c=='+' || c=='-';
}


/**
 read a right-hand side value of an assignment
 return 1 if parsing should continue with the same key
 */
int Glossary::read_value(Glossary::entry_type& res, std::istream & is)
{
    // skip spaces, but do not eat lines
    int c = Tokenizer::get_character(is);
    bool delimited = 0;
    
    std::string k;
    if ( Tokenizer::block_delimiter(c) )
    {
        delimited = true;
        k = Tokenizer::get_block_content(is, 0, Tokenizer::block_delimiter(c));
        c = Tokenizer::get_character(is);
    }
    else
    {
        while ( is_value_char(c) )
        {
            k.push_back(c);
            c = is.get();
        }
    }
    //std::clog << (char)c << "|";
        
    if ( c == EOF || c == '\n' || c == '\r' )
    {
        if ( k.size() || delimited )
            add_value(res, k, true);
        return 0;
    }
    
    if ( c == ';' )
    {
        add_value(res, k, k.size() || delimited);
        return 0;
    }

    if ( c == ',' )
    {
        add_value(res, k, k.size() || delimited);
        return 1;
    }

    if ( c == '%' )
    {
        if ( k.size() || delimited )
            add_value(res, k, true);
        Tokenizer::get_line(is);
        return 0;
    }
    
    if ( c == '\\' )
    {
        if ( k.size() || delimited )
            add_value(res, k, true);
        Tokenizer::get_space(is, true);
        return 1;
    }
    
    is.unget();
    throw InvalidSyntax("Unexpected token `"+std::string(1,c)+"'");
}



/**
 If `no_overwrite == 0`, previous values can be erased without warning,
 If `no_overwrite == 1`, a prexisting symbol cannot be altered, but no exception is thrown
 If `no_overwrite == 2`, an exception is thrown for any duplicate symbol
 */
void Glossary::add_entry(Glossary::entry_type& pair, int no_overwrite)
{
    VLOG("Glossary::ADD " << pair << '\n');
    
    map_type::iterator w = mTerms.find(pair.first);
    
    if ( w == mTerms.end() )
    {
        // for a new key, we accept all values
        rec_type & rec = mTerms[pair.first];
        for ( unsigned v = 0; v < pair.second.size(); ++v )
            rec.push_back(pair.second[v]);
    }
    else
    {
        // for pre-existing keys, we check every values:
        rec_type & rec = w->second;
        for ( unsigned v = 0; v < pair.second.size(); ++v )
        {
            if ( rec.size() <= v )
                rec.push_back(pair.second[v]);
            else
            {
                if ( !rec[v].defined_  ||  !no_overwrite )
                    rec[v] = pair.second[v];
                else if ( pair.second[v].value_ != rec[v].value_  &&  no_overwrite > 1 )
                {
                    std::ostringstream oss;
                    oss << "duplicate value for parameter `" << pair.first << "[" << v << "]':\n";
                    oss << PREF << "previously defined as `" << rec[v].value_ << "'\n";
                    oss << PREF << "new value (ignored)   `" << pair.second[v].value_ << "'\n";
                    throw InvalidSyntax(oss.str());
                }
            }
        }
    }
}


//-------------------------------------------------------------------------------
#pragma mark - Define


/// define one value for the key at specified index: `key[inx]=val`.
void Glossary::define(key_type const& key, unsigned inx, const std::string& val)
{
    map_type::iterator w = mTerms.find(key);
    
    if ( w == mTerms.end() )
    {
        // add new key
        if ( inx > 0 )
            throw InvalidSyntax("index out of range in Glossary::define");
        mTerms[key].push_back(val_type(val, true));

        VLOG1("Glossary::DEFINE " << key << " = |" << val << "|\n");

    }
    else
    {
        rec_type & rec = w->second;
        
        if ( rec.size() > inx )
            rec[inx] = val_type(val, true);
        else if ( rec.size() == inx )
            rec.push_back(val_type(val, true));
        else
            throw InvalidSyntax("index out of range in Glossary::define");
        
        VLOG1("Glossary::DEFINE " << key << "[" << inx << "] = |" << val << "|\n");
    }
}


/**
 This should be equivalent to read('k = rhs')
 */
void Glossary::define(key_type const& k, const std::string& rhs)
{
    define(k, 0, rhs);
}


//------------------------------------------------------------------------------
#pragma mark -

/**
@copydetails Glossary::add_entry
 */

void Glossary::read_entry(std::istream & is, int no_overwrite)
{
    char c = Tokenizer::get_space(is, true);
        
    if ( c == EOF )
        return;
        
    // skip comments:
    if ( c == '%' )
    {
        std::string skip;
        std::getline(is, skip);
        return;
    }
    
    entry_type pair;
    
    int e = read_key(pair, is);
    
    if ( 1 != e )
        throw InvalidParameter("unexpected syntax");

    while ( read_value(pair, is) );
    
    add_entry(pair, no_overwrite);
}


/**
 @copydetails Glossary::add_entry
 */
void Glossary::read(std::istream & is, int no_overwrite)
{
    std::streampos isp;
    std::istream::sentry s(is);
    if ( s ) while ( is.good() )
    {
        isp = is.tellg();
        try
        {
            read_entry(is, no_overwrite);
        }
        catch( Exception& e )
        {
            e << "\n in:\n" << StreamFunc::get_line(is, is.tellg(), PREF);
            throw;
        }
    }
}


void Glossary::read(std::string const& str, int no_overwrite)
{
    VLOG2("Glossary::READ |" << str << "|\n");
    std::istringstream iss(str);
    read(iss, no_overwrite);
}


void Glossary::read_file(const char path[], int no_overwrite)
{
    std::ifstream is(path);
    if ( is.good() )
        read(is, no_overwrite);
    else
        throw InvalidIO("could not open Glossary file");
    is.close();
}


/**
 This is useful to parse the command-line strings given to main().
 
 The following syntax will be accepted:
 FILE.EXT
 and recorded as:
 EXT = FILE.EXT
 
 Strings corresponding to existing directories
 */
void Glossary::read_strings(int argc, char* argv[], int no_overwrite)
{
    for ( int ii = 0; ii < argc; ++ii )
    {
        VLOG("Glossary::ARG |" << argv[ii] << "|\n");
        entry_type pair;
        try
        {
            std::string arg(argv[ii]);
            std::string::size_type spos = arg.rfind("=");
            if ( spos != std::string::npos )
            {
                // a '=' indicates a parameter and its value
                std::istringstream iss(argv[ii]);
                if ( 1 == read_key(pair, iss) )
                {
                    while ( read_value(pair, iss) );
                    add_entry(pair, no_overwrite);
                }
            }
            else
            {
                std::string::size_type spos = arg.rfind(".");
                if ( FilePath::is_dir(arg) )
                {
                    // this is a directory
                    pair.first = "dir";
                    pair.second.push_back(val_type(arg, true));
                }
                else if ( spos != std::string::npos )
                {
                    // a '.' identifies a potential file name
                    pair.first = arg.substr(spos);
                    pair.second.push_back(val_type(arg, true));
                }
                else
                {
                    // anything else is just a string
                    pair.first = arg;
                }
                add_entry(pair, no_overwrite);
            }
        }
        catch( Exception & e )
        {
            e << " in `" << argv[ii] << "'\n";
            throw;
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

std::string Glossary::format_value(std::string const& str)
{
    if ( std::string::npos != str.find(' ') )
        return '(' + str + ')';
    else
        return str;
}


void Glossary::write(std::ostream & os, std::string const& prefix, Glossary::entry_type const& pair)
{
    os << prefix << pair.first;
    if ( pair.second.size() > 0 )
    {
        os << " = " << format_value(pair.second[0].value_);
        for ( unsigned v = 1; v < pair.second.size(); ++v )
            os << ", " << format_value(pair.second[v].value_);
        os << ";";
    }
}

/**
 Write the usage-counter for each value.
 The width of each record will match what is printed by Glossary::write()
 */
void Glossary::write_counts(std::ostream & os, std::string const& prefix, Glossary::entry_type const& pair)
{
    if ( pair.second.size() > 0 )
    {
        os << prefix << std::setw(pair.first.size()) << "used" << " : ";
        os << std::setw(format_value(pair.second[0].value_).size()) << pair.second[0].count_;
        for ( unsigned v = 1; v < pair.second.size(); ++v )
            os << "," << std::setw(format_value(pair.second[v].value_).size()+1) << pair.second[v].count_;
    }
}


void Glossary::write(std::ostream & os, std::string const& prefix) const
{
    for ( map_type::const_iterator n = mTerms.begin(); n != mTerms.end(); ++n )
    {
        write(os, prefix, *n);
        os << std::endl;
    }
}

//------------------------------------------------------------------------------
#pragma mark -

std::istream& operator >> (std::istream & is, Glossary& glos)
{
    glos.read(is);
    return is;
}

std::ostream& operator << (std::ostream& os, Glossary::entry_type const& pair)
{
    Glossary::write(os, "", pair);
    return os;
}

std::ostream& operator << (std::ostream& os, Glossary const& glos)
{
    glos.write(os, "");
    return os;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 @returns:
 - 4 if the parameter was not read
 - 2 if one of the value was not read
 - 1 if some of the value were used multiple times
 - 0 otherwise
 .
 @todo: use color only if the terminal supports it
 */
int Glossary::warnings(std::ostream& os, Glossary::entry_type const& pair, unsigned threshold)
{
    int used = 0, exhausted = 1, overused = 0;
    const rec_type& rec = pair.second;
        
    for ( unsigned v = 0; v < rec.size(); ++v )
    {
        val_type const& val = rec[v];
        if ( val.count_ > 0 )
            used = 1;
        if ( val.count_ == 0 && val.defined_ )
            exhausted = 0;
        else if ( val.count_ > threshold )
            overused = 1;
    }
    
    std::string warn;
    
    if ( !used )
        warn = "this parameter was ignored";
    else if ( !exhausted )
        warn = "a value was unused";
    if ( overused )
        warn = "some value may have been overused";
    
    if ( warn.size() )
    {
        print_magenta(os, "Warning, " + warn + ":\n");
        write(os, PREF, pair);
        os << "\n";
        if ( used )
        {
            write_counts(os, PREF, pair);
            os << "\n";
        }
        os << std::flush;
        
        if ( ! used )
            return 4;
        else if ( ! exhausted )
            return 2;
        return 1;
    }
    return 0;
}


/**
 @returns the number of warnings that were issued
 */

int Glossary::warnings(std::ostream& os, unsigned threshold) const
{
    int res = 0;
    for ( map_type::const_iterator n = mTerms.begin(); n != mTerms.end(); ++n )
        res |= warnings(os, *n, threshold);
    return res;
}

//------------------------------------------------------------------------------

/**
 This copies the string, removing spaces
*/
template <>
void Glossary::set_one(std::string& var, key_type const& key, std::string const& val) const
{
    var = val;
    Tokenizer::trim(var);
    
    VLOG("Glossary::STRING " << key << " = |" << var << "|\n");
}


/**
 This reads a floating point value,
 also accepting 'inf', '+inf' and '-inf' for INFINITY values
 */
template <>
void Glossary::set_one(float& var, key_type const& key, std::string const& val) const
{
/*
    // Infinite values are normally handled by std::strtod()
    if ( val == "inf" || val == "+inf" )
    {
        var = INFINITY;
        return;
    }
    
    if ( val == "-inf" )
    {
        var = -INFINITY;
        return;
    }
*/
    float backup = var;
    char const* str = val.c_str();
    char * end = 0;
    
    var = strtof(str, &end);

    if ( end == str || not_space(end) )
    {
        var = backup;
        throw InvalidSyntax("could not set `"+key+"' from `"+val+"' (expected a scalar value)");
    }
}

/**
 This reads a floating point value,
 also accepting 'inf', '+inf' and '-inf' for INFINITY values
 */
template <>
void Glossary::set_one(double& var, key_type const& key, std::string const& val) const
{
/*
    // Infinite values are normally handled by std::strtod()
    if ( val == "inf" || val == "+inf" )
    {
        var = INFINITY;
        return;
    }
    
    if ( val == "-inf" )
    {
        var = -INFINITY;
        return;
    }
*/
    float backup = var;
    char const* str = val.c_str();
    char * end = 0;
    
    var = strtod(str, &end);
    
    if ( end == str || not_space(end) )
    {
        var = backup;
        throw InvalidSyntax("could not set `"+key+"' from `"+val+"' (expected a scalar value)");
    }
}

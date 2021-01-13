// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by François Nédélec on 1/7/09.

#include "assert_macro.h"
#include "tokenizer.h"
#include "exceptions.h"


// Use the second definition to get some verbose reports:
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG;


//------------------------------------------------------------------------------

/// returns first non-space character in null-terminated C-string
bool Tokenizer::not_space(std::string const& str)
{
    char const* c = str.c_str();
    
    while( *c )
    {
        if ( *c == '\n' )
            break;
        if ( !isspace(*c) )
            return true;
    }
    return false;
}


std::string Tokenizer::get_line(std::istream & is)
{
    std::string res;
    res.reserve(1024);
    std::getline(is, res);
    return res;
}


int Tokenizer::get_space(std::istream & is, bool eat_line)
{
    int c = is.peek();
    while ( isspace(c) )
    {
        if ( c == '\n' && ! eat_line )
            break;
        is.get();
        c = is.peek();
    }
    return c;
}


int Tokenizer::get_character(std::istream & is, bool eat_space, bool eat_line)
{
    int c = 0;
    do {
        c = is.get();
        if ( c == EOF )
            return EOF;
#if ( 0 )
        if ( c == COMMENT_START )
        {
            std::string line;
            std::getline(is, line);
            c = '\n';
        }
#endif
        if ( c == '\n' && ! eat_line )
            break;
    } while ( eat_space && isspace(c) );
    return c;
}

//------------------------------------------------------------------------------
#pragma mark -

std::string get_stuff(std::istream & is, bool (*valid)(int))
{
    std::string res;
    res.reserve(128);
    int c = 0;
    while ( is.good() )
    {
        c = is.peek();
        if ( valid(c) )
            res.push_back(is.get());
        else
            break;
    }
    return res;
}


bool valid_symbol(int c)
{
    return isalnum(c) || c == '_';
}


/**
 get_symbol() reads words:
 - starting with a alpha-character,
 - followed by alpha-numerical characters or '_'
 .
 */
std::string Tokenizer::get_symbol(std::istream & is, bool eat_line)
{
    int c = Tokenizer::get_space(is, eat_line);
    
    if ( !isalpha(c) )
        return "";
    
    std::string res = get_stuff(is, valid_symbol);
    
    VLOG("SYMBOL |" << res << "|\n");

    return res;
}


/**
 get_symbols() reads multiple symbols concatenated with ':'
 */
std::string Tokenizer::get_symbols(std::istream & is, bool eat_line)
{
    int c = Tokenizer::get_space(is, eat_line);
    
    if ( !isalpha(c) )
        return "";
    
    std::string res = get_symbol(is, eat_line);
    while ( is.peek() == ':' )
    {
        res += is.get();
        res += get_symbol(is, false);
    }
    
    return res;
}




bool valid_filename(int c)
{
    return isalnum(c) || c=='_' || c=='-' || c=='/' || c=='\\' || c=='.' || c==':' || c=='*';
}

/**
 Return next token if it looks likes a path to a file name
 */
std::string Tokenizer::get_filename(std::istream & is, bool eat_line)
{
    char c = Tokenizer::get_space(is, eat_line);
    
    if ( !isalpha(c) && c!='.' && c!='/' && c!='\\' && c!='*' )
        return "";
    
    std::string res = get_stuff(is, valid_filename);
    
    VLOG(" FILENAME |" << res << "|\n");
    
    return res;
}


/**
 Read multiple forms of integer:
 - integer-constant: digit-sequence [exponent-part]
 - digit-sequence:   [digit-sequence] digit
 - digit         :   one of [0123456789]
 - exponent-part: ('e' or 'E') [sign] digit-sequence
 */
std::string Tokenizer::get_integer(std::istream & is)
{
    bool accept_expon = true;
    
    char c = is.peek();
    if ( ! is.good() )
        return "";
    
    std::string res;
    
    if ( c == '+' || c == '-' )
        res.push_back(is.get());
    
    while ( !is.eof() )
    {
        c = is.peek();
        if ( isdigit(c) )
            res.push_back(is.get());
        else if (( c == 'e' || c == 'E' ) && accept_expon )
        {
            res.push_back(is.get());
            // accept an optional sign character
            c = is.peek();
            if ( c == '+' || c == '-' )
                res.push_back(is.get());
        }
        else
            break;
    }
    
    VLOG("INTEGER |" << res << "|\n");
    
    return res;
}

/**
 Read an integer, or return false if this was not possible.
 
 Note that the value of `what` will not change, if the input fails.
 */
bool Tokenizer::get_integer(std::istream & is, unsigned& var)
{
    std::streampos isp = is.tellg();
    unsigned backup = var;
    is >> var;
    if ( is.fail() )
    {
        is.clear();
        is.seekg(isp);
        var = backup;
        return false;
    }
    if ( isspace(is.peek()) )
        return true;
    // declare error, if there is no spacing character:
    is.seekg(isp);
    throw InvalidParameter("an integer is expected");
}

/**
 Read an integer, or return false if that was not possible.
 
 Note that the value of `what` will not change, if the input fails.
 */
bool Tokenizer::get_integer(std::istream & is, int& var)
{
    std::streampos isp = is.tellg();
    int backup = var;
    is >> var;
    if ( is.fail() )
    {
        is.clear();
        is.seekg(isp);
        var = backup;
        return false;
    }
    if ( isspace(is.peek()) )
        return true;
    // declare error, if there is no spacing character:
    is.seekg(isp);
    throw InvalidParameter("an integer is expected");
}


/**
 split the string into an integer and the remaining string.
 Any space after the integer is discarded.
 The argument `str` is truncated.
 */
bool Tokenizer::split_integer(std::string& str, unsigned& var)
{
    unsigned backup = var;
    std::istringstream iss(str);
    iss >> var;
    // to succeed, there should also be a space after the number:
    if ( !iss.fail() && iss.get() == ' ' )
    {
        // skip any additional space characters:
        get_space(iss, false);
        str = str.substr(iss.tellg());
        return true;
    }
    // restore value
    var = backup;
    return false;
}


/**
 Read a number specified in the standard (US) format with an optional '.':
 floating-point-constant: 
             [sign] fractional-constant [exponent-part]
             [sign] digit-sequence [exponent-part]
 fractional-constant: digit-sequence.digit-sequence
 digit-sequence:     [digit-sequence] digit
 digit         : one of [0123456789]
 exponent-part : ('e' or 'E') [sign] digit-sequence
 sign          : '+' or '–'
*/
std::string Tokenizer::get_real(std::istream & is)
{
    bool accept_point = true;
    bool accept_expon = true;

    char c = is.peek();
    if ( ! is.good() )
        return "";
    
    std::string res;
    
    if ( c == '+' || c == '-' )
        res.push_back(is.get());
    
    while ( !is.eof() )
    {
        c = is.peek();
        if ( isdigit(c) )
            res.push_back(is.get());
        else if ( c == '.' && accept_point )
        {
            res.push_back(is.get());
            accept_point = false;
        }
        else if (( c == 'e' || c == 'E' ) && accept_expon )
        {
            res.push_back(is.get());
            // only accept integer within exponent
            accept_point = false;
            accept_expon = false;
            // accept an optional sign character
            c = is.peek();
            if ( c == '+' || c == '-' )
                res.push_back(is.get());
        }
        else
            break;
    }
    
    VLOG("REAL |" << res << "|\n");

    return res;
}


bool valid_hexadecimal(int c)
{
    return isxdigit(c) || c=='x';
}


std::string Tokenizer::get_hexadecimal(std::istream & is)
{
    return get_stuff(is, valid_hexadecimal);
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 get_token() reads a block enclosed by '{}', '()' and '""',
 and returns it verbatim with the delimiting characters.
 */

std::string Tokenizer::get_token(std::istream & is, bool eat_line)
{
    Tokenizer::get_space(is, eat_line);
    int c = is.get();
    int d = is.peek();

    if ( c == EOF )
        return "";
    
    if ( d == EOF )
        return std::string(1,c);
   
    if ( Tokenizer::block_delimiter(c) )
        return get_block_content(is, c, block_delimiter(c));

    if ( isalpha(c) )
    {
        is.unget();
        return Tokenizer::get_symbol(is);
    }
    
    if ( c == '0' && d == 'x' )
    {
        is.unget();
        return Tokenizer::get_hexadecimal(is);
    }
    
    if ( isdigit(c) || (( c == '-' || c == '+' ) && isdigit(d)) )
    {
        is.unget();
        return Tokenizer::get_real(is);
    }
    
    VLOG(" ASCII |" << c << "|\n");

    // anything else is one character long:
    return std::string(1,c);
}

//------------------------------------------------------------------------------

char Tokenizer::block_delimiter(char c)
{
    switch(c)
    {
        case '(': return ')';
        case '{': return '}';
        case '[': return ']';
        case '"': return '"';
    }
    return 0;
}

/**
 This will read a block, assuming that opening delimiter has been read already.
 It will read characters until the given closing delimiter `c_out` is found.
 if `c_in` is not zero, the block is returned with `c_in` and `c_out` at the
 first and last positions.
 */
std::string Tokenizer::get_block_content(std::istream & is, char c_in, const char c_out)
{
    assert_true(c_out);
    std::string res;
    res.reserve(16384);
    char c = 0;
    
    if ( c_in )
        res.push_back(c_in);
    is.get(c);
    
    while ( is.good() )
    {
        if ( c == c_out )
        {
            if ( c_in )
                return res+c;
            else
                return res;
        }
        else if ( block_delimiter(c) )
            res.append( get_block_content(is, c, block_delimiter(c)) );
        else if ( c == ')' || c == '}' )
            throw InvalidSyntax("unmatched delimiter '"+std::string(1,c_in)+"'");
#if ( 0 )
        else if ( c == COMMENT_START )
        {
            // Read comments as lines, to inactivate symbols within: ')' and '}'
            std::string line;
            std::getline(is, line);
        }
#endif
        else
            res.push_back(c);
        is.get(c);
    }
    
    throw InvalidSyntax("missing '"+std::string(1,c_out)+"'");
    return "";
}


/**
 This will skip spaces and new-lines until a character is found.
 If this character is equal to `c_in`, then the block is read and returned.
 Otherwise returns empty string "".
 
 @returns content of the block without delimiters
 */
std::string Tokenizer::get_block(std::istream & is, char c_in)
{
    assert_true(c_in);
    
    int c = get_space(is, true);
    
    if ( c == c_in )
    {
        is.get();
        std::string res = get_block_content(is, 0, block_delimiter(c_in));
        VLOG("BLOCK |" << res << "|\n");
        return res;
    }

    return "";
}


std::string Tokenizer::get_block(std::istream & is)
{
    char c = Tokenizer::get_space(is, true);
    
    if ( block_delimiter(c) )
        return get_block_content(is, c, Tokenizer::block_delimiter(c));
    
    return "";
}


std::string Tokenizer::strip_block(std::string const& blok)
{
    int x = blok.size()-1;
    
    if ( x < 1 )
        return blok;
    
    char c = block_delimiter( blok[0] );
    if ( c )
    {
        if ( blok[x] != c )
            throw InvalidSyntax("mismatched enclosing symbols");
        return blok.substr(1, x-1);
    }
    return blok;
}


//------------------------------------------------------------------------------
#pragma mark -


std::string Tokenizer::get_until(std::istream& is, std::string what)
{
    std::string res;
    res.reserve(16384);
    int d = 0;
    char c = 0;
    is.get(c);
    
    while ( is.good() )
    {
        if ( c == what[d] )
        {
            ++d;
            if ( what[d] == '\0' )
                break;
        }
        else
        {
            if ( d == 0 )
            {
                res.push_back(c);
            }
            else
            {
                res.push_back(what[0]);
                if ( d > 1 ) {
                    is.seekg(-d, std::ios_base::cur);
                    d = 0;
                } else {
                    if ( c == what[0] )
                        d = 1;
                    else {
                        res.push_back(c);
                        d = 0;
                    }
                }
            }
        }
        is.get(c);
    }
    //std::clog << "get_until|" << res << std::endl;
    return res;
}


void Tokenizer::trim(std::string& str, const std::string& ws)
{
    std::string::size_type pos = str.find_last_not_of(ws);
    if ( pos != std::string::npos )
    {
        str.erase(pos+1);
        pos = str.find_first_not_of(ws);
        if ( pos != std::string::npos )
            str.erase(0, pos);
    }
    else str.clear();
}




// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef TOKENIZER_H
#define TOKENIZER_H

#include "assert_macro.h"
#include <iostream>
#include <string>

/// elementary tokenizer
/** A Tokenizer is used cut a character stream into words */
namespace Tokenizer
{
    /// return the corresponding closing delimiter character, or 0 if `entry_char` is not a known delimiter
    char block_delimiter(char entry_char);
 
    /// check if any character is not space
    bool not_space(std::string const&);
    
    /// skip space and new-line if `eat_line`==true, return the next character
    int get_space(std::istream & is, bool eat_line);
    
    /// skip upcomming characters for which isspace() is true, and new-line if `eat_line`==true
    int get_character(std::istream & is, bool eat_space=true, bool eat_line=false);
    
    /// read one signed integer, and throw exception if spacing character does not follow
    bool get_integer(std::istream&, int&);
    
    /// read one unsigned integer, and throw exception if spacing character does not follow
    bool get_integer(std::istream&, unsigned&);
    
    /// try to interpret `str` as `UNSIGNED_INT sub`. If successful, `str` is modified to be `sub`
    bool split_integer(std::string& str, unsigned&);

    
    /// read multiple forms of integer numbers
    std::string get_integer(std::istream&);

    /// return next token if it looks like a number
    std::string get_real(std::istream&);

    /// return next token if it looks like a hexadecimal number
    std::string get_hexadecimal(std::istream&);

    /// return next token if it looks like a variable name
    std::string get_symbol(std::istream & is, bool eat_line=false);
    
    /// return next token if it looks like a variable name
    std::string get_symbols(std::istream & is, bool eat_line=false);
    
    /// return next token if it looks like a file name
    std::string get_filename(std::istream & is, bool eat_line=false);
    
    /// return next token
    std::string get_token(std::istream & is, bool eat_line=false);
        
    /// accumulate characters until the next new-line character
    std::string get_line(std::istream & is);

    /// read text until `c_out` is encountered, assuming `c_in` was already read
    std::string get_block_content(std::istream & is, char c_in, char c_out);
    
    /// skip spaces and read a block delimited by `c_in`, or return empty string if `c_in` is not found
    std::string get_block(std::istream & is, char c_in);
    
    /// read a delimited set of characters, return block with delimiters included
    std::string get_block(std::istream & is);
    
    /// remove enclosing parenthesis at the start and at the end of `blok`
    std::string strip_block(std::string const& blok);

    /// read until `what` is found and stop immediately before (`what` is excluded from the returned string)
    std::string get_until(std::istream & is, std::string what);

    /// remove characters present in `ws` from the beggining and at the end of `str`
    void trim(std::string& str, const std::string& ws = " \t\n");
    
}

#endif



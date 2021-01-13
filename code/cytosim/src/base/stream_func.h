// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 30/7/09.


#include <iostream>

/// Simple operations on C++ streams
namespace StreamFunc
{
    
    /// remove non-conventional characters
    void clean_stream(std::ostream &, std::istream&);
    
    /// export lines of `val` that are not identical to `ref`
    void diff_stream(std::ostream &, std::istream& val, std::istream& ref);
    
    /// copy lines that do not start with character `skip`
    void skip_lines(std::ostream &, std::istream&, char skip);

    /// add `insert` before every line, but skip lines that start with `skip`
    void prefix_lines(std::ostream &, std::istream&, const char prefix[], char keep, char skip);

    /// print the line of `istream` indicating the position `pos`, with line number
    void show_line(std::ostream &, std::istream &, std::streampos pos, std::string const& prefix);
    
    /// same as `show_line()`, but output is returned as a string
    std::string get_line(std::istream &, std::streampos, std::string const& prefix);
    
    /// extract the lines located between `start` and `end`, with line numbers
    void show_lines(std::ostream &, std::istream &, std::streampos start, std::streampos end);
    
    /// same as show_lines(), but output is returned as a string
    std::string get_lines(std::istream &, std::streampos start, std::streampos end);
    
    /// return line number corresponding to position
    unsigned line_number(std::istream &, std::streampos);
    
    /// return line number corresponding to position
    unsigned line_number(std::istream &);

    /// replace in `src` all occurences of `fnd` by `rep`
    int  find_and_replace(std::string& src, std::string const& fnd, std::string const& rep);
}


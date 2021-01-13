// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 30/7/09.

#include "stream_func.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cctype>


void StreamFunc::clean_stream(std::ostream & os, std::istream& is)
{
    char c = 0;
    while ( is.good() )
    {
        is.get(c);
        
        // terminate the line for new-line and cariage-return
        if ( c == '\r' )
            os << '\n';
        // all type of spaces are substituted
        else if ( isspace(c) )
            os << ' ';
        // non=printable characters are removed
        else if ( isprint(c) )
            os << c;
        else
            std::cerr << "unprintable ascii "<< (int)c << " found\n";
    }
}


void StreamFunc::diff_stream(std::ostream & os, std::istream& val, std::istream& ref)
{
    std::string val_l, ref_l;
    val.seekg(0);
    ref.seekg(0);
    
    while ( val.good() )
    {
        std::getline(val, val_l);
        std::getline(ref, ref_l);
#if ( 0 )
        // print any line containing '{' or '}' 
        bool par = ( std::string::npos != ref_l.find_first_of("{}") );
        if ( val_l != ref_l || par )
#else
        if ( val_l != ref_l )
#endif
        {
            os << val_l << '\n';
            //os << val_l << "  (" << ref_l << ")\n";
        }
    }
}


void StreamFunc::skip_lines(std::ostream & os, std::istream& is, char skip)
{
    std::string line;
    
    while ( is.good() )
    {
        std::getline(is, line);
        if ( line.size() < 1 )
            continue;
        else if ( line[0] != skip )
            os << line << '\n';
    }
}



void StreamFunc::prefix_lines(std::ostream & os, std::istream& is, const char prefix[],
                              char keep, char skip)
{
    std::string line;
    
    while ( is.good() )
    {
        std::getline(is, line);
        if ( line.size() < 1 )
            continue;
        else if ( line[0] == keep )
            os << line << '\n';
        else if ( line[0] == skip )
            ;
        else
            os << prefix << line << '\n';
    }
}


/**
 The alignment of the vertical bar should match the one in PREF
 */
void print_line(std::ostream & os, int line_nb, std::string const& line)
{
    os << std::setw(7) << line_nb << " | " << line << '\n';
}

/**
 The alignment of the vertical bar should match the one in PREF
 */
void print_line(std::ostream & os, std::string const& prefix, std::string const& line)
{
    os << prefix << " " << line << '\n';
}

/**
 Print one line extracted from `is',
 and indicate the position `pos` with a arrowhead below this line.
 */
void StreamFunc::show_line(std::ostream & os, std::istream & is, std::streampos pos, std::string const& prefix)
{
    if ( !is.good() )
        is.clear();
    
    std::streampos isp, spos = is.tellg();
    if ( pos < 0 ) pos = spos;
    is.seekg(0);

    std::string line;
    
    while ( is.good()  &&  is.tellg() <= pos )
    {
        isp = is.tellg();
        std::getline(is, line);
    }
    
    print_line(os, prefix, line);
    is.clear();
    is.seekg(isp);
    
    line.clear();
    int c = 0;
    while ( is.tellg() < pos )
    {
        c = is.get();
        if ( c == EOF )
            break;
        if ( isspace(c) )
            line.push_back(c);
        else
            line.push_back(' ');
    }
    line.push_back('^');
    print_line(os, prefix, line);

    is.clear();
    is.seekg(spos);
}


std::string StreamFunc::get_line(std::istream & is, std::streampos pos, std::string const& prefix)
{
    std::ostringstream oss;
    show_line(oss, is, pos, prefix);
    return oss.str();
}


/**
 Output enough lines to cover the area specified by [start, end].
 Each line is printed in full and preceded with a line number
 */
void StreamFunc::show_lines(std::ostream & os, std::istream & is,
                            std::streampos start, std::streampos end)
{
    if ( !is.good() )
        is.clear();
    
    std::streampos isp = is.tellg();
    is.seekg(0);
    std::string line;
    
    unsigned int cnt = 0;
    while ( is.good()  &&  is.tellg() <= start  )
    {
        std::getline(is, line);
        ++cnt;
    }

    print_line(os, cnt, line);
    while ( is.good() &&  is.tellg() < end )
    {
        std::getline(is, line);
        ++cnt;
        print_line(os, cnt, line);
    }

    is.clear();
    is.seekg(isp);
}


std::string StreamFunc::get_lines(std::istream & is, std::streampos s, std::streampos e)
{
    std::ostringstream oss;
    show_lines(oss, is, s, e);
    return oss.str();
}


unsigned StreamFunc::line_number(std::istream & is, std::streampos pos)
{
    if ( !is.good() )
        is.clear();
    
    std::streampos spos = is.tellg();
    is.seekg(0);
    
    unsigned cnt = 0;
    std::string line;
    
    while ( is.good()  &&  is.tellg() <= pos )
    {
        std::getline(is, line);
        ++cnt;
    }
    
    is.clear();
    is.seekg(spos);
    return cnt;
}

unsigned StreamFunc::line_number(std::istream & is)
{
    if ( !is.good() )
        is.clear();
    
    return line_number(is, is.tellg());
}


int StreamFunc::find_and_replace(std::string & src,
                                 std::string const& fnd, std::string const& rep)
{
    int num = 0;
    std::string::size_type fLen = fnd.size();
    std::string::size_type rLen = rep.size();
    std::string::size_type pos = src.find(fnd, 0);
    while ( pos != std::string::npos )
    {
        num++;
        src.replace(pos, fLen, rep);
        pos += rLen;
        pos = src.find(fnd, pos);
    }
    return num;
}



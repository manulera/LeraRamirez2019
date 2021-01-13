// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef PARSER_H
#define PARSER_H

#include "interface.h"


/// Cytosim Parser to read and execute config files
/**
 This is where the syntax of the config file is defined
 */
class Parser : private Interface
{
private:
    
    /// disabled default constructor
    Parser();
    
    /// control switch to enable command 'set' (creating a property)
    bool      do_set;
    
    /// control switch to enable command 'change' (change a property)
    bool      do_change;
    
    /// control switch to enable command 'new' and 'delete' (create object)
    bool      do_new;
    
    /// control switch to enable command 'run' (run simulation)
    bool      do_run;
    
    /// control switch to enable command 'write' (write files)
    bool      do_write;
 
    /// position in current stream
    std::streampos spos;
    
public:
    
    /// set the permission of the parser
    Parser(Simul& s, bool allow_set, bool allow_change, bool allow_new, bool allow_run, bool allow_write);
    
    /// destructor
    virtual  ~Parser() {}
    
    //-------------------------------------------------------------------------------
    
    /// parse command \b set
    void      parse_set(std::istream&);
    
    /// parse command \b change
    void      parse_change(std::istream&);
    
    /// parse command \b new
    void      parse_new(std::istream&);
    
    /// parse command \b delete
    void      parse_delete(std::istream&);
    
    /// parse command \b mark
    void      parse_mark(std::istream&);

    /// parse command \b cut
    void      parse_cut(std::istream&);

    /// parse command \b run
    void      parse_run(std::istream&);
    
    /// parse command \b read
    void      parse_read(std::istream&);
    
    /// parse command \b write
    void      parse_report(std::istream&);
    
    /// parse command \b import
    void      parse_import(std::istream&);
    
    /// parse command \b export
    void      parse_export(std::istream&);
    
    /// parse command \b call
    void      parse_call(std::istream&);
    
    /// parse command \b repeat
    void      parse_repeat(std::istream&);

#ifdef NEW_COMMAND_FOR
    /// parse command \b for
    void      parse_for(std::istream&);
#endif
    
    /// parse command \b end
    void      parse_end(std::istream&);

    /// Parse a std::istream, `file_name` is used to build error messages
    void      parse(std::istream&, std::string const& file_name);

    //-------------------------------------------------------------------------------

    /// Open and parse the config file with the given name
    int       readConfig(std::string const& name);

    /// Parse the default config file i.e. calls readConfig(simul.prop->config)
    int       readConfig() { return readConfig(simul.prop->config); }

};

#endif


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef INTERFACE_H
#define INTERFACE_H


#include "simul.h"
#include <iostream>


/// Cytosim Application Programming Interface
/*
 A reduced set of commands to control and simulate
 a system of objects within cytosim.
 */
class Interface
{
private:
    
    /// disabled default constructor
    Interface();

protected:
    
    /// associated Simul
    Simul& simul;
    
public:
    
    /// associates with given Simul
    Interface(Simul& s);
    
    /// destructor
    virtual ~Interface() {}
    
    //-------------------------------------------------------------------------------
    
    /// this is called between commands during the execution process
    /**
     The overwritten version should call simul.relax() to make sure that
     the simulation data structures are coherent.
     It can perform additional things, for example display the simulation world
     */
    virtual void hold() {}
    
    /// Parse a text containing cytosim commands
    /**
     This is used to execute the 'event code' in parse_run(),
     The function must be defined in the derived class Parser
     */
    virtual void parse(std::istream&, std::string const& msg) = 0;
    
    //-------------------------------------------------------------------------------
    
    /// create a new Property of kind `k` from options set in Glossary
    Property*  execute_set(std::string const& kind, std::string const& name, Glossary&);

    /// change values in Property of kind `k` following options specified in Glossary
    Property*  execute_change(std::string const& name, Glossary&);
    
    /// change 'display' (and only this) in corresponding Property
    int        execute_change_all(std::string const& kind, Glossary&);

    /// read the specification of position and orientation of an object
    Isometry   read_position(Glossary&);
    
    /// return position and orientation of an object, with verification of 'placement'
    Isometry   get_placement(Glossary&, int);
    
    /// create 1 object of kind `k` with name `n`, following options in Glossary
    ObjectList execute_new(std::string const& name, Glossary&);
    
    /// create `cnt` objects of kind `k` with name `n`, randomly placed in space (no option)
    void       execute_new(std::string const& name, unsigned cnt);
    
    /// delete `cnt` objects of kind `k` with name `n`, following options in Glossary
    void       execute_delete(std::string const& name, Glossary&);
    
    /// mark `cnt` objects of kind `k` with name `n`, following options in Glossary
    void       execute_mark(std::string const& name, Glossary&);

    /// cut fibers, following different options in Glossary
    void       execute_cut(std::string const& name, Glossary&);
    
    /// import objects from another file
    void       execute_import(std::string const& file, std::string const& what, Glossary&);

    /// export objects from another file
    void       execute_export(std::string& file, std::string const& what, Glossary&);

    /// write output file with object coordinates or information on objects
    void       execute_report(std::string& file, std::string const& what, Glossary&);
    
    /// perform simulation steps
    void       execute_run(bool do_write, Glossary&);
    
    /// execute miscellaneous functions
    void       execute_call(std::string& func, Glossary&);

};

#endif


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FAKE_H
#define FAKE_H

#include "object.h"
#include "organizer.h"
#include "fake_prop.h"
#include "point_exact.h"
#include "solid.h"

//------------------------------------------------------------------------------

///a set of two asters held together by a Solid 
/**
 This object cannot handle the destruction of the Asters
 */

class Fake : public Organizer
{

private:

    /// Property
    FakeProp const* prop;
    
    /// Linking Solid
    Solid *      fkSolid;
    
    /// First Aster
    Aster *      fkAster1;
    
    /// Second Aster
    Aster *      fkAster2;
    
    /// connections
    std::vector<PointExact> asterPoints, solidPoints;
    
    /// find aster
    Aster *    findAster(std::string const&, Simul&);

public:
    
    /// constructor
    Fake(FakeProp const* p) : prop(p) { fkSolid = 0; fkAster1 = 0; fkAster2 = 0; }
 
    /// construct all the dependent Objects of the Organizer
    ObjectList build(Glossary&, Simul&);
 
    /// destructor  
    virtual   ~Fake() { prop = 0; }

    /// perform one Monte-Carlo step
    void       step();
    
    /// add interactions to the Meca
    void       setInteractions(Meca &) const;
    
    /// return pointer to central Solid
    Solid *    solid() const { return fkSolid; }

 
    /// move all associated objects
    void       translate(Vector const& T);
    
    /// rotate all associated objects
    void       rotate(Rotation const& T);

    //------------------------------ read/write --------------------------------
    
    /// get access to display links
    unsigned   getLink(unsigned, Vector&, Vector&) const;

    /// a unique character identifying the class
    static const Tag TAG = 'k';
    
    /// return unique character identifying the class
    Tag        tag() const { return TAG; }
    
    /// return Object Property
    Property const* property() const { return prop; }
    
    /// return PointDisp of Solid
    PointDisp *  disp() const;

 };


#endif


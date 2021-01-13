// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ORGANIZER_H
#define ORGANIZER_H

#include "assert_macro.h"
#include "glossary.h"
#include "object.h"
#include "buddy.h"

class Meca;
class Simul;
class Mecable;
class PointDisp;
class Display;

/// an assemblage of Object
/** 
Organizer contains an Array of pointers of type Mecable*.
These Mecables are organized by Organizer::setInteraction()
which is implemented in the derived classes, eg. Bundle, Aster & Nucleus.
*/
class Organizer: public Object, private Buddy
{

private:

    typedef std::vector<Mecable*> MecableList;
    
    /// list of Objects that are `organized`
    MecableList   mObjects;
    
public:

    /// default constructor
    Organizer() { }
    
    /// destructor
    virtual      ~Organizer();
    
    /// construct all the dependent Objects of the Organizer
    virtual ObjectList    build(Glossary&, Simul&) = 0;

    //--------------------------------------------------------------------------

    /// number of objects currently organized
    unsigned int          nbOrganized() const  { return mObjects.size(); }
    
    /// return Mecable at index `n`
    Mecable *             organized(unsigned n) const { assert_true(n<mObjects.size()); return mObjects[n]; }
    
    /// add Mecable at end of list
    void                  grasp(Mecable *);

    /// add Mecable at index `n`
    void                  grasp(Mecable *, unsigned);

    /// handles the disapearance of one of the organized object
    void                  goodbye(Buddy *);
    
    /// add objects to Simul if they are not already linked
    virtual void          addOrganized(Simul&);
    
    //--------------------------------------------------------------------------

    /// organizers cannot be moved
    bool                  mobile() const { return false; }
    
    /// return the center of gravity
    virtual Vector        position() const;

    /// return the average of all model points
    virtual Vector        positionP(unsigned) const;
    
    /// move all associated objects
    void                  translate(Vector const& T);
    
    /// rotate all associated objects
    void                  rotate(Rotation const& T);

    /// monte-carlo simulation step
    virtual void          step() {}
    
    /// add interactions to the Meca
    virtual void          setInteractions(Meca &) const {}
    
    /// sum the drag coefficient of all objects
    real                  dragCoefficient() const;
    
    
    /// obtain the ends of the link number `inx`, or returns zero if the link does not exist
    virtual unsigned      getLink(unsigned inx, Vector&, Vector&) const { return 0; }
    
    /// display parameters 
    virtual PointDisp *   disp() const { return 0; }
    
    //--------------------------------------------------------------------------
    
    /// a static_cast<> of Node::next()
    Organizer *   next()  const  { return static_cast<Organizer*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Organizer *   prev()  const  { return static_cast<Organizer*>(nPrev); }
    
    //--------------------------------------------------------------------------

    /// read
    void          read(Inputter&, Simul&, Tag);
    
    /// write
    void          write(Outputter&) const;
};



#endif

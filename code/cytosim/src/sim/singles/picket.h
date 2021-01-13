// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef PICKET_H
#define PICKET_H

#include "single.h"
class FiberGrid;

/// a Single attached at a fixed position.
/**
 This Single is fixed at its foot position in absolute space.
 A link is created if the Hand is attached.

 @ingroup SingleGroup
 */
class Picket : public Single
{
public:
    
    /// sPos should never change
    void    beforeDetachment(Hand const*);
       
public:

    /// constructor
    Picket(SingleProp const*, Vector const& = Vector(0,0,0));

    /// destructor
    ~Picket();
    
    ///return the position in space of the object
    Vector  position()           const  { return sPos; }

    /// true if object accepts translations
    bool    mobile()             const  { return true; }
    
    /// translate object's position by the given vector
    void    translate(Vector const& w)  { sPos += w; }
    
    /// true if Single creates an interaction
    bool    hasForce() const { return true; }

    /// tension in the link = stiffness * ( posFoot() - posHand() )
    Vector  force() const;

    /// Monte-Carlo step for a free Single
    void    stepF(const FiberGrid&);
    
    /// Monte-Carlo step for a bound Single
    void    stepA();
    
    /// add interactions to the Meca
    void    setInteractions(Meca &) const;
    
};


#endif

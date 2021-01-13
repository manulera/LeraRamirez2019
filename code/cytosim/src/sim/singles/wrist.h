// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef WRIST_H
#define WRIST_H

#include "single.h"
#include "mecable.h"


/// a Single anchored to a Mecable.
/**
 The Wrist is anchored to Mecable `mBase`:
 - between the DIM+1 consecutive points from 'mPoint' to 'mPoint+3'
 - using interpolation coefficients 'mCoef' on these points.
 .

 @ingroup SingleGroup
 */
class Wrist : public Single
{
protected:
    
    /// Mecable to which Wrist is anchored
    Mecable const* mBase;
    
    // number of point that are being interpolated
    unsigned       mOrder;
    
    /// Index of the point of Mecable to which Wrist is anchored
    unsigned       mPoint;
    
    /// Coefficient on each point
    real           mCoef[4];
    
public:
     
    /// constructor
    Wrist(SingleProp const*, Mecable const*, unsigned point);
    
    /// constructor on two points
    Wrist(SingleProp const*, Mecable const*, unsigned ref, Vector pos);

    /// destructor
    ~Wrist();
    
    //--------------------------------------------------------------------------
    
    /// return the position in space of the object
    Vector  position() const { return posFoot(); }
    
    /// true if object accepts translations
    bool    mobile() const { return false; }
    
    /// translate object's position by the given vector
    void    translate(Vector const& T) { }
    
    /// modulo the position of the grafted
    void    foldPosition(const Modulo * s) { }

    //--------------------------------------------------------------------------
    
    /// Object to which this is attached
    Mecable const* base() const { return mBase; }

    /// the position of what is holding the Hand
    Vector  posFoot() const;
    
    
    /// true if Single creates a link
    bool    hasForce() const { return true; }

    /// force = stiffness * ( posFoot() - posHand() )
    Vector  force() const;
    
    
    /// Monte-Carlo step for a free Single
    void    stepF(const FiberGrid&);
    
    /// Monte-Carlo step for a bound Single
    void    stepA();

    /// add interactions to the Meca
    void    setInteractions(Meca &) const;

    //--------------------------------------------------------------------------
    
    /// the Wrist uses a specific TAG to distinguish itself from the Single
    static const Tag TAG = 'w';
    
    /// return unique character identifying the class
    Tag     tag() const { return TAG; }
    
    /// read from file
    void    read(Inputter&, Simul&, Tag);
    
    /// write to file
    void    write(Outputter&) const;
    
};


#endif

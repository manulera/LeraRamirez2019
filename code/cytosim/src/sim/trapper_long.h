// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef TRAPPER_LONG_H
#define TRAPPER_LONG_H

#include "trapper.h"

class TrapperLong : public Trapper
{
    /// the side (top/bottom) of the interaction
    mutable Torque mArm;
    
    /// used to calculate `mArm`
    static Torque calcArm(const PointInterpolated & pt, Vector const& pos, real len);
    
public:
    
    /// create following the specifications in the CoupleProp
    TrapperLong(TrapperProp const*, Vector const & w = Vector(0,0,0));
    
    /// destructor
    virtual ~TrapperLong();
    
    /// position on the side of fiber1 used for sideInteractions
    Vector  posSide() const;
    
    /// force between hands, essentially: stiffness * ( cHand2->posHand() - cHand1->posHand() )
    Vector  force() const;
    
    /// add interactions to the Meca
    void    setInteractions(Meca &) const;
    
    /// add interactions to the Meca
    void    setInteractionsAF(Meca &) const;
    
    /// add interactions to the Meca
    void    setInteractionsFA(Meca &) const;
    
};


#endif

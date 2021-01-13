// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef trapper_prop_H
#define trapper_prop_H

#include "couple_prop.h"
#include "trapper.h"

/// Additional Property for Trapper
/**
 @ingroup Properties
 */
class TrapperProp : public CoupleProp
{
    
    friend class Trapper;
    friend class TrapperLong;
    
public:
    
    /// The target property
    std::string target_prop_name;
    HandProp * target_prop;
    /// The lattice on which the target single moves got from the properties
    unsigned int target_lat;
    
    /// Rate of trapping singles when they are attached to the microtubule
    real trapping_rate;
    
    /// Rate of untrapping
    real untrapping_rate;
    
    /// Force of the bond between couple and single
    real untrapping_force;
    
    /// Integer indicating whether the single is attached to the hand on the same filament (0), on the other filament (1), or in the center (2)
    bool trap_other_filament;
    
    /// The trapping rate from solution to any trapper (regardless of its binding stage)
    real trapping_rate_solution;
private:
    
    real trapping_rate_dt;
    
    real untrapping_rate_dt;
    
    real inv_untrapping_force;

    /// In this version of the trapper, the trapper stiffness corresponds to the weight of the linker for the TriLink when all three heads are connected. In the current situation, where we connect the three hands to a central dragless point with linkers of equal stiffness. If only two hands are connected (like in a normal couple), stiffness should be the one given by CoupleProp->stiffness, therefore, the stiffness of the linker when the three elements are connected is CoupleProp->stiffness*2/3
    real trap_stiffness;
    
    real trapping_rate_solution_dt;
public:
    
    /// constructor
    TrapperProp(const std::string& n) : CoupleProp(n)  { clear(); }
    
    /// destructor
    ~TrapperProp() { }
    
    /// return a Trapper with this property
    Couple * newCouple(Glossary*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const*);
    
    /// return a carbon copy of object
    Property* clone() const { return new TrapperProp(*this); }
    
    /// write all values
    void write_values(std::ostream &) const;
    
    
};

#endif


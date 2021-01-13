// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef swapper_prop_H
#define swapper_prop_H

#include "couple_prop.h"
#include "swapper.h"

/// Additional Property for Swapper
/**
 @ingroup Properties
 */
class SwapperProp : public CoupleProp
{
    
    friend class Swapper;
    
public:
    
    real swap_rate;
    
private:
    
    real swap_rate_dt;
    
public:
    
    /// constructor
    SwapperProp(const std::string& n) : CoupleProp(n)  { clear(); }
    
    /// destructor
    ~SwapperProp() { }
    
    /// return a Swapper with this property
    Couple * newCouple(Glossary*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const*);
    
    /// return a carbon copy of object
    Property* clone() const { return new SwapperProp(*this); }
    
    /// write all values
    //    void write_values(std::ostream &) const;
    
    
};

#endif


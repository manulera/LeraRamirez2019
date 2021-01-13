// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef flipper_prop_H
#define flipper_prop_H

#include "couple_prop.h"
#include "flipper.h"

/// Additional Property for Flipper
/**
 @ingroup Properties
 */
class FlipperProp : public CoupleProp
{
    
    friend class Flipper;
    
public:
    
    /// Chance of flipping used by the flipper couple subclass
    real flip_rate1_for;
    real flip_rate1_bak;
    real flip_rate2_for;
    real flip_rate2_bak;
    
//    HandProp * hand1_base_prop;
    HandProp * hand1_alter_prop;

//    HandProp * hand2_base_prop;
    HandProp * hand2_alter_prop;
    
//    std::string  hand1_base;
    std::string  hand1_alter;
//    std::string  hand2_base;
    std::string  hand2_alter;
    
    bool hand1_flips;
    bool hand2_flips;
    
private:

    real flip_rate1_for_dt;
    real flip_rate1_bak_dt;
    real flip_rate2_for_dt;
    real flip_rate2_bak_dt;

public:
    
    /// constructor
    FlipperProp(const std::string& n) : CoupleProp(n)  { clear(); }
    
    /// destructor
    ~FlipperProp() { }
    
    /// return a Flipper with this property
    Couple * newCouple(Glossary*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const*);
    
    /// return a carbon copy of object
    Property* clone() const { return new FlipperProp(*this); }
    
    /// write all values
//    void write_values(std::ostream &) const;
    
    ///\todo: propper write_values is missing for this.
};

#endif


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_PROP_H
#define SPACE_PROP_H

#include "real.h"
#include "property.h"
#include "vector.h"

class PointDisp;
class Glossary;
class Space;
class Simul;


#define NEW_DYNAMIC_SPACES


/// Property for Space
/**
 @ingroup Properties
 */
class SpaceProp : public Property
{
    friend class Space;
    
    /// create a new Space
    static Space * newSpace(SpaceProp const*, Glossary&);

public:

    /**
     @defgroup SpacePar Parameters of Space
     @ingroup Parameters
     @{
    */
    
    /// shape followed by parameters (see @ref SpaceGroup)
    std::string  geometry;
    
    /// primitive (rectangle, etc), derived from geometry
    std::string  shape;
    
    /// sizes of the space, or name of input file for polygon
    std::string  dimensions;

    /// display string (see @ref PointDispPar)
    std::string  display;
	
#ifdef NEW_DYNAMIC_SPACES
    /// Viscosity
    real		 viscosity;

    /// Viscosity for rotation
    real		 viscosity_rot;

    /// Set volume
    real         volume;
    
	/// Surface tension
	real		 tension;
#endif
    
    /// @}
    
    /// equal to time_step / viscosity
    real        mobility_dt, mobility_rot_dt;
    
    /// derived variable: a copy of `dimensions`
    std::string  dimensions_old;
    
    /// derived variable: file name including a full path (the first word in `geometry`)
    std::string  shape_spec;
    
    /// derived variable: flag to indicate that `display` has a new value
    bool         display_fresh;

    /// derived variable: parameters extracted from `display`
    PointDisp *  disp;

public:

    /// constructor
    SpaceProp(const std::string& n) : Property(n), disp(0) { clear(); }

    /// destructor
    ~SpaceProp() { }
 
    /// create a new Space
    Space * newSpace(Glossary&) const;
    
    /// identifies the property
    std::string category() const { return "space"; }
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// check and derive parameters
    void complete(Simul const*);
        
    
    /// return a carbon copy of object
    Property* clone() const { return new SpaceProp(*this); }
    
    /// write all values
    void write_values(std::ostream &) const;
    
};

#endif


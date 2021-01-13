// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef HAND_PROP
#define HAND_PROP

#include "real.h"
#include "common.h"
#include "property.h"


class Glossary;
class PointDisp;
class Hand;
class HandMonitor;

/// Property for Hand
/**
 @ingroup Properties
*/
class HandProp : public Property
{
    friend class Hand;
    
public:
      
    /// return one of the Property derived from HandProp
    static HandProp * newProperty(const std::string& n, Glossary&);
    
public:
    
    /**
     @defgroup HandPar Parameters of Hand
     @ingroup Parameters
     @{
     */
    
    /// rate of attachment when the Hand is within `binding_range` (also known as `binding[0]`)
    /**
     According to Leduc et al. PNAS 2004 vol. 101 no. 49 17096-17101
     the molecular binding_rate of kinesin is 4.7 +/- 2.4 /s.
     <em>
     http://dx.doi.org/10.1073/pnas.0406598101 \n
     http://www.pnas.org/content/101/49/17096.abstract
     </em>
     */
    real         binding_rate;
    
    
    /// maximum distance at which the Hand can bind (also known as `binding[1]`)
    real         binding_range;
    
    
    /// can be set to restrict binding to certain type of Fiber
    /**
     The binding to a fiber is allowed only if the keys of the Hand and Fiber match.
     The test is a BITWISE-AND of the two keys:
     @code
     if ( fiber:binding_key & hand:binding_key )
        allowed = true;
     else
        allowed = false;
     @endcode
     */
    unsigned int binding_key;
    
    
    /// detachment rate in the absence of load (also known as `unbinding[0]`)
    /**
     Kramers theory specifies that the detachment rate depends on the force
     in the link:
     @code
     off_rate = RATE * exp( force / FORCE )
     @endcode
     RATE is specified as `unbinding_rate`, and FORCE as `unbinding_force`,
     but one can also directly specify `unbinding = RATE, FORCE`.
     
     Recent measurements:
     <em>
     <b>Examining kinesin processivity within a general gating framework</b>
     Andreasson et al. eLife 2015;4:e07403
     http://dx.doi.org/10.7554/eLife.07403
     </em>
     provide for a conventional kinesin:
     @code
     unbinding_rate = 1 / s
     unbinding_force = 5 pN
     @endcode
     
     (see @ref Stochastic)
     */
    real         unbinding_rate;
    
    
    /// characteristic force of unbinding (also known as `unbinding[1]`)
    /**
     @copydetails unbinding_rate
     */
    real         unbinding_force;
    
    
    /// if true, the Hand can also bind directly to the tip of fibers
    /**
     This parameter determines the binding ability of a Hand that is located within
     `binding_range` of a fiber, but at a position where the orthogonal projection is
     outside the Fiber, either below the Minus end or above the Plus end.
     In this case, the attachement will be at the ends of the fiber, and it is allowed
     only if `bind_also_ends` is true. In other words, with 'bind_also_ends==true', the
     capture regions of the fibers are extended by adding two hemi-spheres at the two
     ends of each fibers, with a radius `binding_range`.
     
     <em>default value = false</em>
     */
    bool         bind_also_ends;
	
	/// if true, the Hand can bind only near the ends of the fibers
	/**
	 This determines that a Hand can only bind near the ends of the fiber.
	 This parameter can be 'none', 'plus_end', 'minus_end' or 'both_ends'.
     Binding is allowed on positions located within a distance 'bind_end_range'
     from the specified end ('bind_end_range' is specified as `bind_only_end[1]`).
	 
	 <em>default value = none</em>
	 */
	FiberEnd     bind_only_end;
	
    
    /// cutoff associated with `bind_only_end` where hand may bind (set as `bind_only_end[1]`)
    real         bind_end_range;

	
	/// if true, only bind fiber tip if no other hand is bound already
	bool         bind_only_free_end;

	
    /// if false, the Hand will detach immediately upon reaching a growing or a static fiber end
    /**
     A Hand may reach the tip of the fiber on which it is bound,
     either by self-movement, or possibly dragged by some other force.
     When this happens, `hold_growing_end` will determine if the Hand
     will detach or not.
     
     <em>default = false</em>
     */
    bool         hold_growing_end;
    
    
    /// if false, the Hand will detach immediately upon reaching a shrinking fiber end
    /**
     A Hand may reach the tip of the fiber on which it is bound,
     of the tip of the fiber may reach a immobile hand because it is disassembling.
     When this happens, `hold_shrinking_end` will determine if the Hand
     will detach or not.
     If `hold_shrinking_end` is true, the hand will be relocated to track the end.

     <em>default = false</em>
     */
    bool         hold_shrinking_end;
    
    
    /// specialization
    /**
     @copydetails HandGroup
     */
    std::string  activity;
    
    
    /// display parameters (see @ref PointDispPar)
    std::string  display;
    
    /** @} */
    
    /// derived variable: inverse of unbinding_force. This is a flag to Kramer
    real unbinding_force_inv;
    
    /// oversampled binding_rate_dt
    real   binding_rate_dt_4;

public:
    
    /// binding_rate_dt = binding_rate * time_step;
    real   binding_rate_dt;
    
    /// binding_range_sqr = binding_range * binding_range;
    real   binding_range_sqr;
    
    /// unbinding_rate_dt = unbinding_rate * time_step;
    real   unbinding_rate_dt;
    
    /// flag to indicate that `display` has a new value
    bool   display_fresh;

    /// the display parameters for this category of Hand
    PointDisp * disp;
    
public:
    
    /// constructor
    HandProp(const std::string& n) : Property(n), disp(0) { clear(); }
    
    /// destructor
    ~HandProp() { }
    
    /// return a Hand with this property
    virtual Hand * newHand(HandMonitor* h) const;

    /// identifies the property
    std::string category() const { return "hand"; }
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    virtual void read(Glossary&);
    
    /// compute values derived from the parameters
    virtual void complete(Simul const*);
    
    /// perform additional tests for the validity of parameters, given the elasticity
    virtual void checkStiffness(real stiff, real len, real mul, real kT) const;
    
    /// Estimate attachment propensity per unit length of fiber
    real  bindingSection(bool) const;
    
    /// write all values
    void write_values(std::ostream &) const;
    
    /// return a carbon copy of object
    Property* clone() const { return new HandProp(*this); }
    
    //Added by manu
    /// logical indicating if the hand should move with gillespie
    bool is_gillespie;
    
    /// 0 if no lattice, lat_val otherwise
    unsigned int lat_val;
    
    /// A value to increase the probability to bind in a position with two MTs. To simulate proteins that are recruited at the overlap by other proteins.
    real overlap_affinity;
    
    /// Complete the property of the hand with some value from the couple, in this case the stiffness
    virtual void complete_from_couple(Simul const* sim, real stiffness){};
    
};


#endif


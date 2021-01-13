// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef CLASSIC_FIBER_PROP
#define CLASSIC_FIBER_PROP

#include "fiber_prop.h"


class Glossary;


/// Enables support for an option to make catastrophe rate dependent on fiber length
//#define NEW_LENGTH_DEPENDENT_CATASTROPHE

/// Enables support for an option that induces catastrophe if the PLUS_END is outside
//#define NEW_CATASTROPHE_OUTSIDE


/// additional Property for ClassicFiber
/**
 @ingroup Properties
 */
class ClassicFiberProp : public FiberProp
{
    friend class ClassicFiber;
    
public:
    
    /**
     @defgroup ClassicFiberPar Parameters of ClassicFiber
     @ingroup Parameters
     Inherits @ref FiberPar.
     The first column of numbers applies to PLUS_END, and the second to MINUS_END
     @{
     */
    
    /// see @ref ClassicFiber
    
    /// Speed of assembly state
    /**
     Antagonistic force decrease assembly rate exponentially if it is directed against the assembly:
     @code
     if ( force < 0 )
         speed = growing_speed * free_polymer * exp( force / growing_force ) + growing_off_speed;
     else
         speed = growing_speed * free_polymer + growing_off_speed;
     @endcode
     
     The parameters are:
     - `growing_speed`, the force-dependent and concentration-dependent assembly rate.
     - `growing_off_speed`, a constant term, normally negative to represent spontaneous disassembly.
     - `growing_force`, the characteristic force
     .
     In this equation, `free_polymer` represents the fraction of free monomers in [0,1].
     Antagonistic force is negative ( force < 0 ) if it is directed against fiber assembly.
     */
    real    growing_speed[2];

    /// Constant term in the growing speed equation
    real    growing_off_speed[2];

    
    /// Characteristic force of assembly state (default=+inf)
    /**
     Antagonistic force decrease assembly rate exponentially.
     */
    real    growing_force[2];
    
    
    /// speed of disassembly state
    /**
     Disassembly occurs always at the specified speed:
     @code
         speed = shrinking_speed;
     @endcode
     */
    real    shrinking_speed[2];
    
    
    /// Rate of stochastic switching from assembly to disassembly
    /**
     The catastrophe rate depends on the growth rate of the corresponding tip,
     which is itself reduced by antagonistic force:
     @code
     catastrophe_rate_real = catastrophe_rate_stalled / ( 1 + coef * growing_speed_real )
     @endcode
     where `growth_speed_real` is calculated as explained in @ref growing_speed,
     and `coef` is set to match the given `catastrophe_rate` in the absence of force:

     @code
     coef = ( catastrophe_rate_stalled/catastrophe_rate - 1.0 ) / growing_speed_unloaded
     growing_speed_unloaded = growing_speed + growing_off_speed;
     @endcode
     
     Note that if `catastrophe_rate_stalled >> catastrophe_rate`, the equation simplies to
     @code
     catastrophe_rate_real = catastrophe_rate * growing_speed_unloaded / growing_speed_real
     @endcode
     */
    real    catastrophe_rate[2];

    /// Rate of catastrophe when the growth is stalled
    /**
     If this parameter is not set, the catastrophe rate will not depend on growth speed.
     */
    real    catastrophe_rate_stalled[2];

#ifdef NEW_CATASTROPHE_OUTSIDE
    
    /// Flag to trigger immediate catastrophe if the PLUS_END is outside
    bool    catastrophe_outside;

#endif
    
#ifdef NEW_LENGTH_DEPENDENT_CATASTROPHE
    
    /// Switch to enable the length-dependent catastrophe rate
    /**
     If this is defined, the catastrophe rate will depend on the length of the fiber:
     @code
     catastrophe_rate_real = catastrophe_rate * length() / catastrophe_length;
     @endcode
     */
    real    catastrophe_length;
    
#endif
    
    /// Rate of stochastic switching from disassembly to assembly
    real    rescue_rate[2];
    
    /// switching rate to the growing state for a fiber shorter than `min_length` (default=0)
    real    rebirth_rate[2];

    /// if `false`, the fiber will be destroyed if it is shorter than `min_length` (default=`false`)
    bool    persistent;
    
    /// @}
    
private:
    
    real    shrinking_speed_dt[2];
    real    growing_speed_dt[2];
    real    growing_off_speed_dt[2];
    real    catastrophe_rate_dt[2];
    real    catastrophe_rate_stalled_dt[2];
    real    catastrophe_coef[2];
    real    rescue_prob[2], rebirth_prob[2];

public:
    
    /// constructor
    ClassicFiberProp(const std::string& n) : FiberProp(n) { clear(); }

    /// destructor
    ~ClassicFiberProp() { }
    
    /// return a Fiber with this property
    Fiber* newFiber() const;
    
    /// set default values
    void clear();
       
    /// set using a Glossary
    void read(Glossary&);
   
    /// check and derive parameter values
    void complete(Simul const*);
    
    /// return a carbon copy of object
    Property* clone() const { return new ClassicFiberProp(*this); }

    /// write
    void write_values(std::ostream &) const;

};

#endif


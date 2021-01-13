// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FIBER_PROP
#define FIBER_PROP

#include "real.h"
#include "property.h"
#include "common.h"
#include "field.h"
#include "sim.h"

class Fiber;
class Glossary;
class FiberDisp;
class SingleProp;
class SingleSet;
class Space;


/// Property for a Fiber
/**
 @ingroup Properties
 */
class FiberProp : public Property
{
    friend class Fiber;
    friend class FiberSet;

public:
    
    /**
     @defgroup FiberPar Parameters of Fiber
     @ingroup Parameters
     These are the parameters for Fiber
     @{
     */
    
    
    /// elastic modulus for bending elasticity
    /**
     The bending elasticity modulus `rigidity` has units of pN.um^2, 
     and is related to the persitence length `Lp` via the Boltzman constant and
     absolute temperature (kT = k_B * T):
     @code
     rigidity = Lp * kT
     @endcode
     
     Many measurments have been made and they agree somewhat:\n
     
     Filament                         |    Lp        | rigidity
     ---------------------------------|--------------|-----------------
     Microtubule                      |   ~ 7200 um  | ~ 30 pN.um^2
     Stabilized Microtubule           |   ~ 2200 um  | ~ 10 pN.um^2
     F-actin                          | ~  9--10 um  | ~ 0.04 pN.um^2
     Phalloidin-stabilized F-actin    | ~ 17--18 um  | ~ 0.08 pN.um^2
     
     <em>
     Flexural rigidity of microtubules and actin filaments measured from thermal fluctuations in shape.\n
     Gittes et al.\n JCB vol. 120 no. 4 923-934 (1993)\n
     http://dx.doi.org/10.1083/jcb.120.4.923 \n
     http://jcb.rupress.org/content/120/4/923
     </em>
     
     <em>
     Flexibility of actin filaments derived from thermal fluctuations.\n
     Isambert, H. et al.\n  J Biol Chem 270, 11437–11444 (1995)\n
     http://www.jbc.org/content/270/19/11437
     </em>
     */
    real         rigidity;
    
    
    /// desired distance between model points
    /**
     `segmentation` is a distance, which affects the precision by which the
     shape of a filament is simulated. Specificially, the number of segments 
     used for a filament of length `L` is the integer `N` that minimizes:
     @code
     fabs( L / N - segmentation )
     @endcode
     
     As a rule of thumb, segmentation should scale with rigidity, depending on 
     the expected magnitude of the forces experienced by the filament:
     @code
     segmentation = sqrt(rigidity/force)
     force ~ rigidity / segmentation^2
     @endcode
     
     Generally, a simulation should not be trusted if any filament contains kinks
     (i.e. if the angle between consecutive segments is greater then 45 degrees).
     In that case, the simulation should be redone with a segmentation divided by 2,
     and the segmentation should be reduced until kinks do not appear.
     */
    real         segmentation;
    
    /// length or initial-length for dynamic fibers
    real         length;
    
    /// Minimum length (this limits the length in some cases)
    real         min_length;
    
    /// Maximum length (this limits the length in some cases)
    real         max_length;

    /// amount of monomer available to make this type of fiber
    /**
     If set, this parameter will limit the total length of all the Fibers of this
     class, by making the assembly rate of the fibers dependent on the amount of
     unused material (ie. 'monomers'):
     @code
     assembly_speed = ( 1 - sum_of_all_fiber_length / total_polymer ) * [...]
     @endcode
     This links the assembly for all the fibers within one class.
     Thus assembly speed decreases linearly with the total amount of polymer in the cell,
     i.e. proportional to the normalized concentration of 'monomers'.
     
     By default `total_polymer = infinite`, and the assembly rate is not reduced.
     */
    real         total_polymer;

    /// effective viscosity (if unspecified, simul:viscosity is used)
    /**
     Set the effective `viscosity` to lower or increase the drag coefficient of a particular class of fibers. This makes it possible for example to reduce the total drag coefficient of an aster.
     If unspecified, the global `simul:viscosity` is used.
     */
    real         viscosity;
    
    /// radius used to calculate mobility
    /**
     These length are used in the formula for the mobility of a cylinder:
     - hydrodynamic_radius[0] corresponds to the radius of the fiber
     - hydrodynamic_radius[1] is a cut-off for the length of the fiber
     .
     */
    real         hydrodynamic_radius[2];

    /// if true, the mobility of a cylinder moving near a plane will be used
    /**
     You can select between two possible formulas to calculate viscous drag coefficient:
     @code
     if ( fiber:surface_effect )
         drag = dragCoefficientSurface();
     else
         drag = dragCoefficientVolume();
     @endcode
     <hr>
     @copydetails Fiber::dragCoefficientVolume
     <hr>
     @copydetails Fiber::dragCoefficientSurface
     */
    bool         surface_effect;
    
    /// distance of fluid between slide and cylinder surface (set as `surface_effect[1]`)
    real         cylinder_height;

    
    /// can be set to control which Hands may bind
    /**
     To decide if a Hand may bind to a Fiber, the two binding_keys are compared:
     @code
     allowed = ( fiber:binding_key & hand:binding_key )
     @endcode
     Attachement is forbiden if the bitwise AND returns false, which is true if the two binding_key do not share any common digit in base 2. For most usage, you would thus use powers of 2 to distinguish fibers:
     - microtubule: binding_key = 1,
     - actin: binding_key = 2,
     - etc.
     .
     More complex combinations can be created by using all the bits of binding_key.
     */
    unsigned int binding_key;

    /// if true, a Lattice is associated to this fiber
    int          lattice;
    
    /// unit length associated with Lattice
    real         lattice_unit;
    
    /// if true, the quantities in the lattice can cut the fiber
    int          lattice_cut_fiber;

    /// flux speed of substance on Lattice (speed<0 is MINUS_END directed)
    real         lattice_flux_speed;
    
    /// loading rate of substance from Field to Lattice
    /**
     This is a binding rate per unit time and per unit length of Fiber.
     Binding is proportional to the concentration of substance in the field.
     */
    real         lattice_binding_rate;

    /// unloading rate of substance from Lattice to Field (unit is 1/second)
    real         lattice_unbinding_rate;
    
    /// flag controlling the forces exerted by Space on fiber points
    /**
     Possible values:
     - `off` (default)
     - `inside`
     - `outside`
     - `surface`
     - `plus_end`
     - `minus_end`
     - `both_ends`
     .
     */
    Confinement  confine;
    
    /// stiffness of confinement (also known as `confine[1]`)
    real         confine_stiffness;
    
    /// name of space used for confinement (also known as `confine[2]`)
    std::string  confine_space;
    
    /// if true, include steric interaction for this object
    /**
     The steric interaction generates a force derived from the potential energy:
     @code
     E = 1/2 k * ( d - d_0 ) ^ 2
     @endcode
     where `d` is the distance between two sections of filament. 
     The force is controlled by two parameters:
     - a stiffness `k`,
     - and equilibrium length `d_0`
     .
     
     This force is repulsive at short range ( d < d_0 ),
     and attractive elsewhere ( d > d_0 ).
     */
    int          steric;
    
    /// radius of repulsive steric interaction (also known as `steric[1]`)
    real         steric_radius;
    
    /// extra radius of attractive steric interaction (also known as `steric[2]`)
    real         steric_range;
    
    /// name of field
    std::string  field;
    
    /// type of glue (interaction between fiber PLUS_END and Space)
    /**
     Parameter fiber:glue is used to create interactions with the boundaries:
     - it creates a Single, everytime a fiber contacts the surface.
     - the Single is deleted if the associated Hand detaches.
     .
    */
    int          glue;
    
    /// name of Single used for glue (set a `glue[1]`)
    std::string  glue_single;
    
    
    /// create a force proportional to the length of the fiber, parallel to the fiber
    /**
     This has unit of force per unit length.
     a positive 'colinear_force' is directed toward the 'plus-end'
     a negative 'colinear_force' is directed toward the 'minus-end'
     */
    real         colinear_force;
    
    /// maximum speed of disassembly due to chewing
    real         max_chewing_speed;
    
    /// specialization
    /**
     @copydetails FiberGroup
     */
    std::string  activity;
    
    /// display string (see @ref FiberDispPar)
    std::string  display;
    
#ifdef BACKWARD_COMPATIBILITY
    /// add a force toward the X-axis
    int  squeeze;
    /// max norm of squeezing force (set as \c squeeze[1])
    real squeeze_force;
    /// range below which squeezing is linear (set as \c squeeze[2])
    real squeeze_range;
#endif
    
    /// @}

    /// derived variable: flag to indicate that `display` has a new value
    bool         display_fresh;
    
    /// derived variable: display
    FiberDisp *  disp;
    
    /// pointer to actual confinement Space, derived from `confine_space`
    Space const* confine_space_ptr;
    
    /// derived variable: pointer to associated Field
    Field *      field_ptr;
    
    /// derived variable: unbinding rates
    real         lattice_unbinding_prob;

protected:
    
    /// maximum speed of shrinkage
    real    max_chewing_speed_dt;

    /// local copy of SimulProp::time_step
    real    time_step;
    
    /// fraction in [0, 1]
    real    free_polymer;
    
    /// total length of fiber for this type
    real    total_length;
    
    /// SingleSet where glue are stored
    SingleSet  * glue_set;
    
    /// SingleProp used for glue
    SingleProp * glue_prop;

public:
    
    /// constructor
    FiberProp(const std::string& n) : Property(n), disp(0) { clear(); }
    
    /// destructor
    ~FiberProp() { }
    
    /// return a non-initialized Fiber with this property
    virtual Fiber* newFiber() const;
    
    /// return a Fiber with this property, initialized
    Fiber* newFiber(Glossary& opt) const;
    
    /// identifies the property
    std::string category() const { return "fiber"; }
    
    /// set default values
    virtual void clear();
       
    /// set using a Glossary
    virtual void read(Glossary&);
   
    /// check and derive parameter values
    virtual void complete(Simul const*);
    
    /// return a carbon copy of object
    Property* clone() const { return new FiberProp(*this); }

    /// write
    virtual void write_values(std::ostream &) const;
    
    
    /// ADDED BY MANU
    bool sweep_digits;
};

#endif


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FIBER_H
#define FIBER_H

#include <set>
#include <stdint.h>
#include "rigid_fiber.h"
#include "fiber_prop.h"
#include "node_list.h"
#include "field.h"
#include "lattice.h"
#include "array.h"
#include "sim.h"


class Single;
class FiberSet;
class FiberLocus;
class FiberBinder;
class LineDisp;


//#define NEW_CHEW_FIBERS


/**
 @todo The type of FiberLattice should be selected at run time,
 by the user via some parameter setting.
 */
#ifdef MANU_LATTICE
typedef Lattice<uint16_t> FiberLattice;
#else
typedef Lattice<real> FiberLattice;
#endif
/// Filament to which FiberBinder may bind
/**
 A Fiber extends the RigidFiber and Filament, definining all the functions
 that are necessary to the simulation.
 
 - FiberProp * prop points to the physical properties (ie. parameters) of the Fiber.
 - frBinders keeps track of all attached FiberBinders.
 - frLocuses provides pointers to the Segments of the Fiber.
 - FiberDisp * disp points to display parameters.
 .
 
 if prop->lattice is true, the Fiber has a Lattice,
 which is used by Digit and derived Hands.
 
 Fibers are stored in a FiberSet.
 */
class Fiber: public RigidFiber
{   
private:
    
    /// Disabled copy constructor
    Fiber(Fiber const&);

    /// class to store the info needed for sever()
    struct SeverPos
    {
        real    abs;
        int     staM;
        int     staP;
        
        SeverPos(real a, int p, int m) { abs=a; staP=p; staM=m; }
        /// sort from PLUS_END to MINUS_END, i.e. with decreasing abscissa
        real operator < (SeverPos const&b) const { return abs > b.abs; }
    };
    
    /// sorted list of future severing positions
    std::set<SeverPos>  futureCuts;

    /// list of attached FiberBinders
    NodeList            frBinders;
    
    /// Associated Lattice
    FiberLattice *      frLattice;
    
    /// array of rods, used in Attachments algorithm
    Array<FiberLocus>   frLocuses;
    
    /// a grafted used to immobilize the Fiber
    Single *            frGlue;
    
    /// time at birth
    real                frBirthTime;
    
protected:
    
    /// stored chewing at the end
    real                frChewM, frChewP;

    /// cut Fiber at point `pti`, return section `[ pti - PLUS_END ]`
    virtual Fiber* severPoint(unsigned int pti);
    
    /// cut fiber at points where consecutive segments make a kink
    void           severKinks();

    
    /// viscous drag coefficient for a cylinder moving close to a surface
    real           dragCoefficientSurface();
    
    /// viscous drag coefficient for a cylinder moving in an infinite volume of fluid
    real           dragCoefficientVolume();
    
public:
    
    /// the Property of this object
    FiberProp const*    prop;
    
    /// the display parameters
    LineDisp mutable*   disp;

    //--------------------------------------------------------------------------

    /// constructor
    Fiber(FiberProp const*);
    
    /// destructor
    virtual ~Fiber();
    
    //--------------------------------------------------------------------------
    
    /// allocate memory for 'nbp' points
    virtual unsigned allocatePoints(unsigned nbp);
    
    /// calculate viscous drag coefficient
    void           setDragCoefficient();
    
    /// prepare for Meca
    void           prepareMecable();
    
    /// add interactions to a Meca
    void           setInteractions(Meca &) const;
    

    /// invert polarity and adjust the abscissa of Hands such that they do not move
    void           flip();
    
    /// remove the portion of size `len` that includes the MINUS_END
    void           cutM(real len);
    
    /// remove the portion of size `len` that includes the PLUS_END
    void           cutP(real len);
    
    /// Cut all segments intersecting the plane defined by <em> n.pos + a = 0 </em>
    void           planarCut(Vector const& n, real a, int stateP, int stateM);
    
    /// cut fiber at distance `abs` from the MINUS_END; returns section `[ abs - PLUS_END ]`
    Fiber *        severM(real abs);

    /// cut fiber at abscissa `abs`; returns section `[ abs - PLUS_END ]`
    Fiber *        severNow(real abs) { return severM(abs-abscissaM()); }

    /// register a cut at abscissa `a` from the ORIGIN, with `m` and `p` the states of the new ends
    void           sever(real a, int p, int m) { futureCuts.insert(SeverPos(a, p, m)); }
    
    /// perform all the cuts registered by sever()
    void           severNow();

    /// register a chewing quantity
    void           chew(const real x, FiberEnd end) { if ( end == PLUS_END ) frChewP += x; else frChewM += x; }
    
    /// call Filament::join(), and transfer Hands (caller should delete `fib`).
    virtual void   join(Fiber * fib);
    
    /// simulation step
    virtual void   step();
    
    /// called if a Fiber tip has elongated or shortened
    void           update();

    /// update all FiberBinder bound to this
    void           updateBinders() const;
    
    //--------------------------------------------------------------------------
    
    /// returns simulation time at which Fiber was created
    real           birthTime() const { return frBirthTime; }

    /// set birth time
    void           birthTime(real t) { frBirthTime = t; }
    
    /// returns current age of the fiber
    real           age() const;

    /// the energy due to rigidity, which grows quadratically with curvature
    real           bendingEnergy() const { return bendingEnergy0()*prop->rigidity; }

    /// FiberLocus representing the segment [p, p+1]
    FiberLocus&    locus(unsigned int p) const;
    
    /// return the abscissa of the closest position to `w` on this Fiber, and set `dis` to the square of the distance
    real           projectPoint(Vector const& w, real & dis) const;
        
    //--------------------------------------------------------------------------
    
    /// return assembly/disassembly state of MINUS_END
    virtual unsigned dynamicStateM() const { return STATE_WHITE; }

    /// return assembly/disassembly state of PLUS_END
    virtual unsigned dynamicStateP() const { return STATE_WHITE; }

    /// return assembly/disassembly state of the FiberEnd
    unsigned         dynamicState(FiberEnd which) const;

    
    /// change state of MINUS_END
    virtual void   setDynamicStateM(unsigned s) {}

    /// change state of PLUS_END
    virtual void   setDynamicStateP(unsigned s) {}

    /// change state of FiberEnd `which` to `s`
    void           setDynamicState(FiberEnd which, unsigned s);
    
    
    /// the amount of freshly assembled polymer during the last time step (this has units of length)
    virtual real   freshAssemblyM() const { return 0; }

    /// the amount of freshly assembled polymer during the last time step (this has units of length)
    virtual real   freshAssemblyP() const { return 0; }

    /// the amount of freshly assembled polymer during the last time step (this has units of length)
    real           freshAssembly(FiberEnd which) const;
    
    
    /// true if the tip `which` has grown in the last time step ( freshAssembly(which) > 0 )
    bool           isGrowing(FiberEnd end) const { return freshAssembly(end) > 0; }
    
    /// true if the tip `which` has shrunk in the last time step ( freshAssembly(which) < 0 )
    bool           isShrinking(FiberEnd end) const { return freshAssembly(end) < 0; }
    
    //--------------------------------------------------------------------------
    
    /// register a new Binder attached to this Fiber
    void           addBinder(FiberBinder*);
    
    /// unregister bound Binder (which has detached)
    void           removeBinder(FiberBinder*);
    
    /// detach all binders
    void           detachBinders();
    
    /// sort binders by order of increasing abscissa
    void           sortBinders();
    
    /// a FiberBinder bound to this fiber (use ->next() to access all other binders)
    FiberBinder*   firstBinder() const;
   
    /// number of attached FiberBinders
    unsigned       nbBinders()   const { return frBinders.size(); }
    
    /// number of attached FiberBinders in a range of abscissa
    unsigned       nbBindersInRange(real abs_min, real abs_max, FiberEnd from) const;
    
    /// number of attached FiberBinders at the specified FiberEnd
    unsigned       nbBindersNearEnd(real len, FiberEnd which) const;
    
    /// a function to count binders using custom criteria
    unsigned       nbBinders(unsigned int (*count)(FiberBinder const&)) const;
    
    //--------------------------------------------------------------------------
    
    /// const pointer to Lattice
    FiberLattice*  lattice() const { return frLattice; }
    
    /// recalculate all lattice sites
    void           resetLattice() const;
    
    /// transfer all lattice substance to the Field
    void           releaseLattice(Field*) const;

    /// transfer a fraction `s` of the lattice substance to the Field
    void           releaseLattice(Field*, real s) const;
    
    /// initialize lattice sites to represent a constant linear concentration
    void           setLattice(real) const;
    
    /// real frLattice
    void           readLattice(Inputter& in);
    
    //--------------------------------------------------------------------------
    
    /// set the box glue for pure pushing
    void           setGlue1(Single* glue, FiberEnd, const Space * space);
    
    /// set the box glue for pure pulling
    void           setGlue2(Single* glue, FiberEnd, const Space * space);
    
    /// set the box glue for pushing and pulling
    void           setGlue3(Single* glue, const Space * space);
    
    /// a setGlue to rule them all
    void           setGlue(Single*& glue, FiberEnd, const Space * space, int mode);
    
    /// create a Single that can be used as glue
    void           makeGlue(Single*& glue);
    
    //--------------------------------------------------------------------------
    
    /// convert pointer to Fiber* is the conversion seems valid; returns 0 otherwise
    static Fiber* toFiber(Object * ptr)
    {
        if ( ptr && ptr->tag()==Fiber::TAG )
            return static_cast<Fiber*>(ptr);
        return 0;
    }

    /// a static_cast<> of Node::next()
    Fiber *  next()  const  { return static_cast<Fiber*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Fiber *  prev()  const  { return static_cast<Fiber*>(nPrev); }

    //--------------------------------------------------------------------------
    
    /// a unique character identifying the class
    static const Tag TAG = 'f';
    
    /// identifies data for dynamic ends of fibers
    static const Tag TAG_DYNAMIC = 'F';
    
    /// identifies FiberLattice data
    static const Tag TAG_LATTICE = 'l';
    
    /// return unique character identifying the class
    Tag         tag() const { return TAG; }
    
    /// return Object Property
    Property const* property() const { return prop; }
    
    ///write to file
    void        write(Outputter&) const;
    
    ///read from file
    void        read(Inputter&, Simul&, Tag);
    
    //---------------------------------------ADDED BY MANU
    bool        causesEntanglement(Vector const &, Vector const & ,const Fiber *) const;
    
    void        sweep();
    
    /// Useful for debugging (added by Manu)
    bool                frBindersBad() const {return frBinders.bad();} ;
    
#ifdef MULTI_LATTICE
    // An array with four positions containing the pointers to the fibers this fibers bind to.
    Fiber *  multi_lattice[4];
    
    unsigned int get_lattice_val(Fiber *);
    
    unsigned int available_lattice();
#endif
#ifdef MANU_LATTICE
    real         measure_cap(unsigned int lat_val);

    FiberBinder *  findHand(unsigned int lat_val, int site) const;
    
    FiberBinder *  findHand(unsigned int lat_val, int site, int fiber_id)  const;
#endif
};

#endif
    

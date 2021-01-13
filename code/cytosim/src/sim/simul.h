// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SIMUL_H
#define SIMUL_H

#include "assert_macro.h"
#include <iostream>
#include <stack>
#include <map>

#include "single_set.h"
#include "couple_set.h"
#include "fiber_set.h"
#include "bead_set.h"
#include "solid_set.h"
#include "sphere_set.h"
#include "organizer_set.h"
#include "field_set.h"
#include "space_set.h"
#include "event_set.h"
#include "fiber_grid.h"
#include "point_grid.h"
#include "simul_prop.h"
#include "property_list.h"
#include "field_values.h"
#include "field.h"
#include "meca.h"

class SimulProp;


/// string defining the start of a frame in the trajectory file
const static char FRAME_TAG[] = "Cytosim ";


/// Simulator class
class Simul
{
private:
    
    /// the first Space defined in the simulation
    Space const*        sSpace;
    
    /// for periodic boundary conditions
    Modulo              sModulo;
    
    /// Meca used to set and integrate the equations of motion of Mecables
    mutable Meca        sMeca;
    
    /// grid used for attachment of Hand to Fiber
    mutable FiberGrid   fiberGrid;
    
public:
    /// grid used for steric interaction between Fiber/Solid/Bead/Sphere
    mutable PointGrid   stericGrid;

    //-------------------------------------------------------------------------------
    
    /// a copy of the properties as they were stored to file
    mutable std::string properties_saved;
    
    /// flag for using precondition on the current time-step
    int                 precond;
    
    /// counter that is used to automatically set the preconditionning option
    unsigned int        precondCounter;
    
    /// stores cpu time to automatically set the preconditionning option
    unsigned long       precondCPU[4];
   
public:

    /// Global cytosim parameters
    SimulProp *      prop;
    
    /// list of all Object Property (does not include the SimulProp)
    PropertyList     properties;
    
    
    /// list of Space in the Simulation
    SpaceSet         spaces;
    
    /// list of Field in the Simulation
    FieldSet         fields;
    
    /// list of Fiber in the Simulation
    FiberSet         fibers;
    
    /// list of Sphere in the Simulation
    SphereSet        spheres;
    
    /// list of Bead in the Simulation
    BeadSet          beads;
    
    /// list of Solid in the Simulation
    SolidSet         solids;
    
    /// list of Single in the Simulation
    SingleSet        singles;
    
    /// list of Couple in the Simulation
    CoupleSet        couples;
    
    /// list of Organizer in the Simulation
    OrganizerSet     organizers;

    /// list of Events in the Simulation
    EventSet         events;
    
    /// list of couples that will undergo Gillespie steps in the simulation, allows for separation of
    /// timescales
    //CoupleSet        couples_gill;
    //-------------------------------------------------------------------------------
    
    /// constructor
    Simul();
    
    /// destructor
    virtual         ~Simul();
        
    //-------------------------------------------------------------------------------
    
    /// add Object to Simulation
    void             add(Object *);

    /// add all Objects from given list
    void             add(ObjectList);

    /// remove Object from Simulation
    void             remove(Object *);

    /// remove all Objects from given list
    void             remove(ObjectList);
    
    /// remove and delete object
    void             erase(Object *);
    
    /// remove and delete all objects from given list
    void             erase(ObjectList);

    /// mark objects from given list
    static void      mark(ObjectList, int);

    /// reset simulation world (clear all sub-lists and variables)
    void             erase();
    
    //-------------------------------------------------------------------------------

    /// total number of objects in the Simulation
    unsigned         nbObjects() const;

    /// time in the simulated world
    real             time()   const;
    
    //-------------------------------------------------------------------------------
   
    /// perform basic initialization; register callbacks
    void             initialize(Glossary&);
    
    /// change 'master' Space
    void             changeSpace(Space const*);

    /// return 'master' Space, usually the first one that was created, with the smallest `identity()`
    Space const*     space() const { return sSpace; }
    
    /// return first Space with given name
    Space const*     findSpace(std::string const& name) const;

    /// call foldPosition() for all objects (implements periodic boundary conditions)
    void             foldPosition() const;
    
    //-------------------------------------------------------------------------------
    
    /// ready the engine for a call to `step()` and `solve()`
    void             prepare();
    
    /// perform one Monte-Carlo step, corresponding to `time_step`
    void             step();
    
    /// this is called after a sequence of `step()` have been done
    void             relax() { singles.relax(); couples.relax(); }
    
    /// call setInteractions(meca) for all objects (this is called before `solve()`
    void             setInteractions(Meca&) const;

    /// simulate the mechanics of the system and move Mecables accordingly, corresponding to `time_step`
    void             solve();
    
    /// calculate the motion of objects, but only in the X-direction
    void             solveX();
    
    /// move every Fiber backward by `shift` (this is a extremely crude model)
    void             solveF(real shift);
    
    /// calculate Forces and Lagrange multipliers on the Mecables, but do not move them
    void             computeForces() const;
    
    /// dump matrix and vector from Meca
    void             dump() const { sMeca.dump(); }
    
private:
    
    /// give an estimate of the cell size of the FiberGrid
    real             estimateFiberGridStep() const;
    
    /// set FiberGrid and StericGrid over the given space
    void             setFiberGrid(Space const*) const;
    
    /// give an estimate for the cell size of the PointGrid used for steric interactions
    real             estimateStericRange() const;
    
    /// initialize the grid for steric interaction (stericGrid)
    void             setStericGrid(Space const*) const;
    
    /// add steric interactions between spheres, solids and fibers to Meca
    void             setStericInteractions(Meca&) const;
    
    //-------------------------------------------------------------------------------
    /// Function used to parse the config file, and to read state from a file:
    //-------------------------------------------------------------------------------

    /// return the ObjectSet corresponding to this Tag in the simulation (used for IO)
    ObjectSet*       findSetT(const Tag);
    
public:
    
    /// convert Object to Mecable* if the conversion seems valid; returns 0 otherwise
    static Mecable* toMecable(Object *);
    
    /// return the ObjectSet corresponding to a class
    ObjectSet*       findSet(const std::string& kind);
    
    /// return the ObjectSet corresponding to a class
    ObjectSet const* findSet(const std::string& kind) const
    {
        return const_cast<Simul*>(this)->findSet(kind);
    }
    
    /// find an object from a string specifying name and inventory number (e.g. 'fiber1')
    Object*         findObject(const std::string& spec) const;
    
    /// read an Object reference and return the corresponding Object (`tag` is set)
    Object*         readReference(Inputter&, Tag& tag);

    /// check if the name corresponds to a property class
    bool            isPropertyClass(const std::string&) const;
    
    /// return existing property of given class and name, or return zero
    Property*       findProperty(const std::string&, const std::string&) const;
    
    /// return existing property of given name, or return zero
    Property*       findProperty(const std::string&) const;

    /// return all existing properties of required class
    PropertyList    findAllProperties(const std::string&) const;

    /// return Property in the requested type
    template < typename T >
    T findProperty(std::string const& kd, std::string const& nm) const
    {
        Property * p = properties.find(kd, nm);
        return static_cast<T>(p);
    }
    
    /// return Property in the requested type
    template < typename T >
    T findProperty(std::string const& kd, unsigned ix) const
    {
        Property * p = properties.find(kd, ix);
        return static_cast<T>(p);
    }
   
    /// create a new property
    Property*       newProperty(const std::string&, const std::string&, Glossary&);
    
    /// export all Properties to speficied file
    void            writeProperties(std::ostream&, bool prune) const;
    
    /// export all Properties to a new file with specified name
    void            writeProperties(char const* filename, bool prune) const;

    //-------------------------------------------------------------------------------
    
    class InputLock;
    
    /// load the properties contained in the standard output property file
    void      loadProperties();
    
    /// read objects from file, and add them to the simulation state
    int       readObjects(Inputter&, ObjectSet* subset);

    /// read objects from a file, and add them to the simulation state
    int       loadObjects(Inputter&, ObjectSet* subset = 0);

    /// load sim-world from the named file
    int       loadObjects(char const* filename);
    
    /// import objects from file, and delete objects that were not referenced in the file
    int       reloadObjects(Inputter&);

    /// write sim-world to specified file
    void      writeObjects(Outputter&) const;
    
    /// write sim-world in binary or text mode, appending to existing file or creating new file
    void      writeObjects(char const* filename, bool binary, bool append) const;
    
    //-------------------------------------------------------------------------------

    /// call `Simul::report0`, adding lines before and after with 'start' and 'end' tags.
    void      report(std::ostream&, std::string const&, Glossary&) const;
    
    /// call one of the report function
    void      report0(std::ostream&, std::string const&, Glossary&) const;

    /// print time
    void      reportTime(std::ostream&) const;
    
    /// give a short inventory of the simulation state, obtained from ObjectSet::report()
    void      reportInventory(std::ostream&) const;
 
    /// print the length and the points of each fiber
    void      reportFiber(std::ostream&, FiberProp const*) const;
    
    /// print the length and the points of each fiber
    void      reportFiber(std::ostream&, std::string const&) const;
    
    /// print the length and the points of each fiber
    void      reportFiber(std::ostream&) const;
    
    /// print the coordinates of the model-points of each fiber
    void      reportFiberPoints(std::ostream&) const;
    
    /// print the mean and standard deviation of model-points of all fibers
    void      reportFiberMoments(std::ostream&) const;

    /// print the coordinates and forces on the model-points of each fiber
    void      reportFiberForces(std::ostream&) const;
    
    /// print radial component of forces experienced by Fibers due to confinement
    void      reportFiberConfinement(std::ostream& out) const;

    /// print the positions and the states of the two ends of each fiber
    void      reportFiberEnds(std::ostream&) const;
    
    /// print average age and standard deviation for each class of fiber
    void      reportFiberAge(std::ostream&) const;

    /// print average length and standard deviation for each class of fiber
    void      reportFiberLengths(std::ostream&) const;
    
    /// print length distribution for each class of fiber
    void      reportFiberLengthDistribution(std::ostream&, Glossary&) const;

    /// print number of kinks in each class of Fiber
    void      reportFiberSegments(std::ostream&) const;
    
    /// print number of fibers according to dynamic state of end
    void      reportFiberDynamic(std::ostream&, FiberEnd) const;
    
    /// print number of fibers according to their dynamic states
    void      reportFiberDynamic(std::ostream&) const;
    
    /// print coordinates of speckles that follow a frozen random sampling
    void      reportFiberSpeckles(std::ostream&, Glossary&) const;
    
    /// print coordinates of points randomly and freshly distributed
    void      reportFiberSamples(std::ostream&, Glossary&) const;

    /// print dynamic states of Fiber
    void      reportFiberStates(std::ostream&) const;
    
    /// print Fiber tensions along certain planes defined in `opt`
    void      reportFiberTension(std::ostream&, Glossary&) const;
    
    /// print sum of all bending energy
    void      reportFiberBendingEnergy(std::ostream&) const;
    
    /// print positions of interection between two fibers
    void      reportFiberIntersections(std::ostream&, Glossary&) const;
    
    /// print interection abscissa between fibers
    void      reportFiberBinders(std::ostream&) const;

    void      reportFiberBindersAll(std::ostream&) const;
    
    /// print interection abscissa between fibers
    void      reportFiberConnectors(std::ostream&, Glossary&) const;
    
    // print the length of the caps
    void      reportFiberCap(std::ostream&, Glossary&) const;
    
    // print the length of the overlap with neighboring fibers
    void      reportFiberOverlap(std::ostream&, Glossary&) const;
    
    /// print the number of fibers
    void      reportFiberCount(std::ostream&, Glossary&) const;
    
    /// print interection abscissa between fibers
    void      reportNetworkBridges(std::ostream&, Glossary&) const;


    /// print Organizer positions
    void      reportOrganizer(std::ostream&) const;

    /// print Aster positions
    void      reportAster(std::ostream&) const;
    
    /// print Bead positions 
    void      reportBeadSingles(std::ostream&) const;

    /// print Bead positions
    void      reportBeadPosition(std::ostream&) const;

    /// print Solid positions 
    void      reportSolid(std::ostream&) const;

    /// print state of Couples 
    void      reportCouple(std::ostream&) const;
    
    /// print state of Couples
    void      reportCoupleHand(std::ostream&) const;
    
    /// print position of Couples 
    void      reportCoupleState(std::ostream&) const;
    
    /// print position of Couples of a certain kind
    void      reportCoupleState(std::ostream&, std::string const&) const;
    
    /// print position of active Couples of a certain kind
    void      reportCoupleActive(std::ostream&, std::string const&) const;
    
    /// print position and forces of doubly-attached Couples
    void      reportCoupleLink(std::ostream&, std::string const&) const;
    
    /// Another mode added by Manu
    void      reportCouplePosition(std::ostream&, std::string const&) const;
    
    /// print position and forces of Couples of a certain kind
    void      reportCoupleForce(std::ostream&, Glossary&) const;
    
    /// print the same as reportCouple
    void      reportFlipper(std::ostream&) const;
    /// print state of Singles
    void      reportSingle(std::ostream&) const;
    
    /// print position of Singles
    void      reportSingleState(std::ostream&) const;
    
    /// print position of Singles of a certain kind
    void      reportSingleState(std::ostream&, std::string const&) const;

    /// print position of Singles
    void      reportSinglePosition(std::ostream&, std::string const&) const;
   
    /// print state of Couples 
    void      reportSphere(std::ostream&) const;

    /// print something about Spaces
    void      reportSpace(std::ostream&) const;
    
    /// print something about Fields
    void      reportField(std::ostream&) const;
#ifdef  TRAP_SINGLES
    /// print the bound/unbound/trapped
    void      reportCoupleTrap(std::ostream&) const;
    
    /// print the bound/unbound/trapped
    void      reportSingleTrap(std::ostream&) const;
#endif
    //-------------------------------------------------------------------------------
    
    /// analyse the network connectivity to identify isolated sub-networks
    void      flagClusters() const;
    
    /// analyse the network connectivity to identify isolated sub-networks
    void      flagClusters(Property const*) const;
    
    /// print size of clusters defined by connections with Couples
    void      reportClusters(std::ostream&, bool) const;
    
    /// estimates if Fibers form a connected ring around the Z-axis
    void      reportRing(std::ostream&) const;

    /// print Aster & Spindle indices
    void      reportIndices(std::ostream&) const;

    /// print number of Fibers pointing left and right that intersect plane YZ at different X positions
    void      reportProfile(std::ostream&) const;

    /// a special print for Romain Gibeaux
    void      reportAshbya(std::ostream&) const;
    
    /// print something
    void      reportCustom(std::ostream&) const;

    //-------------------------------------------------------------------------------
    
    /// custom function
    void      custom0(Glossary&);
    /// custom function
    void      custom1(Glossary&);
    /// custom function
    void      custom2(Glossary&);
    /// custom function
    void      custom3(Glossary&);
    /// custom function
    void      custom4(Glossary&);
    /// custom function
    void      custom5(Glossary&);
    /// custom function
    void      custom6(Glossary&);
    /// custom function
    void      custom7(Glossary&);
    /// custom function
    void      custom8(Glossary&);
    /// custom function
    void      custom9(Glossary&);
    
    //ADDED BY MANU
    /// kill the simulation if the parameters given are not met
    void      spindle_watch();
    
#ifdef TRAP_SINGLES
#if (TRAP_SINGLES==2)
    void      trap_solution();
#endif
#endif

};

#endif


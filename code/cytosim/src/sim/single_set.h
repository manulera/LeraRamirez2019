// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SINGLE_SET_H
#define SINGLE_SET_H

#include "object_set.h"
#include "single_prop.h"
#include "single.h"
#include <deque>

/// Set for Single
/**
 A Single is stored in one of 2 NodeList, depending on its state:
 - fList = free,
 - aList = attached.
 .
 Each list is accessible via its head firstF() and firstA(),
 and subsequent objects obtained with next() are in the same state.
 This way, the state of the Single are known when accessing them.
 
 A Single is automatically transfered to the appropriate list,
 if its Hand binds or unbinds. This is done by the HandMonitor.
 */
class SingleSet: public ObjectSet
{
    
private:
    
    /// List for non-attached Singles (f=free)
    NodeList     fList;
    
    /// List for attached Singles (a=attached)
    NodeList     aList;
    
#ifdef TRAP_SINGLES
    /// Added by Manu list of trapped singles. The reason for keeping them in a list is that at the end of the program when calling ~Simul, it seems that the "unlinked" singles would give a proble. They are not truly unliked though, since they still keep the pointer to the object set, they just dont belong to any list, so if you call them it gives an error.
    NodeList     trapList;
    
#endif
    
    /// a list to hold Couples of one class
    typedef std::deque<Single*> SingleList;
    
    /// an array of CoupleList
    typedef std::vector<SingleList> SingleReserve;
    
    /// uniLists[p] contains the Single with ( property()->index() == p ) that are diffusing
    SingleReserve uniLists;
    
    /// flag to enable couple:fast_diffusion attachment algorithm
    bool          uni;
    
    /// initialize couple:fast_diffusion attachment algorithm
    bool          uniPrepare(PropertyList const& properties);
    
    /// distribute Singles to fiber with average distance `spread`
    void          uniAttach(Array<FiberBinder>&, SingleList&);
    
    /// couple:fast_diffusion attachment algorithm; assumes free Singles are uniformly distributed
    void          uniAttach(FiberSet const&);
    
    /// return Couples in uniLists to the normal lists
    void          uniRelax();

public:
        
    ///creator
    SingleSet(Simul& s) : ObjectSet(s) {}
    
    ///destructor
    virtual      ~SingleSet() {}
    
    //--------------------------

    /// identifies the class
    std::string   title() const { return "single"; }
    
    /// create a new property for class `kind` with given name
    Property *    newProperty(const std::string& kind, const std::string& name, Glossary&) const;
    
    /// create new objects, of class `kind` and type `name`, give the options provided in `opt`
    ObjectList    newObjects(const std::string& name, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *      newObjectT(Tag, unsigned);
    
    //--------------------------

    /// add object
    void          link(Object *);

    /// remove object
    void          unlink(Object *);
    
    /// reassign Single to different sublist following attachement of Hand
    void          relinkA(Single *);
    
    /// reassign Single to different sublist following detachment of Hand
    void          relinkD(Single *);
    
    
    /// create Wrists anchored on given Mecable
    ObjectList    makeWrists(Mecable const*, unsigned, unsigned, std::string&);

    /// return all Wrists anchored on `obj`
    ObjectList    collectWrists(Object * obj) const;
    
    /// remove all Wrists anchored on `obj`
    void          removeWrists(Object * obj);
    
    
    ///returns the first free Single
    Single *      firstF()       const { return static_cast<Single*>(fList.front()); }
    
    ///returns the first bound Single
    Single *      firstA()       const { return static_cast<Single*>(aList.front()); }
#ifdef TRAP_SINGLES
    ///returns the first trapped Single
    Single *      firstTrapped()       const { return static_cast<Single*>(trapList.front()); }
#endif
    /// return pointer to the Object of given ID, or zero if not found
    Single *      findID(ObjectID n) const { return static_cast<Single*>(inventory.get(n)); }
    
    /// collect all objects
    ObjectList    collect() const;
    
    /// collect objects for which func(obj, val) == true
    ObjectList    collect(bool (*func)(Object const*, void const*), void const*) const;

    /// erase all Object and all Property
    void          erase();
    
    /// number of unattached Simgles
    unsigned int  sizeF() const { return fList.size(); }
    
    /// number of attached Singles
    unsigned int  sizeA() const { return aList.size(); }

    
#ifdef TRAP_SINGLES
    /// number of trapped Singles
    unsigned int  sizeTrapped() const { return trapList.size(); }
    /// number of elements
    unsigned int  size()  const { return fList.size() + aList.size() + trapList.size(); }
    
    /// number of elements
    void  checkTrapped()  const;
    
#else
    /// number of elements
    unsigned int  size()  const { return fList.size() + aList.size(); }
#endif
    
    /// mix order of elements
    void          mix();

    /// mark object before import
    void          freeze();
    
    /// delete marked object after import
    void          prune();

    /// unmark objects after import
    void          thaw();
    
    /// prepare for step()
    void          prepare(PropertyList const& properties);
    
    /// Monte-Carlo step
    void          step(FiberSet const& fibers, FiberGrid const&);
    
    /// cleanup at end of simulation period
    void          relax() { uniRelax(); }
    
    
    /// print a summary of the content (nb of objects, class)
    void          report(std::ostream&) const;

    /// write
    void          write(Outputter&) const;
    
    /// modulo the position (periodic boundary conditions)
    void          foldPosition(const Modulo *) const;
    
    /// check internal consistency, returns 0 if everything is OK
    int           bad() const;

#ifdef TRAP_SINGLES
    /// ADDED BY MANU: a method to unlink/relink the object knowing that its trapped/untrapped
    void          untrap_single(Single *);
    
    void          trap_single(Single *);
    

#if (TRAP_SINGLES==2)
    typedef std::deque<HandMonitor * > HaMonList;
    
    HaMonList simulTrapReserve;
    
    void populateTrapList();
    
    void clearTrapReserve(){simulTrapReserve.clear();}
#endif

#endif
    
};


#endif


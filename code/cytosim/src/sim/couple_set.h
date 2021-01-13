// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef COUPLE_SET_H
#define COUPLE_SET_H

#include "object_set.h"
#include "couple.h"
#include "couple_prop.h"
#include <deque>

/// Set for Couple
/**
 A Couple is stored in one of 4 NodeList, depending on its state:
 - ffList = hand1 and hand2 unattached,
 - afList = hand1 attached, hand2 unattached,
 - faList = hand1 unattached, hand2 attached,
 - aaList = hand1 and hand2 attached [the couple makes a `link`].
 .
 The head of each list is accessible via firstFF() and firstFA(), firstAF() and firstAA(),
 and subsequent objects obtained with next() are in the same state.
 This makes iterating through the list of Couple more efficient.
 
 A Couple is automatically transfered to the appropriate list,
 if one of its Hand binds or unbinds. This is done by the HandMonitor.
 */
class CoupleSet: public ObjectSet
{
    
private:
    
    /// list of Couple which are not attached (f=free)
    NodeList    ffList;
    
    /// list of Couple with only one side attached (a=attached, f=free)
    NodeList    afList, faList;
    
    /// list of Couple with both sides attached (a=attached)
    NodeList    aaList;
    
    /// a list to hold Couples of one class
    typedef std::deque<Couple*> CoupleList;
    
    /// an array of CoupleList
    typedef std::vector<CoupleList> CoupleReserve;
    
    /// uniLists[p] contains the Couples with ( property()->index() == p ) that are diffusing
    CoupleReserve uniLists;

    /// flag to enable couple:fast_diffusion attachment algorithm
    bool          uni;
    
    /// gather all Couple with fast_diffusion in dedicated lists
    void          uniCollect();
    
    /// initialize couple:fast_diffusion attachment algorithm
    bool          uniPrepare(PropertyList const& properties);
    
    /// attach Hand1 of Couple on locations specified by first argument
    void          uniAttach1(Array<FiberBinder>&, CoupleList&);
    
    /// attach Hand2 of Couple on locations specified by first argument
    void          uniAttach2(Array<FiberBinder>&, CoupleList&);
    
    /// attach both Hands of `nb` Couple at crossing points specified by first argument
    void          uniAttach12(Array<FiberBinder>&, CoupleList&, int nb);
    
    /// couple:fast_diffusion attachment algorithm; assumes free Couples are uniformly distributed
    void          uniAttach(FiberSet const&);
    
    /// return Couples in uniLists to the normal lists
    void          uniRelax();
    
public:
    
    ///creator
    CoupleSet(Simul& s) : ObjectSet(s), uni(false) {}
    
    ///destructor
    virtual ~CoupleSet() {}
    
    //--------------------------
    
    /// identifies the set
    std::string  title() const { return "couple"; }
    
    /// create a new property for class `kind` with given name
    Property *   newProperty(const std::string& kind, const std::string& name, Glossary&) const;
    
    /// create new objects, of class `kind` and type `name`, give the options provided in `opt`
    ObjectList   newObjects(const std::string& name, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *     newObjectT(Tag, unsigned);

    //--------------------------
    
    /// add object (should be a Couple)
    void         link(Object *);
    
    /// remove object (should be a Couple)
    void         unlink(Object *);

    /// reassign Couple to different sublist following attachement of Hand 1
    void         relinkA1(Couple *);
    /// reassign Couple to different sublist following detachment of Hand 1
    void         relinkD1(Couple *);
    /// reassign Couple to different sublist following attachement of Hand 2
    void         relinkA2(Couple *);
    /// reassign Couple to different sublist following detachment of Hand 2
    void         relinkD2(Couple *);

    /// first unattached Couple
    Couple *     firstFF()     const { return static_cast<Couple*>(ffList.front()); }
    /// first Couple attached by cHand1
    Couple *     firstAF()     const { return static_cast<Couple*>(afList.front()); }
    /// first Couple attached by cHand2
    Couple *     firstFA()     const { return static_cast<Couple*>(faList.front()); }
    /// first Couple attached by both hands
    Couple *     firstAA()     const { return static_cast<Couple*>(aaList.front()); }

    /// last unattached Couple
    Couple *     lastFF()      const { return static_cast<Couple*>(ffList.back()); }
    /// last Couple attached by cHand1
    Couple *     lastAF()      const { return static_cast<Couple*>(afList.back()); }
    /// last Couple attached by cHand2
    Couple *     lastFA()      const { return static_cast<Couple*>(faList.back()); }
    /// last Couple attached by both hands
    Couple *     lastAA()      const { return static_cast<Couple*>(aaList.back()); }

    /// number of free Couples
    unsigned     sizeFF()      const { return ffList.size(); }
    /// number of Couples attached by cHand1 only
    unsigned     sizeAF()      const { return afList.size(); }
    /// number of Couples attached by cHand2 only
    unsigned     sizeFA()      const { return faList.size(); }
    /// number of Couples attached by both hands
    unsigned     sizeAA()      const { return aaList.size(); }
    /// total number of elements
    unsigned     size()        const { return ffList.size() + faList.size() + afList.size() + aaList.size(); }
    
    /// return pointer to the Object of given ID, or zero if not found
    Couple *     findID(ObjectID n)        const { return static_cast<Couple*>(inventory.get(n)); }
    
    /// first Couple in inventory
    Couple *     firstID()                 const { return static_cast<Couple*>(inventory.first()); }
    
    /// next Couple in inventory
    Couple *     nextID(Couple const* obj) const { return static_cast<Couple*>(inventory.next(obj)); }
    
    /// collect all objects
    ObjectList   collect() const;
    
    /// collect objects for which func(this, val) == true
    ObjectList   collect(bool (*func)(Object const*, void const*), void const*) const;

    /// erase all Object and all Property
    void         erase();
    
    /// mix order of elements
    void         mix();

    /// distribute the Couple on the fibers to approximate an equilibrated state
    void         equilibrateSym(FiberSet const&, CoupleList&, CoupleProp const*);

    /// distribute Couples of given class on the fibers to approximate an equilibrated state
    void         equilibrate(FiberSet const&, CoupleList&, CoupleProp const*);
    
    /// distribute all Couple on the fibers to approximate an equilibrated state
    void         equilibrate(FiberSet const&, PropertyList const&);
    
    /// distribute all free Couples on filament intersections
    void         connect(FiberSet const&, PropertyList const&);
    
    /// prepare for step()
    void         prepare(PropertyList const& properties);
    
    /// Monte-Carlo step
    void         step(FiberSet const&, FiberGrid const&);
    
    /// cleanup at end of simulation period
    void         relax() { uniRelax(); }
     
    /// mark object before import
    void         freeze();
    
    /// delete marked object after import
    void         prune();
    
    /// unmark objects after import
    void         thaw();

    /// write
    void         write(Outputter&) const;
    
    /// print a summary of the content (nb of objects, class)
    void         report(std::ostream&) const;

    /// modulo the position (periodic boundary conditions)
    void         foldPosition(const Modulo *) const;
    
    ///debug function
    int          bad() const;
    
    /// Gillespie step
    void         gilles_step(FiberSet const&, FiberGrid const&);
    
    /// Return a pointer to the CoupleSet for the trapping
    CoupleReserve *         getUniLists(){return &uniLists;};
    
#ifdef TRAP_SINGLES
#if (TRAP_SINGLES==2)
    typedef std::deque<HandMonitor * > HaMonList;
    
    HaMonList simulTrapReserve;
    
    HaMonList simulUntrapReserve;
    
    void populateTrapList();
    
    void trap(HaMonList * singles);
    
    void untrap();
    
    void clearTrapReserve(){simulUntrapReserve.clear(); simulTrapReserve.clear();}
#endif
#endif
};


#endif


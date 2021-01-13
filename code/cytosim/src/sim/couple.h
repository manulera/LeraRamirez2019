// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef COUPLE_H
#define COUPLE_H

#include "object.h"
#include "hand_monitor.h"
#include "couple_prop.h"
#include "hand.h"

class Meca;
class Glossary;
class FiberGrid;


/// A set of two Hand linked by an elastic element
/**
 A Couple contains two pointers to Hand:
 - cHand1
 - cHand2
 .
 There are 4 possible states for a Couple:
 - state FF (0): cHand1 and cHand2 are free,
 - state AF (1): cHand1 is bound, cHand2 is free,
 - state FA (2): cHand1 is free, cHand2 is bound,
 - state AA (3): both hands are attached
 .
 The method state() return the state of the Couple in [0-3].

 Generally the Couple behaves according to its state:
 - FF     : the Couple is diffusing and both Hands are trying to bind fibers,
 - AF, FA : the localization is given by the attachement point on the fiber,
 - AA     : the Couple is acting as a Hookean spring between the two fibers.
 .
 
 The default Couple has:
 - a zero resting length (it uses Meca:interLink())
 - no specificity
 .

 @ingroup CoupleGroup
 */
class Couple : public Object, public HandMonitor
{
public:

    /// associated properties
    CoupleProp const* prop;
    
protected:
    
    /// position and position in previous step of complex
    Vector   cPos;
    
    /// first Hand
    Hand    * cHand1;

    /// second Hand
    Hand    * cHand2;
    
    /// specialization of HandMonitor
    bool      allowAttachment(const FiberBinder &);
    /// specialization of HandMonitor
    void      afterAttachment(Hand const*);
    /// specialization of HandMonitor
    void      beforeDetachment(Hand const*);
    /// specialization of HandMonitor
    ObjectID  nucleatorID() const { return Object::identity(); }
    /// specialization of HandMonitor
    Hand *    otherHand(Hand const*) const;
    /// specialization of HandMonitor
    Vector    otherPosition(Hand const*) const;
    /// specialization of HandMonitor
    Vector    otherDirection(Hand const*) const;
    /// specialization of HandMonitor
    real      interactionLength() const;
    
public:
    
    /// create following the specifications in the CoupleProp
    Couple(CoupleProp const*, Vector const & w = Vector(0,0,0));

    /// destructor
    virtual ~Couple();

    /// copy operator
    Couple&  operator=(Couple const&);
    
    //--------------------------------------------------------------------------
    
    /// change the property and update the two Hands
    void           setProperty(CoupleProp *);
    
    /// add interactions to the Meca
    virtual void   setInteractions(Meca &) const;
    
    /// add interactions to the Meca (experimental)
    virtual void   setInteractionsAF(Meca &) const;
    
    /// add interactions to the Meca (experimental)
    virtual void   setInteractionsFA(Meca &) const;
    
    //--------------------------------------------------------------------------
    
    /// the position of the complex, calculated from cPos, cHand1 and cHand2
    virtual Vector position() const;
   
    /// Couple can be displaced only if it is not attached
    virtual bool   mobile()               const { return !cHand1->attached() && !cHand2->attached(); }
    
    /// translate object's position by the given vector
    virtual void   translate(Vector const& w)   { cPos += w; }
    
    /// move object to specified position
    virtual void   setPosition(Vector const& w) { cPos = w; }

    /// modulo the current position vector in the space
    virtual void   foldPosition(const Modulo*);
    
    /// set the position randomly inside prop->confine_space
    void           randomizePosition();
    
    //--------------------------------------------------------------------------
    
    /// activity flag
    virtual bool   active()               const { return true; }
    
    /// true if both Hands are attached
    bool           linking()              const { return cHand1->attached() && cHand2->attached(); }

    /// the state of the Couple in { 0 ... 3 } representing { FF, FA, FA, AA }
    int            state()                const { return cHand1->attached() + 2 * cHand2->attached(); }
    
    ///stiffness of the link ( = prop->stiffness )
    real           stiffness()            const;
    
    /// return one of the Hand that is attached, or zero if both are detached
    Hand *         attachedHand()         const;
    
    /// force between hands, essentially: stiffness * ( cHand2->posHand() - cHand1->posHand() )
    virtual Vector force()                const;
     
    /// cosine of the angle between the two Fibers attached by the hands
    real           cosAngle()             const { return cHand1->dirFiber() * cHand2->dirFiber(); }
   
    /// position on the side of fiber1 used for sideInteractions
    virtual Vector posSide()              const { return cHand1->pos(); }
    
    /// the position of the complex if it is unattached
    Vector         posFree()              const { return cPos; }
   
    //--------------------------------------------------------------------------

    /// simulation step for a free Couple: diffusion
    virtual void   stepFF(const FiberGrid&);
    
    /// simulation step for a Couple attached by Hand1
    virtual void   stepAF(const FiberGrid&);
    
    /// simulation step for a Couple attached by Hand2
    virtual void   stepFA(const FiberGrid&);
    
    /// simulation step for a doubly-attached Couple
    virtual void   stepAA();
    
    //Added by Manu -> steps of gillespie
    real new_gilles_t() const;
    
    /// simulation step for a doubly-attached Couple
    int gillestep();
    
    /// The gillespie time until the next step
    real gilles_t;
    
    
    //--------------------------------------------------------------------------

    /// pointer to Hand1
    Hand const*    hand1()                              const { return cHand1; }
    
    /// true if Hand1 is attached
    bool           attached1()                          const { return cHand1->attached(); }
    
    /// Fiber to which Hand1 is attached, or zero if not attached
    Fiber*         fiber1()                             const { return cHand1->fiber(); }
    
    /// attachment position of Hand1 along fiber (only valid if Hand1 is attached)
    real           abscissa1()                          const { return cHand1->abscissa(); }
    
    /// position of Hand1 when attached (only valid if Hand1 is attached)
    Vector         posHand1()                           const { return cHand1->pos(); }
    
    /// direction of Fiber at attachment point of Hand1 (only valid if Hand1 is attached)
    Vector         dirFiber1()                          const { return cHand1->dirFiber(); }
    
    /// attach Hand1 at given abcissa
    void           attachTo1(Fiber* f, real ab)               { cHand1->attachTo(f, ab); }

    /// attach Hand1 at specified position
    void           attachTo1(Fiber* f, real ab, FiberEnd ref) { cHand1->attachTo(f, ab, ref); }
    
    /// attach Hand1 at the given end
    void           attachToEnd1(Fiber* f, FiberEnd end)       { cHand1->attachToEnd(f, end); }

    /// attach Hand1 at the given FiberBinder
    void           attach1(FiberBinder const& fb)             { cHand1->attach(fb); }
    
    /// move Hand1 to given end
    void           moveToEnd1(FiberEnd end)                   { cHand1->moveToEnd(end); }

    //--------------------------------------------------------------------------

    /// pointer to Hand2
    Hand const*    hand2()                              const { return cHand2; }
    
    /// true if Hand2 is attached
    bool           attached2()                          const { return cHand2->attached(); }
    
    /// Fiber to which Hand2 is attached, or zero if not attached
    Fiber*         fiber2()                             const { return cHand2->fiber(); }
    
    /// attachment position of Hand2 along fiber (only valid if Hand2 is attached)
    real           abscissa2()                          const { return cHand2->abscissa(); }
    
    /// position of Hand2 when attached (only valid if Hand2 is attached)
    Vector         posHand2()                           const { return cHand2->pos(); }
    
    /// direction of Fiber at attachment point of Hand2 (only valid if Hand2 is attached)
    Vector         dirFiber2()                          const { return cHand2->dirFiber(); }
    
    /// attach Hand2 at given abcissa
    void           attachTo2(Fiber* f, real ab)               { cHand2->attachTo(f, ab); }
    
    /// attach Hand2 at specified position
    void           attachTo2(Fiber* f, real ab, FiberEnd ref) { cHand2->attachTo(f, ab, ref); }
    
    /// attach Hand2 at the given end
    void           attachToEnd2(Fiber *f, FiberEnd end)       { cHand2->attachToEnd(f, end); }
    
    /// attach Hand2 at the given FiberBinder
    void           attach2(FiberBinder const& fb)             { cHand2->attach(fb); }
    
    /// move Hand2 to given end
    void           moveToEnd2(FiberEnd end)                   { cHand2->moveToEnd(end); }

    //--------------------------------------------------------------------------

    /// a static_cast<> of Node::next()
    Couple *       next()        const { return static_cast<Couple*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Couple *       prev()        const { return static_cast<Couple*>(nPrev); }
    
    //------------------------------ read/write --------------------------------

    /// a unique character identifying the class
    static const Tag TAG = 'c';
    
    /// return unique character identifying the class
    Tag            tag() const { return TAG; }
    
    /// return Object Property
    Property const* property() const { return prop; }
    
    /// write to file
    void           write(Outputter&) const;
    
    /// read from file
    void           read(Inputter&, Simul&, Tag);
    
    /// return PointDisp of Hand1
    PointDisp *    disp1() const { return cHand1->prop->disp; }
    
    /// return PointDisp of Hand2
    PointDisp *    disp2() const { return cHand2->prop->disp; }
#ifdef TRAP_SINGLES
    /// Display trapped molecule (for trapper)
    virtual int           display_trap(Vector&, Fiber *& , PointDisp*& ){return false;};
    
    void getHands(Hand * &h1, Hand * &h2) const {h1 = cHand1; h2 =cHand2;};
#endif
};


#endif


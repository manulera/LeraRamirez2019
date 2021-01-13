// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "single_set.h"
#include "single_prop.h"
#include "glossary.h"
#include "iowrapper.h"

#include "simul.h"
#include "wrist.h"
#include "wrist_long.h"

//------------------------------------------------------------------------------
/**
 @copydetails SingleGroup
 */
Property* SingleSet::newProperty(const std::string& kd, const std::string& nm, Glossary& opt) const
{
    if ( kd == "single" )
        return new SingleProp(nm);
    else
        return 0;
}

//------------------------------------------------------------------------------

void SingleSet::prepare(PropertyList const& properties)
{
    uni = uniPrepare(properties);
}


void SingleSet::step(FiberSet const& fibers, FiberGrid const& fgrid)
{
    // use alternate attachment strategy:
    if ( uni )
        uniAttach(fibers);

    /*
     ATTENTION: we have multiple lists, and Objects are automatically transfered
     from one list to another if their Hands bind or unbind. We ensure here that
     step() is called exactly once for each object. THe code relies on the fact
     that a transfered node would be linked at the start of the new list.
     We start always at the node, which was first before any transfer could occur.
     */
    
    /*
     MSG(9,"CoupleSet::stepCouple entry : F %5i AF %5i FA %5i AA %5i\n",
     ffList.size(), afList.size(), faList.size(), bList.size());
     */
    
    Single *const fHead = firstF();
    Single * obj, * nxt;
    
    obj = firstA();
    while ( obj )
    {
        nxt = obj->next();
        obj->stepA();
        if ( ! nxt ) break;
        obj = nxt->next();
        nxt->stepA();
    }
    
    obj = fHead;
    while ( obj )
    {
        nxt = obj->next();
        obj->stepF(fgrid);
        if ( ! nxt ) break;
        obj = nxt->next();
        nxt->stepF(fgrid);
    }
}


//------------------------------------------------------------------------------
#pragma mark -

Object * SingleSet::newObjectT(const Tag tag, unsigned idx)
{
    SingleProp * p = simul.findProperty<SingleProp*>("single", idx);
    if ( p == 0 )
        throw InvalidIO("no single class defined with id "+sMath::repr(idx));

    if ( tag == Wrist::TAG )
        return p->newWrist(0, 0);
    
    if ( tag == Single::TAG )
        return p->newSingle();

    return 0;
}

/**
 @addtogroup SingleGroup
 
 A newly created Single can be anchored to a Mecable:
 @code
 new single NAME {
   base = OBJECT, POINT
 }
 @endcode
 
 where:
 - OBJECT is the concatenation of the class name with the serial number of the object:
     - 'bead1' for the first bead
     - 'bead2' for the second...
     - 'bead0' designates the last bead made,
     - 'bead-1' is the penultimate one, etc.
     .
 - POINT designates a point on this object:
     - point0 = first point
     - point1 = second point...
     .
 .

 You can attach a Single to a fiber:
 @code
 new single protein
 {
    attach = FIBER, REAL, REFERENCE
 }
 @endcode
 
 where:
 - FIBER designates the fiber:
     - `fiber1` of `fiber2` correspond to fibers directly
     - `first` or `last` to the oldest and youngest fiber
     - `fiber-1` the penultimate, etc.
     .
 - REAL is the abscissa of the attachment point.
   If the abscissa is not specified, and random position along
   along the fiber will be selected.
 - REFERENCE can be `minus_end`, `center` or `plus_end` (default = `origin`).
   This defines from which position the abscissa is measured.
 .
 */
ObjectList SingleSet::newObjects(const std::string& name, Glossary& opt)
{
    Property * p = simul.properties.find_or_die("single", name);
    SingleProp * sp = static_cast<SingleProp*>(p);
        
    Single * obj = 0;
    std::string str;
    if ( opt.set(str, "base") )
    {
        Object * omec = simul.findObject(str);
        if ( omec == 0 )
            throw InvalidParameter("Could not find Mecable specified in single:base");
        Mecable * mec = Simul::toMecable(omec);
        if ( mec == 0 )
            throw InvalidParameter("Object specified in single:base is not a Mecable");
        // get index of point in second argument
        unsigned ip = 0;
        if ( opt.set(str, "base", 1) )
        {
            if ( str.compare(0,5,"point") == 0 )
                ip = strtol(str.c_str()+5, 0, 10);
            else
                ip = strtol(str.c_str(), 0, 10);
        }
        //std::clog << "base = " << str << "|" << io << "|" << ip << std::endl;
        if ( ip >= mec->nbPoints() )
            throw InvalidParameter("point index single:base[1] is out of range");
        
        obj = sp->newWrist(mec, ip);
    }
    else
        obj = sp->newSingle();

 
    ObjectList res;
    res.push_back(obj);
    
    /*
     This provides a way for the user to attach the Single to an existing fiber
     */
    std::string spe;
    if ( opt.set(spe, "attach") )
    {
        ObjectList fibs = simul.fibers.findObjects(spe);
        if ( fibs.empty() )
            throw InvalidParameter("Could not find Fiber in single::attach");
        
        Fiber* fib = Fiber::toFiber(fibs.pick_one(RNG));
        if ( fib == 0 )
            throw InvalidParameter("Unexpected Fiber specification in single::attach");
        
        real abs = 0;
        if ( opt.set(abs, "attach", 1) )
        {
            FiberEnd ref = ORIGIN;
            opt.set(ref, "attach", KeyList<FiberEnd>("plus_end", PLUS_END, "minus_end", MINUS_END, "center", CENTER), 2);
            abs = fib->abscissaFrom(abs, ref);
        }
        else
        {
            // abscissa is set randomly:
            abs = RNG.real_uniform(fib->abscissaM(), fib->abscissaP());
        }
        
        if ( fib->betweenMP(abs) )
            obj->attachTo(fib, abs);
        else
            throw InvalidParameter("single::attach[1] (abscissa) is out of range");
    }

    return res;
}


//------------------------------------------------------------------------------
#pragma mark -


void SingleSet::relinkA(Single * obj)
{
    assert_true( obj->attached() );
    
    fList.pop(obj);
    aList.push_front(obj);
}


void SingleSet::relinkD(Single * obj)
{
    assert_true( obj->attached() );
    
    aList.pop(obj);
    fList.push_front(obj);
}


void SingleSet::link(Object * obj)
{
    assert_true( obj->objset() == 0 );
    assert_true( obj->tag()==Single::TAG || obj->tag()==Wrist::TAG );

    obj->objset(this);
#ifndef TRAP_SINGLES
    if ( static_cast<Single*>(obj)->attached() )
        aList.push_front(obj);
    else
        fList.push_front(obj);
#else
    Single * s = static_cast<Single*>(obj);
    if (s->trapped())
        trapList.push_front(obj);
    else if ( s->attached() )
        aList.push_front(obj);
    else
        fList.push_front(obj);

#endif
}


void SingleSet::unlink(Object * obj)
{
    assert_true( obj->objset() == this );
    
    obj->objset(0);
#ifndef TRAP_SINGLES
    if ( static_cast<Single*>(obj)->attached() )
        aList.pop(obj);
    else
        fList.pop(obj);
#else

    Single * s = static_cast<Single*>(obj);
    if (s->trapped())
        trapList.pop(obj);
    else if ( s->attached() )
        aList.pop(obj);
    else
        fList.pop(obj);
#endif
    
}


void SingleSet::foldPosition(const Modulo * s) const
{
    Single * obj;
    for ( obj=firstF(); obj; obj=obj->next() )  obj->foldPosition(s);
    for ( obj=firstA(); obj; obj=obj->next() )  obj->foldPosition(s);
#ifdef TRAP_SINGLES
    for ( obj=firstTrapped(); obj; obj=obj->next() )  obj->foldPosition(s);
#endif
}


void SingleSet::mix()
{
    aList.mix(RNG);
    fList.mix(RNG);
#ifdef TRAP_SINGLES
    // This is unnecesary, but maybe we want to use this at some point?
    trapList.mix(RNG);
#endif
}


void SingleSet::erase()
{
    relax();
    ObjectSet::erase(fList);
    ObjectSet::erase(aList);
#ifdef TRAP_SINGLES
    ObjectSet::erase(trapList);
#endif
    inventory.clear();
}


void SingleSet::freeze()
{
    relax();
    ObjectSet::flag(aList, 1);
    ObjectSet::flag(fList, 1);
#ifdef TRAP_SINGLES
    ObjectSet::flag(trapList,1);
#endif
}


void SingleSet::prune()
{
    ObjectSet::prune(aList, 1);
    ObjectSet::prune(fList, 1);
#ifdef TRAP_SINGLES
    ObjectSet::prune(trapList,1);
#endif
}


void SingleSet::thaw()
{
    ObjectSet::flag(aList, 0);
    ObjectSet::flag(fList, 0);
#ifdef TRAP_SINGLES
    ObjectSet::flag(trapList,0);
#endif
}


void SingleSet::write(Outputter & out) const
{
    if ( sizeA() > 0 )
    {
        out.put_line("#section single A");
        ObjectSet::write(aList, out);
    }
    if ( sizeF() > 0 )
    {
        out.put_line("#section single F");
        ObjectSet::write(fList, out);
    }
#ifdef TRAP_SINGLES
    if (sizeTrapped()>0)
    {
        out.put_line("#section single Trapped");
        ObjectSet::write(trapList, out);
    }
#endif
}


int SingleSet::bad() const
{
    int code = 0;
    Single * ghi;
    code = fList.bad();
    if ( code ) return code;
    for ( ghi = firstF(); ghi ; ghi=ghi->next() )
    {
        if ( ghi->attached() ) return 100;
    }
    
    code = aList.bad();
    if ( code ) return code;
    for ( ghi = firstA();  ghi ; ghi=ghi->next() )
    {
        if ( !ghi->attached() ) return 101;
    }

    
    return code;
}

//------------------------------------------------------------------------------
#pragma mark -


void SingleSet::report(std::ostream& os) const
{
    if ( size() > 0 )
    {
        os << title() << "\n";
        PropertyList plist = simul.properties.find_all(title());
        for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
        {
            SingleProp * sp = static_cast<SingleProp*>(*ip);
            ObjectList olist = collect(match_property, sp);
            os << std::setw(10) << olist.size() << " " << sp->name();
            os << " ( " << sp->hand << " )\n";
        }
        if ( plist.size() > 1 )
            os << std::setw(10) << size() << " total\n";
    }
}


ObjectList SingleSet::collect() const
{
    ObjectList res = ObjectSet::collect(fList);
    res.append( ObjectSet::collect(aList) );
#ifdef TRAP_SINGLES
    res.append(ObjectSet::collect(trapList));
#endif
    return res;
}


ObjectList SingleSet::collect(bool (*func)(Object const*, void const*), void const* arg) const
{
    ObjectList res = ObjectSet::collect(fList, func, arg);
    res.append( ObjectSet::collect(aList, func, arg) );
#ifdef TRAP_SINGLES
    res.append( ObjectSet::collect(trapList, func, arg) );
#endif
    return res;
}


//------------------------------------------------------------------------------
#pragma mark - Wrists

/**
 This will create Wrists with `obj` as Base, following the specifications given in `str`.
 These Wrists will be anchored on points `fip` to `fip+nbp-1` of `obj`.
 
 The syntax understood for `str` is as follows:
 @code 
 [INTEGER] [NAME_OF_SINGLE] [each]
 @endcode
 
 The first optional integer specifies the number of Singles to be attached.
 If 'each' is specified, this number is multiplied by the number of point `nbp`,
 and every point receives the same number of Singles.
 
 This is used to decorate Solid and Sphere
 */
ObjectList SingleSet::makeWrists(Mecable const* obj, unsigned fip, unsigned nbp, std::string& arg)
{
    ObjectList res;
    unsigned num = 1;

    std::istringstream iss(arg);
    iss >> num;
    
    if ( iss.fail() )
    {
        num = 1;
        iss.clear();
    }
    
    if ( num == 0 || nbp == 0 )
        return res;
    
    std::string str, mod;
    iss >> str >> mod;
    
    SingleProp * sip = simul.findProperty<SingleProp*>("single", str);
    if ( sip == 0 )
        throw InvalidParameter("could not find fiber:attach single `"+str+"'");

    if ( mod == "each" )
    {
        for ( unsigned u = 0; u < num; ++u )
        {
            for ( unsigned i = 0; i < nbp; ++i )
                res.push_back(sip->newWrist(obj, fip+i));
        }
    }
    else
    {
        for ( unsigned u = 0; u < num; ++u )
        {
            unsigned i = RNG.pint(nbp);
            res.push_back(sip->newWrist(obj, fip+i));
        }
    }
    
    return res;
}


ObjectList SingleSet::collectWrists(Object * foot) const
{
    ObjectList res;
    
    for ( Single * s=firstF(); s; s=s->next() )
        if ( s->base() == foot )
            res.push_back(s);
    
    for ( Single * s=firstA(); s; s=s->next() )
        if ( s->base() == foot )
            res.push_back(s);
    
    return res;
}


void SingleSet::removeWrists(Object * foot)
{
    Single *nxt, *obj;
    
    obj = firstF();
    while ( obj )
    {
        nxt = obj->next();
        if ( obj->base() == foot )
            delete(obj);
        obj = nxt;
    }

    obj = firstA();
    while ( obj )
    {
        nxt = obj->next();
        if ( obj->base() == foot )
            delete(obj);
        obj = nxt;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Fast Diffusion


/**
 Distribute Singles on the sites specified in `loc`.
 */
void SingleSet::uniAttach(Array<FiberBinder>& loc, SingleList& reserve)
{
    Single * s = reserve.front();
    
    for ( Array<FiberBinder>::iterator i = loc.begin(); i < loc.end(); ++i )
    {
        Hand const* h = s->hand();
        
        if ( h->attachmentAllowed(*i) )
        {
            Vector pos = i->pos();
            
            Space const* spc = i->fiber()->prop->confine_space_ptr;

            // Only attach if position is within the confining Space:
            if ( spc && spc->outside(pos) )
                continue;

            if ( s->prop->fast_diffusion == 3 )
            {
                if ( ! spc )
                    continue;
                // Only attach if position is near the edge of the Space:
                Vector prj;
                spc->project(pos, prj);
                if ( distanceSqr(pos, prj) > h->prop->binding_range_sqr )
                    continue;
                // Single will be placed on the edge of the Space:
                pos = prj;
            }
            else
            {
                /*
                 place the Single in the line perpendicular to the attachment point,
                 at a random distance within the range of attachment of the Hand:
                 */
                pos += i->dirFiber().randOrthoB(h->prop->binding_range);
            }

            reserve.pop_front();
            s->setPosition(pos);
            s->attach(*i);
            link(s);
            
            if ( reserve.empty() )
                return;
            s = reserve.front();
        }
    }
}



/**
 Implements a Monte-Carlo approach for attachments of free Single, assumming that
 diffusion is sufficiently fast to maintain a uniform spatial distribution,
 and that the distribution of fibers is more-or-less uniform such that the
 attachments are distributed randomly along the fibers.
 
 Diffusing (free) Single are removed from the standard list, and thus the
 random walk that is used for simulating diffusion will be skipped,
 as well as the detection of neighboring fibers done for attachments.
 
 Algorithm:
 - Remove diffusing Single from the simulation, transfering them to a 'reserve'.
 - Estimate the distance between binding sites occuring in one time-step, from:
    - the total length of fibers,
    - the volume of the Space,
    - the binding parameters of the relevant Hand.
    .
 - Attach Singles from the reserve, at random positions along the Fibers
 .
 
 Note: there is a similar feature for Couple
 */
void SingleSet::uniAttach(FiberSet const& fibers)
{
    // transfer free Single that fast-diffuse to the reserve
    Single * obj = firstF(), * nxt;
    while ( obj )
    {
        nxt = obj->next();
        SingleProp const* p = static_cast<SingleProp const*>(obj->property());
        if ( p->fast_diffusion )
        {
            unlink(obj);
            assert_true(!obj->attached());
            assert_true((size_t)p->index() < uniLists.size());
            uniLists[p->index()].push_front(obj);
        }
        obj = nxt;
    }
    
    Array<FiberBinder> loc(1024);
    
    // uniform attachment for the reserves:
    for ( SingleReserve::iterator cr = uniLists.begin(); cr < uniLists.end(); ++cr )
    {
        SingleList & reserve = *cr;
        if ( !reserve.empty() )
        {
            Single * obj = reserve.front();
            SingleProp const* p = static_cast<SingleProp const*>(obj->property());
            
            const real vol = p->spaceVolume();
            const int cnt = reserve.size();
            
            if ( p->fast_diffusion == 2 )
            {
                real dis = vol / ( cnt * p->hand_prop->bindingSection(0) );
                fibers.newFiberSitesP(loc, dis);
            }
            else
            {
                real dis = vol / ( cnt * p->hand_prop->bindingSection(1) );
                fibers.uniFiberSites(loc, dis);
            }
            
            uniAttach(loc, reserve);
        }
    }
}


/**
 
 Return true if at least one single:fast_diffusion is true,
 and in this case allocate uniLists.
 
 The Volume of the Space is assumed to remain constant until the next uniPrepare()
 */
bool SingleSet::uniPrepare(PropertyList const& properties)
{
    unsigned inx = 0;
    bool res = false;
    
    PropertyList plist = properties.find_all("single");
    
    for ( PropertyList::const_iterator n = plist.begin(); n != plist.end(); ++n )
    {
        SingleProp const* p = static_cast<SingleProp const*>(*n);
        if ( p->fast_diffusion )
            res = true;
        
        if ( p->index() > inx )
            inx = p->index();
    }
    
    if ( res )
        uniLists.resize(inx+1);
    
    return res;
}


/**
 empty uniLists, reversing all Singles in the normal lists.
 This is useful if ( single:fast_diffusion == true )
 */
void SingleSet::uniRelax()
{
    for ( SingleReserve::iterator res = uniLists.begin(); res != uniLists.end(); ++res )
    {
        SingleList& reserve = *res;
        while( ! reserve.empty() )
        {
            Single * s = reserve.front();
            assert_true(!s->attached());
            s->randomizePosition();
            reserve.pop_front();
            link(s);
        }
    }
}
#ifdef TRAP_SINGLES

void SingleSet::trap_single(Single * s)
{
    
    if ( s->attached() )
        aList.pop(s);
    else
        fList.pop(s);
    trapList.push_front(s);
}

void SingleSet::untrap_single(Single * s)
{
    
    trapList.pop(s);
    if ( s->attached() )
        aList.push_front(s);
    else
        fList.push_front(s);
}


#if (TRAP_SINGLES==2)
void SingleSet::populateTrapList()
{
    // Here only the free ones should be included
    Single * obj = firstF(), * nxt;
    while ( obj )
    {
        assert_false(obj->trapped());
        nxt = obj->next();
        simulTrapReserve.push_back(obj);
        obj = nxt;
    }
    
}

void SingleSet::checkTrapped() const
{
    Single * obj = firstTrapped(), * nxt;
    while ( obj )
    {
        nxt = obj->next();
        if (!obj->trapped()) {
            std::cout<<"Untrapped single found in trapped"<<std::endl;
        }
        obj = nxt;
    }
}

#endif
#endif



// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "couple_set.h"
#include "couple_prop.h"
#include "fork_prop.h"
#include "crosslink_prop.h"
#include "shackle_prop.h"
#include "bridge_prop.h"
#include "duo_prop.h"
#include "glossary.h"
#include "simul.h"
#include "swapper_prop.h"
#include "flipper_prop.h"
#ifdef TRAP_SINGLES
#include "trapper_prop.h"
#endif




/**
 @defgroup CoupleGroup Couple and Derived Activities
 @ingroup ObjectGroup
 @ingroup NewObject
 @brief A Couple contains two Hand, and can thus crosslink two Fibers.

 The plain Couple may crosslink two Fiber irrespective of their configuration.
 Derived classes implement specificity, angular stiffness, etc.
 
 List of classes accessible by specifying `couple:activity`.

 `activity`    |   Classes               | Parameters           |  Property
 --------------|-------------------------|----------------------|------------
 `diffuse`     | Couple CoupleLong       | @ref CouplePar       | CoupleProp
 `crosslink`   | Crosslink CrosslinkLong | @ref CrosslinkPar    | CrosslinkProp
 `bridge`      | Bridge                  | @ref BridgePar       | BridgeProp
 `duo`         | Duo  DuoLong            | @ref DuoPar          | DuoProp
 `slide`       | Shackle ShackleLong     | @ref ShacklePar      | ShackleProp
 `fork`        | Fork                    | @ref ForkPar         | ForkProp

 Example:
 @code
 set couple complex
 {
   hand1 = kinesin
   hand2 = kinesin
   stiffness = 100
   diffusion = 10
   activity = crosslink
   length = 0.02
 }
 @endcode
 */

Property* CoupleSet::newProperty(const std::string& kd, const std::string& nm, Glossary& opt) const
{
    CoupleProp * sp = 0;
    if ( kd == "couple" )
    {
        std::string a;
        if ( opt.peek(a, "activity") )
        {
            if ( a == "fork" )
                sp = new ForkProp(nm);
            else if ( a == "crosslink" )
                sp = new CrosslinkProp(nm);
            else if ( a == "bridge" )
                sp = new BridgeProp(nm);
            else if (a == "flip")
                sp = new FlipperProp(nm);
#ifdef TRAP_SINGLES
            else if (a == "trap")
                sp = new TrapperProp(nm);
#endif
            else if (a == "swap")
                sp = new SwapperProp(nm);
            else if ( a == "duo" )
                sp = new DuoProp(nm);
            else if ( a == "slide" )
                sp = new ShackleProp(nm);
            else if ( a == "diffuse" )
                sp = new CoupleProp(nm);
            else 
                throw InvalidParameter("unknown single:activity `"+a+"'");
        }
        if ( sp == 0 )
            sp = new CoupleProp(nm);
    }
    return sp;
}

//------------------------------------------------------------------------------

void CoupleSet::prepare(PropertyList const& properties)
{
    uni = uniPrepare(properties);
}


void CoupleSet::step(FiberSet const& fibers, FiberGrid const& fgrid)
{
    // use alternative attachment strategy:
    if ( uni )
    {
        uniCollect();
        uniAttach(fibers);
    }

    /*
     ATTENTION: We ensure here that step() is called exactly once for each object.
     The Couples are stored in multiple lists, and are automatically transfered
     from one list to another one if their Hands bind or unbind.
     The code relies on the fact that a Couple will be moved to the start of the
     list to which it is transfered. By proceeding always from the node, which was
     first before any transfer could occur, we process each Couple only once.
     Moreover, we get the 'next' in the list always before calling 'step()', because
     'step()' may transfer the node to another list, changing the value of 'next()'
     */
    
    /*
     MSG(9,"CoupleSet::step : FF %5i AF %5i FA %5i AA %5i\n",
         ffList.size(), afList.size(), faList.size(), bList.size());
     */
    
    Couple *const ffHead = firstFF();
    Couple *const afHead = firstAF();
    Couple *const faHead = firstFA();
    
    bool const faOdd = faList.size() % 2;
    bool const afOdd = afList.size() % 2;
    bool const ffOdd = ffList.size() % 2;

    Couple * obj, * nxt;
    
    obj = firstAA();
    // this loop is unrolled, processing objects 2 by 2:
    if ( aaList.size() % 2 )
    {
        nxt = obj->next();
        obj->stepAA();
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        obj->stepAA();
        obj = nxt->next();
        nxt->stepAA();
    }
    
    obj = faHead;
    // this loop is unrolled, processing objects 2 by 2:
    if ( faOdd )
    {
        nxt = obj->next();
        obj->stepFA(fgrid);
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        obj->stepFA(fgrid);
        obj = nxt->next();
        nxt->stepFA(fgrid);
    }

    obj = afHead;
    // this loop is unrolled, processing objects 2 by 2:
    if ( afOdd )
    {
        nxt = obj->next();
        obj->stepAF(fgrid);
        obj = nxt;
    }
    while ( obj )
    {

        nxt = obj->next();
        obj->stepAF(fgrid);
        obj = nxt->next();
        nxt->stepAF(fgrid);
    }
    
    //std::clog << "CoupleSet::step : FF " << ffList.size() << " head " << ffHead << std::endl;

    obj = ffHead;
    // this loop is unrolled, processing objects 2 by 2:
    if ( ffOdd )
    {
        nxt = obj->next();
        obj->stepFF(fgrid);
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        obj->stepFF(fgrid);
        obj = nxt->next();
        nxt->stepFF(fgrid);
    }
    // added by manu
#ifdef  GILLESPIE_DIGITS
    gilles_step(fibers, fgrid);
#endif
}


//------------------------------------------------------------------------------
#pragma mark -

Object * CoupleSet::newObjectT(const Tag tag, unsigned idx)
{
    if ( tag == Couple::TAG )
    {
        CoupleProp * p = simul.findProperty<CoupleProp*>("couple", idx);
        if ( p == 0 )
            throw InvalidIO("no couple class defined with id "+sMath::repr(idx));
        return p->newCouple();
    }
    return 0;
}


/**
 @addtogroup CoupleGroup

 You can attach the hands of a Couple:
 @code
 new couple protein
 {
    attach1 = FIBER, REAL, REFERENCE
    attach2 = FIBER, REAL, REFERENCE
 }
 @endcode
 
 where:
 - FIBER designates the fiber:
     - `fiber1` of `fiber2` correspond to fibers directly
     - `first` or `last` to the oldest and youngest fiber
     - `last-1` the penultimate, etc.
     .
 - REAL is the abscissa of the attachment point.
   If the abscissa is not specified, and random position along
   along the fiber will be selected.
 - REFERENCE can be `minus_end`, `center` or `plus_end` (default = `origin`).
   This defines from which position the abscissa is measured.
 .
 
 */
ObjectList CoupleSet::newObjects(const std::string& name, Glossary& opt)
{
    Property * p = simul.properties.find_or_die(name);
    Couple * obj = static_cast<CoupleProp*>(p)->newCouple(&opt);
    
    ObjectList res;
    res.push_back(obj);
        
    std::string spe;
    /*
     This provides a way for the user to attach hand1:
     */
    if ( opt.set(spe, "attach1") )
    {
        ObjectList fibs = simul.fibers.findObjects(spe);
        if ( fibs.empty() )
            throw InvalidParameter("Could not find Fiber in couple::attach1");
        
        Fiber* fib = Fiber::toFiber(fibs.pick_one(RNG));
        if ( fib == 0 )
            throw InvalidParameter("Unexpected Fiber specification in couple::attach1");
        
        real abs = 0;
        if ( opt.set(abs, "attach1", 1) )
        {
            FiberEnd ref = ORIGIN;
            opt.set(ref, "attach1", KeyList<FiberEnd>("plus_end", PLUS_END, "minus_end", MINUS_END, "center", CENTER), 2);
            abs = fib->abscissaFrom(abs, ref);
        }
        else
        {
            // abscissa is set randomly:
            abs = RNG.real_uniform(fib->abscissaM(), fib->abscissaP());
        }
        
        if ( fib->betweenMP(abs) )
            obj->attachTo1(fib, abs);
        else
            throw InvalidParameter("couple:attach1[1] (abscissa) is out of range");
    }
    
    /*
     This provides a way for the user to attach hand2:
     */
    if ( opt.set(spe, "attach2") )
    {
        ObjectList fibs = simul.fibers.findObjects(spe);
        if ( fibs.empty() )
            throw InvalidParameter("Could not find Fiber in couple::attach2");
        
        Fiber* fib = Fiber::toFiber(fibs.pick_one(RNG));
        if ( fib == 0 )
            throw InvalidParameter("Unexpected Fiber specification in couple::attach2");
        
        real abs = 0;
        if ( opt.set(abs, "attach2", 1) )
        {
            FiberEnd ref = ORIGIN;
            opt.set(ref, "attach2", KeyList<FiberEnd>("plus_end", PLUS_END, "minus_end", MINUS_END, "center", CENTER), 2);
            abs = fib->abscissaFrom(abs, ref);
        }
        else
        {
            // abscissa is set randomly:
            abs = RNG.real_uniform(fib->abscissaM(), fib->abscissaP());
        }
        
        if ( fib->betweenMP(abs) )
            obj->attachTo2(fib, abs);
        else
            throw InvalidParameter("couple:attach2[1] (abscissa) is out of range");
    }
    
    return res;
}

//------------------------------------------------------------------------------
#pragma mark -


void CoupleSet::relinkA1(Couple * obj)
{
    assert_true( obj->attached1() );

    if ( obj->attached2() )
    {
        faList.pop(obj);
        aaList.push_front(obj);
    }
    else
    {
        ffList.pop(obj);
        afList.push_front(obj);
    }
}


void CoupleSet::relinkD1(Couple * obj)
{
    assert_true( obj->attached1() );
    
    if ( obj->attached2() )
    {
        aaList.pop(obj);
        faList.push_front(obj);
    }
    else
    {
        afList.pop(obj);
        ffList.push_front(obj);
    }
}


void CoupleSet::relinkA2(Couple * obj)
{
    assert_true( obj->attached2() );

    if ( obj->attached1() )
    {
        afList.pop(obj);
        aaList.push_front(obj);
    }
    else
    {
        ffList.pop(obj);
        faList.push_front(obj);
    }
}


void CoupleSet::relinkD2(Couple * obj)
{
    assert_true( obj->attached2() );

    if ( obj->attached1() )
    {
        aaList.pop(obj);
        afList.push_front(obj);
    }
    else
    {
        faList.pop(obj);
        ffList.push_front(obj);
    }
}


void CoupleSet::link(Object * obj)
{
    assert_true( obj->objset() == 0 );
    assert_true( obj->tag() == Couple::TAG );
    
    obj->objset(this);
    
    if ( static_cast<Couple*>(obj)->attached1() )
    {
        if ( static_cast<Couple*>(obj)->attached2() )
            aaList.push_front(obj);
        else
            afList.push_front(obj);
    }
    else
    {
        if ( static_cast<Couple*>(obj)->attached2() )
            faList.push_front(obj);
        else
            ffList.push_front(obj);
    }
}


void CoupleSet::unlink(Object * obj)
{
    assert_true( obj->objset() == this );
    
    obj->objset(0);

    if ( static_cast<Couple*>(obj)->attached1() )
    {
        if ( static_cast<Couple*>(obj)->attached2() )
            aaList.pop(obj);
        else
            afList.pop(obj);
    }
    else
    {
        if ( static_cast<Couple*>(obj)->attached2() )
            faList.pop(obj);
        else
            ffList.pop(obj);
    }
}



void CoupleSet::foldPosition(const Modulo * s) const
{
    Couple * cx;
    for ( cx=firstAA(); cx; cx=cx->next() )  cx->foldPosition(s);
    for ( cx=firstFA(); cx; cx=cx->next() )  cx->foldPosition(s);
    for ( cx=firstAF(); cx; cx=cx->next() )  cx->foldPosition(s);
    for ( cx=firstFF(); cx; cx=cx->next() )  cx->foldPosition(s);
}


void CoupleSet::mix()
{
    ffList.mix(RNG);
    afList.mix(RNG);
    faList.mix(RNG);
    aaList.mix(RNG);
}


void CoupleSet::erase()
{
    relax();
    ObjectSet::erase(aaList);
    ObjectSet::erase(faList);
    ObjectSet::erase(afList);
    ObjectSet::erase(ffList);
    inventory.clear();
}

void CoupleSet::freeze()
{
    relax();
    ObjectSet::flag(aaList, 1);
    ObjectSet::flag(faList, 1);
    ObjectSet::flag(afList, 1);
    ObjectSet::flag(ffList, 1);
}


void CoupleSet::prune()
{
    ObjectSet::prune(aaList, 1);
    ObjectSet::prune(faList, 1);
    ObjectSet::prune(afList, 1);
    ObjectSet::prune(ffList, 1);
}


void CoupleSet::thaw()
{
    ObjectSet::flag(aaList, 0);
    ObjectSet::flag(faList, 0);
    ObjectSet::flag(afList, 0);
    ObjectSet::flag(ffList, 0);
}


void CoupleSet::write(Outputter & out) const
{
    if ( sizeAA() > 0 )
    {
        out.put_line("#section couple AA");
        ObjectSet::write(aaList, out);
    }
    if ( sizeAF() > 0 )
    {
        out.put_line("#section couple AF");
        ObjectSet::write(afList, out);
    }
    if ( sizeFA() > 0 )
    {
        out.put_line("#section couple FA");
        ObjectSet::write(faList, out);
    }
    if ( sizeFF() > 0 )
    {
        out.put_line("#section couple FF");
        ObjectSet::write(ffList, out);
    }
}


void CoupleSet::report(std::ostream& os) const
{
    if ( size() > 0 )
    {
        os << title() << "\n";
        PropertyList plist = simul.properties.find_all(title());
        for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
        {
            CoupleProp * cp = static_cast<CoupleProp*>(*ip);
            ObjectList olist = collect(match_property, cp);
            os << std::setw(10) << olist.size() << " " << cp->name();
            os << " ( " << cp->hand1 << " | " << cp->hand2 << " )\n";
        }
        if ( plist.size() > 1 )
            os << std::setw(10) << size() << " total\n";
    }
}


ObjectList CoupleSet::collect() const
{
    ObjectList res = ObjectSet::collect(ffList);
    res.append( ObjectSet::collect(afList) );
    res.append( ObjectSet::collect(faList) );
    res.append( ObjectSet::collect(aaList) );
    return res;
}


ObjectList CoupleSet::collect(bool (*func)(Object const*, void const*), void const* arg) const
{
    ObjectList res = ObjectSet::collect(ffList, func, arg);
    res.append( ObjectSet::collect(afList, func, arg) );
    res.append( ObjectSet::collect(faList, func, arg) );
    res.append( ObjectSet::collect(aaList, func, arg) );
    return res;
}



int CoupleSet::bad() const
{
    int code = 0;
    Couple * obj;
    code = ffList.bad();
    if ( code ) return 100+code;
    for ( obj=firstFF(); obj ; obj = obj->next() )
    {
        if ( obj->attached1() || obj->attached2() )
            return 100;
    }
    
    code = afList.bad();
    if ( code ) return 200+code;
    for ( obj=firstAF(); obj ; obj = obj->next() )
    {
        if ( !obj->attached1() || obj->attached2() )
            return 200;
    }
    
    code = faList.bad();
    if ( code ) return 300+code;
    for ( obj=firstFA(); obj ; obj = obj->next() )
    {
        if ( obj->attached1() || !obj->attached2() )
            return 300;
    }
    
    code = aaList.bad();
    if ( code ) return 400+code;
    for ( obj=firstAA(); obj ; obj = obj->next() )
    {
        if ( !obj->attached1() || !obj->attached2() )
            return 400;
    }
    return code;
}


//------------------------------------------------------------------------------
#pragma mark - Fast Diffusion


/**
Distribute Hand1 of Couples on the sites specified in `loc`.
 */
void CoupleSet::uniAttach1(Array<FiberBinder>& loc, CoupleList& reserve)
{
    for ( Array<FiberBinder>::iterator i = loc.begin(); i < loc.end(); ++i )
    {
        if ( reserve.empty() )
            return;
        Couple * c = reserve.front();
        if ( c->hand1()->attachmentAllowed(*i) )
        {
            reserve.pop_front();
            c->attach1(*i);
            link(c);
        }
    }
}


/**
 Distribute Hand2 of Couples on the sites specified in `loc`.
 */
void CoupleSet::uniAttach2(Array<FiberBinder>& loc, CoupleList& reserve)
{
    for ( Array<FiberBinder>::iterator i = loc.begin(); i < loc.end(); ++i )
    {
        if ( reserve.empty() )
            return;
        Couple * c = reserve.front();
        if ( c->hand2()->attachmentAllowed(*i) )
        {
            reserve.pop_front();
            c->attach2(*i);
            link(c);
        }
    }
}



/**
 Distribute Couples on crossing points specified in `loc`.
 `loc` contains positions on the fibers corresponding to crossing points
 as returned by FiberSet::allIntersections()
 */
void CoupleSet::uniAttach12(Array<FiberBinder>& loc, CoupleList& reserve, int nb)
{
    int nbc = loc.size() / 2;
    for ( int n = 0; n < nb; ++n )
    {
        if ( reserve.empty() )
            return;
        Couple * c = reserve.front();
        reserve.pop_front();
        int p = RNG.pint(nbc);
        c->attach1(loc[2*p  ]);
        c->attach2(loc[2*p+1]);
        link(c);
    }
}


/**
 Implements a Monte-Carlo approach for attachments of free Couple, assumming that
 diffusion is sufficiently fast to maintain a uniform spatial distribution,
 and that the distribution of fibers is more-or-less uniform such that the
 attachments are distributed randomly along the fibers.
 
 Diffusing (free) Couple are removed from the standard list, and thus the
 random walk that is used for simulating diffusion will be skipped,
 as well as the detection of neighboring fibers done for attachments.
 The attachment of already attached Couple is unchanged.
 
 Algorithm:
 - Remove diffusing Single from the simulation, transfering them to a 'reserve'.
 - Estimate the distance between binding sites occuring in one time-step, from:
    - the total length of fibers,
    - the volume of the Space,
    - the binding parameters of the relevant Hand.
    .
 - Attach Singles from the reserve, at random positions along the Fibers
 .
 
 Note: there is a similar feature for Single
 */
void CoupleSet::uniAttach(FiberSet const& fibers)
{
    // preallocate array:
    Array<FiberBinder> loc(1024);
    
#if ( 0 )
    
    // this performs a basic verification of fibers.uniFiberSites()
    real dis = 1;
    unsigned cnt = 1<<10;
    real avg = 0;
    real var = 0;
    for ( unsigned i = 0; i < cnt; ++i )
    {
        fibers.uniFiberSites(loc, dis);
        real s = loc.size();
        avg += s;
        var += s*s;
    }
    avg /= cnt;
    var = var/cnt - avg * avg;
    printf("UNI-FIBER-SITES(1)  avg = %9.2f   var = %9.2f\n", avg, var);
    
#endif
    
    // uniform attachment for reserved couples:
    for ( CoupleReserve::iterator cr = uniLists.begin(); cr < uniLists.end(); ++cr )
    {
        CoupleList & reserve = *cr;
        
        if ( reserve.empty() )
            continue;
        
        Couple * obj = reserve.front();
        CoupleProp const * p = obj->prop;
        
        const real vol = p->spaceVolume();
        const int cnt = reserve.size();

        if ( p->fast_diffusion == 2 )
        {
            real dis = 2 * vol / ( cnt * p->hand1_prop->bindingSection(0) );
            fibers.newFiberSitesP(loc, dis);
        }
        else
        {
            real dis = 2 * vol / ( cnt * p->hand1_prop->bindingSection(1) );
            fibers.uniFiberSites(loc, dis);
        }
        
        uniAttach1(loc, reserve);
        
        if ( reserve.empty() || p->trans_activated )
            continue;
        
        // if ( couple:trans_activated == true ), Hand2 cannot bind
        if ( p->fast_diffusion == 2 )
        {
            real dis = 2 * vol / ( cnt * p->hand2_prop->bindingSection(0) );
            fibers.newFiberSitesP(loc, dis);
        }
        else
        {
            real dis = 2 * vol / ( cnt * p->hand2_prop->bindingSection(1) );
            fibers.uniFiberSites(loc, dis);
        }
        
        uniAttach2(loc, reserve);
    }
}


/**
 
 Return true if at least one couple:fast_diffusion is true,
 and in this case allocate uniLists.
 
 The Volume of the Space is assumed to remain constant until the next uniPrepare() 
 */
bool CoupleSet::uniPrepare(PropertyList const& properties)
{
    unsigned inx = 0;
    bool res = false;
    
    PropertyList plist = properties.find_all("couple");
    
    for ( PropertyList::const_iterator n = plist.begin(); n != plist.end(); ++n )
    {
        CoupleProp const * p = static_cast<CoupleProp const*>(*n);
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
 Transfer free complex that fast-diffuse to the reserve lists
*/
void CoupleSet::uniCollect()
{
    Couple * obj = firstFF(), * nxt;
    while ( obj )
    {
        nxt = obj->next();
        CoupleProp const* p = static_cast<CoupleProp const*>(obj->property());
        if ( p->fast_diffusion )
        {
            unlink(obj);
            assert_true((size_t)p->index() < uniLists.size());
            uniLists[p->index()].push_back(obj);
        }
        obj = nxt;
    }
}


/**
 empty uniLists, reversing all Couples in the normal lists.
 This is useful if ( couple:fast_diffusion == true )
 */
void CoupleSet::uniRelax()
{
    for ( CoupleReserve::iterator res = uniLists.begin(); res != uniLists.end(); ++res )
    {
        CoupleList& reserve = *res;
        while( ! reserve.empty() )
        {
            Couple * c = reserve.front();
            assert_true(!c->attached1() && !c->attached2());
            c->randomizePosition();
            link(c);
            reserve.pop_front();
        }
    }
}


//------------------------------------------------------------------------------
# pragma mark - Equilibrate


void CoupleSet::equilibrateSym(FiberSet const& fibers, CoupleList& reserve, CoupleProp const* cop)
{
    if ( cop->hand1_prop != cop->hand2_prop )
        throw InvalidParameter("Cannot equilibrate heterogeneous Couple");
    
    if ( cop->trans_activated )
        throw InvalidParameter("Cannot equilibrate trans_activated Couple");

    const real space_volume = cop->spaceVolume();
    const real total_length = fibers.totalLength();

    const real binding_rate = cop->hand1_prop->binding_rate;
    const real binding_range = cop->hand1_prop->binding_range;
    const real unbinding_rate = cop->hand1_prop->unbinding_rate;
    
    // get all crosspoints:
    Array<FiberBinder> loc(1024);
    fibers.allIntersections(loc, binding_range);
    const unsigned nb_crossings = loc.size() / 2;
    //const real nb_crossings = sqr(total_length) / ( M_PI * space_volume );

    const real ratio_fibs = 2 * total_length * binding_range / space_volume;
    const real ratio_cros = 4 * M_PI * nb_crossings * sqr(binding_range) / space_volume;
    
    real bind = binding_rate / unbinding_rate;
    real BsG = bind / 2;
    real AsF = ( ratio_fibs - ratio_cros ) * bind;
    real GsF = ratio_cros * bind;
    
    real popF = reserve.size() / ( 1 + AsF + GsF + BsG * GsF );
    real popA = AsF * popF;
    real popG = GsF * popF;
    real popB = BsG * popG;
    
    if ( 0 )
    {
        printf("Couple::equilibrate %s (sym):\n", cop->name().c_str());
        printf("     total %lu\n", reserve.size());
        const real nb_fiber = fibers.size();
        const real fiber_length = total_length / nb_fiber;
        const real nbc = nb_fiber * ( nb_fiber - 1 ) * sqr(fiber_length) / ( M_PI * space_volume );
        //const real nbc = sqr(total_length) / ( M_PI * space_volume );
        printf("     nb_crossings predicted  %9.2f   true %9i\n", nbc, nb_crossings);
        printf("     F %9.2f A %9.2f G %9.2f B %9.2f\n", popF, popA, popG, popB);
    }
    
    // create doubly-attached Couples at the crossing positions:
    uniAttach12(loc, reserve, RNG.poisson(popB));
    
    real dis = 2 * total_length / ( popA + popG );

    if ( !reserve.empty() )
    {
        fibers.uniFiberSites(loc, dis);
        uniAttach1(loc, reserve);
    }
    
    if ( !reserve.empty() )
    {
        fibers.uniFiberSites(loc, dis);
        uniAttach2(loc, reserve);
    }
}


/**
 This attempts to create a configuration of Couple that is close to the equilibrium
 that would be reached, after a sufficient time is given for binding and unbinding.
 This assumes that the configuration of filaments does not change, and also that
 it is random, in particular without bundles. The motion of the motor is also ignored.
 */
void CoupleSet::equilibrate(FiberSet const& fibers, CoupleList& reserve, CoupleProp const* cop)
{
    if ( cop->trans_activated )
        throw InvalidParameter("Cannot equilibrate trans_activated Couple");
    
    const real space_volume = cop->spaceVolume();
    const real total_length = fibers.totalLength();

    const real binding_rate1 = cop->hand1_prop->binding_rate;
    const real binding_range1 = cop->hand1_prop->binding_range;
    const real unbinding_rate1 = cop->hand1_prop->unbinding_rate;

    const real binding_rate2 = cop->hand2_prop->binding_rate;
    const real binding_range2 = cop->hand2_prop->binding_range;
    const real unbinding_rate2 = cop->hand2_prop->unbinding_rate;

    
    // get all crosspoints:
    Array<FiberBinder> loc(1024);
    fibers.allIntersections(loc, std::max(binding_range1, binding_range2));
    const unsigned nb_crossings = loc.size() / 2;
    
    const real ratio_fibs1 = 2 * total_length * binding_range1 / space_volume;
    const real ratio_fibs2 = 2 * total_length * binding_range2 / space_volume;
    const real ratio_cros1 = 4 * M_PI * nb_crossings * sqr(binding_range1) / space_volume;
    const real ratio_cros2 = 4 * M_PI * nb_crossings * sqr(binding_range2) / space_volume;
    
    real BsG1 = binding_rate1 / unbinding_rate1;
    real BsG2 = binding_rate2 / unbinding_rate2;
    real A1sF = ( ratio_fibs1 - ratio_cros1 ) * BsG1 / 2;
    real A2sF = ( ratio_fibs2 - ratio_cros2 ) * BsG2 / 2;
    real G1sF = ratio_cros1 * BsG1 / 2;
    real G2sF = ratio_cros2 * BsG2 / 2;
    
    // the two should be equal
    real BsF = 0.5 * ( BsG1 * G1sF + BsG2 * G2sF );

    real popF = reserve.size() / ( 1.0 + A1sF + A2sF + G1sF + G2sF + BsF );
    real popA1 = A1sF * popF;
    real popA2 = A2sF * popF;
    real popG1 = G1sF * popF;
    real popG2 = G2sF * popF;
    real popB = BsF * popF;

    if ( 0 )
    {
        printf("Couple::equilibrate %s:\n", cop->name().c_str());
        printf("     total %lu\n", reserve.size());
        const real nb_fiber = fibers.size();
        const real fiber_length = total_length / nb_fiber;
        const real nbc = nb_fiber * ( nb_fiber - 1 ) * sqr(fiber_length) / ( M_PI * space_volume );
        //const real nbc = sqr(total_length) / ( M_PI * space_volume );
        printf("     nb_crossings predicted  %9.2f   true %9i\n", nbc, nb_crossings);
        printf("     F %9.2f A %9.2f G %9.2f B %9.2f\n", popF, popA1+popA2, popG1+popG2, popB);
    }

    // create doubly-attached Couples at the crossing positions:
    uniAttach12(loc, reserve, RNG.poisson(popB));
    
    if ( !reserve.empty() )
    {
        const real dis = total_length / ( popA1 + popG1 );
        fibers.uniFiberSites(loc, dis);
        uniAttach1(loc, reserve);
    }
    
    if ( !reserve.empty() )
    {
        const real dis = total_length / ( popA2 + popG2 );
        fibers.uniFiberSites(loc, dis);
        uniAttach2(loc, reserve);
    }
}


void CoupleSet::equilibrate(FiberSet const& fibers, PropertyList const& properties)
{
    PropertyList plist = properties.find_all("couple");

    for ( PropertyList::const_iterator n = plist.begin(); n != plist.end(); ++n )
    {
        CoupleProp const* cop = static_cast<CoupleProp const*>(*n);
        
        if ( !cop->trans_activated )
        {
            CoupleList list;
            
            // collect all Couple of this kind:
            Couple * c = firstFF(), * nxt;
            while ( c )
            {
                nxt = c->next();
                if ( c->property() == cop )
                {
                    unlink(c);
                    list.push_back(c);
                }
                c = nxt;
            }
            if ( list.size() > 0 )
            {
                equilibrate(fibers, list, cop);
                
                // release all collected Couple
                while ( !list.empty() )
                {
                    link(list.front());
                    list.pop_front();
                }
            }
        }
    }
    //printf("Couple::equilibrate    FF %u FA %u AF %u AA %u\n", sizeFF(), sizeFA(), sizeAF(), sizeAA());
}


/**
 This takes all the Free Couple and attach them at the intersection points of the network of filaments
 */
void CoupleSet::connect(FiberSet const& fibers, PropertyList const& properties)
{
    PropertyList plist = properties.find_all("couple");
    
    if ( plist.empty() )
        return;

    real range = 0;
    for ( PropertyList::const_iterator n = plist.begin(); n != plist.end(); ++n )
    {
        CoupleProp const* cop = static_cast<CoupleProp const*>(*n);
        range = std::max(range, cop->hand1_prop->binding_range);
        range = std::max(range, cop->hand2_prop->binding_range);
    }
    
    if ( range <= 0 )
        throw InvalidParameter("hand:binding_range should be > 0!");
    
    // get all crossings:
    Array<FiberBinder> loc(1024);
    fibers.allIntersections(loc, range);
    const unsigned nb_crossings = loc.size() / 2;
    
    //std::clog << nb_crossings << " intersections at range " << range << "\n";

    if ( nb_crossings > 0 )
    {
        Couple * c = firstFF(), * nxt;
        while ( c )
        {
            nxt = c->next();
            int p = RNG.pint(nb_crossings);
            c->attach1(loc[2*p  ]);
            c->attach2(loc[2*p+1]);
            c = nxt;
        }
    }
}



#pragma mark - Gillespie steps


typedef std::set<Couple *, bool(*)( Couple * , Couple *)> GillesList;

bool comp_gillest(Couple * a, Couple * b)
{
    return a->gilles_t < b->gilles_t;
}


// class called objectlist-> container class can be used for couple pointers


void separate_gilles(const NodeList & ori, GillesList & gilles, const real max_time)
{
    
    if (ori.size())
    {
        Couple * obj = static_cast<Couple*>(ori.front());
        Couple * nxt;
        while ( obj )
        {
            nxt = obj->next();
            if (obj->prop->is_gillespie)
            {
                // Calculate the firing time:
                obj->gilles_t = obj->new_gilles_t();
                // If the firing time is bigger than the dt, no point in including it
                if (obj->gilles_t < max_time)
                {
                    gilles.insert(obj);
                }
            }
            obj = nxt;
        }
    }
    return;
}

void CoupleSet::gilles_step(FiberSet const& fibers, FiberGrid const& fgrid)
{
    const real max_time = simul.prop->time_step;
    //Added by manu
    // The problem here would be that the sum of the rates would be small compared to the dt, causing the number of steps to be zero. However, this is not the case where you would you use the gillespie, so in principle it should not matter.
    // store the couples which should do a gillespie step in these lists, then run the gillespie_step
    GillesList sorted_list(&comp_gillest);
    GillesList::iterator it;

    Couple * obj;
    real t;

    // Collect all
    separate_gilles(afList, sorted_list, max_time);
    separate_gilles(faList, sorted_list, max_time);
    separate_gilles(aaList, sorted_list, max_time);

    while (!sorted_list.empty())
    {
        // The first element is the one that should react first
        it = sorted_list.begin();
        obj = *it;
        sorted_list.erase(it);
        t = obj->gilles_t;
        // If this is true, the couple has detached by reaching the end of the filament
        if (!obj->gillestep())
        {
            //std::clog<<t << " fire "<<obj->identity()<<std::endl;
            obj->gilles_t = t + obj->new_gilles_t();
            // add back only if the new time of gillespie is smaller than the time step
            if (obj->gilles_t < max_time)
                sorted_list.insert(obj);
        }
    }
    
    return;
}

#ifdef TRAP_SINGLES
#if (TRAP_SINGLES==2)
void CoupleSet::populateTrapList()
{
    // Any couple can be trapped, whether they are attached or not
    Couple * obj = firstFF(), * nxt;
    while ( obj )
    {
        nxt = obj->next();
        if (obj->prop->activity=="trap")
        {
            if (obj->trapped())
                simulUntrapReserve.push_back(obj);
            else
                simulTrapReserve.push_back(obj);
        }
        obj = nxt;
    }
    
    obj = firstAF();
    while ( obj )
    {
        nxt = obj->next();
        if (obj->prop->activity=="trap" && !obj->trapped())
            simulTrapReserve.push_back(obj);
        obj = nxt;
    }
    
    obj = firstFA();
    while ( obj )
    {
        nxt = obj->next();
        if (obj->prop->activity=="trap" && !obj->trapped())
            simulTrapReserve.push_back(obj);
        obj = nxt;
    }

    obj = firstAA();
    while ( obj )
    {
        nxt = obj->next();
        if (obj->prop->activity=="trap" && !obj->trapped())
            simulTrapReserve.push_back(obj);
        obj = nxt;
    }
}

void CoupleSet::trap(HaMonList * singles)
{
    if (!simulTrapReserve.size() || !singles->size())
        return;
    HandMonitor * coup = simulTrapReserve.front();
    real trap_rate = coup->trap_rate();
    int nb_trappers = simulTrapReserve.size();
    int nb_singles =singles->size();
    // Given the rate, get a poissonian number
    int nb_trappings = RNG.poisson(trap_rate*nb_trappers*nb_singles);
    // The number of trapping events cannot be bigger than the available molecules
    if (nb_trappings>nb_singles)
        nb_trappings = nb_singles;
    if (nb_trappings>nb_trappers)
        nb_trappings = nb_trappers;
    // In principle, the list is randomized, so we just go along the list.
    HandMonitor * couple_i = 0, * single_i = 0;
    for (int i=0; i<nb_trappings; i++)
    {
        couple_i = simulTrapReserve.front();
        single_i = singles->front();
        assert_false(couple_i->trapped());
        couple_i->trap(single_i);
        
        simulTrapReserve.pop_front();
        singles->pop_front();
        
    }


}

void CoupleSet::untrap()
{
    if (!simulUntrapReserve.size())
        return;
    HandMonitor * coup = simulUntrapReserve.front();
    int nb_trapped = simulUntrapReserve.size();
    real untrap_rate = coup->untrap_rate();
    int nb_untrapping = RNG.poisson(nb_trapped*untrap_rate);
    
    if (nb_untrapping>nb_trapped)
        nb_untrapping = nb_trapped;
    HandMonitor * couple_i = 0;
    for (int i=0; i<nb_untrapping; i++)
    {
        couple_i = simulUntrapReserve.front();
        couple_i->untrap();
        simulUntrapReserve.pop_front();
    }
}
#endif
#endif


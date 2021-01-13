// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "sim.h"
#include "simul.h"
#include "meca.h"
#include "exceptions.h"
#include "hand_prop.h"
#include "simul_prop.h"
#include "backtrace.h"
#include "modulo.h"

extern Modulo const* modulo;

#include "simul_step.cc"
#include "simul_file.cc"
#include "simul_custom.cc"
#include "simul_report.cc"
#include "simul_solve.cc"
#include "simul_solvex.cc"
#include "simul_solvef.cc"

#include "nucleus.h"
#include "aster.h"
#include "bundle.h"
#include "fake.h"
#include "wrist.h"
#include "space_strip.h"
#include "space_periodic.h"
#include "space_cylinderP.h"
#include "fiber.h"
#include "event.h"

#include <csignal>

//---------------------------  global variables/functions ---------------------

bool functionKey[17] = { 0 };


void out_of_memory_handler()
{
    fprintf(stderr, "* * * * *\n");
    fprintf(stderr, "Cytosim: memory allocation failed\n");
    fprintf(stderr, "* * * * *\n");
    print_backtrace(stderr);
    std::exit(1);
}

void termination_handler()
{
    fprintf(stderr, "* * * * *\n");
    fprintf(stderr, "Cytosim: uncaught exception\n");
    fprintf(stderr, "* * * * *\n");
    print_backtrace(stderr);
    std::abort();
}

void fpe_handler(int sig)
{
    std::cout << std::endl;
    fprintf(stderr, "* * * * *\n");
    psignal(sig, "Cytosim");
    fprintf(stderr, "* * * * *\n");
    print_backtrace(stderr);
    std::exit(sig);
}

//------------------------------------------------------------------------------
#pragma mark -

Simul::Simul()
: prop(0), spaces(*this), fields(*this),
fibers(*this), spheres(*this), beads(*this), solids(*this),
singles(*this), couples(*this), organizers(*this), events(*this)
{
    sSpace        = 0;
    precond       = 0;
    precondCPU[0] = 0;
    precondCPU[1] = 0;
    precondCPU[2] = 0;
    precondCPU[3] = 0;
    precondCounter = 0;
    
    prop = new SimulProp("undefined");
}

Simul::~Simul()
{
    erase();
    
    if ( prop )
    {
        delete(prop);
        prop = 0;
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 set current Space to `spc`. (spc==NULL is a valid argument).
*/
void Simul::changeSpace(Space const* spc)
{
    if ( spc != sSpace )
    {
        sSpace = spc;
        
#if ( 0 )
        if ( spc )
            std::clog << "Simul::changeSpace(" << spc->prop->name() << ")" << std::endl;
        else
            std::clog << "Simul::changeSpace(NULL)" << std::endl;
#endif
    }
    
    modulo = 0;

    if ( sSpace )
    {
        assert_true(spc->prop);

        sSpace->setModulo(&sModulo);
        
        if ( sModulo.isPeriodic() )
            modulo = &sModulo;
    }
}


const Space* Simul::findSpace(std::string const& str) const
{
    if ( str == "first" )
        return static_cast<Space*>(spaces.inventory.first());

    if ( str == "last" )
        return static_cast<Space*>(spaces.inventory.last());
    
    Property * sp = properties.find("space", str);
    
    if ( sp )
        return spaces.first(static_cast<SpaceProp*>(sp));
    else
        return 0;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 This will initialize the simulation by registering callbacks.
 You should still call Simul::prepare() before calling step()
 */
void Simul::initialize(Glossary & glos)
{
    // Register a function to be called if operator new fails:
    std::set_new_handler(out_of_memory_handler);
    
    // Register a function to be called upon abortion:
    std::set_terminate(termination_handler);
    
    // Register a function to be called for Floating point exceptions:
    if ( signal(SIGFPE, fpe_handler) == SIG_ERR )
        std::cerr << "Could not register SIGFPE handler\n";
    
    // read parameters, and complete
    prop->read(glos);
}

//------------------------------------------------------------------------------

real Simul::time() const
{
    return prop->time;
}

void Simul::erase()
{
    relax();
    organizers.erase();
    fibers.erase();
    spheres.erase();
    beads.erase();
    solids.erase();
    singles.erase();
    couples.erase();
    fields.erase();
    spaces.erase();
    events.erase();
    
    // destroy all properties, except the SimulProp:
    properties.erase();
 
    prop->time = 0;
    sSpace     = 0;
    modulo     = 0;
}


unsigned Simul::nbObjects() const
{
    return  (  organizers.size()
             + singles.size()
             + couples.size()
             + fibers.size()
             + beads.size()
             + solids.size()
             + spheres.size()
             + spaces.size()
             + fields.size() );
}


void Simul::foldPosition() const
{
    if ( modulo )
    {
        fibers.foldPosition(modulo);
        beads.foldPosition(modulo);
        solids.foldPosition(modulo);
        spheres.foldPosition(modulo);
        singles.foldPosition(modulo);
        couples.foldPosition(modulo);
        organizers.foldPosition(modulo);
    }
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Convert pointer to Mecable*
 Note: We could have a virtual function in Object:
 Mecable * Object::toMecable() { return 0; }
 */
Mecable* Simul::toMecable(Object * obj)
{
    if ( obj )
    switch( obj->tag() )
    {
        case  Fiber::TAG:  return static_cast<Mecable*>(obj);
        case   Bead::TAG:  return static_cast<Mecable*>(obj);
        case  Solid::TAG:  return static_cast<Mecable*>(obj);
        case Sphere::TAG:  return static_cast<Mecable*>(obj);
    }
    return 0;
}

/**
 Find an object from a Human-friendly representation, such as
 fiber1
 single1
 */
Object * Simul::findObject(const std::string& arg) const
{
    // extract name part:
    size_t n = 0;
    while ( arg.size() > n && isalpha(arg[n]) )
        ++n;
    // extract serial number:
    long num = strtol(arg.c_str()+n, 0, 10);
    std::string str = arg.substr(0, n);
    // find ObjectSet:
    ObjectSet const* set = findSet(str);
    if ( set == 0 )
        throw InvalidIO("unknown class `"+str+"'");
    return set->findObject(num);
}


void Simul::add(Object * w)
{
    //std::clog << " Simul::add(" << w->reference() << ")" << std::endl;
    assert_true(w);
    ObjectSet * set = findSetT(w->tag());
    set->add(w);
}


void Simul::add(ObjectList objs)
{
    //std::clog << " Simul::add("<< objs.size() <<" objects):" << std::endl;
    for ( Object ** oi = objs.begin(); oi < objs.end(); ++oi )
    {
        Object * obj = *oi;
        if ( obj )
            add(obj);
    }
}


void Simul::remove(Object * w)
{
    assert_true( w->objset() );
    w->objset()->remove(w);
}


void Simul::remove(ObjectList objs)
{
    //std::clog << " Simul::remove("<< objs.size() <<" objects):" << std::endl;
    for ( Object ** oi = objs.begin(); oi < objs.end(); ++oi )
    {
        Object * obj = *oi;
        if ( obj )
            remove(obj);
    }
}


void Simul::erase(Object * w)
{
    //std::clog << "Simul::erase " << w->reference() << std::endl;
    remove(w);
    delete(w);
}


void Simul::erase(ObjectList objs)
{
    //std::clog << " Simul::erase(" << objs.size() << " objects):\n";
    for ( Object ** oi = objs.begin(); oi < objs.end(); ++oi )
    {
        Object * obj = *oi;
        if ( obj )
        {
            //std::clog << " Simul::erase(" << obj << ")\n";
            remove(obj);
            delete(obj);
        }
    }
}


void Simul::mark(ObjectList objs, int mrk)
{
    //std::clog << " Simul::erase("<< objs.size() <<" objects):" << std::endl;
    for ( Object ** oi = objs.begin(); oi < objs.end(); ++oi )
        (*oi)->mark(mrk);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 This is used primarily to parse the configuration file,
 using full class name
 */
ObjectSet * Simul::findSet(const std::string& kind)
{
    //std::clog << "findSet("<<kind<<")"<<std::endl;
    if ( kind == "space" )        return &spaces;
    if ( kind == "field" )        return &fields;
    if ( kind == "fiber" )        return &fibers;
    if ( kind == "bead" )         return &beads;
    if ( kind == "solid" )        return &solids;
    if ( kind == "sphere" )       return &spheres;
    if ( kind == "single" )       return &singles;
    if ( kind == "couple" )       return &couples;
    if ( kind == "aster" )        return &organizers;
    if ( kind == "bundle" )       return &organizers;
    if ( kind == "nucleus" )      return &organizers;
    if ( kind == "fake" )         return &organizers;
    if ( kind == "event" )        return &events;
#ifdef BACKWARD_COMPATIBILITY
    if ( kind == "tubule" )       return &fibers;
#endif
    return 0;
}


/**
 This is used primarily to read the binary trajectory file,
 using a single character to refer to each class in Cytosim
 */
ObjectSet * Simul::findSetT(const Tag tag)
{
    switch( tag )
    {
        case        Couple::TAG:    return &couples;
        case        Single::TAG:    return &singles;
        case         Wrist::TAG:    return &singles;
        case         Fiber::TAG:    return &fibers;
        case Fiber::TAG_DYNAMIC:    return &fibers;
        case Fiber::TAG_LATTICE:    return &fibers;
        case          Bead::TAG:    return &beads;
        case         Solid::TAG:    return &solids;
        case        Sphere::TAG:    return &spheres;
        case         Field::TAG:    return &fields;
        case         Space::TAG:    return &spaces;
        case       Nucleus::TAG:    return &organizers;
        case        Bundle::TAG:    return &organizers;
        case         Aster::TAG:    return &organizers;
        case          Fake::TAG:    return &organizers;
        case         Event::TAG:    return &events;
#ifdef BACKWARD_COMPATIBILITY
        case 'm':                   return &fibers;
#endif
    }
    return 0;
}

//------------------------------------------------------------------------------
#pragma mark -

bool Simul::isPropertyClass(const std::string& kind) const
{
    if ( kind == "simul" )
        return true;
    
    if ( kind == "hand" )
        return true;
    
    return ( 0 != findSet(kind) );
}


Property* Simul::findProperty(const std::string& kd, const std::string& nm) const
{
    if ( kd == "simul" )
        return prop;

    return properties.find(kd, nm);
}


Property* Simul::findProperty(const std::string& nm) const
{
    if ( nm == prop->name() )
        return prop;
    
    return properties.find(nm);
}


PropertyList Simul::findAllProperties(const std::string& kd) const
{
    if ( kd == "simul" )
    {
        PropertyList list;
        list.push_back(prop);
        return list;
    }
    
    return properties.find_all(kd);
}


/**
 @defgroup ObjectGroup List of objects
 
 The command `set simul` will define the global parameters.
 The `simul` is automatically created, and you cannot use 'new simul'.

 Objects       | Base Class    | Parameters     
 --------------|---------------|----------------
 `simul`       |  Simul        | @ref SimulPar  
 
 
 These objects cannot move:
 
 Class Name    | Base Class    | Parameters       |  Specialization
 --------------|---------------|------------------|-------------------
 `space`       |  Space        | @ref SpacePar    | @ref SpaceGroup
 `field`       |  FieldBase    | @ref FieldPar    | -                 
 `event`       |  Event        | @ref EventPar    | -
 
 
 `Mecables` can move or deform, and come in 4 basic forms:
 
 Class Name    | Base Class    | Parameters       |  Specialization
 --------------|---------------|------------------|-------------------
 `fiber`       |  Fiber        | @ref FiberPar    | @ref FiberGroup
 `bead`        |  Bead         | @ref BeadPar     | -
 `solid`       |  Solid        | @ref BeadPar     | -
 `sphere`      |  Sphere       | @ref SpherePar   | -

 
 A `Hand` is an object that can bind to fiber, but it can only be used
 as a sub-part of `Single` or `Couple`.

 Class Name    | Base Class    | Parameters       |  Specialization
 --------------|---------------|------------------|-------------------
 `hand`        |  Hand         | @ref HandPar     | @ref HandGroup
 
 
 `Single` and `Couple` contain one or two `Hand` respectively:

 Class Name    | Base Class    | Parameters       |  Specialization
 --------------|---------------|------------------|-------------------
 `single`      |  Single       | @ref SinglePar   | @ref SingleGroup
 `couple`      |  Couple       | @ref CouplePar   | @ref CoupleGroup
 
 
 
 The `Organizers` describe composite objects build from multiple Mecables:
 
 Organizers    | Base Class    | Parameters       
 --------------|---------------|------------------
 `aster`       |  Aster        | @ref AsterPar    
 `bundle`      |  Bundle       | @ref BundlePar   
 `nucleus`     |  Nucleus      | @ref NucleusPar  
 `fake`        |  Fake         | @ref FakePar     
 .
 
 */
Property* Simul::newProperty(const std::string& kd, const std::string& nm, Glossary& glos)
{
    if ( kd == "simul" )
    {
        assert_true(prop);
        if ( prop->is_named("undefined") )
            prop->rename(nm);
        return prop;
    }
    
    if ( isPropertyClass(nm) )
        throw InvalidSyntax("`"+nm+"' is a reserved keyword");
    
    Property * p = findProperty(nm);
    
    if ( p )
        throw InvalidSyntax("property `"+nm+"' is already defined");
    
    if ( kd == "hand" )
    {
        p = HandProp::newProperty(nm, glos);
        properties.deposit(p);
    }
    else
    {
        ObjectSet * set = findSet(kd);
        
        if ( set == 0 )
            throw InvalidSyntax("unknown class `"+kd+"'");
        
        p = set->newProperty(kd, nm, glos);
        properties.deposit(p);
    }
    
    return p;
}



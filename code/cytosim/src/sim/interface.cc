// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "interface.h"
#include "glossary.h"
#include "filepath.h"
#include "tictoc.h"
#include <fstream>

extern Random RNG;

// Use the second definition to get some verbose reports:
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG;

//------------------------------------------------------------------------------

Interface::Interface(Simul& s)
: simul(s)
{
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 Property::complete() is called after a property is set.
 This ensures that inconsistencies are detected as early as possible.
 
 The drawback is that we cannot support cross-dependencies (A needs B and vice-versa).
 If that is necessary, we could:
 - call complete() for all Properties, after the parsing process is complete.
 - remove any check for the existence of invoked properties, in which case 
 error would be detected only when objects are created later.
 */
Property* Interface::execute_set(std::string const& kind, std::string const& name, Glossary& def)
{
    VLOG("-SET " << kind << " `" << name << "'\n");
    
    Property* pp = simul.newProperty(kind, name, def);
    
    if ( pp == 0 )
        throw InvalidSyntax("failed to create property of class `"+kind+"'");
    
    pp->read(def);
    pp->complete(&simul);
    
    return pp;
}


Property * Interface::execute_change(std::string const& name, Glossary& def)
{
    // in this form, 'name' designs the parameter
    Property * pp = simul.findProperty(name);
    
    if ( pp == 0 )
        throw InvalidSyntax("unknown property `"+name+"'");
    
    pp->read(def);
    pp->complete(&simul);
    
    return pp;
}


int Interface::execute_change_all(std::string const& kind, Glossary& def)
{
    int cnt = 0;
    PropertyList list = simul.findAllProperties(kind);
    
    if ( list.size() == 0 )
        throw InvalidSyntax("could not find any property `"+kind+"'");
    
    for ( PropertyList::iterator n = list.begin(); n != list.end(); ++n )
    {
        ++cnt;
        (*n)->read(def);
        (*n)->complete(&simul);
    }
    
    return cnt;
}


//------------------------------------------------------------------------------
#pragma mark -
#include <stream_func.h>

Isometry Interface::read_position(Glossary& opt)
{
    Isometry res;
    std::string str;
    
    Space const* spc = simul.space();
    
    // Space specified as second argument to 'position'
    if ( opt.set(str, "position", 1) )
        spc = simul.findSpace(str);
    
    // Position
    res.vec.zero();
    if ( opt.set(str, "position") )
    {
        std::istringstream iss(str);
        res.vec = Movable::readPosition(iss, spc);

        if ( Tokenizer::get_symbol(iss) == "to" )
        {
            int mx = 0;
            real alpha = 0;
            if ( opt.set(mx, "instance") && mx > 0 && opt.set(alpha, "instance", 1) )
                alpha /= mx;
            else
                alpha = RNG.preal();
            Vector vec = Movable::readPrimitive(iss, spc);
            res.vec += ( vec - res.vec ) * alpha;
        }
    }
    else if ( spc )
        res.vec = spc->randomPlace();
    
    // Rotation applied before the translation
    if ( opt.set(str, "orientation") )
    {
        std::istringstream iss(str);
        res.rot = Movable::readRotation(iss, res.vec, spc);
    }
    
    else if ( opt.set(str, "direction") )
    {
        std::istringstream iss(str);
        Vector vec = Movable::readDirection(iss, res.vec, spc);
        res.rot = Rotation::rotationToVector(vec, RNG);
    }
    else
        res.rot = Rotation::randomRotation(RNG);
    
    // Second rotation applied after the translation
    if ( opt.set(str, "orientation", 1) )
    {
        std::istringstream iss(str);
        Rotation rot = Movable::readRotation(iss, res.vec, spc);
        res.rotate(rot);
    }
    
    return res;
}


enum PlacementType { NONE, ANYWHERE, INSIDE, EDGE, OUTSIDE, ALL_INSIDE };


/**
 @code
 new INTEGER CLASS NAME
 {
   position = POSITION
   placement = PLACEMENT, SPACE_NAME
   nb_trials = INTEGER
 }
 @endcode
 
 
 PLACEMENT can be:
 - if placement = `inside` (default), it tries to find a place inside the Space
 - if placement = `anywhere`, the position is returned
 - if placement = `outside`, the object is created only if it is outside the Space
 - if placement = `surface`, the position is projected on the edge of current Space
 .
 
 By default, the specifications are relative to the last Space that was defined,
 but a different space can be specified as second argument of PLACEMENT.
 
 You can set the density of objects by setting `nb_trials=1`:
 @code
 new 100 single grafted
 {
   position = ( rectangle 10 10 )
   nb_trials = 1
 }
 @endcode
 In this way an object will be created only if its randomly chosen position falls
 inside the Space, and the density will be exactly what is specified from the 
 `position` range (here 100/10*10 = 1).
 */
Isometry Interface::get_placement(Glossary& opt, int placement)
{
    int n = 0, nb_trials = 1<<13;
    std::string str;
    
    opt.set(nb_trials, "nb_trials");

    const Space* spc = simul.space();
    
    if ( opt.set(str, "placement", 1) )
        spc = simul.findSpace(str);

    Isometry iso = read_position(opt);

    while ( ++n < nb_trials )
    {
        if ( spc == 0 || placement == ANYWHERE )
            return iso;
        
        if ( placement == EDGE )
        {
            iso.vec = spc->project(iso.vec);
            return iso;
        }
        
        if ( spc->inside(iso.vec) )
        {
            if ( placement == INSIDE || placement == ALL_INSIDE )
                return iso;
        }
        else
        {
            if ( placement == OUTSIDE )
                return iso;
        }
        
        // generate a new position:
        iso = read_position(opt);
    }
    
    MSG.warning("placement failed\n");
    iso.zero();
    return iso;
}



/**
 This would normally create one object of class 'kind' and type 'name'.
 */
ObjectList Interface::execute_new(std::string const& name, Glossary& opt)
{
    Property * pp = simul.properties.find(name);
    
    if ( pp == 0 )
        throw InvalidSyntax("unknown property `"+name+"'");
    
    ObjectSet * set = simul.findSet(pp->category());

    if ( set == 0 )
        throw InvalidSyntax("could not determine class of `"+name+"'");

    ObjectList res;
    
    do {
        
        res = set->newObjects(name, opt);
        
        //if ( res.count(0) ) std::clog << "cytosim: empty slots in newObjects()\n";
        
        res.remove_pack(0);
        
        PlacementType placement = INSIDE;
        
        opt.set(placement, "placement",
                KeyList<PlacementType>("none",       NONE,
                                       "anywhere",   ANYWHERE,
                                       "inside",     INSIDE,
                                       "all_inside", ALL_INSIDE,
                                       "outside",    OUTSIDE,
                                       "surface",    EDGE));
        
        if ( placement != NONE )
        {
            Isometry iso = get_placement(opt, placement);

            ObjectSet::moveObjects(res, iso);
            
            // special case for which we check all model points:
            if ( placement == ALL_INSIDE )
            {
                for ( Object ** oi = res.begin(); oi < res.end(); ++oi )
                {
                    if ( (*oi)->tag() == Fiber::TAG || (*oi)->tag() == Solid::TAG )
                    {
                        PointSet * ps = static_cast<PointSet*>(*oi);
                        if ( !ps->allInside(simul.space()) )
                        {
                            res.destroy();
                            break;
                        }
                    }
                }
            }
        }
        
    } while ( res.empty() );
    
    // optionally mark the objects:
    int mk = 0;
    if ( opt.set(mk, "mark") )
    {
        for ( Object ** oi = res.begin(); oi < res.end(); ++oi )
            (*oi)->mark(mk);
    }
    
    /* 
     Because the objects in ObjectList are not necessarily all of the same class,
     we call simul.add() rather than directly set->add()
     */
    simul.add(res);
    
    //hold();

    VLOG("-NEW `" << name << "' made " << res.size() << " objects\n");
    
    return res;
}


//------------------------------------------------------------------------------
/**
 Creates `cnt` objects of class `kind` and type `name`.
 The objects are placed at random position in a random orientation within the current Space.
 
 This is meant to be faster than calling execute_new(set, kind, name, opt) 
 `cnt` times.
 */
void Interface::execute_new(std::string const& name, unsigned cnt)
{
    Property * pp = simul.properties.find(name);
    
    if ( pp == 0 )
        throw InvalidSyntax("unknown property `"+name+"'");
    
    ObjectSet * set = simul.findSet(pp->category());
    
    if ( set == 0 )
        throw InvalidSyntax("could not determine class of `"+name+"'");

    Glossary opt;

    for ( unsigned n = 0; n < cnt; ++n )
    {
        ObjectList objs = set->newObjects(name, opt);
        
        if ( simul.space() )
        {
            Isometry iso(simul.space()->randomPlace(), Rotation::randomRotation(RNG));
            ObjectSet::moveObjects(objs, iso);
        }
        
        /* 
         Because the objects in ObjectList are not necessarily all of the same class,
         we call simul.add() rather than directly set->add()
         */
        simul.add(objs);
    }
    
    VLOG("-NEW " << cnt << "`" << name << "' objects\n");

    //hold();
}

//------------------------------------------------------------------------------
#pragma mark -

/// holds a set of criteria used to select Objects
class SelectionCriteria
{
public:

    unsigned     mrk;
    int          st;
    int          st1;
    int          st2;
    void  const* prp;
    Space const* ins;
    Space const* ous;

    /// initialize
    SelectionCriteria()
    {
        mrk = 0;
        st  = -1;
        st1 = -1;
        st2 = -1;
        prp = 0;
        ins = 0;
        ous = 0;
    }
    
    void set(Simul& simul, std::string const& name, Glossary& opt)
    {
        if ( name != "*" )
            prp = simul.properties.find_or_die(name);
        
        std::string str;
        if ( opt.set(str, "position") )
        {
            Space const* spc = 0;
            std::string spn;
            if ( opt.set(spn, "position", 1) )
                spc = simul.findSpace(spn);
            else
                spc = simul.space();
            if ( spc == 0 )
                throw InvalidSyntax("unknown Space `"+spn+"'");
            
            if ( str == "inside" )
                ins = spc;
            else if ( str == "outside" )
                ous = spc;
            else
                throw InvalidSyntax("unknown specification `"+str+"'");
        }
        
        opt.set(mrk, "mark");
        opt.set(st1, "state")    || opt.set(st1, "state1") || opt.set(st1, "stateP");
        opt.set(st2, "state", 1) || opt.set(st2, "state2") || opt.set(st2, "stateM");
    }
    
    /// return `true` if given object fulfills all the conditions specified
    bool match(Object const* obj) const
    {
        if ( mrk > 0 && obj->mark() != mrk )
            return false;
        if ( ins && ins->outside(obj->position()) )
            return false;
        if ( ous && ous->inside(obj->position()) )
            return false;
        if ( prp && obj->property() != prp )
            return false;
        if ( st1 >= 0 )
        {
            if ( obj->tag()==Single::TAG && static_cast<Single const*>(obj)->attached() != st1 )
                return false;
            if ( obj->tag()==Couple::TAG && static_cast<Couple const*>(obj)->attached1() != st1 )
                return false;
            if ( obj->tag()==Fiber::TAG && static_cast<Fiber const*>(obj)->dynamicStateP() != st1 )
                return false;
        }
        if ( st2 >= 0 )
        {
            if ( obj->tag()==Single::TAG )
                throw InvalidParameter("to select Single, 'state[1]' is irrelevant");
            if ( obj->tag()==Couple::TAG && static_cast<Couple const*>(obj)->attached2() != st2 )
                return false;
            if ( obj->tag()==Fiber::TAG && static_cast<Fiber const*>(obj)->dynamicStateM() != st2 )
                return false;
        }
        return true;
    }
};


bool match_criteria(Object const* obj, void const* val)
{
    return static_cast<SelectionCriteria const*>(val)->match(obj);
}


void Interface::execute_delete(std::string const& name, Glossary& opt)
{
    Property * pp = simul.properties.find(name);
    
    if ( pp == 0 )
        throw InvalidSyntax("unknown property `"+name+"'");
    
    ObjectSet * set = simul.findSet(pp->category());
    
    if ( set == 0 )
        throw InvalidSyntax("could not determine class of `"+name+"'");
    
    SelectionCriteria cri;
    cri.set(simul, name, opt);

    ObjectList objs = set->collect(match_criteria, &cri);
    
    //if `opt:nb_objects` is specified, limit the list to a random subset
    unsigned cnt;
    if ( opt.set(cnt, "nb_objects")  &&  cnt < objs.size() )
    {
        objs.mix(RNG);
        objs.truncate(cnt);
    }
    
    //std::clog << " deleting " << objs.size() << " objects" << std::endl;
    simul.erase(objs);
}



void Interface::execute_mark(std::string const& name, Glossary& opt)
{
    Property * pp = simul.properties.find(name);
    
    if ( pp == 0 )
        throw InvalidSyntax("unknown property `"+name+"'");
        
    ObjectSet * set = simul.findSet(pp->category());
    
    if ( set == 0 )
        throw InvalidSyntax("could not determine class of `"+name+"'");

    int mrk;
    if ( ! opt.set(mrk, "mark") )
        throw InvalidParameter("mark must be specified for command `mark'");
    opt.clear("mark");
    
    SelectionCriteria cri;
    cri.set(simul, name, opt);
    
    ObjectList objs = set->collect(match_criteria, &cri);
    
    //if `opt:nb_objects` is specified, limit the list to a random subset
    unsigned cnt;
    if ( opt.set(cnt, "nb_objects")  &&  cnt < objs.size() )
    {
        objs.mix(RNG);
        objs.truncate(cnt);
    }
    
    simul.mark(objs, mrk);
}



void Interface::execute_cut(std::string const& name, Glossary& opt)
{
    Property * pp = simul.properties.find(name);
    
    if ( pp == 0 )
        throw InvalidSyntax("unknown property `"+name+"'");

    if ( pp->category() != "fiber" )
        throw InvalidSyntax("only `cut fiber' is supported");
    
    SelectionCriteria cri;
    cri.set(simul, name, opt);
    
    Vector n(1,0,0);
    real a = 0;
    
    opt.set(n, "plane");
    opt.set(a, "plane", 1);
    
    int stateP = STATE_RED, stateM = STATE_GREEN;
    opt.set(stateP, "new_end_state");
    opt.set(stateM, "new_end_state", 1);
    
    VLOG("-CUT PLANE (" << n << ").x = " << -a << "\n");
    
    ObjectList objs = simul.fibers.collect(match_criteria, &cri);
    
    simul.fibers.planarCut(objs, n, a, stateP, stateM);
}

//------------------------------------------------------------------------------
#pragma mark -

void reportCPUtime(int frame, real sec)
{
    static clock_t clock;
    static double cum = 0;
    
    static int hour = -1;
    int h = TicToc::hours_today();
    if ( hour != h )
    {
        hour = h;
        char date[26];
        TicToc::get_date(date, sizeof(date));
        MSG("%s\n", date);
    }
    
    char cpu[64];
    clock = TicToc::processor_time(cpu, sizeof(cpu), clock, cum);
    
    MSG("F%-6i  %7.2fs   CPU %s\n", frame, sec, cpu);
    MSG.flush();
}


/**
 Perform simulation steps, and write frames to files.
 Currently, only 'run simul *' is supported.
 
 @code
 run [NB_STEPS] simul *
 {
    nb_steps   = INTEGER
    time_span  = REAL
    solve      = SOLVE_MODE
    event      = RATE, ( CODE )
    nb_frames  = INTEGER, ( CODE )
    prune      = BOOL
    flux_speed = NEGATIVE REAL
 }
 @endcode
 
 The optional specification [NB_STEPS] enables the short syntax:
 @code
 run NB_STEPS simul *
 @endcode
 
 Parameter    | Default |   Description
 -------------|---------|-----------------------------------------------------------
 `nb_steps`   |  1      | number of simulation steps
 `time_span`  |  -      | when specified, `nb_steps` is set to `ceil(time_span/time_step)`
 `solve`      |  `all`  | Define the type of method used for the mechanics
 `event`      |  `none` | custom code executed stochastically with prescribed rate
 `nb_frames`  |  0      | number of states written to trajectory file
 `prune`      |  `true` | Print only parameters that are different from default
 
 
 The parameter `solve` can be used to select alternative mechanical engines.
 The monte-carlo part of the simulation is always done, including
 fiber assembly dynamics, binding/unbinding and diffusion of molecules.
 
 Value        |  Result
 -------------|---------------------------------------------------------------------
 `all`        | The mechanics is solved fully (default).
 `off`        | Objects are immobile.
 `horizontal` | Objects can move in the X-direction. The mechanics is solved partly.
 `flux`       | Fibers are translated at `flux_speed` according to their orientation.
 
 If set, `event` defines an event occuring at a rate specified by the positive real `RATE`.
 The action is defined by CODE, a string enclosed with parenthesis containing cytosim commands.
 This code will be executed at stochastic times with the specified rate.
 
 Example:
 @code
 event = 10, ( new fiber actin { position=(rectangle 1 6); length=0.1; } )
 @endcode
 
 Calling `run` will not output the initial state, but this can be done with a separate command:
 @code
 export objects objects.cmo { append = 0 }
 
 run 1000 simul *
 {
    nb_frames = 10
 }
 @endcode
 */
void Interface::execute_run(bool do_write, Glossary& opt)
{
    unsigned     nb_steps   = 1;
    unsigned     nb_frames  = 0;
    std::string  code;
    int          solve      = 1;
    bool         prune      = true;
    bool         binary     = true;
    real         event_rate = 0;
    std::string  event_code;
    real         flux_speed = 0;

    opt.set(nb_steps,   "nb_steps");
    opt.set(nb_frames,  "nb_frames");
    bool has_code = opt.set(code, "nb_frames", 1);
    opt.set(event_rate, "event");
    opt.set(event_code, "event", 1);
    opt.set(solve,      "solve", KeyList<int>("off",0, "on", 1, "horizontal", 2, "flux", 3));
    opt.set(prune,      "prune");
    opt.set(binary,     "binary");
    
    if ( solve == 3 )
    {
        opt.set(flux_speed, "flux_speed");
        if ( flux_speed > 0 )
            throw InvalidParameter("simul:flux_speed should be <= 0");
        flux_speed *= simul.prop->time_step;
    }
    
    // user can specify the final time in seconds:
    real span;
    if ( opt.set(span, "time_span") || opt.set(span, "time") )
        nb_steps = ceil(span/simul.prop->time_step);
    
    unsigned int  frame = 1;
    real          delta = nb_steps;
    unsigned long check = nb_steps;
    
    VLOG("-RUN START " << nb_steps << '\n');

    if ( nb_frames <= 0 )
        do_write = false;
    
    if ( do_write )
    {
        simul.writeProperties(0, prune);
        if ( simul.prop->clear_trajectory )
        {
            simul.writeObjects(0, binary, false);
            simul.prop->clear_trajectory = false;
        }
        if ( has_code )
        {
            std::istringstream iss(code);
            parse(iss, ", while executing 'run' code");
        }

        delta = real(nb_steps) / real(nb_frames);
        check = (int)delta;
    }
    
    simul.prepare();
    
    // Gillespie countdown timer for next event:
    real etime = RNG.exponential();
    // decrement of Gillespie timer corresponding to one time-step:
    real event_rate_dt = event_rate * simul.prop->time_step;
    
    unsigned n = 0;
    while ( 1 )
    {
        if ( n >= check )
        {
            if ( do_write )
            {
                simul.relax();
                simul.writeObjects(0, binary, true);
                reportCPUtime(frame, simul.time());
                if ( has_code )
                {
                    std::istringstream iss(code);
                    parse(iss, ", while executing `run' code");
                }
            }
            if ( n >= nb_steps )
                break;
            check = (int)( ++frame * delta );
        }
        
        simul.step();
        
        if ( solve == 1 )
            simul.solve();
        else if ( solve == 2 )
            simul.solveX();
        else if ( solve == 3 )
            simul.solveF(flux_speed);

        hold();
        
        etime -= event_rate_dt;
        while ( etime < 0 )
        {
            simul.relax();
            VLOG("-EVENT\n");
            std::istringstream iss(event_code);
            parse(iss, ", while executing `event' code");
            etime += RNG.exponential();
        }
        
        ++n;
    }
    
    simul.relax();
    
    VLOG("-RUN COMPLETED\n");
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Import a simulation snapshot from a trajectory file
 
 The frame to be imported can be specified as an option: `frame=INTEGER`:
 @code
 import objects sim_objects.cmo { frame = 10 }
 @endcode
 
 By default, this will replace the simulation state by the one loaded from file.
 To add the file objects to the simulation without deleting the current world,
 you should specify `append=1`:
 
 @code
 import objects sim_objects.cmo { append = 1 }
 @endcode
 */
void Interface::execute_import(std::string const& file, std::string const& what, Glossary& opt)
{
    ObjectSet * subset = 0;
    
    if ( what != "objects" && what != "*" )
    {
        subset = simul.findSet(what);
        if ( subset == 0 )
            throw InvalidIO("unexpected class specified for import");
    }

    Inputter in(file.c_str(), "rb");
    in.vectorSize(DIM);

    if ( ! in.good() )
        throw InvalidIO("Could not open file `"+file+"'");
    
    bool append = false;
    unsigned cnt = 0, frm = 0;

    opt.set(frm, "frame");
    opt.set(append, "append");

    VLOG("-IMPORT frame " << frm << " from " << file << '\n');

    while ( in.good() )
    {
        if ( append )
        {
            real t = simul.prop->time;
            simul.loadObjects(in, subset);
            simul.prop->time = t;
        }
        else
            simul.reloadObjects(in);
        if ( cnt >= frm )
            break;
        ++cnt;
    }
    
    if ( cnt < frm )
        throw InvalidIO("Could not import requested frame");
    
#if ( 0 )
    //unfinished code to mark imported objects
    int mrk;
    if ( opt.set(mrk, "mark") )
    {
         simul.mark(objs, mrk);
    }
#endif
    
    // set time
    real t;
    if ( opt.set(t, "time") )
        simul.prop->time = t;
}


/**
 see Parser::parse_export
 */
void Interface::execute_export(std::string& file, std::string const& what, Glossary& opt)
{
    bool append = true;
    bool binary = true;
    
    opt.set(append, "append");
    opt.set(binary, "binary");

    VLOG("-EXPORT " << what << " to " << file << '\n');
    
    if ( what == "objects" || what == "*" )
    {
        // a '*' designates the usual file filename for output:
        if ( file == "*" )
            file = simul.prop->trajectory_file;

        simul.writeObjects(file.c_str(), binary, append);
    }
    else if ( what == "properties" )
    {
        // a '*' designates the usual file name for output:
        if ( file == "*" )
            file = simul.prop->property_file;
        
        simul.writeProperties(file.c_str(), false);
    }
    else
        throw InvalidIO("only `objects' or `properties' can be exported");
}


/**
 see Parser::parse_report
 */
void Interface::execute_report(std::string& file, std::string const& what, Glossary& opt)
{
    std::string str;
    VLOG("-WRITE " << what << " to " << file << '\n');
    
    // a '*' designates the C-standard output:
    if ( file == "*" )
    {
        simul.report(std::cout, what, opt);
    }
    else
    {
        bool append = true;
        opt.set(append, "append");
        std::ofstream out(file.c_str(), append ? std::ios_base::app : std::ios_base::out);
        simul.report(out, what, opt);
        out.close();
    }
}



void Interface::execute_call(std::string& str, Glossary& opt)
{
    if ( str == "equilibrate" )
        simul.couples.equilibrate(simul.fibers, simul.properties);
    else if ( str == "connect" )
        simul.couples.connect(simul.fibers, simul.properties);
    else if ( str == "custom0" )
        simul.custom0(opt);
    else if ( str == "custom1" )
        simul.custom1(opt);
    else if ( str == "custom2" )
        simul.custom2(opt);
    else if ( str == "custom3" )
        simul.custom3(opt);
    else if ( str == "custom4" )
        simul.custom4(opt);
    else if ( str == "custom5" )
        simul.custom5(opt);
    else if ( str == "custom6" )
        simul.custom6(opt);
    else if ( str == "custom7" )
        simul.custom7(opt);
    else if ( str == "custom8" )
        simul.custom8(opt);
    else if ( str == "custom9" )
        simul.custom9(opt);
    else
        throw InvalidSyntax("called unknown command");
}



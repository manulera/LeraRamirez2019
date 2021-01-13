// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "glossary.h"
#include "iowrapper.h"
#include "key_list.h"
#include "smath.h"
#include "aster.h"
#include <iostream>
#include <numeric>
#include <set>


/// width of columns in formatted output, in number of characters
int column_width = 10;


/// this macro defines the beginning of a line of comment
#define COM "\n%" << std::setw(column_width-1)

/// this macro defines the beginning of a new line of data
#define LIN  '\n' << std::setw(column_width)

/// this macro define the separator used between values
#define SEP   ' ' << std::setw(column_width-1)


/// pad string by adding white-space on the right up to size 'n * column_width - p'
std::string ljust(std::string const& str, unsigned n, unsigned p=0)
{
    size_t s = n * column_width - p;
    if ( str.size() < s )
        return str + std::string(s-str.size(), ' ');
    else
        return str;
}

/// pad string by adding white-space on the left up to size 'n * column_width - p'
std::string rjust(std::string const& str, unsigned n, unsigned p=0)
{
    size_t s = n * column_width - p;
    if ( str.size() < s )
        return std::string(s-str.size(), ' ') + str;
    else
        return str;
}

/// repeat string DIM times with X, Y, Z endings as appropriate
std::string repeatXYZ(std::string const& str, size_t s = column_width)
{
    std::string res(" " + rjust(str+"X", 1, 1));
#if ( DIM > 1 )
    res += " " + rjust(str+"Y", 1, 1);
#endif
#if ( DIM > 2 )
    res += " " + rjust(str+"Z", 1, 1);
#endif
    return res;
}

/// remove any 's' at the end of the argument
void remove_plural(std::string & str)
{
    if ( str.size() > 2  &&  str.at(str.size()-1) == 's' )
        str.resize(str.size()-1);
}


#pragma mark - Report


/** 
 @copydetails Simul::report0
 */
void Simul::report(std::ostream& out, std::string const& arg, Glossary& opt) const
{
    int p = 4;
    opt.set(p, "precision");
    out.precision(p);
    opt.set(column_width, "column") || opt.set(column_width, "column_width");

    out << "\n% start   " << prop->time;
    out << "\n% time    " << prop->time;
    try {
        report0(out, arg, opt);
        out << "\n% end\n";
    }
    catch( Exception & e )
    {
        out << "\n% error: " << e.what();
        out << "\n% end\n";
        throw;
    }
}


/**
 
 WHAT            |   output
 ----------------|------------------------------------------------------
 `bead`          | Position of beads
 `couple`        | Number and state of couples
 `fiber`         | Length and position of the ends of fibers
 `single`        | Number and state of singles
 `solid`         | Position of center and first point of solids
 `sphere`        | Position of center and first point of spheres
 `organizer`     | Position of the center of asters and other organizers
 `field`         | Total quantity of substance in field and Lattices
 `time`          | Time
 `inventory`     | summary list of objects
 `property`      | All object properties
 `parameter`     | All object properties
 
 WHAT                 |   output
 ---------------------|--------------------------------------------------------------------
 `fiber:age`          | Average age of fibers
 `fiber:length`       | Average length and standard deviation of fibers
 `fiber:distribution` | length distribution of fiber lengths (option: `max` and `interval`)
 `fiber:dynamic`      | Number of fiber classified by PLUS_END Dynamic state
 `fiber:point`        | coordinates of model points of all fibers
 `fiber:speckle`      | coordinates of points randomly distributed along all fibers (option: `interval`)
 `fiber:sample`       | coordinates of points newly distributed along all fibers (option: `interval`)
 `fiber:segment`      | information about lengths of segments, number of kinks
 `fiber:end`          | Positions and dynamic states of all fiber ends
 `fiber:force`        | Position of model points and Forces acting on model points
 `fiber:tension`      | Internal stress along fibers
 `fiber:cluster`      | Clusters made of fibers connected by Couples
 `fiber:energy`       | Report fiber's elastic bending energy
 `fiber:binder`       | Report positions of bridging hands along each fiber

 WHAT                    |   output
 ------------------------|--------------------------------------------------------------------
 `bead:all`              | Position of beads
 `bead:single`           | Number of Beads with no single attached, 1 single attached etc.
 `spindle:indice`        | Two scalar indices that caracterize the organization of fibers
 `spindle:profile`       | Number of right- and left-pointing fiber as a function of position
 `single:all`            | Position and force of singles
 `single:NAME`           | Position and force of singles of class NAME
 `single:position:NAME`  | Position and force of singles of class NAME
 `couple:all`            | Position of couples
 `couple:NAME`           | Position of couples of class NAME
 `couple:hands`          | Composition of couples
 
 */
void Simul::report0(std::ostream& out, std::string const& arg, Glossary& opt) const
{
    std::string who = arg, what, which;
    
    // split the argument string into 3 parts separated by ':':
    std::string::size_type pos = arg.find(':');
    if ( pos != std::string::npos )
    {
        who  = arg.substr(0, pos);
        what = arg.substr(pos+1);
        std::string::size_type pas = what.find(':');
        if ( pas != std::string::npos )
        {
            which = what.substr(pas+1);
            what.resize(pas);
        }
    }
    
    // allow for approximate English:
    remove_plural(who);
    remove_plural(what);

    //std::clog << "report("<< what << "|" << who << "|" << which << ")\n";

    if ( who == "fiber" )
    {
        if ( !which.empty() )
            return reportFiber(out, which);

        if ( what.empty() )
            return reportFiber(out);
        if ( what == "end" )
            return reportFiberEnds(out);
        if ( what == "point" )
            return reportFiberPoints(out);
        if ( what == "moment" )
            return reportFiberMoments(out);
        if ( what == "speckle" )
            return reportFiberSpeckles(out, opt);
        if ( what == "sample" )
            return reportFiberSamples(out, opt);
        if ( what == "segment" )
            return reportFiberSegments(out);
        if ( what == "length" )
            return reportFiberLengths(out);
        if ( what == "distribution" )
            return reportFiberLengthDistribution(out, opt);
        if ( what == "tension" )
            return reportFiberTension(out, opt);
        if ( what == "energy" )
            return reportFiberBendingEnergy(out);
        if ( what == "dynamic" )
            return reportFiberDynamic(out);
        if ( what == "force" )
            return reportFiberForces(out);
        if ( what == "confinement" )
            return reportFiberConfinement(out);
        if ( what == "cluster" )
            return reportClusters(out, 1);
        if ( what == "age" )
            return reportFiberAge(out);
        if ( what == "intersection" )
            return reportFiberIntersections(out, opt);
        if ( what == "binder" )
            return reportFiberBinders(out);
        if ( what == "allbinder" )
            return reportFiberBindersAll(out);
        if ( what == "connector" )
            return reportFiberConnectors(out, opt);
        if (what == "cap")
            return reportFiberCap(out, opt);
        if (what == "overlap")
            return reportFiberOverlap(out, opt);
        if (what == "count")
            return reportFiberCount(out, opt);
        throw InvalidSyntax("I only know fiber: end, point, moment, speckle, sample, segment, dynamic, length, distribution, tension, force, cluster, age, energy, binder");
    }
    if ( who == "bead" )
    {
        if ( what.empty() )
            return reportBeadPosition(out);
        if ( what == "single" )
            return reportBeadSingles(out);
        if ( what == "position" )
            return reportBeadPosition(out);
        throw InvalidSyntax("I only know bead: position, single");
    }
    if ( who == "solid" )
    {
        if ( what.empty() )
            return reportSolid(out);
        throw InvalidSyntax("I only know `solid'");
    }
    if ( who == "space" )
    {
        if ( what.empty() )
            return reportSpace(out);
        throw InvalidSyntax("I only know `space'");
    }
    if ( who == "sphere" )
    {
        if ( what.empty() )
            return reportSphere(out);
        throw InvalidSyntax("I only know `sphere'");
    }
    if ( who == "single" )
    {
        if ( what.empty() )
            return reportSingle(out);
        if ( what == "state" || what == "force" )
        {
            if ( which.empty() )
                return reportSingleState(out);
            else
                return reportSingleState(out, which);
        }
        if ( what == "position" )
            return reportSinglePosition(out, which);
        if (what == "trap")
            return reportSingleTrap(out);
        throw InvalidSyntax("I only know single: state, force, NAME");
    }
    if ( who == "flipper")
    {
        if (what.empty())
            return reportFlipper(out);
        throw InvalidSyntax("I only know flipper");
    }
    
    if ( who == "couple" )
    {
        if ( what.empty() )
            return reportCouple(out);
        else if ( what == "state" )
        {
            if ( which.empty() )
                return reportCoupleState(out);
            else
                return reportCoupleState(out, which);
        }
        else if ( what == "link" )
            return reportCoupleLink(out, which);
        else if ( what == "position" )
            return reportCouplePosition(out, which);
        else if ( what == "force" )
            return reportCoupleForce(out, opt);
        else if ( what == "active" )
            return reportCoupleActive(out, which);
        else if ( what == "hand" )
            return reportCoupleHand(out);
#ifdef TRAP_SINGLES
        else if ( what == "trap")
            return reportCoupleTrap(out);
#endif
        throw InvalidSyntax("I only know couple: all, force, NAME");
    }
    if ( who == "organizer" )
    {
        if ( what.empty() )
            return reportOrganizer(out);
        throw InvalidSyntax("I only know `organizer'");
    }
    if ( who == "aster" )
    {
        if ( what.empty() )
            return reportAster(out);
        throw InvalidSyntax("I only know `aster'");
    }
    if ( who == "field" )
    {
        return reportField(out);
    }
    if ( who == "time" )
    {
        if ( what.empty() )
            return reportTime(out);
        throw InvalidSyntax("I only know `time'");
    }
    if ( who == "inventory" )
    {
        if ( what.empty() )
            return reportInventory(out);
        throw InvalidSyntax("I only know `inventory'");
    }
    if ( who == "property" || who == "parameter" )
    {
        if ( what.empty() )
            return writeProperties(out, false);
        else
        {
            Property * p = findProperty(what);
            if ( p == 0 )
                throw InvalidSyntax("unknown Property `"+what+"'");
            p->write(out);
            return;
        }
    }
    if ( who == "spindle" )
    {
        if ( what == "indice" )
            return reportIndices(out);
        if ( what == "profile" )
            return reportProfile(out);
        throw InvalidSyntax("I only know spindle: indices, profile");
    }
    if ( who == "network" && what == "bridge" )
        return reportNetworkBridges(out, opt);

    if ( who == "ring" )
    {
        return reportRing(out);
    }
    if ( who == "ashbya" )
    {
        return reportAshbya(out);
    }
    if ( who == "custom" )
        return reportCustom(out);

    throw InvalidSyntax("I do not know how to report `"+who+"'");
}

//------------------------------------------------------------------------------
#pragma mark - Fiber

/**
 Export average length and standard-deviation for each class of fiber
 */
void Simul::reportFiberAge(std::ostream& out) const
{
    out << COM << ljust("class", 2, 1) << SEP << "count" << SEP << "avg_age";
    out << SEP << "std_dev" << SEP << "min_age" << SEP << "max_age";
    
    unsigned cnt;
    real avg, dev, mn, mx;
    
    PropertyList plist = properties.find_all("fiber");
    
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        FiberProp * fp = static_cast<FiberProp*>(*ip);
        ObjectList objs = fibers.collect(fp);
        fibers.infoBirthtime(objs, cnt, avg, dev, mn, mx);
        real now = prop->time;
        out << LIN << ljust(fp->name(), 2);
        out << SEP << cnt;
        out << SEP << avg;
        out << SEP << dev;
        out << SEP << now-mx;
        out << SEP << now-mn;
    }
}

/**
 Export average length and standard-deviation for each class of fiber
 */
void Simul::reportFiberLengths(std::ostream& out) const
{
    out << COM << ljust("class", 2, 1) << SEP << "count" << SEP << "avg_len" << SEP << "std_dev";
    out << SEP << "min_len" << SEP << "max_len" << SEP << "total";

    unsigned cnt;
    real avg, dev, mn, mx;
    
    PropertyList plist = properties.find_all("fiber");
    
    std::streamsize p = out.precision();
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        FiberProp * fp = static_cast<FiberProp*>(*ip);
        
        ObjectList objs = fibers.collect(fp);
        fibers.infoLength(objs, cnt, avg, dev, mn, mx);
        
        out << LIN << ljust(fp->name(), 2);
        out << SEP << cnt;
        out.precision(3);
        out << SEP << std::fixed << avg;
        out << SEP << std::fixed << dev;
        out << SEP << std::fixed << mn;
        out << SEP << std::fixed << mx;
        out.precision(1);
        out << SEP << std::fixed << avg*cnt;
    }
    out.precision(p);
}


/**
 Export average length and standard-deviation for each class of fiber
 */
void Simul::reportFiberLengthDistribution(std::ostream& out, Glossary & opt) const
{
    real del = 1;
    unsigned nbin = 32;
    opt.set(del, "interval");
    opt.set(nbin, "bins");
    
    unsigned * cnt = new unsigned[nbin+1];
    
    PropertyList plist = properties.find_all("fiber");
    
    std::streamsize p = out.precision();
    out.precision(4);

    if ( 1 )
    {
        out << "\n% length_distribution (`scale` indicates the center of each bin)";
        out << LIN << ljust("scale", 2);
        for ( unsigned u = 0; u <= nbin; ++u )
            out << " " << std::setw(5) << del * ( u + 0.5 );
    }
    
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        FiberProp * fp = static_cast<FiberProp*>(*ip);
        
        for ( unsigned u = 0; u <= nbin; ++u )
            cnt[u] = 0;
        
        for ( Fiber * obj=fibers.first(); obj; obj=obj->next() )
        {
            if ( obj->prop == fp )
            {
                unsigned u = floor( obj->length() / del );
                if ( u < nbin )
                    ++cnt[u];
                else
                    ++cnt[nbin];
            }
        }

        out << LIN << ljust(fp->name(), 2);
        for ( unsigned u = 0; u <= nbin; ++u )
            out << " " << std::setw(5) << cnt[u];
    }
    out.precision(p);
    delete[] cnt;
}


void Simul::reportFiberSegments(std::ostream& out) const
{
    out << COM << ljust("class", 2, 1) << SEP << "nb_fibers" << SEP << "nb_joints";
    out << SEP << "nb_kinks" << SEP << "min_seg" << SEP << "max_seg";
    
    PropertyList plist = properties.find_all("fiber");
    
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        FiberProp * fp = static_cast<FiberProp*>(*ip);
        
        unsigned cnt, joints, kinks;
        real mn = 0, mx = 0;
        
        ObjectList objs = fibers.collect(fp);
        fibers.infoSegments(objs, cnt, joints, kinks, mn, mx);
        out << LIN << ljust(fp->name(), 2);
        out << SEP << cnt;
        out << SEP << joints;
        out << SEP << kinks;
        out << SEP << mn;
        out << SEP << mx;
    }
}


/**
 Export number of fiber, classified according to dynamic state of one end
 */
void Simul::reportFiberDynamic(std::ostream& out, FiberEnd end) const
{
    const int MAX = 5;
    int cnt[MAX] = { 0 };
    int sum = 0;
    
    for ( Fiber * obj=fibers.first(); obj; obj=obj->next() )
    {
        ++sum;
        unsigned s = obj->dynamicState(end);
        if ( s < MAX )
            ++cnt[s];
    }

    if ( end == PLUS_END )
        out << LIN << ljust("plus_end", 1);
    else if ( end == MINUS_END )
        out << LIN << ljust("minus_end", 1);
 
    out << SEP << sum;
    for ( int ii = 0; ii < MAX; ++ii )
        out << SEP << cnt[ii];
}

/**
 Export number of fiber, classified according to dynamic state of one end
 */
void Simul::reportFiberDynamic(std::ostream& out) const
{
    out << COM << "fiber_end" << SEP << "total" << SEP << "static";
    out << SEP << "red" << SEP << "orange" << SEP << "yellow" << SEP << "green";
    reportFiberDynamic(out, PLUS_END);
    reportFiberDynamic(out, MINUS_END);
}


/**
 Export length, position and directions at center of fibers
 */
void Simul::reportFiber(std::ostream& out, FiberProp const* fp) const
{
    out << COM << "class" << SEP << "identity" << SEP << "length";
    out << repeatXYZ("position") << repeatXYZ("direction");
    out << "end-to-end" << SEP << "cosinus";

    for ( Fiber * obj=fibers.first(); obj; obj=obj->next() )
    {
        if ( obj->prop == fp )
        {
            out << LIN << obj->prop->index();
            out << SEP << obj->identity();
            out << SEP << obj->length();
            out << SEP << obj->posEnd(CENTER);
            out << SEP << obj->dirEnd(CENTER);
            out << SEP << (obj->posEndM()-obj->posEndP()).norm();
            out << SEP << obj->dirEndM() * obj->dirEndP();
        }
    }
}


/**
 Export length, position and directions at center of fibers
 */
void Simul::reportFiber(std::ostream& out, std::string const& which) const
{
    Property * fp = properties.find("fiber", which);
    
    if ( !fp )
        throw InvalidSyntax("unknown fiber class `"+which+"'");

    reportFiber(out, static_cast<FiberProp*>(fp));
}



/**
 Export length, position and directions at center of fibers
 */
void Simul::reportFiber(std::ostream& out) const
{
    PropertyList plist = properties.find_all("fiber");
    
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        FiberProp const* fp = static_cast<FiberProp*>(*ip);
        out << COM << " fiber class " + sMath::repr(fp->index()) + " is " + fp->name();
        reportFiber(out, fp);
    }
}


/**
 Export dynamic state, positions and directions of fiber at both ends
 */
void Simul::reportFiberEnds(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << SEP << "length";
    out << SEP << "stateM";
    for (int i=0;i<DIM;i++)
        out << SEP << "positionM";
    for (int i=0;i<DIM;i++)
        out << SEP << "directM";
    out << SEP << "stateP";
    for (int i=0;i<DIM;i++)
        out << SEP << "positionP";
    for (int i=0;i<DIM;i++)
        out << SEP << "directP";
        
    for ( Fiber * obj=fibers.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->prop->index();
        out << SEP << obj->identity();
        out << SEP << obj->length();
        out << SEP << obj->dynamicStateM();
        out << SEP << obj->posEndM();
        out << SEP << obj->dirEndM();
        out << SEP << obj->dynamicStateP();
        out << SEP << obj->posEndP();
        out << SEP << obj->dirEndP();
    }
    out << std::endl;
}


/**
 Export Fiber-number, position of model points
 */
void Simul::reportFiberPoints(std::ostream& out) const
{
    out << COM << "identity" << repeatXYZ("position");

    // we print Fibers in the order of the inventory:
    Fiber * fib = fibers.firstID();

    while ( fib )
    {
        out << COM << " fiber " << fib->reference();
        
        for ( unsigned p = 0; p < fib->nbPoints(); ++p )
        {
            out << LIN << fib->identity();
            out << SEP << fib->posP(p);
        }
        
        fib = fibers.nextID(fib);
    }
}


/// Helper class to calculate moments of a cloud of points
class Accumulator
{
public:
    real sum;
    real avg[3];
    real var[9];
    
    void reset()
    {
        sum = 0;
        for ( int i = 0; i < 3; ++i ) avg[i] = 0;
        for ( int i = 0; i < 9; ++i ) var[i] = 0;
    }
    
    Accumulator()
    {
        reset();
    }

    void add(real w, Vector const& p)
    {
        sum += w;
        avg[0] += w * p.XX;
        var[0] += w * p.XX * p.XX;
#if ( DIM > 1 )
        avg[1] += w * p.YY;
        var[1] += w * p.YY * p.XX;
        var[4] += w * p.YY * p.YY;
#endif
#if ( DIM > 2 )
        avg[2] += w * p.ZZ;
        var[2] += w * p.ZZ * p.XX;
        var[5] += w * p.ZZ * p.YY;
        var[8] += w * p.ZZ * p.ZZ;
#endif
    }
    
    void add(Vector const& p)
    {
        sum += 1;
        avg[0] += p.XX;
        var[0] += p.XX * p.XX;
#if ( DIM > 1 )
        avg[1] += p.YY;
        var[1] += p.YY * p.XX;
        var[4] += p.YY * p.YY;
#endif
#if ( DIM > 2 )
        avg[2] += p.ZZ;
        var[2] += p.ZZ * p.XX;
        var[5] += p.ZZ * p.YY;
        var[8] += p.ZZ * p.ZZ;
#endif
    }
    
    void subtract_mean()
    {
        //Remove the mean:
        avg[0] /= sum;
        var[0] = var[0]/sum - avg[0] * avg[0];
#if ( DIM > 1 )
        avg[1] /= sum;
        var[1] = var[1]/sum - avg[1] * avg[0];
        var[4] = var[4]/sum - avg[1] * avg[1];
#endif
#if ( DIM > 2 )
        avg[2] /= sum;
        var[2] = var[2]/sum - avg[2] * avg[0];
        var[5] = var[5]/sum - avg[2] * avg[1];
        var[8] = var[8]/sum - avg[2] * avg[2];
#endif
    }
    
    void print_doc(std::ostream& out) const
    {
        out << LIN << "cnt";
        out << SEP << "avgX" << SEP << "avgY" << SEP << "avgZ";
        out << SEP << "varX" << SEP << "varY" << SEP << "varZ";
        out << SEP << "var_sum";
    }
    
    void print(std::ostream& out, bool mode)
    {
        if ( mode )
            out << LIN << (int) sum;
        else
            out << SEP << sum;
        out << SEP << avg[0];
        out << SEP << avg[1];
        out << SEP << avg[2];
        out << SEP << var[0];
        out << SEP << var[4];
        out << SEP << var[8];
        out << SEP << var[0] + var[4] + var[8];
    }
};


/**
 Export first and second-order moments of model points
 */
void Simul::reportFiberMoments(std::ostream& out) const
{
    out << COM << ljust("class", 2, 1) << SEP << "sum";
    out << SEP << "avgX" << SEP << "avgY" << SEP << "avgZ";
    out << SEP << "varX" << SEP << "varY" << SEP << "varZ" << SEP << "var_sum";
    out << std::fixed;
    
    Accumulator accum;
    
    PropertyList plist = properties.find_all("fiber");
    
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        FiberProp * fp = static_cast<FiberProp*>(*ip);
        
        accum.reset();
       
        for ( Fiber * fib=fibers.first(); fib; fib=fib->next() )
        {
            if ( fib->prop == fp )
            {
                const real w = fib->segmentation();
                accum.add(0.5*w, fib->posEndM());
                for ( unsigned n = 1; n < fib->lastPoint(); ++n )
                    accum.add(w, fib->posP(n));
                accum.add(0.5*w, fib->posEndP());
            }
        }
        
        accum.subtract_mean();
        out << LIN << ljust(fp->name(), 2);
        accum.print(out, 0);
    }
}


/**
 Export Fiber-number, position of model points
 */
void Simul::reportFiberForces(std::ostream& out) const
{
    computeForces();

    out << LIN << "identity" << repeatXYZ("position") << repeatXYZ("force") << SEP << "tension";
    
    // we print Fibers in the order of the inventory:
    Fiber * fib = fibers.firstID();
    while ( fib )
    {
        out << "\n% fiber " << fib->reference();
            
        for ( unsigned p = 0; p < fib->nbPoints(); ++p )
        {
            out << LIN << fib->identity();
            out << SEP << fib->posP(p);
            out << SEP << fib->netForce(p);
            if ( p == fib->lastPoint() )
                out << SEP << 0.0;
            else
                out << SEP << fib->tension(p);
        }
        
        fib = fibers.nextID(fib);
    }
}


/**
 Export total magnitude of force exerted by Fiber on the confinement
 */
void Simul::reportFiberConfinement(std::ostream& out) const
{
    out << COM << " radial_force";
    real sum = 0;
    for ( Fiber * fib=fibers.first(); fib; fib=fib->next() )
    {
        Space const* spc = findSpace(fib->prop->confine_space);
        const real stiff = fib->prop->confine_stiffness;
        
        for ( unsigned p = 0; p < fib->nbPoints(); ++p )
        {
            Vector w, pos = fib->posP(p);
            if ( spc->outside(pos) )
            {
                spc->project(pos, w);
#if ( DIM > 1 )
                Vector dir = Vector(pos.XX, pos.YY, 0).normalized();
                sum += stiff * ( ( pos - w ) * dir );
#endif
            }
        }
    }
    out << LIN << sum;
}



/**
 Export positions of points taken randomly along all fibers,
 but that remain static with respect to the lattice of each fiber,
 during the life-time of this fiber.
 
 This is meant to simulate the `speckle microscopy` that is obtained
 in microcscopy with a low amount of fluorescent-monomers.
 
 The distance between the speckles follows an exponential distribution
 with an average defined by the parameter `spread`.
 */
void Simul::reportFiberSpeckles(std::ostream& out, Glossary& opt) const
{
    const real tiny = 0x1p-32;
    real spread = 1;
    if ( opt.set(spread, "density") )
        spread = 1.0 / spread;
    else
        opt.set(spread, "interval");

    Fiber * fib = fibers.first();
    while ( fib )
    {
        out << "\n% fiber " << fib->reference();
        
        // generate speckles below the origin of abscissa
        if ( fib->abscissaM() < 0 )
        {
            uint32_t z = fib->signature();
            real a = spread * log(z*tiny);
            while ( a > fib->abscissaP() )
            {
                z = lcrng2(z);
                a += spread * log(z*tiny);
            }
            while ( a >= fib->abscissaM() )
            {
                out << '\n' << fib->pos(a);
                z = lcrng2(z);
                a += spread * log(z*tiny);
            }
        }
        // generate speckles above the origin of abscissa
        if ( fib->abscissaP() > 0 )
        {
            uint32_t z = ~fib->signature();
            real a = -spread * log(z*tiny);
            while ( a < fib->abscissaM() )
            {
                z = lcrng1(z);
                a -= spread * log(z*tiny);
            }
            while ( a <= fib->abscissaP() )
            {
                out << '\n' << fib->pos(a);
                z = lcrng1(z);
                a -= spread * log(z*tiny);
            }
        }
        
        fib = fib->next();
    }
}

/**
 Export positions of points taken randomly along all fibers,
 changing the distribution at every time.
 */
void Simul::reportFiberSamples(std::ostream& out, Glossary& opt) const
{
    real spread = 1;
    if ( opt.set(spread, "density") )
        spread = 1.0 / spread;
    else
        opt.set(spread, "interval");
    
    Array<FiberBinder> loc(1024);
    fibers.uniFiberSites(loc, spread);
    
    Fiber * ofib = 0;
    for ( Array<FiberBinder>::iterator i = loc.begin(); i < loc.end(); ++i )
    {
        FiberBinder & site = *i;
        
        if ( ofib != site.fiber() )
        {
            out << "\n% fiber " << site.fiber()->reference();
            ofib = site.fiber();
        }
        
        out << LIN << i->pos();
    }
}


/**
 Export fiber elastic bending energy
 */
void Simul::reportFiberBendingEnergy(std::ostream& out) const
{
    out << COM << ljust("class",2,1) << SEP << "count" << SEP << "sum" << SEP << "avg" << SEP << "dev";
    
    unsigned cnt;
    real avg, dev;
    
    PropertyList plist = properties.find_all("fiber");
    
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        FiberProp * fp = static_cast<FiberProp*>(*ip);
        ObjectList objs = fibers.collect(fp);
        fibers.infoBendingEnergy(objs, cnt, avg, dev);
        if ( cnt > 0 )
        {
            out << LIN << ljust(fp->name(), 2);
            out << SEP << cnt;
            out << SEP << avg*cnt;
            out << SEP << avg;
            out << SEP << dev;
        }
    }
}


/**
 Sum of the internal tensions from fiber segments that intersect a plane
 specified in `opt`.
 The plane is defined by <em> n.pos + a = 0 </em>
 @code
 plane = NORMAL, SCALAR
 @endcode
 */
void Simul::reportFiberTension(std::ostream& out, Glossary& opt) const
{
    computeForces();
    
    out << COM << "count" << SEP << "force";
    for ( int d = 1; d < DIM; ++d )
        out << SEP << "count" << SEP << "force";

    real a = 0;
    Vector n(1,0,0);
    real ten = 0;
    unsigned cnt = 0;
    if ( opt.value_is("plane", 0, "all") )
    {
        fibers.infoTension(cnt, ten, Vector(1,0,0), 0);
        out << LIN << cnt;
        out << SEP << ten;
#if ( DIM > 1 )
        fibers.infoTension(cnt, ten, Vector(0,1,0), 0);
        out << SEP << cnt;
        out << SEP << ten;
#endif
#if ( DIM > 2 )
        fibers.infoTension(cnt, ten, Vector(0,0,1), 0);
        out << SEP << cnt;
        out << SEP << ten;
#endif
        return;
    }
    if ( opt.set(n, "plane") )
    {
        opt.set(a, "plane", 1);
        fibers.infoTension(cnt, ten, n, a);
        out << "\n% plane (" << n << ").pos + " << a << " = 0";
    }
    else
    {
        fibers.infoTension(cnt, ten);
    }
    
    out << LIN << cnt;
    out << SEP << ten;
}


void Simul::reportFiberIntersections(std::ostream& out, Glossary& opt) const
{
    int details = 2;
    real up = 0;
    opt.set(up, "distance");
    opt.set(details, "details");
    
    const real mds = up * up;
    real abs1, abs2, dis;
    
    if ( details == 2 )
    {
        out << COM << "id1" << SEP << "abs1";
        out << SEP << "id2" << SEP << "abs2" << SEP << repeatXYZ("position");
    }
    Accumulator accum;
    
    Fiber * fib = fibers.firstID();
    
    while ( fib )
    {
        unsigned cnt = 0;
        Fiber * fob = fibers.nextID(fib);
        while ( fob )
        {
            for ( unsigned ii = 0; ii < fib->nbSegments(); ++ii )
            {
                FiberLocus & loc1 = fib->locus(ii);
                for ( unsigned jj = 0; jj < fob->nbSegments(); ++jj )
                {
                    FiberLocus & loc2 = fob->locus(jj);
                    if ( loc1.shortestDistance(loc2, abs1, abs2, dis) )
                    {
                        if ( dis <= mds )
                        {
                            ++cnt;
                            Vector pos1 = loc1.pos(abs1/loc1.len());
                            //Vector pos2 = loc2.pos(abs2/loc2.len());
                            if ( details == 2 )
                            {
                                out << LIN << fib->identity();
                                out << SEP << abs1 + loc1.abscissa1();
                                out << SEP << fob->identity();
                                out << SEP << abs2 + loc2.abscissa1();
                                out << SEP << pos1;
                            }
                            accum.add(pos1);
                        }
                    }
                }
            }
            fob = fibers.nextID(fob);
        }
        if ( cnt && details >= 1 )
        {
            out << COM << "total";
            out << SEP << fib->identity();
            out << SEP << cnt;
        }
        fib = fibers.nextID(fib);
    }
    accum.subtract_mean();
    accum.print_doc(out);
    accum.print(out, 1);
}


void Simul::reportFiberBinders(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << SEP << "abs1" << SEP << "abs2" << SEP << "...";
    for ( Fiber * fib=fibers.first(); fib; fib=fib->next() )
    {
        fib->sortBinders();
        out << LIN << fib->prop->index();
        out << SEP << fib->identity();
        for ( FiberBinder * fb = fib->firstBinder(); fb; fb = fb->next() )
        {
            Hand * h = static_cast<Hand*>(fb)->otherHand();
            if ( h && h->attached() )
                out << SEP << fb->abscissa();
        }
    }
}

void Simul::reportFiberBindersAll(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << SEP << "abs1" << SEP << "abs2" << SEP << "...";\
    
    for ( Fiber * fib=fibers.first(); fib; fib=fib->next() )
    {
        fib->sortBinders();
        out << LIN << fib->prop->index();
        out << SEP << fib->identity();
        for ( FiberBinder * fb = fib->firstBinder(); fb; fb = fb->next() )
        {
            Hand * h = static_cast<Hand*>(fb);
            if (h->otherHand())
            {
                out << SEP << fb->abscissa();
            }
            
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Networks


/// accessory class to analyse the connections in a network of fibers
struct Connector
{
    real a;
    int  f;
    int  g;
    int  h;
    real s;
    Connector(real as, int fs) { a = as; f = fs; g = -1; h = -1; s = 0; }
    Connector(real as, int fs, int gs) { a = as; f = fs; g = gs; h = -1; s = 0; }
    Connector(real as, int fs, int gs, int hs) { a = as; f = fs; g = gs; h = hs; s = 0; }
    Connector(real as, int fs, int gs, int hs, real ss) { a = as; f = fs; g = gs; h = hs; s = ss; }
    bool operator < (const Connector& r) const { return a < r.a; }
};

bool comp_real(const real& a, const real& b)
{
    return a < b;
}

real sum_real(const real& a, const real& b)
{
    return a + b;
}

/**
 This is the older version
 */
void Simul::reportFiberConnectors(std::ostream& out, Glossary& opt) const
{
    int details = 2;
    opt.set(details, "details");

    if ( details > 1 )
    {
        out << "\n%    class  identity      abs1    fiber1     hand1      dist";
        out << "      abs2    fiber2     hand2      dist ...";
    }
    else
    {
        out << "\n% fiber connectors";
    }
    
    // used to calculate the size of the network from the position of connectors
    Accumulator accum;
    
    typedef std::vector<Connector> clist_t;
    typedef std::map<ObjectID, clist_t> map_t;
    
    map_t map;

    Fiber * fib = fibers.firstID();
    
    while ( fib )
    {
        map.clear();
        // check all connecting Hands and record abscissa, depending on the fiber that is linked
        for ( FiberBinder * fb = fib->firstBinder(); fb; fb = fb->next() )
        {
            Hand * h1 = static_cast<Hand*>(fb);
            Hand * h2 = h1->otherHand();
            if ( h2 && h2->attached() )
            {
                int f2 = h2->fiber()->identity();
                map[f2].push_back(Connector(fb->abscissa(), h1->prop->index()));
            }
        }
        if ( map.size() )
        {
            if ( details > 1 )
            {
                out << LIN << fib->prop->index() << SEP << fib->identity();
            }
            
            clist_t list;
            // average all the abscissa linking to the same fiber:
            for ( map_t::const_iterator mi = map.begin(); mi != map.end(); ++mi )
            {
                clist_t const& sec = mi->second;
                // average abscissa:
                //real a = std::accumulate(sec.begin(), sec.end(), (real)0.0, sum_real) / sec.size();
                real a = 0.0;
                // number of connector of each type
                int c1 = 0, c2 = 0;
                for ( clist_t::const_iterator ci = sec.begin(); ci != sec.end(); ++ci )
                {
                    a += ci->a;
                    if ( ci->f == 1 ) ++c1;
                    if ( ci->f == 2 ) ++c2;
                 }
                a /= sec.size();
                list.push_back(Connector(a, mi->first, c1, c2, sec.size()-c1-c2));
            }
            // sort the list in increasing abscissa
            std::sort(list.begin(), list.end());
            
            clist_t::const_iterator p = list.begin();
            for ( clist_t::const_iterator c = list.begin(); c != list.end(); ++c )
            {
                if ( details > 1 )
                {
                    out << SEP << c->a;
                    out << SEP << c->f;
                    out << SEP << sMath::repr(c->g)+"+"+sMath::repr(c->h);
                    // calculate direct distance to previous point of intersection:
                    if ( c != p )
                        out << SEP << ( fib->pos(p->a) - fib->pos(c->a) ).norm();
                    else
                        out << SEP << 0;
                }
                p = c;
                accum.add(fib->pos(p->a));
            }
            if ( details > 0 )
            {
                out << COM << "total";
                out << SEP << fib->prop->index();
                out << SEP << fib->identity();
                out << SEP << list.size();
            }
        }
        fib = fibers.nextID(fib);
    }
    accum.subtract_mean();
    accum.print_doc(out);
    accum.print(out, 1);
}



#include "motor_prop.h"

real hand_speed(HandProp const* hp)
{
    if ( hp->activity == "move" )
        return static_cast<MotorProp const*>(hp)->unloaded_speed;
    return 0;
}

/**
 F. Nedelec, 18/08/2017
 */
void Simul::reportNetworkBridges(std::ostream& out, Glossary& opt) const
{
    int details = 0;
    opt.set(details, "details");

    out << COM << "length" << SEP << "speed" << SEP << "type1" << SEP << "type2";

    typedef std::vector<Connector> clist_t;
    typedef std::map<ObjectID, clist_t> map_t;
    
    map_t map;
    
    HandProp const* hp1 = static_cast<HandProp*>(properties.find_or_die("hand", 1));
    HandProp const* hp2 = static_cast<HandProp*>(properties.find_or_die("hand", 2));
    
    const real speedh1 = hand_speed(hp1);
    const real speedh2 = hand_speed(hp2);

    Fiber * fib = fibers.firstID();
    
    while ( fib )
    {
        map.clear();
        // check all connecting Hands and record abscissa, depending on the fiber that is linked
        for ( FiberBinder * fb = fib->firstBinder(); fb; fb = fb->next() )
        {
            Hand * h1 = static_cast<Hand*>(fb);
            Hand * h2 = h1->otherHand();
            if ( h2 && h2->attached() )
            {
                int f2 = h2->fiber()->identity();
                if ( h1->prop == hp1 )
                    map[f2].push_back(Connector(fb->abscissa(), 1));
                else if ( h1->prop == hp2 )
                    map[f2].push_back(Connector(fb->abscissa(), 2));
                else
                    out << COM << "report network:bridge can only handle 2 hand types";
            }
        }
        if ( map.size() )
        {
            clist_t list;
            // average all the abscissa linking to the same fiber:
            for ( map_t::const_iterator mi = map.begin(); mi != map.end(); ++mi )
            {
                clist_t const& sec = mi->second;
                // average abscissa:
                real a = 0.0;
                // number of connector of each type
                int c1 = 0, c2 = 0;
                for ( clist_t::const_iterator ci = sec.begin(); ci != sec.end(); ++ci )
                {
                    a += ci->a;
                    if ( ci->f == 1 ) ++c1;
                    if ( ci->f == 2 ) ++c2;
                }
                a /= sec.size();
                real speed = 0;
                if ( c1 > 0 && c2 > 0 )
                    speed = std::min(speedh1, speedh2);
                else if ( c1 > 0 )
                    speed = speedh1;
                else if ( c2 > 0 )
                    speed = speedh2;
                list.push_back(Connector(a, mi->first, c1, c2, speed));
            }
            // sort the list in increasing abscissa
            std::sort(list.begin(), list.end());
            
            if ( details > 0 )
            {
                out << COM << " connectors on fiber f" << fib->identity() << ":";
                // print all connector attachment positions:
                out << COM << "abscissa" << SEP << "fiber_id" << SEP << "speed" << SEP << "type";
                for ( clist_t::const_iterator c = list.begin(); c != list.end(); ++c )
                {
                    out << LIN << c->a;
                    out << SEP << c->f;
                    out << SEP << c->s;
                    out << SEP << sMath::repr(c->g)+"+"+sMath::repr(c->h);
                }
            }
            if ( list.size() > 1 )
            {
                // print all bridges
                if ( details > 0 )
                    out << COM << "length" << SEP << "speed" << SEP << "type1" << SEP << "type2";
#if ( 1 )
                for ( clist_t::const_iterator p = list.begin(); p != list.end(); ++p )
                for ( clist_t::const_iterator c = p+1; c != list.end(); ++c )
#else
                for ( clist_t::const_iterator p = list.begin(), c = p+1; c != list.end(); ++p, ++c )
#endif
                {
                    out << LIN << c->a - p->a;
                    out << SEP << c->s - p->s;
                    out << SEP << sMath::repr(p->g)+"+"+sMath::repr(p->h);
                    out << SEP << sMath::repr(c->g)+"+"+sMath::repr(c->h);
                }
            }
        }
        fib = fibers.nextID(fib);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Objects


void Simul::reportTime(std::ostream& out) const
{
    out << LIN << prop->time;
}


void Simul::reportInventory(std::ostream& out) const
{
    //out << COM << " properties:";
    //properties.write_names(out, "");
    out << COM << " objects:\n";
    spaces.report(out);
    fields.report(out);
    fibers.report(out);
    spheres.report(out);
    beads.report(out);
    solids.report(out);
    singles.report(out);
    couples.report(out);
    organizers.report(out);
    events.report(out);
}


/**
 Export position of all organizers
 */
void Simul::reportOrganizer(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("position");

    for ( Organizer * obj=organizers.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->property()->index();
        out << SEP << obj->identity();
        out << SEP << obj->position();
        out << SEP << obj->nbOrganized();
    }
}


/**
 Export position of Asters
 */
void Simul::reportAster(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("position");
    
    for ( Organizer * obj=organizers.first(); obj; obj=obj->next() )
    {
        if ( obj->tag() == Aster::TAG )
        {
            out << LIN << obj->property()->index();
            out << SEP << obj->identity();
            out << SEP << obj->position();
        }
    }
}


/**
 Export position of Beads
 */
void Simul::reportBeadPosition(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("position");
    
    for ( Bead * obj=beads.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->prop->index();
        out << SEP << obj->identity();
        out << SEP << obj->position();
    }
}


/**
 Export number of beads classified as a function of
 the number of grafted Single that are attached to Fibers
 */
void Simul::reportBeadSingles(std::ostream& out) const
{
    out << COM << "identity" << "amount(nb_attached_hands)";
    
    std::map<ObjectID, int> cnt;
    
    for ( Single * sig=singles.firstA(); sig; sig=sig->next() )
    {
        Mecable const* mec = sig->base();
        if ( mec && mec->tag() == Bead::TAG )
            ++cnt[ mec->identity() ];
    }

    const int max = 12;
    int nb[max] = { 0 };
    for ( Bead * obj=beads.first(); obj; obj=obj->next() )
        ++nb[ cnt[obj->identity()] ];
    
    for ( int c = 0; c < max; ++c )
        out << " " << std::setw(3) << nb[c];
}


/**
 Export position of Solids
 */
void Simul::reportSolid(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << repeatXYZ("centroid");
    out << repeatXYZ("point0") << repeatXYZ("point1");
    
    for ( Solid * obj=solids.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->prop->index();
        out << SEP << obj->identity();
        out << SEP << obj->centroid();
        out << SEP << obj->posP(0);
        if ( obj->nbPoints() > 1 )
            out << SEP << obj->posP(1);

        if ( modulo ) 
        {
            Vector pos = obj->centroid();
            modulo->fold(pos);
            out << SEP << pos;
        }
    }
}


/**
 Report position of Sphere
 */
void Simul::reportSphere(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity";
    out << repeatXYZ("point0") << repeatXYZ("point1");
  
    for ( Sphere * obj=spheres.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->prop->index();
        out << SEP << obj->identity();
        out << SEP << obj->posP(0);
        if ( obj->nbPoints() > 1 )
            out << SEP << obj->posP(1);
    }
}



/**
 Report something about Space (incomplete)
 */
void Simul::reportSpace(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity";
    
    for ( Space * obj=spaces.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->prop->name();
        out << SEP << obj->identity();
        for ( unsigned d = 0; d < Space::DMAX; ++d )
            out << SEP << std::fixed << obj->length(d);
    }
}


/**
 Report quantity of substance in Field
 */
void Simul::reportField(std::ostream& out) const
{
    out << COM << ljust("class", 2, 1);
    out << SEP << "total" << SEP << "avg" << SEP << "min" << SEP << "max";
    
    // report total substance in each Field
    for ( Field * obj=fields.first(); obj; obj=obj->next() )
    {
        real vol = obj->FieldGrid::cellVolume();
        Field::value_type s, n, x;
        obj->infoValues(s, n, x);
        out << LIN << ljust(obj->prop->name(), 2);
        out << SEP << s;
        out << SEP << s / vol;
        out << SEP << n / vol;
        out << SEP << x / vol;
    }
    
    // report substance on Fiber Lattices
    real len, sm, mn, mx;
    fibers.infoLattice(len, sm, mn, mx);
    out << LIN << ljust("fiber:lattice", 2);
    out << SEP << sm;
    out << SEP << sm / len;
    out << SEP << mn / len;
    out << SEP << mx / len;
}


//------------------------------------------------------------------------------
#pragma mark - Single

void writeF(std::ostream& out, Single * obj)
{
    assert_true( !obj->attached() );
    out << LIN << obj->prop->index();
    out << SEP << obj->identity();
    out << SEP << obj->position();
    out << SEP << Vector(0,0,0);
    out << SEP << "0";
    out << SEP << "nan";
    out << SEP << "0";
}

void writeA(std::ostream& out, Single * obj, Simul const* simul)
{
    assert_true( obj->attached() );
    out << LIN << obj->prop->index();
    out << SEP << obj->identity();
    out << SEP << obj->position();
    out << SEP << obj->force();
    Fiber const* fib = obj->fiber();
    out << SEP << fib->identity();
    out << SEP << obj->abscissa();
    Organizer * o = simul->organizers.findOrganizer(fib);
    if ( o )
        out << SEP << static_cast<Object*>(o)->identity();
    else
        out << SEP << "0";
}


/**
 Export details of Singles
 */
void Simul::reportSingleState(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity";
    out << repeatXYZ("position") << repeatXYZ("force");
    out << "     fiber  abscissa     aster";
    
    for ( Single * obj=singles.firstF(); obj ; obj=obj->next() )
        writeF(out, obj);
    
    for ( Single * obj=singles.firstA(); obj ; obj=obj->next() )
        writeA(out, obj, this);
}


/**
 Export details of Singles of a certain kind
 */
void Simul::reportSingleState(std::ostream& out, std::string const& which) const
{
    Property * prop = properties.find_or_die("single", which);
    
    out << COM << "class" << SEP << "identity";
    out << repeatXYZ("position") << repeatXYZ("force");
    out << SEP << "fiber" << SEP << "abscissa" << SEP << "aster";
    
    for ( Single * obj=singles.firstF(); obj ; obj=obj->next() )
        if ( obj->prop == prop )
            writeF(out, obj);
    
    for ( Single * obj=singles.firstA(); obj ; obj=obj->next() )
        if ( obj->prop == prop )
            writeA(out, obj, this);
}


/**
 Export details of attached Singles
 */
void Simul::reportSinglePosition(std::ostream& out, std::string const& which) const
{
    Property * prop = 0;
    
    if ( which.size() )
        prop = properties.find_or_die("single", which);

    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("hand") << repeatXYZ("foot");
    out << SEP << "fiber" << SEP << "abscissa";

    for ( Single * obj=singles.firstA(); obj ; obj=obj->next() )
    {
        if ( !prop  ||  obj->prop == prop )
        {
            out << LIN << obj->prop->index();
            out << SEP << obj->identity();
            out << SEP << obj->posHand();
            out << SEP << obj->posFoot();
            out << SEP << obj->fiber()->identity();
            out << SEP << obj->abscissa();
        }
    }
}


/**
 Export number of Single in each state
 */
void Simul::reportSingle(std::ostream& out) const
{
    PropertyList plist = properties.find_all("single");
    
    const unsigned mx = 128;
    
    int nb[mx] = { 0 }, free[mx] = { 0 }, bound[mx] = { 0 };
    
    for ( Single * si = singles.firstF(); si ; si = si->next() )
    {
        unsigned ix = si->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            ++free[ix];
        }
    }
    
    for ( Single * si = singles.firstA(); si ; si = si->next() )
    {
        unsigned ix = si->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            ++bound[ix];
        }
    }
    
    if ( 1 )
    {
        out << COM << ljust("single", 2, 1);
        out << SEP << "total";
        out << SEP << "free";
        out << SEP << "bound";
    }
    
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        out << LIN << ljust((*ip)->name(), 2);
        unsigned ix = (*ip)->index();
        if ( ix < mx )
        {
            out << SEP << nb[ix];
            out << SEP << free[ix];
            out << SEP << bound[ix];
        }
        else
            out << SEP << " out-of-range ";
    }
}


//------------------------------------------------------------------------------
#pragma mark - Couple


void write(std::ostream& out, Couple * obj)
{
    out << LIN << obj->prop->index();
    out << SEP << obj->identity();
    out << SEP << obj->active();
    out << SEP << obj->position();
    if (obj->linking())
        out << SEP << obj->force();
    else
        out << SEP << Vector(0,0,0);
    
    Fiber const* fib = obj->fiber1();
    if ( fib )
    {
        out << SEP << fib->identity();
        out << SEP << obj->abscissa1();
    }
    else
    {
        out << SEP << "0";
        out << SEP << "nan";
    }

    fib = obj->fiber2();
    if ( fib )
    {
        out << SEP << fib->identity();
        out << SEP << obj->abscissa2();
    }
    else
    {
        out << SEP << "0";
        out << SEP << "nan";
    }
    
}

/**
 Export position of Couples
 */
void Simul::reportCoupleState(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << SEP << "active" << repeatXYZ("position")<<repeatXYZ("force");
    out << SEP << "fiber1" << SEP << "abscissa1" << SEP << "fiber2" << SEP << "abscissa2";

    for ( Couple * obj=couples.firstFF(); obj ; obj=obj->next() )
        write(out, obj);
    
    for ( Couple * obj=couples.firstAF(); obj ; obj=obj->next() )
        write(out, obj);
    
    for ( Couple * obj=couples.firstFA(); obj ; obj=obj->next() )
        write(out, obj);
    
    for ( Couple * obj=couples.firstAA(); obj ; obj=obj->next() )
        write(out, obj);
}

/**
 Export position of Couples of a certain kind
 */
void Simul::reportCoupleState(std::ostream& out, std::string const& which) const
{
    Property * prop = properties.find_or_die("couple", which);
    
    out << COM << "class" << SEP << "identity" << SEP << "active" << repeatXYZ("position")<<repeatXYZ("force");
    out << SEP << "fiber1" << SEP << "abscissa1" << SEP << "fiber2" << SEP << "abscissa2";
    
    for ( Couple * obj=couples.firstFF(); obj ; obj=obj->next() )
        if ( obj->prop == prop )
            write(out, obj);
    
    for ( Couple * obj=couples.firstAF(); obj ; obj=obj->next() )
        if ( obj->prop == prop )
            write(out, obj);
    
    for ( Couple * obj=couples.firstFA(); obj ; obj=obj->next() )
        if ( obj->prop == prop )
            write(out, obj);
    
    for ( Couple * obj=couples.firstAA(); obj ; obj=obj->next() )
        if ( obj->prop == prop )
            write(out, obj);
}


/**
 Export position of active Couples of a certain kind
 */
void Simul::reportCoupleActive(std::ostream& out, std::string const& which) const
{
    Property * prop = properties.find_or_die("couple", which);
    
    out << "\n% state" << repeatXYZ("position");
    
    for ( Couple * obj=couples.firstFF(); obj ; obj=obj->next() )
        if ( obj->active()  &&  obj->prop == prop )
            out << "\n 0 " << obj->position();
   
    for ( Couple * obj=couples.firstAF(); obj ; obj=obj->next() )
        if ( obj->prop == prop )
            out << "\n 1 " << obj->position();
    
    for ( Couple * obj=couples.firstFA(); obj ; obj=obj->next() )
        if ( obj->prop == prop )
            out << "\n 2 " << obj->position();
    
    for ( Couple * obj=couples.firstAA(); obj ; obj=obj->next() )
        if ( obj->prop == prop )
            out << "\n 3 " << obj->position();
}


/**
 Export position and force of Couples that are bound to 2 filaments
 */
void Simul::reportCoupleLink(std::ostream& out, std::string const& which) const
{
    Property * prop = 0;
    
    if ( which.size() )
        prop = properties.find_or_die("couple", which);
   
    out << COM << "class" << SEP << "identity";
    out << SEP << "fiber1" << SEP << "abscissa1";
    out << SEP << "fiber2" << SEP << "abscissa2";
    out << SEP << "force_nrm" << SEP << "cos_angle";
#if ( 0 )
    out << SEP << "pos_hand1" << SEP << "pos_hand2";
#endif
    
    for ( Couple * obj=couples.firstAA(); obj ; obj=obj->next() )
    {
        if ( !prop  ||  obj->prop == prop )
        {
            out << LIN << obj->prop->index();
            out << SEP << obj->identity();
            
            Fiber const* fib = obj->fiber1();
            out << SEP << fib->identity();
            out << SEP << obj->abscissa1();
            
            fib = obj->fiber2();
            out << SEP << fib->identity();
            out << SEP << obj->abscissa2();
            
            out << SEP << obj->force().norm();
            out << SEP << obj->dirFiber1() * obj->dirFiber2();
#if ( 0 )
            out << SEP << obj->posHand1();
            out << SEP << obj->posHand2();
#endif
        }
    }
}
void Simul::reportCouplePosition(std::ostream& out, std::string const& which) const
{
    Property * prop = 0;
    
    if ( which.size() )
        prop = properties.find_or_die("couple", which);
    
    out << COM << "class" << SEP << "identity";
#if ( 1 )
    out << SEP << "pos_hand1" << SEP << "pos_hand2";
#endif
    
    for ( Couple * obj=couples.firstAA(); obj ; obj=obj->next() )
    {
        if ( !prop  ||  obj->prop == prop )
        {
            out << LIN << obj->prop->index();
            out << SEP << obj->identity();
#if ( 1 )
            out << SEP << obj->posHand1();
            out << SEP << obj->posHand2();
#endif
        }
    }
}

/**
 Export histogram of Couples forces
 */
void Simul::reportCoupleForce(std::ostream& out, Glossary& opt) const
{
    const unsigned imax = 8;
    const unsigned xbin = 128;
    unsigned nbin = 64;
    real delta = 0.5;
    
    opt.set(delta,  "interval");
    opt.set(nbin, "bins");
    if ( nbin > xbin )
        nbin = xbin;

    int  cnt[imax][xbin+1];
    
    // reset counts:
    for ( unsigned ii = 0; ii <  imax; ++ii )
    for ( unsigned jj = 0; jj <= xbin; ++jj )
        cnt[ii][jj] = 0;
    
    // accumulate counts:
    for ( Couple * cxi=couples.firstAA(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < imax )
        {
            int f = (int)( cxi->force().norm() / delta );
            if ( f < xbin )
                ++cnt[ix][f];
            else
                ++cnt[ix][xbin];
        }
    }
    
    if ( 1 )
    {
        out << COM << "force_distribution" << " (`scale` indicates the center of each bin)";
        out << LIN << ljust("scale", 2);
        for ( unsigned u = 0; u <= nbin; ++u )
            out << " " << std::setw(5) << delta * ( u + 0.5 );
    }
    
    for ( unsigned ii = 0; ii < imax; ++ii )
    {
        unsigned sum = 0;
        for ( unsigned jj = 0; jj < nbin; ++jj )
            sum += cnt[ii][jj];
        if ( sum )
        {
            Property const* ip = properties.find("couple", ii);
            out << LIN << ljust(ip->name(), 2);
            for ( unsigned jj = 0; jj <= nbin; ++jj )
                out << ' ' << std::setw(5) << cnt[ii][jj];
        }
    }
}


/**
 Export number of Couples in each state
 */
void Simul::reportCouple(std::ostream& out) const
{
    PropertyList plist = properties.find_all("couple");
    
    const unsigned mx = 128;
    int nb[mx] = { 0 }, act[mx] = { 0 }, cnt[mx][4];
    
    //reset counts:
    for ( unsigned ii = 0; ii < mx; ++ii )
    {
        cnt[ii][0] = 0;
        cnt[ii][1] = 0;
        cnt[ii][2] = 0;
        cnt[ii][3] = 0;
    }
    
    for ( Couple * cxi=couples.firstFF(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if ( cxi->active() ) ++act[ix];
            ++(cnt[ix][0]);
        }
    }
    
    for ( Couple * cxi=couples.firstAF(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if ( cxi->active() ) ++act[ix];
            ++(cnt[ix][1]);
        }
    }
    for ( Couple * cxi=couples.firstFA(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if ( cxi->active() ) ++act[ix];
            ++(cnt[ix][2]);
        }
    }
    
    for ( Couple * cxi=couples.firstAA(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if ( cxi->active() ) ++act[ix];
            ++(cnt[ix][3]);
        }
    }
    
    if ( 1 )
    {
        out << COM << ljust("couple", 2, 1);
        out << SEP << "total";
        out << SEP << "active";
        out << SEP << "FF";
        out << SEP << "AF";
        out << SEP << "FA";
        out << SEP << "AA";
    }
    
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        out << LIN << ljust((*ip)->name(), 2);
        unsigned ix = (*ip)->index();
        if ( ix < mx )
        {
            out << SEP << nb[ix];
            out << SEP << act[ix];
            for ( unsigned int d = 0; d < 4; ++d )
                out << SEP << cnt[ix][d];
        }
        else
            out << SEP << "out-of-range";
    }
}
void Simul::reportFlipper(std::ostream& out) const
{
    PropertyList plist = properties.find_all("couple");
    
    const unsigned mx = 128;
    int nb[mx] = { 0 }, act[mx] = { 0 }, cnt[mx][7];
    
    //reset counts:
    for ( unsigned ii = 0; ii < mx; ++ii )
    {
        cnt[ii][0] = 0;
        cnt[ii][1] = 0;
        cnt[ii][2] = 0;
        cnt[ii][3] = 0;
        cnt[ii][4] = 0;
        cnt[ii][5] = 0;
        cnt[ii][6] = 0;
        
    }
    
    for ( Couple * cxi=couples.firstFF(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if ( cxi->active() ) ++act[ix];
            ++(cnt[ix][0]);
        }
    }
    
    for ( Couple * cxi=couples.firstAF(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if ( cxi->active() ) ++act[ix];
            ++(cnt[ix][1]);
        }
    }
    for ( Couple * cxi=couples.firstFA(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if ( cxi->active() ) ++act[ix];
            ++(cnt[ix][2]);
        }
    }
    
    for ( Couple * cxi=couples.firstAA(); cxi ; cxi = cxi->next() )
    {
        //Only for the flippers
        if (cxi->prop->activity!="flip")
            continue;
//        Flipper * c = nullptr;
        unsigned ix = cxi->prop->index();
        unsigned indexer = 3;
        if ( ix < mx )
        {
            ++nb[ix];
            // Add 1 if hand 1 is flipped
            indexer+=(cxi->hand1()->prop!=cxi->prop->hand1_prop);
            // Add 2 if hand 2 is flipped
            indexer+=(cxi->hand2()->prop!=cxi->prop->hand2_prop)*2;
            if ( cxi->active() ) ++act[ix];
            ++(cnt[ix][indexer]);
        }
    }
    
    if ( 1 )
    {
        out << COM << ljust("couple", 2, 1);
        out << SEP << "total";
        out << SEP << "active";
        out << SEP << "FF";
        out << SEP << "AF";
        out << SEP << "FA";
        out << SEP << "A1A1";
        out << SEP << "A2A1";
        out << SEP << "A1A2";
        out << SEP << "A2A2";
    }
    
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        out << LIN << ljust((*ip)->name(), 2);
        unsigned ix = (*ip)->index();
        if ( ix < mx )
        {
            out << SEP << nb[ix];
            out << SEP << act[ix];
            for ( unsigned int d = 0; d < 7; ++d )
                out << SEP << cnt[ix][d];
        }
        else
            out << SEP << "out-of-range";
    }
}



/**
 Export composition of each type of Couple
 */
void Simul::reportCoupleHand(std::ostream& out) const
{
    PropertyList plist = properties.find_all("couple");
    
    out << COM << ljust("class", 2, 1);
    out << SEP << rjust("hand1", 2, 1) << SEP << rjust("hand2", 2, 1);

    bool mixed = false;
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        CoupleProp * p = static_cast<CoupleProp*>(*ip);
        out << LIN << ljust(p->name(), 2);
        out << SEP << rjust(p->hand1_prop->name(), 2, 1);
        out << SEP << rjust(p->hand2_prop->name(), 2, 1);
        if ( p->index() != p->hand1_prop->index() )
            mixed = true;
        if ( p->index() != p->hand2_prop->index() )
            mixed = true;
    }
    if ( mixed )
        out << "\nERROR: MIXED UP INDICES";
}

//------------------------------------------------------------------------------
#pragma mark - Clusters

/**
 Substitute the values of fiber:flag() such that `a` and `b`
 values are replaced by the smallest value between `a` and `b`.
 */
void reFleck(FiberSet const& set, int a, int b)
{
    // swap to make sure a is smallest of (a,b)
    if ( b < a )
        std::swap(a, b);

    for ( Fiber* fib = set.first(); fib; fib=fib->next() )
    {
        if ( fib->flag() == b )
            fib->flag(a);
    }
}


/**
 Set Fiber::flag() to indicated Fibers that are connected by Couple.
 
 The clusters are defined by the Couple that are bridging Fibers:
 Two fibers are in the same cluster if there is a Couple connecting them,
 of if they can be indirectly connected in this way via other Fibers.
 
 This analysis can be useful to identify mechanically isolated sub-networks
 in the simulation. 
 The result can be visualized in `play` with the option fiber:coloring=4, 
 and it can also be printed with the tool `report fiber:cluster`
 */
void Simul::flagClusters() const
{
    // set unique flag() for each fiber:
    for ( Fiber * fib = fibers.first(); fib; fib=fib->next() )
        fib->flag(fib->identity());
    
    // equalize flag() when fibers are connected by a Couple:
    for ( Couple * cx = couples.firstAA(); cx ; cx=cx->next() )
    {
        if ( cx->fiber1()->flag() != cx->fiber2()->flag() )
            reFleck(fibers, cx->fiber1()->flag(), cx->fiber2()->flag());
    }
}


/**
 Perform the same analysis as flagClusters(), 
 but only considering Couple of the given Property
 */
void Simul::flagClusters(Property const* cop) const
{
    // set unique flag() for each fiber:
    for ( Fiber * fib = fibers.first(); fib; fib=fib->next() )
        fib->flag(fib->identity());

    // equalize flag() when fibers are connected by a Couple:
    for ( Couple * cx = couples.firstAA(); cx ; cx=cx->next() )
    {
        if ( cx->prop == cop  &&  cx->fiber1()->flag() != cx->fiber2()->flag() )
            reFleck(fibers, cx->fiber1()->flag(), cx->fiber2()->flag());
    }
}



/// class to store info about a Cluster
struct ClusterInfo
{
    unsigned id;
    unsigned nb;
    ClusterInfo(int i, int n) { id=i; nb=n; }
    real operator < (ClusterInfo const&b) const { return nb > b.nb; }
};


/**
 Export size of clusters found by Simul::flagClusters()
 Clusters are ordered in decreasing size.
 */
void Simul::reportClusters(std::ostream& out, bool details) const
{
    flagClusters();

    typedef std::set<ObjectID> list_t;
    typedef std::map<unsigned, list_t> map_t;
    
    map_t map;
    
    for ( Fiber * fib = fibers.first(); fib; fib=fib->next() )
        map[fib->flag()].insert(fib->identity());

    // the std::set will keep its elements ordered:
    std::set<ClusterInfo> clusters;
    
    for ( map_t::const_iterator c = map.begin(); c != map.end(); ++c )
        if ( c->second.size() > 1 )
            clusters.insert(ClusterInfo(c->first, c->second.size()));
    
    // this will output clusters ordered by decreasing size:
    out << COM << "cluster" << SEP << "nb_fibers :" << SEP << "fiber_id";
    for ( std::set<ClusterInfo>::const_iterator cc = clusters.begin(); cc != clusters.end(); ++cc )
    {
        out << LIN << cc->id;
        out << SEP << cc->nb << " :";
        if ( details )
        {
            list_t & list = map[cc->id];
            for ( list_t::const_iterator f = list.begin(); f != list.end(); ++f )
                out << " " << (*f);
        }
    }
}


/**
 Evaluates if the Fiber distribution makes a connected ring around the Z-axis
 FJN 8.07.2017, for Blood Platelet project
 */
void Simul::reportRing(std::ostream& out) const
{
    flagClusters();
    
    typedef std::map<ObjectID, unsigned> map_t;
    typedef std::pair<ObjectID, unsigned> pair_t;
    
    map_t ring, sec;
    
    for ( Fiber const* fib = fibers.first(); fib; fib=fib->next() )
    {
        ring.insert(pair_t(fib->flag(), 0));
        ++ring[fib->flag()];
    }
    
    // rotate plane around the Z-axis and find intersecting fibers
    for ( real ang = 0; ang < 2*M_PI; ang += M_PI/360 )
    {
        Vector nor( cosf(ang), sinf(ang), 0);
        Vector dir(-sinf(ang), cosf(ang), 0);
        
        sec = ring;
        ring.clear();
        for ( Fiber const* fib=fibers.first(); fib; fib=fib->next() )
        {
            for ( unsigned s = 0; s < fib->nbSegments(); ++s )
            {
                // check that fiber intersect with plane:
                real abs = fib->planarIntersect(s, nor, 0);
                if ( 0 <= abs  &&  abs < 1 )
                {
                    // check that intersection is located on 'dir' side of Z-axis:
                    real H = fib->interpolatePoints(s,s+1,abs) * dir;
                    if ( H > 0 )
                    {
                        // check that cluster was already present before:
                        map_t::const_iterator f = sec.find(fib->flag());
                        if ( f != sec.end() )
                        {
                            ring.insert(pair_t(fib->flag(), 0));
                            ++ring[fib->flag()];
                        }
                    }
                }
            }
        }
        
        // retain the minimum number of fiber:
        for ( map_t::iterator cc = ring.begin(); cc != ring.end(); ++cc )
        {
            unsigned n = sec.find(cc->first)->second;
            if ( n < cc->second )
                cc->second = n;
        }
    }
    
    // output cluster and its minimum size:
    out << COM << "min_cnt" << SEP << "fiber_ids";
    for ( map_t::const_iterator cc = ring.begin(); cc != ring.end(); ++cc )
    {
        ObjectID f = cc->first;
        out << LIN << cc->second;
        
        for ( Fiber const* fib = fibers.first(); fib; fib=fib->next() )
        {
            if ( fib->flag() == f )
                out << " " << std::setw(5) << fib->identity();
        }
    }
    
    if ( ring.empty() )
        out << LIN << "0";
}



//------------------------------------------------------------------------------
#pragma mark - Misc

/**
 Export indices calculated by FiberSet::infoSpindle
 */
void Simul::reportIndices(std::ostream& out) const
{
    out << COM << "amount" << SEP << "radial" << SEP << "polar";
    real ixa, ixp;
    fibers.infoSpindle(ixa, ixp, Vector(1,0,0), 0, 30, 1);
    out << LIN << fibers.size();
    out << SEP << ixa;
    out << SEP << ixp;
}


/**
 Export number of Fibers pointing left and right,
 that intersect a plane parallel to YZ.
 The planes are distributed regularly every 0.5 um along the X-axis.
 */
void Simul::reportProfile(std::ostream& out) const
{
    out << COM << "position" << SEP << "left-pointing" << SEP << "right-pointing";
    Vector n(1,0,0);
    real m = 40, dm = 0.5;
    int nr, nl;
    for ( real p = -m ; p <= m ; p += dm )
    {
        fibers.infoPlane(nr, nl, n, -p);
        out << LIN << p;
        out << SEP << nl;
        out << SEP << nr;
    }
}


/**
 l'angle entre un vecteur 1 (centre du noyau --> SPB)
 et un vecteur 2 (axe de l'hyphe; gauche --> droite = sens du flow).
 */
void Simul::reportAshbya(std::ostream& out) const
{
    out << "\n% class id point_0, vector_1, angle";
    for ( Solid * obj=solids.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->prop->index();
        out << SEP << obj->identity();
        out << SEP << obj->posP(0);
        if ( obj->nbPoints() > 1 )
        {
            Vector vec = obj->diffPoints(0).normalized();
            Vector dir(1,0,0);
            out << SEP << vec;
            out << SEP << acos( vec * dir );
        }
    }
}


/**
 Export end-to-end distance of Fiber
 */
void Simul::reportCustom(std::ostream& out) const
{
    for ( Fiber * obj=fibers.first(); obj; obj=obj->next() )
    {
        Vector ee = obj->posEndP() - obj->posEndM();
        out << ee.norm() << " ";
    }
}



//ADDED BY MANU

void Simul::spindle_watch()
{
    //First, check that the spindle has grown at least a bit
    if  (prop->watch_dist_solids && time()>prop->watch_dist_solids_t)
    {
        Solid * s1 = solids.first();
        Solid * s2 = s1->next();
        if (s1->centroid().distance(s2->centroid())<prop->watch_dist_solids)
        {
            std::cout<< "Killed: distance between solids smaller than " << prop->watch_dist_solids << " at time "<<time() <<std::endl;
            exit(666);
        }
    }
    
    // Check that there are less MTs than x
    if (prop->watch_mt_nb && time()>prop->watch_mt_nb_t)
    {
        int count = 0;
        Fiber * f = fibers.first();
        while(f)
        {
            count++;
            f = f->next();
        }
        if (count>prop->watch_mt_nb)
        {
            std::cout<< "Killed: more than " << prop->watch_mt_nb << "fibers at time "<<time() <<std::endl;
            exit(666);
        }
    }
    
    // Check that there are x double bound couples of type y
    if (prop->watch_AAcouples && time()>prop->watch_AAcouples_t)
    {
        
        int count = 0;
        Couple * cou = couples.firstAA();

        while(cou)
        {
            if (cou->prop->name()==prop->watch_AAcouples_name)
                count += (cou->attached1() && cou->attached2());
            cou = cou->next();
        }
        if (count<prop->watch_AAcouples)
        {
            std::cout<< "Killed: less than " <<prop->watch_AAcouples <<" AA couples at time "<<time() <<std::endl;
            exit(666);
        }
    }
}

void Simul::reportFiberCap(std::ostream& out,Glossary& opt) const
{
    out << COM << "class" << SEP << "identity" << SEP << "length";
    out << SEP << "cap_length";
    unsigned int lat_val;
    opt.set(lat_val, "lat_val");
    for ( Fiber * obj=fibers.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->prop->index();
        out << SEP << obj->identity();
        out << SEP << obj->length();
        out << SEP << obj->measure_cap(lat_val);
    }
    out << std::endl;
}

void Simul::reportFiberOverlap(std::ostream& out,Glossary& opt) const
{
    out << COM << "class" << SEP << "identity" << SEP << "length";
    out << SEP << "overlaps";
    unsigned int lat_val;
    opt.set(lat_val, "lat_val");
    for ( Fiber * fib=fibers.first(); fib; fib=fib->next() )
    {
        for ( FiberBinder * fb = fib->firstBinder(); fb; fb = fb->next() )
        {
            out << LIN << fib->prop->index();
            out << SEP << fib->identity();
            out << SEP << fib->length();
            
            // A dictionary to count
            std::map<int, int> dict_max;
            std::map<int, int> dict_min;
            
            Hand * h = static_cast<Hand*>(fb)->otherHand();
            // Verify that the hand is bridging
            if ( h && h->attached() )
            {
                Fiber * other_fib = h->fiber();
                // Only count the pairs formed with fibers with higher identity, not to count the same overlap twice
                
                if (other_fib->identity()>fib->identity()) {
                    int id = other_fib->identity();
                    // If it exist, check whether its bigger or smaller than the max
                    if (dict_max.count(id))
                    {
                        if (h->abscissa()>dict_max[id])
                            dict_max[id]=h->abscissa();
                        if (h->abscissa()<dict_min[id])
                            dict_min[id]=h->abscissa();
                    }
                    else
                    {
                        dict_max[id] = h->abscissa();
                        dict_min[id] = h->abscissa();
                    }
                }
            }
        }

    }
    out << std::endl;
}

void Simul::reportFiberCount(std::ostream& out,Glossary& opt) const
{
    out << COM << "Number of fibers";
    unsigned int class_index;
    opt.set(class_index, "class_index");
    int cnt = 0;
    for ( Fiber * fib=fibers.first(); fib; fib=fib->next() )
    {
        if (fib->prop->index()==class_index)
            cnt+=1;
    }
    out << LIN << cnt;
    out << std::endl;
}

void Simul::reportCoupleTrap(std::ostream& out) const
{

    PropertyList plist = properties.find_all("couple");
    
    const unsigned mx = 128;
    int nb[mx] = { 0 }, act[mx] = { 0 }, cnt[mx][4], cnt_trapped[mx][4], cnt_trapped_att[mx][4];
    
    //reset counts:
    for ( unsigned ii = 0; ii < mx; ++ii )
    {
        cnt[ii][0] = 0;
        cnt[ii][1] = 0;
        cnt[ii][2] = 0;
        cnt[ii][3] = 0;
        cnt_trapped[ii][0] = 0;
        cnt_trapped[ii][1] = 0;
        cnt_trapped[ii][2] = 0;
        cnt_trapped[ii][3] = 0;
        cnt_trapped_att[ii][0] = 0;
        cnt_trapped_att[ii][1] = 0;
        cnt_trapped_att[ii][2] = 0;
        cnt_trapped_att[ii][3] = 0;
        
    }
    Hand * h1 = 0, * dummy = 0;
    for ( Couple * cxi=couples.firstFF(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if ( cxi->active() ) ++act[ix];
            if (cxi->trapped())
            {
                ++(cnt_trapped[ix][0]);
                cxi->trapped_haMon->getHands(h1, dummy);
                if (h1->attached())
                {
                    ++(cnt_trapped_att[ix][0]);
                }
            }
            ++(cnt[ix][0]);
        }
    }
    
    for ( Couple * cxi=couples.firstAF(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if ( cxi->active() ) ++act[ix];
            if (cxi->trapped())
            {
                ++(cnt_trapped[ix][1]);
                cxi->trapped_haMon->getHands(h1, dummy);
                if (h1->attached())
                {
                    ++(cnt_trapped_att[ix][1]);
                }
            }
            ++(cnt[ix][1]);
        }
    }
    for ( Couple * cxi=couples.firstFA(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if ( cxi->active() ) ++act[ix];
            if (cxi->trapped())
            {
                ++(cnt_trapped[ix][2]);
                cxi->trapped_haMon->getHands(h1, dummy);
                if (h1->attached())
                {
                    ++(cnt_trapped_att[ix][2]);
                }
            }
            ++(cnt[ix][2]);
        }
    }
    
    for ( Couple * cxi=couples.firstAA(); cxi ; cxi = cxi->next() )
    {
        unsigned ix = cxi->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if ( cxi->active() ) ++act[ix];
            if (cxi->trapped())
            {
                ++(cnt_trapped[ix][3]);
                cxi->trapped_haMon->getHands(h1, dummy);
                if (h1->attached())
                {
                    ++(cnt_trapped_att[ix][3]);
                }
            }
            ++(cnt[ix][3]);
        }
    }
    
    if ( 1 )
    {
        out << COM << ljust("couple", 2, 1);
        out << SEP << "total";
        out << SEP << "active";
        out << SEP << "FF";
        out << SEP << "AF";
        out << SEP << "FA";
        out << SEP << "AA";
    }
    
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        out << LIN << ljust((*ip)->name(), 2);
        unsigned ix = (*ip)->index();
        if ( ix < mx )
        {
            out << SEP << nb[ix];
            out << SEP << act[ix];
            for ( unsigned int d = 0; d < 4; ++d )
                out << SEP << cnt[ix][d];
        }
        else
            out << SEP << "out-of-range";

        out << LIN << ljust(" -->trapped", 2);
        ix = (*ip)->index();
        if ( ix < mx )
        {
            out << SEP << nb[ix];
            out << SEP << act[ix];
            for ( unsigned int d = 0; d < 4; ++d )
                out << SEP << cnt_trapped[ix][d];
        }
        else
            out << SEP << "out-of-range";
        
        out << LIN << ljust(" -->trap_attached", 2);
        ix = (*ip)->index();
        if ( ix < mx )
        {
            out << SEP << nb[ix];
            out << SEP << act[ix];
            for ( unsigned int d = 0; d < 4; ++d )
                out << SEP << cnt_trapped_att[ix][d];
        }
        else
            out << SEP << "out-of-range";
        
    }


}
void Simul::reportSingleTrap(std::ostream& out) const
{
    PropertyList plist = properties.find_all("single");
    
    const unsigned mx = 128;
    
    int nb[mx] = { 0 }, free[mx] = { 0 }, bound[mx] = { 0 },free_trapped[mx] = {0}, bound_trapped[mx] = {0};
    
    for ( Single * si = singles.firstF(); si ; si = si->next() )
    {
        unsigned ix = si->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            ++free[ix];
        }
    }
    
    for ( Single * si = singles.firstA(); si ; si = si->next() )
    {
        unsigned ix = si->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            ++bound[ix];
        }
    }
    for ( Single * si = singles.firstTrapped(); si ; si = si->next() )
    {
        unsigned ix = si->prop->index();
        if ( ix < mx )
        {
            ++nb[ix];
            if (si->attached())
                ++bound_trapped[ix];
            else
                ++free_trapped[ix];
        }
    }
    
    
    if ( 1 )
    {
        out << COM << ljust("single", 2, 1);
        out << SEP << "total";
        out << SEP << "free";
        out << SEP << "bound";
    }
    
    for ( PropertyList::iterator ip = plist.begin(); ip < plist.end(); ++ip )
    {
        out << LIN << ljust((*ip)->name(), 2);
        unsigned ix = (*ip)->index();
        if ( ix < mx )
        {
            out << SEP << nb[ix];
            out << SEP << free[ix]+free_trapped[ix];
            out << SEP << bound[ix]+bound_trapped[ix];
        }
        else
            out << SEP << " out-of-range ";
        out << LIN << ljust(" -->trapped", 2);
        ix = (*ip)->index();
        if ( ix < mx )
        {
            out << SEP << nb[ix];
            out << SEP << free_trapped[ix];
            out << SEP << bound_trapped[ix];
        }
        else
            out << SEP << " out-of-range ";
    }

}


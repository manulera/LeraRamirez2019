// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "sim.h"
#include <fstream>
#include <unistd.h>
#include "filepath.h"
#include "tictoc.h"
#include "parser.h"

/// Current format version number used for writing object-files.
/**
 This is related to Inputter::formatID()
 */
const int currentFormatID = 47;

//History of changes of currentFormatID:

/*
 47: 13/02/2017 Wrist and Aster built on a local reference frame of Solid
 46: 23/10/2015 GrowingFiber writes dynamic states, changed ClassicFiber:write()
 45: 18/09/2015 Indices of Properties are numbers starting from one, and not zero
 44: 25/07/2015 the dynamic state of fiber ends is stored with a separate TAG
 43: 24/04/2014 number of dimension in Space stored in 16 bits
 42: 09/11/2013 All fibers store end_state on 32 bits
     08/12/2012 FRAME_TAG was changed from "#frame " to "#Cytosim "
 41: 17/11/2012 Space stores its shape as a string in objects.cmo
 40: 11/09/2012 Aster format simplified
 39: 11/07/2012 Object::mark is stored on 32 bits instead of 16 previously
 38: 03/05/2012 Fiber stores length instead of segment-length
 37: 30/04/2012 Couple::Duo stores its activity state
 36: 22/02/2012 All Spaces store their dimensions in objects.cmo
 35: 15/09/2011 Some Spaces stores their dimensions in objects.cmo  
 34: 20/12/2010 Moved Fiber::mark to Object::mark
 33: 29/04/2010 Added Fiber::mark
 32: 15/04/2010 Space became an Object
 31: 01/04/2010 Fiber became a Mecable
 30: The Tag were reduced to 1 char: saves space & simplifies code
     26/05/2009 started cytosim-PI: a revolution!
 27: 22/03/2008 new Fiber::write(), called in Tubule::write()
 26: 03/11/2007 Hand do not record haEnd flag
 24: 14/12/2006 started cytosim 3, lots of changes
 23: 10/12/2005 new Solid
 22: modified Sphere
 21: modified Sphere
 20: 12/07/2004
 19: introduced different kind of Single
*/


//------------------------------------------------------------------------------
#pragma mark - Read Objects

/**
 We do not allow property()->index() of an Object to change during import from a file.
 However, there is no structural reason that prevent this in the code.
 If necessary, it should be possible to remove this limitation.
 
 The Object is not modified
 */
Object * Simul::readReference(Inputter & in, Tag & tag)
{
    tag = in.get_char();
    
    if ( tag == EOF )
        throw InvalidIO("unexpected end of file");

    // skip spaces in text mode:
    if ( !in.binary() )
    {
        while ( isspace(tag) || tag == 0 )
            tag = in.get_char();
    }
    
    // Object::TAG is the 'void' reference
    if ( tag == Object::TAG )
        return 0;
    
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 32 )
    {
        ObjectID n = isupper(tag) ? in.readUInt32() : in.readUInt16();
        if ( n == 0 )
            return 0;
        ObjectSet const* set = findSetT(tolower(tag));
        if ( set == 0 )
            throw InvalidIO("unknown object tag in Simul::find()");
        Object * w =  set->findID(n);
        if ( w == 0 )
            throw InvalidIO("Unknown object referenced (old style)");
        return w;
    }
#endif
    
    char pretag = 0;
    if ( tag == '$' )
    {
        pretag = tag;
        tag = in.get_char();
    }
    
    if ( !isalpha(tag) )
        throw InvalidIO("`"+std::string(1,tag)+"' is not a valid class TAG");

    const ObjectSet * set = findSetT(tag);
    
    if ( set == 0 )
        throw InvalidIO("`"+std::string(1,tag)+"' is not a known class TAG");
    
    unsigned ix = 0, mk = 0;
    ObjectID nb = 0;
    
    Object::readReference(in, ix, nb, mk, pretag);
    
    if ( nb == 0 )
        return 0;
    
    Object * res = set->findID(nb);
    
    if ( res == 0 )
    {
        throw InvalidIO("Unknown object `"+((char)tag+sMath::repr(nb))+"' referenced");
    }
    
    assert_true( res->identity() == nb );
    
    if ( res->property() == 0 || res->property()->index() != ix )
        throw InvalidIO("The property of a `"+res->property()->category()+"' should not change!");
    
    return res;
}


//------------------------------------------------------------------------------
#pragma mark - Read Objects

/// InputLock is a helper class used to import a cytosim state from a file
class Simul::InputLock
{
private:
    Simul * sim;
    bool    era;

public:
    
    /// change the mandate to 'erase' objects
    void erase(bool e) { era = e; }
    
    InputLock(Simul * s, bool e)
    : sim(s), era(e)
    {
        //MSG(21, "Simul::InputLock created with %i objects\n", sim->nbObjects());
        sim->couples.freeze();
        sim->singles.freeze();
        sim->fibers.freeze();
        sim->beads.freeze();
        sim->solids.freeze();
        sim->spheres.freeze();
        sim->organizers.freeze();
        sim->fields.freeze();
        sim->spaces.freeze();
        sim->events.freeze();
    }
    ~InputLock()
    {
        /*
         Attention: The order of the thaw() below is important:
         destroying a Fiber will detach any motor attached to it,
         and thus automatically move them to the 'unattached' list,
         as if they had been updated from reading the file.
         Destroying couples and singles before the fibers avoid this problem.
         */
        if ( era )
        {
            sim->events.prune();
            sim->organizers.prune();
            sim->couples.prune();
            sim->singles.prune();
            sim->beads.prune();
            sim->solids.prune();
            sim->spheres.prune();
            sim->fibers.prune();
            sim->spaces.prune();
            sim->fields.prune();
        }
        else
        {
            sim->events.thaw();
            sim->organizers.thaw();
            sim->couples.thaw();
            sim->singles.thaw();
            sim->beads.thaw();
            sim->solids.thaw();
            sim->spheres.thaw();
            sim->fibers.thaw();
            sim->spaces.thaw();
            sim->fields.thaw();
        }
        //MSG(21, "Simul::InputLock deleted with %i objects\n", sim->nbObjects());
    }
};


/**
 This will update the current state to make it identical to what has been saved in the file.
 
 Before reading, all objects are marked with flag().
 Every object encountered in the file is unmarked as it is updated.
 
 When the read is complete, the objects that are still marked are deleted.
 In this way the new state reflects exactly the system that was stored on file.
 
 @returns
 - 0 = success
 - 1 = EOF
 .
 */
int Simul::reloadObjects(Inputter & in)
{
    // set flag to erase any object that was not updated
    InputLock lock(this, true);

    // if an error occurs, we do not erase objects:
    if ( loadObjects(in) )
        lock.erase(false);

    return in.eof();
}


/**
 Read Objects from a file:
 update the ones that are present in the simulation world,
 and add the ones that are new.
 @returns
 - 0 = success
 - 1 = EOF
 .
 */
int Simul::loadObjects(Inputter & in, ObjectSet* subset)
{
    if ( in.eof() )
        return 1;
    
    if ( ! in.good() )
        throw InvalidIO("invalid file in Simul::loadObjects()");
    
    int res = 0;
    
    in.lock();
    try
    {
        res = readObjects(in, subset);
        //std::clog << "loadObjects returns " << res << std::endl;
    }
    catch(Exception & e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        in.unlock();
        throw;
    }
    
    in.unlock();
    return res;
}


/**
 Create an Inputter 'in' and call 'loadObjects(in)'
 */
int Simul::loadObjects(char const* filename)
{
    Inputter in(filename, true);
    in.vectorSize(DIM);

    if ( ! in.good() )
        throw InvalidIO("Could not open specified file for reading");
    
    return loadObjects(in);
}


//------------------------------------------------------------------------------

/**
 Read Objects from a file, and add new ones to the simulation
 The Inputter should be locked in a multithreaded application
 
 @returns
 - 0 : success
 - 1 : EOF
 - 2 : the file does not appear to be a valid cytosim archive
 
 if `subset` is not zero, only objects from this class will be imported
  */
int Simul::readObjects(Inputter & in, ObjectSet* subset)
{
    bool recognized = false;
    char c = 0, tag, pretag;
    ObjectSet * objset = 0;
    std::string section;
    std::string line;
    
    while ( in.good() )
    {
        // find the next entry: a newline followed by a ascii character
        line.clear();
        while ( c != '\n' )
        {
            c = in.get_char();
            if ( c == EOF )
                return 1;
            line.push_back(c);
        }
        
        if ( line.size() > 2 )
        {
            std::clog << "Skipped " << line.size() << " characters" << std::endl;
        }

        do
            tag = in.get_char();
        while ( tag == '\n' );

        if ( tag == EOF )
            return 1;
        
        //std::clog << "TAG |" << (char)tag << "|\n";
        
#ifdef BACKWARD_COMPATIBILITY
        // Compatibility with older format (before 2010)
        if ( in.formatID() < 32 )
        {
            ObjectSet * set = findSetT(tolower(tag));
            if ( set )
            {
                ObjectID n = isupper(tag) ? in.readUInt32() : in.readUInt16();
                if ( n == 0 )
                    throw InvalidIO("invalid (null) Object reference");
                Object * obj = set->findID(n);
                if ( obj )
                {
                    if ( tag!='i'  &&  ( tag!='m' || in.formatID()!=31 ))
                        in.readUInt16();
                    obj->read(in, *this, tag);
                    obj->flag(0);
                }
                else
                {
                    int pi = 0;
                    if ( tag!='i'  &&  ( tag!='m' || in.formatID()!=31 ))
                        pi = in.readUInt16();
                    obj = set->newObjectT(tolower(tag), pi);
                    obj->identity(n);
                    obj->read(in, *this, tag);
                    set->add(obj);
                }
                continue;
            }
        }
#endif
        
        if ( tag == '$' )
        {
            pretag = '$';
            tag = in.get_char();
            if ( tag == EOF )
                return 1;
        }
        else
            pretag = 0;

        if ( isalpha(tag) )
        {
            try
            {
                if ( objset )
                {
                    // check that we are in the correct ObjectSet:
                    assert_true( objset == findSetT(tag) );

                    Object * obj = objset->readObject(in, tag, pretag);
                    obj->flag(0);

                    if ( !obj->linked() )
                    {
                        if ( subset==0 || subset==objset )
                            objset->add(obj);
                        else
                            delete(obj);
                    }
                }
#ifdef BACKWARD_COMPATIBILITY
                else
                {
                    ObjectSet * set = findSetT(tag);
                    if ( set )
                    {
                        Object * obj = set->readObject(in, tag, pretag);
                        obj->flag(0);
                        
                        if ( !obj->linked() )
                            set->add(obj);
                    }
                }
#endif
            }
            catch( Exception & e )
            {
                if ( section.size() )
                {
                    std::cerr << "Error in section " << section << ": " << e.what() << std::endl;
                    if ( objset )
                        in.skip_until("#section ");
                }
                else
                    std::cerr << "Error : " << e.what() << std::endl;

            }
            continue;
        }

        //check for meta-data, contained in lines starting with '#'
        if ( tag == '#' )
        {
            in.get_line(line);
            //std::clog << "------|" << line << "|" << std::endl;
            
            //detect frame start
            if ( 0 == line.compare(0, strlen(FRAME_TAG)-1, FRAME_TAG) )
            {
                if ( recognized )
                    return 2;
                recognized = true;
                continue;
            }
            
#ifdef BACKWARD_COMPATIBILITY
            //detect frame start
            if ( 0 == line.compare(0, 6, "frame ") )
            {
                if ( recognized )
                    return 2;
                recognized = true;
                continue;
            }
#endif
            
            //detect section headers
            if ( 0 == line.compare(0, 8, "section ") )
            {
                std::string sub = line.substr(8);
                section = sub.substr(0, sub.find(' '));
                //std::clog << " section |" << section << "| " << std::endl;
#if ( 0 )
                // this skips loading Couple that are not bridging:
                if ( section == "couple" && sub != "couple AA" )
                {
                    in.skip_until("#section ");
                    std::clog << "skipped " << sub << std::endl;
                }
#endif
                objset = findSet(section);
                continue;
            }
            

            //detect info line
            if ( 0 == line.compare(0, 5, "time ") )
            {
                /* 
                 Parse the line, which should be formatted as follows:
                 #time 14.000000, dim 2, format 47
                 */
                std::stringstream ss(line);
                while ( ss.good() )
                {
                    std::string str;
                    ss >> str;
                    if ( str == "time" )
                        ss >> prop->time;
                    if ( str == "dim" )
                    {
                        int d;
                        ss >> d;
                        if ( ! ss.fail() && d != in.vectorSize() )
                        {
                            std::cerr << "Warning: Mismatch between file ("<<d<<"D) and executable ("<<DIM<<"D)\n";
                            in.vectorSize(d);
                        }
                    }
                    if ( str == "format" )
                    {
                        int f;
                        ss >> f;
                        if ( ! ss.fail() )
                        {
                            in.formatID(f);
                            if ( f != currentFormatID )
                                MSG(20, "Cytosim is reading data format %d\n", f);
                        }
                    }
                }
                continue;
            }
            
            //binary signature
            if ( 0 == line.compare(0, 7, "binary ") )
            {
                in.setEndianess(line.substr(7).c_str());
                continue;
            }
           
            //detect the mark at the end of the frame
            if ( 0 == line.compare(0, 10, "end frame ") )
                return 0;
            
            //detect the mark at the end of the frame
            if ( 0 == line.compare(0, 12, "end cytosim ") )
                return 0;
            
            continue;
        }
        else {
            //finally, we just skip the character
            if ( isprint(tag ) )
                std::cerr << "Warning: unknown TAG `" << tag << "'\n";
            else
                std::cerr << "Warning: unknown TAG \\" << (int)tag << "\n";
            c = tag;
        }
    }
    return 2;
}


//------------------------------------------------------------------------------
#pragma mark - Write Objects


void Simul::writeObjects(Outputter & out) const
{
    if ( ! out.good() )
        throw InvalidIO("output file is invalid");
    
    //std::clog << "Writing frame " << frame_index << std::endl;

    char date[26] = { 0 };
    TicToc::get_date(date, sizeof(date));
    
    // lock file:
    out.lock();
    // write a line identifying a new frame:
    fprintf(out, "\n\n#%s %i  %s", FRAME_TAG, getpid(), date);
    
    // record the simulated time:
    fprintf(out, "\n#time %.6f, dim %i, format %i", prop->time, DIM, currentFormatID);
    
    // identify the file as binary, with its endianess:
    if ( out.binary() )
        out.writeEndianess("\n#binary ");
    
    /*
     An object should be written after any other objects that it refers to.
     For example, Aster is written after Fiber, Couple after Fiber...
     This makes it easier to reconstruct the state during input.
     */
    
    spaces.write(out);
    fields.write(out);
    fibers.write(out);
    solids.write(out);
    beads.write(out);
    spheres.write(out);
    singles.write(out);
    couples.write(out);
    organizers.write(out);
    events.write(out);
    
    out.put_line("#section end");
    out.put_line("#end cytosim");
    fprintf(out, " %s\n\n", date);
    out.unlock();
}



/**
 This appends the current state to a trajectory file.
 Normally, this is objects.cmo in the current directory.
 
 If the file does not exist, it is created de novo.
*/
void Simul::writeObjects(char const* name, bool binary, bool append) const
{
    // use default file name if 'name' is empty or not provided
    if ( name == 0 || *name==0 )
        name = prop->trajectory_file.c_str();
    
    try
    {
        Outputter out(name, append, binary);
        writeObjects(out);
    }
    catch( InvalidIO & e )
    {
        std::cerr << "Error writing trajectory file: " << e.what() << '\n';
    }
}


//------------------------------------------------------------------------------
#pragma mark - Write Properties


/**
 The order of the output is important, since properties may depend
 on each other (eg. SingleProp and CoupleProp use HandProp).
 Luckily, there is no circular dependency in Cytosim at the moment.
 
 Thus we simply follow the order in which properties were defined,
 and which is the order in which properties appear in the PropertyList.
 */

void Simul::writeProperties(std::ostream& os, const bool prune) const
{
    //std::clog << "Writing properties" << std::endl;
    char date[26];
    TicToc::get_date(date, sizeof(date));
    os << "% Cytosim property file" << std::endl;
    os << "% " << date << std::endl;
    os << "% pid " << getpid() << '\n';

    prop->write(os, prune);
    properties.write(os, prune);
}


/**
 At the first call, this will write all properties to file, 
 and save a copy of what was written to a string `properties_saved`.
 
 The next time this is called, the properties will be compared to the string,
 and the file will be rewritten only if there is a difference.
 */
void Simul::writeProperties(char const* name, bool prune) const
{
    std::ostringstream oss;
    writeProperties(oss, prune);
    if ( oss.str() != properties_saved )
    {
        properties_saved = oss.str();

        // use default file name if 'name' is empty or not provided
        if ( name == 0 || *name==0 )
            name = prop->property_file.c_str();
        
        std::ofstream os(name);
        //this should be equivalent to: writeProperties(os, prune);
        os << properties_saved << std::endl;
        os.close();
        //std::clog << "Writing properties at frame " << currFrame() << std::endl;
    }
}


void Simul::loadProperties()
{
    if ( Parser(*this, 1, 1, 0, 0, 0).readConfig(prop->property_file) )
        std::cerr << "Error : File `" << prop->property_file << "' not found";
}

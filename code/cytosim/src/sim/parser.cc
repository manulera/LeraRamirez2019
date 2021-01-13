// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "sim.h"
#include "parser.h"
#include "glossary.h"
#include "tokenizer.h"
#include "filepath.h"
#include "stream_func.h"
#include <fstream>


// Use the second definition to get some verbose reports:
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG;


//------------------------------------------------------------------------------
/**
 The permission of the parser are:
 - allow_change: existing Property or Object can be modified
 - allow_set: new Properties can be created
 - allow_new: new Object can be created
 - allow_write: can write to disc
 .
 */
Parser::Parser(Simul& s, bool ds, bool dc, bool dn, bool dr, bool dw)
: Interface(s), do_set(ds), do_change(dc), do_new(dn), do_run(dr), do_write(dw)
{
    spos = 0;
}


//------------------------------------------------------------------------------
#pragma mark - Parse

/**
 Create a new Property, which is a set of parameters associated with a class.
 
 @code
 set CLASS NAME
 {
   PARAMETER =  VALUE
   ...
 }
 @endcode
 
 CLASS should be one of the predefined object (see @ref ObjectGroup).\n
 NAME should be a string, starting with a letter, followed by alphanumeric characters.
 The underscore character is also allowed, as in `fiber_23`\n
 
 It is also possible to use 'set' to change a parameter of an existing Property:
 
 @code
 set NAME
 {
   PARAMETER =  VALUE
   ...
 }
 @endcode

 or: 
 
 @code
 set NAME PARAMETER VALUE
 @endcode
 */

void Parser::parse_set(std::istream & is)
{
    std::string kind = Tokenizer::get_symbol(is);
    bool create_property = simul.isPropertyClass(kind);
    
#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    unsigned inx = 0;
    Tokenizer::get_integer(is, inx);
#endif

    std::string name = Tokenizer::get_symbol(is);

#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( is.peek() == ':' )
    {
        is.get();
        name = Tokenizer::get_symbol(is);
        if ( kind == "simul" )
        {
            kind = simul.prop->name();
            Tokenizer::get_token(is);
        }
        else
            kind = Tokenizer::get_token(is);
        create_property = false;
    }
#endif
    
    std::string blok = Tokenizer::get_token(is, true);

    if ( blok.empty() )
        throw InvalidSyntax("missing/empty value block");
    
    Glossary opt;
    Property * pp = 0;
    
    if ( create_property )
    {
        // in this form, 'set' defines a new Property:
        if ( do_set )
        {
            VLOG("+SET |" << kind << "|" << name << "|\n");
            opt.read(Tokenizer::strip_block(blok));
            pp = execute_set(kind, name, opt);
            
            int i;
            if ( opt.set(i, "property_index") )
            {
                if ( i != pp->index() )
                    throw InvalidSyntax("Property index missmatch");
            }
        }
    }
    else
    {
        if ( name.empty() )
            throw InvalidSyntax("unexpected syntax");

        // in this form, 'set' changes the value of an existing Property:
        opt.define(name, Tokenizer::strip_block(blok));
        
        if ( do_change )
        {
            VLOG(" SET |" << kind << "|" << name << "|\n");
            execute_change(kind, opt);
        }
        else if ( name == "display" )
        {
            pp = simul.findProperty(kind);
            if ( pp )
            {
                VLOG("-SET |" << kind << "|" << name << "|\n");
                pp->read(opt);
            }
        }
    }

    if ( pp && opt.warnings(std::cerr) )
    {
        std::cerr << "in\n";
        StreamFunc::show_lines(std::cerr, is, spos, is.tellg());
    }
}

//------------------------------------------------------------------------------
/**
 Change the value of one (or more) parameters for property `NAME`.

 @code
 change NAME
 {
   PARAMETER = VALUE
   ...
 }
 @endcode
 
 Short syntax:
 
 @code
 change NAME PARAMETER VALUE
 @endcode
 
The NAME should have been defined previously with the command `set`.
It is also possible to change all properties of a particular class:
 
 @code
 change all CLASS
 {
   PARAMETER = VALUE
   ...
 }
 @endcode

 Examples:
 
 @code
  change simul { viscosity = 0.5; }
  change simul display { back_color=red; }
  change actin { rigidity = 1; }
  change microtubule rigidity 20
  change all fiber { confine = inside, 10; }
  change all fiber display { color = white; }
 @endcode

 */

void Parser::parse_change(std::istream & is)
{
    bool change_all = false;
    
    std::string name = Tokenizer::get_symbol(is);
    
    if ( name == "all" )
    {
        change_all = true;
        name = Tokenizer::get_symbol(is);
    }
    
    std::string para = Tokenizer::get_symbol(is);

#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( is.peek() == ':' )
    {
        is.get();
        para = Tokenizer::get_symbol(is);
        std::string str = Tokenizer::get_token(is);
        if ( str == "*" )
            change_all = true;
        else
            name = str;
    }
    if ( is.peek() == '*' )
    {
        is.get();
        change_all = true;
    }
    Property * p = simul.findProperty(para);
    if ( p  &&  name == p->category() )
    {
        name = para;
        para = "";
    }
#endif

    std::string blok = Tokenizer::get_token(is, true);
    
    if ( blok.empty() )
        throw InvalidSyntax("missing/empty value block");
    
    if ( do_change )
    {
        Glossary opt;
        
        if ( para.size() > 0 )
            opt.define(para, Tokenizer::strip_block(blok));
        else
            opt.read(Tokenizer::strip_block(blok));
        
        int cnt = 1;
        if ( change_all )
        {
            VLOG("-CHG ALL |" << name << "|" << para << "|\n");
            cnt = execute_change_all(name, opt);
        }
        else
        {
            VLOG("-CHG |" << name << "|" << para << "|\n");
            execute_change(name, opt);
        }
        if ( opt.warnings(std::cerr, cnt) )
        {
            std::cerr << "in" << std::endl;
            StreamFunc::show_lines(std::cerr, is, spos, is.tellg());
        }
    }
}

//------------------------------------------------------------------------------
/**
 The command `new` creates one or more objects with given specifications:
 
 @code
 new [MULTIPLICITY] NAME
 {
   position         = POSITION, [SPACE]
   direction        = DIRECTION
   orientation      = ROTATION, [ROTATION]
   mark             = INTEGER
   required         = INTEGER
 }
 @endcode
 
 The NAME should have been defined previously with the command `set`.\n

 The other parameters are:
 
 Parameter        | type      | Default |   Description
 -----------------|-----------|---------|---------------------------------------------------------------------
 MULTIPLICITY     | INTEGER   |   1     | number of objects.
 `position`       | POSITION  | random  | initial position of the object.
 `orientation`    | ROTATION  | random  | a rotation specified with respect to the object's center of gravity.
 `orientation[1]` | ROTATION  | none    | a rotation specified around the origin.
 `direction`      | DIRECTION | random  | specifies the direction of a fiber.
 `mark`           | INTEGER   |   0     | specifies a mark to be given to all objects created.
 `required`       | INTEGER   |   0     | minimum number of objects that should be created.
 

 Note that `position` only applies to movable objects, and `orientation` will have an effect only on objects that can be rotated. In addition, `position[1]` and `orientation[1]` are relevant only if `(MULTIPLICITY > 1)`, and do not apply to the first object.\n
 
 
 Short syntax:
 
 @code
 new [MULTIPLICITY] NAME ( POSITION )
 @endcode

 Shorter syntax:
 
 @code
 new [MULTIPLICITY] NAME
 @endcode

*/

void Parser::parse_new(std::istream & is)
{
    Glossary opt;
    unsigned cnt = 1;
    Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
    
#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( simul.isPropertyClass(name) )
        name = Tokenizer::get_symbol(is);
#endif

    // Syntax sugar: () specify only position
    std::string blok = Tokenizer::get_block(is, '(');
    if ( blok.empty() )
    {
        blok = Tokenizer::get_block(is, '{');
        opt.read(blok);
    }
    else {
        opt.define("position", 0, blok);
    }
    
    if ( do_new  &&  cnt > 0 )
    {
        if ( opt.nb_keys() == 0 )
        {
            VLOG("-NEW |" << name << "|\n");
            execute_new(name, cnt);
        }
        else
        {
            VLOG("+NEW |" << name << "|\n");
            unsigned nb_objects = simul.nbObjects();
            
            // syntax sugar, to specify the position of the Fiber ends
            if ( opt.has_key("position_ends") )
            {
                Vector a, b;
                if ( opt.set(a, "position_ends") && opt.set(b, "position_ends", 1) )
                {
                    opt.define("length",      0, (a-b).norm());
                    opt.define("position",    0, (a+b)*0.5);
                    opt.define("orientation", 0, (b-a).normalized());
                }
            }

#ifdef BACKWARD_COMPATIBILITY
            // post_translation was replaced by keyword 'to'
            Vector vec;
            if ( opt.set(vec, "post_translation") )
            {
                opt.define("instance", 0, cnt-1);
                for ( unsigned n = 0; n < cnt; ++n )
                {
                    opt.define("instance", 1, n);
                    ObjectList objs = execute_new(name, opt);
                    ObjectSet::translateObjects(objs, n*vec);
                }
                opt.clear("instance");
            }
            else
#endif
            {
                opt.define("instance", 0, cnt-1);
                
                for ( unsigned n = 0; n < cnt; ++n )
                {
                    opt.define("instance", 1, n);
                    //std::clog << "instance " << n << " of " << cnt << ":\n";
                    ObjectList objs = execute_new(name, opt);
                }
                
                opt.clear("instance");
            }
            
            int created = simul.nbObjects() - nb_objects;
            
            int required = 0;
            if ( opt.set(required, "required")  &&  created < required )
            {
                std::cerr << "created  = " << created << std::endl;
                std::cerr << "required = " << required << std::endl;
                throw InvalidParameter("could not create enough `"+name+"'");
            }
            
            if ( opt.warnings(std::cerr, -1) )
            {
                std::cerr << "in" << std::endl;
                StreamFunc::show_lines(std::cerr, is, spos, is.tellg());
            }
        }
    }
}

//------------------------------------------------------------------------------
/**
 Delete objects:

 @code
 delete [MULTIPLICITY] NAME
 {
    mark      = INTEGER
    position  = [inside|outside], SPACE
    state     = [0|1], [0|1]
 }
 @endcode
 
 NAME can be '*', and the parameters (mark, position, state) are all optional.
 All specified conditions must be fulfilled (this is a logical AND).
 The parameter `state` refers to bound/unbound state of Hands for Single and Couple,
 and to dynamic state for Fibers:
 - for Single, `state[0]` refers to the Hand: 0=free, 1=attached.
 - for Couple, `state[0]` refers to the first Hand: 0=free, 1=attached,
           and `state[1]` refers to the second Hand: 0=free, 1=attached.
 - for Fibers, `state[0]` refers to the Dynanic state of the PLUS end,
           and `state[1]` refers to the Dynanic state of the MINUS end.
 .
 
 To delete all objects of specified NAME:
 @code
 delete NAME
 @endcode
 
 To delete at most CNT objects of class NAME:
 @code
 delete CNT NAME
 @endcode
 
 To delete all objects with a specified mark:
 @code
 delete NAME
 {
    mark = INTEGER
 }
 @endcode

 To delete all objects within a Space:
 @code
 delete NAME
 {
    position = inside, SPACE
 }
 @endcode
 
 The SPACE must be the name of an existing Space.
 Only 'inside' and 'outside' are valid specifications.

 To delete all Couple called NAME that are not bound:
 @code
 delete NAME { state = 0, 0; }
 @endcode

 To delete all Couple TYPE that are not crosslinking, use two calls:
 @code
 delete NAME { state1 = 0; }
 delete NAME { state2 = 0; }
 @endcode
*/

void Parser::parse_delete(std::istream & is)
{
    unsigned cnt = 0;
    bool has_cnt = Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( simul.isPropertyClass(name) )
    {
        name = Tokenizer::get_symbol(is);
    }
#endif
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_new )
    {
        Glossary opt(blok);

        if ( has_cnt )
            opt.define("nb_objects", 0, cnt);

        execute_delete(name, opt);
        
        if ( opt.warnings(std::cerr) )
        {
            std::cerr << "in" << std::endl;
            StreamFunc::show_lines(std::cerr, is, spos, is.tellg());
        }
    }
}


/**
 Mark objects:
 
 @code
 mark [MULTIPLICITY] NAME
 {
   mark       = INTEGER
   position   = POSITION
 }
 @endcode
 
 NAME can be '*', and the parameter position is optional.
 The syntax is the same as for command `delete`.
 */

void Parser::parse_mark(std::istream & is)
{
    unsigned cnt = 0;
    bool has_cnt = Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( simul.isPropertyClass(name) )
        name = Tokenizer::get_symbol(is);
#endif
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_new )
    {
        Glossary opt(blok);

        if ( has_cnt )
            opt.define("nb_objects", 0, cnt);

        execute_mark(name, opt);
        
        if ( opt.warnings(std::cerr) )
        {
            std::cerr << "in" << std::endl;
            StreamFunc::show_lines(std::cerr, is, spos, is.tellg());
        }
    }
}

//------------------------------------------------------------------------------
/**
 Cut all fibers that intersect a given plane.
 
 @code
 cut NAME
 {
    plane = VECTOR, REAL
 }
 @endcode
 
 NAME can be '*' to cut all fibers.
 The plane is specified by a normal vector `n` (VECTOR) and a scalar `a` (REAL).
 The plane is defined by <em> n.pos + a = 0 </em>
 */

void Parser::parse_cut(std::istream & is)
{    
    std::string name = Tokenizer::get_token(is);
#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( simul.isPropertyClass(name) )
        name = Tokenizer::get_symbol(is);
#endif
    std::string blok = Tokenizer::get_block(is, '{');

    if ( blok.empty() )
        throw InvalidSyntax("missing block after `cut'");
    
    if ( do_run )
    {
        Glossary opt(blok);
        execute_cut(name, opt);
    }
}


//------------------------------------------------------------------------------

/**
 @copydetails Interface::execute_run
 */
void Parser::parse_run(std::istream & is)
{
    unsigned cnt = 0;
    bool has_cnt = Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( name == "simul" )
        name = Tokenizer::get_token(is);
#endif
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_run )
    {
        Glossary opt(blok);

        if ( has_cnt )
        {
            if ( opt.has_key("nb_steps") )
                throw InvalidSyntax("the number of steps was specified twice");
            opt.define("nb_steps", 0, cnt);
        }
        
        execute_run(do_write, opt);
        
        if ( opt.warnings(std::cerr) )
        {
            std::cerr << "in" << std::endl;
            StreamFunc::show_lines(std::cerr, is, spos, is.tellg());
        }
    }
}

//------------------------------------------------------------------------------
/**
 Read and execute another config file.
 
 @code
 read FILE_NAME
 {
   required = BOOL
 }
 @endcode
 
 By default, `required = 1`, and execution will terminate if the file is not found.
 If `required=0`, the file will be executed if it is found, but execution will continue
 in any case.
 
 \todo: able to specify do_set and do_new for command 'include' 
*/

void Parser::parse_read(std::istream & is)
{
    bool required = true;
    std::string file = Tokenizer::get_filename(is);
    
    if ( file.empty() )
        throw InvalidSyntax("missing/invalid file name after 'read'");
    
    std::string blok = Tokenizer::get_block(is, '{');
    if ( ! blok.empty() )
    {
        Glossary opt(blok);
        opt.set(required, "required");
        if ( opt.warnings(std::cerr) )
        {
            std::cerr << "in" << std::endl;
            StreamFunc::show_lines(std::cerr, is, spos, is.tellg());
        }
    }
    
    std::ifstream fis(file.c_str());
    if ( ! fis.fail() )
    {
        VLOG("-READ "<< file << "\n");
        parse(fis, ", while reading `"+file+"'");
    }
    else
    {
        if ( required )
            throw InvalidSyntax("could not open file `"+file+"'");
        else
            MSG.warning("could not open file `%s'\n", file.c_str());
    }
}

//------------------------------------------------------------------------------
/**
 Import a simulation snapshot from a trajectory file
 
 @code
   import WHAT FILE_NAME
   {
     append = BOOL
     frame = INTEGER
   }
 @endcode
 
 The frame to be imported can be specified as an option: `frame=INTEGER`:
 @code
 import * my_file.cmo { frame = 10 }
 @endcode
 
 By default, this will replace the simulation state by the one loaded from file.
 To add the file objects to the simulation without deleting any of the current 
 object, you should specify `append = 1`:
 
 @code
 import * my_file.cmo { append = 1 }
 @endcode
 
 Finally instead of importing all the objects from the file, one can restrict
 the import to a desired class:
 
 @code
 import fiber my_file.cmo { append = 1 }
 @endcode
 
 Note that the simulation time will be changed to the one specified in the file,
 but this behavior can be changed by specifying the time:
 @code
 change simul * { time = 0 }
 @endcode
 */

void Parser::parse_import(std::istream & is)
{
    std::string what = Tokenizer::get_token(is);
    std::string file = Tokenizer::get_filename(is);
    
    if ( file.empty() )
        throw InvalidSyntax("missing/invalid file name after 'import'");
    
    std::string blok = Tokenizer::get_block(is, '{');
    Glossary opt(blok);
    
    if ( do_new )
    {
        execute_import(file, what, opt);
        if ( opt.warnings(std::cerr) )
        {
            std::cerr << "in" << std::endl;
            StreamFunc::show_lines(std::cerr, is, spos, is.tellg());
        }
    }
}

/**
 Export state to file. The general syntax is:
 
 @code
 export WHAT FILE_NAME
 {
   append = BOOL
   binary = BOOL
 }
 @endcode
 
 WHAT must be ``objects`` of ``properties``, and by default, both `binary` 
 and `append` are `true`. If `*` is specified instead of a file name,
 the current trajectory file will be used.
 
 
 Short syntax:
 @code
 export objects FILE_NAME
 @endcode
 
 
 Examples:
 
 @code
 export objects sim_objects.cmo { append=0 }
 export properties properties.txt
 @endcode
 
 Attention: this command is disabled for `play`.
 */

void Parser::parse_export(std::istream & is)
{
    std::string what = Tokenizer::get_token(is);
    std::string file = Tokenizer::get_filename(is);
    
    if ( file.empty() )
        throw InvalidSyntax("missing/invalid file name after 'export'");
    
    std::string blok = Tokenizer::get_block(is, '{');
    Glossary opt(blok);
    
    if ( do_write )
    {
        //what = Tokenizer::strip_block(what);
        execute_export(file, what, opt);
        if ( opt.warnings(std::cerr) )
        {
            std::cerr << "in" << std::endl;
            StreamFunc::show_lines(std::cerr, is, spos, is.tellg());
        }
    }
}


/**
 Export formatted data to file. The general syntax is:
 
 @code
 report WHAT FILE_NAME
 {
   append = BOOL
 }
 @endcode
 
 Short syntax:
 @code
 report WHAT FILE_NAME
 @endcode
 

 WHAT should be a valid argument to `report`:
 @copydetails Simul::report
 
 If `*` is specified instead of a file name, the report is sent to the standard output.
 
 Examples:
 
 @code
 report parameters parameters.cmo { append=0 }
 report fiber:length fibers_length.txt
 @endcode
 
 Note that this command is disabled for `play`.
 */

void Parser::parse_report(std::istream & is)
{
    std::string what = Tokenizer::get_symbols(is);
    std::string file = Tokenizer::get_filename(is);

    if ( file.empty() )
        throw InvalidSyntax("missing/invalid file name after 'report'");
    
    std::string blok = Tokenizer::get_block(is, '{');
    Glossary opt(blok);
    
    if ( do_write )
    {
        execute_report(file, what, opt);

        if ( opt.warnings(std::cerr) )
        {
            std::cerr << "in" << std::endl;
            StreamFunc::show_lines(std::cerr, is, spos, is.tellg());
        }
    }
}

//------------------------------------------------------------------------------

/**
 Call custom function
 
 @code
 call FUNCTION_NAME { OPTIONS }
 @endcode
 
 FUNCTION_NAME should be `equilibrate`, `custom0`, `custom1`, ... `custom9`.

 Note: The Simul::custom() functions need to be modified, to do something!
 */
void Parser::parse_call(std::istream & is)
{
    std::string str = Tokenizer::get_symbol(is);
    
    if ( str.empty() )
        throw InvalidSyntax("missing function name after 'call'");
    
    std::string blok = Tokenizer::get_block(is, '{');
    Glossary opt(blok);
    
    if ( do_run )
        execute_call(str, opt);
    
    if ( opt.warnings(std::cerr) )
    {
        std::cerr << "in" << std::endl;
        StreamFunc::show_lines(std::cerr, is, spos, is.tellg());
    }
}

//------------------------------------------------------------------------------
/**
 Repeat specified code.
 
 @code
 repeat INTEGER { CODE }
 @endcode
 
 */

void Parser::parse_repeat(std::istream & is)
{
    unsigned cnt = 1;
    
    if ( ! Tokenizer::get_integer(is, cnt) )
        throw InvalidSyntax("missing integer number after 'repeat'");

    std::string code = Tokenizer::get_block(is, '{');
    
    for ( unsigned c = 0; c < cnt; ++c )
    {
        //it is best to use a fresh stream for each instance:
        std::istringstream iss(code);
        parse(iss, ", while executing `repeat'");
    }
}


#ifdef NEW_COMMAND_FOR
/**
 Repeat specified code.
 
 @code
 for VAR=INTEGER:INTEGER { CODE }
 @endcode
 
 The two integers are the first and last iterations counters.
 Any occurence of VAR is replaced before the code is executed.
 
 Example:
 @code
 for CNT=1:10 {
   new 10 fiber tube { length = CNT }
 }
 @endcode
 
 NOTE: This code is a hack, and it might be removed in the future!
 */
void Parser::parse_for(std::istream & is)
{
    unsigned int start = 0, end = 1;
    
    std::string var = Tokenizer::get_symbol(is);
    
    std::string s = Tokenizer::get_token(is);
    if ( s != "=" )
        throw InvalidSyntax("missing '=' in command 'for'");
    
    if ( ! Tokenizer::get_integer(is, start) )
        throw InvalidSyntax("missing number after 'repeat'");

    s = Tokenizer::get_token(is);
    if ( s != ":" )
        throw InvalidSyntax("missing ':' in command 'for'");
    
    if ( ! Tokenizer::get_integer(is, end) )
        throw InvalidSyntax("missing number after 'repeat'");
    
    std::string code = Tokenizer::get_block(is, '{');
    
    for ( unsigned int c = start; c < end; ++c )
    {
        std::string blok = code;
        // substitute Variable name for this iteration:
        StreamFunc::find_and_replace(blok, var, sMath::repr(c));
        // it is best to use a fresh stream for each instance:
        std::istringstream iss(blok);
        parse(iss, ", while executing `for'");
        hold();
    }
}

#endif

//------------------------------------------------------------------------------

/**
 Terminates execution
 
 @code
 end
 @endcode
 */
void Parser::parse_end(std::istream & is)
{
    if ( do_run )
        throw Exception("terminating program at command 'end'");

    /*
    std::string str = Tokenizer::get_symbol(is);
    
    if ( str == "if" )
    {
        str = Tokenizer::get_token(is);
        ABORT_NOW("unfinished code");
    }
     */
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 The configuration file consist in a succession of commands, that mostly
 follow the following syntax:
 @code
 COMMAND CLASS NAME
 {
    PARAMETERS
 }
 @endcode
 
 Essential commands:
 
 Command        |   Effect
 ---------------|---------------------------------------------------------
 `set`          | Create a new Property, and set its parameter values
 `new`          | Create objects of a certain Property
 `run`          | Simulate the system to advance in time
 
 Additional commands:
 
 Command        |   Effect
 ---------------|---------------------------------------------------------
 `change`       | Change parameter values in an existing Property
 `read`         | Read another file and excutes the commands it contains
 `delete`       | Delete objects from the simulation
 `import`       | Import objects from a trajectory file
 `export`       | Export all simulated objects to a file
 `report`       | generate file or text with formatted information

 Other commands:
 
 Command        |   Effect
 ---------------|---------------------------------------------------------
 `mark`         | Mark objects
 `cut`          | Cut fibers along a plane
 `repeat`       | Execute code a number of times
 `for`          | Execute code a number of times (disabled)
 `end`          | Terminates simulation
 `call`         | Call a custom function

 */
void Parser::parse(std::istream & is, std::string const& msg)
{
    std::streampos fpos;
    std::string tok;
    
    try {
        while ( is.good() )
        {
            do {
                spos = is.tellg();
                fpos = spos;
                tok = Tokenizer::get_token(is);
                if ( is.fail() ) return;
            } while ( tok.length() < 1 || isspace(tok[0]) );
            
            
            // skip matlab-style comments
            if ( tok[0] == '%' )
            {
                if ( '{' == is.peek() )
                    Tokenizer::get_block_content(is, is.get(), '}');
                else
                    Tokenizer::get_line(is);
                continue;
            }

#ifdef BACKWARD_COMPATIBILITY
            /*
             skip C-style comments:
             - single-line comment start with '/' and '/'
             - multi-line comments start with '/' and '*'
            */
            if ( tok[0] == '/' )
            {
                int c = is.get();
                if ( '/' == c )
                    Tokenizer::get_line(is);
                else if ( '*' == c )
                    Tokenizer::get_until(is, "*/");
                else
                    throw InvalidSyntax("unexpected token after / `"+tok+"'");
                continue;
            }
#endif
            
            // StreamFunc::show_line(std::cout, is, fpos, "PARSE ");
            
            if ( tok == "set" )
                parse_set(is);
            else if ( tok == "change" )
                parse_change(is);
            else if ( tok == "new" )
                parse_new(is);
            else if ( tok == "delete" )
                parse_delete(is);
            else if ( tok == "mark" )
                parse_mark(is);
            else if ( tok == "run" )
                parse_run(is);
            else if ( tok == "read" )
                parse_read(is);
#ifdef BACKWARD_COMPATIBILITY
            else if ( tok == "include" )
                parse_read(is);
#endif
            else if ( tok == "cut" )
                parse_cut(is);
            else if ( tok == "report" )
                parse_report(is);
#ifdef BACKWARD_COMPATIBILITY
            else if ( tok == "write" )
                parse_report(is);
#endif
            else if ( tok == "import" )
                parse_import(is);
            else if ( tok == "export" )
                parse_export(is);
            else if ( tok == "call" )
                parse_call(is);
            else if ( tok == "repeat" )
                parse_repeat(is);
            else if ( tok == "skip" )
                Tokenizer::get_block(is, '{');
#ifdef NEW_COMMAND_FOR
            else if ( tok == "for" )
                parse_for(is);
#endif
            else if ( tok == "end" )
                return;//parse_end(is);
            else if ( tok == ";" )
                continue;
            else if ( tok == "dump" )
            {
                if ( do_write && do_run )
                    simul.dump();
            }
            else {
                throw InvalidSyntax("unknown command `"+tok+"'");
            }
            
            hold();
        }
    }
    catch( Exception & e )
    {
        e << msg + '\n';
        e << StreamFunc::get_lines(is, fpos, is.tellg());
        throw;
    }
}



int Parser::readConfig(std::string const& file)
{
    std::ifstream is(file.c_str());
    
    if ( ! is.good() )
        return 1;
    
    VLOG("-----------  Cytosim reads " << file << " (");
    VLOG(" set=" << do_set << "  change=" << do_change << "  new=" << do_new);
    VLOG("  run=" << do_run << "  write=" << do_write << " )\n");
    
    parse(is, ", while reading `"+file+"':");
    return 0;
}



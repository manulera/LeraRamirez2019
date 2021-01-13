// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "field_set.h"
#include "field_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"

//------------------------------------------------------------------------------

void FieldSet::prepare()
{
    for ( Field * f=first(); f; f=f->next() )
    {
        assert_true( f->hasField() );
        f->prepare();
    }
}


void FieldSet::step()
{
    for ( Field * f=first(); f; f=f->next() )
    {
        if ( f->hasField() )
        {
            PRINT_ONCE("!!!! Field is active\n");
            f->step(simul.fibers);
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

Property* FieldSet::newProperty(const std::string& kd, const std::string& nm, Glossary&) const
{
    if ( kd == "field" )
        return new FieldProp(nm);
    return 0;
}


Object * FieldSet::newObjectT(const Tag tag, unsigned idx)
{
    Field * obj = 0;
    if ( tag == Field::TAG )
    {
        FieldProp * p = simul.findProperty<FieldProp*>("field", idx);
        if ( p == 0 )
            throw InvalidIO("no field class defined with id "+sMath::repr(idx));
        obj = new Field(p);
        //the field is not initialized, because it should be set by FieldBase::read
    }
    return obj;
}


/**
 @ingroup NewObject

 Specify the initial value of the Field:
 
 @code
 new field NAME
 {
    value = 0
 }
 @endcode
 
 \todo: read the value of the field from a file, at initialization
 */
ObjectList FieldSet::newObjects(const std::string& name, Glossary& opt)
{
    Property * p = simul.properties.find_or_die("field", name);
    FieldProp * fp = static_cast<FieldProp*>(p);
        
    Field * obj = new Field(fp);
        
    // initialize field:
    obj->setField();
        
    // an initial concentration can be specified:
    Field::value_type val = 0;
    if ( opt.set(val, "value") || opt.set(val, "initial_value") )
    {
        std::string str;
        if ( opt.set(str, "value", 1) )
        {
            Space const* spc = simul.findSpace(str);
            if ( spc == 0 )
                spc = obj->prop->confine_space_ptr;
            obj->setConcentration(spc, val, 0);
        }
        else
        {
            obj->setConcentration(val);
        }
    }

    ObjectList res;
    res.push_back(obj);
    return res;
}

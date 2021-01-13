// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "assert_macro.h"
#include "nucleus.h"
#include "exceptions.h"
#include "sphere_prop.h"
#include "bundle_prop.h"
#include "point_exact.h"
#include "fiber_set.h"
#include "glossary.h"
#include "bundle.h"
#include "simul.h"
#include "meca.h"

extern Random RNG;


void Nucleus::step()
{
}


void Nucleus::setInteractions(Meca & meca) const
{
    Sphere * sph = sphere();
    
    if ( sph )
    {
        unsigned nix = sphere()->nbPoints() - Sphere::nbRefPts;
        
        for ( unsigned ix = 0; ix < nix; ++ix )
        {
            const Fiber * fib = fiber(ix);
            if ( fib )
                meca.interLink(PointExact(sph, ix+Sphere::nbRefPts),
                               fib->exactEnd(prop->focus),
                               prop->stiffness );
        }
    }
}



//------------------------------------------------------------------------------
ObjectList Nucleus::build(Glossary& opt, Simul& simul)
{
    assert_true(prop);
    ObjectList res;
    
    real rad = -1;
    if ( !opt.set(rad, "radius" ) || rad <= 0 )
        throw InvalidParameter("nucleus:radius should be specified and > 0");
   
    Sphere * sph = new Sphere(prop->sphere_prop, rad);
    grasp(sph, 0);
    res.push_back(sph);
    
    // get the center of the sphere
    Vector c = sph->posP(0);
    
    if ( prop->nb_fibers > 0 )
    {        
        // create points and clamps and add fiber attached to them
        for ( unsigned ii = 0; ii < prop->nb_fibers; ++ii )
        {
            ObjectList objs = simul.fibers.newObjects(prop->fibers, opt);
            if ( objs.size() )
            {
                Fiber * fib = Fiber::toFiber(objs[0]);
                Vector pos = c + Vector::randU(rad);
                Vector dir = Vector::randU();
                fib->setStraight(pos, dir, prop->focus);
                sph->addPoint(pos);
                res.append(objs);
                grasp(fib);
            }
        }
    }
    
    if ( prop->nb_bundles > 0 )
    {
        Rotation rotation;
        // add bundles        
        const real len = 0.5 * prop->bundle_prop->overlap;
        for ( unsigned int ii=0; ii < prop->nb_bundles; ++ii  )
        {
            rotation = Rotation::randomRotation(RNG);
            //a random position on the sphere:
            Vector pos = rotation * Vector(0,rad,0);
            //a direction tangent to the sphere:
            Vector dir = rotation * Vector(1,0,0);
            
            Bundle * bu = new Bundle(prop->bundle_prop);
            ObjectList objs = bu->build(opt, simul);
            res.append(objs);
            res.push_back(bu);
            
            //position the bundle correctly:
            bu->rotate(rotation);
            bu->translate(pos);
            
            sph->addPoint( c + (pos-len*dir).normalized(rad) );
            grasp(bu->organized(0));
            
            sph->addPoint( c + (pos+len*dir).normalized(rad) );
            grasp(bu->organized(1));
        }        
    }
    
    return res;
}

//------------------------------------------------------------------------------


/**
 This sets the ends of the link number `inx`
 or returns zero if the link does not exist
 */
unsigned Nucleus::getLink(unsigned inx, Vector& pos1, Vector& pos2) const
{
    if ( sphere() && ( inx + Sphere::nbRefPts < sphere()->nbPoints() ))
    {
        pos1 = sphere()->posP(inx+Sphere::nbRefPts);
        
        Fiber const* fib = Fiber::toFiber(organized(inx+1));
        if ( fib )
            pos2 = fib->posEnd(prop->focus);
        else
            pos2 = pos1;
        return 1;
    }
    
    return 0;
}


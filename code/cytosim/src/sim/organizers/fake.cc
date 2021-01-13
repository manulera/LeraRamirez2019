// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "glossary.h"
#include "messages.h"
#include "point_exact.h"
#include "simul.h"
#include "fake.h"
#include "aster.h"
#include "meca.h"



void Fake::step()
{
}


void Fake::setInteractions(Meca & meca) const
{
    assert_true( linked() );
    assert_true( asterPoints.size() == solidPoints.size() );
    
    for (unsigned n = 0; n < asterPoints.size(); ++n )
        meca.interLink(asterPoints[n], solidPoints[n], prop->stiffness);
}


void Fake::translate(Vector const& T)
{
    if ( fkSolid  &&  fkSolid->mobile() )   fkSolid->translate(T);
    if ( fkAster1 && fkAster1->mobile() )  fkAster1->translate(T);
    if ( fkAster2 && fkAster2->mobile() )  fkAster2->translate(T);
}


void Fake::rotate(Rotation const& T)
{
    if ( fkSolid  &&  fkSolid->mobile() )   fkSolid->rotate(T);
    if ( fkAster1 && fkAster1->mobile() )  fkAster1->rotate(T);
    if ( fkAster2 && fkAster2->mobile() )  fkAster2->rotate(T);
}


Aster * Fake::findAster(std::string const& spec, Simul& simul)
{
    ObjectList objs = simul.organizers.findObjects(spec);
    
    if ( objs.size() != 1 )
        throw InvalidParameter("could not find fake:aster1 `"+spec+"'");
    
    // we found an Aster already made
    if ( objs[0]->tag() != Aster::TAG )
        throw InvalidParameter("object `"+spec+"' is not an Aster");
    
    return static_cast<Aster*>(objs[0]);
}


ObjectList Fake::build(Glossary& opt, Simul& simul)
{
    ObjectList res;
    real rad = 0;
    if ( ! opt.set(rad, "radius") ||  rad <= 0 )
        throw InvalidParameter("fake:radius must be specified and > 0");

    assert_true(prop);
    
    // find the Aster specified:
    std::string spec;
    if ( !opt.set(spec, "aster1") )
        throw InvalidParameter("fake:aster1 must be specified");
    Aster * ax = findAster(spec, simul);
    
    if ( !opt.set(spec, "aster2") )
        throw InvalidParameter("fake:aster2 must be specified");
    Aster * bx = findAster(spec, simul);
   
    Solid * as = ax->solid();
    Solid * bs = bx->solid();
    Solid * so = new Solid(as->prop);
    
    Vector apos = ax->position();
    Vector bpos = bx->position();
    
    // define two orthogonal directions:
    Vector dir1, dir2;
#if ( DIM == 3 )
    Vector dir0 = ( apos - bpos ).normalized();
    dir0.orthonormal(dir1, dir2);
#else
    dir1 = ( apos - bpos ).orthogonal(1);
    dir2 = dir1;
#endif
    
    asterPoints.clear();
    solidPoints.clear();

    for ( int d = 1; d < DIM; ++d )
    {
        for ( int s = -1; s < 2; s += 2 )
        {
            Vector x = s * ( d == 2 ? dir2 : dir1 );
            solidPoints.push_back(PointExact(so, so->addSphere(apos+x, rad)));
            asterPoints.push_back(PointExact(as, as->addPoint(apos+x)));
        
            solidPoints.push_back(PointExact(so, so->addSphere(bpos+x, rad)));
            asterPoints.push_back(PointExact(bs, bs->addPoint(bpos+x)));
        }
    }
    
    so->fixShape();
    as->fixShape();
    bs->fixShape();

    /*
     // print for debugging:
     so->write(std::clog, true);
     as->write(std::clog, true);
     bs->write(std::clog, true);
    */
    
    grasp(so);
    res.push_back(so);

    // remember the Objects
    fkSolid = so;
    fkAster1 = ax;
    fkAster2 = bx;
    
    return res;
}

//------------------------------------------------------------------------------

PointDisp * Fake::disp() const
{
    if ( fkSolid )
        return fkSolid->prop->disp;
    return 0;
}


/**
 This sets the ends of the link number `inx`
 or returns zero if the link does not exist
 */
unsigned Fake::getLink(unsigned inx, Vector& pos1, Vector& pos2) const
{
    if ( inx < asterPoints.size() )
    {
        pos1 = asterPoints[inx].pos();
        pos2 = solidPoints[inx].pos();
        return 1;
    }
    
    return 0;
}


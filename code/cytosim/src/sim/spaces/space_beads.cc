// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_beads.h"
#include "bead_prop.h"
#include "object_set.h"
#include "simul.h"



/**
 The parameters BeadSet and BeadProp define the inner volume for this Space
*/
SpaceBeads::SpaceBeads(const SpaceProp* p)
: Space(p)
{
    for ( int d = 0; d < 3; ++d )
    {
        bbMin[d] = 0;
        bbMax[d] = 0;
    }
}

/**
 refresh the list of Beads 
 */
void SpaceBeads::resize()
{
    if ( objset() )
    {
        Simul const& sim = objset()->simul;
        BeadProp * bip = sim.findProperty<BeadProp*>("bead", prop->shape_spec);

        mBeads.clear();
        if ( bip == 0 )
        {
            for ( Bead * bd = sim.beads.first(); bd; bd=bd->next() )
                mBeads.push_back(bd);
        }
        else
        {
            for ( Bead * bd = sim.beads.first(); bd; bd=bd->next() )
            {
                if ( bd->property() == bip )
                    mBeads.push_back(bd);
            }
        }
#if ( 1 )
        static unsigned nb = 0;
        if ( nb != mBeads.size() )
        {
            nb = mBeads.size();
            std::clog << "SpaceBeads has " << nb << " beads\n";
        }
#endif
    }
}

/**
 Calculate a box in which the box are entirely inside
 */
void SpaceBeads::setBoundaries()
{
    BeadList::iterator bi = mBeads.begin();
    const BeadList::iterator end = mBeads.end();
    if ( bi )
    {
        real rad = (*bi)->radius();
        Vector pos = (*bi)->position();
        for ( int d = 0; d < DIM; ++d )
        {
            bbMin[d] = pos[d]-rad;
            bbMax[d] = pos[d]+rad;
        }
        while ( ++bi < end )
        {
            real rad = (*bi)->radius();
            Vector pos = (*bi)->position();
            for ( int d = 0; d < DIM; ++d )
            {
                if ( pos[d]-rad < bbMin[d] ) bbMin[d] = pos[d]-rad;
                if ( pos[d]+rad > bbMax[d] ) bbMax[d] = pos[d]+rad;
            }
        }
    }
    //std::clog << "SpaceBeads::bounding " << bbMin[0] << "  " << bbMax[0] << std::endl;
}


void SpaceBeads::step()
{
    resize();
    setBoundaries();
}


void SpaceBeads::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(bbMin[0], bbMin[1], bbMin[2]);
    sup.set(bbMax[0], bbMax[1], bbMax[2]);
}


/**
 Calculates the sum of the bead's volumes,
 thus assuming that the beads do no overlap
 */
real SpaceBeads::volume() const
{
    real res = 0;
    for ( BeadList::iterator oi = mBeads.begin(); oi < mBeads.end(); ++oi )
        res += (*oi)->volume();
    return res;
}


bool SpaceBeads::inside( const real point[] ) const
{
    for ( int d = 0; d < DIM; ++d )
    {
        if ( point[d] < bbMin[d] ) return false;
        if ( bbMax[d] < point[d] ) return false;
    }
    
    Vector pos(point);
    
    for ( BeadList::iterator oi = mBeads.begin(); oi < mBeads.end(); ++oi )
        if ( pos.distanceSqr((*oi)->position()) < (*oi)->radiusSqr() )
            return true;
    
    return false;
}


void SpaceBeads::project( const real point[], real proj[] ) const
{
    Vector pos(point);
    
    BeadList::iterator oi = mBeads.begin();
    Bead * closest = *oi;
    real dmin = fabs(pos.distance(closest->position()) - closest->radius());

    while ( ++oi < mBeads.end() )
    {
        real d = fabs(pos.distance((*oi)->position()) - (*oi)->radius());
        if ( d < dmin )
        {
            dmin    = d;
            closest = *oi;
        }
    }
    
    Vector bp = ( pos - closest->position() ).normalized(closest->radius());
    (closest->position()+bp).put(proj);
}


//------------------------------------------------------------------------------

void SpaceBeads::setInteraction(Vector const& pos, PointExact const& pe, Meca & meca, real stiff) const
{
    ABORT_NOW("SpaceBeads is incomplete");
}


void SpaceBeads::setInteraction(Vector const& pos, PointExact const& pe, real rad, Meca & meca, real stiff) const
{
    ABORT_NOW("SpaceBeads is incomplete");
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

bool SpaceBeads::display() const
{
    return true;
}



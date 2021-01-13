// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "point_set.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "point_interpolated.h"
#include "space.h"
#include "modulo.h"
#include "random.h"
#include "cblas.h"

extern Random RNG;

//------------------------------------------------------------------------------
void PointSet::psConstructor()
{
    psPos       = 0;
    psFor       = 0;
    psAllocated = 0;
}


PointSet::PointSet()
{
    psConstructor();
}

//------------------------------------------------------------------------------

PointSet::PointSet(const PointSet & o)
{
    psConstructor();
    allocatePoints( o.nbPoints() );
    pSize = o.pSize;
    for ( unsigned int p = 0; p < DIM*pSize; ++p )
        psPos[p] = o.psPos[p];
}


PointSet& PointSet::operator =(const PointSet& o)
{
    allocatePoints( o.nbPoints() );
    pSize = o.pSize;
    for ( unsigned int p = 0; p < DIM*pSize; ++p )
        psPos[p] = o.psPos[p];
    return *this;
}

//------------------------------------------------------------------------------
/** allocate(size) ensures that the set can hold `size` points
it returns the size if new memory was allocated 
*/
unsigned PointSet::allocatePoints(const unsigned nbp)
{
    if ( psAllocated < nbp )
    {
        // Keep memory aligned to 64 bytes:
        const unsigned chunk = 64 / sizeof(real);
        // make a multiple of chunk to align memory:
        const unsigned size = ( nbp + chunk - 1 ) & ~( chunk -1 );
        //std::clog << "PointSet::allocatePoints(" << nbp << ") allocates " << size << std::endl;

        real* mem = new real[DIM*size];
        
        if ( psPos )
        {
            //copy the current position to the new array
            for ( unsigned int p = 0; p < DIM*psAllocated; ++p )
                mem[p] = psPos[p];
            delete[] psPos;
        }
        psPos = mem;
        psAllocated = size;
        return size;
    }
    
    return 0;
}


void PointSet::deallocatePoints()
{
    if ( psPos )
    {
        delete[] psPos;
        psPos = 0;
    }
    psFor = 0;
    psAllocated = 0;
    pSize = 0;
}

//------------------------------------------------------------------------------
#pragma mark - Modifying points

unsigned PointSet::addPoint( Vector const& vec )
{
    allocatePoints(pSize+1);
    unsigned p = pSize++;
    vec.put(psPos+DIM*p);
    return p;
}


void PointSet::removePoints(const unsigned inx, const unsigned nbp)
{
    assert_true( inx + nbp <= pSize );
    
    pSize -= nbp;
    
    //move part of the array down, to erase 'nbp' points from index 'inx'
    for ( unsigned ii = DIM*inx; ii < DIM*pSize; ++ii )
        psPos[ii] = psPos[ii+DIM*nbp];
}


void PointSet::shiftPoints(const unsigned inx, const unsigned nbp)
{
    allocatePoints(pSize+nbp);
    
    //move part of the array up, making space for 'nbp' points from index 'inx'
    for ( unsigned ii = DIM*inx; ii < DIM*pSize; ++ii )
        psPos[ii+DIM*nbp] = psPos[ii];

    pSize += nbp;
}

//------------------------------------------------------------------------------
/**
 shifts ending-part of the array to indices starting at 0.
*/
void PointSet::truncateM(const unsigned int p)
{
    assert_true( p < pSize - 1 );

    unsigned int np = pSize - p;
    
    for ( unsigned int ii = 0; ii < DIM*np; ++ii )
        psPos[ii] = psPos[ii+DIM*p];
    
    pSize = np;
}

/**
 erase higher indices of array
*/
void PointSet::truncateP(const unsigned int p)
{
    assert_true( p < pSize );
    assert_true( p > 0 );
    
    pSize = p+1;
}

//------------------------------------------------------------------------------

void PointSet::resetPoints()
{
    if ( psPos )
    {
        for ( unsigned int p = 0; p < DIM*psAllocated; ++p )
            psPos[p] = 0;
    }
}


void PointSet::addNoise(const real amount)
{
    for ( unsigned int p = 0; p < DIM*pSize; ++p )
        psPos[p] += amount * RNG.sreal();
}


void PointSet::translate(Vector const& T)
{
    for ( unsigned p = 0; p < DIM*pSize; p += DIM )
    {
        psPos[p  ] += T.XX;
#if ( DIM > 1 )
        psPos[p+1] += T.YY;
#endif
#if ( DIM > 2 )
        psPos[p+2] += T.ZZ;
#endif
    }
}


void PointSet::rotate(Rotation const& T)
{
    for ( unsigned p = 0; p < DIM*pSize; p += DIM)
        T.vecMul( &psPos[p] );
}


//------------------------------------------------------------------------------
#pragma mark - Export/Inport


void PointSet::putPoints(real * x) const
{
    blas_xcopy(DIM*nbPoints(), psPos, 1, x, 1);
}


void PointSet::getPoints(const real * x)
{
    blas_xcopy(DIM*nbPoints(), x, 1, psPos, 1);
}


Vector PointSet::netForce(const unsigned p) const
{
    assert_true( p < pSize );
    if ( psFor )
        return Vector(psFor+DIM*p);
    else
        return Vector(0,0,0);
}

//------------------------------------------------------------------------------
/**
 Returns the center of gravity of all points
 */
Vector PointSet::position() const
{
    Vector sum(0,0,0);
    for ( unsigned p = 0; p < pSize; ++p )
        sum += posP(p);
    if ( pSize > 1 )
        return sum / real(pSize);
    return sum;
}

/**
 Calculate first and second moment of point distribution:
 - avg = sum( P ) / nb_points
 - sec = sum( P * P );
 .
 if 'sub = true', the average is substracted from 'sec'
 */
void PointSet::calculateMomentum(Vector& avg, Vector& sec, bool sub)
{
    avg.zero();
    sec.zero();
    
    // calculate first and second moments:
    for ( unsigned p = 0; p < pSize; ++p )
    {
        avg += posP(p);
        sec += posP(p).e_squared();
/*
        real const* pp = psPos + DIM*p;
        avg.XX += pp[0];
        sec.XX += pp[0] * pp[0];
#if ( DIM > 1 )
        avg.YY += pp[1];
        sec.YY += pp[1] * pp[1];
#endif
#if ( DIM > 2 )
        avg.ZZ += pp[2];
        sec.ZZ += pp[2] * pp[2];
#endif
 */
    }
    
    if ( pSize > 1 )
    {
        avg /= pSize;
        sec /= pSize;
    }
    
    if ( sub )
        sec -= avg.e_squared();
}


void PointSet::foldPosition(const Modulo * s)
{
    Vector off = s->offset(position());
    if ( !off.null() )
        translate(-off);
}


bool PointSet::allInside(const Space * spc) const
{
    bool res = true;
    for ( unsigned ii = 0; ii < nbPoints(); ++ii )
    {
        if ( spc->outside(posP(ii)) )
        {
            res = false;
            break;
        }
    }
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Read/write



void PointSet::write(Outputter& out) const
{
    out.writeUInt16(pSize);
    for ( unsigned int p = 0; p < pSize ; ++p )
        out.writeFloatVector(psPos+DIM*p, DIM, '\n');
}



void PointSet::read(Inputter& in, Simul&, Tag)
{
    try
    {
        unsigned nb = in.readUInt16();
        allocatePoints(nb);
    
        //we reset the point for a clean start:
        resetPoints();
        
        pSize = nb;
#if ( 1 )
        for ( unsigned int p = 0; p < nb ; ++p )
            in.readFloatVector(psPos+DIM*p, DIM);
#else
        in.readFloatVector(psPos, nb, DIM);
#endif
    }
    catch( Exception & e ) {
        e << ", in PointSet::read()";
        clearPoints();
        throw;
    }
}


void PointSet::write(std::ostream& os) const
{
    os << "new mecable " << reference() << "\n{\n";
    os << " nb_points = " << nbPoints() << std::endl;
    for ( unsigned pp = 0; pp < nbPoints() ; ++pp )
    {
        os << " point" << pp << " = " << posP(pp) << std::endl;
    }
    os << "}" << std::endl;
}


std::ostream& operator << (std::ostream& os, PointSet const& obj)
{
    obj.write(os);
    return os;
}



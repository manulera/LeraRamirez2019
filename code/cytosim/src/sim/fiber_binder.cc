// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "fiber_binder.h"
#include "fiber_locus.h"
#include "iowrapper.h"
#include "simul.h"
#include "sim.h"


FiberBinder::FiberBinder(Fiber* f, real a)
: fbFiber(f), fbAbs(a)
{
    assert_true(f);
    inter = f->interpolate(a);
}



void FiberBinder::locate(Fiber * f, real a)
{
    assert_true(f);
    assert_true(fbFiber==0);
    assert_true(f->abscissaM() <= a + REAL_EPSILON);
    assert_true(a <= f->abscissaP() + REAL_EPSILON);
    
    fbAbs   = a;
    fbFiber = f;
    
    f->addBinder(this);
    updateBinder();
}


void FiberBinder::delocate()
{
    assert_true( fbFiber );
    fbFiber->removeBinder(this);
    fbFiber = 0;
}


void FiberBinder::relocate(Fiber* f)
{
    if ( f != fbFiber )
    {
        if ( fbFiber )
            fbFiber->removeBinder(this);
        fbFiber = f;
        f->addBinder(this);
        updateBinder();
    }
}


void FiberBinder::relocate(real a)
{
    fbAbs = a;
    updateBinder();
}


void FiberBinder::relocate(Fiber* f, real a)
{
    if ( f != fbFiber )
    {
        if ( fbFiber )
            fbFiber->removeBinder(this);
        fbFiber = f;
        f->addBinder(this);
    }
    fbAbs = a;
    updateBinder();
}


void FiberBinder::moveToEndM()
{
    assert_true(fbFiber);
    fbAbs = fbFiber->abscissaM();
    inter = fbFiber->interpolateEndM();
}


void FiberBinder::moveToEndP()
{
    assert_true(fbFiber);
    fbAbs = fbFiber->abscissaP();
    inter = fbFiber->interpolateEndP();
}


void FiberBinder::moveToEnd(const FiberEnd end)
{
    assert_true(fbFiber);
    assert_true(end==PLUS_END || end==MINUS_END);
    
    if ( end == PLUS_END )
        moveToEndP();
    else
        moveToEndM();
}

//------------------------------------------------------------------------------
#pragma mark -


FiberEnd FiberBinder::nearestEnd() const
{
    assert_true(fbFiber);
    if ( fbAbs > fbFiber->abscissaC() )
        return PLUS_END;
    else
        return MINUS_END;
}


real  FiberBinder::abscissaFrom(const FiberEnd from) const
{
    assert_true(fbFiber);
    switch( from )
    {
        case MINUS_END:  return fbAbs - fbFiber->abscissaM();
        case PLUS_END:   return fbFiber->abscissaP() - fbAbs;
        case ORIGIN:     return fbAbs;
        case CENTER:     return fbAbs - fbFiber->abscissaC();
        default:         ABORT_NOW("invalid argument value");
    }
    return 0;
}


//------------------------------------------------------------------------------
#pragma mark -

void FiberBinder::write(Outputter& out) const
{
    out.writeSoftSpace();
    if ( fbFiber )
    {
        checkAbscissa();
        fbFiber->writeReference(out);
        out.writeFloat(fbAbs);
    }
    else {
        Object::writeNullReference(out);
    }
}


void FiberBinder::read(Inputter& in, Simul& sim)
{
    Tag tag = 0;
    Object * w = sim.readReference(in, tag);

    if ( w )
    {
        //std::clog << "FiberBinder::read() " << (char)tag << std::endl;
        Fiber * newfib = static_cast<Fiber*>(w);

        if ( tag == Fiber::TAG )
        {
            fbAbs  = in.readFloat();
        }
        else if ( tag == Fiber::TAG_LATTICE )
        {
            fbAbs  = in.readFloat();
            //skip information of site:
            in.readUInt32();
        }
#ifdef BACKWARD_COMPATIBILITY
        else if ( tag == 'm' )
        {
            fbAbs  = in.readFloat();
        } 
#endif
        else
        {
            ///\todo: we should allow binder to refer to any Mecable
            throw InvalidIO("FiberBinder should be bound to a Fiber!");
        }
        
        // link the FiberBinder as in attach():
        if ( newfib != fbFiber )
        {
            if ( fbFiber )
                fbFiber->removeBinder(this);
            fbFiber = newfib;
            fbFiber->addBinder(this);
        }
        updateBinder();
        checkAbscissa();
    }
    else
    {
        if ( fbFiber )
            FiberBinder::delocate();
    }
}

//------------------------------------------------------------------------------
#pragma mark -


std::ostream& operator << (std::ostream& os, FiberBinder const& obj)
{
    if ( obj.fiber() )
        os << "(" << obj.fiber()->reference() << " abs " << obj.abscissa() << ")";
    else
        os << "(null)";
    
    return os;
}


void FiberBinder::checkAbscissa() const
{
    assert_true(fbFiber);
    if ( fbAbs < fbFiber->abscissaM() - 1e-3 )
        MSG.warning("FiberBinder:abscissa < fiber:abscissa(MINUS_END) :  %e\n", fbFiber->abscissaM()-fbAbs );
    
    if ( fbAbs > fbFiber->abscissaP() + 1e-3 )
        MSG.warning("FiberBinder:abscissa > fiber:abscissa(PLUS_END)  :  %e\n", fbAbs-fbFiber->abscissaP() );
}


int FiberBinder::bad() const
{
    if ( fbFiber != inter.mecable() )
    {
        std::cerr << "Interpolation mismatch " << fbFiber << " " << inter.mecable() << std::endl;
        return 7;
    }
    
    if ( fbFiber->betweenMP(fbAbs) )
    {
        const real e = fbAbs - abscissaInter();
        
        //std::clog << "Interpolation " << std::scientific << e << std::endl;
        if ( fabs(e) > 1e-6 )
        {
            std::cerr << "Interpolation error is " << std::scientific << e << "\n";
            std::cerr << " abscissa:\n";
            std::cerr << "    binder       " << fbAbs << "\n";
            std::cerr << "    interpolated " << abscissaInter() << "\n";
            PointInterpolated pi = fbFiber->interpolate(fbAbs);
            std::cerr << "    updated      " << fbFiber->abscissaPoint(pi.point1()+pi.coef1()) << "\n";
            return 8;
        }
    }
    return 0;
}



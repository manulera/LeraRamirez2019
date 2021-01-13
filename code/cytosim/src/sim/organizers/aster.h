// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ASTER_H
#define ASTER_H

#include "object.h"
#include "organizer.h"
#include "aster_prop.h"
#include "solid.h"
#include "fiber.h"

class Glossary;

class AsterLink
{
    friend class Aster;
    
private:
    
    /// index of reference point
    unsigned ref;
    
    /// interpolation coefficient for Fiber end
    real coef1[4];
    
    /// interpolation coefficient for Fiber side
    real coef2[4];
    
    /// distance between the two anchoring points
    real len;
    
    /// true if the interpolation corresponds exactly to point 'ref'
    bool degenerate;
    
    /// only used for backward compatibility
    unsigned alt;
    
public:
    
    /// constructor
    AsterLink()
    {
        reset();
    }
    
    void reset()
    {
        ref = 0;
        len = 0;
        for ( int i = 0; i < 4; ++i )
        {
            coef1[i]=0.0;
            coef2[i]=0.0;
        }
        degenerate = 0;
        alt = 0;
    }
    
    void set(Vector const& A, Vector const& B)
    {
        len = ( A - B ).norm();
        
        coef1[1] = A.XX;
        coef2[1] = B.XX;
#if ( DIM == 1 )
        coef1[2] = 0;
        coef2[2] = 0;
        coef1[3] = 0;
        coef2[3] = 0;
        coef1[0] = 1.0 - A.XX;
        coef2[0] = 1.0 - B.XX;
#elif ( DIM == 2 )
        coef1[2] = A.YY;
        coef2[2] = B.YY;
        coef1[3] = 0;
        coef2[3] = 0;
        coef1[0] = 1.0 - A.XX - A.YY;
        coef2[0] = 1.0 - B.XX - B.YY;
#elif ( DIM == 3 )
        coef1[2] = A.YY;
        coef2[2] = B.YY;
        coef1[3] = A.ZZ;
        coef2[3] = B.ZZ;
        coef1[0] = 1.0 - A.XX - A.YY - A.ZZ;
        coef2[0] = 1.0 - B.XX - B.YY - B.ZZ;
#endif
        degenerate = ( fabs(coef1[1]) + fabs(coef1[2]) + fabs(coef1[2]) < REAL_EPSILON);
        alt = 0;
    }
    
    void write(Outputter& out) const
    {
        out.writeUInt16(ref);
        for ( int d = 1; d < 4; ++d )
            out.writeFloat(coef1[d]);
        for ( int d = 1; d < 4; ++d )
            out.writeFloat(coef2[d]);
    }
    
    void read(Inputter& in)
    {
        ref = in.readUInt16();
        for ( int d = 1; d < 4; ++d )
            coef1[d] = in.readFloat();
        for ( int d = 1; d < 4; ++d )
            coef2[d] = in.readFloat();
        coef1[0] = 1.0 - coef1[1] - coef1[2] - coef1[3];
        coef2[0] = 1.0 - coef2[1] - coef2[2] - coef2[3];
        len = (Vector3(coef1+1)-Vector3(coef2+1)).norm();
    }
    
    void write(std::ostream& out) const
    {
        const unsigned w = 12;
        out << std::setw(w) << coef1[0];
        for ( int d = 1; d < 4; ++d )
            out << " " << std::setw(w) << coef1[d];
        out << "     " << std::setw(w) << coef2[0];
        for ( int d = 1; d < 4; ++d )
            out << " " << std::setw(w) << coef2[d];
        out << "\n";
    }
};


/// A radial configuration of Fiber(s) built around a Solid
/**
 The parameters are defined in AsterProp.
 
 Each Fiber is attached to the Solid:
 - at the end of the Fiber
 - at a secondary point that is tied to the Fiber at some distance from this end.
 .
 This anchors the Fiber to the Solid, both in position and direction.
 The stiffness of the links is defined in AsterProp::stiffness, and can be adjusted independently.
 .
 
 */
class Aster : public Organizer
{
private:
    
    /// scale of local reference frame
    real       asRadius;
    
    /// store the coefficients needed to make the links between Solid and Fiber
    Array<AsterLink> asLinks;

    /// create and configure the Solid
    ObjectList makeSolid(Glossary& opt, Simul&, unsigned& origin);

    /// create a Fiber for position 'inx'
    ObjectList makeFiber(unsigned inx, Glossary& opt, Simul&);

    /// define the anchor points of Fibers
    void       placeFibers(Glossary& opt, unsigned origin);

    /// define the attachment position of fiber 'inx'
    void       anchorFiber(unsigned inx, Vector const&, Vector const&, unsigned origin);
    
    /// Property
    AsterProp const* prop;
    
public:
    
    /// constructor
    Aster(AsterProp const* p) : prop(p) { asRadius = 0; }
    
    /// destructor
    virtual      ~Aster();
    
    /// construct all the dependent Objects of the Organizer
    ObjectList    build(Glossary&, Simul&);
    
    /// return the scaffolding Solid
    Solid *       solid() const { return static_cast<Solid*>(organized(0)); }
    
    /// return the center of the Solid
    Vector        position() const { return solid()->posP(0); }
    
    /// return Fiber `n`
    Fiber *       fiber(unsigned n) const { return Fiber::toFiber(organized(n+1)); }
    
    /// perform one Monte-Carlo step
    void          step();
    
    /// add interactions to the Meca
    void          setInteractions(Meca &) const;
    

    
    /// position of first clamp for Fiber n
    Vector        posLink1(unsigned n) const;
    
    /// position of second clamp for Fiber n
    Vector        posLink2(unsigned n) const;

    /// position of attachment point on Fiber corresponding to second link
    Vector        posFiber2(unsigned n) const;
    
    /// obtain positions of a link
    unsigned      getLink(unsigned, Vector&, Vector&) const;
    
    /// a unique character identifying the class
    static const Tag TAG = 'a';
    
    /// return unique character identifying the class
    Tag           tag() const { return TAG; }
    
    /// return Object Property
    Property const* property() const { return prop; }

    /// read from IO
    void          read(Inputter&, Simul&, Tag);
    
    /// write to IO
    void          write(Outputter&) const;
    
    /// return PointDisp of Solid
    PointDisp *   disp() const { if ( solid() ) return solid()->prop->disp; return 0; }
  
};


#endif


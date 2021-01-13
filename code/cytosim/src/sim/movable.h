// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MOVABLE_H
#define MOVABLE_H

#include "vector.h"
#include "isometry.h"

class Space;
class Modulo;
class Glossary;
class Simul;


/// Can be moved and rotated in space
/**
This provides a common interface, to move and rotation object in space.
The actual operations need to be implemented by redefining the virtual functions:
 
 - mobile()  should return true
 - position() must be implemented
 - translate() must be implemented
 - rotate() may be implemented if the object has an orientation
 .
 To support periodic boundary conditions, foldPosition() should be defined.
 */
class Movable
{
    
public:
    
    /// read a position specified with primitives, such as 'circle 5', etc.
    static Vector readPrimitive(std::istream&, const Space*);
    
    /// read a position in space
    static Vector readPosition(std::istream&, const Space*);

    /// read an orientation, and return a normalized vector
    static Vector readDirection(std::istream&, const Vector&, const Space*);

    /// read a rotation specified in stream, at position `pos`
    static Rotation readRotation(std::istream&, const Vector&, const Space*);
    
public:
    
    /// constructor
    Movable() {}
    
    /// destructor
    virtual ~Movable() {}
    
    
    /// true if object accepts translations (default=false)
    virtual bool      mobile()  const { return false; }
    
    /// return the position in space of the object
    virtual Vector    position()  const { return Vector(0,0,0); }
    
    /// move object to specified position
    virtual void      setPosition(Vector const&);
    
    /// move the object ( position += given vector )
    virtual void      translate(Vector const&);
    
    /// rotate the object around the origin of coordinates
    virtual void      rotate(Rotation const&);
    
    /// rotate the object around its current position
    virtual void      revolve(Rotation const&);
    
    
    /// perform modulo for periodic boundary conditions
    /** This brings the object to the centered mirror image defined by Modulo */
    virtual void      foldPosition(Modulo const*) {}
    
};


#endif


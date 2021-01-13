// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MODULO_H
#define MODULO_H

#include "real.h"
#include "vector.h"

/// Modulo implements periodic boundary conditions
/**
 This class is used to apply periodic boundaries conditions in one or 
 in multiple dimensions in space.

 The methods are not virtual to avoid the calling overload in C++.
 */
class Modulo
{
private:
    
    /// bitfield indicating the dimensions that are periodic
    int   mMode;

    /// half-period in each dimension
    real  mSize[DIM];
    
    /// adjust 'x' to canonical representation within periodicity 'p':
    static void fold(real& x, const real p)
    {
        while ( x >  p ) x -= p+p;
        while ( x < -p ) x += p+p;
    }

public:
    
    /// constructor
    Modulo() { mMode = 0; }

    /// destructor
    ~Modulo() {}
    
    /// disable periodicity in all dimensions
    void disable() { mMode = 0; }
    
    /// enable periodicity in dimension 'd'
    void enable(int d, real size) { mMode |= 1<<d; mSize[d] = size; }
    
    /// true if direction `d` has periodic boundaries
    bool isPeriodic() const { return mMode; }

    /// true if direction `d` has periodic boundaries
    bool isPeriodic(int d) const { return mMode & (1<<d); }
    
    /// set vector `off` as the i-th direction of periodicity
    const Vector periodicity(int d) const;
    
    /// set `pos` to its periodic representation closest to the origin
    void         fold(Vector & pos) const;
    
    /// remove periodic repeats in `pos`, to bring it closest to `ref`
    void         fold(Vector & pos, const Vector& ref) const;
    
    /// calculate translation necessary to bring `pos` closest to origin
    const Vector offset(const Vector& pos) const;
    
    /// calculate translation necessary to bring `pos` closest to origin
    const Vector foldOffset(Vector& pos) const;
    
};

#endif



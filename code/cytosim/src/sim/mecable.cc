// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "mecable.h"
#include "exceptions.h"
#include "cblas.h"
#include "iowrapper.h"
#include "organizer.h"


/**
 mIndex is not initialized, since it must be set by Meca
 */
Mecable::Mecable() :
pSize(0), pBlock(0), pBlockAllocated(0), pBlockSize(0), pBlockUse(false)
{
}

/**
 
 The first time, we allocate exactly what is demanded, but if new allocation 
 is required again, we allocate with some margin, because it means that this
 object is probably growing.

 */
bool Mecable::allocateBlock(unsigned size)
{
    const size_t chunk = 64 / sizeof(real);
    
    //std::clog << "Mecable::allocateBlock " << size << "  " << pBlockSize << "\n";
    
    pBlockSize = size;
    
    if ( size > pBlockAllocated )
    {
        if ( pBlock )
            delete[] pBlock;
 
        size_t bum = ( size + chunk - 1 ) & ~( chunk -1 );
        
        pBlock = new real[ bum * bum ];
        pBlockAllocated = bum;
        
        //blas_xzero(bum*bum, pBlock);
       
        return true;
    }
    return false;
}


Mecable::~Mecable()
{
    if ( pBlock )
        delete[] pBlock;
    pBlock = 0;
}


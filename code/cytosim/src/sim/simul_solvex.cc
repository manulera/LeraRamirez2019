// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "simul.h"
#include "array.h"
#include "mecable.h"
#include "matsparsesym1.h"
#include "bicgstab.h"
#include "cblas.h"
extern Random RNG;


/// Solves the motion of Objects along the X axis
/**
 Used to solve the motion of Mecables in 1D, along the X axis.
 Each Mecable is represented by only one coordinate X.
 */
class Meca1D
{
    size_t allocated;
public:
   
    Array<Mecable *> objs;   ///< list of mobile objects

    real * vPTS;         ///< position of the points
    real * vFOR;         ///< base points of forces
    real * vMOB;         ///< the mobility of the objects
    real * vRHS;         ///< right-hand side term in the equation
    real * vTMP;         ///< intermediate of calculus

    MatrixSparseSymmetric1  mB;    ///< isotropic symmetric part of the dynamic

    Meca1D()
    {
        allocated = 0;
        vPTS = 0;
        vFOR = 0;
        vMOB = 0;
        vRHS = 0;
        vTMP = 0;
    }
    
    void clear()
    {
        objs.clear();
    }
    
    void add(Fiber * fib)
    {
        objs.push_back(fib);
    }

    void prepare(real time_step, real kT)
    {
        if ( size() > allocated )
        {
            // Keep memory aligned to 64 bytes:
            const size_t chunk = 64 / sizeof(real);
            // make a multiple of chunk to align memory:
            allocated = ( size() + chunk - 1 ) & ~( chunk -1 );
            
            if ( vFOR )  delete[] vFOR;
            if ( vPTS )  delete[] vPTS;
            if ( vTMP )  delete[] vTMP;
            if ( vMOB )  delete[] vMOB;
            if ( vRHS )  delete[] vRHS;
            
            vFOR = new real[ allocated ];
            vPTS = new real[ allocated ];
            vTMP = new real[ allocated ];
            vMOB = new real[ allocated ];
            vRHS = new real[ allocated ];
        }
        
        mB.resize(objs.size());
        mB.makeZero();

        blas_xzero(objs.size(), vFOR);
        blas_xzero(objs.size(), vRHS);

        unsigned ii = 0;
        for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci, ++ii )
        {
            Mecable * mec = *mci;
            mec -> matIndex(ii);
            mec -> setDragCoefficient();
            // Put the x coordinate of the origin in vPTS[ii]
            vPTS[ii] = mec->posPoint(0).XX;
            vMOB[ii] = time_step / mec->dragCoefficient();
        }
    }
    
    void addClamp(Matrix::index_type ii, real w, real dx)
    {
        mB(ii, ii) -= w;
        vFOR[ii]   += w * dx;
    }
    
    void addLink(Matrix::index_type ii, Matrix::index_type jj, real w, real dx)
    {
        mB(ii, ii) -= w;
        mB(ii, jj) += w;
        mB(jj, jj) -= w;
        vFOR[ii] += w * dx;
        vFOR[jj] -= w * dx;
    }
    
    real brownian(real kT)
    {
        real res = INFINITY;
        for ( unsigned ii = 0; ii < objs.size(); ++ii )
        {
            real b = sqrt( 2 * kT * vMOB[ii] );
            vRHS[ii] = vMOB[ii] * vFOR[ii] + b * RNG.gauss();
            if ( b < res )
                res = b;
        }
        return res;
    }
    
    bool solve(real precision)
    {
        mB.prepareForMultiply();
        static LinearSolvers::Allocator allocator;
        LinearSolvers::Monitor monitor(size(), precision);
        LinearSolvers::BCGS(*this, vRHS, vPTS, monitor, allocator);
        return monitor.converged();
    }
    
    void moveMecables()
    {
        unsigned ii = 0;
        for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci, ++ii )
        {
            Mecable * mec = *mci;
            // Move the Mecable along the X direction as calculated
            mec -> translate(Vector(vPTS[ii], 0, 0));
        }
    }
    
    //----------------------------- Implementation of LinearSolvers::LinearOperator
    
    unsigned size() const { return objs.size(); }
    
    /// Y <- X - mob * A * X
    void multiply(const real * X, real * Y) const
    {
        assert_true( X != Y  &&  X != vTMP  &&  Y != vTMP );
        
        blas_xzero(objs.size(), vTMP);
        mB.vecMulAdd(X, vTMP);
        
        for( unsigned ii = 0; ii < objs.size(); ++ii )
            Y[ii] = X[ii] - vMOB[ii] * vTMP[ii];
    }
};


//==========================================================================
//                                SOLVE
//==========================================================================


void Simul::solveX()
{
    static Meca1D meca1D;

    //-----initialize-----
    
    meca1D.clear();
    
    for(Fiber * fib = fibers.first(); fib; fib=fib->next())
        meca1D.add(fib);
    
    if ( meca1D.size() == 0 ) return;

    meca1D.prepare(prop->time_step, prop->kT);
    
    //-----set matrix-----

    for ( Couple * co = couples.firstAA(); co ; co=co->next() )
    {
        Hand const* h1 = co->hand1();
        Hand const* h2 = co->hand2();
        
        const Matrix::index_type  i1 = h1->fiber()->matIndex();
        const Matrix::index_type  i2 = h2->fiber()->matIndex();
        assert_true( i1 != i2 );
        
        meca1D.addLink(i1, i2, co->stiffness(), h2->pos().XX - h1->pos().XX);
    }
    
    for ( Single * gh = singles.firstA(); gh ; gh=gh->next() )
    {
        Hand const* h = gh->hand();
        const Matrix::index_type ii = h->fiber()->matIndex();
        
        meca1D.addClamp(ii, gh->prop->stiffness, gh->position().XX - h->pos().XX);
    }
    
    //-----resolution-----

    real estimate = meca1D.brownian(prop->kT);
    
    if ( meca1D.solve(prop->tolerance * estimate) )
    {
        meca1D.moveMecables();
    }
    else
    {
        std::cerr << "Simul::solveX() did not converge!\n";
    }
}   



// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 projections performed with explicit matrices
 This is a slow method, but it can be useful to compare with other methods
*/
//#include "vecprint.h"


void RigidFiber::buildProjection()
{
    //reset all variables for the projections:
    rfAllocated  = 0;
    mtP          = 0;
    mtDiffP      = 0;
    mtJJtiJ      = 0;
    mtJJtiJforce = 0;
}


void RigidFiber::allocateProjection(const unsigned int nbp)
{
    if ( rfAllocated < nbp )
    {
        //std::clog << reference() << "allocateProjection(" << nbp << ")\n";
        if ( mtP )          delete[] mtP;
        if ( mtDiffP )      delete[] mtDiffP;
        if ( mtJJtiJ )      delete[] mtJJtiJ;
        rfAllocated  = nbp;
        mtP          = new real[DIM*nbp*DIM*nbp];
        mtDiffP      = new real[DIM*nbp*DIM*nbp];
        mtJJtiJforce = new real[nbp];
        mtJJtiJ      = new real[DIM*nbp*nbp];
    }
}


void RigidFiber::destroyProjection()
{
    //std::clog << reference() << "destroyProjection\n";
    if ( mtP )          delete[] mtP;
    if ( mtDiffP )      delete[] mtDiffP;
    if ( mtJJtiJforce ) delete[] mtJJtiJforce;
    if ( mtJJtiJ )      delete[] mtJJtiJ;
    mtP          = 0;
    mtDiffP      = 0;
    mtJJtiJforce = 0;
    mtJJtiJ      = 0;
}


/*
 Computes the projection matrix
 @code
 P = I - J' ( J J' )^-1 J
 @endcode
 that is associated with the length constraints:
 @code
 | point(p+1) - point(p) |^2 = lambda^2
 @endcode
 */
void RigidFiber::makeProjection()
{
    const unsigned int nbc = nbSegments();             //number of constraints
    const unsigned int nbv = DIM * nbPoints();         //number of variables
    assert_true( nbc > 0 );
    
    assert_true( rfAllocated >= nbPoints() );
    
    //----- we allocate the arrays needed:
    real* J   = new real[nbv*nbc];
    real* JJt = new real[2*nbc];
    
    //------------compute the projection matrix
    real x;
    Vector v, w, dv, dw;
    int ofs, jj, kk;
    int info=0;
    
    blas_xzero(nbv*nbc, J);
    
    //set up the Jacobian matrix J and the diagonals of J * Jt
    w  = posP(0);
    for ( jj = 0; jj < nbc ; ++jj )
    {
        //set J:
        for ( ofs = 0; ofs < DIM ; ++ofs )
        {
            kk = DIM * jj + ofs;
            x = psPos[kk+DIM] - psPos[kk];
            J[jj+nbc*kk      ] = -x;
            J[jj+nbc*(kk+DIM)] =  x;
        }
        
        //set the diagonal and off-diagonal term of JJt:
        v  = w;
        w  = posP(jj+1);
        dv = dw;
        dw = w - v;
        JJt[jj] = 2 * dw.normSqr();     //diagonal term
        JJt[jj+nbc] = - dv * dw;        //off-diagonal term
    }
    
    // JJtiJ <- J
    blas_xcopy( nbc * nbv, J, 1, mtJJtiJ, 1 );
    
    // JJtiJ <- inv( JJt ) * J
    lapack_xptsv(nbc, nbv, JJt, JJt+nbc+1, mtJJtiJ, nbc, &info);
    if ( info ) ABORT_NOW("lapack_ptsv() failed");
    
    // mtP <-  - Jt * JJtiJ
    blas_xgemm('T', 'N', nbv, nbv, nbc, -1., J, nbc,
               mtJJtiJ, nbc, 0., mtP, nbv );
    
    // mtP <- mtP + I
    for ( jj = 0; jj < nbv*nbv; jj += nbv+1 )
        mtP[jj] += 1.0;
    
    //printf(" m%lu\n", name); VecPrint::matPrint( nbv, nbv, mtP );
    delete[] J;
    delete[] JJt;
}


void RigidFiber::projectForces( const real* X, real s, real* Y, real* ) const
{
    const unsigned nbv = DIM * nbPoints();
    blas_xsymv('U', nbv, s, mtP, nbv, X, 1, 0.0, Y, 1);
}


void RigidFiber::computeTensions(const real* force)
{
    const unsigned nbs = nbSegments();
    const unsigned nbv = DIM * nbPoints();
    
    // calculate the lagrangian coefficients:
    blas_xgemv('N', nbs, nbv, 1., mtJJtiJ, nbs, force, 1, 0., rfLag, 1);
}


//------------------------------------------------------------------------------


void RigidFiber::makeProjectionDiff(const real* force) const
{
    const unsigned nbs = nbSegments();             //number of constraints
    const unsigned nbv = DIM * nbPoints();         //number of variables
    assert_true( nbs > 0 );
    
    // calculate the lagrangian coefficients:
    blas_xgemv('N', nbs, nbv, 1., mtJJtiJ, nbs, force, 1, 0., rfLag, 1);
    
    //printf("Lagrange: "); sMath::vecPrint(std::clog, nbc, rfLag);
    
    // select expensive forces ( lagrangian > 0 )
    for ( int ii = 0; ii < nbs; ++ii )
    {
        if ( rfLag[ii] > 0 )
            mtJJtiJforce[ii] = rfLag[ii];
        else
            mtJJtiJforce[ii] = 0;
    }
    
    //printf("diffP ");VecPrint::vecPrint(std::clog, nbs, mtJJtiJforce);
    
    //set up the first term in the derivative of J with respect to variable x[ii]
    //set up term  P * (DJ)t (JJti) J force:
    for ( int jj = 0; jj < nbv; ++jj )
    {
        real* coljj = mtDiffP + nbv * jj;
        blas_xzero(nbv, coljj);
        int lin = jj / DIM;
        if ( lin > 0 ) {
            coljj[jj-DIM] = +mtJJtiJforce[lin-1];
            coljj[jj    ] = -mtJJtiJforce[lin-1];
        }
        if ( lin < nbs ) {
            coljj[jj    ] += -mtJJtiJforce[lin];
            coljj[jj+DIM]  = +mtJJtiJforce[lin];
        }
    }

    /*
     The final matrix is symmetric, for any force,
     as can be seen from the above relations to set its columns
     */
    //printf("projectionDiff\n");
    //VecPrint::matPrint(std::clog, nbv, nbv, mtDiffP);
}


void RigidFiber::addProjectionDiff( const real* X, real* Y ) const
{
    unsigned int nbv = DIM * nbPoints();
    blas_xsymv('U', nbv, 1.0, mtDiffP, nbv, X, 1, 1.0, Y, 1);
}


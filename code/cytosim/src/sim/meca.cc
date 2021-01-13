// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 * ------------------------------------------------------------------------------
 *                   -- Meca is the heart of cytosim --
 * ------------------------------------------------------------------------------
 *              solving the equations of motion for the PointSet
 *                   implicit integration, sparse matrix
 * ------------------------------------------------------------------------------
 * @todo See if Lagrangian dynamics could work better than constrainted dynamics
 * @todo Implement the PARDISO sparse matrix format
 * @todo Check if IDR(s) would perform better than BCGS
 * @todo Calculate the preconditionner in single precision
 * ------------------------------------------------------------------------------
 */

#include <fstream>

#include "meca.h"
#include "mecable.h"
#include "messages.h"
#include "cblas.h"
#include "clapack.h"
#include "bicgstab.h"
#include "exceptions.h"
#include "vecprint.h"

#include "meca_inter.cc"

/**
 Add correction term to the constrainted dynamics
 The effect is to stabilize fibers under traction, at some modest CPU cost.
*/
//#define ADD_PROJECTION_DIFF


/**
 The forces are usually:
 @code
 force = vBAS + ( mB + mC + mR + mdiffP ) * vPTS
 @endcode
 where mR represent the bending elasticity of fibers.
 
 The term RIGIDITY_IN_MATRIX Toggles between:
 1- Including mR into mB with Mecable::addRigidityUpper()
    The matrix type for mB should be able to cope with many near-diagonal terms

 2- Calculating mR on the fly with Mecable::addRigidity()
    This option is usually faster.
 
 RIGIDITY_IN_MATRIX should not be defined
 */
//#define RIGIDITY_IN_MATRIX


/// this define will enable explicit integration (should be off)
//#define EXPLICIT_INTEGRATION


/// this is used to check the validity of the results
#define not_a_number(x) ((x) != (x))


/// use this to generate code for validation
#define DEBUG_MECA 0


//------------------------------------------------------------------------------
#pragma mark - Allocate

Meca::Meca()
: objs(32)
{
    nbPts = 0;
    largestBlock = 0;
    allocated = 0;
    vPTS = 0;
    vSOL = 0;
    vBAS = 0;
    vRND = 0;
    vRHS = 0;
    vFOR = 0;
    vTMP = 0;
    use_mC = false;
    displayInteractions = false;
}


void allocate_vector(unsigned s, real *& ptr, bool reset)
{
    if ( ptr )
        delete [] ptr;
    ptr = new real[s];
    if ( reset )
        blas_xzero(s, ptr);
}


void Meca::allocate(unsigned alc)
{
    //allocate the vectors
    if ( alc > allocated )
    {
        // Keep memory aligned to 64 bytes:
        const size_t chunk = 64 / sizeof(real);
        // make a multiple of chunk to align pointers:
        allocated = ( alc + chunk - 1 ) & ~( chunk -1 );
        
        allocate_vector(DIM*allocated, vPTS, 1);
        allocate_vector(DIM*allocated, vSOL, 1);
        allocate_vector(DIM*allocated, vBAS, 0);
        allocate_vector(DIM*allocated, vRND, 1);
        allocate_vector(DIM*allocated, vRHS, 1);
        allocate_vector(DIM*allocated, vFOR, 1);
        allocate_vector(DIM*allocated, vTMP, 0);
    }
}

void Meca::release()
{
    //std::clog << "Meca::release()\n";
    if ( vPTS ) delete[] vPTS;
    if ( vSOL ) delete[] vSOL;
    if ( vBAS ) delete[] vBAS;
    if ( vRND ) delete[] vRND;
    if ( vRHS ) delete[] vRHS;
    if ( vFOR ) delete[] vFOR;
    if ( vTMP ) delete[] vTMP;
    vPTS = 0;
    vSOL = 0;
    vBAS = 0;
    vRND = 0;
    vRHS = 0;
    vFOR = 0;
    vTMP = 0;
}

void Meca::clear()
{
    objs.clear();
    nbPts = 0;
    largestBlock = 0;
    use_mC = false;
}

/**
 Register a new Mecable, sets its position in the matrix
 and update Meca::nbPts and Meca::largestBlock
 */
void Meca::add(Mecable* o)
{
    objs.push_back(o);
    
    o->matIndex(nbPts);
    
    const unsigned int n = o->nbPoints();
    
    nbPts += n;
    
    if ( largestBlock < n )
        largestBlock = n;
}


//------------------------------------------------------------------------------
#pragma mark - Multiply

/**
 Compute the linear part of the forces.
 The forces in a system with coordinates X are:
 @code
 forces = (mB+mC)*X + vBAS
 @endcode
 
 This will calculate the part that is variable with 'X', and add it to 'Y':
 @code
 Y = Y + (mB+mC)*X
 @endcode
 */
void Meca::addLinearForces(const real* X, real* Y) const
{    
    // PARALLELIZABLE CODE after modifications
    // multiplication by mB and mC can be done in parallel, summing afterward

    // Y <- Y + mB * X
#if ( DIM == 1 )
    mB.vecMulAdd( X, Y );
#elif ( DIM == 2 )
    mB.vecMulAddIso2D( X, Y );
#elif ( DIM == 3 )
    mB.vecMulAddIso3D( X, Y );
#endif
    
    // Y <- Y + mC * X
    if ( use_mC )
        mC.vecMulAdd( X, Y );
}


void Meca::addRigidity(const real* X, real* Y) const
{
#if ( DIM > 1 ) && !defined(RIGIDITY_IN_MATRIX)
    // PARALLELIZABLE LOOP
    for ( Mecable ** mec = objs.begin(); mec < objs.end(); ++mec )
    {
        const index_type inx = DIM * (*mec)->matIndex();
        (*mec)->addRigidity( X+inx, Y+inx );
    }
#endif
}


/**
 calculate the matrix product needed for the conjugate gradient algorithm
 @code
 Y = X - time_step * P ( mB + mC + P' ) * X;
 @endcode
*/
void Meca::multiply( const real* X, real* Y ) const
{
    // vTMP is used for temporary storage !
    assert_true( X != Y  &&  X != vTMP  &&  Y != vTMP );

    // vTMP <= Forces = ( mB + mC ) * X
    blas_xzero(DIM*nbPts, vTMP);
    addLinearForces(X, vTMP);
    
#if ( 0 )
    // Bypass the Projection:
    // Y <- Y - time_step * X
    addRigidity(X, vTMP);
    blas_xcopy(DIM*nbPts, X, 1, Y, 1);
    blas_xaxpy(DIM*nbPts, -time_step, vTMP, 1, Y, 1);
    return;
#endif
    
    // PARALLELIZABLE LOOP
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable const * mec = *mci;
        const index_type inx = DIM * mec->matIndex();
#if ( DIM > 1 ) && !defined(RIGIDITY_IN_MATRIX)
        mec->addRigidity(X+inx, vTMP+inx);
#endif
#ifdef ADD_PROJECTION_DIFF
        mec->addProjectionDiff(X+inx, vTMP+inx);
#endif
        mec->setSpeedsFromForces(vTMP+inx, -time_step, Y+inx);

        // Y <- Y + X
        blas_xaxpy(DIM*mec->nbPoints(), 1.0, X+inx, 1, Y+inx, 1);
    }
}

//------------------------------------------------------------------------------
#pragma mark - Helper functions


/**
 Fill-in matrix 'dst' as the duplicate of 'src', for each 'DIM' dimension.
 For 'DIM==1', this makes a simple copy
 For 'DIM==2', 'src' is copied twice, into odd indices, and into even indices.
 For 'DIM==3', three copies of 'src' are made into 'dst'.
 
 Both 'src' and 'dst' must be symmetrix square matrices. 
 The size of 'dst' is DIM times the size of 'src'.
 Only the upper diagonal of 'src' is specified.
 Matrix Y is specified in full.
 */

void duplicate_matrix(unsigned siz, real const* src, real * dst)
{
    unsigned ddd = DIM * siz;
    
    blas_xzero(ddd*ddd, dst);
    
    for ( unsigned ii = 0; ii < siz; ++ii )
    {
        real xx = src[ii + siz * ii];
        
        unsigned kk = ( ddd+1 ) * DIM * ii;
        for ( unsigned d = 0; d < DIM; ++d, kk += ddd+1 )
            dst[kk] = xx;
        
        for ( unsigned jj = ii+1; jj < siz; ++jj )
        {
            xx = src[ii + siz * jj];
            kk = DIM * ( ii + ddd * jj );
            unsigned ll = DIM * ( jj + ddd * ii );
            for ( unsigned d = 0; d < DIM; ++d, kk += ddd+1, ll += ddd+1 )
            {
                dst[kk] = xx;
                dst[ll] = xx;
            }
        }
    }
    
#if ( 0 )
    std::clog << "\nOriginal:\n";
    VecPrint::matPrint(std::clog, siz, siz, src);
    std::clog << "Duplicated:\n";
    VecPrint::matPrint(std::clog, ddd, ddd, dst);
#endif
}


/**
 This should symmetrize a matrix, and also copy the terms
 that are within the first subspace `X` into the other dimensions
 */
void expand_matrix(unsigned siz, real * mat)
{
#if ( 0 )
    std::clog << "\nOriginal:\n";
    VecPrint::matPrint(std::clog, siz, siz, mat, siz);
#endif
    
    for ( unsigned jj = 0; jj < siz; jj += DIM  )
    {
        for ( unsigned ii = 0; ii < jj; ii += DIM  )
        {
            real val = mat[ii+siz*jj];
            // expand term in other dimensions:
#if ( DIM >= 2 )
            mat[ii+1+siz*(jj+1)] = val;
#endif
#if ( DIM >= 3 )
            mat[ii+2+siz*(jj+2)] = val;
#endif
            
            // symmetrize matrix:
            mat[jj  +siz* ii   ] = val;
#if ( DIM >= 2 )
            mat[jj+1+siz*(ii+1)] = val;
#endif
#if ( DIM >= 3 )
            mat[jj+2+siz*(ii+2)] = val;
#endif
        }
        // expand diagonal term in other dimensions:
        real val = mat[jj+siz*jj];
#if ( DIM >= 2 )
        mat[jj+1+siz*(jj+1)] = val;
#endif
#if ( DIM >= 3 )
        mat[jj+2+siz*(jj+2)] = val;
#endif
    }

#if ( 0 )
    std::clog << "Expanded:\n";
    VecPrint::matPrint(std::clog, siz, siz, mat, siz);
#endif
}



/**
 Set to zero all the terms that are not within 'diag' from the diagonal.
 With 'diag==0' the entire matrix is set to zero.
 With 'diag==1', only the diagonal is kept.
 With 'diag==2', the matrix is made tri-diagonal
 etc.
 */
void truncate_matrix(unsigned siz, real* mat, unsigned diag)
{
#if ( 0 )
    std::clog << "\nOriginal:\n";
    VecPrint::matPrint(std::clog, siz, siz, mat);
#endif
    
    for ( unsigned ii = 0; ii < siz; ++ii )
    {
        real * col = mat + siz * ii;
        for ( unsigned jj = 0; jj+diag < ii; ++jj )
            col[jj] = 0.0;
        
        for ( unsigned kk = ii+diag+1; kk < siz; ++kk )
            col[kk] = 0.0;
    }
    
#if ( 0 )
    std::clog << "Truncated:\n";
    VecPrint::matPrint(std::clog, siz, siz, mat);
#endif
}


/// sim(element^2) / sum(diagonal^2)
real off_diagonal_norm(int siz, real * mat)
{
    real all = 0;
    for ( int k = 0; k < siz*siz; ++k )
        all += mat[k] * mat[k];

    real dia = 0;
    for ( int k = 0; k < siz*siz; k+=siz+1 )
        dia += mat[k] * mat[k];

    return sqrt( ( all - dia ) / dia );
}


/// set all values between '-val' and 'val' to zero
void threshold_matrix(int siz, real * mat, real val)
{
    for ( int k = 0; k < siz*siz; ++k )
    {
        if ( fabs(mat[k]) < val )
            mat[k] = 0.0;
    }
}


/// set 'mat' of order `siz` to `diag * I`
void diagonal_matrix(int siz, real * mat, real diag)
{
    for ( int k = 0; k < siz*siz; ++k )
        mat[k] = 0.0;
    for ( int k = 0; k < siz*siz; k+=siz+1 )
        mat[k] = diag;
}


/// erase all off-diagonal terms in `mat` of order `siz`
void make_diagonal(int siz, real * mat)
{
    for ( int j = 0; j < siz; ++j )
    {
        real * col = mat + j * siz;
        int i;
        for ( i = 0; i < j; ++i )
            col[i] = 0.0;
        for ( i = j+1; i < siz; ++i )
            col[i] = 0.0;
    }
}


/**
 Convert a full matrix into a LAPACK banded matrix data with 'diag' diagonals.
 `src` is a square matrix of side 'siz'
 
 `dst` is of size ldd * siz, with ldd > ku+kl+1
 */
void banded_matrix(int siz, real const* src, int ku, int kl, real * dst, int ldd)
{
    assert_true( ldd > ku + kl + 1 );
#if ( 0 )
    std::clog << "\noriginal:\n";
    VecPrint::matPrint(std::clog, siz, siz, src);
#endif
    
    for ( int jj = 0; jj < siz; ++jj )
    {
        real * col = dst + ldd * jj;
        // number of value not set before the diagonal:
        int off = ldd - kl - 1 - ( jj < ku ? jj : ku );
        // index offset between 'src' and 'col':
        int shift = ( siz + 1 ) * jj - ( ldd - kl - 1 );
        // last index in 'col' to be set:
        int e = ldd - ( jj + kl + 1 > siz ? jj + kl + 1 - siz : 0 );
        //std::clog << " ldd " << ldd << " off " << off << " shift " << shift << " e " << e << "\n";
        
        int ii;
        for ( ii = 0; ii < off; ++ii )
            col[ii] = 0;
        for ( ; ii < e; ++ii )
            col[ii] = src[ii+shift];
        for ( ; ii < ldd; ++ii )
            col[ii] = 0;
    }
    
#if ( 0 )
    std::clog << "banded:\n";
    VecPrint::matPrint(std::clog, ldd, siz, dst);
#endif
}


/*
 uses power iterations to estimate the largest eigenvalue of `mat * tam - I`
 @returns an estimate of the largest eigenvalue
 The precision of the estimate is low: 10%
 */
real largest_eigenvalue(int siz, real * mat, real * tam, real * vec, real * tmp)
{
    assert_true(siz > 0);
    real oge, eig = blas_xnrm2(siz,vec,1);
    
    int n;
    for ( n = 0; n < siz; n += 2 )
    {
        blas_xgemv('N', siz, siz, 1.0/eig, mat, siz, vec, 1,  0.0,     tmp, 1);
        blas_xgemv('N', siz, siz, 1.0,     tam, siz, tmp, 1, -1.0/eig, vec, 1);
        oge = blas_xnrm2(siz, vec, 1);
        //VecPrint::vecPrint(std::clog, std::min(16, siz), vec, 3);
        //std::clog << "        iter " << std::setw(3) << n << " eigen " << std::setw(9) << oge << "\n";
        blas_xgemv('N', siz, siz, 1.0/oge, mat, siz, vec, 1,  0.0,     tmp, 1);
        blas_xgemv('N', siz, siz, 1.0,     tam, siz, tmp, 1, -1.0/oge, vec, 1);
        eig = blas_xnrm2(siz, vec, 1);
        if ( fabs( oge - eig ) < 0.1 * fabs( eig ) )
            break;
        //VecPrint::vecPrint(std::clog, std::min(16, siz), vec, 3);
        //std::clog << "        iter " << std::setw(3) << n+1 << " eigen " << std::setw(9) << eig << "   " << std::setw(9) << oge << "\n";
    }
    //std::clog << "  power iter " << std::setw(3) << n << " eigen " << std::setw(9) << eig << "   " << std::setw(9) << oge << "\n";
    
    return std::max(eig, oge);
}


/// return LAPACK's optimal work size, for a matrix of size 'ord'
unsigned lapack_work_size(const int ord)
{
    real tmpA, tmpW;
    int tmpi, info = 0;
    lapack_xgetri(ord, &tmpA, ord, &tmpi, &tmpW, -1, &info);
    if ( info == 0 )
    {
        //std::clog << "Lapack::dgetri optimal size is " << work_size << std::endl;
        return (unsigned)tmpW;
    }
    return 1024;
}

//------------------------------------------------------------------------------
#pragma mark - Precondition


/**
 Return the diagonal block of the full matrix corresponding to an Object,
 which is:
 @code
 I - time_step * P ( B + C + P' ).
 @endcode
 This block is square but not symmetric!
*/
void Meca::extractBlock(real* res, const Mecable * mec, real* tmp, real* vec) const
{
    const unsigned ps = mec->nbPoints();
    const unsigned bs = DIM * ps;
    
    blas_xzero(bs*bs, tmp);
    
#if ( DIM > 1 ) && ! defined(RIGIDITY_IN_MATRIX)
    // set the Rigidity terms:
    mec->addRigidityUpper(tmp);
    //std::clog<<"Rigidity block " << mec->reference() << "\n";
    //VecPrint::matPrint(std::clog, bs, bs, tmp, bs);
#endif
    
    mB.addTriangularBlock(tmp, bs, mec->matIndex(), ps, DIM);
    
    expand_matrix(bs, tmp);
    
    if ( use_mC )
        mC.addDiagonalBlock(tmp, bs, DIM*mec->matIndex(), bs);
    
#if ( 0 )
    std::clog<<"mB+mC block:\n";
    VecPrint::matPrint(std::clog, bs, bs, tmp, bs);
#endif
    
#ifdef ADD_PROJECTION_DIFF
    // Include the corrections P' in preconditioner, vector by vector.
    blas_xzero(bs, vec);
    for ( unsigned ii = 0; ii < bs; ++ii )
    {
        vec[ii] = 1.0;
        mec->addProjectionDiff(vec, tmp+bs*ii);
        vec[ii] = 0.0;
    }
#if ( 0 )
    std::clog<<"dynamic with P'\n";
    VecPrint::matPrint(std::clog, bs, bs, tmp, bs);
#endif
#endif
    
#if ( 0 )
    // Bypass the projection (never enable this!)
    for ( unsigned ii = 0; ii < bs; ++ii )
    {
        for ( int u = 0; u < bs; ++u )
            res[bs*ii+u] = -time_step * tmp[bs*ii+u];
        res[bs*ii+ii] += 1.0;
    }
    return;
#endif

    //compute the projection, by applying it to each column vector:
    /*
     This is slow and could be parallelized by having setSpeedFromForces()
     accept multiple vectors as arguments
     */
    for ( unsigned ii = 0; ii < bs; ++ii )
    {
        mec->setSpeedsFromForces(tmp+bs*ii, -time_step, res+bs*ii);
        // add Identity matrix:
        res[bs*ii+ii] += 1.0;
    }
}



/**
 This version builds the diagonal block directly from Meca:multiply().
 This is a very slow method that calls 'multiply()' n-times, where
 'n' is the size of the block. It does not use any of the Meca vectors.
 
 It should used for validation only.
*/
void Meca::extractBlockS(real* res, const Mecable* mec) const
{
    const unsigned sz = DIM * nbPoints();
    const unsigned bs = DIM * mec->nbPoints();
    const unsigned off= DIM * mec->matIndex();
    
    assert_true( off+bs <= sz );
    real * vec = new real[sz];
    real * tmp = new real[sz];
    
    blas_xzero(sz, vec);
    //blas_xzero(bs*bs, res);
    
    for ( unsigned jj=0; jj<bs; ++jj )
    {
        vec[jj+off] = 1.0;
        multiply(vec, tmp);
        vec[jj+off] = 0.0;
        blas_xcopy(bs, tmp+off, 1, res+jj*bs, 1);
    }
    
    delete [] vec;
    delete [] tmp;
}
 

/**
 @todo We could call lapack_xgetrf() here, and lapack_xgetrs() in Meca::precondition()
 Arrays 'tmp' and 'wrk' should be of size (nb_points*DIM)^2 or more
 'mode' should be {1, 2 or 3}
 */
int Meca::computePreconditionner(Mecable* mec, int mode, int* ipiv, real* tmp, real* wrk, int wrksize)
{
    assert_true( ipiv && tmp && wrk );
    
    unsigned bs = DIM * mec->nbPoints();
    bool may_keep = false;
    
    if ( mec->blockSize() == bs )
        may_keep = ( mode == 2 );
    else
        mec->allocateBlock(bs);
    
    real* blk = mec->block();
    real* vec = vTMP + DIM * mec->matIndex();

    if ( may_keep )
    {
        // extract diagonal matrix block corresponding to this Mecable:
        extractBlock(wrk, mec, tmp, vec);
        
        // chose initial vector for power iteration
        blas_xcopy(bs, vRHS+DIM*mec->matIndex(), 1, vec, 1);
        
        real eig = largest_eigenvalue(bs, wrk, blk, vec, tmp);

#if ( 0 )
        // matrix-matrix multiplication requires O(N^3) operations:
        blas_xgemm('N','N',bs,bs,bs,1.0,blk,bs,wrk,bs,0.0,tmp,bs);
        
        if ( 0 )
        {
            unsigned s = std::min(10U, bs);
            std::clog<<"*block " << mec->reference() << " size " << bs << "\n";
            VecPrint::matPrint(std::clog, s, s, wrk, bs, 3);
        }

        //printf("\n precond * block:\n"); VecPrint::matPrint(std::clog,bs,bs,tmp, bs);
        
        for ( int k=0; k < bs*bs; k += 1+bs )
            tmp[k] -= 1.0;
        
        real err = blas_xnrm2(bs*bs,tmp,1) / bs;
        std::clog << std::setw(4) << bs << "  NORM " << std::setw(9) << err << "   EIGEN " << std::setw(9) << eig << "\n";
#endif
 
        //std::clog << "mec " << std::setw(4) << mec->identity();
        if ( eig < 1.0 )
        {
            //std::clog << "    keep     || 1 - precond * block || = " << eig << "\n";
            return 0;
        }
        else
        {
            //std::clog << "    RENEW    || 1 - precond * block || = " << eig << "\n";
            blas_xcopy(bs*bs, wrk, 1, blk, 1);
        }
    }
    else
    {
        // extract diagonal matrix block corresponding to this Mecable:
        extractBlock(blk, mec, tmp, vec);
    }
    
#if ( 0 )
    // DEBUG: compare with block extracted using different method:
    extractBlockS(wrk, mec);
    real err = max_diff(bs*bs, blk, wrk);
    std::clog<<" norm( block - block' ) = " << err << "\n";

    if ( err > REAL_EPSILON )
    {
        unsigned s = std::min(16U, bs);
        std::clog<<"*block " << mec->reference() << " size " << bs << "\n";
        VecPrint::matPrint(std::clog, s, s, blk, bs, 3);
        std::clog<<" block'\n";
        VecPrint::matPrint(std::clog, s, s, wrk, bs, 3);
        //mB.printSparse(std::clog);
        exit(1);
    }
#endif
    

    int info = 0;

    if ( mode <= 2 )
    {
        /**
         We invert here the full block using LAPACK general matrix functions
         the workload scales as SIZE ^ 3, as a function of the size of the block
         */
        //TEST truncate_matrix(bs, blk, 3*DIM-1);
        
        // LU factorization:
        lapack_xgetf2(bs, bs, blk, bs, ipiv, &info);
        
        if ( info ) return 2;      //failed to factorize matrix !!!
        
        // calculate matrix inverse from factorization:
        lapack_xgetri(bs, blk, bs, ipiv, wrk, wrksize, &info);
        
        //TEST threshold_matrix(bs, blk, 0);
    }
    else if ( mode == 3 )
    {
        /**
         We compute here a preconditionner containing 2*DIM+1 diagonals.
         This should contain most of the terms due to bending elasticity of fibers,
         and the preconditionner is quite good therefore.
         This preconditionner is easier to compute, since inverting a
         banded-matrix scales as K * SIZE ^ 2, for a matrix of SIZE with K digonals
         
         F. Nedelec, 18.02.2017
         */
        const int diag = 3 * DIM - 1;  // bs - 1
        const int ldd = 3 * diag + 1;
        
        real * banded = new real[ldd*bs];
        
        banded_matrix(bs, blk, diag, diag, banded, ldd);
        
        // LU factorization:
        lapack_xgbtf2(bs, bs, diag, diag, banded, ldd, ipiv, &info);
        
        if ( info ) return 2;      //failed to factorize matrix !!!
        
        diagonal_matrix(bs, blk, 1.0);
        
        lapack_xgbtrs('N', bs, diag, diag, bs, banded, ldd, ipiv, blk, bs, &info);
        
        delete[] banded;
    }
    
    if ( info ) return 3;      //failed to invert matrix !!!
    
#if ( 0 )
    // DEBUG: print preconditionner block for visual inspection:
    //if ( mec->identity() == 30 )
    {
        unsigned s = std::min(16U, bs);
        std::clog<<"*inverse " << mec->reference() << "\n";
        VecPrint::matPrint(std::clog, s, s, blk, bs, 3);
    }
#endif
    
    if ( 0 ) // mec->identity() == 1 )
    {
        /*
         Multiply here the preconditionner block with the dynamic block to
         check if we recover the identity matrix
         */
        extractBlock(wrk, mec, tmp, vec);
        
        // printf("\nblock2:\n"); VecPrint::matPrint(std::clog,bs,bs,blk2);
        blas_xgemm('N','N',bs,bs,bs,1.0,blk,bs,wrk,bs,0.0,tmp,bs);
        
        //printf("\n precond * block:\n"); VecPrint::matPrint(std::clog,bs,bs,res);
        
        for ( int k=0; k < bs*bs; k += 1+bs )
            tmp[k] -= 1.0;
        
        real err = blas_xnrm2(bs*bs,tmp,1) / bs;
        
        std::clog << "mec" << std::setw(10) << mec->identity();
        std::clog << "    mode " << mode << " NORM2( 1 - precond * block ) = " << err << "\n";
    }
    
    return 0;
}


// Compute the preconditionner block by block
/**
 The code can be parallelized here:
 - allocate temporary memory for each thread
 - distribute block calculation to different threads
 */
int Meca::computePreconditionner(int mode)
{
    const unsigned vecsize = DIM * largestBlock;
    const unsigned matsize = vecsize * vecsize;
    
    unsigned wrksize = lapack_work_size(vecsize);
    if ( wrksize < matsize )
        wrksize = matsize;
    
    // allocate memory:
    int* ipiv = new int[vecsize];
    real* mat = new real[matsize];
    real* wrk = new real[wrksize];
    
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable * mec = *mci;
        assert_true( mec->nbPoints() <= largestBlock );
        int res = computePreconditionner(mec, mode, ipiv, mat, wrk, wrksize);
        if ( res )
            std::clog << "Warning: could not compute preconditionner for " << mec->reference() << "\n";
        mec->useBlock(res==0);
    }

    delete[] ipiv;
    delete[] wrk;
    delete[] mat;
    return 0;
}


void Meca::precondition(const real* X, real* Y) const
{
    // PARALLELIZABLE LOOP
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable const* mec = *mci;
        const unsigned bs = DIM * mec->nbPoints();
        const index_type inx = DIM * mec->matIndex();
        if ( mec->useBlock() )
        {
            // we multiply by the diagonal block that was calculated
            blas_xgemv('N', bs, bs, 1.0, mec->block(), bs, X+inx, 1, 0.0, Y+inx, 1);
        }
        else
        {
            // we effectively 'multiply' by the Identity Matrix
            blas_xcopy(bs, X+inx, 1, Y+inx, 1);
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Solve


/**
 Allocate and reset matrices and vectors necessary for Meca::solve(),
 copy coordinates of Mecables into vPTS[]
 */
void Meca::prepare(SimulProp const* prop)
{
#ifndef NDEBUG
    // verify all the Mecable indices and matrix size:
    unsigned int n = 0;
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        assert_true( (*mci)->matIndex() == n );
        n += (*mci)->nbPoints();
    }
    assert_true( n == nbPts );
#endif

    allocate(nbPts);

    //allocate the sparse matrices:
    mB.resize(nbPts);
    mC.resize(DIM*nbPts);
    
    //reset matrices:
    mB.makeZero();
    mC.makeZero();
    
    // reset base:
    blas_xzero(DIM*nbPts, vBAS);
    
    // get global time step
    time_step = prop->time_step;
    
    // PARALLELIZABLE LOOP
    // import coordinates of mecables:
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable * mec = *mci;
        mec->putPoints(vPTS+DIM*mec->matIndex());
        mec->prepareMecable();
#if ( DIM > 1 ) && defined(RIGIDITY_IN_MATRIX)
        //include the rigidity terms in matrix mB
        mec->addRigidityMatrix(mB, mec->matIndex(), 1);
#endif
    }
}


/**
 Prepare matrices mB and mC for multiplication
 This should be called after setInteractions()
 */
void Meca::prepareMatrices()
{
    mB.prepareForMultiply();
    
    if ( mC.nonZero() )
    {
        use_mC = true;
        mC.prepareForMultiply();
    }
    else
        use_mC = false;
}


/**
 Calculates the force in the objects, that can be accessed by Mecable::netForce()
 and calculate the speed of the objects in vRHS, in the abscence of Thermal motion,
 ie. the motion is purely due to external forces.

 This also sets the Lagrange multipliers for the Fiber.
 
 The function will not change the position of the Mecables.
 */
void Meca::computeForces()
{
    prepareMatrices();
    
    // calculate forces in vFOR, but without adding Brownian noise:
    blas_xcopy(DIM*nbPts, vBAS, 1, vFOR, 1);
    addLinearForces(vPTS, vFOR);
    
    // return forces to Mecable, and calculate resulting speed:
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable * mec = *mci;
        real * f = vFOR + DIM * mec->matIndex();
        mec->getForces(f);
        mec->computeTensions(f);
    }
}



/**
 This solves the equation:
 @code
 (Xnew - Xold)/time_step = P*force(X) + BrownianNoise
 @endcode
 
 Explicit integration is:
 @code
 Xnew = Xold + time_step * P * force + BrownianNoise
 @endcode
 
 Implicit integration using the linearized force(X) = A.X + B:
 @code
 ( I - time_step*P*A ) ( Xnew - Xold ) = time_step * P * force + BrownianNoise
 @endcode
 
 where, in all cases,
 @code
 force = A * Xold + B
 BrownianNoise = sqrt(2*kT*time_step*mobility) * Gaussian(0,1)
 @endcode
 
 Implicit integration enables a bigger time_step to be used.
 
 Vector 'vRHS' is the right-hand-side of the system, and 'vSOL' the solution:
 
 @code
 vRHS = time_step * P * force + BrownianNoise
 vSOL = Xnew - Xold
 @endcode
 */
unsigned Meca::solve(SimulProp const* prop, const int precondition)
{
    assert_true( time_step == prop->time_step );

    if ( nbPts == 0 )
        return 0;
    
    prepareMatrices();
    
    // calculate forces before constraints in vFOR:
    blas_xcopy(DIM*nbPts, vBAS, 1, vFOR, 1);
    addLinearForces(vPTS, vFOR);
    addRigidity(vPTS, vFOR);

    //----------------------------------------------------------------
 
    /* 
     Fill `vRND` with Gaussian random numbers 
     This operation can be done in parallel, in a separate thread
     */
    RNG.gauss_set(vRND, DIM*nbPts);
    //blas_xzero(DIM*nbPts, vRND);
    
    /*
     Add Brownian motions to 'vFOR', and calculate vRHS by multiplying by mobilities.
     As Brownian terms are added, we record the magnitude of the typical smallest
     scalar contribution in `noiseLevel`. The dynamics will later be solved with 
     a residual that is proportional to this level:
     SimulProp::tolerance * noiseLevel
     As long as SimulProp::tolerance is smaller than 1, this should allow for a
     level of numerical error is small with respect to the Brownian noise in
     the system, and the results should be physically appropriate.
     */
    
    real noiseLevel = INFINITY;
    
    //add the Brownian contribution
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable const * mec = *mci;
        const index_type inx = DIM * mec->matIndex();

        real n = mec->addBrownianForces(vFOR+inx, vRND+inx, prop->kT/time_step);
        
        if ( n < noiseLevel )
            noiseLevel = n;
        
        // calculate the right-hand-side of the system in vRHS:
        mec->setSpeedsFromForces(vFOR+inx, time_step, vRHS+inx, true);
        
#ifdef ADD_PROJECTION_DIFF
        mec->makeProjectionDiff(vFOR+DIM*mec->matIndex());
#endif
    }
    noiseLevel *= time_step;
    
#ifdef NEW_CYTOPLASMIC_FLOW
    /**
     Includes a constant fluid flow displacing all the objects along
     */
    if ( prop->flow.norm() > REAL_EPSILON )
    {
        PRINT_ONCE("NEW_CYTOPLASMIC_FLOW code enabled\n");
        Vector flow_dt = prop->flow * time_step;
        
        const real *const end = vRHS + DIM*nbPts;
        
        for ( real * mx = vRHS; mx < end; mx += DIM )
            flow_dt.add_to(mx);
    }
#endif
    
#ifdef EXPLICIT_INTEGRATION
    /*
     This implements the forward Euler integration, for testing purposes.
     The result is very inefficient, since we have built the entire stiffness matrix,
     which is not necessary for this simple explicit scheme.
     */
    blas_xaxpy(DIM*nbPts, 1.0, vRHS, 1, vPTS, 1);
    
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable * mec = *mci;
        mec->getPoints(vPTS+DIM*mec->matIndex());
        mec->getForces(vFOR+DIM*mec->matIndex());
    }
    return 1;
#endif
    
    /*
     Choose the initial guess for the solution of the system (Xnew - Xold):
     we could use the solution at the previous step, or a vector of zeros.
     Using the previous solution could be advantageous if the speed were 
     somehow continuous. However, the system is without inertia. In addition,
     objects are considered in a random order to build the linear system, such
     that the blocks from two consecutive iterations do not match.
     From this, using zero for the initial guess seems safer:
     */
    blas_xzero(DIM*nbPts, vSOL);

    /*
     We now solve the system MAT * vSOL = vRHS  by an iterative method:
     the tolerance is in scaled to the contribution of Brownian
     motions contained in vRHS, assuming that the error is equally spread 
     along all degrees of freedom, this should work for tolerance << 1
     here a printf() can be used to check that the estimate is correct:
    */ 
     //printf("noiseLeveld = %8.2e   variance(vRHS) / estimate = %8.4f\n", 
     //       noiseLevel, blas_xnrm2(DIM*nbPts, vRHS, 1) / (noiseLevel * sqrt(DIM*nbPts)) );
    
    // TEST: the tolerance to solve the system should be such that the solution
    // found does not depend on the initial guess.
    
    real residual;
    
    if ( noiseLevel > 0 )
        residual = noiseLevel * prop->tolerance;
    else
        residual = prop->tolerance;

    /*
     One can use static memory by making the Allocator static
     */
    static LinearSolvers::Allocator allocator;

    /*
     With exact arithmetic, biConjugate Gradient should converge at most
     in a number of iterations equal to the size of the linear system.
     This is the max limit that is set here to the number of iterations:
     */
    LinearSolvers::Monitor monitor(DIM*nbPts, residual);

    //------- call the iterative solver:
    //std::clog << "Solve: " << DIM*nbPts << "  " << residual_ask << std::endl;

    if ( precondition  &&  0 == computePreconditionner(precondition) )
        LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
    else
        LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator);
    
    //------- in case the solver did not converge, we try other methods:
    
    if ( !monitor.converged() )
    {
        MSG("Solver failed: precond %i flag %i, nb_iter %3i residual %.2e\n", 
            precondition, monitor.flag(), monitor.iterations(), monitor.residual());
        
        //---we try vRHS as a different initial seed:
        blas_xcopy(DIM*nbPts, vRHS, 1, vSOL, 1);
        
        //---reset tolerance and iteration counters:
        monitor.reset();
        
        //---try the same method again:
        if ( precondition )
            LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
        else
            LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator);
        
        //---check again for convergence:
        if ( monitor.converged() )
            MSG("Solver rescued by changing seed: nb_iter %3i residual %.2e\n", monitor.iterations(), monitor.residual());
        else 
        {
            //---use zero as an initial guess:
            blas_xzero(DIM*nbPts, vSOL);
            
            //---reset tolerance and iteration counters:
            monitor.reset();
            
            if ( precondition )
            {
                // try with the other preconditionner:
                if ( 0 == computePreconditionner(3-precondition) )
                    LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
                else // try without preconditionning:
                    LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator);
            }
            else
            {
                // try with a preconditioner
                if ( 0 == computePreconditionner(1) )
                    LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
                else
                    MSG("Failed to compute precondionner");
            }
            
            //---check again for convergence:
            if ( monitor.converged() )
                MSG("Solver rescued by changing precond: nb_iter %3i residual %.2e\n", monitor.iterations(), monitor.residual());
            else {
                //no solver could converge... this is really bad!
                //we could still try to change the initial guess, to recover convergence
                MSG("Solver dead nb_iter %i residual %.2e\n", monitor.iterations(), monitor.residual());
                throw Exception("convergence failure in solver");
                return monitor.iterations();
            }
        }
    }
    
#ifndef NDEBUG
    
    //check validity of the data:
    for( unsigned ii = 0; ii < DIM*nbPts; ++ii )
    {
        if ( not_a_number(vSOL[ii]) )
        {
            fprintf(stderr, "Meca::solve detected NaN\n");
            abort();
        }
    }
    
#endif

    //add the solution of the system (=dPTS) to the points coordinates
    blas_xaxpy(DIM*nbPts, 1., vSOL, 1, vPTS, 1);

    /*
     Re-calculate forces once the objects have moved:
     Bending elasticity of fibers, an internal forces, is not included,
     and Brownian terms in vFOR are also excluded.
     In this way the result returned to the fibers does not sum-up to zero,
     and is appropriate for example to calculate the effect of force on assembly.
     */
    blas_xcopy(DIM*nbPts, vBAS, 1, vFOR, 1);
    addLinearForces(vPTS, vFOR);


    // export new coordinates to Mecables:
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable * mec = *mci;
        mec->getForces(vFOR+DIM*mec->matIndex());
        mec->getPoints(vPTS+DIM*mec->matIndex());
    }
    
    // report on the matrix type and size, sparsity, and the number of iterations
    if ( prop->verbose )
    {
        std::stringstream oss;
        oss << "Meca " << DIM << "x" << nbPts;
        oss << " block " << largestBlock;
        oss << " " << mB.what();
        if ( use_mC ) oss << " " << mC.what();
        oss << " prec " << precondition;
        oss << " iter " << monitor.iterations();
        oss << " resid " << monitor.residual() << "\n";
        if ( prop->verbose > 2 )
            std::clog << oss.str();
        else
            MSG.put_line(oss.str(), 0);
    }
    return monitor.iterations();
}






//------------------------------------------------------------------------------
#pragma mark - Dump

/**
 Extract the full matrix associated with matVect, in `mat` of size `sz`.
 */
void Meca::getSystem(int sz, real * mat) const
{
    if ( sz != size() )
        throw InvalidIO("unexpected matrix dimension");
    
    real * src = new real[sz];
    real * res = new real[sz];
    
    blas_xzero(sz, src);
    blas_xzero(sz, res);
    
    for ( unsigned ii = 0; ii < sz; ++ii )
    {
        src[ii] = 1.0;
        multiply(src, res);
        blas_xcopy(sz, res, 1, mat+ii*sz, 1);
        src[ii] = 0.0;
    }
    
    delete [] res;
    delete [] src;
}


/**
 Save the full matrix associated with multiply(), in binary format
 */
void Meca::dumpSystem(FILE * file) const
{
    const unsigned sz = size();
    
    real * src = new real[sz];
    real * res = new real[sz];
    
    blas_xzero(sz, src);
    
    for ( unsigned ii = 0; ii < sz; ++ii )
    {
        src[ii] = 1.0;
        multiply(src, res);
        fwrite(res, sizeof(real), sz, file);
        src[ii] = 0.0;
    }
    
    delete [] res;
    delete [] src;
}


/**
 Save the elasticity matrix
 */
void Meca::dumpElasticity(FILE * file) const
{
    const unsigned sz = size();
    
    real * src = new real[sz];
    real * res = new real[sz];
    
    blas_xzero(sz, src);
    
    for ( unsigned ii = 0; ii < sz; ++ii )
    {
        src[ii] = 1.0;
        
        blas_xzero(sz, res);
        addLinearForces(src, res);
        addRigidity(src, res);

#ifdef ADD_PROJECTION_DIFF
        for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
        {
            Mecable const * mec = *mci;
            const index_type inx = DIM * mec->matIndex();
            mec->addProjectionDiff(src+inx, res+inx);
        }
#endif
        
        fwrite(res, sizeof(real), sz, file);
        src[ii] = 0.0;
    }
    
    delete [] res;
    delete [] src;
}


/**
 Save the projection matrix multiplied by the mobility
 */
void Meca::dumpProjection(FILE * file) const
{
    const unsigned sz = size();
    
    real * src = new real[sz];
    real * res = new real[sz];
    
    blas_xzero(sz, src);
    
    for ( unsigned ii = 0; ii < sz; ++ii )
    {
        src[ii] = 1.0;
        
        blas_xzero(sz, res);
        
        for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
        {
            Mecable const * mec = *mci;
            const index_type inx = DIM * mec->matIndex();
            mec->setSpeedsFromForces(src+inx, 1.0, res+inx);
        }
        
        fwrite(res, sizeof(real), sz, file);
        src[ii] = 0.0;
    }
    
    delete [] res;
    delete [] src;
}


/**
 Save matrix associated with the preconditionner, in binary format
 */
void Meca::dumpPreconditionner(FILE * file) const
{
    const unsigned sz = size();
    
    real * src = new real[sz];
    real * res = new real[sz];
    
    blas_xzero(sz, src);
    
    for ( unsigned ii = 0; ii < sz; ++ii )
    {
        src[ii] = 1.0;
        precondition(src, res);
        fwrite(res, sizeof(real), sz, file);
        src[ii] = 0.0;
    }
    
    delete [] res;
    delete [] src;
}


void Meca::dumpObjectID(FILE * file) const
{
    unsigned ms = 0;
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable const * mec = *mci;
        if ( ms < mec->nbPoints() )
            ms = mec->nbPoints();
    }
    
    real * vec = new real[ms];
    
    for ( unsigned ii = 0; ii < objs.size(); ++ii )
    {
        unsigned n = objs[ii]->nbPoints();
        for ( unsigned p=0; p < n; ++p )
            vec[p] = ii;
        for ( int d = 0; d < DIM; ++ d )
            fwrite(vec, sizeof(real), n, file);
    }
    
    delete[] vec;
}


void Meca::dumpDrag(FILE * file) const
{
    unsigned ms = 0;
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable const * mec = *mci;
        if ( ms < mec->nbPoints() )
            ms = mec->nbPoints();
    }
    
    real * vec = new real[ms];
    
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable const * mec = *mci;
        const real drag = mec->dragCoefficient() / mec->nbPoints();
        for ( unsigned p=0; p < mec->nbPoints(); ++p )
            vec[p] = drag;
        for ( int d = 0; d < DIM; ++ d )
            fwrite(vec, sizeof(real), mec->nbPoints(), file);
    }
    
    delete[] vec;
}


/**
 This dump the total matrix and some vectors in binary files.
 
 This matlab code should read the output:
 @code
 ord = load('cytosim_ord.txt');
 obj = fread(fopen('cytosim_obj.bin'), ord, 'double');
 drg = fread(fopen('cytosim_drg.bin'), ord, 'double');
 sys = fread(fopen('cytosim_sys.bin'), [ord, ord], 'double');
 ela = fread(fopen('cytosim_ela.bin'), [ord, ord], 'double');
 prj = fread(fopen('cytosim_prj.bin'), [ord, ord], 'double');
 con = fread(fopen('cytosim_con.bin'), [ord, ord], 'double');
 pts = fread(fopen('cytosim_pts.bin'), ord, 'double');
 rhs = fread(fopen('cytosim_rhs.bin'), ord, 'double');
 sol = fread(fopen('cytosim_sol.bin'), ord, 'double');
 @endcode
 
 Display the matrices:
 @code
 imshow(abs(sys))
 imshow(abs(ela))
 @endcode
 
 You can then compare the results with matlab's own iterative method,
 and compare the result using a scatter plot:
 
 @code
 x = bicgstab(sys, rhs, 0.001, 100);
 plot(x, sol, '.');
 @endcode
 */
void Meca::dump() const
{
    FILE * f = fopen("cytosim_ord.txt", "w");
    fprintf(f, "%i\n", size());
    fclose(f);
    
    f = fopen("cytosim_drg.bin", "wb");
    dumpDrag(f);
    fclose(f);
    
    f = fopen("cytosim_obj.bin", "wb");
    dumpObjectID(f);
    fclose(f);
    
    f = fopen("cytosim_rhs.bin", "wb");
    fwrite(vRHS, sizeof(real), size(), f);
    fclose(f);
    
    f = fopen("cytosim_sol.bin", "wb");
    fwrite(vSOL, sizeof(real), size(), f);
    fclose(f);
    
    f = fopen("cytosim_pts.bin", "wb");
    fwrite(vPTS, sizeof(real), size(), f);
    fclose(f);
    
    f = fopen("cytosim_sys.bin", "wb");
    dumpSystem(f);
    fclose(f);
    
    f = fopen("cytosim_ela.bin", "wb");
    dumpElasticity(f);
    fclose(f);
    
    f = fopen("cytosim_prj.bin", "wb");
    dumpProjection(f);
    fclose(f);
    
    f = fopen("cytosim_con.bin", "wb");
    dumpPreconditionner(f);
    fclose(f);
    
    std::clog << " Dumped a system of size " << size() << std::endl;
}


/**
 output of matrices in a text-based sparse format
 */
void Meca::dumpSparse()
{
#ifndef RIGIDITY_IN_MATRIX
    std::clog << "dump is incorrect if RIGIDITY_IN_MATRIX is not defined";
    return;
#endif
    std::clog << "dumping matrices for a system of size " << nbPts << std::endl;
    
    unsigned ms = 0;
    for ( Mecable ** mci = objs.begin(); mci < objs.end(); ++mci )
    {
        Mecable const * mec = *mci;
        if ( ms < mec->nbPoints() )
            ms = mec->nbPoints();
    }
    
    FILE * f = fopen("d_drg.bin", "wb");
    dumpDrag(f);
    fclose(f);
    
    f = fopen("d_obj.bin", "wb");
    dumpObjectID(f);
    fclose(f);
    
    std::ofstream os("d_sol.txt");
    VecPrint::vecDump(os, DIM*nbPts, vPTS);
    os.close();
    
    os.open("d_rhs.txt");
    VecPrint::vecDump(os, DIM*nbPts, vRHS);
    os.close();
    
    os.open("d_matB.txt");
    mB.printSparse(os);
    os.close();
    
    os.open("d_matC.txt");
    mC.printSparse(os);
    os.close();
    
    
    real* tmp1 = new real[DIM*ms];
    real* tmp2 = new real[DIM*ms*DIM*ms];
    
    os.open("diagonal.txt");
    
    for ( unsigned int oo = 0; oo < objs.size(); ++oo )
    {
        unsigned int bs = DIM * objs[oo]->nbPoints();
        extractBlockS(tmp2, objs[oo]);
        VecPrint::matSparsePrintOffset(os, bs, bs, tmp2, bs, DIM*objs[oo]->matIndex() );
    }
    os.close();
    
    delete [] tmp1;
    delete [] tmp2;
}


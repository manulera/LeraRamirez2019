// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 18/12/07.

#include "field.h"
#include "fiber_binder.h"
#include "fiber_set.h"
#include "cblas.h"
#include "sim.h"


template < >
real FieldBase<FieldScalar>::display_amp = 0;
template < >
Space const* FieldBase<FieldScalar>::display_spc = 0;


extern Random RNG;

/**
 Initialize the diffusion matrix using periodic boundary conditions
 if the underlying space is peridic
 */
template < >
void FieldBase<FieldScalar>::prepareDiffusion(real theta)
{
    const unsigned nbc = FieldGrid::nbCells();
    
    fiDiffusionMatrix.resize(nbc);
    fiDiffusionMatrix.makeZero();
    
    for ( unsigned c = 0; c < nbc; ++c )
    {
        for ( int d = 0; d < DIM; ++d )
        {
            unsigned n = FieldGrid::next(c, d);
                
            if ( n != c )
            {
                fiDiffusionMatrix(c, n) += theta;
                fiDiffusionMatrix(c, c) -= theta;
                fiDiffusionMatrix(n, n) -= theta;
            }
        }
    }
    fiDiffusionMatrix.prepareForMultiply();
    
    /*
    std::clog << "tight Field has diffusion matrix with ";
    std::clog << fiDiffusionMatrix.nbElements() << " elements" << std::endl;
     */
}


/**
 Initialize the diffusion matrix.
 Diffusion is allowed between neighboring cells that are in the same domain:
 @code
 ( domain[c] > 0 ) && ( domain[c] == domain[n] )
 @endcode
 */
template < >
void FieldBase<FieldScalar>::prepareDiffusion(real theta, unsigned char * domain)
{
    const unsigned nbc = FieldGrid::nbCells();
    
    fiDiffusionMatrix.resize(nbc);
    fiDiffusionMatrix.makeZero();
    
    for ( unsigned c = 0; c < nbc; ++c )
    {
        if ( domain[c] )
        {
            for ( int d = 0; d < DIM; ++d )
            {
                unsigned n = c + FieldGrid::stride(d);
                
                if ( n < nbc  &&  domain[c] == domain[n] )
                {
                    fiDiffusionMatrix(c, n) += theta;
                    fiDiffusionMatrix(c, c) -= theta;
                    fiDiffusionMatrix(n, n) -= theta;
                }
            }
        }
    }
    fiDiffusionMatrix.prepareForMultiply();
    
    /*
     std::clog << "Field has diffusion matrix with ";
     std::clog << fiDiffusionMatrix.nbElements() << " elements" << std::endl;
     */
}







/**
 Initialize Field to be ready for step()
 */
template < > 
void FieldBase<FieldScalar>::prepare()
{
    const Space * spc = prop->confine_space_ptr;

    if ( spc == 0 )
        throw InvalidParameter("A Space must be created before the field");

    const unsigned nbc = FieldGrid::nbCells();
    assert_true( nbc > 0 );
    
    if ( fiTMP )
        delete [] fiTMP;
    fiTMP = new real[nbc];
    fiTMPSize = nbc;

    if ( prop->diffusion > 0 )
    {
        real theta = prop->diffusion * prop->time_step / ( prop->step * prop->step );

        if ( DIM == 1 || prop->periodic )
            prepareDiffusion(theta);
        else
        {
            unsigned char * domain = new unsigned char[nbc];
            
            // determine which cell is inside the space:
#if ( 1 )
            for ( unsigned c = 0; c < nbc; ++c )
            {
                Vector pos;
                setPositionFromIndex(pos, c, 0.5);
                domain[c] = spc->inside(pos);
            }
#else
            // extended covered area:
            const real range = 2 * cellWidth();
            for ( unsigned c = 0; c < nbc; ++c )
            {
                Vector pos;
                setPositionFromIndex(pos, c, 0.5);
                domain[c] = ! spc->allOutside(pos, range);
            }
#endif
            
            prepareDiffusion(theta, domain);
            
            delete[] domain;
        }
    }
}



#pragma mark -


template < >
void FieldBase<FieldScalar>::diffuseX(real * field, real c)
{
    const int nbc = FieldGrid::nbCells();

    unsigned nx = FieldGrid::dim(0);
    unsigned nyz = nbc / nx;
    unsigned stride = FieldGrid::stride(1);
    
    real * a = new real[nyz];
    real * b = new real[nyz];
    
    // diffusion in X-direction:
    blas_xzero(nyz, a);
    for ( unsigned x = 1; x < nx; ++x )
    {
        real * h = field + x - 1;
        real * n = field + x;
        // b = n - h
        blas_xcopy(nyz,  n, stride, b, 1);
        blas_xaxpy(nyz, -1, h, stride, b, 1);
        // a = a - b
        blas_xaxpy(nyz, -1, b, 1, a, 1);
        // h = h - c * a
        blas_xaxpy(nyz, -c, a, 1, h, stride);
        // swap a and b
        real * t = a;
        a = b;
        b = t;
    }
    real * h = field + nx - 1;
    blas_xaxpy(nyz, -c, a, 1, h, stride);
    
    if ( prop->periodic )
    {
        real * n = field;
        blas_xcopy(nyz,  n, stride, b, 1);
        blas_xaxpy(nyz, -1, h, stride, b, 1);
        blas_xaxpy(nyz,  c, b, 1, h, stride);
        blas_xaxpy(nyz, -c, b, 1, n, stride);
    }
    
    delete[] a;
    delete[] b;
}


template < >
void FieldBase<FieldScalar>::laplacian(const real* field, real * mat) const
{
    const int nbc = FieldGrid::nbCells();
    
    const real sc = 2 * DIM;

    for ( int c = 0; c < nbc; ++c )
        mat[c] = sc * field[c];
    
    const int nx = FieldGrid::dim(0);
#if ( 1 )
    // derivative in the X-direction:
    const int nyz = nbc / nx;
    for ( int xx = 1; xx < nx; ++xx )
    {
        blas_xaxpy(nyz, -1, field+xx-1, nx, mat+xx  , nx);
        blas_xaxpy(nyz, -1, field+xx  , nx, mat+xx-1, nx);
    }
    // index of last valid X index:
    int xx = FieldGrid::dim(0) - 1;
    
    if ( prop->periodic )
    {
        blas_xaxpy(nyz, -1, field+xx, nx, mat   , nx);
        blas_xaxpy(nyz, -1, field   , nx, mat+xx, nx);
    }
    else
    {
        blas_xaxpy(nyz, -1, field   , nx, mat   , nx);
        blas_xaxpy(nyz, -1, field+xx, nx, mat+xx, nx);
    }
#endif
    
#if ( DIM == 2 )
    // derivative in the Y-direction:
    blas_xaxpy(nbc-nx, -1, field,    1, mat+nx, 1);
    blas_xaxpy(nbc-nx, -1, field+nx, 1, mat,    1);
    
    int yy = FieldGrid::dim(1) - 1;
    if ( prop->periodic )
    {
        blas_xaxpy(nx, -1, field+nx*yy, 1, mat      , 1);
        blas_xaxpy(nx, -1, field      , 1, mat+nx*yy, 1);
    }
    else
    {
        blas_xaxpy(nx, -1, field      , 1, mat      , 1);
        blas_xaxpy(nx, -1, field+nx*yy, 1, mat+nx*yy, 1);
    }
#endif

#if ( DIM == 3 )
    // derivative in the Y-direction:
    const int sz = FieldGrid::stride(2);
    for ( int yy = 1; yy < FieldGrid::dim(1); ++yy )
    for ( int zz = 0; zz < FieldGrid::dim(2); ++zz )
    {
        blas_xaxpy(nx, -1, field+nx*(yy-1)+sz*zz, 1, mat+nx*(yy  )+sz*zz, 1);
        blas_xaxpy(nx, -1, field+nx*(yy  )+sz*zz, 1, mat+nx*(yy-1)+sz*zz, 1);
    }
    int yy = FieldGrid::dim(1) - 1;
    
    if ( prop->periodic )
    {
        for ( int zz = 0; zz < FieldGrid::dim(2); ++zz )
        {
            blas_xaxpy(nx, -1, field+nx*yy+sz*zz, 1, mat      +sz*zz, 1);
            blas_xaxpy(nx, -1, field      +sz*zz, 1, mat+nx*yy+sz*zz, 1);
        }
    }
    else
    {
        for ( int zz = 0; zz < FieldGrid::dim(2); ++zz )
        {
            blas_xaxpy(nx, -1, field      +sz*zz, 1, mat      +sz*zz, 1);
            blas_xaxpy(nx, -1, field+nx*yy+sz*zz, 1, mat+nx*yy+sz*zz, 1);
        }
    }
#endif

#if ( DIM == 3 )
    // derivative in the Z-direction:
    const int nxy = nbc / FieldGrid::dim(2);
    assert_true( nxy == sz );
    blas_xaxpy(nbc-nxy, -1, field,     1, mat+nxy, 1);
    blas_xaxpy(nbc-nxy, -1, field+nxy, 1, mat,     1);
    int zz = FieldGrid::dim(2) - 1;
    
    if ( prop->periodic )
    {
        blas_xaxpy(nxy, -1, field+sz*zz, 1, mat      , 1);
        blas_xaxpy(nxy, -1, field      , 1, mat+sz*zz, 1);
    }
    else
    {
        blas_xaxpy(nxy, -1, field      , 1, mat      , 1);
        blas_xaxpy(nxy, -1, field+sz*zz, 1, mat+sz*zz, 1);
    }
#endif
}



template < >
void FieldBase<FieldScalar>::setEdgesX(real * field, real val)
{
    const int nbc = FieldGrid::nbCells();
    
    // set X-edges:
    const int nx = FieldGrid::dim(0);
    
    real * lastf = field + nx - 1;
    for ( int xx = 0; xx < nbc; xx += nx )
    {
        field[xx] = val;
        lastf[xx] = val;
    }
}


template < >
void FieldBase<FieldScalar>::setEdgesY(real * field, real val)
{
#if ( DIM > 1 )
    const int nbc = FieldGrid::nbCells();
    const int nx = FieldGrid::dim(0);
#endif
    
#if ( DIM == 2 )
    // set Y-edges:
    real * lastf = field + nbc - nx;
    for ( int xx = 0; xx < nx; ++xx )
    {
        field[xx] = val;
        lastf[xx] = val;
    }
#endif
    
#if ( DIM == 3 )
    // set Y-edges:
    const int nz = FieldGrid::dim(2);
    const int nxy = nbc / nz;
    
    real * lastf = field + nxy - nx;
    for ( int zz = 0; zz < nz; ++zz )
    {
        for ( int xx = 0; xx < nx; ++xx )
        {
            field[xx+zz*nxy] = val;
            lastf[xx+zz*nxy] = val;
        }
    }
#endif
}


template < >
void FieldBase<FieldScalar>::setEdgesZ(real * field, real val)
{
#if ( DIM == 3 )
    const int nbc = FieldGrid::nbCells();
    
    const int nz = FieldGrid::dim(2);
    const int nxy = nbc / nz;
    
    real * lastf = field + nxy * ( nz - 1 );
    for ( int xy = 0; xy < nxy; ++xy )
    {
        field[xy] = val;
        lastf[xy] = val;
    }
#endif
}



/**
 //\todo implement Crank-Nicholson for diffusion
 */
template < >
void FieldBase<FieldScalar>::step(FiberSet& fibers)
{
    assert_true( prop );
    
    // we cast FieldScalar to floating-point type :
    assert_true( sizeof(FieldScalar) == sizeof(real) );
    real * field = reinterpret_cast<real*>(cell_data());
    const int nbc = FieldGrid::nbCells();
    
    real * dup = fiTMP;

    // decay:
    if ( prop->decay_rate > 0 )
    {
        // field = field * exp( - decay_rate * dt ):
        blas_xscal(nbc, prop->decay_frac, field, 1);
    }

    // full grid diffusion:
    if ( prop->full_diffusion > 0 )
    {
        real c = prop->full_diffusion * prop->time_step / ( prop->step * prop->step );

#if ( DIM > 1 )
        laplacian(field, dup);
        blas_xaxpy(nbc, -c, dup, 1, field, 1);
#else
        diffuseX(field, c);
#endif
    }

    // diffusion:
    if ( prop->diffusion > 0 )
    {
        assert_true( fiTMP );
        assert_true( fiTMPSize == nbc );
        assert_true( fiDiffusionMatrix.size() == nbc );

        // dup = field:
        blas_xcopy(nbc, field, 1, dup, 1);
        
        // field = field + fiDiffusionMatrix * dup:
        fiDiffusionMatrix.vecMulAdd(dup, field);
    }

    if ( prop->boundary_condition & 1 )
        setEdgesX(field, prop->boundary_value * cellVolume());
    
#if ( DIM > 1 )
    if ( prop->boundary_condition & 2 )
        setEdgesY(field, prop->boundary_value * cellVolume());
#endif
    
#if ( DIM == 3 )
    if ( prop->boundary_condition & 4 )
        setEdgesZ(field, prop->boundary_value * cellVolume());
#endif
    
#if ( 0 )

    Array<FiberBinder> loc(1024);
    
    // binding along Fibers
    /**
     This functionality is replaced by fiber:lattice_binding_rate,
     in which a different field can be specified for each type of fiber
     */
    if ( prop->bind_fibers > 0 )
    {
        // we want roughly one point per cell:
        const real spread = 0.5 * cellWidth();
        // each point represents a Fiber chunk of length 'spread':
        const real rate = prop->bind_fibers * spread / cellVolume();
        // fraction of the cell content that will bind:
        const real frac = 1 - exp( -rate * prop->time_step );
        
        fibers.uniFiberSites(loc, spread);

        for ( Array<FiberBinder>::iterator i = loc.begin(); i < loc.end(); ++i )
        {
            FiberBinder & site = *i;
            Fiber* fib = site.fiber();
            FiberLattice* lat = fib->lattice();
            
            value_type & cell = FieldGrid::cell(site.pos());
            
            // amount to be transfered:
            real mass = cell * frac;
            
            cell -= mass;
            lat->site(site.abscissa()) += mass;
        }
        
        //std::clog << "field sum = " << sumValues() << std::endl;
    }
    
    // instantaneous transport along Fibers
    if ( prop->transport_strength > 0 )
    {
        const real spread = 0.5 * cellWidth();
        const real rate = prop->transport_strength * spread / cellVolume();
        const real frac = 1 - exp( -rate * prop->time_step );
        
        if ( frac >= 0.5 )
            throw InvalidParameter("field:transport_strength is too high");
        
        fibers.uniFiberSites(loc, spread);
        for ( Array<FiberBinder>::iterator i = loc.begin(); i < loc.end(); ++i )
        {
            FiberBinder & site = *i;
            
            // abscissa for exit point of transport:
            real abs = site.abscissa() + RNG.exponential(prop->transport_length);

            // find index of cell:
            value_type cell = FieldGrid::cell(site.pos());
            
            // amount to be transfered:
            real mass = cell * frac;
            
            // transport:
            cell -= mass;
            field[FieldGrid::index(site.fiber()->pos(abs))] += mass;
        }
    }
    
    // direct cutting of fiber
    // this is deprecated in favor of fiber:lattice_cut_fiber
    if ( prop->cut_fibers )
    {
        PRINT_ONCE("!!!! Field severs fibers\n");
        const real spread = 0.5 / prop->time_step;
        const real fac = spread * prop->time_step / cellVolume();
        
        fibers.uniFiberSites(loc, spread);
        for ( Array<FiberBinder>::iterator i = loc.begin(); i < loc.end(); ++i )
        {
            FiberBinder & site = *i;
            real val = field[FieldGrid::index(site.pos())];
            if ( prop->cut_fibers == 2 )
                val = val * val / cellVolume();
            if ( RNG.test_not( exp(-fac*val) ) )
                site.fiber()->sever(site.abscissa(), STATE_RED, STATE_GREEN);
        }
    }
    
    if ( prop->chew_fibers )
    {
        PRINT_ONCE("!!!! Field chews PLUS_END\n");
        const real fac = -prop->time_step / cellVolume();
        for ( Fiber * fib = fibers.first(); fib ; fib = fib->next() )
            fib->growP(fac*cell(fib->posEndP()));
    }
#endif
}


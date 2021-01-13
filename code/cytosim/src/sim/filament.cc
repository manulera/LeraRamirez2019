// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "assert_macro.h"
#include "filament.h"
#include "iowrapper.h"
#include "messages.h"
#include "point_exact.h"
#include "point_interpolated.h"
#include "fiber_binder.h"
#include "exceptions.h"
#include "clapack.h"
#include "modulo.h"

extern Random RNG;
extern Modulo const* modulo;

/**
 This returns N+1, where N is the integer that minimizes:
 fabs( length / N - segmentation ),
 */
unsigned Filament::bestNumberOfPoints(const real ratio)
{
    unsigned n = (int)ratio;
    
    if ( (2*n+1)*ratio > 2*n*(n+1) )
        return n+2;
    
    return n+1;
}


Filament::Filament()
{
    normal.set(0, 0, 0);
    fnCut           = 0;
    fnSegmentation  = 0;
    fnAbscissa      = 0;
#ifdef CURVATURE_DEPENDENT_SEGMENTATION
    fnCutError      = 0;
    ///giving a different 'seed' desynchronize the adjustSegmentation() operations
    fnCutErrorIndex = RNG.pint(RECUT_PERIOD);
#endif
    needUpdate      = false;
}



//------------------------------------------------------------------------------
#pragma mark -

void Filament::setStraight(Vector const& pos, Vector const& dir)
{
    assert_true( dir.norm() > 0.1 );
    // 'dir' is normalized for safety:
    Vector dpts = dir * ( fnCut / dir.norm() );
    //
    for ( unsigned int p = 0 ; p < nbPoints(); ++p )
        setPoint( p, pos + p * dpts );
}


void Filament::setStraight(Vector const& pos, Vector const& dir, const FiberEnd ref)
{
    switch( ref )
    {
        case MINUS_END:
            setStraight( pos, dir );
            break;
            
        case PLUS_END:
            setStraight( pos + dir*length(), -dir );
            break;
            
        case CENTER:
            setStraight( pos - 0.5*dir*length(), dir );
            break;
            
        default:
            ABORT_NOW("invalid argument `ref`");
    }
}


void Filament::setStraight(Vector const& pos, Vector const& dir, real len, const FiberEnd ref)
{
    assert_true( fnSegmentation > REAL_EPSILON );

    if ( len <= 0 )
        throw InvalidParameter("fiber:length must be > 0");

    int nbp = bestNumberOfPoints(len/fnSegmentation);
    assert_true( nbp > 1 );
    
    fnCut = len / real(nbp-1);
    setNbPoints(nbp);
    
    setStraight(pos, dir, ref);
    postUpdate();
}


real Filament::trueLength(const real* pts, unsigned n_pts)
{
    real len = 0;
    Vector a(pts), b;
    for ( unsigned n = 1; n < n_pts; ++n )
    {
        b.get(pts+DIM*n);
        len += (b-a).norm();
        a = b;
    }
    return len;
}

/**
 This will set the Fiber with `np` points unless `np == 0`, in which case
 the number of points will be set automatically from fnSegmentation.
 pts[] should be of size DIM * n_pts and contain coordinates.

 The given set of points do not need to be equally distributed.
 The MINUS_END and PLUS_END will be set to the first and last points in `pts[]`,
 and intermediate points will be interpolated at regular intervals on `pts[]`.
 
 The length of the resulting fiber will be roughly equal to the sum of all segment lengths.
 However, the length of the segments will only be approximately equal to each other,
 and reshape() should be called to equalize them if necessary.
 */
void Filament::setShape(const real pts[], unsigned n_pts, unsigned np)
{
    assert_true(n_pts > 1);
    Vector a(pts), b;
    
    //calculate the total length
    real len = trueLength(pts, n_pts);
    
    if ( np == 0 )
    {
        assert_true( fnSegmentation > REAL_EPSILON );
        np = bestNumberOfPoints(len/fnSegmentation);
    }
    fnCut = len / real(np-1);
    setNbPoints(np);
    
    a.get(pts);
    b.get(pts+DIM);
    setPoint(0, a);
    
    len = (b-a).norm();
    real h = 0;
    unsigned p = 1;
    --np;
    
    for ( unsigned n = 1; n < np; ++n )
    {
        h += fnCut;

        while ( h > len )
        {
            h -= len;
            a = b;
            ++p;
            assert_true(p<n_pts);
            b.get(pts+DIM*p);
            len = (b-a).norm();
        }
        
        setPoint(n, a+(h/len)*(b-a));
    }
    b.get(pts+DIM*n_pts-DIM);
    setPoint(np, b);
    postUpdate();
}


//===================================================================
#pragma mark -

/*
 This deals with Fiber having one segment only,
 for which the procedure is trivial
 */
void Filament::reshape_two(const real* src, real* dst, real cut)
{
#if ( DIM == 1 )
        Vector dif(src[1]-src[0], 0);
#elif ( DIM == 2 )
        Vector dif(src[2]-src[0], src[3]-src[1]);
#else
        Vector dif(src[3]-src[0], src[4]-src[1], src[5]-src[2]);
#endif
        dif *= 0.5 * ( 1 - cut/dif.norm() );
        
        dst[0    ] = src[0    ] + dif.XX;
        dst[  DIM] = src[  DIM] - dif.XX;
#if ( DIM > 1 )
        dst[1    ] = src[1    ] + dif.YY;
        dst[1+DIM] = src[1+DIM] - dif.YY;
#endif
#if ( DIM > 2 )
        dst[2    ] = src[2    ] + dif.ZZ;
        dst[2+DIM] = src[2+DIM] - dif.ZZ;
#endif
}



/**
 Provide a reasonable initial guess of the coefficients
 */
void Filament::reshape_guess(real * sca, const unsigned ns, const Vector* dif, real cut,
                               real * ang, real * tmp)
{
    sca[0] = 0;
    for ( unsigned pp = 1; pp < ns; ++pp )
    {
        sca[pp] = 0;
        ang[pp] = ( dif[pp-1] * dif[pp] ) / sqrt( dif[pp-1].normSqr() * dif[pp].normSqr() );
    }
    
    for ( unsigned pp = 0; pp < ns; ++pp )
    {
        printf("\npp %2i  ", pp);
        
        real L = 1;
        for ( int qq = pp; qq > 0; --qq )
        {
            tmp[qq] = L;
            L *= ang[qq];
        }
        tmp[0] = L;
        printf("  LL %8.4f   ", L);
        for ( unsigned qq = 1; qq <= pp; ++qq )
        {
            L = tmp[qq] + L * ang[qq];
            tmp[qq] = L;
        }
        
        real R = 1;
        tmp[pp] = R;
        for ( unsigned qq = pp+1; qq < ns; ++qq )
        {
            R *= ang[qq];
            tmp[qq] = R;
        }
        printf("   RR %8.4f   ", R);
        for ( unsigned qq = ns-2; qq >= pp; --qq )
        {
            R = tmp[qq] + R * ang[qq+1];
            tmp[qq] = R;
        }
        
        real e = ( dif[pp].norm() - cut );
        real s = e * ( R * L ) / ( L + R );
        printf("\n L %8.4f  R %8.4f  s %8.4f  ", L, R, s);
        
        sca[pp] += s;
        
        L = s / L;
        for ( unsigned qq = 0; qq < pp; ++qq )
            sca[qq] += L * tmp[qq];
        
        R = s / R;
        for ( unsigned qq = pp+1; qq < ns; ++qq )
            sca[qq] += R * tmp[qq];
        
        printf("   %i sca  ", pp);
        for ( unsigned n = 0; n < ns; ++n )
            printf("%+6.4f ", sca[n]);
    }
    
    for ( unsigned pp = 0; pp < ns; ++pp )
        sca[pp] /= dif[pp].norm();
    
#if ( 1 )
    printf("\n---> sca  ");
    for ( unsigned n = 0; n < ns; ++n )
        printf("%+6.4f ", sca[n]);
    printf("\n");
#endif
}


/**
 Apply correction
 */
void Filament::reshape_apply(const unsigned ns, const real* src, real* dst,
                             const real * sca, const Vector* dif)
{
    assert_true( ns > 1 );
    Vector d, e = sca[0] * dif[0];
    dst[0] = src[0] + e.XX;
#if ( DIM > 1 )
    dst[1] = src[1] + e.YY;
#endif
#if ( DIM > 2 )
    dst[2] = src[2] + e.ZZ;
#endif
    for ( unsigned pp = 1; pp < ns; ++pp )
    {
        d = sca[pp] * dif[pp];
        dst[DIM*pp  ] = src[DIM*pp  ] + d.XX - e.XX;
#if ( DIM > 1 )
        dst[DIM*pp+1] = src[DIM*pp+1] + d.YY - e.YY;
#endif
#if ( DIM > 2 )
        dst[DIM*pp+2] = src[DIM*pp+2] + d.ZZ - e.ZZ;
#endif
        e = d;
    }
    dst[DIM*ns  ] = src[DIM*ns  ] - d.XX;
#if ( DIM > 1 )
    dst[DIM*ns+1] = src[DIM*ns+1] - d.YY;
#endif
#if ( DIM > 2 )
    dst[DIM*ns+2] = src[DIM*ns+2] - d.ZZ;
#endif
}


/**
 Shorten segments to restore their length to 'cut'.
 We use a multidimensional Newton's method, to find iteratively the scalar
 coefficients that define the amount of displacement of each point.
 
 X[i] = vector of position
 We note 'dif' the differences between consecutive points:  dif[i] = X[i+1] - X[i]
 
 Given one scalar per segment: sca[i], the point is displaced as:
 Y[i] = X[i] + sca[i] * dif[i] - sca[i-1] * dif[i-1]
 except for the first and last points, for which there is only one term:
 Y[0] = X[0] + sca[  0] * dif[  0]
 Y[L] = X[L] - sca[L-1] * dif[L-1]
 
 We want 'sca' to restore the length of segments:
 ( Y[i+1] - Y[i] )^2 = cut^2
 
 i.e. 'sca' should fulfill a set of equalities F[i] = 0, with:
 F[i] = ( Y[i+1] - Y[i] )^2 - cut^2
 
 Method: use all zeros as first guess for 'sca', and apply a multidimensional
 Newton's method to iteratively refine the guess.
 
 In practice, we calculate `sca_next` from `sca` using the relationship:
 J(sca) * ( sca_next - sca ) = -F(sca)
 
 Where J is the Jacobian matrix: J[i,j] = dF[i] / dX[j]
 
 For this problem, J is square and tri-diagonal but not symmetric,
 and must be recalculated at each iteration.

 FJN, Strasbourg, 22 Feb 2015
 */
#if ( 1 )
int Filament::reshape_it(const unsigned ns, const real* src, real* dst, real cut)
{
    assert_true( ns > 1 );
    int info = 0;
    unsigned ii = 0;
    real err = 0;
    const real alphaSqr = cut * cut;
    
    Vector * dif = new Vector[ns];
    
    // Keep memory aligned to 64 bytes:
    const size_t chunk = 64 / sizeof(real);
    // make a multiple of chunk to align memory to 32 bytes:
    size_t chk = ( ns + chunk - 1 ) & ~( chunk -1 );

    //std::clog << "fiber::reshape_it allocates " << chk << std::endl;
    real * mem = new real[chk*5];
    real * sca = mem;
    real * val = mem+ns;
    real * dia = mem+ns*2;
    real * low = mem+ns*3;
    real * upe = mem+ns*4;

    // calculate differences:
    {
    assert_true( sizeof(Vector) == DIM*sizeof(real) );
    real * dif_ = dif->data();
    const int end = DIM * ns;
    for ( int p = 0; p < end; ++p )
        dif_[p] = src[p+DIM] - src[p];
    }
        
#if ( 1 )
    /*
     Perform here the first iteration of Newton's method
     the formula is the same as below, with all `sca` equal to zero,
     and thus 'vec == dif'
     The system is symmetric, and we can use a faster factorization
     */
    val[0] = dif[0].normSqr() - alphaSqr;
    dia[0] = 2 * dif[0].normSqr();
    for ( unsigned pp = 1; pp < ns; ++pp )
    {
        real n = dif[pp].normSqr();
        val[pp] = n - alphaSqr;
        low[pp] = -( dif[pp] * dif[pp-1] );
        dia[pp] = 2 * n;
    }
    
    lapack_xpttrf(ns, dia, low+1, &info);
    if ( info ) {
        std::cerr << " LAPACK dpttrf failed " << info << std::endl;
        goto finish;
    }
    lapack_xptts2(ns, 1, dia, low+1, val, ns);
    
    for ( unsigned pp = 0; pp < ns; ++pp )
    {
        sca[pp] = 0.5 * val[pp];
        err += fabs(val[pp]);
    }
#else
    // start with a naive guess
    err = INFINITY;
    reshape_guess(sca, ns, dif, cut, low, upe);
#endif
#if ( 0 )
    printf("\n --- sca ");
    for ( unsigned pp = 0; pp < ns; ++pp )
        printf(" %+8.6f", sca[pp]);
#endif

    while ( err > 1e-12 )
    {
        assert_true( ns > 1 );
        // set the matrix elements and RHS of system,
        // calculating 'vec' on the fly
        Vector vec0 = (1-2*sca[0])*dif[0] + sca[1]*dif[1];
        val[0] = vec0.normSqr() - alphaSqr;
        dia[0] = -2 * ( vec0 * dif[0] );
        upe[0] = vec0 * dif[1];
        int pp = 1;
        while ( pp+1 < ns )
        {
            vec0 = sca[pp-1]*dif[pp-1] + (1-2*sca[pp])*dif[pp] + sca[pp+1]*dif[pp+1];
            val[pp] = vec0.normSqr() - alphaSqr;
            low[pp] = vec0 * dif[pp-1];
            dia[pp] = -2 * ( vec0 * dif[pp] );
            upe[pp] = vec0 * dif[pp+1];
            ++pp;
        }
        assert_true( pp == ns-1 );
        vec0 = sca[pp-1]*dif[pp-1] + (1-2*sca[pp])*dif[pp];
        val[pp] = vec0.normSqr() - alphaSqr;
        low[pp] = vec0 * dif[pp-1];
        dia[pp] = -2 * ( vec0 * dif[pp] );

#if ( 0 )
        real sum = 0;
        for ( unsigned pp = 0; pp < ns; ++pp )
            sum += val[pp];
        printf("\n   %i sum %8.5f", ii, sum);
#endif
#if ( 0 )
        printf("\n   %i val  ", ii);
        for ( unsigned pp = 0; pp < ns; ++pp )
            printf("%+6.4f ", val[pp]);
        printf("\n   %i upe  ", ii);
        for ( unsigned pp = 0; pp+1 < ns; ++pp )
            printf("%+6.4f ", upe[pp]);
        printf("\n   %i dia  ", ii);
        for ( unsigned pp = 0; pp < ns; ++pp )
            printf("%+6.4f ", dia[pp]);
        printf("\n   %i low  ", ii);
        for ( unsigned pp = 1; pp < ns; ++pp )
            printf("%+6.4f ", low[pp]);
#endif
        
        lapack_xgtsv(ns, 1, low+1, dia, upe, val, ns, &info);
        if ( info )
        {
            std::cerr << " LAPACK dgtsv failed " << info << std::endl;
            goto finish;
        }

        // calculate residual error
        err = 0;
        for ( unsigned pp = 0; pp < ns; ++pp )
        {
            sca[pp] -= 0.5 * val[pp];
            err += fabs(val[pp]);
        }
        if ( ++ii > 31 )
        {
            info = 1;
            goto finish;
        }
        
#if ( 0 )
        printf("\n %3i sca ", ii);
        for ( unsigned pp = 0; pp < ns; ++pp )
            printf(" %+8.6f", sca[pp]);
#endif
    }
    
#if ( 0 )
    printf("\n%2i err %e", ii, err);
    printf("\n%2i sca  ", ii);
    for ( unsigned pp = 0; pp < ns; ++pp )
        printf("%+6.4f ", sca[pp]);
    printf("\n");
#endif
    
    //apply corrections:
    reshape_apply(ns, src, dst, sca, dif);
    
finish:
    delete [] mem;
    delete [] dif;
    
    return info;
}

#else

int Filament::reshape_it(const unsigned ns, const real* src, real* dst, real cut)
{
    assert_true( ns > 1 );
    int info = 0;
    const real alphaSqr = cut * cut;
    
    Vector * dif = new Vector[ns];
    Vector * vec = new Vector[ns];
    real * sca = new real[ns];
    real * val = new real[ns];
    
    real * dia = new real[ns];
    real * low = new real[ns];
    real * upe = new real[ns];
    
    // calculate differences
    for ( unsigned pp = 0; pp < ns; ++pp )
    {
        dif[pp] = diffPoints(src, pp);
        sca[pp] = 0;
    }
    
    real err = 0;
    unsigned ii = 0;
    do {
#if ( 0 )
        printf("\n   %i sca  ", ii);
        for ( unsigned pp = 0; pp < ns; ++pp )
            printf("%+6.4f ", sca[pp]);
#endif

        // calculate all values of 'vec'
        vec[0] = (1-2*sca[0])*dif[0] + sca[1]*dif[1];
        for ( unsigned pp = 1; pp+1 < ns; ++pp )
            vec[pp] = sca[pp-1]*dif[pp-1] + (1-2*sca[pp])*dif[pp] + sca[pp+1]*dif[pp+1];
        vec[ns-1] = sca[ns-2]*dif[ns-2] + (1-2*sca[ns-1])*dif[ns-1];
        
        // calculate the matrix elements and RHS of system
        val[0] = vec[0].normSqr() - alphaSqr;
        dia[0] = -2 * ( vec[0] * dif[0] );
        for ( unsigned pp = 1; pp < ns; ++pp )
        {
            val[pp] = vec[pp].normSqr() - alphaSqr;
            low[pp] = vec[pp] * dif[pp-1];
            dia[pp] = -2 * ( vec[pp] * dif[pp] );
            upe[pp-1] = vec[pp-1] * dif[pp];
        }
        
#if ( 0 )
        printf("\n   %i val  ", ii);
        for ( unsigned pp = 0; pp < ns; ++pp )
            printf("%+6.4f ", val[pp]);
        printf("\n   %i upe  ", ii);
        for ( unsigned pp = 0; pp+1 < ns; ++pp )
            printf("%+6.4f ", upe[pp]);
        printf("\n   %i dia  ", ii);
        for ( unsigned pp = 0; pp < ns; ++pp )
            printf("%+6.4f ", dia[pp]);
        printf("\n   %i low  ", ii);
        for ( unsigned pp = 1; pp < ns; ++pp )
            printf("%+6.4f ", low[pp]);
#endif
        
        lapack_xgtsv(ns, 1, low+1, dia, upe, val, ns, &info);
        if ( info )
        {
            std::cerr << " LAPACK dgtsv failed " << info << std::endl;
            goto finish;
        }
        
        err = 0;
        for ( unsigned pp = 0; pp < ns; ++pp )
        {
            sca[pp] += -0.5 * val[pp];
            err += fabs(val[pp]);
        }
        if ( ++ii > 32 )
        {
            info = 1;
            goto finish;
        }
    } while ( err > 0.0001 );

    
#if ( 0 )
    printf("\n%2i err %e", ii, err);
    printf("\n%2i sca  ", ii);
    for ( unsigned pp = 0; pp < ns; ++pp )
        printf("%+6.4f ", sca[pp]);
    printf("\n");
#endif
    
    //apply corrections:
    {
        Vector d, e = sca[0] * dif[0];
        dst[0] = src[0] + e.XX;
#if ( DIM > 1 )
        dst[1] = src[1] + e.YY;
#endif
#if ( DIM > 2 )
        dst[2] = src[2] + e.ZZ;
#endif
        for ( unsigned pp = 1; pp < ns; ++pp )
        {
            d = sca[pp] * dif[pp];
            dst[DIM*pp  ] = src[DIM*pp  ] + d.XX - e.XX;
#if ( DIM > 1 )
            dst[DIM*pp+1] = src[DIM*pp+1] + d.YY - e.YY;
#endif
#if ( DIM > 2 )
            dst[DIM*pp+2] = src[DIM*pp+2] + d.ZZ - e.ZZ;
#endif
            e = d;
        }
        dst[DIM*ns+0] = src[DIM*ns+0] - d.XX;
#if ( DIM > 1 )
        dst[DIM*ns+1] = src[DIM*ns+1] - d.YY;
#endif
#if ( DIM > 2 )
        dst[DIM*ns+2] = src[DIM*ns+2] - d.ZZ;
#endif
    }
    
finish:
    delete [] sca;
    delete [] dif;
    delete [] vec;
    delete [] val;
    delete [] dia;
    delete [] low;
    delete [] upe;
    
    return info;
}
#endif

/**
 The response of this method to a sudden perpendicular force is not ideal:
 For example, a force applied to the bottom of a vertical fibers leads
 to a 'L' configuration after one step of `solve()`.
 reshape() reduces the bottom leg of the 'L', by translating the entire vertical portion
 of the fiber, irrespective of the length of this section.
 */

#if ( 1 )   // 1 = optimized version of Filament::reshape_sure()

/**
 Move the model-points relative to each other, such that when this is done,
 all segments have the same distance `fnCut` ( =segmentation() ).
 This is operation does not change the center of gravity of the fiber.

 
 NOTE: if two consecutive points overlap, there is no unique way to
 restore the constraints! We do nothing in that case, because most 
 likely, the Brownian motion will push the points appart soon.
 */

void Filament::reshape_sure(const unsigned ns, real* vec, real cut)
{
    Vector dp(0,0,0), sum(0,0,0);
    Vector seg = diffPoints(vec, 0);
    real   dis = seg.norm();
    
    // translation needed to restore first segment
    if ( dis > REAL_EPSILON )
        dp = ( cut/dis - 1.0 ) * seg;
    
    for ( unsigned pp = 1; pp < ns; ++pp )
    {
        seg = diffPoints(vec, pp);
        dis = seg.norm();
        
        //move the left point by dp:
        dp.add_to(vec+DIM*pp);
        //update the uniform motion of the points:
        sum += dp;
        
        //add to the translation needed to restore this segment
        if ( dis > REAL_EPSILON )
            dp += ( cut/dis - 1.0 ) * seg;
    }
    
    //move the last point by dy[]:
    dp.add_to(vec+DIM*ns);
    
    // calculte a uniform motion to conserve the center of gravity:
    sum = ( sum + dp ) * ( -1.0 / ( ns + 1 ) );
    
    //translate the entire fiber uniformly:
    for ( unsigned pp = 0; pp <= ns; ++pp )
        sum.add_to(vec+DIM*pp);
}

#else

// ------------  old ( unoptimal ) version:
/**
 Move the model-points relative to each other, such that when this is done,
 all segments have the same distance segmentation() = fnCut.
 This is operation does not change the center of gravity of the fiber.
 */

void Filament::reshape_sure(const unsigned ns, real* vec, real cut)
{
    Vector dp, sum(0,0,0);
    
    for ( unsigned pp = 1; pp <= ns; ++pp )
    {
        dp       = diffPoints(vec, pp-1);
        real dis = dp.norm();
        if ( dis > REAL_EPSILON )
        {
            dp  *= ( cut/dis - 1.0 );
            for ( unsigned qq = pp; qq <= ns; ++qq )
                dp.add_to(vec+DIM*qq);
            sum += ( 1 + ns - pp ) * dp;
        }
    }
    
    sum *= ( -1.0 / (1+ns) );
    for ( unsigned pp = 0; pp <= ns; ++pp )
        sum.add_to(vec+DIM*pp);
}

#endif


void Filament::reshape()
{
    assert_true( nbPoints() > 1 );
#if ( DIM > 1 )
    if ( nbPoints() == 2 )
        reshape_two(psPos, psPos, fnCut);
    else if ( reshape_it(nbSegments(), psPos, psPos, fnCut) )
#endif
        reshape_sure(nbSegments(), psPos, fnCut);
}


void Filament::getPoints(const real * x)
{
#if ( DIM == 1 )
    PointSet::getPoints(x);
    reshape_sure(nbSegments(), psPos, fnCut);
#else
    if ( nbPoints() == 2 )
        reshape_two(x, psPos, fnCut);
    else if ( reshape_it(nbSegments(), x, psPos, fnCut) )
    {
        PointSet::getPoints(x);
        reshape_sure(nbSegments(), psPos, fnCut);
        //std::cerr << "Note that a crude method was used to reshape " << reference() << std::endl;
    }
#endif
    //dump(std::cerr);
}


/**
 Flip all the points. We do not change fnAscissa,
 and the abscissa of center thus stays as it is:
*/
void Filament::flip()
{
    unsigned ii = 0;
    unsigned jj = lastPoint();
    
    while ( ii < jj )
    {
        Vector P(psPos+DIM*ii);
        Vector Q(psPos+DIM*jj);
        Q.put(psPos+DIM*ii);
        P.put(psPos+DIM*jj);
        ++ii;
        --jj;
    }
}


//========================================================================
//=====================GROWING/SHRINKING==================================
//========================================================================
#pragma mark -

/**
 The argument 'dlen' can be positive or negative:
 - dlen > 0 : growth,
 - dlen < 0 : shrinkage
 .
 
 Note: This works nicely only if `dlen` is small compared to segmentation().
 For large decrease in length, use cutM().
 */
void Filament::growM(const real dlen)
{
    assert_true( length() + dlen > 0 );
    real a = -dlen / length();
    
    if ( dlen > 0 )
    {
        unsigned p = 0, n = nbSegments();
        Vector dp0 = diffPoints(0), dp1;
        movePoint(p, ( a * n ) * dp0);
        ++p;
        --n;
        
        if ( n > 0  &&  ( n & 1 ) )
        {
            dp1 = diffPoints(p);
            movePoint(p, ( a * n ) * dp0);
            dp0 = dp1;
            ++p;
            --n;
        }
        
        while ( n > 1 )
        {
            //assert_true( 0 == (p & 1) );
            dp1 = diffPoints(p);
            movePoint(p, ( a * n ) * dp0);
            ++p; --n;
            //assert_true( 1 == (p & 1) );
            dp0 = diffPoints(p);
            movePoint(p, ( a * n ) * dp1);
            ++p; --n;
        }
    }
    else if ( dlen < 0 )
    {
        for ( unsigned p = 0, n = nbSegments(); n > 0; ++p, --n )
            movePoint(p, ( a * n ) * diffPoints(p));
    }
    
    fnCut += dlen / real( nbSegments() );
    fnAbscissa -= dlen;
    postUpdate();
}

/**
 This extends the fiber by adding one segment at the MINUS_END.
 Thus `segmentation()` is not changed, and the existing points are not displaced.
 */
void Filament::addSegmentM()
{
    unsigned int pp = 1+nbPoints();
    setNbPoints(pp);
    
    pp *= DIM;
    while ( --pp >= DIM )
        psPos[pp] = psPos[pp-DIM];
    
    for ( pp = 0; pp < DIM; ++pp )
        psPos[pp] += psPos[pp] - psPos[pp+2*DIM];
    
    fnAbscissa -= fnCut;
    postUpdate();
}


/**
 The Fiber length is reduced by `dlen` ( which must be >= 0 ).
 The portion of size `dlen` near the MINUS_END is removed,
 the (fewer) model-points are recalculated.
 
 Note: after cutM(), the distance between the points is not exactly
 equal to segmentation(). This is true only if the fiber is straight.
 */
void Filament::cutM(const real dlen)
{
    real len = length();
    assert_true( 0 <= dlen );
    assert_true( dlen < len );
    
    const unsigned nbp = bestNumberOfPoints((len-dlen)/fnSegmentation);
    const real cut = (len-dlen) / (nbp-1);
    real* tmp = new real[DIM*nbp];

    // calculate intermediate points into tmp[]:
    for ( unsigned pp=0; pp+1 < nbp; ++pp )
    {
        Vector w = interpolateM(dlen+pp*cut).pos();
        w.put(tmp+DIM*pp);
    }
    
    // copy the position of plus-end into tmp[]:
    const unsigned lp = lastPoint();
    for ( unsigned d = 0 ; d < DIM; ++d )
        tmp[DIM*(nbp-1)+d] = psPos[DIM*lp+d];
    
    setNbPoints(nbp);
    
    // copy calculated points to psPos[]
    for ( unsigned pp = 0; pp < DIM*nbp; ++pp )
        psPos[pp] = tmp[pp];
    
    delete [] tmp;
    fnAbscissa += dlen;
    fnCut = cut;
    postUpdate();
}


/**
 The argument 'dlen' can be positive or negative:
 - dlen > 0 : growth,
 - dlen < 0 : shrinkage
 .
 
 Note: This works nicely only if `dlen` is small compared to segmentation().
 For large decrease in length, use cutP().
 */
void Filament::growP(const real dlen)
{
    assert_true( length() + dlen > 0 );
    real a = dlen / length();
    
    if ( dlen > 0 )
    {
        unsigned p = lastPoint();
        Vector dp0 = diffPoints(p-1), dp1;
        movePoint(p, ( a * p ) * dp0);
        --p;
        
        if ( p > 0  &&  ( p & 1 ) )
        {
            dp1 = diffPoints(p-1);
            movePoint(p, ( a * p ) * dp0);
            dp0 = dp1;
            --p;
        }
        
        while ( p > 1 )
        {
            //assert_true( 0 == (p & 1) );
            dp1 = diffPoints(p-1);
            movePoint(p, ( a * p ) * dp0);
            --p;
            //assert_true( 1 == (p & 1) );
            dp0 = diffPoints(p-1);
            movePoint(p, ( a * p ) * dp1);
            --p;
        }
    }
    else if ( dlen < 0 )
    {
        for ( unsigned p = lastPoint() ; p > 0 ; --p )
            movePoint(p, ( a * p ) * diffPoints(p-1));
    }
    
    fnCut += dlen / real( nbSegments() );
    postUpdate();
}


/**
 This extends the fiber by adding one segment at the PLUS_END.
 Thus `segmentation()` is not changed, and the existing points are not displaced.
 */
void Filament::addSegmentP()
{
    unsigned pp = nbPoints();
    setNbPoints(pp+1);
    
    real * psp = psPos + pp * DIM;
    for ( unsigned int dd = 0; dd < DIM; ++dd )
        psp[dd] = 2 * psp[dd-DIM] - psp[dd-2*DIM];
    
    postUpdate();
}


/**
 The Fiber length is reduced by `dlen` ( which must be >= 0 ).
 The portion of size `dlen` near the PLUS_END is removed,
 and the fewer model-points are recalculated.

 Note: after cutP(), the distance between the points is not exactly
 equal to segmentation(). This is true only if the fiber is straight.
*/
void Filament::cutP(const real dlen)
{
    real len = length();
    assert_true( 0 <= dlen );
    assert_true( dlen < len );
    
    const unsigned nbp = bestNumberOfPoints((len-dlen)/fnSegmentation);
    const real cut = (len-dlen) / (nbp-1);
    real* tmp = new real[DIM*nbp];
    
    // calculate intermediate points into tmp[]:
    for ( unsigned pp = 1; pp < nbp; ++pp )
    {
        Vector w = interpolateM(pp*cut).pos();
        w.put(tmp+DIM*pp);
    }
    
    setNbPoints(nbp);
    
    // copy calculated points to psPos[]
    // point at minus-end has not changed
    for ( unsigned pp = DIM; pp < DIM*nbp; ++pp )
        psPos[pp] = tmp[pp];
    
    delete [] tmp;
    fnCut = cut;
    postUpdate();
}

//------------------------------------------------------------------------------

void Filament::grow(FiberEnd end, const real dlen)
{
    if ( end == PLUS_END )
        growP(dlen);
    else if ( end == MINUS_END )
        growM(dlen);
}


void Filament::adjustLength(real len, FiberEnd ref)
{
    assert_true( len > 0 );
    
    if ( ref == PLUS_END )
    {
        if ( len < length() )
            cutP(length()-len);
        else
            growP(len-length());
    }
    else if ( ref == MINUS_END )
    {
        if ( len < length() )
            cutM(length()-len);
        else
            growM(len-length());
    }
}


void Filament::truncateM(const unsigned int p)
{
    PointSet::truncateM(p);
    fnAbscissa = abscissaPoint(p);
    postUpdate();
}


void Filament::truncateP(const unsigned int p)
{
    PointSet::truncateP(p);
    postUpdate();
}


/**
 `fib` is attached at the PLUS_END of `*this`
 
 The model-point are reinterpolated linearly, and the length of the
 segments will not fullfil the constraints of segmentation.
 If this is a problem, Filament::reshape() should be called.
 
 `fib` should usually be destroyed afterward.
 */
void Filament::join(Filament const* fib)
{
    const real len1 = length();
    const real lenT = len1 + fib->length();
    const unsigned nbr = bestNumberOfPoints(lenT/fnSegmentation) - 1;
    const real cut = lenT / real(nbr);
    
    // save position of PLUS_END:
    Vector ppe = fib->posEndP();
    
    real* tmp = new real[DIM*nbr];
    
    // calculate new points into tmp[]:
    for ( unsigned pp = 1; pp < nbr; ++pp )
    {
        Vector w;
        if ( pp*cut < len1 )
            w = interpolateM(pp*cut).pos();
        else
            w = fib->interpolateM(pp*cut-len1).pos();
        
        w.put(tmp+DIM*pp);
    }
    
    setNbPoints(nbr+1);
    
    // copy point back in place:
    for ( unsigned int pp = DIM; pp < DIM*nbr; ++pp )
        psPos[pp] = tmp[pp];
    
    ppe.put(psPos+DIM*nbr);
    
    delete [] tmp;
    fnCut = cut;
    postUpdate();
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 Returns the minimum and maximum distance between consecutive points
 */
void Filament::minMaxSegments(real& mn, real& mx) const
{
    real r = diffPoints(0).norm();
    mn = r;
    mx = r;
    for ( unsigned n = 1; n < lastPoint(); ++n )
    {
        real r = diffPoints(n).norm();
        if ( r > mx )
            mx = r;
        if ( r < mn )
            mn = r;
    }
}

/**
 Returns the average and variances of segment length
 */
void Filament::infoSegments(real& avg, real& var) const
{
    avg = 0;
    var = 0;
    unsigned cnt = nbSegments();
    for ( unsigned n = 0; n < cnt; ++n )
    {
        real r = diffPoints(n).norm();
        avg += r;
        var += r*r;
    }
    var = var - avg * avg / cnt;
    avg /= cnt;
}

/**
 curvature is the inverse of the radius of the circle containing three 
 consecutive model points A, B, C. It is calculated from the relations:
 @code
 cos(angle) = scalar product AB*BC
 sin(angle) = sqrt( 1 - cos(angle)^2 )
 Diameter-of-curvature = AC / sin(angle)
 Curvature = 2.0/Diameter
 @endcode
 */
real Filament::curvature(unsigned p) const
{
    assert_true( 0 < p && p < lastPoint() );
    Vector ab = diffPoints(p-1);
    Vector bc = diffPoints(p);
    real cs = ab * bc;
    real si = sqrt( 1 - ( cs * cs ) / ( ab.normSqr() * bc.normSqr() ) );
    return 2 * si / (ab+bc).norm();
}


/**
 The bending energy is an integral:
 1/2 * rigidity * sum( curvature^2 ds ),
 where s is curvilinear abscissa
 
 The normalized bending energy is an integral:
 1/2 * sum( curvature^2 ds )
 
 if theta = angle between two consecutive segment,
 Curvature = 1/R = 2 * sin(angle/2) / segmentation
 and
 2 * sin(angle/2) = sqrt(2-2*cos(angle))
 hence:
 0.5 * Curvature^2 = ( 1 - cos(angle) ) / ( segmentation * segmentation )
 and finaly:
 sum[ 0.5 * Curvature^2 * segmentation ] = sum[ 1 - cos(angle) ] / segmentation
 */
real Filament::bendingEnergy0() const
{
    real e = 0;
    real s = 1 / ( fnCut * fnCut );
    
    Vector dp0 = diffPoints(0);
    for ( unsigned p = 1; p < lastPoint() ; ++p )
    {
        Vector dp1 = diffPoints(p);
        e += 1 - s * ( dp0 * dp1 );
        dp0 = dp1;
        //std::clog << "curvature " << p << " " << curvature(p) << std::endl;
    }
    
    // we correct the result, because we only considered (nbPoints()-2) junctions,
    // and thus only a fraction of the total length of the fiber
    
    if ( nbPoints() > 2 )
        e *= ( nbPoints() - 1 ) / ( ( nbPoints() - 2 ) * segmentation() );
    
    return e;
}


real Filament::minCosinus() const
{
    real result;
    Vector dir1, dir2;
    
    unsigned ps = nbSegments() % 2;
    if ( ps )
    {
        dir1   = diffPoints(0);
        result = fnCut * fnCut;
    }
    else
    {
        dir1   = diffPoints(1);
        result = diffPoints(0) * dir1;
        ps = 2;
    }
    
    for ( ; ps < nbSegments(); ps+=2 )
    {
        dir2 = diffPoints(ps);
        real s = dir1 * dir2;
        if ( s < result ) result = s;
        dir1 = diffPoints(ps+1);
        real t = dir1 * dir2;
        if ( t < result ) result = t;
    }
    
    return result / ( fnCut * fnCut );
}


/**
 Returns the minimum and maximum distance between consecutive points
 */
unsigned Filament::nbKinks(real threshold) const
{
    threshold *= fnCut * fnCut;
    unsigned res = 0;
    Vector d = diffPoints(0);
    
    for ( unsigned n = 1; n < lastPoint(); ++n )
    {
        Vector r = diffPoints(n);
        if ( d * r < threshold )
            ++res;
        d = r;
    }
    return res;
}

/**
 This calculates the intersection between the support line of segment `s`,
 and the plane defined by <em> n.pos + a = 0 </em>

 @return scalar `x` specifying the intersection with the support line:
 - `x = 0` if intersection occurs at point 's'
 - `x in ]0, 1[` for intersections that are within the segment boundaries
 - `x = 1` if intersection occurs at point 's+1'
 - `x = INFINITY` if the segment is parallel to the plane
 .
 
 The abscissa of the intersection is `abscissaPoint(s+a)`.
 The position of the cut is `interpolatePoints(s, s+1, a)`
 */

real Filament::planarIntersect(unsigned s, Vector const& n, const real a) const
{
    assert_true( s < nbSegments() );
    
    real sca = diffPoints(s) * n;
    
    // if segment is parallel to plane, there is no intersection:
    if ( -REAL_EPSILON < sca  &&  sca < REAL_EPSILON )
        return INFINITY;
    
    Vector pos = posP(s);
    
    if ( modulo )
        modulo->fold(pos);
    
    return - ( pos * n + a ) / sca;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Note: Unless the Fiber is straight, the segments will not be exactly of length `fnCut`
 after the reinterpolation, and calling reshape() may be necessary.
 
 @todo 2d-order interpolation in Filament::resegment()
 */
void Filament::resegment(unsigned nps)
{
    assert_true( nps > 1 );
    
    //now nps is the number of segment
    --nps;
    
    real cut = length() / real(nps);
    
    // calculate new intermediate points in tmp[]:
    real* tmp = new real[DIM*nps];
    Vector a = posP(0), b = posP(1);
    
    real h = 0;
    unsigned p = 1;
    
    for ( unsigned n = 1; n < nps; ++n )
    {
        h += cut;
        
        while ( h > fnCut )
        {
            h -= fnCut;
            a = b;
            ++p;
            assert_true(p<nbPoints());
            b.get(psPos+DIM*p);
        }
        
        Vector w = a + ( h / fnCut ) * ( b - a );
        w.put(tmp+DIM*n);
    }
    
    // save index of PLUS_END
    p = DIM*lastPoint();
    
    // resize array:
    setNbPoints(nps+1);

    // move coordinates of last point
    for ( unsigned d = 0; d < DIM; ++d )
        psPos[DIM*nps+d] = psPos[p+d];

    // copy calculated coordinates back into psPos
    for ( unsigned d = DIM; d < DIM*nps; ++d )
        psPos[d] = tmp[d];

    delete [] tmp;
    fnCut = cut;
    reshape();
}


#ifndef CURVATURE_DEPENDENT_SEGMENTATION

/**
 A fiber is segmented as a function of its length.
 The number of segments `nseg` is the one that minimizes:
 @code
 abs( length/nseg - fnSegmentation )
 @endcode
 
 It is such that:
 @code
 Filament::segmentation() < 4/3 * FiberProp::segmentation
 Filament::segmentation() > 2/3 * FiberProp::segmentation
 @endcode
 */
void Filament::adjustSegmentation()
{
    assert_true( fnSegmentation > REAL_EPSILON );
    
    unsigned best = bestNumberOfPoints(length()/fnSegmentation);
    
    if ( best != nbPoints() )
    {
#if ( 1 )
        resegment(best);
#else
        unsigned np = nbPoints();
        // copy current points in temporary array:
        real* tmp = new real[DIM*np];
        for ( int n = 0; n < DIM*np; ++n )
            tmp[n] = psPos[n];
        // re-interpolate:
        setShape(tmp, np, best);
        delete[] tmp;
#endif
    }
}


#else

static const real RECUT_PRECISION = 0.05;

/**
 A fiber is segmented as a function of its curvature.
 Only one new segmentation is done every time step at maximum,
 to allow the solver to equilibrate the new model-points.
 */
void Filament::adjustSegmentation()
{
    PRINT_ONCE("adjustSegmentation() is using the fiber curvature\n");
    
    const int upLimit = 8;
    real rl = length();
    
    assert_true( fnSegmentation > REAL_EPSILON );
    //one segment for very short tubes
    if ( rl <= fnSegmentation )
    {
        if ( nbPoints() > 2 )
        {
            resegment(2);
            fnCutErrorIndex = 0;
        }
        return;
    }
    //at least two segments if length exceeds 2*fnSegmentation
    if ( nbPoints() < 3 )
    {
        if ( rl > fnSegmentation )
        {
            resegment(3);
            fnCutErrorIndex = 0;
        }
        return;
    }
    //the segment-length should not exceed 4*fnSegmentation
    if ( fnCut >= upLimit*fnSegmentation )
    {
        resegment(2*nbPoints()-1);
        fnCutErrorIndex = 0;
        return;
    }
    
    // accumulate the error in variable fnCutError, 
    // The error is in angle-square: 1-cos(angle) ~ (angle)^2
    if ( fnCutErrorIndex == 0 )
        fnCutError  =  1.0 - minCosinus();
    else
        fnCutError +=  1.0 - minCosinus();
    
    // after accumulation of RECUT_PERIOD time steps, we check the result
    if ( ++fnCutErrorIndex >= RECUT_PERIOD )
    {
        //we scale the error accumulated over the last steps
        fnCutError *= fnCut / ( RECUT_PRECISION * RECUT_PERIOD );
        
        //cut more finely if the error is large
        if ( fnCutError > 1.0 )
        {
            if ( fnCut > fnSegmentation )
                resegment(2*nbPoints()-1);
        }
        else if ( fnCutError < 0.1 )
        {
            //we want the error that we will have with fewer points to be small
            //at constant curvature, the error scales like (length-of-rods)^2
            //we expect the error to raise by 4x, if we divide the number of rods by 2.
            //to be safe we use 0.125, instead of 1/4, i.e. another factor 2.
            if (( nbPoints() > 3 ) && ( fnCut < (upLimit/2)*fnSegmentation ))
                resegment((nbPoints()+1)/2);
        }
        
        //reset the counter and error accumulator
        fnCutErrorIndex = 0;
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark -

/**
 return the abscissa with respect to the ORIGIN.
 */
real Filament::abscissaEnd(const FiberEnd end) const
{
    switch( end )
    {
        case ORIGIN:    return 0;
        case MINUS_END: return fnAbscissa;
        case PLUS_END:  return fnAbscissa + fnCut * nbSegments();
        case CENTER:    return fnAbscissa + 0.5 * fnCut * nbSegments();
        default:        ABORT_NOW("invalid argument value"); return 0;
    }
}


/**
 convert the abscissa that is specified from the given reference,
 to the standard abscissa measured from the ORIGIN.
 
 ATTENTION: the direction is inverted if `ref = PLUS_END`
 */
real Filament::abscissaFrom(const real ab, const FiberEnd ref) const
{
    switch( ref )
    {
        case ORIGIN:     return ab;
        case MINUS_END:  return ab + fnAbscissa;
        case CENTER:     return ab + fnAbscissa + 0.5 * fnCut * nbSegments();
        case PLUS_END:   return fnAbscissa + fnCut * nbSegments() - ab;
        default:         ABORT_NOW("invalid argument value"); return 0;
    }
}


/**
 The Fiber is partitionned by this function in three regions:
 - a MINUS_END part, which is of length 'lambda'
 - a PLUS_END part, also of length 'lambda'
 - and a NO_END section in between
 .
 A Fiber shorter than 2*len does not have a central region,
 and is composed of PLUS_END and MINUS_END parts of equal size.
 */    
FiberEnd Filament::whichEndDomain(const real ab, const real lambda) const
{
    const real abs = ab - fnAbscissa;
    const real len = length();
    
    if ( 2 * abs > len )
    {
        if ( abs >= len - lambda )
            return PLUS_END;
    }
    else
    {
        if ( abs <= lambda )
            return MINUS_END;
    }
    return NO_END;
}


//------------------------------------------------------------------------------
#pragma mark -


PointExact Filament::exactEnd(const FiberEnd which) const
{
    if ( which == MINUS_END )
        return PointExact(this, 0);
    else
    {
        assert_true( which == PLUS_END );
        return PointExact(this, lastPoint());
    }
}


PointInterpolated Filament::interpolateEnd(const FiberEnd which) const
{
    if ( which == MINUS_END )
        return interpolateEndM();
    else
    {
        assert_true( which == PLUS_END );
        return interpolateEndP();
    }
}


PointInterpolated Filament::interpolateCenter() const
{
    unsigned int n = lastPoint() / 2;
    if ( 2*n == lastPoint() )
        return PointInterpolated(this, n, n+1, 0);
    else
        return PointInterpolated(this, n, n+1, 0.5);
}


/**
 Same as interpolate(), but the abscissa `ab' is taken from the MINUS_END of the Fiber.
 */
PointInterpolated Filament::interpolateM(const real ab) const
{
    if ( ab <= 0 )
        return PointInterpolated(this, 0, 1, 0.0);
    
    double n, co = modf( ab / fnCut, &n );
    unsigned rd = (unsigned)n;
    unsigned rx = rd + 1;
    
    //beyond the last point, we interpolate the PLUS_END
    if ( rx < nbPoints() )
        return PointInterpolated(this, rd, rx, co);
    else
        return PointInterpolated(this, nbPoints()-2, nbPoints()-1, 1.0);
}


/**
 Convert abscissa `ab' into a PointInterpolated = ( a model-point `r' + a coefficient `a' ).
 The corresponding point X = P(r) * (1-a) + P(r+1) * a:
 - `r' is an integer: 0 <= r < lastPoint(),
 - `a' is a positive real coefficient: 0 <= a <= 1
 .
 In this function, the abscissa `ab` is taken from the ORIGIN of the Fiber.
 */
PointInterpolated Filament::interpolate(const real abo) const
{
    real ab = abo - fnAbscissa;
    
    if ( ab <= 0 )
        return PointInterpolated(this, 0, 1, 0.0);
    
    double n, co = modf( ab / fnCut, &n );
    unsigned rd = (unsigned)n;
    unsigned rx = rd + 1;
    
    //beyond the last point, we interpolate the PLUS_END
    if ( rx < nbPoints() )
        return PointInterpolated(this, rd, rx, co);
    else
        return PointInterpolated(this, nbPoints()-2, nbPoints()-1, 1.0);
}


PointInterpolated Filament::interpolate(const real ab, const FiberEnd from) const
{
    switch( from )
    {
        case ORIGIN:
            return interpolate(ab);
            
        case MINUS_END:
            return interpolateM(ab);
            
        case CENTER:
            return interpolateM(ab + 0.5*fnCut*nbSegments());
            
        case PLUS_END:  //this is counted from the plus towards the minus end
            return interpolateM(fnCut*nbSegments() - ab);
        
        default:
            ABORT_NOW("invalid argument value");
    }
    return interpolate(0);
}

//------------------------------------------------------------------------------
#pragma mark -

Vector Filament::posEnd(FiberEnd which) const
{
    if ( which == MINUS_END )
        return posP(0);
    else if ( which == PLUS_END )
        return posP(lastPoint());
    else
        return interpolate(0, which).pos();
}

#if ( DIM > 1 )
Vector Filament::posM(const real ab) const
{
    //return MINUS_END
    if ( ab <= 0 )
        return posP(0);
    
    double n, co = modf( ab / fnCut, &n );
    unsigned int rd = (unsigned int)n;
    
    //beyond the last point, we return the position of the PLUS_END
    if ( rd < lastPoint() )
        return interpolatePoints(rd, rd+1, co);
    else
        return posP(lastPoint());
}
#endif

Vector Filament::pos(const real ab, const FiberEnd from) const
{
    return interpolate(ab, from).pos();
}


//------------------------------------------------------------------------------
#pragma mark -

Vector Filament::dir(const real ab) const
{
    return dirPoint(interpolate(ab).point1());
}

Vector Filament::dir(const real ab, const FiberEnd from) const
{
    return dirPoint(interpolate(ab, from).point1());
}


Vector Filament::dirEnd(const FiberEnd which) const
{
    if ( which == MINUS_END )
        return dirPoint(0);
    else if ( which == PLUS_END )
        return dirPoint(lastSegment());
    else
        return dirPoint(interpolate(0, which).point1());
}


/**
 The returned value is negative when the force antagonizes elongation,
 and this is true at both ends. 
 */
real Filament::projectedForceEnd(const FiberEnd which) const
{
    if ( which == PLUS_END )
        return projectedForceEndP();
    else
    {
        assert_true( which == MINUS_END );
        return projectedForceEndM();
    }
}


real Filament::angleXY(const unsigned p) const
{
    Vector d = diffPoints(p);
#if ( DIM == 3 )
    //return atan2(sqrt(d.YY*d.YY+d.ZZ*d.ZZ), d.XX);
    return atan2(d.YY, d.XX);
#elif ( DIM == 2 )
    return atan2(d.YY, d.XX);
#else
    return d.XX > 0 ? 0 : M_PI;
#endif
}


//------------------------------------------------------------------------------
#pragma mark -


void Filament::checkLength() const
{
    real L = trueLength();
    if ( std::abs( L - length() ) > 0.1 )
        std::clog << "Warning: length of " << reference() << " is " << L << " but " << length() << " was expected\n";
}


/**
 Prints info on the length of Segments, which can be useful for debugging
 */
void Filament::dump(std::ostream& os) const
{
    os << "Fiber " << std::setw(7) << reference();
    os << "  " << std::left << std::setw(6) << fnCut << " {";
    
    if ( 1 )
    {
        real mn, mx;
        minMaxSegments(mn, mx);
        real p = 100 * ( mx - mn ) / mn;
        os.precision(4);
        os << " " << mn << " + " << std::setw(7) << std::fixed << p << " %";
    }
    else
    {
        for ( unsigned pp = 0; pp < lastPoint(); ++pp )
        {
            real p = 100 * ( diffPoints(pp).norm() / fnCut - 1.0 );
            os << " " << std::left << std::setw(5) << p;
        }
    }
    os << " }" << std::endl;
}


void Filament::write(Outputter& out) const
{
    out.writeUInt32(signature());
    out.writeFloat(length());
    out.writeFloat(fnSegmentation);
    out.writeFloat(fnAbscissa);
    PointSet::write(out);
}


/**
 The fiber will be re-segmented if its current desired segmentation 
 does not match the one stored in the file.
 */
void Filament::read(Inputter & in, Simul& sim, Tag tag)
{
#if ( 0 )
    fpos_t pos;
    in.get_pos(pos);
    std::clog << "  reading Filament at " << pos << std::endl;
#endif

    unsigned s = in.readUInt32();
    if ( s ) signature(s);
    
    real len   = in.readFloat();
    real seg   = in.readFloat();
    fnAbscissa = in.readFloat();
    
    if ( len <= 0 )
        throw InvalidIO("invalid (negative) fiber length");

    if ( len > 1e6 )
        throw InvalidIO("excessive fiber length");
    
    if ( seg <= 1e-6 || seg > 1e6 )
        throw InvalidIO("invalid fiber segmentation");

    PointSet::read(in, sim, tag);
    
    if ( nbPoints() < 2 )
        throw InvalidIO("invalid fiber with 0 or 1 point");

#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() <= 37 )
        fnCut = len;
    else
#endif
        fnCut = len / nbSegments();

    postUpdate();
    
    /// resegment if the sementation parameter has changed:
    if ( fnSegmentation != seg )
        adjustSegmentation();
    
    //PointSet::write(std::cerr);
    
    // verify the length and segmentation:
    if ( in.vectorSize() == DIM )
    {
        checkLength();
        real mn, mx;
        minMaxSegments(mn, mx);
        real p = ( mx - mn ) / mn;
        if ( p > 0.01 )
        {
            real p = ( mx - mn ) / mn;
            std::cerr << "Warning: Fiber has non-uniform segments in [ " << mn << ", " << mx << " ] dev = " << p << std::endl;
        }
    }
}


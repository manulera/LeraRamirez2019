// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


#include <cmath>
#include "real.h"
#include <cstdio>
#include "assert_macro.h"

/**
 Calculate the projection P = (pX, pY) of the point W = (wX, wY) on the ellipse that
 is aligned with the X and Y axis, and has radii (lenX, lenY).
 
 Method:
 
 The normal to the ellipse at position P is N = ( pX / lenX^2, pY / lenY^2 ),
 and we can thus write W = P + h * N, which leads to:
 @code
 pX = wX * lenX^2 / ( lenX^2 + h );
 pY = wY * lenY^2 / ( lenY^2 + h );
 @endcode
 if wX and wY are not both null.
 
 Moreover, the projection should be on the ellipse and thus `h` should be a zero of:
 @code
 F(h) = ( pX / lenX )^2 + ( pY / lenY )^2 - 1
 @endcode
 We follow Newton's rule to find the root of F(h), and use the formula above to
 calculate the projection.
 */
void projectEllipse(real&   pX, real&  pY,
                    real    wX, real   wY,
                    real  lenX, real lenY)
{
    assert_true( lenX > REAL_EPSILON );
    assert_true( lenY > REAL_EPSILON );
    
    // handle the pathological cases:
    if ( wX == 0 )
    {
        pX = 0;
        pY = ( wY > 0 ) ? lenY : -lenY;
        return;
    }
    if ( wY == 0 )
    {
        pX = ( wX > 0 ) ? lenX : -lenX;
        pY = 0;
        return;
    }
    
    real aa = lenX * lenX;
    real bb = lenY * lenY;
    
    real h = 0, h_old;
    
    // we derive a lower limit for 'h' from  pX^2 + pY^2 > max(lenX,lenY)^2
    real RR = ( bb < aa ) ? aa : bb;
    // 'hmin' is the minimum value that 'h' could have
    real hmin = sqrt( ( wX*wX*aa*aa + wY*wY*bb*bb ) / RR ) - RR;
    
    // we derive another lower limit for 'h' from  |pX| < lenX
    h = ( fabs(wX) - lenX ) * lenX;
    if ( h > hmin )
        hmin = h;

    // we derive another lower limit for 'h' from  |pY| < lenY
    h = ( fabs(wY) - lenY ) * lenY;
    if ( h > hmin )
        hmin = h;

    // if the point is outside, then 'h' should be positive:
    if ( wX*wX/aa + wY*wY/bb > 1 )
    {
        if ( hmin < 0 )
            hmin = 0;
    }
    
    h = hmin;

    //fprintf(stderr, " <<< %+.10f  %+.10f    hmin %+10.4f", wX, wY, hmin);

    // follow Newton's iteration to find the root
    unsigned cnt = 0;
    do {
        real aah = aa + h;
        real bbh = bb + h;
#if ( 0 )
        pX = wX * aa / aah;
        pY = wY * bb / bbh;
         
        real pXX = pX * pX / aa;
        real pYY = pY * pY / bb;
#else
        real waX = wX / aah;
        real waY = wY / bbh;
        
        real pXX = waX * waX * aa;
        real pYY = waY * waY * bb;
#endif
        h_old = h;
        
        real F    = 1 - ( pXX         + pYY       );
        real dF   = 2 * ( pXX / aah   + pYY / bbh );
#if ( 0 )
        real ddF  =   - pXX / ( aah * aah ) - pYY / ( bbh * bbh );  // * 2
        //dddF =    + pXX / ( aah * aah * aah ) + pYY / ( bbh * bbh * bbh );  // * 6
        //fprintf(stderr, "       %+.10f   %+.10f   %+.10f\n", F, dF, ddF);

        // Halley's method convergence is cubic in general
        h -= ( F * dF ) / ( dF * dF - F * ddF );
#else
        // Newtons' method
        h -= F / dF;
#endif
        
        //fprintf(stderr, "  %i : h %+f  F %+20.16f  dF %+20.16f  dh %e\n", cnt, h, F, dF, h-h_old);
        
        if ( h < hmin )
        {
            h = 0.5 * ( h_old + hmin );
            continue;
        }
        
#if ( 0 )
        if ( cnt > 16 )
            fprintf(stderr, "projectEllipse fails %u :  h %+f  F %+e  dh %e\n", cnt, h, F, h-h_old);
#endif

        if ( ++cnt > 20 )
            break;
        
    } while ( h > h_old );

    // calculate the projection from h
    pX = wX * aa / ( aa + h );
    pY = wY * bb / ( bb + h );
    
#if ( 0 )
    // verify that projection is on ellipse:
    real F = 1 - ( pX*pX/aa + pY*pY/bb );
    fprintf(stderr, " %2i  >>> h %12.8f  F  %+e\n", cnt, h, F);
#endif
}






/**
 Calculates the projection P = (pX, pY, pZ) of the point W = (wX, wY, wZ) on the ellipse that
 is aligned with the X and Y axis, and has radii (lenX, lenY, lenZ).
 
 Method:
 
 The normal to the ellipse at position P is N = ( pX / lenX^2, pY / lenY^2, pZ / lenZ^2 ),
 and we can thus write W = P + h * N, for some scalar `h' which leads to:
 @code
 pX = wX / ( 1 + h / lenX^2 );
 pY = wY / ( 1 + h / lenY^2 );
 pZ = wZ / ( 1 + h / lenZ^2 );
 @endcode
 
 Moreover, the projection should be on the ellipse and thus `h` should be a zero of:
 @code
 F(h) = ( pX / lenX )^2 + ( pY / lenY )^2 + ( pZ / lenZ )^2 - 1
 @endcode
 We follow Newton's rule to find the root of F(h), and use the formula above to
 calculate the projection.
 */
void projectEllipsoid(real  p[3],
                      const real w[3],
                      const real len[3])
{
    assert_true( len[0] > REAL_EPSILON );
    assert_true( len[1] > REAL_EPSILON );
    assert_true( len[2] > REAL_EPSILON );
    
    // handle the pathological cases:
    if ( w[0] == 0 )
    {
        p[0] = 0;
        projectEllipse(p[1], p[2], w[1], w[2], len[1], len[2]);
        return;
    }
    if ( w[1] == 0 )
    {
        p[1] = 0;
        projectEllipse(p[0], p[2], w[0], w[2], len[0], len[2]);
        return;
    }
    if ( w[2] == 0 )
    {
        p[2] = 0;
        projectEllipse(p[0], p[1], w[0], w[1], len[0], len[1]);
        return;
    }

    real aa = len[0] * len[0];
    real bb = len[1] * len[1];
    real cc = len[2] * len[2];
    
    real h = 0, h_old;

    // we derive a lower limit for 'h' from  pX^2 + pY^2 + pZ^2 < max(lenX,lenY,lenZ)^2
    real RR = ( bb < aa ) ? ( cc < aa ? aa : cc ) : ( cc < bb ? bb : cc );
    // 'hmin' is the minimum value that 'h' can have
    real hmin = sqrt( ( w[0]*w[0]*aa*aa + w[1]*w[1]*bb*bb + w[2]*w[2]*cc*cc ) / RR ) - RR;

    // we derive another lower limit for 'h' from  |pX| < lenX
    h = ( fabs(w[0]) - len[0] ) * len[0];
    if ( h > hmin )
        hmin = h;

    // we derive another lower limit for 'h' from  |pY| < lenY
    h = ( fabs(w[1]) - len[1] ) * len[1];
    if ( h > hmin )
        hmin = h;
    
    // we derive another lower limit for 'h' from  |pZ| < lenZ
    h = ( fabs(w[2]) - len[2] ) * len[2];
    if ( h > hmin )
        hmin = h;

    if ( w[0]*w[0]/aa + w[1]*w[1]/bb + w[2]*w[2]/cc > 1 )
    {
        // if the point is outside, then 'h' should be positive:
        if ( hmin < 0 )
            hmin = 0;
    }

    h = hmin;
    //fprintf(stderr, "----- h %+f\n", h);

    /*
     Follow Newton's iteration to find the largest root.
     We start with h>0, and h should only increase
     */
    unsigned cnt = 0;
    do {
        real aah = aa + h;
        real bbh = bb + h;
        real cch = cc + h;
#if ( 0 )
        real pX = w[0] * aa / aah;
        real pY = w[1] * bb / bbh;
        real pZ = w[2] * cc / cch;
        
        real pXX = pX * pX / aa;
        real pYY = pY * pY / bb;
        real pZZ = pZ * pZ / cc;
#else
        real waX = w[0] / aah;
        real waY = w[1] / bbh;
        real waZ = w[2] / cch;
        
        real pXX = waX * waX * aa;
        real pYY = waY * waY * bb;
        real pZZ = waZ * waZ * cc;
#endif
        h_old = h;

        real   F = 1 - ( pXX         + pYY         + pZZ       );
        real  dF = 2 * ( pXX / aah   + pYY / bbh   + pZZ / cch );
#if ( 0 )
        real ddF =   - pXX/(aah*aah) - pYY/(bbh*bbh) - pZZ/(cch*cch);  // * 2
        
        // Halley's method convergence is cubic in general
        h -= ( F * dF ) / ( dF * dF - F * ddF );
#else
        // Newton's method
        h -= F / dF;
#endif
        
        //fprintf(stderr, "  %i : h %+f  F %+e dh %+.20f\n", cnt, h_old, F, h-h_old);
        //fprintf(stderr, "       %+.10f   %+.10f   %+.10f   %+.10f\n", F, F/dF, ddF/dF, dddF/dF);

        if ( h < hmin )
        {
            h = 0.5 * ( h_old + hmin );
            continue;
        }

#if ( 0 )
        if ( cnt > 16 )
        {
            fprintf(stderr, "projectEllipsoid fails %u :  h %+f  F %.6e dh %.6e\n", cnt, h_old, F, h-h_old);
            //fprintf(stderr, "    pos  %+.10f     %+.10f       %+.10f\n", w[0], w[1], w[2]);
            //fprintf(stderr, "    F    %+.10f  dF %+.10f   ddF %+.10f\n", F, dF, ddF);
        }
#endif

        if ( ++cnt > 20 )
            break;
        
    } while ( h > h_old );

    // calculate the projection from h
    p[0] = w[0] * aa / ( aa + h );
    p[1] = w[1] * bb / ( bb + h );
    p[2] = w[2] * cc / ( cc + h );
    
#if ( 0 )
    // verify that projection is on ellipse
    real F = 1 - ( p[0]*p[0]/aa + p[1]*p[1]/bb + p[2]*p[2]/cc );
    fprintf(stderr, " %2i  >>> h %12.8f  F  %+e\n", cnt, h, F);
#endif
}





#if ( 0 )

// Halley's method convergence is cubic in general
dh = -( F * dF ) / ( dF * dF - 0.5 * F * ddF );

// 4th order method:
dh = -F * ( dF * dF - 0.5 * F * ddF ) / ( dF * dF * dF - F * dF * ddF + dddF * F * F / 6 );

// modified 4th order method:
dh = -F / dF * ( 1 + 0.5 * F * ddF / ( dF * dF ) + F * F * ( 3 * ddF * ddF - dF * dddF ) / ( 6 * dF * dF * dF * dF ));

// or equivalently:
real R = F / dF;
real T = ddF / dF;
dh = -R * ( 1 + 0.5 * R * T + R * R * ( 0.5 * T * T - dddF / ( 6 * dF ) ));

#endif


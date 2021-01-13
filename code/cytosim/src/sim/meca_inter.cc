// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "sim.h"
#include "space.h"
#include "simul.h"
#include "modulo.h"
#include "point_exact.h"
#include "point_interpolated.h"
//#include "vecprint.h"


extern Modulo const* modulo;

/// enable to address the block of MatrixSparseSymmetricBLock directly
#define USE_MATRIX_BLOCK

//------------------------------------------------------------------------------
#pragma mark - Display

/**
 Option to allow the user to display the Meca-links online in 'play',
 by pressing the function key 'F1'.  
 This requires two calls to Meca::setInteractions()
 This option is normally OFF.
 */
//#define DISPLAY_INTERACTIONS



#ifdef DISPLAY_INTERACTIONS

  #include "opengl.h"
  #include "gle.h"
  #include "gle_color_list.h"

/// line with to display interactions
const GLfloat lw = 3.0;

/// Display a 2-vector interaction
void displayInteraction(const Vector & a, const Vector & b)
{
    glLineWidth(lw);
    glBegin(GL_LINES);
    gle::gleVertex(a);
    gle::gleVertex(b);
    glEnd();
}

/// Display a 2-vector interaction
void displayInteraction(const Vector & a, Vector b, real len)
{
    if ( modulo ) modulo->fold(b, a);
    Vector dx = ( b - a );
    dx *= 0.5 * ( 1.0 - len / dx.norm() );
    glLineWidth(lw);
    glBegin(GL_LINES);
    gle::gleVertex(a);
    gle::gleVertex(a+dx);
    gle::gleVertex(b-dx);
    gle::gleVertex(b);
    glEnd();
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(1, 0x3333);
    glBegin(GL_LINES);
    gle::gleVertex(a+dx);
    gle::gleVertex(b-dx);
    glEnd();
    glDisable(GL_LINE_STIPPLE);
    glPointSize(10);
    glBegin(GL_POINTS);
    gle::gleVertex(a);
    gle::gleVertex(b);
    glEnd();
}

/// Display a 3-vector interaction
void displayInteraction(const Vector & a, Vector b, Vector c)
{
    if ( modulo )
    {
        modulo->fold(b, a);
        modulo->fold(c, a);
    }
    glLineWidth(lw);
    glLineStipple(1, 0x7310);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    gle::gleVertex(a);
    gle::gleVertex(b);
    glEnd();
    glDisable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    gle::gleVertex(b);
    gle::gleVertex(c);
    glEnd();
    glPointSize(10);
    glBegin(GL_POINTS);
    gle::gleVertex(b);
    glEnd();
}


/// Display a 4-vector interaction
void displayInteraction(const Vector & a, const Vector & b,
                        const Vector & c, const Vector & d)
{
    glLineWidth(lw);
    glLineStipple(1, 0x7171);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    gle::gleVertex(a);
    gle::gleVertex(b);
    gle::gleVertex(c);
    gle::gleVertex(d);
    glEnd();
    glDisable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    gle::gleVertex(b);
    gle::gleVertex(c);
    glEnd();
    glPointSize(10);
    glBegin(GL_POINTS);
    gle::gleVertex(b);
    gle::gleVertex(c);
    glEnd();
}

#endif

//------------------------------------------------------------------------------
#pragma mark - any_equal

/// true if two values are equal
bool any_equal(const Meca::index_type a, const Meca::index_type b)
{
    if ( a == b ) return true;
    return false;
}


/// true if any two values are equal
bool any_equal(const Meca::index_type a, const Meca::index_type b,
               const Meca::index_type c)
{
    if ( a == b ) return true;
    if ( a == c ) return true;
    if ( b == c ) return true;
    return false;
}


/// true if any two values are equal
bool any_equal(const Meca::index_type a, const Meca::index_type b,
               const Meca::index_type c, const Meca::index_type d)
{
    if ( a == b ) return true;
    if ( a == c ) return true;
    if ( a == d ) return true;
    if ( b == c ) return true;
    if ( b == d ) return true;
    if ( c == d ) return true;
    return false;
}


//------------------------------------------------------------------------------
#pragma mark - Force
//------------------------------------------------------------------------------

/**
 Add constant force at index 'inx'
 */
void Meca::addForce(const index_type inx, const Vector & force)
{
    assert_true( inx < DIM*nbPts );
    force.add_to(vBAS+inx);
}


/**
 Add constant force to `pte`
 */
void Meca::addForce(const PointExact & pte, const Vector & force)
{
    const index_type inx = DIM * pte.matIndex();
    force.add_to(vBAS+inx);
}


/**
Add constant force to `pti`
 */
void Meca::addForce(const PointInterpolated & pti, const Vector & force)
{
    const index_type inx1 = DIM * pti.matIndex1();
    const index_type inx2 = DIM * pti.matIndex2();
    
    force.add_to(pti.coef2(), vBAS+inx1);
    force.add_to(pti.coef1(), vBAS+inx2);
}


//------------------------------------------------------------------------------
#pragma mark - Torque
//------------------------------------------------------------------------------

/**
 Add constant torque:
 
 @code
 force = torque ^ position
 @endcode

 This is explicit and all contributions go in the force vector vBAS[]
*/
void Meca::addTorque(const PointInterpolated & pti, const Torque & torque)
{
    const index_type inx1 = DIM * pti.matIndex1();
    const index_type inx2 = DIM * pti.matIndex2();
    
    Vector d = pti.diff();
    Vector f = cross(torque/d.normSqr(), d);
    
    f.sub_to(vBAS+inx1);
    f.add_to(vBAS+inx2);
}


/**
 Add an explicit torque to constrain a segment in direction `dir`,
 with a given weight:
 
 @code
 torque = weigth * ( segment.normalized() ^ dir )
 force = torque ^ position
 @endcode

 This assumes norm(dir) == 1
 This is explicit and all contributions go in the force vector vBAS[]
*/
void Meca::addTorqueClamp(const PointInterpolated & pti,
                          const Vector & dir,
                          const real weight)
{
    assert_true( weight >= 0 );
    const index_type inx1 = DIM * pti.matIndex1();
    const index_type inx2 = DIM * pti.matIndex2();
    
    Vector d = pti.diff();
    real n = d.normSqr();

    Torque Tq = cross(d, dir);

#if ( DIM == 3 )
    
    const real Tn = Tq.norm();

#else
    
    const real Tn = std::abs(Tq);

#endif
    
    // calculate current angle between segments in [0, pi]:
    const real angle = atan2(Tn, d*dir);
    
    /**
     Scale torque to make it proportional to angle:
     we multiply the vector Tq by angle / sin(angle),
     knowing that Tn = Tq.norm() = sin(angle) * sqrt(n)
     
     To have a Torque proportional to sin(angle), use:
     real nn = weight / ( n * sqrt(n) );
     */
    real nn = weight * angle / ( n * Tn );
    
    Vector f = cross(Tq * nn, d);
    
    f.sub_to(vBAS+inx1);
    f.add_to(vBAS+inx2);
}


/**
 Add an explicit torque to bring two segments parallel to each other,
 with a given weight:
 
 @code
 torque = weigth * ( dirA ^ dirB )
 forceA =  torque ^ dirA
 forceB = -torque ^ dirB
 @endcode

 This is explicit and all contributions go in the force vector vBAS[]
 */
void Meca::interTorque(const PointInterpolated & pta,
                       const PointInterpolated & ptb,
                       const real weight)
{
    assert_true( weight >= 0 );

    const index_type ia1 = DIM * pta.matIndex1();
    const index_type ia2 = DIM * pta.matIndex2();

    const index_type ib1 = DIM * ptb.matIndex1();
    const index_type ib2 = DIM * ptb.matIndex2();

    if ( any_equal(ia1, ia2, ib1, ib2) )
        return;
    
    Vector da = pta.diff();
    Vector db = ptb.diff();

    real na = da.normSqr();
    real nb = db.normSqr();

    Torque Tq = cross(da, db);
    
#if ( DIM == 3 )
    
    const real Tn = Tq.norm();
    
#else
    
    const real Tn = std::abs(Tq);
    
#endif
    
    // calculate current angle between segments in [0, pi]:
    const real angle = atan2(Tn, da*db);
    
    /**
     Scale torque to make it proportional to angle:
     we multiply the vector Tq by angle / sin(angle),
     knowing that Tn = Tq.norm() = sin(angle) * sqrt( na * nb )
     
     To have a Torque proportional to sin(angle), use:
     real nn = sqrt( na * nb ); or nn = Tn / angle
     na = weight / ( na * nn );
     nb = weight / ( nb * nn );
     */
    na = weight * angle / ( na * Tn );
    nb = weight * angle / ( nb * Tn );
    
    Vector fa = cross(Tq * na, da);
    Vector fb = cross(db, Tq * nb);

    fa.sub_to(vBAS+ia1);
    fa.add_to(vBAS+ia2);
    fb.sub_to(vBAS+ib1);
    fb.add_to(vBAS+ib2);
}



/**
 Add an explicit torque to bring two segments parallel to each other,
 making an angle defined by (cosinus, sinus), and with a given weight:
 
 @code
 torque = weigth * ( dirA ^ dirB.rotated(angle) )
 forceA =  torque ^ dirA
 forceB = -torque ^ dirB
 @endcode
 
 The direction of `ptb` is rotated around `axis` defined as ( dirA ^ dirB ).
 
 This is explicit and all contributions go in the force vector vBAS[]
 */
void Meca::interTorque(const PointInterpolated & pta,
                       const PointInterpolated & ptb,
                       const real cosinus, const real sinus,
                       const real weight)
{
    assert_true( weight >= 0 );

    const index_type ia1 = DIM * pta.matIndex1();
    const index_type ia2 = DIM * pta.matIndex2();
    
    const index_type ib1 = DIM * ptb.matIndex1();
    const index_type ib2 = DIM * ptb.matIndex2();
    
    if ( any_equal(ia1, ia2, ib1, ib2) )
        return;
    
    Vector da = pta.diff();
    Vector db = ptb.diff();

    real na = da.normSqr();
    real nb = db.normSqr();
    
#if ( DIM == 3 )
    
    /*
     in 3D the axis of torque is perpendicular to both `da` and `db`,
     ad the angle is only defined between 0 and PI,
     */
    Vector axis = cross(db, da).normalized( sinus > 0 ? 1 : -1 );
    
    // rotate vector `db` around `arm` by angle specified as (cosinus, sinus):
    Vector rot = cosinus * db + sinus * cross(axis, db);
    
#elif ( DIM == 2 )

    // this correspond to the Z-direction, up or down:
    real dir = cross(da, db) > 0 ? 1 : -1;

    // rotate vector `db` by angle defined by (cosinus, sinus) around Z
    Vector rot( db.XX*cosinus + db.YY*sinus*dir, db.YY*cosinus - db.XX*sinus*dir );
    
#else
    
    // this is meaningless but makes compilation possible
    Vector rot(0, 0);
    
    throw InvalidParameter("Meca::interTorque is meaningless in 1D");

#endif

    // calculate torque by vector-product:
    Torque Tq = cross(da, rot);
    
#if ( DIM == 3 )
    
    const real Tn = Tq.norm();
    
#else
    
    const real Tn = std::abs(Tq);
    
#endif
    
    // calculate current angle between segments in [0, pi]:
    const real angle = atan2(Tn, da*rot);
    
    /**
     Scale torque to make it proportional to angle:
     we multiply the vector Tq by angle / sin(angle),
     but knowing that Tn = Tq.norm() = sin(angle) * sqrt( na * nb )
     
     To have a Torque proportional to sin(angle), use:
     real nn = sqrt( na * nb ); or nn = Tn / angle;
     na = weight / ( na * nn );
     nb = weight / ( nb * nn );
     */
    na = weight * angle / ( na * Tn );
    nb = weight * angle / ( nb * Tn );
 
    // forces are divided appropriately to match desired torque:
    Vector fa = cross(Tq * na, da);
    Vector fb = cross(db, Tq * nb);
    
    vBAS[ia1  ] -= fa.XX;
    vBAS[ia2  ] += fa.XX;
    vBAS[ib1  ] -= fb.XX;
    vBAS[ib2  ] += fb.XX;
#if ( DIM > 1 )
    vBAS[ia1+1] -= fa.YY;
    vBAS[ia2+1] += fa.YY;
    vBAS[ib1+1] -= fb.YY;
    vBAS[ib2+1] += fb.YY;
#endif
#if ( DIM > 2 )
    vBAS[ia1+2] -= fa.ZZ;
    vBAS[ia2+2] += fa.ZZ;
    vBAS[ib1+2] -= fb.ZZ;
    vBAS[ib2+2] += fb.ZZ;
#endif
}



#if (DIM == 2)
/**
 Update Meca to include torque between segment A-B and C-D containing pt1 and pt2.
 Implicit version with linearized force 2D
 Angle is between AB and CD. Force is along normal N_A and N_C pointing to the other filament
 L_AB and L_CD is the length of the segments AB and CD
 force_A = torque_weight * ( Delta angle ) * N_A/L_AB =-force_B
 force_C = torque_weight * ( Delta angle ) * N_C/L_CD =-force_D
 Delta_angle is the difference between actual angle and resting angle between AB and CD
 
 Antonio Politi, 2013
 */
void Meca::interTorque2D(const PointInterpolated & pt1,
                         const PointInterpolated & pt2,
                         const real cosinus, const real sinus,
                         const real weight)
{
    assert_true( weight >= 0 );
    if ( pt1.overlapping(pt2) )
        return;
    
    //index in the matrix mC:
    const index_type index[] = { DIM*pt1.matIndex1(),DIM*pt1.matIndex1()+1, DIM*pt1.matIndex2(),
        DIM*pt1.matIndex2()+1, DIM*pt2.matIndex1(), DIM*pt2.matIndex1()+1,  DIM*pt2.matIndex2(),
        DIM*pt2.matIndex2()+1 };
    
    //Vectors and points of torque
    Vector ab = pt1.diff();
    Vector cd = pt2.diff();
    Vector a = pt1.pos1();
    Vector b = pt1.pos2();
    Vector c = pt2.pos1();
    Vector d = pt2.pos2();
    const real coord[]={a.XX, a.YY, b.XX, b.YY, c.XX, c.YY, d.XX, d.YY};
    //Helping vector this vector is at torque_angle from cd.
    //Therefore in resting state angle difference between ab and ce is zero. This vector is used to compute the strength of torque
    Vector ce;
    ce.XX =  cd.XX*cosinus + cd.YY*sinus;
    ce.YY = -cd.XX*sinus   + cd.YY*cosinus;
    //normalize
    const real abn = ab.norm();
    const real abnS= ab.normSqr();
    const real cdn = cd.norm();
    const real cdnS= cd.normSqr();
    if (abn < REAL_EPSILON || cdn < REAL_EPSILON ) return;
    
    //normalize the vectors
    ab /= abn; cd /= cdn; ce /= cdn;
    
    //Coordinates of normal vectors yielding the direction of the force
    //fa = torque_weight*dangle*(h[0], h[1]) = torque_weight*dangle*na/la
    const real h[]={ ab.YY/abn, -ab.XX/abn, -ab.YY/abn, ab.XX/abn, -cd.YY/cdn, cd.XX/cdn, cd.YY/cdn, -cd.XX/cdn };
    
    //dangle = angle - torque_angle
    //real dangle = atan2( cross(ab, ce), ab * ce );
    real dangle = atan2( ab.XX*ce.YY - ab.YY*ce.XX, ab * ce );
    //Computation of the jacobian for the linearization
    //M = d_x f = M1 + M2
    //M1 = w1/l normal d_x dangle
    //M2 = w2 * dangle  d_x normal/l
    real w1 = weight;
    real w2 = weight*dangle;
    
    
    //Matrix M1 with k*hxh (outer product) this yieald a matrix stored with its lower triangular part in m. The -w1 is because ab = b-a
    real m[36] = { 0 };
    //blas_xspr('U', 8, -w1, h, 1, m);
    blas_xspr('L', 8, -w1, h, 1, m);
    
    
    
    //Matrix M2
    real Da = w2*( -2*ab.XX*ab.YY )/abnS;
    real da = w2*( ab.XX*ab.XX-ab.YY*ab.YY )/abnS;
    real Dc = w2*( -2*cd.XX*cd.YY )/cdnS;
    real dc = w2*(  cd.XX*cd.XX-cd.YY*cd.YY )/cdnS;
    real entrya[] = {-Da, -da, Da,  da}; //={d(na_x/la)/dxa, d(na_x/la)/dya, d(na_x/l)/dxb, ...}
    real entryc[] = { Dc,  dc,  -Dc,  -dc};//={ d(nc_x/lc)/dxc, d(nc_x/lc)/dyc, d(nc_x/l)/dxd, ...}
    int shifta = 0;
    int shiftc= 26;
    int mm;
    
    //Add second part of matrix.
    //The pos(-1, jj) accounts for the different signs of the matrix
    for ( int jj=0; jj <  4; ++jj) {
        for ( int ii=jj ; ii < 4; ++ii ) {
            m[ii + shifta] += pow(-1,jj)*entrya[ii-jj];
            m[ii + shiftc] += pow(-1,jj)*entryc[ii-jj];
        }
        shifta += 7 - jj;
        shiftc += 3 - jj;
    }
    
    
    //very Cumbersome!!!
    //Entries for Matrix mC  and vector vBAS
    //vBAS = fa - M*P0
    for ( int ii = 0; ii < 8; ++ii )
    {
        vBAS[index[ii]] += w2*h[ii];
        for (int jj = 0; jj < 8; ++jj) {
            if (jj < ii)
                mm = jj*(7.5-0.5*jj)+ii;
            else {
                mm = ii*(7.5-0.5*ii)+jj;
                mC( index[ii], index[jj] ) += m[mm];
            }
            vBAS[index[ii]] -= m[mm]*coord[jj];
        }
    }
}
#endif


//------------------------------------------------------------------------------
#pragma mark - Links
//------------------------------------------------------------------------------


/**
 Update Meca to include an interaction between A and B
 The force is linear with a zero resting length:
 @code
 force_A = weight * ( B - A )
 force_B = weight * ( A - B )
 @endcode
 
 In practice, Meca::interLink() will update the matrix mB,
 adding `weight` at the indices corresponding to `A` and `B`.
 
 Note: with modulo, the position of the fibers may be shifted in space,
 and a correction is necessary to make the force calculation correct:
 
 @code
 force_A = weight * ( B - A - offset )
 force_B = weight * ( A - B + offset )
 @endcode

 Here 'offset' is a multiple of the space periodicity, corresponding to B-A:
 offset = modulo->offset( A - B )

 In practice, Meca::interLink() will update the vector vBAS[]:
 @code
 vBAS[A] += weight * offset;
 vBAS[B] -= weight * offset;
 @endcode

 In principle, what goes to vBAS[] with modulo can be derived
 simply by multiplying the matrix block by 'offset'.
 */

void Meca::interLink(const PointExact & pta,
                     const PointExact & ptb, 
                     const real weight)
{
    assert_true( weight >= 0 );
    
    const index_type ia = pta.matIndex();
    const index_type ib = ptb.matIndex();

    if ( any_equal(ia, ib) )
        return;

    mB(ia, ia) -= weight;
    mB(ia, ib) += weight;
    mB(ib, ib) -= weight;
    
    if ( modulo )
    {
        Vector off = modulo->offset( pta.pos() - ptb.pos() );
        if ( !off.null() )
        {
            off.add_to(weight, vBAS+DIM*ia);
            off.sub_to(weight, vBAS+DIM*ib);
        }
    }
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(),ptb.pos());
    }
#endif
}


/**
Update Meca to include an interaction between A and B
 The force is linear with a zero resting length:
 @code
 force_A = weight * ( B - A )
 force_B = weight * ( A - B )
 @endcode
 */

void Meca::interLink(const PointInterpolated & pta, 
                     const PointExact & ptb, 
                     const real weight)
{
    assert_true( weight >= 0 );
    
    //the index of the points in the matrix mB:
    const index_type inx[] = { pta.matIndex1(), pta.matIndex2(), ptb.matIndex() };
    
    if ( any_equal(inx[0], inx[1], inx[2]) )
        return;

    //coefficients on the points:
    const real cc[] = { pta.coef2(), pta.coef1(), -1.0 };
    const real ww[] = { weight*cc[0], weight*cc[1], -weight };
    
    for ( int kk = 0;  kk < 3; ++kk )
    for ( int ll = kk; ll < 3; ++ll )
        mB(inx[kk], inx[ll]) -= cc[kk] * ww[ll];
    
    if ( modulo )
    {
        Vector off = modulo->offset( pta.pos() - ptb.pos() );
        if ( !off.null() )
        {
            off.add_to(ww[0], vBAS+DIM*inx[0]);
            off.add_to(ww[1], vBAS+DIM*inx[1]);
            off.add_to(ww[2], vBAS+DIM*inx[2]);
        }
    }
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(),ptb.pos());
    }
#endif
}


/**
Update Meca to include an interaction between A and B
 The force is linear with a zero resting length:
 @code
 force_A = weight * ( B - A )
 force_B = weight * ( A - B )
 @endcode
 */

void Meca::interLink(const PointInterpolated & pta,
                     const PointInterpolated & ptb,
                     const real weight)
{
    assert_true( weight >= 0 );
    
    //the index of the points in the matrix mB:
    const index_type inx[] = { pta.matIndex1(), pta.matIndex2(), ptb.matIndex1(), ptb.matIndex2() };
    
    if ( any_equal(inx[0], inx[1], inx[2], inx[3]) )
        return;
    
    //interpolation coefficients:
    const real cc[] = { pta.coef2(), pta.coef1(),  -ptb.coef2(), -ptb.coef1() };
    const real ww[] = { weight*cc[0], weight*cc[1], weight*cc[2], weight*cc[3] };

    if ( any_equal(inx[0], inx[1], inx[2], inx[3]) )
        return;
    
#if ( 1 )
    
    mB(inx[0], inx[0]) -= ww[0] * cc[0];
    mB(inx[1], inx[0]) -= ww[1] * cc[0];
    mB(inx[2], inx[0]) -= ww[2] * cc[0];
    mB(inx[3], inx[0]) -= ww[3] * cc[0];
    
    mB(inx[1], inx[1]) -= ww[1] * cc[1];
    mB(inx[2], inx[1]) -= ww[2] * cc[1];
    mB(inx[3], inx[1]) -= ww[3] * cc[1];
    
    mB(inx[2], inx[2]) -= ww[2] * cc[2];
    mB(inx[3], inx[2]) -= ww[3] * cc[2];
    
    mB(inx[3], inx[3]) -= ww[3] * cc[3];

#else
    
    for ( int jj = 0 ; jj < 4; ++jj )
    for ( int ii = jj; ii < 4; ++ii )
        mB(inx[ii], inx[jj]) -= cc[jj] * ww[ii];

#endif
    
    if ( modulo )
    {
        Vector off = modulo->offset( pta.pos() - ptb.pos() );
        if ( !off.null() )
        {
            off.add_to(ww[0], vBAS+DIM*inx[0]);
            off.add_to(ww[1], vBAS+DIM*inx[1]);
            off.add_to(ww[2], vBAS+DIM*inx[2]);
            off.add_to(ww[3], vBAS+DIM*inx[3]);
        }
    }
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(),ptb.pos());
    }
#endif
    
}

/**
 Update Meca to include an interaction between A and B
 The force is linear with a zero resting length:
 @code
 force_A = weight * ( B - A )
 force_B = weight * ( A - B )
 @endcode
 Point A is a PointExact.
 Point B in interpolated over 2 model-points of a Mecable, at index 'off'.
 Diagonal and lower elements of mB are set.
 */
void Meca::interLink2(const PointExact & pte, const unsigned off,
                      const index_type pts[2], const real coef[2],
                      const real weight)
{
    assert_true( weight >= 0 );
    
    index_type inx[] = { pte.matIndex(), off+pts[0], off+pts[1] };
    
    const real cc[] = {    -1.0,       coef[0],      coef[1] };
    const real ww[] = { -weight,  weight*cc[1], weight*cc[2] };
    
    assert_small(coef[0]+coef[1]-1.0);
    //assert_false( any_equal(inx[0], inx[1], inx[2]) );
    
    mB(inx[0], inx[0]) += ww[0];
    mB(inx[1], inx[0]) += ww[1];
    mB(inx[2], inx[0]) += ww[2];
    
    mB(inx[1], inx[1]) -= ww[1] * cc[1];
    mB(inx[2], inx[1]) -= ww[2] * cc[1];
    
    mB(inx[2], inx[2]) -= ww[2] * cc[2];
    
    if ( modulo )
        throw Exception("interLink(...) is not usable with periodic boundary conditions");
}

/**
 Update Meca to include an interaction between A and B
 The force is linear with a zero resting length:
 @code
 force_A = weight * ( B - A )
 force_B = weight * ( A - B )
 @endcode
 Point A is a PointInterpolated.
 Point B in interpolated over 2 model-points of a Mecable, at index 'off'.
 Diagonal and lower elements of mB are set.
 */
void Meca::interLink2(const PointInterpolated & pti, const unsigned off,
                      const index_type pts[2], const real coef[2],
                      const real weight)
{
    assert_true( weight >= 0 );
    
    index_type inx[] = { pti.matIndex1(), pti.matIndex2(), off+pts[0], off+pts[1] };
    
    const real cc[] = { -pti.coef2(), -pti.coef1(),      coef[0],      coef[1] };
    const real ww[] = { weight*cc[0], weight*cc[1], weight*cc[2], weight*cc[3] };

    assert_small(coef[0]+coef[1]-1.0);
    //assert_false( any_equal(inx[0], inx[1], inx[2], inx[3]) );
    
    mB(inx[0], inx[0]) -= ww[0] * cc[0];
    mB(inx[1], inx[0]) -= ww[1] * cc[0];
    mB(inx[2], inx[0]) -= ww[2] * cc[0];
    mB(inx[3], inx[0]) -= ww[3] * cc[0];

    mB(inx[1], inx[1]) -= ww[1] * cc[1];
    mB(inx[2], inx[1]) -= ww[2] * cc[1];
    mB(inx[3], inx[1]) -= ww[3] * cc[1];

    mB(inx[2], inx[2]) -= ww[2] * cc[2];
    mB(inx[3], inx[2]) -= ww[3] * cc[2];

    mB(inx[3], inx[3]) -= ww[3] * cc[3];
    
    if ( modulo )
        throw Exception("interLink(...) is not usable with periodic boundary conditions");
}


/**
 Update Meca to include an interaction between A and B
 The force is linear with a zero resting length:
 @code
 force_A = weight * ( B - A )
 force_B = weight * ( A - B )
 @endcode
 Point A is a PointExact.
 Point B in interpolated over 3 model-points of a Mecable, at index 'off'.
 Diagonal and lower elements of mB are set.
 */
void Meca::interLink3(const PointExact & pte, const unsigned off,
                      const index_type pts[3], const real coef[3],
                      const real weight)
{
    assert_true( weight >= 0 );
    
    index_type inx[] = { pte.matIndex(), off+pts[0], off+pts[1], off+pts[2] };
    
    const real cc[] = {    -1.0,        coef[0],        coef[1],        coef[2] };
    const real ww[] = { -weight, weight*coef[0], weight*coef[1], weight*coef[2] };

    assert_small(coef[0]+coef[1]+coef[2]-1.0);
    
    mB(inx[0], inx[0]) += ww[0];
    mB(inx[1], inx[0]) += ww[1];
    mB(inx[2], inx[0]) += ww[2];
    mB(inx[3], inx[0]) += ww[3];
    
    mB(inx[1], inx[1]) -= ww[1] * cc[1];
    mB(inx[2], inx[1]) -= ww[2] * cc[1];
    mB(inx[3], inx[1]) -= ww[3] * cc[1];
    
    mB(inx[2], inx[2]) -= ww[2] * cc[2];
    mB(inx[3], inx[2]) -= ww[3] * cc[2];
    
    mB(inx[3], inx[3]) -= ww[3] * cc[3];
    
    if ( modulo )
        throw Exception("interLink(...) is not usable with periodic boundary conditions");
}



/**
 Update Meca to include an interaction between A and B
 The force is linear with a zero resting length:
 @code
 force_A = weight * ( B - A )
 force_B = weight * ( A - B )
 @endcode
 Point A is a PointInterpolated.
 Point B in interpolated over 4 model-points of a Mecable, at index 'off'.
 Diagonal and lower elements of mB are set.
*/
void Meca::interLink3(const PointInterpolated & pti, const unsigned off,
                      const index_type pts[3], const real coef[3],
                      const real weight)
{
    assert_true( weight >= 0 );
    
    index_type inx[] = { pti.matIndex1(), pti.matIndex2(), off+pts[0], off+pts[1], off+pts[2] };
    
    const real cc[] = { -pti.coef2(), -pti.coef1(),        coef[0],        coef[1],        coef[2] };
    const real ww[] = { weight*cc[0], weight*cc[1], weight*coef[0], weight*coef[1], weight*coef[2] };
    
    assert_small(coef[0]+coef[1]+coef[2]-1.0);
    
    mB(inx[0], inx[0]) -= ww[0] * cc[0];
    mB(inx[1], inx[0]) -= ww[1] * cc[0];
    mB(inx[2], inx[0]) -= ww[2] * cc[0];
    mB(inx[3], inx[0]) -= ww[3] * cc[0];
    mB(inx[4], inx[0]) -= ww[4] * cc[0];
    
    mB(inx[1], inx[1]) -= ww[1] * cc[1];
    mB(inx[2], inx[1]) -= ww[2] * cc[1];
    mB(inx[3], inx[1]) -= ww[3] * cc[1];
    mB(inx[4], inx[1]) -= ww[4] * cc[1];
    
    mB(inx[2], inx[2]) -= ww[2] * cc[2];
    mB(inx[3], inx[2]) -= ww[3] * cc[2];
    mB(inx[4], inx[2]) -= ww[4] * cc[2];
    
    mB(inx[3], inx[3]) -= ww[3] * cc[3];
    mB(inx[4], inx[3]) -= ww[4] * cc[3];
    
    mB(inx[4], inx[4]) -= ww[4] * cc[4];
    
    if ( modulo )
        throw Exception("interLink(...) is not usable with periodic boundary conditions");
}

/**
 Update Meca to include an interaction between A and B
 The force is linear with a zero resting length:
 @code
 force_A = weight * ( B - A )
 force_B = weight * ( A - B )
 @endcode
 Point A is a PointExact.
 Point B in interpolated over 4 model-points of a Mecable, at index 'off'.
 Diagonal and lower elements of mB are set.
 */
void Meca::interLink4(const PointExact & pte, const unsigned off,
                      const index_type pts[4], const real coef[4],
                      const real weight)
{
    assert_true( weight >= 0 );
    
    index_type inx[] = { pte.matIndex(), off+pts[0], off+pts[1], off+pts[2], off+pts[3] };
    
    const real cc[] = {    -1.0,        coef[0],        coef[1],        coef[2],        coef[3] };
    const real ww[] = { -weight, weight*coef[0], weight*coef[1], weight*coef[2], weight*coef[3] };
    
    assert_small(coef[0]+coef[1]+coef[2]+coef[3]-1.0);
    
    mB(inx[0], inx[0]) += ww[0];
    mB(inx[1], inx[0]) += ww[1];
    mB(inx[2], inx[0]) += ww[2];
    mB(inx[3], inx[0]) += ww[3];
    mB(inx[4], inx[0]) += ww[4];
    
    mB(inx[1], inx[1]) -= ww[1] * cc[1];
    mB(inx[2], inx[1]) -= ww[2] * cc[1];
    mB(inx[3], inx[1]) -= ww[3] * cc[1];
    mB(inx[4], inx[1]) -= ww[4] * cc[1];
    
    mB(inx[2], inx[2]) -= ww[2] * cc[2];
    mB(inx[3], inx[2]) -= ww[3] * cc[2];
    mB(inx[4], inx[2]) -= ww[4] * cc[2];
    
    mB(inx[3], inx[3]) -= ww[3] * cc[3];
    mB(inx[4], inx[3]) -= ww[4] * cc[3];
    
    mB(inx[4], inx[4]) -= ww[4] * cc[4];
    
    if ( modulo )
        throw Exception("interLink(...) is not usable with periodic boundary conditions");
}


/**
 Update Meca to include an interaction between A and B
 The force is linear with a zero resting length:
 @code
 force_A = weight * ( B - A )
 force_B = weight * ( A - B )
 @endcode
 Point A is a PointInterpolated.
 Point B in interpolated over 4 model-points of a Mecable, at index 'off'.
 Diagonal and lower elements of mB are set.
*/
void Meca::interLink4(const PointInterpolated & pti, const unsigned off,
                      const index_type pts[4], const real coef[4],
                      const real weight)
{
    assert_true( weight >= 0 );
    
    index_type inx[] = { pti.matIndex1(), pti.matIndex2(), off+pts[0], off+pts[1], off+pts[2], off+pts[3] };
    
    const real cc[] = { -pti.coef2(), -pti.coef1(),        coef[0],        coef[1],        coef[2],        coef[3] };
    const real ww[] = { weight*cc[0], weight*cc[1], weight*coef[0], weight*coef[1], weight*coef[2], weight*coef[3] };
    
    assert_small(coef[0]+coef[1]+coef[2]+coef[3]-1.0);
    
    mB(inx[0], inx[0]) -= ww[0] * cc[0];
    mB(inx[1], inx[0]) -= ww[1] * cc[0];
    mB(inx[2], inx[0]) -= ww[2] * cc[0];
    mB(inx[3], inx[0]) -= ww[3] * cc[0];
    mB(inx[4], inx[0]) -= ww[4] * cc[0];
    mB(inx[5], inx[0]) -= ww[5] * cc[0];
    
    mB(inx[1], inx[1]) -= ww[1] * cc[1];
    mB(inx[2], inx[1]) -= ww[2] * cc[1];
    mB(inx[3], inx[1]) -= ww[3] * cc[1];
    mB(inx[4], inx[1]) -= ww[4] * cc[1];
    mB(inx[5], inx[1]) -= ww[5] * cc[1];
    
    mB(inx[2], inx[2]) -= ww[2] * cc[2];
    mB(inx[3], inx[2]) -= ww[3] * cc[2];
    mB(inx[4], inx[2]) -= ww[4] * cc[2];
    mB(inx[5], inx[2]) -= ww[5] * cc[2];
    
    mB(inx[3], inx[3]) -= ww[3] * cc[3];
    mB(inx[4], inx[3]) -= ww[4] * cc[3];
    mB(inx[5], inx[3]) -= ww[5] * cc[3];
    
    mB(inx[4], inx[4]) -= ww[4] * cc[4];
    mB(inx[5], inx[4]) -= ww[5] * cc[4];

    mB(inx[5], inx[5]) -= ww[5] * cc[5];

    if ( modulo )
        throw Exception("interLink(...) is not usable with periodic boundary conditions");
}

//------------------------------------------------------------------------------
#pragma mark - Long Links
//------------------------------------------------------------------------------


void setInterLongLink1(real T[DIM*DIM], Vector const& ab)
{
    T[0      ] = ab.XX * ab.XX;
#if ( DIM > 1 )
    T[1      ] = ab.XX * ab.YY;
    T[1+DIM  ] = ab.YY * ab.YY;
#endif
#if ( DIM > 2 )
    T[2      ] = ab.XX * ab.ZZ;
    T[2+DIM  ] = ab.YY * ab.ZZ;
    T[2+DIM*2] = ab.ZZ * ab.ZZ;
#endif
}

void setInterLongLink2(real T[DIM*DIM], Vector const& ab, const real len)
{
    T[0      ] = 1.0 + len * ( ab.XX * ab.XX - 1.0 );
#if ( DIM > 1 )
    T[1      ] = len * ab.XX * ab.YY;
    T[1+DIM  ] = 1.0 + len * ( ab.YY * ab.YY - 1.0 );
#endif
#if ( DIM > 2 )
    T[2      ] = len * ab.XX * ab.ZZ;
    T[2+DIM  ] = len * ab.YY * ab.ZZ;
    T[2+DIM*2] = 1.0 + len * ( ab.ZZ * ab.ZZ - 1.0 );
#endif
}

/**
 Update Meca to include an interaction between A and B,
 The force is affine with non-zero resting length: 
 @code
 force_A = weight * ( B - A ) * ( length / |AB| - 1 )
 force_B = weight * ( A - B ) * ( length / |AB| - 1 )
 @endcode
 */

void Meca::interLongLink(const PointExact & pta, 
                         const PointExact & ptb,
                         real len, 
                         const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    const index_type ia = DIM * pta.matIndex();  // coef is +weight
    const index_type ib = DIM * ptb.matIndex();  // coef is -weight

    if ( any_equal(ia, ib) )
        return;

#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(), ptb.pos(), len);
    }
#endif
    
    Vector offset, ab = ptb.pos() - pta.pos();

    if ( modulo )
        offset = modulo->foldOffset(ab);
    
    const real abn = ab.norm();
    if ( abn < REAL_EPSILON )
        return;
    ab /= abn;
    
    ab.add_to(-weight*len, vBAS+ia);
    ab.add_to( weight*len, vBAS+ib);
    
    len /= abn;
    
    /* To stabilize the matrix with compression, we remove the negative eigenvalues
        This is done by using len = 1 in the formula if len > 1.0. */
    const bool cooked = ( len > 1.0 );
    
    real T[DIM*DIM];

    if ( cooked )
        setInterLongLink1(T, ab);
    else
        setInterLongLink2(T, ab, len);
    
    for ( int x = 0; x < DIM; ++x )
    for ( int y = x; y < DIM; ++y )
        T[y+DIM*x] *= weight;

#ifdef USE_MATRIX_BLOCK
    mC.block(ia, ia).sub_lower(T);
    mC.block(ia, ib).add_sym(T);
    mC.block(ib, ib).sub_lower(T);
#else
    for ( int x = 0; x < DIM; ++x )
    {
        mC(ia+x, ia+x) -= T[x+DIM*x];
        mC(ia+x, ib+x) += T[x+DIM*x];
        mC(ib+x, ib+x) -= T[x+DIM*x];
        for ( int y = x+1; y < DIM; ++y )
        {
            mC(ia+x, ia+y) -= T[y+DIM*x];
            mC(ia+x, ib+y) += T[y+DIM*x];
            mC(ib+x, ia+y) += T[y+DIM*x];
            mC(ib+x, ib+y) -= T[y+DIM*x];
        }
    }
#endif
    
    if ( modulo && !offset.null() )
    {
        real s = offset * ab;  //scalar product
        Vector v;
        if ( cooked )
            v = ( weight * s ) * ab;
        else
            v = ( weight * len * s ) * ab + ( weight * ( 1.0 - len )) * offset;
        v.sub_to(vBAS+ia);
        v.add_to(vBAS+ib);
    }
}


/**
 Update Meca to include an interaction between A and B,
 The force is affine with non-zero resting length: 
 @code
 force_A = weight * ( B - A ) * ( len / |AB| - 1 )
 force_B = weight * ( A - B ) * ( len / |AB| - 1 )
 @endcode
 */

void Meca::interLongLink(const PointInterpolated & pta, 
                         const PointExact & pte, 
                         real len,
                         const real weight )
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    //index in the matrix mC:
    const index_type inx[] = { DIM*pta.matIndex1(), DIM*pta.matIndex2(), DIM*pte.matIndex() };
    
    if ( any_equal(inx[0], inx[1], inx[2]) )
        return;

    //force coefficients on the points:
    const real  c[] = { pta.coef2(), pta.coef1(), -1.0 };
    const real cw[] = { weight*c[0], weight*c[1], -weight };

#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(), pte.pos(), len);
    }
#endif
    
    Vector offset, ab = pte.pos() - pta.pos();

    if ( modulo )
        offset = modulo->foldOffset(ab);

    const real abn = ab.norm();
    if ( abn < REAL_EPSILON ) return;
    ab /= abn;
    
    for ( int ii = 0; ii < 3; ++ii )
        ab.add_to(-len*cw[ii], vBAS+inx[ii]);
    
    len /= abn;
    
    /* To stabilize the matrix with compression, we remove the negative eigenvalues
        This is done by using len = 1 in the formula if len > 1.0. */
    const bool cooked = ( len > 1.0 );
    
    real T[DIM*DIM];
    
    if ( cooked )
        setInterLongLink1(T, ab);
    else
        setInterLongLink2(T, ab, len);

#ifdef USE_MATRIX_BLOCK
    for ( int kk = 0; kk < 3; ++kk )
    {
        mC.block(inx[kk], inx[kk]).add_lower(-c[kk] * cw[kk], T);
        for ( int ll = kk+1; ll < 3; ++ll )
            mC.block(inx[kk], inx[ll]).add_sym(-c[kk] * cw[ll], T);
    }
#else
    for ( int x = 0; x < DIM; ++x )
    {
        for ( int kk =  0; kk < 3; ++kk )
        for ( int ll = kk; ll < 3; ++ll )
            mC(inx[kk]+x, inx[ll]+x) -= c[kk] * cw[ll] * T[x+DIM*x];
        
        for ( int y = x+1; y < DIM; ++y )
        {
            for ( int kk = 0; kk < 3; ++kk )
            for ( int ll = 0; ll < 3; ++ll )
                mC(inx[kk]+y, inx[ll]+x) -= c[kk] * cw[ll] * T[y+DIM*x];
        }
    }
#endif
    
    if ( modulo && !offset.null() )
    {
        real s = offset * ab;  //scalar product
        Vector vec;
        if ( cooked )
            vec = s * ab;
        else
            vec = ( len * s ) * ab + ( 1.0 - len ) * offset;
        for ( int jj = 0; jj < 3; ++jj )
            vec.add_to(-cw[jj], vBAS+inx[jj]);
    }
}



/**
 Update Meca to include an interaction between A and B,
 The force is affine with non-zero resting length: 
 @code
 force_A = weight * ( B - A ) * ( len / |AB| - 1 )
 force_B = weight * ( A - B ) * ( len / |AB| - 1 )
 @endcode
 */

void Meca::interLongLink(const PointInterpolated & pta, 
                         const PointInterpolated & ptb, 
                         real len, 
                         const real weight )
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    //index in the matrix mC:
    const index_type inx[] = { DIM*pta.matIndex1(), DIM*pta.matIndex2(),
                               DIM*ptb.matIndex1(), DIM*ptb.matIndex2() };

    if ( any_equal(inx[0], inx[1], inx[2], inx[3]) )
        return;
    
    //force coefficients on the points:
    const real  c[] = { pta.coef2(), pta.coef1(), -ptb.coef2(), -ptb.coef1() };
    const real cw[] = { weight*c[0], weight*c[1], weight*c[2], weight*c[3] };
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(),ptb.pos(), len);
    }
#endif
    Vector offset, ab = ptb.pos() - pta.pos();

    if ( modulo )
        offset = modulo->foldOffset(ab);

    const real abn = ab.norm();
    if ( abn < REAL_EPSILON ) return;
    ab /= abn;
    
    for ( int ii = 0; ii < 4; ++ii )
        ab.add_to(-cw[ii]*len, vBAS+inx[ii]);

    len /= abn;
    
    /* To stabilize the matrix with compression, we remove the negative eigenvalues
        This is done by using len = 1 in the formula if len > 1.0. */
    const bool cooked = ( len > 1.0 );
    
    real T[DIM*DIM];
    
    if ( cooked )
        setInterLongLink1(T, ab);
    else
        setInterLongLink2(T, ab, len);
    
#ifdef USE_MATRIX_BLOCK
    for ( int kk = 0; kk < 4; ++kk )
    {
        mC.block(inx[kk], inx[kk]).add_lower(-c[kk] * cw[kk], T);
        for ( int ll = kk+1; ll < 4; ++ll )
            mC.block(inx[kk], inx[ll]).add_sym(-c[kk] * cw[ll], T);
    }
#else
    for ( int x = 0; x < DIM; ++x )
    {
        for ( int kk =  0; kk < 4; ++kk )
        for ( int ll = kk; ll < 4; ++ll )
            mC(inx[kk]+x, inx[ll]+x) -= c[kk] * cw[ll] * T[x+DIM*x];

        for ( int y = x+1; y < DIM; ++y )
        {
            for ( int kk = 0; kk < 4; ++kk )
            for ( int ll = 0; ll < 4; ++ll )
                mC(inx[kk]+y, inx[ll]+x) -= c[kk] * cw[ll] * T[y+DIM*x];
        }
    }
#endif
    
    if ( modulo && !offset.null() )
    {
        real s = offset * ab;  //scalar product
        if ( cooked )
        {
            for ( int jj = 0; jj < 4; ++jj )
                ab.add_to(-s*cw[jj], vBAS+inx[jj]);
        }
        else
        {
            Vector v = ( len * s ) * ab + ( 1.0 - len ) * offset;
            for ( int jj = 0; jj < 4; ++jj )
                v.add_to(-cw[jj], vBAS+inx[jj]);
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Side Links
//------------------------------------------------------------------------------

/**
 Update Meca to include an interaction between A and B,
 but it is taken between B and a point S located on the side of A:
 S = A + len * N,
 where N is a vector of unit norm that is orthogonal to the fiber in A.
 S is therefore linearly related to the two model points on the sides of an.
 The force is linear of zero resting length:
 @code
 force_S = weight * ( S - B )
 force_B = weight * ( B - S )
 @endcode
 */


#if ( DIM == 2 )

void Meca::interSideLink2D(const PointInterpolated & pta, 
                           const PointExact & ptb, 
                           const real arm,
                           const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    index_type ia1 = pta.matIndex1();
    index_type ia2 = pta.matIndex2();
    index_type ib  = ptb.matIndex();

    if ( any_equal(ia1, ia2, ib) )
        return;

    //force coefficients on the points:
    const real ca1 = pta.coef2();
    const real ca2 = pta.coef1();
    const real eps = arm / pta.len();
    
    const real ca1w = weight * ca1;
    const real ca2w = weight * ca2;
    const real epsw = weight * eps;
    const real epsepsw = eps * epsw;
    
    //we put the isotropic terms in mB
    mB(ia1, ia1) -=  ca1w * ca1 + epsepsw;
    mB(ia1, ia2) -=  ca1w * ca2 - epsepsw;
    mB(ia2, ia2) -=  ca2w * ca2 + epsepsw;
    
    mB(ib,  ib) -=  weight;
    mB(ia1, ib) +=  ca1w;
    mB(ia2, ib) +=  ca2w;
    
    //index in the matrix mC:
    ia1 *= DIM;
    ia2 *= DIM;
    ib  *= DIM;
    
    mC(ia1  , ia2+1) += epsw;
    mC(ia1+1, ia2  ) -= epsw;
    
    mC(ia1  , ib+1 ) -= epsw;
    mC(ia1+1, ib   ) += epsw;
    mC(ia2,   ib+1 ) += epsw;
    mC(ia2+1, ib   ) -= epsw;
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(), pta.pos()+cross(arm,pta.dir()), ptb.pos());
    }
#endif
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptb.pos() - pta.pos() );
        
        real offx = off.XX;
        if ( offx )
        {
            vBAS[ia1  ] += ca1w * offx;
            vBAS[ia1+1] += epsw * offx;
            vBAS[ia2  ] += ca2w * offx;
            vBAS[ia2+1] -= epsw * offx;
            vBAS[ib   ] += offx;
        }
        real offy = off.YY;
        if ( offy )
        {
            vBAS[ia1  ] -= epsw * offy;
            vBAS[ia1+1] += ca1w * offy;
            vBAS[ia2  ] += epsw * offy;
            vBAS[ia2+1] += ca2w * offy;
            vBAS[ib +1] += offy;
        }
    }
}

#elif ( DIM == 3 )


/**
 This is experimental and should not be used.
 
 This creates a interaction between `ptb`, 
 and a point S which is on the side of `pta`.
 @code
 S = pos_a + arm ^ dir_a
 @endcode
 arm must be perpendicular to link
 */
void Meca::interSideLink3D(const PointInterpolated & pta, 
                           const PointExact & ptb, 
                           const Vector & arm,
                           const real weight)
{
    assert_true( weight >= 0 );

    // indices to mC:
    const index_type ia1 = pta.matIndex1();
    const index_type ia2 = pta.matIndex2();
    const index_type ib  = ptb.matIndex();
    
    if ( any_equal(ia1, ia2, ib) )
        return;
    
    const index_type inx[6] = { DIM*ia1, DIM*ia1+1, DIM*ia1+2, DIM*ia2, DIM*ia2+1, DIM*ia2+2 };

    real a = pta.coef2();
    real b = pta.coef1();
    real s = 1.0 / pta.len();
    
    real ex = s * arm.XX;
    real ey = s * arm.YY;
    real ez = s * arm.ZZ;
    
    /* The transfer matrix transforms the two PointExact in pta,
     to the side point S:
     S = aa * pt1 + bb * pt2 + arm ^ ( pt2 - pt1 ).normalized
     
     It was generated in Maxima:
     MVP: matrix([0, -ez, ey], [ez, 0, -ex], [-ey, ex, 0]);
     MD: addcol(-ident(3), ident(3));
     MC: addcol(aa*ident(3), bb*ident(3));
     T: MC+MVP.MD;
     */
    const real T[18] = {
        a,  ez, -ey,   b, -ez,  ey,
      -ez,   a,  ex,  ez,   b, -ex,
       ey, -ex,   a, -ey,  ex,   b
    };
    
    real a2 = a * a, ab = a * b;
    real b2 = b * b;
    
    real exx = ex * ex, exy = ex*ey, exz = ex*ez;
    real eyy = ey * ey, eyz = ey*ez;
    real ezz = ez * ez;
    
    // TT = transpose(T) * T is symmetric, and thus we only set half of it:
    /* Maxima code:
     TT: expand(transpose(T) . T);
     */
    real TT[36] = {
        eyy+ezz+a2,  0,           0,           0,           0,           0,
        -exy,        exx+ezz+a2,  0,           0,           0,           0,
        -exz,       -eyz,         exx+eyy+a2,  0,           0,           0,
        -ezz-eyy+ab, ez+exy,      exz-ey,      eyy+ezz+b2,  0,           0,
        -ez+exy,    -ezz-exx+ab,  eyz+ex,     -exy,         exx+ezz+b2,  0,
        exz+ey,      eyz-ex,     -eyy-exx+ab, -exz,        -eyz,         exx+eyy+b2
    };
    
    // we project to bring all forces in the plane perpendicular to 'arm'
    //real sca = 1.0 / arm.norm();
    //real an = a * sca;
    //real bn = b * sca;
    // Maxima code: matrix([ex, ey, ez]) . T;
    //real TP[9] = { an*ex, an*ey, an*ez, bn*ex, bn*ey, bn*ez, -ex, -ey, -ez };    
    //blas_xgemm('N','N', 6, 1, 3, sca, T, 6, arm, 3, 0.0, TP, 6);
    
    //blas_xsyrk('U','N', 6, 1, weight, TP, 6, -weight, TT, 6);
    blas_xscal(36, -weight, TT, 1); 
    
    for ( int ii=0; ii<6; ++ii )
        for ( int jj=ii; jj<6; ++jj )
            mC(inx[ii], inx[jj]) += TT[ii+6*jj];
    
    //mB(ia1, ib) += -a * weight;
    //mB(ia2, ib) += -b * weight;
    mB(ib, ib) -= weight;
    
    for ( int ii=0; ii<6; ++ii )
        for ( int jj=0; jj<3; ++jj )
            mC(inx[ii], DIM*ib+jj) += weight * T[ii+6*jj];
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(), pta.pos()+cross(arm,pta.dir()), ptb.pos());
    }
#endif
    
    if ( modulo )
        throw Exception("interSideLink3D is not usable with periodic boundary conditions");
}


void Meca::interSideLinkS(const PointInterpolated & pta, 
                          const PointExact & ptb, 
                          const Vector & arm,
                          const real len,
                          const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len > REAL_EPSILON );

    const index_type inx[] = { DIM*pta.matIndex1(), DIM*pta.matIndex2(),
                               DIM*ptb.matIndex() };

    if ( any_equal(inx[0], inx[1], inx[2]) )
        return;

    // vector 'a' is parallel to first Fiber
    Vector a = pta.dir();
    Vector b = arm / len;
    
    // we can set directly the interaction coefficient matrix:
    const real xx = a.XX*a.XX+b.XX*b.XX, xy = a.XX*a.YY+b.XX*b.YY, xz = a.XX*a.ZZ+b.XX*b.ZZ;
    const real yy = a.YY*a.YY+b.YY*b.YY, yz = a.YY*a.ZZ+b.YY*b.ZZ, zz = a.ZZ*a.ZZ+b.ZZ*b.ZZ;
    const real T[9] = { xx, xy, xz, xy, yy, yz, xz, yz, zz };
    
    // we set directly the transformed offset vector:
    const real RB[3] = { arm.XX, arm.YY, arm.ZZ };
    
    // weights and indices:
    const real  c[3] = { pta.coef2(), pta.coef1(), -1.0 };
    const real cw[3] = { -weight*c[0], -weight*c[1], weight };
    
    // fill the matrix mC
    for ( int ii=0; ii<3; ++ii )
    {
        for ( int x = 0; x < DIM; ++x )
            vBAS[inx[ii]+x] += cw[ii] * RB[x];
        
        const real g = c[ii] * cw[ii];
#ifdef USE_MATRIX_BLOCK
        mC.block(inx[ii], inx[ii]).add_lower(g, T);
#else
        for ( int x = 0; x < DIM; ++x )
        for ( int y = x; y < DIM; ++y )
            mC(inx[ii]+y, inx[ii]+x) += g * T[y+3*x];
#endif

        for ( int jj=ii+1; jj<3; ++jj )
        {
            const real h = c[ii] * cw[jj];
#ifdef USE_MATRIX_BLOCK
            mC.block(inx[ii], inx[jj]).add(h, T);
#else
            for ( int x = 0; x < DIM; ++x )
            for ( int y = 0; y < DIM; ++y )
                mC(inx[ii]+y, inx[jj]+x) += h * T[y+3*x];
#endif
        }
    }
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(), pta.pos()+arm, ptb.pos());
    }
#endif
    
    if ( modulo )
        throw Exception("interSideLinkS is not usable with periodic boundary conditions");
}
#endif


void Meca::interSideLink(const PointInterpolated & pta, 
                         const PointExact & ptb, 
                         const real len,
                         const real weight )
{
#if ( DIM == 1 )
    
    throw Exception("Meca::interSideLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    real arm = len * RNG.sign_exc( cross(pta.diff(), ptb.pos()-pta.pos()));
    interSideLink2D(pta, ptb, arm, weight);

#elif ( DIM == 3 )
    
    // 'arm' is perpendicular to A-Fiber and link:
    Vector a   = pta.diff();
    Vector as  = ptb.pos() - pta.pos();
    Vector arm = as - ( ( as * a ) / a.normSqr() ) * a;
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        interSideLinkS(pta, ptb, arm * (len / n), len, weight);
    
#endif
}


#if ( DIM == 2 )

void Meca::interSideLink2D(const PointInterpolated & pta, 
                           const PointInterpolated & ptb, 
                           const real arm,
                           const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    index_type  ia1 = pta.matIndex1(), ia2 = pta.matIndex2();
    index_type  ib1 = ptb.matIndex1(), ib2 = ptb.matIndex2();

    if ( any_equal(ia1, ia2, ib1, ib2) )
        return;

    const real ca1 =  pta.coef2(), ca2 =  pta.coef1();
    const real cb1 = -ptb.coef2(), cb2 = -ptb.coef1();
    const real eps = arm / pta.len();
    
    const real w = -weight;
    const real ca1w = ca1 * w, ca2w = ca2 * w;
    const real cb1w = cb1 * w, cb2w = cb2 * w;
    const real epsw = eps * w, epsepsw = eps * epsw;
    
    // isotropic terms go in matrix mB:
    mB(ia1, ia1) += ca1w * ca1 + epsepsw;
    mB(ia1, ia2) += ca1w * ca2 - epsepsw;
    mB(ia2, ia2) += ca2w * ca2 + epsepsw;
    
    mB(ib1, ib1) += cb1w * cb1;
    mB(ib1, ib2) += cb1w * cb2;
    mB(ib2, ib2) += cb2w * cb2;
    
    mB(ia1, ib1) += ca1w * cb1;
    mB(ia1, ib2) += ca1w * cb2;
    mB(ia2, ib1) += ca2w * cb1;
    mB(ia2, ib2) += ca2w * cb2;
    
    // update indices to address matrix mC:
    ia1 *= DIM;
    ia2 *= DIM;
    ib1 *= DIM;
    ib2 *= DIM;
    
    // anistropic terms go in the matrix mC:
    mC(ia1  , ia2+1) -= epsw;
    mC(ia1+1, ia2  ) += epsw;
    
    const real epscb1w = epsw * cb1;
    const real epscb2w = epsw * cb2;
    
    mC(ia1, ib1+1) -=  epscb1w;
    mC(ia1, ib2+1) -=  epscb2w;
    
    mC(ia1+1, ib1) +=  epscb1w;
    mC(ia1+1, ib2) +=  epscb2w;
    
    mC(ia2, ib1+1) +=  epscb1w;
    mC(ia2, ib2+1) +=  epscb2w;
    
    mC(ia2+1, ib1) -=  epscb1w;
    mC(ia2+1, ib2) -=  epscb2w;
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptb.pos() - pta.pos() );
        
        real offx = off.XX;
        if ( offx )
        {
            vBAS[ia1  ] += ca1w * offx;
            vBAS[ia1+1] += epsw * offx;
            vBAS[ia2  ] += ca2w * offx;
            vBAS[ia2+1] -= epsw * offx;
            vBAS[ib1  ] += cb1w * offx;
            vBAS[ib2  ] += cb2w * offx;
        }
        real offy = off.YY;
        if ( offy )
        {
            vBAS[ia1  ] -= epsw * offy;
            vBAS[ia1+1] += ca1w * offy;
            vBAS[ia2  ] += epsw * offy;
            vBAS[ia2+1] += ca2w * offy;
            vBAS[ib1+1] += cb1w * offy;
            vBAS[ib2+1] += cb2w * offy;
        }
    }
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(), pta.pos()+cross(arm,pta.diff()), ptb.pos());
    }
#endif
}

#elif ( DIM == 3 )

void Meca::interSideLinkS(const PointInterpolated & pta, 
                          const PointInterpolated & ptb,
                          const Vector & arm,
                          const real len,
                          const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len > REAL_EPSILON );

    const index_type inx[] = { DIM*pta.matIndex1(), DIM*pta.matIndex2(),
                               DIM*ptb.matIndex1(), DIM*ptb.matIndex2() };
    
    if ( any_equal(inx[0], inx[1], inx[2], inx[3]) )
        return;

    Vector a = pta.dir();
    Vector b = arm / len;
    // Vector c = cross(a, b);
    
    // we can set directly the interaction coefficient matrix:
    const real xx = a.XX*a.XX+b.XX*b.XX, xy = a.XX*a.YY+b.XX*b.YY, xz = a.XX*a.ZZ+b.XX*b.ZZ;
    const real yy = a.YY*a.YY+b.YY*b.YY, yz = a.YY*a.ZZ+b.YY*b.ZZ, zz = a.ZZ*a.ZZ+b.ZZ*b.ZZ;
    const real T[9] = { xx, xy, xz, xy, yy, yz, xz, yz, zz };
    
    // we set directly the transformed offset vector:
    const real RB[3] = { arm.XX, arm.YY, arm.ZZ };
    
    // weights and indices:
    const real  c[4] = { pta.coef2(), pta.coef1(), -ptb.coef2(), -ptb.coef1() };
    const real cw[4] = { -weight*c[0], -weight*c[1], -weight*c[2], -weight*c[3] };
    
    // fill the matrix mC
    for ( int ii=0; ii<4; ++ii )
    {
        for ( int x = 0; x < DIM; ++x )
            vBAS[inx[ii]+x] += cw[ii] * RB[x];
        
        const real g = c[ii] * cw[ii];
#ifdef USE_MATRIX_BLOCK
        mC.block(inx[ii], inx[ii]).add_lower(g, T);
#else
        for ( int x = 0; x < DIM; ++x )
        for ( int y = x; y < DIM; ++y )
            mC(inx[ii]+y, inx[ii]+x) += g * T[y+3*x];
#endif
        
        for ( int jj=ii+1; jj<4; ++jj )
        {
            const real h = c[ii] * cw[jj];
#ifdef USE_MATRIX_BLOCK
            mC.block(inx[ii], inx[jj]).add(h, T);
#else
            for ( int x = 0; x < DIM; ++x )
            for ( int y = 0; y < DIM; ++y )
                mC(inx[ii]+y, inx[jj]+x) += h * T[y+3*x];
#endif
            
        }
    }
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(), pta.pos()+arm, ptb.pos());
    }
#endif
    
    if ( modulo )
        throw Exception("interSideLinkS is not usable with periodic boundary conditions");
}

#endif


/**
 Update Meca to include an interaction between A and B,
 Which is taken between B and a point S located on the side of A:
 S = A + len * N,
 where N is a normalized vector orthogonal to the fiber in an.
 S is linearly related to the two model points on the sides of A, P1 and P2
 In 3D S is choosen in the plane of P1, P2 and B.
 The force is linear of zero resting length:
 @code
 force_S = weight * ( S - B )
 force_B = weight * ( B - S )
 @endcode
 */

void Meca::interSideLink(const PointInterpolated & pta, 
                         const PointInterpolated & ptb, 
                         const real len,
                         const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::interSideLink is meaningless in 1D");

#elif ( DIM == 2 )
    
    real arm = len * RNG.sign_exc( cross(pta.diff(), ptb.pos()-pta.pos()) );
    interSideLink2D(pta, ptb, arm, weight);
    
#elif ( DIM == 3 )
    
    // 'arm' is perpendicular to A-Fiber and link:
    Vector a   = pta.diff();
    Vector as  = ptb.pos() - pta.pos();
    Vector arm = as - ( ( as * a ) / a.normSqr() ) * a;
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        interSideLinkS(pta, ptb, arm * (len / n), len, weight);
    
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Side Side Links
//------------------------------------------------------------------------------

#if ( DIM == 2 )

void Meca::interSideSideLink2D(const PointInterpolated & pta,
                               const PointInterpolated & ptb, 
                               const real len,
                               const real weight,
                               int side1, int side2 )
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );
    
    //index in the matrix mB:
    index_type ia1 = pta.matIndex1(), ia2 = pta.matIndex2();
    index_type ib1 = ptb.matIndex1(), ib2 = ptb.matIndex2();
    
    if ( any_equal(ia1, ia2, ib1, ib2) )
        return;

    const real ca1 =  pta.coef2(), ca2 =  pta.coef1();
    const real cb1 = -ptb.coef2(), cb2 = -ptb.coef1();
    
    const real ee1 = side1 * len / ( 2 * pta.len() );
    const real ee2 = side2 * len / ( 2 * ptb.len() );
    
    const real w = -weight;
    const real ca1w = ca1 * w, ca2w = ca2 * w;
    const real cb1w = cb1 * w, cb2w = cb2 * w;
   
    const real ee1w = ee1 * w, ee1ee1w = ee1 * ee1w;
    const real ee2w = ee2 * w, ee2ee2w = ee2 * ee2w;
    const real ee1ee2w = ee1 * ee2w;
    
    //we put the isotropic terms in mB
    mB(ia1, ia1) +=  ca1w * ca1 + ee1ee1w;
    mB(ia1, ia2) +=  ca1w * ca2 - ee1ee1w;
    mB(ia2, ia2) +=  ca2w * ca2 + ee1ee1w;
    
    mB(ib1, ib1) +=  cb1w * cb1 + ee2ee2w;
    mB(ib1, ib2) +=  cb1w * cb2 - ee2ee2w;
    mB(ib2, ib2) +=  cb2w * cb2 + ee2ee2w;
    
    mB(ia1, ib1) +=  ca1w * cb1 - ee1ee2w;
    mB(ia1, ib2) +=  ca1w * cb2 + ee1ee2w;
    mB(ia2, ib1) +=  ca2w * cb1 + ee1ee2w;
    mB(ia2, ib2) +=  ca2w * cb2 - ee1ee2w;
    
    //index in the matrix mC:
    ia1 *= DIM;
    ia2 *= DIM;
    ib1 *= DIM;
    ib2 *= DIM;
    
    mC(ia1  , ia2+1) -= ee1w;
    mC(ia1+1, ia2  ) += ee1w;
    
    mC(ib1  , ib2+1) -= ee2w;
    mC(ib1+1, ib2  ) += ee2w;
    
    const real ee1cb1w = ee1w * cb1;
    const real ee1cb2w = ee1w * cb2;
    const real ee2ca1w = ee2w * ca1;
    const real ee2ca2w = ee2w * ca2;
    
    mC(ia1, ib1+1) -=  ee2ca1w + ee1cb1w;
    mC(ia1, ib2+1) +=  ee2ca1w - ee1cb2w;
    
    mC(ia1+1, ib1) +=  ee2ca1w + ee1cb1w;
    mC(ia1+1, ib2) -=  ee2ca1w - ee1cb2w;
    
    mC(ia2, ib1+1) -=  ee2ca2w - ee1cb1w;
    mC(ia2, ib2+1) +=  ee2ca2w + ee1cb2w;
    
    mC(ia2+1, ib1) +=  ee2ca2w - ee1cb1w;
    mC(ia2+1, ib2) -=  ee2ca2w + ee1cb2w;
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(), pta.pos()+cross(ee1,pta.diff()),
                           ptb.pos()+cross(ee2,ptb.diff()), ptb.pos());
    }
#endif
    
    if ( modulo )
        throw Exception("interSideSideLink2D is not usable with periodic boundary conditions");
}

#endif


/**
 Update Meca to include an interaction between A and B,
 but the links are maded between SA and SB which are located
 on the side of A and B, respectively:
 @code
 SA = A + len * N_A,
 SB = B + len * N_B,
 @endcode
 N_X is a normalized vector orthogonal to the fiber carrying X, in X:
 The force is linear of zero resting length,
 @code
 force_SA = weight * ( SA - SB )
 force_SB = weight * ( SB - SA )
 @endcode
 */

void Meca::interSideSideLink(const PointInterpolated & pta,
                             const PointInterpolated & ptb,
                             const real len,
                             const real weight )
{
#if ( DIM == 1 )
    
    throw Exception("Meca::interSideSideLink meaningless in 1D");
    
#elif ( DIM == 2 )
    
    Vector dir = ptb.pos() - pta.pos();
    int side1 = RNG.sign_exc( cross(pta.diff(), dir) );
    int side2 = RNG.sign_exc( cross(dir, ptb.diff()) );
    interSideSideLink2D(pta, ptb, len, weight, side1, side2);
    
#elif ( DIM == 3 )
    
    throw Exception("Meca::interSideSideLink was not implemented in 3D");
    
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Sliding Links
//------------------------------------------------------------------------------

/**
 Update Meca to include an interaction between A and B,
 The force is linear of zero resting length, but is anisotropic:
 The component of the force parallel to the fiber in A is removed

 If T is the normalized direction of the fiber in A:
 @code
 force_A = weight * ( 1 - T T' ) ( A - B )
 force_B = weight * ( 1 - T T' ) ( B - A )
 @endcode
 */

void Meca::interSlidingLink(const PointInterpolated & pta,
                            const PointExact & ptb,
                            const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const index_type ia1 = DIM * pta.matIndex1();
    const index_type ia2 = DIM * pta.matIndex2();
    const index_type ib  = DIM * ptb.matIndex();
    
    if ( any_equal(ia1, ia2, ib) )
        return;
    
    //force coefficients on the points:
    const real A = pta.coef2();
    const real B = pta.coef1();
    const real AA = A * A, AB = A * B, BB = B * B;
    
    Vector dir = pta.dir();
    
    // on points (a, b, e), (ab) being the PointInterpolated, and e the PointExact,
    // P is the projection on the plane perpendicular to (ab): P.v= (v - (T.v)T/normSqr(T))
    // the interaction is  -weigth * transpose(bb, aa, -1) * P * ( bb, aa, -1 )
    // we set only the upper part of this symmetric matrix:
    
    for ( int dd = 0; dd < DIM; ++dd )
    {
        for ( int ee = dd; ee < DIM; ++ee )
        {
            real P = weight * ((dd==ee) - dir[dd]*dir[ee] );
            
            mC(ia1+dd, ia1+ee) -= P * AA;
            mC(ia1+dd, ia2+ee) -= P * AB;
            mC(ia2+dd, ia2+ee) -= P * BB;
            mC(ia1+dd, ib +ee) += P * A;
            mC(ia2+dd, ib +ee) += P * B;
            mC(ib +dd, ib +ee) -= P;
            if ( dd != ee )
            {
                mC(ia1+ee, ia2+dd) -= P * AB;
                mC(ia1+ee, ib +dd) += P * A;
                mC(ia2+ee, ib +dd) += P * B;
            }
        }
    }
    
    if ( modulo )
    {
        Vector off = modulo->offset( pta.pos() - ptb.pos() );
        if ( !off.null() )
        {
            real s = dir * off; //scalar product
            Vector v = weight * off - ( weight * s ) * dir;
            v.sub_to(vBAS+ib);
            v.add_to(A, vBAS+ia1);
            v.add_to(B, vBAS+ia2);
        }
    }
}


/**
Update Meca to include an interaction between A and B,
 The force is linear of zero resting length, but is anisotropic:
 The component of the force parallel to the fiber in A is removed
 
 If T is the normalized direction of the fiber in A:
 @code
 force_A = weight * ( 1 - T T' ) ( A - B )
 force_B = weight * ( 1 - T T' ) ( B - A )
 @endcode
 */

void Meca::interSlidingLink(const PointInterpolated & pta,
                            const PointInterpolated & ptb, 
                            const real weight)
{
    assert_true( weight >= 0 );
    
    //the index of the points in the matrix mB:
    const index_type inx[] = { DIM*pta.matIndex1(), DIM*pta.matIndex2(),
                               DIM*ptb.matIndex1(), DIM*ptb.matIndex2() };
    
    if ( any_equal(inx[0], inx[1], inx[2], inx[3]) )
        return;
    
    //interpolation coefficients
    const real  c[4] = { pta.coef2(), pta.coef1(), -ptb.coef2(), -ptb.coef1() };
    const real cw[4] = { -weight*c[0], -weight*c[1], -weight*c[2], -weight*c[3] };
    
    // on points (a, b, e), (ab) being the PointInterpolated, and e the PointExact,
    // P is the projection on the plane perpendicular to (ab): P.v= (v - (T.v)T/normSqr(T))
    // the interaction is  -wh' * P * h
    // we set only the upper part of this symmetric matrix:
    
    Vector dir = pta.dir();

#if ( DIM == 1 )
    
    throw Exception("Meca::interSlidingLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    real dd[] = { dir.XX*dir.XX, dir.XX*dir.YY, dir.YY*dir.YY };

    for ( int jj = 0; jj < 4; ++jj )
    {
        const real g = c[jj] * cw[jj];
        mC(inx[jj],   inx[jj]  ) += ( 1.0 - dd[0] ) * g;
        mC(inx[jj]+1, inx[jj]+1) += ( 1.0 - dd[2] ) * g;
        mC(inx[jj],   inx[jj]+1) += (     - dd[1] ) * g;
        
        for ( int ii = jj+1; ii < 4; ++ii )
        {
            const real h = c[jj] * cw[ii];
            mC(inx[ii],   inx[jj]  ) += ( 1.0 - dd[0] ) * h;
            mC(inx[ii]+1, inx[jj]+1) += ( 1.0 - dd[2] ) * h;
            mC(inx[ii],   inx[jj]+1) += (     - dd[1] ) * h;
            mC(inx[ii]+1, inx[jj]  ) += (     - dd[1] ) * h;
        }
    }
    
#elif ( DIM == 3 )
    
    real dd[] = { dir.XX*dir.XX, dir.XX*dir.YY, dir.XX*dir.ZZ, dir.YY*dir.YY, dir.YY*dir.ZZ, dir.ZZ*dir.ZZ };

    for ( int jj = 0; jj < 4; ++jj )
    {
        const real g = c[jj] * cw[jj];
        mC(inx[jj],   inx[jj]  ) += ( 1.0 - dd[0] ) * g;
        mC(inx[jj]+1, inx[jj]+1) += ( 1.0 - dd[3] ) * g;
        mC(inx[jj]+2, inx[jj]+2) += ( 1.0 - dd[5] ) * g;
        mC(inx[jj],   inx[jj]+1) += (     - dd[1] ) * g;
        mC(inx[jj],   inx[jj]+2) += (     - dd[2] ) * g;
        mC(inx[jj]+1, inx[jj]+2) += (     - dd[4] ) * g;
        
        for ( int ii = jj+1; ii < 4; ++ii )
        {
            const real h = c[jj] * cw[ii];
            mC(inx[ii],   inx[jj]  ) += ( 1.0 - dd[0] ) * h;
            mC(inx[ii]+1, inx[jj]+1) += ( 1.0 - dd[3] ) * h;
            mC(inx[ii]+2, inx[jj]+2) += ( 1.0 - dd[5] ) * h;
            mC(inx[ii],   inx[jj]+1) += (     - dd[1] ) * h;
            mC(inx[ii],   inx[jj]+2) += (     - dd[2] ) * h;
            mC(inx[ii]+1, inx[jj]+2) += (     - dd[4] ) * h;
            mC(inx[ii]+1, inx[jj]  ) += (     - dd[1] ) * h;
            mC(inx[ii]+2, inx[jj]  ) += (     - dd[2] ) * h;
            mC(inx[ii]+2, inx[jj]+1) += (     - dd[4] ) * h;
        }
    }
#endif
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptb.pos() - pta.pos() );
        if ( !off.null() )
        {
            Vector vec = off - ( dir * off ) * dir;
            for ( int ii = 0; ii < 4; ++ii )
                vec.add_to(cw[ii], vBAS+inx[ii]);
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Side Sliding Links
//------------------------------------------------------------------------------

#if ( DIM == 2 )

void mat4print(const real M[4], const real alpha)
{
    fprintf(stderr, "%9.3f %9.3f\n", alpha*M[0], alpha*M[2]);
    fprintf(stderr, "%9.3f %9.3f\n", alpha*M[1], alpha*M[3]);
}
///calculates M = A*B
void mat4Mul(real M[4], const real A[4], const real B[4])
{
    M[0] = A[0]*B[0] + A[2]*B[1];
    M[1] = A[1]*B[0] + A[3]*B[1];
    M[2] = A[0]*B[2] + A[2]*B[3];
    M[3] = A[1]*B[2] + A[3]*B[3];
}
///calculates M = At*B
void mat4MulT(real M[4], const real A[4], const real B[4])
{
    M[0] = A[0]*B[0] + A[1]*B[1];
    M[1] = A[2]*B[0] + A[3]*B[1];
    M[2] = A[0]*B[2] + A[1]*B[3];
    M[3] = A[2]*B[2] + A[3]*B[3];
}


void Meca::interSideSlidingLink2D(const PointInterpolated & pta,
                                  const PointExact & pte, 
                                  const real arm,
                                  const real weight)
{
    assert_true( weight >= 0 );

    // indices
    const index_type inx[] = { DIM*pta.matIndex1(), DIM*pta.matIndex2(), DIM*pte.matIndex() };
    
    if ( any_equal(inx[0], inx[1], inx[2]) )
        return;
    
    Vector dir = pta.dir();
    const real aa = pta.coef2();
    const real bb = pta.coef1();
    const real ee = arm / pta.len();
    
    // the projection matrix (symmetric):
    const real P[4] = { 1.0-dir.XX*dir.XX,    -dir.XX*dir.YY,
                           -dir.XX*dir.YY, 1.0-dir.YY*dir.YY };

#ifdef USE_MATRIX_BLOCK
    /*
     We use block operations to set the matrix block by block:
     | A'PA  A'PB  A'P |
     | B'PA  B'PB  B'P |
     |   PA    PB    P |
     */
    
    // anti-symmetric matrix blocks:
    const real A[4] = { -aa,  ee, -ee, -aa };
    const real B[4] = { -bb, -ee,  ee, -bb };
    
    if ( inx[2] > inx[1] )
    {
        // in this case, we directly set the lower side
        real PA[4], PB[4], M[4];

        mat4Mul(PA,P,A);
        mat4Mul(PB,P,B);
    
        mat4MulT(M, A, PA);
        mC.block(inx[0], inx[0]).add_lower(-weight, M);
    
        mat4MulT(M, B, PA);
        mC.block(inx[1], inx[0]).add(-weight, M);
    
        mat4MulT(M, B, PB);
        mC.block(inx[1], inx[1]).add_lower(-weight, M);
        
        mC.block(inx[2], inx[0]).add(-weight, PA);
        mC.block(inx[2], inx[1]).add(-weight, PB);
    }
    else
    {
        real AP[4], BP[4], M[4];

        mat4MulT(AP,A,P);
        mat4MulT(BP,B,P);
        
        mat4Mul(M, AP, A);
        mC.block(inx[0], inx[0]).add_lower(-weight, M);
        
        mat4Mul(M, BP, A);
        mC.block(inx[1], inx[0]).add(-weight, M);
        
        mat4Mul(M, BP, B);
        mC.block(inx[1], inx[1]).add_lower(-weight, M);

        // in this case, the two blocks are transposed
        mC.block(inx[0], inx[2]).add(-weight, AP);
        mC.block(inx[1], inx[2]).add(-weight, BP);
    }

    mC.block(inx[2], inx[2]).add_lower(-weight, P);
    
    if ( modulo )
    {
        Vector off = modulo->offset( pta.pos() - pte.pos() );
        if ( !off.null() )
        {
            Vector vec = weight * ( ( off * dir ) * dir - off );
            vBAS[inx[0]  ] += A[0] * vec.XX + A[1] * vec.YY;
            vBAS[inx[0]+1] += A[2] * vec.XX + A[3] * vec.YY;
            vBAS[inx[1]  ] += B[0] * vec.XX + B[1] * vec.YY;
            vBAS[inx[1]+1] += B[2] * vec.XX + B[3] * vec.YY;
            vBAS[inx[2]  ] += vec.XX;
            vBAS[inx[2]+1] += vec.YY;
        }
    }
#else
    
    // matrix of coefficients 2x6:
    real T[2*6] = { aa, -ee, ee, aa, bb, ee,
                   -ee,  bb, -1,  0,  0, -1 };

    // here we brutally multiply matrices
    real PT[2*6], TPT[6*6];
    blas_xgemm('N','N', 2, 6, 2, -weight, P, 2, T, 2, 0.0, PT, 2);
    blas_xgemm('T','N', 6, 6, 2, 1.0, T, 2, PT, 2, 0.0, TPT, 6);
    
    //printf("\n"); VecPrint::matPrint(std::clog, 6, 6, TPT, 6);
    
    const index_type ixx[] = { inx[0], inx[0]+1, inx[1], inx[1]+1, inx[2], inx[2]+1 };
    
    for ( int x=0; x<6; ++x )
    for ( int y=x; y<6; ++y )
        mC(ixx[y], ixx[x]) += TPT[y+6*x];
    
    if ( modulo )
    {
        Vector off = modulo->offset( pta.pos() - pte.pos() );
        if ( !off.null() )
        {
            for ( int ii=0; ii<6; ++ii )
                vBAS[ixx[ii]] += TPT[ii+6*4] * off.XX + TPT[ii+6*5] * off.YY;
        }
     }
#endif
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(),pta.pos()+cross(arm, pta.dir()),pte.pos());
    }
#endif
}



/**
 Alternative method in which we add an offset to vBAS 
 */
void Meca::interSideSlidingLinkS(const PointInterpolated & pta,
                                 const PointExact & pte,
                                 const real arm,
                                 const real weight)
{
    assert_true( weight >= 0 );
    
    // indices
    const index_type inx[] = { DIM*pta.matIndex1(), DIM*pta.matIndex2(), DIM*pte.matIndex() };
    
    if ( any_equal(inx[0], inx[1], inx[2]) )
        return;
    
    // choose vector 'a' parallel to the first Fiber:
    Vector a = pta.dir();
    
    const int side = RNG.sign_exc(arm);
    const real bXX = -side*a.YY, bYY = side*a.XX;
    const real T[4] = { bXX*bXX, bXX*bYY, bXX*bYY, bYY*bYY };
    
    // we set directly the transformed offset vector:
    const real RB[2] = { -arm*a.YY, arm*a.XX };
    
    // weights and indices:
    const real  c[3] = { pta.coef2(), pta.coef1(), -1.0 };
    const real cw[3] = { -weight*c[0], -weight*c[1], -weight*c[2] };
    
    // fill the matrix mC
    for ( int ii=0; ii<3; ++ii )
    {
        vBAS[inx[ii]  ] += cw[ii] * RB[0];
        vBAS[inx[ii]+1] += cw[ii] * RB[1];
        
        const real g = c[ii] * cw[ii];
        mC(inx[ii]  , inx[ii]  ) += g * T[0];
        mC(inx[ii]  , inx[ii]+1) += g * T[2];
        mC(inx[ii]+1, inx[ii]+1) += g * T[3];
        
        for ( int jj=ii+1; jj<3; ++jj )
        {
            const real h = c[ii] * cw[jj];
            mC(inx[ii]  , inx[jj]  ) += h * T[0];
            mC(inx[ii]  , inx[jj]+1) += h * T[2];
            mC(inx[ii]+1, inx[jj]  ) += h * T[1];
            mC(inx[ii]+1, inx[jj]+1) += h * T[3];
        }
    }
    
    if ( modulo )
        throw Exception("interSideSlidingLinkS is not usable with periodic boundary conditions");    
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(),pta.pos()+cross(arm, pta.dir()),pte.pos());
    }
#endif
}


#elif ( DIM == 3 )

/**
 Vector 'arm' must be parallel to the link and orthogonal to 'pta'
 */

void Meca::interSideSlidingLinkS(const PointInterpolated & pta,
                                 const PointExact & pte,
                                 const Vector & dir,
                                 const real len,
                                 const real weight)
{    
    assert_true( weight >= 0 );
    assert_true( len > REAL_EPSILON );
    
    const index_type inx[] = { DIM*pta.matIndex1(), DIM*pta.matIndex2(),
                               DIM*pte.matIndex() };

    if ( any_equal(inx[0], inx[1], inx[2]) )
        return;

    // choose vector 'a' parallel to the first Fiber:
    const Vector a = pta.dir();
    // vector 'b' aligned with the link:
    const Vector b = dir;
    
    /*      
     Without tangential force, a 'long link' is in the perpendicular direction.
     In the local reference frame, the matrix of interaction coefficients would be:
     real T[9] = { 0, 0, 0, 0, -weight, 0, 0, 0, 0 };
     we could transform it with a change-of-coordinates matrix R:
     Vector c = cross(a, b);
     real R[9] = { a.XX, a.YY, a.ZZ, b.XX, b.YY, b.ZZ, c.XX, c.YY, c.ZZ };
     real TR[3*3];
     blas_xgemm('N','T', 3, 3, 3, 1.0, T, 3, R, 3, 0.0, TR, 3);
     blas_xgemm('N','N', 3, 3, 3, 1.0, R, 3, TR, 3, 0.0, T, 3);
     equivalently, we can set directly the interaction coefficient matrix: 
     */
    
    const real xx=b.XX*b.XX, xy=b.XX*b.YY, xz=b.XX*b.ZZ;
    const real yy=b.YY*b.YY, yz=b.YY*b.ZZ, zz=b.ZZ*b.ZZ;
    const real T[9] = { xx, xy, xz, xy, yy, yz, xz, yz, zz };
    
    // we set directly the transformed offset vector:
    const real RB[3] = { len*b.XX, len*b.YY, len*b.ZZ };
    
    // weights and indices:
    const real  c[3] = { pta.coef2(), pta.coef1(), -1.0 };
    const real cw[3] = { -weight*c[0], -weight*c[1], -weight*c[2] };
    
    // fill the matrix mC
    for ( int ii=0; ii<3; ++ii )
    {
        for ( int x = 0; x < DIM; ++x )
            vBAS[inx[ii]+x] += cw[ii] * RB[x];
        
        const real g = c[ii] * cw[ii];
#ifdef USE_MATRIX_BLOCK
        mC.block(inx[ii], inx[ii]).add_lower(g, T);
#else
        for ( int x = 0; x < DIM; ++x )
        for ( int y = x; y < DIM; ++y )
            mC(inx[ii]+y, inx[ii]+x) += g * T[y+3*x];
#endif

        for ( int jj=ii+1; jj<3; ++jj )
        {
            const real h = c[ii] * cw[jj];
#ifdef USE_MATRIX_BLOCK
            mC.block(inx[ii], inx[jj]).add(h, T);
#else
            for ( int x = 0; x < DIM; ++x )
            for ( int y = 0; y < DIM; ++y )
                mC(inx[ii]+y, inx[jj]+x) += h * T[y+3*x];
#endif
        }
    }
        
    if ( modulo )
    {
        Vector off = modulo->offset( pte.pos() - pta.pos() );
        if ( !off.null() )
        {
            Vector vec = ( b * off ) * b;
            for ( int ii = 0; ii < 3; ++ii )
                vec.add_to(cw[ii], vBAS+inx[ii]);
        }
    }
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(),pta.pos()+Vector(RB),pte.pos());
    }
#endif
    
}
#endif

/**
 Update Meca to include an interaction between A and B,
 This is a combination of a SideLink with a Sliding Link:
 The force is linear of zero resting length, but it is taken between B,
 and another point S located on the side of A:
 S = A + len * N,
 where N is a normalized vector orthogonal to the fiber in A, in the direction of B.
 In addition, the tangential part of the force is removed.
 
 If T is the normalized direction of the fiber in A:
 @code
 force_S = weight * ( 1 - T T' ) ( S - B )
 force_B = weight * ( 1 - T T' ) ( B - S )
 @endcode
 */
void Meca::interSideSlidingLink(const PointInterpolated & pta, 
                                const PointExact & ptb, 
                                const real len,
                                const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::interSideLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    Vector as = ptb.pos()-pta.pos();
    if ( modulo )
        modulo->fold(as);
    real arm  = len * RNG.sign_exc( cross(pta.diff(), as) );
    interSideSlidingLink2D(pta, ptb, arm, weight);
    
#elif ( DIM == 3 )
    
    // 'arm' is perpendicular to Fiber and parallel to link:
    Vector a   = pta.diff();
    Vector as  = ptb.pos() - pta.pos();
    if ( modulo ) 
        modulo->fold(as);
    Vector arm = as - ( ( as * a ) / a.normSqr() ) * a;
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        interSideSlidingLinkS(pta, ptb, arm / n, len, weight);
    
#endif
}



#if ( DIM == 2 )


void Meca::interSideSlidingLink2D(const PointInterpolated & pta,
                                  const PointInterpolated & ptb, 
                                  const real arm,
                                  const real weight)
{
    assert_true( weight >= 0 );

    const index_type inx[] = { DIM*pta.matIndex1(),  DIM*pta.matIndex1()+1,
                               DIM*pta.matIndex2(),  DIM*pta.matIndex2()+1,
                               DIM*ptb.matIndex1(),  DIM*ptb.matIndex1()+1,
                               DIM*ptb.matIndex2(),  DIM*ptb.matIndex2()+1 };

    if ( any_equal(inx[0], inx[2], inx[4], inx[6]) )
        return;

    Vector dir = pta.dir();
    const real A1 =  pta.coef2(), A2 =  pta.coef1();
    const real B1 = -ptb.coef2(), B2 = -ptb.coef1();
    
    const real ee = arm / pta.len();

    //this is done the 'hard' way by multiplying all matrices
    //coefficient matrix:
    real T[2*8] = { A1, -ee, ee, A1, A2, ee, -ee,  A2,
                    B1,   0,  0, B1, B2,  0,   0,  B2 };
    
    //the projection matrix:
    const real P[4] = { 1-dir.XX*dir.XX, -dir.XX*dir.YY, -dir.XX*dir.YY, 1-dir.YY*dir.YY };
    
    real PT[2*8], TPT[8*8];
    blas_xgemm('N','N', 2, 8, 2, -weight, P, 2, T, 2, 0.0, PT, 2);
    blas_xgemm('T','N', 8, 8, 2, 1.0, T, 2, PT, 2, 0.0, TPT, 8);
    
    //printf("\n");  VecPrint::matPrint(8,8, TPT);
    
    for ( int ii=0; ii<8; ++ii )
        for ( int jj=ii; jj<8; ++jj )
            mC(inx[ii], inx[jj]) += TPT[ii+8*jj];
    
    if ( modulo )
    {
        Vector off = modulo->offset( pta.pos() - ptb.pos() );
        if ( !off.null() )
        {
            for ( int ii=0; ii<8; ++ii )
            {
                vBAS[inx[ii]] += TPT[ii+8*4] * off.XX;
                vBAS[inx[ii]] += TPT[ii+8*5] * off.YY;
                vBAS[inx[ii]] += TPT[ii+8*6] * off.XX;
                vBAS[inx[ii]] += TPT[ii+8*7] * off.YY;
            }
        }
    }
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(),pta.pos()+cross(arm, pta.dir()),ptb.pos());
    }
#endif
}

    
void Meca::interSideSlidingLinkS(const PointInterpolated & pta,
                                 const PointInterpolated & ptb,
                                 const real arm,
                                 const real weight)
{
    assert_true( weight >= 0 );
    
    const index_type inx[] = { DIM*pta.matIndex1(), DIM*pta.matIndex2(),
                               DIM*ptb.matIndex1(), DIM*ptb.matIndex2() };
    
    if ( any_equal(inx[0], inx[1], inx[2], inx[3]) )
        return;

    // vector 'a' parallel to the first Fiber:
    Vector a = pta.dir();
    
    // vector 'b' perpendicular to 'a', and aligned with the link:
    const int side = RNG.sign_exc(arm);
    const real bXX = -side*a.YY, bYY = side*a.XX;
    const real T[4] = { bXX*bXX, bXX*bYY, bXX*bYY, bYY*bYY };
    
    // we set directly the transformed offset vector:
    const real RB[2] = { -arm*a.YY, arm*a.XX };

    // weights and indices:
    const real  c[4] = { pta.coef2(), pta.coef1(), -ptb.coef2(), -ptb.coef1() };
    const real cw[4] = { -weight*c[0], -weight*c[1], -weight*c[2], -weight*c[3] };
    
    // fill the matrix mC
    for ( int ii=0; ii<4; ++ii )
    {
        vBAS[inx[ii]  ] += cw[ii] * RB[0];
        vBAS[inx[ii]+1] += cw[ii] * RB[1];
        
        const real g = c[ii] * cw[ii];
        mC(inx[ii]  , inx[ii]  ) += g * T[0];
        mC(inx[ii]  , inx[ii]+1) += g * T[2];
        mC(inx[ii]+1, inx[ii]+1) += g * T[3];
        
        for ( int jj=ii+1; jj<4; ++jj )
        {
            const real h = c[ii] * cw[jj];
            mC(inx[ii]  , inx[jj]  ) += h * T[0];
            mC(inx[ii]  , inx[jj]+1) += h * T[2];
            mC(inx[ii]+1, inx[jj]  ) += h * T[1];
            mC(inx[ii]+1, inx[jj]+1) += h * T[3];
        }
    }
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptb.pos() - pta.pos() );
        if ( !off.null() )
        {
            Vector vec = off - ( a * off ) * a;
            for ( int ii = 0; ii < 4; ++ii )
                vec.add_to(cw[ii], vBAS+inx[ii]);
        }
    }
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(),pta.pos()+cross(arm, pta.dir()),ptb.pos());
    }
#endif
}


#elif ( DIM == 3 )

    /**
     Vector 'arm' must be parallel to the link and orthogonal to 'pta'
     */
void Meca::interSideSlidingLinkS(const PointInterpolated & pta,
                                 const PointInterpolated & ptb, 
                                 const Vector & dir,
                                 const real arm,
                                 const real weight)
{
    assert_true( weight >= 0 );
    assert_true( arm > REAL_EPSILON );
    
    const index_type inx[] = { DIM*pta.matIndex1(), DIM*pta.matIndex2(),
                               DIM*ptb.matIndex1(), DIM*ptb.matIndex2() };

    if ( any_equal(inx[0], inx[1], inx[2], inx[3]) )
        return;

    // choose vector 'a' parallel to the first Fiber:
    const Vector a = pta.dir();
    // make vector 'b' perpendicular to 'a', and aligned with the link:
    const Vector b = dir;
    
    /*
     Without tangential force, a 'long link' is in the perpendicular direction.
     In the local reference frame, the matrix of interaction coefficients would be:
     real T[9] = { 0, 0, 0, 0, -weight, 0, 0, 0, 0 };
     we could transform it with a change-of-coordinates matrix R:
     Vector c = cross(a, b);
     real R[9] = { a.XX, a.YY, a.ZZ, b.XX, b.YY, b.ZZ, c.XX, c.YY, c.ZZ };
     real TR[3*3];
     blas_xgemm('N','T', 3, 3, 3, 1.0, T, 3, R, 3, 0.0, TR, 3);
     blas_xgemm('N','N', 3, 3, 3, 1.0, R, 3, TR, 3, 0.0, T, 3);
     equivalently, we can set directly the interaction coefficient matrix: 
     */
    
    const real xx = b.XX*b.XX, xy = b.XX*b.YY, xz = b.XX*b.ZZ;
    const real yy = b.YY*b.YY, yz = b.YY*b.ZZ, zz = b.ZZ*b.ZZ;
    const real T[9] = { xx, xy, xz, xy, yy, yz, xz, yz, zz };
    
    // we set directly the transformed offset vector:
    const real RB[3] = { arm*b.XX, arm*b.YY, arm*b.ZZ };
    
    // weights and indices:
    const real  c[4] = { pta.coef2(), pta.coef1(), -ptb.coef2(), -ptb.coef1() };
    const real cw[4] = { -weight*c[0], -weight*c[1], -weight*c[2], -weight*c[3] };
    
    // fill the matrix mC
    for ( int ii=0; ii<4; ++ii )
    {
        for ( int x = 0; x < DIM; ++x )
            vBAS[inx[ii]+x] += cw[ii] * RB[x];
        
        const real g = c[ii] * cw[ii];
#ifdef USE_MATRIX_BLOCK
        mC.block(inx[ii], inx[ii]).add_lower(g, T);
#else
        for ( int x = 0; x < DIM; ++x )
        for ( int y = x; y < DIM; ++y )
            mC(inx[ii]+y, inx[ii]+x) += g * T[y+3*x];
#endif

        for ( int jj=ii+1; jj<4; ++jj )
        {
            const real h = c[ii] * cw[jj];
#ifdef USE_MATRIX_BLOCK
            mC.block(inx[ii], inx[jj]).add(h, T);
#else
            for ( int x = 0; x < DIM; ++x )
            for ( int y = 0; y < DIM; ++y )
                mC(inx[ii]+y, inx[jj]+x) += h * T[y+3*x];
#endif

        }
    }

    if ( modulo )
    {
        Vector off = modulo->offset( ptb.pos() - pta.pos() );
        if ( !off.null() )
        {
            Vector vec = ( b * off ) * b;
            for ( int ii = 0; ii < 4; ++ii )
                vec.add_to(cw[ii], vBAS+inx[ii]);
         }
    }
    
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        //gle::bright_color(pta.mecable()->signature()).load();
        glColor3f(0,0,1);
        displayInteraction(pta.pos(),pta.pos()+Vector(RB),ptb.pos());
    }
#endif
    
}

#endif

/**
 Update Meca to include an interaction between A and B,
 This is a combination of Side- and Sliding Links:
 The force is linear of zero resting length, but it is taken between B
 and another point S which is located on the side of A:
 S = A + len * N,
 where N is a normalized vector orthogonal to the fiber in A, in the direction of B.
 In addition, the part of the force tangential to A is removed.
 
 If T is the normalized direction of the fiber in A:
 @code
 force_S = weight * ( 1 - T T' ) ( S - B )
 force_B = weight * ( 1 - T T' ) ( B - S )
 @endcode
 */

void Meca::interSideSlidingLink(const PointInterpolated & pta, 
                                const PointInterpolated & ptb, 
                                const real len,
                                const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::interSideSlidingLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    real arm = len * RNG.sign_exc( cross(pta.diff(), ptb.pos()-pta.pos()) );
    interSideSlidingLink2D(pta, ptb, arm, weight);
    
#elif ( DIM == 3 )
    
    // 'arm' is perpendicular to A-Fiber and parallel to link:
    Vector a   = pta.diff();
    Vector as  = ptb.pos() - pta.pos();
    if ( modulo ) 
        modulo->fold(as);
    Vector arm = as - ( ( as * a ) / a.normSqr() ) * a;
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        interSideSlidingLinkS(pta, ptb, arm / n, len, weight);
    
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Point Clamps
//------------------------------------------------------------------------------

/**
 Update Meca to include a link between a point A and a fixed position G
 The force is linear:
 @code
 force_A = weight * ( G - A );
 @endcode
 There is no counter-force in G, since G is immobile.
 */

void Meca::addPointClamp(PointExact const& pta,
                         Vector pos,
                         const real weight)
{
    assert_true( weight >= 0 );
    const index_type inx = pta.matIndex();
    
    mB(inx, inx) -= weight;
    
    if ( modulo )
        modulo->fold(pos, pta.pos());

    pos.add_to(weight, vBAS+inx);
}



/**
 Update Meca to include a link between a point A and a fixed position G
 The force is linear:  
 @code
 force_A = weight * ( G - A );
 @endcode
 The point G is not associated to a Mecable, and there is no counter-force in G.
 */

void Meca::addPointClamp(PointInterpolated const& pti,
                         Vector pos,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    index_type inx1 = pti.matIndex1();
    index_type inx2 = pti.matIndex2();
    
    const real c1 = pti.coef2();
    const real c2 = pti.coef1();
    
    assert_true( inx1 != inx2 );
    assert_true( 0 <= c1  &&  c1 <= 1 );
    assert_true( 0 <= c2  &&  c2 <= 1 );

    const real c2w = weight * c2;
    const real c1w = weight * c1;
    
    mB(inx1, inx1) -=  c1w * c1;
    mB(inx1, inx2) -=  c2w * c1;
    mB(inx2, inx2) -=  c2w * c2;
    
    inx1 *= DIM;
    inx2 *= DIM;
    
    if ( modulo )
        modulo->fold(pos, pti.pos());
    
    pos.add_to(c1w, vBAS+inx1);
    pos.add_to(c2w, vBAS+inx2);
}


//------------------------------------------------------------------------------
#pragma mark - Long Point Clamps
//------------------------------------------------------------------------------

/**
Update Meca to include a non-zero resting force between a point `P` and a position `G`
 The force is affine with non-zero resting length: 
 @code
 force = weight * ( G - P ) * ( len / |PG| - 1 )
 @endcode
 There is no force on G, which is an immobile position.
 */

void Meca::addLongPointClamp(Vector const& pos,
                             PointExact const& pte,
                             Vector const& center,
                             real len,
                             const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );
    const index_type inx = DIM * pte.matIndex();
    
    Vector axis = pos - center;
    real axis_n = axis.norm();
    
    if ( axis_n > REAL_EPSILON )
        axis /= axis_n;
    else
    {
        axis   = Vector::randU(len/2);
        axis_n = len/2;
    }
    
    if ( len < axis_n )
    {
#if ( DIM > 1 )
        len /= axis_n;
        real facX = weight * len * ( axis_n + axis * center );
        real facC = weight * ( 1.0 - len );
        
        mC(inx  , inx  ) += weight * ( len * ( 1.0 - axis.XX * axis.XX ) - 1.0 );
        mC(inx+1, inx  ) -= weight * len * axis.XX * axis.YY;
        mC(inx+1, inx+1) += weight * ( len * ( 1.0 - axis.YY * axis.YY ) - 1.0 );
        
        vBAS[inx  ] += facX * axis.XX + facC * center.XX;
        vBAS[inx+1] += facX * axis.YY + facC * center.YY;
#endif
#if ( DIM == 3 )
        mC(inx+2, inx  ) -= weight * len * axis.XX * axis.ZZ;
        mC(inx+2, inx+1) -= weight * len * axis.YY * axis.ZZ;
        mC(inx+2, inx+2) += weight * ( len * ( 1.0 - axis.ZZ * axis.ZZ ) - 1.0 );

        vBAS[inx+2] += facX * axis.ZZ + facC * center.ZZ;
#endif
    }
    else
    {
#if ( DIM > 1 )
        real facX = weight * ( len + axis * center );

        mC(inx  , inx  ) -= weight * axis.XX * axis.XX;
        mC(inx+1, inx  ) -= weight * axis.XX * axis.YY;
        mC(inx+1, inx+1) -= weight * axis.YY * axis.YY;
        
        vBAS[inx  ] += facX * axis.XX;
        vBAS[inx+1] += facX * axis.YY;
#endif
#if ( DIM == 3 )
        mC(inx+2, inx  ) -= weight * axis.XX * axis.ZZ;
        mC(inx+2, inx+1) -= weight * axis.YY * axis.ZZ;
        mC(inx+2, inx+2) -= weight * axis.ZZ * axis.ZZ;
        
        vBAS[inx+2] += facX * axis.ZZ;
#endif
    }
    
    if ( modulo )
        throw Exception("addLongPointClamp is not usable with periodic boundary conditions");
}


void Meca::addLongPointClamp(PointExact const& pte,
                          Vector const& center,
                          real len,
                          const real weight)
{
    addLongPointClamp(pte.pos(), pte, center, len, weight);
}


void Meca::addLongPointClamp(PointInterpolated const& pti,
                          Vector const& center,
                          real len,
                          const real weight)
{
    // interpolate on the two flanking model points using coefficients:
    addLongPointClamp(pti.pos(), pti.exact1(), center, len, weight*pti.coef2());
    addLongPointClamp(pti.pos(), pti.exact2(), center, len, weight*pti.coef1());
}




/**
 Update Meca to include a non-zero resting force between a point `P` and a position `G`
 The force is affine with non-zero resting length:
 @code
 force = weight * ( G - P ) * ( len / |PG| - 1 )
 @endcode
 The force is only in the YZ plane.
 */

void Meca::addLongPointClampYZ(const PointExact & pte,
                               real  len,
                               const real weight)
{
    assert_true( weight >= 0 );
    const index_type inx = DIM * pte.matIndex();
    
#if ( DIM == 2 )
    
    if ( pte.pos().YY > 0 )
    {
        mC(inx+1, inx+1) -= weight;
        vBAS[inx+1]      += weight * len;
    }
    else
    {
        mC(inx+1, inx+1) -= weight;
        vBAS[inx+1]      -= weight * len;
    }

#elif ( DIM == 3 )

    Vector pos = pte.pos();
    real axis_n = pos.normYZ();
    if ( axis_n < REAL_EPSILON )
        return;
    
    Vector axis(0, pos.YY/axis_n, pos.ZZ/axis_n);
    
    real facX;
    
    if ( len < axis_n )
    {
        len /= axis_n;
        facX = weight * len * axis_n;
        
        mC(inx+1, inx+1) += weight * ( len * ( 1.0 - axis.YY * axis.YY ) - 1.0 );
        mC(inx+1, inx+2) -= weight * len * axis.YY * axis.ZZ;
        mC(inx+2, inx+2) += weight * ( len * ( 1.0 - axis.ZZ * axis.ZZ ) - 1.0 );
    }
    else
    {
        facX = weight * len;

        mC(inx+1, inx+1) -= weight * axis.YY * axis.YY;
        mC(inx+1, inx+2) -= weight * axis.YY * axis.ZZ;
        mC(inx+2, inx+2) -= weight * axis.ZZ * axis.ZZ;
    }
    
    vBAS[inx+1] += facX * axis[1];
    vBAS[inx+2] += facX * axis[2];
    
#endif
}


/**
 Update Meca to include a non-zero resting force between a point `P` and a position `G`
 The force is affine with non-zero resting length:
 @code
 force = weight * ( G - P ) * ( len / |PG| - 1 )
 @endcode
 The force is constrained in the XY plane.
 */

void Meca::addLongPointClampXY(const PointExact & pte,
                               real  len,
                               const real weight)
{
    assert_true( weight >= 0 );
    
#if ( DIM > 1 )

    const index_type inx = DIM * pte.matIndex();
    Vector pos = pte.pos();
    real axis_n = pos.normXY();
    if ( axis_n < REAL_EPSILON )
        return;
    
    Vector axis(pos.XX/axis_n, pos.YY/axis_n, 0);

    real facX;
    
    if ( len < axis_n )
    {
        len /= axis_n;
        facX = weight * len * axis_n;
        
        mC(inx  , inx  ) += weight * ( len * ( 1.0 - axis.XX * axis.XX ) - 1.0 );
        mC(inx  , inx+1) -= weight * len * axis.XX * axis.YY;
        mC(inx+1, inx+1) += weight * ( len * ( 1.0 - axis.YY * axis.YY ) - 1.0 );
    }
    else
    {
        facX = weight * len;
        
        mC(inx  , inx  ) -= weight * axis.XX * axis.XX;
        mC(inx  , inx+1) -= weight * axis.XX * axis.YY;
        mC(inx+1, inx+1) -= weight * axis.YY * axis.YY;
    }
    
    vBAS[inx  ] += facX * axis.XX;
    vBAS[inx+1] += facX * axis.YY;
    
#endif
}



//------------------------------------------------------------------------------
#pragma mark - Side Point Clamps
//------------------------------------------------------------------------------


#if ( DIM == 2 )

void Meca::addSidePointClamp2D(PointInterpolated const& pta,
                               Vector const& pos,
                               const real arm,
                               const real weight)
{
    //force coefficients on the points:
    const real A = pta.coef2(),  wA = weight * A;
    const real B = pta.coef1(),  wB = weight * B;
    
    const real E = arm / pta.len();
    const real wE = weight * E;
    const real wEE = weight * E * E;
    
    //index in the matrix mB:
    index_type inx1 = pta.matIndex1();
    index_type inx2 = pta.matIndex2();
    
    //we put the isotropic terms in mB
    mB(inx1, inx1) -=  wA * A + wEE;
    mB(inx1, inx2) -=  wA * B - wEE;
    mB(inx2, inx2) -=  wB * B + wEE;
    
    //index in the matrix mC:
    inx1 *= DIM;
    inx2 *= DIM;
    
    mC(inx1  , inx2+1) += wE;
    mC(inx1+1, inx2  ) -= wE;
    
    //it seems to works also fine without the term in eew* below:
    vBAS[inx1  ] += wA * pos.XX - wE * pos.YY;
    vBAS[inx1+1] += wA * pos.YY + wE * pos.XX;
    vBAS[inx2  ] += wB * pos.XX + wE * pos.YY;
    vBAS[inx2+1] += wB * pos.YY - wE * pos.XX;

#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(), pta.pos()+cross(arm, pta.dir()), pos);
    }
#endif
    
    if ( modulo )
        throw Exception("addSidePointClamp2D is not usable with periodic boundary conditions");
}

#elif ( DIM == 3 )

/**
 A link of stiffness weight, between offset_point on the side of pta, and the fixed position `g`.

 This uses the vector product x -> arm ^ x to offset the point on which the link is acting:
 offset_point = fiber_point + arm ^ fiber_dir,
 with fiber_point = pta.pos() and fiber_dir = pta.diff().normalized.
 
 arm must be perpendicular to link ( g - pta.pos() )

 F. Nedelec, March 2011
 */
void Meca::addSidePointClamp3D(PointInterpolated const& pta,
                               Vector const& pos,
                               Vector const& arm,
                               real const weight)
{
    real aa = pta.coef2();
    real bb = pta.coef1();
    
    real s = 1.0 / pta.len();

    real ex = s * arm.XX;
    real ey = s * arm.YY;
    real ez = s * arm.ZZ;
    
    // indices to mC:
    const index_type inx1 = DIM * pta.matIndex1();
    const index_type inx2 = DIM * pta.matIndex2();
    const index_type inx[6] = { inx1, inx1+1, inx1+2, inx2, inx2+1, inx2+2 };
    
    /* The transfer matrix transforms the two PointExact in pta,
     to the side point S:
     S = aa * pt1 + bb * pt2 + arm ^ ( pt2 - pt1 ).normalized
     
     It was generated in Maxima:
     MVP: matrix([0, -ez, ey], [ez, 0, -ex], [-ey, ex, 0]);
     MD: addcol(-ident(3), ident(3));
     MC: addcol(aa*ident(3), bb*ident(3));
     T: MC+MVP.MD;
     */
    const real T[18] = {
         aa,  ez, -ey,  bb, -ez,  ey,
        -ez,  aa,  ex,  ez,  bb, -ex,
         ey, -ex,  aa, -ey,  ex,  bb
    };
    
#if ( 0 )
    
    real TT[36];
    // TT = transpose(T) * T
    blas_xsyrk('U','N', 6, 3, 1.0, T, 6, 0.0, TT, 6);
    
#else
    
    real a2 = aa * aa;
    real b2 = bb * bb;
    real ab = aa * bb;
    
    real exx = ex * ex, exy = ex*ey, exz = ex*ez;
    real eyy = ey * ey, eyz = ey*ez;
    real ezz = ez * ez;
    
    // TT = transpose(T) * T is symmetric, and thus we only set half of it:
    /* Maxima code:
    TT: expand(transpose(T) . T);
     */
    real TT[36] = {
        eyy+ezz+a2,  0,           0,           0,           0,           0,
        -exy,        exx+ezz+a2,  0,           0,           0,           0,
        -exz,       -eyz,         exx+eyy+a2,  0,           0,           0,
        -ezz-eyy+ab, ez+exy,      exz-ey,      eyy+ezz+b2,  0,           0,
        -ez+exy,    -ezz-exx+ab,  eyz+ex,     -exy,         exx+ezz+b2,  0,
        exz+ey,      eyz-ex,     -eyy-exx+ab, -exz,        -eyz,         exx+eyy+b2
    };
    
#endif
    
    // we project to bring all forces in the plane perpendicular to 'arm'
    real sca = 1.0 / arm.norm();
    real aan = aa * sca;
    real bbn = bb * sca;
    real TP[6] = { aan*ex, aan*ey, aan*ez, bbn*ex, bbn*ey, bbn*ez };    
    //blas_xgemm('N','N', 6, 1, 3, sca, T, 6, arm, 3, 0.0, TP, 6);

    blas_xsyrk('U','N', 6, 1, weight, TP, 6, -weight, TT, 6);
    
    for ( int ii=0; ii<6; ++ii )
    {
        for ( int jj=ii; jj<6; ++jj )
            mC(inx[ii], inx[jj]) += TT[ii+6*jj];
    }
    
    // { gx, gy, gz } is the projection of `pos` in the plane perpendicular to 'arm'
    real ws = ( arm * pos ) * sca * sca;
    real gx = weight * ( pos.XX - ws * arm.XX );
    real gy = weight * ( pos.YY - ws * arm.YY );
    real gz = weight * ( pos.ZZ - ws * arm.ZZ );
    
    for ( int ii=0; ii<6; ++ii )
        vBAS[inx[ii]] += T[ii] * gx + T[ii+6] * gy + T[ii+12] * gz;
                 
#ifdef DISPLAY_INTERACTIONS
    if ( displayInteractions )
    {
        gle::bright_color(pta.mecable()->signature()).load();
        displayInteraction(pta.pos(), pta.pos()+cross(arm, pta.dir()), pos);
    }
#endif
    
    if ( modulo )
        throw Exception("addSidePointClamp3D is not usable with periodic boundary conditions");
}

#endif  

/**
 Update Meca to include a connection between A and a fixed position G.
 The force is of zero resting length, but it is taken between B
 and another point S which is located on the side of A:
 @code
 S = A + len * N,
 force_S = weight * ( G - S )
 @endcode
 There is no counter-force in G, since G is immobile.
 */

void Meca::addSidePointClamp(PointInterpolated const& pta,
                             Vector const& pos,
                             const real len,
                             const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSidePointClamp is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    // 'arm' is a vector in the Z direction
    real arm = len * RNG.sign_exc( cross(pta.diff(), pos-pta.pos()));
    addSidePointClamp2D(pta, pos, arm, weight);
   
#elif ( DIM == 3 )
    
    // 'arm' perpendicular to link and fiber is obtained by vector product:
    Vector arm = cross( pta.pos()-pos, pta.diff() );
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        addSidePointClamp3D(pta, pos, arm * ( len / n ), weight);

#endif  
}

//------------------------------------------------------------------------------
#pragma mark - Line & Plane Clamps
//------------------------------------------------------------------------------

/**
 Add interaction between `pta` and the line defined by `G` and tangent vector `dir`.
 The force is linear with position with a stiffness `weight`, and its
 components parallel to `dir` are removed corresponding to a frictionless line:
 
 @code
 matrix P = 1 - dir (x) dir'
 force(A) = weight * P * ( G - A )
 @endcode

 The vector `dir` should be of norm = 1.
 */

void Meca::addLineClamp(const PointExact & pta,
                        const Vector & g,
                        const Vector & dir,
                        const real weight )
{
    assert_true( weight >= 0 );
    
    const index_type inx = DIM * pta.matIndex();
    const real gdir = g * dir;

    
    mC(inx  , inx  ) += ( dir.XX * dir.XX - 1 ) * weight;
    vBAS[inx  ] += ( g.XX - gdir * dir.XX ) * weight;
    
#if ( DIM >= 2 )
    mC(inx+1, inx+1) += ( dir.YY * dir.YY - 1 ) * weight;
    mC(inx  , inx+1) +=   dir.XX * dir.YY * weight;
    vBAS[inx+1] += ( g.YY - gdir * dir.YY ) * weight;
#endif
    
#if ( DIM >= 3 )
    mC(inx+2, inx+2) += ( dir.ZZ * dir.ZZ - 1 ) * weight;
    mC(inx+1, inx+2) +=   dir.YY * dir.ZZ * weight;
    mC(inx  , inx+2) +=   dir.XX * dir.ZZ * weight;
    vBAS[inx+2] += ( g.ZZ - gdir * dir.ZZ ) * weight;
#endif
}



/**
 Add interaction between `pta` and the line defined by `G` and tangent vector `dir`.
 The force is linear with position with a stiffness `weight`, and its
 components parallel to `dir` are removed corresponding to a frictionless line:
 
 @code
 matrix P = 1 - dir (x) dir'
 force(A) = weight * P * ( G - A )
 @endcode

 The vector `dir` should be of norm = 1.
 */

void Meca::addLineClamp(const PointInterpolated & pta,
                        const Vector & pos,
                        const Vector & dir,
                        const real weight )
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const index_type inx1 = DIM * pta.matIndex1();
    const index_type inx2 = DIM * pta.matIndex2();

    //force coefficients on the points:
    const real A = pta.coef2();
    const real B = pta.coef1();
    const real BB = B * B, AA = A * A, AB = A * B;
    const real gdir = pos * dir;

    //build the upper part of the projection on ga:
    for ( unsigned dd = 0; dd < DIM; ++dd )
    {
        real P = weight * ( dir[dd] * dir[dd] - 1 );
        
        mC(inx1+dd, inx1+dd) += AA * P;
        mC(inx1+dd, inx2+dd) += AB * P;
        mC(inx2+dd, inx2+dd) += BB * P;
        
        for ( unsigned ee = dd+1; ee < DIM; ++ee )
        {
            real P = weight * dir[dd] * dir[ee];
            
            mC(inx1+dd, inx1+ee) += AA * P;
            mC(inx1+dd, inx2+ee) += AB * P;
            mC(inx1+ee, inx2+dd) += AB * P;
            mC(inx2+dd, inx2+ee) += BB * P;
        }
        
        //add the constant term:
        real pg = weight * ( pos[dd] - gdir * dir[dd] );
        vBAS[inx1+dd] += pg * A;
        vBAS[inx2+dd] += pg * B;
    }
}



/**
 Add interaction between `pta` and the plane defined by `G` and the normal `dir`.
 The force is linear and the components parallel to the plane are removed,
 corresponding to an interaction with a frictionless plane:
 @code
 matrix P = dir (x) dir'
 force(A) = weight * P * ( G - A )
 @endcode
 
 The vector `dir` should be of norm = 1, or alternatively
 `weight` can be the true weigth divided by |dir|^2.
 */

void Meca::addPlaneClamp(const PointExact & pta, 
                         const Vector & pos,
                         const Vector & dir,
                         const real weight )
{
    assert_true( weight >= 0 );
    
    const index_type inx = DIM * pta.matIndex();
    const real pr = ( pos * dir ) * weight;
    
    mC(inx  , inx  ) -= dir.XX * dir.XX * weight;
    vBAS[inx  ] += pr * dir.XX;
    
#if ( DIM >= 2 )
    mC(inx+1, inx+1) -= dir.YY * dir.YY * weight;
    mC(inx+1, inx  ) -= dir.XX * dir.YY * weight;
    vBAS[inx+1] += pr * dir.YY;
#endif
    
#if ( DIM >= 3 )
    mC(inx+2, inx+2) -= dir.ZZ * dir.ZZ * weight;
    mC(inx+2, inx+1) -= dir.YY * dir.ZZ * weight;
    mC(inx+2, inx  ) -= dir.XX * dir.ZZ * weight;
    vBAS[inx+2] += pr * dir.ZZ;
#endif
}


/**
 Add interaction between `pta` and the plane defined by `G` and the normal `dir`.
 The force is linear and the perpendicular forces are removed, to create a frictionless plane:
 @code
 matrix P = dir (x) dir'
 force(A) = weight * P * ( G - A )
 @endcode
 
 The vector `dir` should be of norm = 1, or alternatively 
 `weight` can be the true weigth divided by |dir|^2.
 */

void Meca::addPlaneClamp(const PointInterpolated & pta,
                         const Vector & pos,
                         const Vector & dir,
                         const real weight )
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const index_type inx1 = DIM * pta.matIndex1();
    const index_type inx2 = DIM * pta.matIndex2();

    //force coefficients on the points:
    const real A = pta.coef2();
    const real B = pta.coef1();
    const real BB = B * B, AA = A * A, AB = A * B;
    
    //build the upper part of the projection on ga:
    const real pr = weight * ( pos * dir );
    for ( unsigned dd = 0; dd < DIM; ++dd )
    {
        real P = weight * dir[dd] * dir[dd];
        
        mC(inx1+dd, inx1+dd) -= AA * P;
        mC(inx1+dd, inx2+dd) -= AB * P;
        mC(inx2+dd, inx2+dd) -= BB * P;

        for ( unsigned ee = dd+1; ee < DIM; ++ee )
        {
            real P = weight * dir[dd] * dir[ee];
            
            mC(inx1+dd, inx1+ee) -= AA * P;
            mC(inx1+dd, inx2+ee) -= AB * P;
            mC(inx1+ee, inx2+dd) -= AB * P;
            mC(inx2+dd, inx2+ee) -= BB * P;
        }

        //add the constant term:
        vBAS[inx1+dd] += pr * A * dir[dd];
        vBAS[inx2+dd] += pr * B * dir[dd];
    }
}


//------------------------------------------------------------------------------
#pragma mark - Experimental
//------------------------------------------------------------------------------

/**
 Do not use this function!
 
 If weigth>0, this creates an attractive force that varies like 1/R^3
 */
void Meca::interCoulomb( const PointExact & pta, const PointExact & ptb, real weight )
{
    Vector ab = ptb.pos() - pta.pos();
    real abnSqr = ab.normSqr(), abn=sqrt(abnSqr);
    
    const index_type inxA = DIM * pta.matIndex();
    const index_type inxB = DIM * ptb.matIndex();
    
    if ( abn < REAL_EPSILON ) return;
    ab /= abn;
    
    real abn3 = weight / abnSqr;
    real abn5 = weight / ( abnSqr * abn );
    
    for ( int ii = 0; ii < DIM; ++ii )
    {
        vBAS[inxA+ii] -= 3 * abn3 * ab[ii];
        vBAS[inxB+ii] += 3 * abn3 * ab[ii];
    }
    
    for ( int ii = 0; ii < DIM; ++ii )
    {
        for ( int jj = ii; jj < DIM; ++jj )
        {
            real m = abn5 * ( (ii==jj) - 3 * ab[ii] * ab[jj] );
            
            mC(inxA+ii, inxA+jj) -= m;
            mC(inxB+ii, inxB+jj) -= m;
            
            mC(inxA+ii, inxB+jj) += m;
            if ( ii != jj )
                mC(inxA+jj, inxB+ii) += m;
        }
    }
}


// A simple version where stiffness is always the same
void  Meca::interSimpleTriLink(PointInterpolated const& pt1, PointInterpolated const& pt2,PointInterpolated const& pt3, real w)
{
    interLink(pt1, pt2, w);
    interLink(pt1, pt3, w);
    interLink(pt2, pt3, w);
}

void  Meca::interTriLink(PointInterpolated const& pt1, const real w1,
                         PointInterpolated const& pt2, const real w2,
                         PointInterpolated const& pt3, const real w3)
{
    /* we replace `x` obtained from the first equation into the second equation,
     to derive the stiffness of each pair of points */
    const real sum = w1 + w2 + w3;
    assert_true( sum > REAL_EPSILON );
    interLink(pt1, pt2, w1*w2/sum);
    interLink(pt1, pt3, w1*w3/sum);
    interLink(pt2, pt3, w2*w3/sum);
}

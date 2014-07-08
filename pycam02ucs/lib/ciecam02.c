/**
 * Copyright (c) 2003 Billy Biggs <vektor@dumbterm.net>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * Version 0.3
 *  - Further cleanups and a function to return all of J,C,h,Q,M,s.
 * Version 0.2
 *  - Cleanup, added missing functions.
 * Version 0.1
 *  - Initial release.
 */

#include <math.h>
#include "ciecam02.h"

static double calculate_fl_from_la_ciecam02( double la )
{
    double la5 = la * 5.0;
    double k = 1.0 / (la5 + 1.0);

    /* Calculate k^4. */
    k = k * k;
    k = k * k;

    return (0.2 * k * la5) + (0.1 * (1.0 - k)
                                  * (1.0 - k)
                                  * pow(la5, 1.0 / 3.0));
}

/**
 *              [  0.7328  0.4296  -0.1624 ]
 *    M_CAT02 = [ -0.7036  1.6975   0.0061 ]
 *              [  0.0030  0.0136   0.9834 ]
 *
 *              [  1.096124 -0.278869 0.182745 ]
 * M^-1_CAT02 = [  0.454369  0.473533 0.072098 ]
 *              [ -0.009628 -0.005698 1.015326 ]
 */
static void xyz_to_cat02( double *r, double *g, double *b,
                          double x, double y, double z )
{
    *r = ( 0.7328 * x) + (0.4296 * y) - (0.1624 * z);
    *g = (-0.7036 * x) + (1.6975 * y) + (0.0061 * z);
    *b = ( 0.0030 * x) + (0.0136 * y) + (0.9834 * z);
}

static void cat02_to_xyz( double *x, double *y, double *z,
                          double r, double g, double b )
{
    *x = ( 1.096124 * r) - (0.278869 * g) + (0.182745 * b);
    *y = ( 0.454369 * r) + (0.473533 * g) + (0.072098 * b);
    *z = (-0.009628 * r) - (0.005698 * g) + (1.015326 * b);
}

static void hpe_to_xyz( double *x, double *y, double *z,
                        double r, double g, double b )
{
    *x = (1.91020 * r) - (1.11212 * g) + (0.20191 * b);
    *y = (0.37095 * r) + (0.62905 * g) - (0.00001 * b);
    *z = b;
}

static void cat02_to_hpe( double *rh, double *gh, double *bh,
                          double r, double g, double b )
{
    *rh = ( 0.7409792 * r) + (0.2180250 * g) + (0.0410058 * b);
    *gh = ( 0.2853532 * r) + (0.6242014 * g) + (0.0904454 * b);
    *bh = (-0.0096280 * r) - (0.0056980 * g) + (1.0153260 * b);
}

/**
 * Theoretically, D ranges from
 *     0 = no adaptation to the adopted white point,
 *  to 1 = complete adaptation to the adopted white point.
 * In practice, the minimum D value will not be less than 0.65 for a
 * dark surround and exponentially converges to 1 for average surrounds
 * with increasingly large values of L_A.
 *
 * L_A is the luminance of the adapting field in cd/m^2.
 */
static double d_factor( double f, double la )
{
    return f * (1.0 - ((1.0 / 3.6) * exp((-la - 42.0) / 92.0)));
}

static double nonlinear_adaptation( double c, double fl )
{
    double p = pow( (fl * c) / 100.0, 0.42 );
    return ((400.0 * p) / (27.13 + p)) + 0.1;
}

static double inverse_nonlinear_adaptation( double c, double fl )
{
    return (100.0 / fl) * pow( (27.13 * fabs( c - 0.1 )) / (400.0 - fabs( c - 0.1 )), 1.0 / 0.42 );
}

static double achromatic_response_to_white( double x, double y, double z,
                                            double d, double fl, double nbb )
{
    double r, g, b;
    double rc, gc, bc;
    double rp, gp, bp;
    double rpa, gpa, bpa;

    xyz_to_cat02( &r, &g, &b, x, y, z );

    rc = r * (((y * d) / r) + (1.0 - d));
    gc = g * (((y * d) / g) + (1.0 - d));
    bc = b * (((y * d) / b) + (1.0 - d));

    cat02_to_hpe( &rp, &gp, &bp, rc, gc, bc );

    rpa = nonlinear_adaptation( rp, fl );
    gpa = nonlinear_adaptation( gp, fl );
    bpa = nonlinear_adaptation( bp, fl );

    return ((2.0 * rpa) + gpa + ((1.0 / 20.0) * bpa) - 0.305) * nbb;
}

void xyz2jch_ciecam02( double *J, double *C, double *h,
                       double x, double y, double z,
                       double xw, double yw, double zw,
                       double yb, double la,
                       double f, double c, double nc )
{
    double r, g, b;
    double rw, gw, bw;
    double rc, gc, bc;
    double rp, gp, bp;
    double rpa, gpa, bpa;
    double a, ca, cb;
    double fl, d;
    double n, nbb, ncb;
    double aw;
    double e, t;
    double cz;

    xyz_to_cat02( &r, &g, &b, x, y, z );
    xyz_to_cat02( &rw, &gw, &bw, xw, yw, zw );

    n = yb / yw;
    d = d_factor( f, la );
    fl = calculate_fl_from_la_ciecam02( la );
    nbb = ncb = 0.725 * pow( 1.0 / n, 0.2 );
    cz = 1.48 + sqrt( n );
    aw = achromatic_response_to_white( xw, yw, zw, d, fl, nbb );

    rc = r * (((yw * d) / rw) + (1.0 - d));
    gc = g * (((yw * d) / gw) + (1.0 - d));
    bc = b * (((yw * d) / bw) + (1.0 - d));

    cat02_to_hpe( &rp, &gp, &bp, rc, gc, bc );

    rpa = nonlinear_adaptation( rp, fl );
    gpa = nonlinear_adaptation( gp, fl );
    bpa = nonlinear_adaptation( bp, fl );

    ca = rpa - ((12.0 * gpa) / 11.0) + (bpa / 11.0);
    cb = (1.0 / 9.0) * (rpa + gpa - (2.0 * bpa));

    *h = (180.0 / M_PI) * atan2( cb, ca );
    if( *h < 0.0 ) *h += 360.0;

    a = ((2.0 * rpa) + gpa + ((1.0 / 20.0) * bpa) - 0.305) * nbb;

    *J = 100.0 * pow( a / aw, c * cz );

    e = ((12500.0 / 13.0) * nc * ncb) * (cos( ((*h * M_PI) / 180.0) + 2.0 ) + 3.8);
    t = (e * sqrt( (ca * ca) + (cb * cb) )) / (rpa + gpa + ((21.0 / 20.0) * bpa));

    *C = pow( t, 0.9 ) * sqrt( *J / 100.0 )
                       * pow( 1.64 - pow( 0.29, n ), 0.73 );
}

void xyz2jchqms_ciecam02( double *J, double *C, double *h,
                          double *Q, double *M, double *s,
                          double x, double y, double z,
                          double xw, double yw, double zw,
                          double yb, double la,
                          double f, double c, double nc )
{
    double r, g, b;
    double rw, gw, bw;
    double rc, gc, bc;
    double rp, gp, bp;
    double rpa, gpa, bpa;
    double a, ca, cb;
    double fl, d;
    double n, nbb, ncb;
    double aw;
    double e, t;
    double cz;
    double myh, myj, myc, myq, mym, mys;

    xyz_to_cat02( &r, &g, &b, x, y, z );
    xyz_to_cat02( &rw, &gw, &bw, xw, yw, zw );

    n = yb / yw;
    d = d_factor( f, la );
    fl = calculate_fl_from_la_ciecam02( la );
    nbb = ncb = 0.725 * pow( 1.0 / n, 0.2 );
    cz = 1.48 + sqrt( n );
    aw = achromatic_response_to_white( xw, yw, zw, d, fl, nbb );

    rc = r * (((yw * d) / rw) + (1.0 - d));
    gc = g * (((yw * d) / gw) + (1.0 - d));
    bc = b * (((yw * d) / bw) + (1.0 - d));

    cat02_to_hpe( &rp, &gp, &bp, rc, gc, bc );

    rpa = nonlinear_adaptation( rp, fl );
    gpa = nonlinear_adaptation( gp, fl );
    bpa = nonlinear_adaptation( bp, fl );

    ca = rpa - ((12.0 * gpa) / 11.0) + (bpa / 11.0);
    cb = (1.0 / 9.0) * (rpa + gpa - (2.0 * bpa));

    myh = (180.0 / M_PI) * atan2( cb, ca );
    if( myh < 0.0 ) myh += 360.0;

    a = ((2.0 * rpa) + gpa + ((1.0 / 20.0) * bpa) - 0.305) * nbb;

    myj = 100.0 * pow( a / aw, c * cz );

    e = ((12500.0 / 13.0) * nc * ncb) * (cos( ((myh * M_PI) / 180.0) + 2.0 ) + 3.8);
    t = (e * sqrt( (ca * ca) + (cb * cb) )) / (rpa + gpa + ((21.0 / 20.0) * bpa));

    myc = pow( t, 0.9 ) * sqrt( myj / 100.0 )
                        * pow( 1.64 - pow( 0.29, n ), 0.73 );

    myq = ( 4.0 / c ) * sqrt( myj / 100.0 ) * ( aw + 4.0 ) * pow( fl, 0.25 );
    mym = myc * pow( fl, 0.25 );
    mys = 100.0 * sqrt( mym / myq );

    if( J ) *J = myj;
    if( C ) *C = myc;
    if( h ) *h = myh;
    if( Q ) *Q = myq;
    if( M ) *M = mym;
    if( s ) *s = mys;
}


static void Aab_to_rgb( double *r, double *g, double *b, double A, double aa,
                        double bb, double nbb )
{
    double x = (A / nbb) + 0.305;

    /*       c1              c2               c3       */
    *r = (0.32787 * x) + (0.32145 * aa) + (0.20527 * bb);
    /*       c1              c4               c5       */
    *g = (0.32787 * x) - (0.63507 * aa) - (0.18603 * bb);
    /*       c1              c6               c7       */
    *b = (0.32787 * x) - (0.15681 * aa) - (4.49038 * bb);
}

static void calculate_ab( double *aa, double *bb, double h, double e, double t,
                          double nbb, double a )
{
    double hrad = (h * M_PI) / 180.0;
    double sinh = sin( hrad );
    double cosh = cos( hrad );
    double x = (a / nbb) + 0.305;

    if( fabs( sinh ) >= fabs( cosh ) ) {
        *bb =    ((0.32787 * x) * (2.0 + (21.0 / 20.0))) /
                 ((e / (t * sinh)) -
                  ((0.32145 - 0.63507 - ((21.0 / 20.0) * 0.15681)) * (cosh / sinh)) -
                   (0.20527 - 0.18603 - ((21.0 / 20.0) * 4.49038)));
        *aa = (*bb * cosh) / sinh;
    } else {
        *aa =    ((0.32787 * x) * (2.0 + (21.0 / 20.0))) /
                 ((e / (t * cosh)) -
                   (0.32145 - 0.63507 - ((21.0 / 20.0) * 0.15681)) -
                  ((0.20527 - 0.18603 - ((21.0 / 20.0) * 4.49038)) * (sinh / cosh)));
        *bb = (*aa * sinh) / cosh;
    }
}

void jch2xyz_ciecam02( double *x, double *y, double *z,
                       double J, double C, double h,
                       double xw, double yw, double zw,
                       double yb, double la,
                       double f, double c, double nc )
{
    double r, g, b;
    double rc, gc, bc;
    double rp, gp, bp;
    double rpa, gpa, bpa;
    double rw, gw, bw;
    double fl, d;
    double n, nbb, ncb;
    double a, ca, cb;
    double aw;
    double e, t;
    double cz;

    xyz_to_cat02( &rw, &gw, &bw, xw, yw, zw );

    n = yb / yw;
    d = d_factor( f, la );
    fl = calculate_fl_from_la_ciecam02( la );
    nbb = ncb = 0.725 * pow( 1.0 / n, 0.2 );
    cz = 1.48 + sqrt( n );
    aw = achromatic_response_to_white( xw, yw, zw, d, fl, nbb );

    e = ((12500.0 / 13.0) * nc * ncb) * (cos( ((h * M_PI) / 180.0) + 2.0 ) + 3.8);
    a = pow( J / 100.0, 1.0 / (c * cz) ) * aw;
    t = pow( C / (sqrt( J / 100) * pow( 1.64 - pow( 0.29, n ), 0.73 )), 10.0 / 9.0 );

    calculate_ab( &ca, &cb, h, e, t, nbb, a );
    Aab_to_rgb( &rpa, &gpa, &bpa, a, ca, cb, nbb );

    rp = inverse_nonlinear_adaptation( rpa, fl );
    gp = inverse_nonlinear_adaptation( gpa, fl );
    bp = inverse_nonlinear_adaptation( bpa, fl );

    hpe_to_xyz( x, y, z, rp, gp, bp );
    xyz_to_cat02( &rc, &gc, &bc, *x, *y, *z );

    r = rc / (((yw * d) / rw) + (1.0 - d));
    g = gc / (((yw * d) / gw) + (1.0 - d));
    b = bc / (((yw * d) / bw) + (1.0 - d));

    cat02_to_xyz( x, y, z, r, g, b );
}


#include "fit-utils.h"

#include "datatypes.h"
#include "niget-common.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <libprocess/gwyprocess.h>
#include <libgwyddion/gwynlfit.h>

#ifndef NOFORTRAN
#ifdef _MSC_VER /* If compiling with MSVC, we expect Intel Fortran to be used, which generates different function names */
extern void ODR_WRAP_POWER_3PARAM(gint *n, gint *m, const gint *np, gint *nq, gdouble *x, gdouble *y, gdouble *beta, const gint *ifix,
                                  gchar *fnm, gint *lenfnm, gdouble *ux, gdouble *uy, gdouble *std_beta, gdouble *covariance,
                                  gdouble *ressum,  gint *info);
extern void ODR_WRAP_POWER_2PARAM(gint *n, gint *m, const gint *np, gint *nq, gdouble *x, gdouble *y, gdouble *beta, const gint *ifix,
                                  gchar *fnm, gint *lenfnm, gdouble *ux, gdouble *uy, gdouble *std_beta, gdouble *covariance,
                                  gdouble *ressum,  gint *info);
extern void ODR_WRAP_HERTZ_RADIAL(gint *n, gint *m, const gint *np, gint *nq, gdouble *x, gdouble *y, gdouble *beta, const gint *ifix,
                                  gchar *fnm, gint *lenfnm, gdouble *ux, gdouble *uy, gdouble *std_beta, gdouble *covariance,
                                  gdouble *ressum,  gint *info);
extern void ODR_WRAP_LINEAR(gint *n, gint *m, const gint *np, gint *nq, gdouble *x, gdouble *y, gdouble *beta, const gint *ifix,
                            gchar *fnm, gint *lenfnm, gdouble *ux, gdouble *uy, gdouble *std_beta, gdouble *covariance,
                            gdouble *ressum,  gint *info);
#else
extern void odr_wrap_power_3param_(gint *n, gint *m, const gint *np, gint *nq, gdouble *x, gdouble *y, gdouble *beta, const gint *ifix,
                                   gchar *fnm, gint *lenfnm, gdouble *ux, gdouble *uy, gdouble *std_beta, gdouble *covariance,
                                   gdouble *ressum,  gint *info);
extern void odr_wrap_power_2param_(gint *n, gint *m, const gint *np, gint *nq, gdouble *x, gdouble *y, gdouble *beta, const gint *ifix,
                                   gchar *fnm, gint *lenfnm, gdouble *ux, gdouble *uy, gdouble *std_beta, gdouble *covariance,
                                   gdouble *ressum,  gint *info);
extern void odr_wrap_hertz_radial_(gint *n, gint *m, const gint *np, gint *nq, gdouble *x, gdouble *y, gdouble *beta, const gint *ifix,
                                   gchar *fnm, gint *lenfnm, gdouble *ux, gdouble *uy, gdouble *std_beta, gdouble *covariance,
                                   gdouble *ressum,  gint *info);
extern void odr_wrap_linear_(gint *n, gint *m, const gint *np, gint *nq, gdouble *x, gdouble *y, gdouble *beta, const gint *ifix,
                             gchar *fnm, gint *lenfnm, gdouble *ux, gdouble *uy, gdouble *std_beta, gdouble *covariance,
                             gdouble *ressum,  gint *info);
#endif
#endif


#define NParam 3
static gdouble fun(gdouble x, G_GNUC_UNUSED gint nparam, gdouble *param, G_GNUC_UNUSED gpointer user_data, gboolean *fres);

void ordinary_least_squares_line_fit__(GwyDataLine *x, GwyDataLine *y, LinReg *reg, gdouble ux, gdouble uy, gboolean correct, gboolean rescale)
{
    gint i, n;
    double k, q, p, r, c;
    GwyDataLine *w, *z;
    gdouble xmin, xmax, ymin, ymax;

    if (rescale) {

        xmax = gwy_data_line_get_max(x);
        ymax = gwy_data_line_get_max(y);

        xmin = gwy_data_line_get_min(x);
        ymin = gwy_data_line_get_min(y);

        k = 2 / (xmax - xmin);
        p = 2 / (ymax - ymin);

        q =  -(xmin + xmax) / (xmax - xmin);
        r =  -(ymin + ymax) / (ymax - ymin);

        w = gwy_data_line_new_alike(x, FALSE);
        z = gwy_data_line_new_alike(y, FALSE);

        n = x->res;

        for (i = 0; i < n; i++) {
            w->data[i] = k * x->data[i] + q;
            z->data[i] = p * y->data[i] + r;
        }

        ordinary_least_squares_line_fit(w, z, reg);

        reg->intercept = (reg->intercept + reg->slope * q - r) / p;
        reg->slope *= k / p;

        g_object_unref(w);
        g_object_unref(z);
    }
    else {
        ordinary_least_squares_line_fit(x, y, reg);
    }

    if (correct) {
        c = least_squares_line_get_correction(x, ux);
        reg->slope *= c;
        reg->intercept = gwy_data_line_get_avg(y) - reg->slope * gwy_data_line_get_avg(x);
    }
}

void ordinary_least_squares_line_fit(GwyDataLine *x, GwyDataLine *y, LinReg *reg)
{
    gdouble sx, sy, sxx, sxy;
    gint i, n;
    gdouble a, b;

    n = x->res;

    sx = 0;
    sy = 0;
    sxx = 0;
    sxy = 0;

    for (i = 0; i < n; i++) {
        sx += x->data[i];
        sy += y->data[i];
        sxx += x->data[i] * x->data[i];
        sxy += x->data[i] * y->data[i];
    }

    a = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    b = (sxx * sy - sx * sxy) / (n * sxx - sx * sx);

    reg->slope = a;
    reg->intercept = b;
}

void ordinary_least_squares_line_fit_unc_const(GwyDataLine *x, GwyDataLine *y, LinReg *reg, gdouble ux, gdouble uy)
{
    gdouble sx, sy, sxx, sxy, syy;
    gint i, n;
    gdouble a, b ;
    gdouble W;
    gdouble uax, uay, ubx, uby, uabx, uaby;

    n = x->res;

    sx = 0;
    sy = 0;
    sxx = 0;
    sxy = 0;
    syy = 0;

    for (i = 0; i < n; i++) {
        sx += x->data[i];
        sy += y->data[i];
        sxx += x->data[i] * x->data[i];
        sxy += x->data[i] * y->data[i];
        syy += y->data[i] * y->data[i];
    }

    W = n * sxx - sx * sx;

    a = (n * sxy - sx * sy) / W;
    /* no correction for regression attenuation */
    b = (sy - a * sx) / n;

    reg->slope = a;
    reg->intercept = b;

    /* uncertainty components */
    uax = n * ux * ux / W * (n * syy - sy * sy);
    uay = n * uy * uy / W;
    ubx = uax * sx * sx / n / n + ux * ux * a * a / n;
    uby = uay * sx * sx / n / n + uy * uy / n;
    uabx = -uax * sx / n;
    uaby = -uay * sy / n;

    reg->unc_slope_x =  sqrt(uax);
    reg->unc_slope_y =  sqrt(uay);
    reg->unc_slope = sqrt(uax + uay);

    reg->unc_intercept_x = sqrt(ubx);
    reg->unc_intercept_y = sqrt(uby);
    reg->unc_intercept = sqrt(ubx + uby);

    reg->cov_slope_intercept_x = uabx;
    reg->cov_slope_intercept_y =  uaby;
    reg->cov_slope_intercept = uabx + uaby;
}

void ordinary_least_squares_line_fit_unc_general(GwyDataLine *x, GwyDataLine *y, LinReg *reg, GwyDataLine *ux, GwyDataLine *uy)
{
    gdouble sx, sy, sxx, sxy;
    gint i, n;
    gdouble a, b;
    gdouble  c, W;
    gdouble *axi, *bxi, *ayi, *byi;
    gdouble uax, uay, ubx, uby, uabx, uaby;
    gdouble sux2;

    n = x->res;

    sx = 0;
    sy = 0;
    sxx = 0;
    sxy = 0;

    for (i = 0; i < n; i++) {
        sx += x->data[i];
        sy += y->data[i];
        sxx += x->data[i] * x->data[i];
        sxy += x->data[i] * y->data[i];
    }

    W = n * sxx - sx * sx;

    a = (n * sxy - sx * sy) / W;

    /* correction for regression attenuation */
    sux2 = 0;

    for (i = 0; i < n; i++) {
        sux2 += ux->data[i] * ux->data[i];
    }

    c = 1 + sux2 / sxx;
    a *= c;

    b = (sy - a * sx) / n;

    reg->slope = a;
    reg->intercept = b;

    /* derivative of a and b with respect to x[i] and y[i] */
    axi = (gdouble *)g_malloc(n * sizeof(gdouble));
    ayi = (gdouble *)g_malloc(n * sizeof(gdouble));
    bxi = (gdouble *)g_malloc(n * sizeof(gdouble));
    byi = (gdouble *)g_malloc(n * sizeof(gdouble));

    for (i = 0; i < n; i++) {
        axi[i] = (n * y->data[i] - sy) / W - 2 / W / W * (n * sxy - sx * sy) * (n * x->data[i] - sx);
        ayi[i] = (n * x->data[i] - sx) / W;
        axi[i] *= c;
        ayi[i] *= c;
        bxi[i] = (-axi[i] * sx - a) / n;
        byi[i] = (1 - ayi[i] * sx) / n;
    }

    uax = 0;
    uay = 0;
    ubx = 0;
    uby = 0;
    uabx = 0;
    uaby = 0;

    /* propagate uncertainties
     * u(a;x)^2 = sum_i (da/dx[i])^2 ux[i]^2
     * u(a;y)^2 = sum_i (da/dy[i])^2 uy[i]^2
     * u(a)^2 = u(a;x)^2 + u(a;y)^2
     *
     * u(b;x)^2 = sum_i (db/dx[i])^2 ux[i]^2
     * u(b;y)^2 = sum_i (db/dy[i])^2 uy[i]^2
     * u(b)^2 = u(b;x)^2 + u(b;y)^2
     *
     * cov(a,b;x) = sum_i (da/dx[i])(db/dx[i]) ux[i]^2)
     * cov(a,b;y) = sum_i (da/dy[i])(db/dy[i]) uy[i]^2)
     * cov(a,b) = cov(a,b;x)+ cov(a,b;y)
    */
    for (i = 0; i < n; i++) {
        uax += axi[i] * axi[i] * ux->data[i] * ux->data[i];
        ubx += bxi[i] * bxi[i] * ux->data[i] * ux->data[i];
        uay += ayi[i] * ayi[i] * uy->data[i] * uy->data[i];
        uby += byi[i] * byi[i] * uy->data[i] * uy->data[i];
        uabx += axi[i] * bxi[i] * ux->data[i] * ux->data[i];
        uaby += ayi[i] * byi[i] * uy->data[i] * uy->data[i];
    }

    reg->unc_slope_x =  sqrt(uax);
    reg->unc_slope_y =  sqrt(uay);
    reg->unc_slope = sqrt(uax + uay);

    reg->unc_intercept_x = sqrt(ubx);
    reg->unc_intercept_y = sqrt(uby);
    reg->unc_intercept = sqrt(ubx + uby);

    reg->cov_slope_intercept_x = uabx;
    reg->cov_slope_intercept_y =  uaby;
    reg->cov_slope_intercept = uabx + uaby;

    g_free(axi);
    g_free(ayi);
    g_free(bxi);
    g_free(byi);
}

/* correction for regression attenuation */
gdouble least_squares_line_get_correction(GwyDataLine *x, gdouble ux)
{
    gdouble xavg, sxx;
    gint i, n;

    xavg = gwy_data_line_get_avg(x);

    n = x->res;

    sxx = 0;

    for (i = 0; i < n; i++) {
        sxx += (x->data[i] - xavg) * (x->data[i] - xavg);
    }

    sxx /= n - 1;

    return 1 + ux * ux / sxx;
}

void total_least_squares_line_fit__(GwyDataLine *x, GwyDataLine *y, LinReg *reg, gdouble ux, gdouble uy, gboolean correct, gboolean rescale, gboolean delta)
{
    gint i, n;
    double k, q, p, r, c;
    GwyDataLine *w, *z;
    gdouble xmin, xmax, ymin, ymax;

    if (delta) {
        rescale = FALSE;    //delta and rescaling is redundant
    }

    if (rescale) {

        xmax = gwy_data_line_get_max(x);
        ymax = gwy_data_line_get_max(y);

        xmin = gwy_data_line_get_min(x);
        ymin = gwy_data_line_get_min(y);

        k = 2 / (xmax - xmin);
        p = 2 / (ymax - ymin);

        q =  -(xmin + xmax) / (xmax - xmin);
        r =  -(ymin + ymax) / (ymax - ymin);

        w = gwy_data_line_new_alike(x, FALSE);
        z = gwy_data_line_new_alike(y, FALSE);

        n = x->res;

        for (i = 0; i < n; i++) {
            w->data[i] = k * x->data[i] + q;
            z->data[i] = p * y->data[i] + r;
        }

        total_least_squares_line_fit_delta(w, z, reg, 1);

        reg->intercept = (reg->intercept + reg->slope * q - r) / p;
        reg->slope *= k / p;

        g_object_unref(w);
        g_object_unref(z);
    }
    else {
        if (delta && ux > 0) {
            total_least_squares_line_fit_delta(x, y, reg, uy * uy / ux / ux);
        }
        else {
            total_least_squares_line_fit_delta(x, y, reg, 1);
        }
    }

    if (correct) {
        c = least_squares_line_get_correction(x, ux);
        reg->slope *= c;
        reg->intercept = gwy_data_line_get_avg(y) - reg->slope * gwy_data_line_get_avg(x);
    }
}

void total_least_squares_line_fit(GwyDataLine *x, GwyDataLine *y, LinReg *reg)
{
    gdouble xavg, yavg, sxx, syy, sxy;
    gint i, n;
    gdouble a, b;

    xavg = gwy_data_line_get_avg(x);
    yavg = gwy_data_line_get_avg(y);

    n = x->res;

    sxx = 0;
    syy = 0;
    sxy = 0;

    for (i = 0; i < n; i++) {
        sxx += (x->data[i] - xavg) * (x->data[i] - xavg);
        sxy += (x->data[i] - xavg) * (y->data[i] - yavg);
        syy += (y->data[i] - yavg) * (y->data[i] - yavg);
    }

    sxy /= n;
    sxx /= n;
    syy /= n;

    reg->corr = sxx * syy;


    a = (syy - sxx + sqrt((syy - sxx) * (syy - sxx) + 4 * sxy * sxy)) / (2 * sxy);
    b = yavg - a * xavg;

    reg->slope = a;
    reg->intercept = b;
}

void total_least_squares_line_fit_delta(GwyDataLine *x, GwyDataLine *y, LinReg *reg, gdouble delta)
{
    gdouble xavg, yavg, sxx, syy, sxy;
    gint i, n;
    gdouble a, b;
    gdouble *xstar;
    gdouble SSres;

    xavg = gwy_data_line_get_avg(x);
    yavg = gwy_data_line_get_avg(y);

    n = x->res;

    sxx = 0;
    syy = 0;
    sxy = 0;

    for (i = 0; i < n; i++) {
        sxx += (x->data[i] - xavg) * (x->data[i] - xavg);
        sxy += (x->data[i] - xavg) * (y->data[i] - yavg);
        syy += (y->data[i] - yavg) * (y->data[i] - yavg);
    }

    sxy /= n;
    sxx /= n;
    syy /= n;

    a = (syy - delta * sxx + sqrt((syy - delta * sxx) * (syy - delta * sxx) + 4 * delta * sxy * sxy)) / (2 * sxy);
    b = yavg - a * xavg;

    reg->slope = a;
    reg->intercept = b;

    xstar = (gdouble *)g_malloc(n * sizeof(gdouble));

    for (i = 0; i < n; i++) {
        xstar[i] = x->data[i] + a / (a * a + delta) * (y->data[i] - b - a * x->data[i]);
    }

    SSres = 0;

    for (i = 0; i < n; i++) {
        SSres += (y->data[i] - a * xstar[i] - b) * (y->data[i] - a * xstar[i] - b) + delta * (x->data[i] - xstar[i]) * (x->data[i] - xstar[i]);
    }

    g_free(xstar);

    reg->chi2 =  SSres / (syy * n + delta * n * sxx);
    reg->SSres = SSres;
    reg->R2 = 1 - reg->chi2;
    reg->R2adj = 1 - (1 - reg->R2) * (n - 1) / (n - 2 - 1); // p=2

    // R2 , R2adj
    // SStot = syy *n + delta * n*sxx;
    // SSres = [ sum(yi-ax^*i -b)^2 +delta sum(xi-x^i*)^2
    // set sigma_y^2 = syy*n
    // R2 = 1 -SSres/SStot = 1-  [sum (yi-ax^*i -b)^2 + delta sum(xi-x^*i)^2] / (syy *n +delta *n*sxx) = 1-chi2
    // R2adj = 1 - (1-R2) (n-1)/(n-p-1), p=2
}

void total_least_squares_line_fit_const(GwyDataLine *x, GwyDataLine *y, LinReg *reg, gdouble ux, gdouble uy)
{
    gdouble xavg, yavg, sxx, syy, sxy;
    gint i, n;
    gdouble a, b;
    gdouble uax, uay, ubx, uby, uabx, uaby;
    gdouble asxx, asxy, asyy;
    gdouble sdet;
    gdouble delta;
    gdouble *xstar;
    gdouble SSres;

    xavg = gwy_data_line_get_avg(x);
    yavg = gwy_data_line_get_avg(y);

    n = x->res;

    sxx = 0;
    syy = 0;
    sxy = 0;

    for (i = 0; i < n; i++) {
        sxx += (x->data[i] - xavg) * (x->data[i] - xavg);
        sxy += (x->data[i] - xavg) * (y->data[i] - yavg);
        syy += (y->data[i] - yavg) * (y->data[i] - yavg);
    }

    sxy /= n;
    sxx /= n;
    syy /= n;

    delta = uy * uy / ux / ux;
    a = (syy - delta * sxx + sqrt((syy - delta * sxx) * (syy - delta * sxx) + 4 * delta * sxy * sxy)) / (2 * sxy);
    b = yavg - a * xavg;

    reg->slope = a;
    reg->intercept = b;

    /* derivative of a with respect to sxx, sxy, syy */
    sdet = sqrt((syy - delta * sxx) * (syy - delta * sxx) + 4 * delta * sxy * sxy);
    asxx = -1 / (2 * sxy) * (1 + (syy - delta * sxx) / sdet);
    asyy = 1 / (2 * sxy) * (1 + (syy - delta * sxx) / sdet);
    asxy = -1 / (2 * sxy * sxy) * (syy - delta * sxx + (syy - delta * sxx) * (syy - delta * sxx) / sdet);

    /* uncertainty components */
    uax = (asxx * asxx * 4 / n * sxx + asxy * asxy * syy / n + 2 * asxy * asxx * 2 * sxy / n) * ux * ux;
    uay = (asyy * asyy * 4 / n * syy + asxy * asxy * sxx / n + 2 * asxy * asyy * 2 * sxy / n) * uy * uy;
    ubx = uax * xavg * xavg + ux * ux * a * a / n;
    uby = uy * uy / n + uay * xavg * xavg;
    uabx = -xavg * uax;
    uaby = -xavg * uay;

    reg->unc_slope_x =  sqrt(uax);
    reg->unc_slope_y =  sqrt(uay);
    reg->unc_slope = sqrt(uax + uay);

    reg->unc_intercept_x = sqrt(ubx);
    reg->unc_intercept_y = sqrt(uby);
    reg->unc_intercept = sqrt(ubx + uby);

    reg->cov_slope_intercept_x = uabx;
    reg->cov_slope_intercept_y =  uaby;
    reg->cov_slope_intercept = uabx + uaby;

    // R2 upravit
    delta = 1;
    xstar = (gdouble *)g_malloc(n * sizeof(gdouble));

    for (i = 0; i < n; i++) {
        xstar[i] = x->data[i] + a / (a * a + delta) * (y->data[i] - b - a * x->data[i]);
    }

    SSres = 0;

    for (i = 0; i < n; i++) {
        SSres += (y->data[i] - a * xstar[i] - b) * (y->data[i] - a * xstar[i] - b) + delta * (x->data[i] - xstar[i]) * (x->data[i] - xstar[i]);
    }

    g_free(xstar);

    reg->chi2 =  SSres / (syy * n + delta * n * sxx);
    reg->R2 = 1 - reg->chi2;
    reg->R2adj = 1 - (1 - reg->R2) * (n - 1) / (n - 2 - 1); // p=2
}

// TODO uses Deming fit with delta = 1, but uncertainties for ux, uy varying
void total_least_squares_line_fit_nonconst(GwyDataLine *x, GwyDataLine *y, LinReg *reg, GwyDataLine *ux, GwyDataLine *uy)
{
    gdouble xavg, yavg, sxx, syy, sxy;
    gint i, n;
    gdouble *axi, *ayi, *bxi, *byi;
    gdouble *xstar;
    gdouble a, b;
    gdouble uax, uay, ubx, uby, uabx, uaby;
    gdouble asxx, asxy, asyy;
    gdouble sdet;

    xavg = gwy_data_line_get_avg(x);
    yavg = gwy_data_line_get_avg(y);

    n = x->res;

    sxx = 0;
    syy = 0;
    sxy = 0;

    for (i = 0; i < n; i++) {
        sxx += (x->data[i] - xavg) * (x->data[i] - xavg);
        sxy += (x->data[i] - xavg) * (y->data[i] - yavg);
        syy += (y->data[i] - yavg) * (y->data[i] - yavg);
    }

    sxy /= n;
    sxx /= n;
    syy /= n;

    reg->corr = sxx * syy;

    a = (syy - sxx + sqrt((syy - sxx) * (syy - sxx) + 4 * sxy * sxy)) / (2 * sxy);
    b = yavg - a * xavg;

    reg->slope = a;
    reg->intercept = b;

    xstar = (gdouble *)malloc(n * sizeof(gdouble));

    for (i = 0; i < n; i++) {
        xstar[i] = x->data[i] + a / (a * a + 1) * (y->data[i] - b - a * x->data[i]);
    }

    /* derivative of a and b with respect to x[i] and y[i] */
    axi = (gdouble *)malloc(n * sizeof(gdouble));
    ayi = (gdouble *)malloc(n * sizeof(gdouble));
    bxi = (gdouble *)malloc(n * sizeof(gdouble));
    byi = (gdouble *)malloc(n * sizeof(gdouble));

    /* derivative of a with respect to sxx, sxy, syy */
    sdet = sqrt((syy - sxx) * (syy - sxx) + 4 * sxy * sxy);
    asxx = -1 / (2 * sxy) * (1 + (syy - sxx) / sdet);
    asyy = 1 / (2 * sxy) * (1 + (syy - sxx) / sdet);
    asxy = -1 / (2 * sxy * sxy) * (syy - sxx + (syy - sxx) * (syy - sxx) / sdet);

    for (i = 0; i < n; i++) {
        axi[i] =  asxy * (y->data[i] - yavg) / n + asxx * 2 * (x->data[i] - xavg) / n;

        ayi[i] =  asxy * (x->data[i] - xavg) / n + asyy * 2 * (y->data[i] - yavg) / n ;

        bxi[i] = -axi[i] * xavg - a / n;
        byi[i] = 1. / n - ayi[i] * xavg ;
    }

    uax = 0;
    uay = 0;
    ubx = 0;
    uby = 0;
    uabx = 0;
    uaby = 0;

    /* propagate uncertainties
     * u(a;x)^2 = sum_i (da/dx[i])^2 ux[i]^2
     * u(b;x)^2 = sum_i (db/dx[i])^2 ux[i]^2
     * u(a)^2 = u(a;x)^2 + u(a;y)^2
     *
     * u(a;y)^2 = sum_i (da/dy[i])^2 uy[i]^2
     * u(b;y)^2 = sum_i (db/dy[i])^2 uy[i]^2
     * u(b)^2 = u(b;x)^2 + u(b;y)^2
     *
     * cov(a,b;x) = sum_i (da/dx[i])(db/dx[i]) ux[i]^2)
     * cov(a,b;y) = sum_i (da/dy[i])(db/dy[i]) uy[i]^2)
     * cov(a,b) = cov(a,b;x)+ cov(a,b;y)
    */
    for (i = 0; i < n; i++) {
        uax += axi[i] * axi[i] * ux->data[i] * ux->data[i];
        ubx += bxi[i] * bxi[i] * ux->data[i] * ux->data[i];
        uay += ayi[i] * ayi[i] * uy->data[i] * uy->data[i];
        uby += byi[i] * byi[i] * uy->data[i] * uy->data[i];
        uabx += axi[i] * bxi[i] * ux->data[i] * ux->data[i];
        uaby += ayi[i] * byi[i] * uy->data[i] * uy->data[i];
    }

    reg->unc_slope_x =  sqrt(uax);
    reg->unc_slope_y =  sqrt(uay);
    reg->unc_slope = sqrt(uax + uay);

    reg->unc_intercept_x = sqrt(ubx);
    reg->unc_intercept_y = sqrt(uby);
    reg->unc_intercept = sqrt(ubx + uby);

    reg->cov_slope_intercept_x = uabx;
    reg->cov_slope_intercept_y =  uaby;
    reg->cov_slope_intercept = uabx + uaby;

    reg->auxtoberemoved_axi = (gdouble *)malloc(n * sizeof(gdouble));

    for (i = 0; i < n; i++) {
        reg->auxtoberemoved_axi[i] = axi[i];
    }

    free(axi);
    free(ayi);
    free(bxi);
    free(byi);
    free(xstar);
}

void least_squares_32power_fit_two(GwyDataLine *x, GwyDataLine *y, LinReg *reg) //TODO zkontrolovat teorii
{
    gdouble sx32y, sx3, sy, sx32;
    GwyDataLine *x32;
    gdouble a, b;
    gint i, n;

    x32 = gwy_data_line_new_alike(x, FALSE);
    n = x->res;

    for (i = 0; i < n; i++) {
        //temporary to make it work with MonteCarlo which creates negative forces TODO
        x32->data[i] = ((x->data[i] > 0) ? pow(x->data[i], 1.5) : 0);
    }

    sx32 = 0;
    sx3 = 0;
    sx32y = 0;
    sy = 0;

    for (i = 0; i < n; i++) {
        sx32 += x32->data[i] ;
        sx3 += x->data[i] * x->data[i] * x->data[i];
        sx32y += x32->data[i] * y->data[i];
        sy += y->data[i];
    }

    a = (n * sx32y - sx32 * sy) / (n * sx3 - sx32 * sx32);
    b = (sy - a * sx32) / n;

    reg->slope = a;
    reg->intercept = b;

    g_object_unref(x32);
}

void least_squares_32power_fit(GwyDataLine *x, GwyDataLine *y, gdouble ux, gdouble uy, LinReg *reg) //TODO zkontrolovat teorii
{
    gdouble sx32y, sx3, sx52y, sx4, sxy2;
    gdouble sy2, my;
    GwyDataLine *x32;
    gdouble a, W;
    gint i, n;
    gint p = 1;

    x32 = gwy_data_line_new_alike(x, FALSE);
    n = x->res;

    for (i = 0; i < n; i++) {
        //temporary to make it work with MonteCarlo which creates negative forces TODO
        x32->data[i] = ((x->data[i] > 0) ? pow(x->data[i], 1.5) : 0);
    }

    sx32y = 0;
    sx3 = 0;

    for (i = 0; i < n; i++) {
        sx32y += x32->data[i] * y->data[i];
        sx3 += x->data[i] * x->data[i] * x->data[i];
    }

    a = sx32y / sx3;

    W = 0;

    for (i = 0; i < n; i++) {
        W += (a * x32->data[i] - y->data[i]) * (a * x32->data[i] - y->data[i]);
    }

    if (verbose) {
        g_print("Least squares 32 power fit:\n");
        g_print("  least squares sum %g \n", W);
        g_print("  weighted least squares sum %g \n", W / uy / uy);
    }

    reg->chi = W / n;

    /* find variance of y */
    /* mean value */
    my = gwy_data_line_get_avg(y);

    sy2 = 0;

    for (i = 0; i < n; i++) {
        sy2 += (y->data[i] - my) * (y->data[i] - my);
    }

    /*definice z wikipedie */
    reg->SSres = W / uy / uy; /* squared sum of residuals */
    reg->R2 = 1 - W / sy2; /*R - squared */
    reg->R2adj = 1 - (1 - reg->R2) * (n - 1) / (n - p - 1); /* adjusted R-squared, p- number of parameters fitted, here p=1 */
    reg->chi2 = 1 - reg->R2;

    sx52y = 0;
    sx4 = 0;
    sxy2 = 0;

    for (i = 0; i < n; i++) {
        sx52y += x32->data[i] * x->data[i] * y->data[i];
        sx4 += x->data[i] * x->data[i] * x->data[i] * x->data[i];
        sxy2 += x->data[i] * y->data[i] * y->data[i];
    }

    reg->slope = a;
    reg->intercept = 0;
    reg->unc_slope_y = sqrt(1 / sx3);
    reg->unc_slope_x = 1 / sx3 / sx3 * 3 / 2 * sqrt(sx3 * sx3 * sxy2 - 4 * sx3 * sx52y * sx32y + 4 * sx32y * sx32y * sx4);
    reg->unc_intercept = 0;

    g_object_unref(x32);
}

void nonlinear_least_squares_power_fit(GwyDataLine *x, GwyDataLine *y, gdouble *param)
{
    GwyNLFitFunc ff;
    GwyNLFitter *fitter;
    gboolean fixed[3];
    gpointer user_data;

    fixed[0] = FALSE;
    fixed[1] = TRUE;
    fixed[2] = FALSE;

    ff = (GwyNLFitFunc)fun;

    fitter = gwy_math_nlfit_new(ff, gwy_math_nlfit_derive);

    user_data = NULL; /* avoid warning from using uninitialized pointer */
    gwy_math_nlfit_fit_full(fitter, gwy_data_line_get_res(x), gwy_data_line_get_data_const(x), gwy_data_line_get_data_const(y),
                            NULL, NParam, param, fixed, NULL, user_data);

    if (fitter) {
        gwy_math_nlfit_free(fitter);
    }
}

#ifndef NOFORTRAN
/* old, should be changed for newer */
gint odr_run_fit_old(GwyDataLine *data_x, GwyDataLine *data_y, gdouble beta[3], gint ifix[3], gchar *fnm, gdouble ux, gdouble uy, gdouble std_beta[3], gdouble covariance[3], gdouble *SSres, gdouble *R2, gdouble *R2adj, gdouble *chi2, gboolean rescale)
{
    gint n, m, np, nq;
    gdouble *y, *x;
    gdouble xmax, xmin, ymax, ymin;
    GwyDataLine *rx, *ry;
    gint info;
    gdouble qx, qy;
    gint infofnm;
    gdouble ressum, sxx, syy, xavg, yavg;
    gint p, i;

    if (verbose) {
        g_print("ODR run fit:\n");
    }

    n = gwy_data_line_get_res(data_x);
    m = 1; /* 1-dimensional */
    np = 3; /* # of function parameters */
    nq = 1; /* # of responses per observation */

    p = np; /* # of fitted function parameters */

    gwy_data_line_invert(data_y, TRUE, FALSE);  //// DATA WERE INVERTED FOR SOME REASON...
    y = gwy_data_line_get_data(data_y); // pro nq = 1 to muze byt

    gwy_data_line_invert(data_x, TRUE, FALSE);  //// DATA WERE INVERTED FOR SOME REASON...
    x = gwy_data_line_get_data(data_x); // pro m = 1 to muze byt

    if (rescale) {
        /*rescale data */
        xmax = gwy_data_line_get_max(data_x);
        ymax = gwy_data_line_get_max(data_y);
        xmin = gwy_data_line_get_min(data_x);
        ymin = gwy_data_line_get_min(data_y);

        qx = xmax;
        qy = ymax;

        if (verbose) {
            g_print("  xmax %g xmin %g ymax %g ymin %g \n", xmax, xmin, ymax, ymin);
        }

        rx = gwy_data_line_duplicate(data_x);
        gwy_data_line_multiply(rx, 1. / qx);
        x = gwy_data_line_get_data(rx); // pro m = 1 to muze byt

        ry = gwy_data_line_duplicate(data_y);
        gwy_data_line_multiply(ry, 1. / qy);
        y = gwy_data_line_get_data(ry); // pro nq = 1 to muze byt

        beta[0] *= pow(qx, beta[2]) / qy;
        beta[1] /= qx;

        if (verbose) {
            g_print("  qx %g qy %g \n", qx, qy);
        }

        /*
        for (i=0; i< n; i++)
            printf(" i %d  x %g y %g rx %g ry %g \n", i, data_x->data[i], data_y->data[i], rx->data[i], ry->data[i]);

        printf("rxmax %g rxmin %g rymax %g rymin %g \n", gwy_data_line_get_max(rx), gwy_data_line_get_min(rx), gwy_data_line_get_max(ry), gwy_data_line_get_min(ry));
        for (i=0; i< n; i++)
            printf(" i %d  x %g y %g\n", i, x[i], y[i]);
            */
    }

    //    __odrpack95_wrap_MOD_odr_wrap(&n, &m, &np, &nq, x, y, beta); // GNU COMPILER
    //    odrpack95_wrap_mp_odr_wrap_(&n, &m, &np, &nq, x, y, beta); //INTEL COMPILER

    if (fnm == NULL) {
        infofnm = 0;
    }
    else {
        infofnm = strlen(fnm);
    }

    ressum = 0;

#ifdef _MSC_VER
    ODR_WRAP_POWER_3PARAM(&n, &m, &np, &nq, x, y, beta, ifix, fnm, &infofnm, &ux, &uy, std_beta, covariance, &ressum, &info); // NOT MODULE
#else
    odr_wrap_power_3param_(&n, &m, &np, &nq, x, y, beta, ifix, fnm, &infofnm, &ux, &uy, std_beta, covariance, &ressum, &info); // NOT MODULE
#endif

    p = 3;

    for (i = 0; i < 3 ; i++)
        if (ifix[i] == 0) {
            p--;
        }

    //printf(" p %d \n", p);

    xavg = 0;
    yavg = 0;

    for (i = 0; i < n; i++) {
        xavg += x[i] ;
        yavg += y[i] ;
    }

    xavg /= n;
    yavg /= n;

    sxx = 0;
    syy = 0;

    for (i = 0; i < n; i++) {
        sxx += (x[i] - xavg) * (x[i] - xavg);
        syy += (y[i] - yavg) * (y[i] - yavg);
    }

    /*
    printf(" ressum %g \n", ressum);
    printf(" syy %g yavg %g  \n", syy, yavg);
    printf(" syy/uy/uy %g sxx/ux/ux %g \n", syy/uy/uy, sxx/ux/ux);
    printf(" syy/uy/uy + sxx/ux/ux %g \n", syy/uy/uy+ sxx/ux/ux);
    */
    *SSres = ressum;
    //printf("ressum/(syy/uy/uy +sxx/ux/ux) %g \n", ressum/(syy/uy/uy + sxx/ux/ux));
    //    *R2 = 1 - ressum / (syy / uy / uy + sxx / ux / ux);
    *R2adj = 1 - (1 - *R2) * (n - 1) / (n - p - 1);
    *chi2 = 1 - (*R2);

    if (rescale) {
        g_object_unref(rx);
        g_object_unref(ry);

        beta[0] /= pow(qx, beta[2]) / qy;
        beta[1] *= qx;

        std_beta[0] /= pow(qx, beta[2]) / qy;
        std_beta[1] *= qx;

        covariance[0] *= qx;
        covariance[1] /= pow(qx, beta[2]) / qy;
        covariance[2] /= pow(qx, beta[2]) / qy;
        covariance[2] *= qx;
    }

    return info;
}
#endif

#ifndef NOFORTRAN
gint odr_run_fit(GwyDataLine *data_x, GwyDataLine *data_y, const gint  np, gdouble *beta, const gint *ifix, gchar *fnm, gdouble ux, gdouble uy, gdouble *std_beta, gdouble *covariance, gdouble *SSres, gdouble *R2, gdouble *R2adj, gdouble *chi2, gint what_to_fit)
{
    gint n, m, nq;
    gdouble *y, *x;
    gint info;
    gint infofnm;
    gdouble ressum, sxx, syy, xavg, yavg;
    gint p, i;

    n = gwy_data_line_get_res(data_x);
    m = 1; /* 1-dimensional */
    nq = 1; /* # of responses per observation */

    p = np; /* # of fitted function parameters */

    gwy_data_line_invert(data_y, TRUE, FALSE);  //// DATA WERE INVERTED FOR SOME REASON...
    y = gwy_data_line_get_data(data_y); // pro nq = 1 to muze byt

    gwy_data_line_invert(data_x, TRUE, FALSE);  //// DATA WERE INVERTED FOR SOME REASON...
    x = gwy_data_line_get_data(data_x); // pro m = 1 to muze byt

    if (fnm == NULL) {
        infofnm = 0;
    }
    else {
        infofnm = strlen(fnm);
    }

    ressum = 0;

    if (verbose) {
        g_print("in odrwrap fit-utils np %d \n", np);
        g_print("ux %g uy %g \n", ux, uy);
    }

#ifdef _MSC_VER

    switch (what_to_fit) {
    case POWER_3PARAM:
        ODR_WRAP_POWER_3PARAM(&n, &m, &np, &nq, x, y, beta, ifix, fnm, &infofnm, &ux, &uy, std_beta,  covariance, &ressum, &info); // NOT MODULE
        break;

    case HERTZ_RADIAL:
        ODR_WRAP_HERTZ_RADIAL(&n, &m, &np, &nq, x, y, beta, ifix, fnm, &infofnm, &ux, &uy, std_beta,  covariance, &ressum, &info); // NOT MODULE
        break;

    case POWER_2PARAM:
        ODR_WRAP_POWER_2PARAM(&n, &m, &np, &nq, x, y, beta, ifix, fnm, &infofnm, &ux, &uy, std_beta,  covariance, &ressum, &info); // NOT MODULE
        break;

    case LINEAR:
        ODR_WRAP_LINEAR(&n, &m, &np, &nq, x, y, beta, ifix, fnm, &infofnm, &ux, &uy, std_beta,  covariance, &ressum, &info); // NOT MODULE
        break;

    default:
        g_printerr("Should not get here. \n");
        break;
    }

#else

    switch (what_to_fit) {
    case POWER_3PARAM:
        odr_wrap_power_3param_(&n, &m, &np, &nq, x, y, beta, ifix, fnm, &infofnm, &ux, &uy, std_beta, covariance, &ressum, &info); // NOT MODULE
        break;

    case HERTZ_RADIAL:
        odr_wrap_hertz_radial_(&n, &m, &np, &nq, x, y, beta, ifix, fnm, &infofnm, &ux, &uy, std_beta, covariance, &ressum, &info); // NOT MODULE
        break;

    case POWER_2PARAM:
        odr_wrap_power_2param_(&n, &m, &np, &nq, x, y, beta, ifix, fnm, &infofnm, &ux, &uy, std_beta, covariance, &ressum, &info); // NOT MODULE
        break;

    case LINEAR:
        odr_wrap_linear_(&n, &m, &np, &nq, x, y, beta, ifix, fnm, &infofnm, &ux, &uy, std_beta, covariance, &ressum, &info); // NOT MODULE
        break;

    default:
        g_printerr("Should not get here. \n");
        break;
    }

#endif

    p = np;

    for (i = 0; i < np ; i++)
        if (ifix[i] == 0) {
            p--;
        }

    xavg = 0;
    yavg = 0;

    for (i = 0; i < n; i++) {
        xavg += x[i] ;
        yavg += y[i] ;
    }

    xavg /= n;
    yavg /= n;

    sxx = 0;
    syy = 0;

    for (i = 0; i < n; i++) {
        sxx += (x[i] - xavg) * (x[i] - xavg);
        syy += (y[i] - yavg) * (y[i] - yavg);
    }

    *SSres = ressum;
    *R2 = 1 - ressum / (syy / uy / uy + sxx / ux / ux);
    *R2adj = 1 - (1 - *R2) * (n - 1) / (n - p - 1);
    *chi2 = 1 - (*R2);

    return info;
}
#endif

/* static functions */

static gdouble fun(gdouble x, G_GNUC_UNUSED gint nparam, gdouble *param, G_GNUC_UNUSED gpointer user_data, gboolean *fres)
{
    *fres = TRUE;

    // param[0]  a
    // param[1]  hp
    // param[2]  m

    return  param[0] * pow(x - param[1], param[2]);
}

gint filter_negative_data(GwyDataLine *x, GwyDataLine *y)
{
    gint n, i;
    gint nskip = 0;
    gdouble *p, *q, *pn, *qn;

    n = x->res;

    p = x->data;
    pn = x->data;
    q = y->data;
    qn = y->data;

    for (i = 0; i < n; i++) {
        if (*p <= 0) {
            p++;
            q++;
            nskip++;
        }
        else {
            *pn = *p;
            *qn = *q;
            p++;
            q++;
            pn++;
            qn++;
        }
    }

    gwy_data_line_resize(x, 0, n - nskip);
    gwy_data_line_resize(y, 0, n - nskip);

    return nskip;
}

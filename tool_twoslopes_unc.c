#ifndef NOFORTRAN

#include "tool_twoslopes_unc.h"
#include "tool_twoslopes.h"

#include "datatypes.h"
#include "fddata.h"
#include "fit-utils.h"
#include "file-utils.h"
#include "niget-common.h"
#include "mc-utils.h"

#include <math.h>
#include <string.h>


static const gchar *logfnmload = "fit.log.slopes.load.unc";
static const gchar *logfnmunload = "fit.log.slopes.unload.unc";

//static const gchar *mclogfnmload = "fit.log.slopes.load.mc";
//static const gchar *mclogfnmunload = "fit.log.slopes.unload.mc";

static void init_shifted_Slopesdata(const Slopesdata *src, Slopesdata *dest, gdouble dh);
static gboolean nparam_to_string(gint i, gchar *str, gint nstr);

void init_SlopesUncdata(SlopesUncdata *slopesunc, Instdata *instdata)
{
    slopesunc->instdata = *instdata;

    slopesunc->Nmc = N_MC;

    slopesunc->umuh = 0;
    slopesunc->umuF = 0;
    slopesunc->unuh = 0;
    slopesunc->unuF = 0;

    slopesunc->uSuuh = 0;
    slopesunc->uSuuF = 0;
    slopesunc->uSluh = 0;
    slopesunc->uSluF = 0;

    slopesunc->uAuh = 0;
    slopesunc->uAuF = 0;
    slopesunc->uA = 0;

    slopesunc->uHituh = 0;
    slopesunc->uHituF = 0;
    slopesunc->uHit = 0;

    slopesunc->uEruh = 0;
    slopesunc->uEruF = 0;
    slopesunc->uEr = 0;

    slopesunc->uEituh = 0;
    slopesunc->uEituF = 0;
    slopesunc->uEitunu = 0;
    slopesunc->uEitunui = 0;
    slopesunc->uEituEi = 0;
    slopesunc->uEit = 0;

    slopesunc->Ec = NULL;
    slopesunc->Hc = NULL;
    slopesunc->Ac = NULL;
}

static void init_shifted_Slopesdata(const Slopesdata *src, Slopesdata *dest, gdouble dh)
{
    init_Slopesdata(dest);

    // shift fitting range  in h-mode by dh
    // don't shift range for F-mode

    dest->Finputload = src->Finputload;

    dest->loadrange[0] = src->loadrange[0] - dh;
    dest->loadrange[1] = src->loadrange[1] - dh;
    dest->loadfrom = dest->loadrange[0];
    dest->loadto = dest->loadrange[1];

    dest->loadrange_pct_Fmax[0] = src->loadrange_pct_Fmax[0];
    dest->loadrange_pct_Fmax[1] = src->loadrange_pct_Fmax[1];
    dest->loadfrom_pct_Fmax = src->loadfrom_pct_Fmax;
    dest->loadto_pct_Fmax = src->loadto_pct_Fmax;

    dest->loadifixb[0] = src->loadifixb[0];
    dest->loadifixb[1] = src->loadifixb[1];
    dest->loadifixb[2] = src->loadifixb[2];

    dest->Finputunload = src->Finputunload;

    dest->unloadrange[0] = src->unloadrange[0] - dh;
    dest->unloadrange[1] = src->unloadrange[1] - dh;
    dest->unloadfrom = dest->unloadrange[0];
    dest->unloadto = dest->unloadrange[1];

    dest->unloadrange_pct_Fmax[0] = src->unloadrange_pct_Fmax[0];
    dest->unloadrange_pct_Fmax[1] = src->unloadrange_pct_Fmax[1];
    dest->unloadfrom_pct_Fmax = src->unloadfrom_pct_Fmax;
    dest->unloadto_pct_Fmax = src->unloadto_pct_Fmax;

    g_free(dest->logfnmload);
    g_free(dest->logfnmunload);
    dest->logfnmload = create_log_filename_path(logfnmload);
    dest->logfnmunload = create_log_filename_path(logfnmunload);

    // copy beta
    dest->beta = src->beta;
}

void slopes_propagate_uncertainties_calc(const Slopesdata *sldata, SlopesUncdata *unc, const FDdata *fd)
{
    gdouble uh, uF;
    gdouble mh, mF, epsh, epsF ;
    gdouble umh, umF, uhph, uhpF, um, uhp;
    gdouble unh, unF, uh0h, uh0F, un, uh0;
    gdouble uSuh, uSuF, uSu;
    gdouble uSlh, uSlF, uSl;
    gdouble  depsdm, depsdmh, depsdmF;
    gdouble uepsh, uepsF;
    gdouble ueps;
    gdouble beta[3], std_beta[3], covariance[3];
    gint istart, iend, ndata;
    GwyDataLine *x, *y;
    gdouble ux, uy;
    gdouble n, Sl, h0, m, Su, hp, eps, hmax, Fmax, Aphc, Hit, Er, Eit, nu, nui, Ei;
    gdouble uEi, unu, unui;
    gdouble covnh0h, covnh0F, covnh0;
    gdouble covmhph, covmhpF, covmhp;
    gdouble K, uKuh, uKuF, uK;
    gdouble covKSuh, covKSuF, covKSu;
    gdouble covSumh, covSumF, covSum;
    gdouble slbeta;
    gchar *unc_logfnmload, *unc_logfnmunload;
    gdouble SSres, R2, R2adj, chi2;
    gint np;

    //make copies to simplify code
    uh = unc->instdata.sigma_h;
    uF = unc->instdata.sigma_F;
    n = sldata->n;
    h0 = sldata->h0;
    Sl = sldata->Sload;
    m = sldata->m;
    hp = sldata->hp;
    Su = sldata->Sunload;
    eps = sldata->eps;
    slbeta = sldata->beta;
    hmax = fd->hmax;
    Fmax = fd->Fmax;
    Aphc = sldata->Aphc;
    Hit = sldata->Hit;
    Eit = sldata->Eit;
    Er = sldata->Er;
    Ei = unc->instdata.Ei;
    nui = unc->instdata.nui;
    nu = unc->instdata.nu;
    uEi = unc->instdata.uEi;
    unui = unc->instdata.unui;
    unu = unc->instdata.unu;

    unc_logfnmload =  create_log_filename_path(logfnmload);
    unc_logfnmunload =  create_log_filename_path(logfnmunload);

    np = 3;

    if (VERBOSE) {
        printf(" \n\nPropagate uncertainties\n");
    }

    /*   u(n)    */

    /* prepare odr fit */

    beta[0] = sldata->gamma;
    beta[1] = sldata->h0;
    beta[2] = sldata->n;

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (sldata->Finputload) {
        range_to_indices(sldata->loadrange_pct_Fmax[0] * 0.01 * fd->Fmax, sldata->loadrange_pct_Fmax[1] * 0.01 * fd->Fmax, fd->Fload, FALSE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(sldata->loadrange[0], sldata->loadrange[1], fd->hload, FALSE, &istart, &iend, &ndata);
    }

    x = gwy_data_line_part_extract(fd->hload, istart, ndata);
    y = gwy_data_line_part_extract(fd->Fload, istart, ndata);

    /* fit with only uh */

    unc_logfnmload =  create_log_filename_path(logfnmload);
    unc_logfnmunload =  create_log_filename_path(logfnmunload);

    np = 3;

    if (VERBOSE) {
        printf(" \n\nPropagate uncertainties\n");
    }

    /*   u(n)    */

    /* prepare odr fit */

    beta[0] = sldata->gamma;
    beta[1] = sldata->h0;
    beta[2] = sldata->n;

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (sldata->Finputload) {
        range_to_indices(sldata->loadrange_pct_Fmax[0] * 0.01 * fd->Fmax, sldata->loadrange_pct_Fmax[1] * 0.01 * fd->Fmax, fd->Fload, FALSE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(sldata->loadrange[0], sldata->loadrange[1], fd->hload, FALSE, &istart, &iend, &ndata);
    }

    x = gwy_data_line_part_extract(fd->hload, istart, ndata);
    y = gwy_data_line_part_extract(fd->Fload, istart, ndata);

    /* fit with only uh */
    ux = uh;
    uy = 0;
    beta[0] = sldata->gamma;
    beta[1] = sldata->h0;
    beta[2] = sldata->n;
    odr_run_fit(x, y, np, beta, sldata->loadifixb, unc_logfnmload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);
    unh = std_beta[2];
    uh0h = std_beta[1];
    covnh0h = covariance[0];

    if (verbose)
        g_print(" fit with only uh:\n  gamma %g n %g h0 %g \n  u(n) %g u(h0) %g corr(n,h0) %g \n",
                beta[0], beta[1], beta[2], unh, uh0h, covnh0h / unh / uh0h);

    /* fit with only uF */
    ux = 0;
    uy = uF;
    odr_run_fit(x, y, np, beta, sldata->loadifixb, unc_logfnmload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);
    unF = std_beta[2];
    uh0F = std_beta[1];
    covnh0F = covariance[0];

    if (verbose)
        g_print(" fit with only uF:\n  gamma %g n %g h0 %g \n  u(n) %g u(h0) %g corr(n,h0) %g \n",
                beta[0], beta[1], beta[2], unF, uh0F, covnh0F / unF / uh0F);

    unc->unuh = unh;
    unc->unuF = unF;

    /*   u(m)    */

    /* prepare odr fit */

    beta[0] = sldata->alpha;
    beta[1] = sldata->hp;
    beta[2] = sldata->m;

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (sldata->Finputunload) {
        range_to_indices(sldata->unloadrange_pct_Fmax[0] * 0.01 * fd->Fmax, sldata->unloadrange_pct_Fmax[1] * 0.01 * fd->Fmax, fd->Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(sldata->unloadrange[0], sldata->unloadrange[1], fd->hunload, TRUE, &istart, &iend, &ndata);
    }

    x = gwy_data_line_part_extract(fd->hunload, istart, ndata);
    y = gwy_data_line_part_extract(fd->Funload, istart, ndata);

    /* fit with only uh */
    ux = uh;
    uy = 0;
    beta[0] = sldata->alpha;
    beta[1] = sldata->hp;
    beta[2] = sldata->m;
    odr_run_fit(x, y, np, beta, sldata->unloadifixb, unc_logfnmunload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);
    mh = beta[2];
    umh = std_beta[2];
    uhph = std_beta[1];
    covmhph = covariance[0];
    unc->umuh = umh;

    if (verbose)
        g_print(" fit with only uh:\n  alpha %g m %g hp %g \n  u(m) %g u(hp) %g corr(m,hp) %g \n",
                beta[0], beta[1], beta[2], umh, uhph, covmhph / umh / uhph);

    /*   u(eps)    */
    epsh = epsilon(mh);
    depsdmh = der_epsilon(mh);
    uepsh = depsdmh * umh;

    if (verbose) {
        g_print(" uh fit \n");
        g_print("depsdm %g   \n", depsdmh);
        g_print(" u(eps;h) %g urel(eps;h)  %g \n", uepsh, uepsh / epsh);
    }

    /* fit with only uF */
    ux = 0;
    uy = uF;
    odr_run_fit(x, y, np, beta, sldata->unloadifixb, unc_logfnmunload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);
    mF = beta[2];
    umF = std_beta[2];
    uhpF = std_beta[1];
    covmhpF = covariance[0];

    if (verbose)
        g_print(" fit with only uF:\n  alpha %g m %g hp %g \n  u(m) %g u(hp) %g corr(m,hp) %g \n",
                beta[0], beta[1], beta[2], umF, uhpF, covmhpF / umF / uhpF);

    unc->umuF = umF;

    /*   u(eps)    */
    epsF = epsilon(mF);
    depsdmF = der_epsilon(mF);
    uepsF = depsdmF * umF;

    if (verbose) {
        g_print(" uF fit \n");
        g_print("depsdm %g   \n", depsdmF);
        g_print(" u(eps;F)  %g urel(eps;F) %g \n", uepsF, uepsF / epsF);
    }

    /*    u(Sl)   */
    uSlh = sqrt(unh * unh / n / n +
                uh * uh / (hmax - h0) / (hmax - h0) +
                uh0h * uh0h / (hmax - h0) / (hmax - h0) +
                2 * covnh0h / n / (hmax - h0)
               ) * Sl;

    uSlF = sqrt(unF * unF / n / n +
                uF * uF / Fmax / Fmax  +
                uh0F * uh0F / (hmax - h0) / (hmax - h0) +
                2 * covnh0F / n / (hmax - h0)
               ) * Sl;
    unc->uSluh = uSlh;
    unc->uSluF = uSlF;
    /* u(Sl;h) = Sl *sqrt [ (unh/n)^2 + (uh/(hmax-h0))^2 + (uh0h/(hmax-h0))^2  + 2 cov(n,h0;h)/(n*(hmax-h0))]
     * u(Sl;F) = Sl *sqrt [ (unF/n)^2 + (uFmax/Fmax)^2   + (uh0F/(hmax-h0))^2  + 2 cov(n,h0;F)/(n*(hmax-h0))]
     */

    if (verbose) {
        g_print(" u(Sl;h) %g u(Sl;F) %g  urel(Sl;h) %g urel(Sl;F) %g \n", uSlh, uSlF, uSlh / Sl, uSlF / Sl);
    }

    /*    u(Su)   */
    uSuh = sqrt(umh * umh / m / m +
                uh * uh / (hmax - hp) / (hmax - hp) +
                uhph * uhph / (hmax - hp) / (hmax - hp) +
                2 * covmhph / m / (hmax - hp)
               ) * Su;

    uSuF = sqrt(umF * umF / m / m +
                uF * uF / Fmax / Fmax  +
                uhpF * uhpF / (hmax - hp) / (hmax - hp) +
                2 * covmhpF / m / (hmax - hp)
               ) * Su;

    unc->uSuuh = uSuh;
    unc->uSuuF = uSuF;
    /* u(Su;h) = Su *sqrt [ (umh/m)^2 + (uh/(hmax-hp))^2 + (uhph/(hmax-hp))^2  + 2 cov(m,hp;h)/(m*(hmax-hp))]
     * u(Su;F) = Su *sqrt [ (umF/m)^2 + (uFmax/Fmax)^2   + (uhpF/(hmax-hp))^2  + 2 cov(m,hp;F)/(m*(hmax-hp))]
     */

    covSumh =  m * Fmax / (hmax - hp) / (hmax - hp) * covmhph + Fmax / (hmax - hp) * umh * umh;
    covSumF =  m * Fmax / (hmax - hp) / (hmax - hp) * covmhpF + Fmax / (hmax - hp) * umF * umF;

    /* cov(Su,m;h) = m Fmax/(hmax-hp)^2  cov(hp,m;h) + Fmax/(hmax-hp) u(m;h)^2
     * cov(Su,m;F) = m Fmax/(hmax-hp)^2  cov(hp,m;F) + Fmax/(hmax-hp) u(m;F)^2
     * */

    if (verbose) {
        g_print(" u(Su;h) %g u(Su;F) %g  urel(Su;h) %g urel(Su;F) %g \n", uSuh, uSuF, uSuh / Su, uSuF / Su);
        g_print(" cov(Su,m;h) %g cov(Su,m;F) %g \n", covSumh, covSumF);
    }

    /*    u(K)   */
    K = (2 * Su - Sl * eps * slbeta) / Su / Sl;

    uKuh = sqrt(4 * uSlh * uSlh / Sl / Sl / Sl / Sl + eps * slbeta * eps * slbeta / Su / Su / Su / Su * uSuh * uSuh + slbeta * slbeta / Su / Su * uepsh * uepsh + slbeta * eps / Su / Su * depsdmh * covSumh);
    uKuF = sqrt(4 * uSlF * uSlF / Sl / Sl / Sl / Sl + eps * slbeta * eps * slbeta / Su / Su / Su / Su * uSuF * uSuF + slbeta * slbeta / Su / Su * uepsF * uepsF + slbeta * eps / Su / Su * depsdmF * covSumF);

    covKSuh = eps * slbeta / Su / Su * uSuh * uSuh - slbeta / Su * depsdmh * covSumh;
    covKSuF = eps * slbeta / Su / Su * uSuh * uSuF - slbeta / Su * depsdmF * covSumF;
    /*  K = 2 /Sl - eps*BETA/Su
     *  u(K;h) = sqrt [ (-2/S1^2 u(S1;h))^2 + ( eps*BETA/Su^2 u(Su;h))^2  + ( BETA/Su u(eps;h))^2 + BETA eps/Su^2 depsdm(h) cov(S,m;h) ]
     *  u(K;F) = sqrt [ (-2/S1^2 u(S1;F))^2 + ( eps*BETA/Su^2 u(Su;F))^2  + ( BETA/Su u(eps;F))^2 + BETA eps/Su^2 depsdm(F) cov(S,m;F) ]
     *  cov(K,Su; h) = eps*BETA/Su^2 u(Su;h)^2 - BETA/Su deps/dm cov(Su,m;h)
     *  cov(K,Su; F) = eps*BETA/Su^2 u(Su;F)^2 - BETA/Su deps/dm cov(Su,m;F)
     *  */

    /*    u(A)   */
    unc->uAuh = 2 * Aphc * uKuh / K;
    unc->uAuF = 2 * Aphc * sqrt(uKuF * uKuF / K / K + uF * uF / Fmax / Fmax);
    /* u(A; h) =  2A *(uK;h)/K ]
       u(A; F) =  2A * sqrt [ (uFmax/Fmax)^2 + ((uK;F)/K)^2 ] */

    if (verbose) {
        g_print(" u(A;h) %g u(A;F) %g  urel(A;h) %g urel(A;F) %g \n", unc->uAuh, unc->uAuF, unc->uAuh / Aphc, unc->uAuF / Aphc);
    }

    /*     u(H_IT)   */
    unc->uHituh = 2 * Hit * uKuh / K;
    unc->uHituF = Hit * sqrt(uF * uF / Fmax / Fmax + 4 * uKuF * uKuF / K / K);
    /* u(H; h ) = 2 H_IT u(K;h)/ K
     * u(H; F)  = H_IT sqrt [ (uFmax/Fmax)^2 + (2 u(K;F)/K)^2 ] */

    if (verbose) {
        g_print(" u(Hit;h) %g u(Hit;F) %g  urel(Hit;h) %g urel(Hit;F) %g \n", unc->uHituh, unc->uHituF, unc->uHituh / Hit, unc->uHituF / Hit);
    }

    /*     u(E_r)   */
    unc->uEruh = Er * sqrt(uKuh * uKuh / K / K + uSuh * uSuh / Su / Su - covKSuh / K / Su);
    unc->uEruF = Er * sqrt(uF * uF / Fmax / Fmax + uKuF * uKuF / K / K + uSuF * uSuF / Su / Su - covKSuF / K / Su);
    /* u(Er;h) = Er sqrt [  [u(K;h)/ K]^2 + [u(Su;h)/Su]^2 - cov(K,Su;h)/K/Su ]
     * u(Er;F) = Er sqrt [ (uFmax/Fmax)^2 +  [u(K;h)/ K]^2 + [u(Su;h)/Su]^2 - cov(K,Su;F)/K/Su] */

    if (verbose) {
        g_print(" u(Er;h) %g u(Er;F) %g  urel(Er;h) %g urel(Er;F) %g \n", unc->uEruh, unc->uEruF, unc->uEruh / Er, unc->uEruF / Er);
    }

    /*     u(E_IT)   */
    unc->uEituh = Eit * Eit / Er / Er / (1 - nu * nu) * unc->uEruh; // Eit^2/Er^2 /(1-nu^2) u(Er;h)
    unc->uEituF = Eit * Eit / Er / Er / (1 - nu * nu) * unc->uEruF; // Eit^2/Er^2 /(1-nu^2) u(Er;F)

    if (verbose) {
        g_print("u(Eit;h) %g  u(Eit;F) %g  urel(Eit;h) %g  urel(Eit;F) %g\n", unc->uEituh, unc->uEituF, unc->uEituh / Eit, unc->uEituF / Eit);
    }

    // total uncertainties
    // should take um as given by odr_wrap, otherwise doublecounting
    /* fit with both uh and uF */
    ux = uh;
    uy = uF;
    beta[0] = sldata->gamma;
    beta[1] = sldata->h0;
    beta[2] = sldata->n;

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (sldata->Finputload) {
        range_to_indices(sldata->loadrange_pct_Fmax[0] * 0.01 * fd->Fmax, sldata->loadrange_pct_Fmax[1] * 0.01 * fd->Fmax, fd->Fload, FALSE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(sldata->loadrange[0], sldata->loadrange[1], fd->hload, FALSE, &istart, &iend, &ndata);
    }

    x = gwy_data_line_part_extract(fd->hload, istart, ndata);
    y = gwy_data_line_part_extract(fd->Fload, istart, ndata);

    odr_run_fit(x, y, np, beta, sldata->loadifixb, unc_logfnmload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);
    un = std_beta[2];
    uh0 = std_beta[1];
    covnh0 = covariance[0];

    unc->un = un;

    if (verbose)
        g_print(" fit with both uh and uF:\n  gamma %g h0 %g n %g \n  u(n) %g u(h0) %g corr(n,h0) %g \n",
                beta[0], beta[1], beta[2], un, uh0, covnh0 / un / uh0);

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (sldata->Finputunload) {
        range_to_indices(sldata->unloadrange_pct_Fmax[0] * 0.01 * fd->Fmax, sldata->unloadrange_pct_Fmax[1] * 0.01 * fd->Fmax, fd->Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(sldata->unloadrange[0], sldata->unloadrange[1], fd->hunload, TRUE, &istart, &iend, &ndata);
    }

    x = gwy_data_line_part_extract(fd->hunload, istart, ndata);
    y = gwy_data_line_part_extract(fd->Funload, istart, ndata);

    beta[0] = sldata->alpha;
    beta[1] = sldata->hp;
    beta[2] = sldata->m;
    odr_run_fit(x, y, np, beta, sldata->unloadifixb, unc_logfnmunload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);
    um = std_beta[2];
    uhp = std_beta[1];
    covmhp = covariance[0];

    unc->um = um;

    if (verbose)
        g_print(" fit with both uh and uF:\n  alpha %g hp %g m %g \n  u(m) %g u(hp) %g corr(m,hp) %g \n",
                beta[0], beta[1], beta[2], um, uhp, covmhp / um / uhp);

    /*   u(eps)    */
    depsdm = der_epsilon(m);
    ueps = depsdm * um;

    if (verbose) {
        g_print("depsdm %g  \n", depsdm);
        g_print(" u(eps) %g    urel(eps) %g  \n", ueps, ueps / eps);
    }

    /*    u(Sl)   */
    uSl = sqrt(un * un / n / n  + uh * uh / (hmax - h0) / (hmax - h0) +  uF * uF / Fmax / Fmax +
               uh0 * uh0 / (hmax - h0) / (hmax - h0) + 2 * covnh0 / n / (hmax - h0)) * Sl;

    unc->uSl = uSl;
    /* u(Sl) = Sl *sqrt [ (umh/m)^2 +  (uhmax/(hmax-hp))^2 + (uFmax/Fmax)^2 +  (uh0/(hmax-h0))^2 + 2 cov(m,h0)/(n*(hmax-h0)) ]
     */

    if (verbose) {
        g_print(" u(Sl) %g    urel(Sl) %g \n", uSl, uSl / Sl);
    }

    /*    u(Su)   */
    uSu = sqrt(um * um / m / m  + uh * uh / (hmax - hp) / (hmax - hp) +  uF * uF / Fmax / Fmax +
               uhp * uhp / (hmax - hp) / (hmax - hp) + 2 * covmhp / m / (hmax - hp)) * Su;
    unc->uSu = uSu;
    /* u(Su) = Su *sqrt [ (umh/m)^2 +  (uhmax/(hmax-hp))^2 + (uFmax/Fmax)^2 +  (uhp/(hmax-hp))^2 + 2 cov(m,hp)/(m*(hmax-hp)) ]
     */

    covSum = m * Fmax / (hmax - hp) / (hmax - hp) * covmhp + Fmax / (hmax - hp) * um * um;

    /* cov(Su,m) = m Fmax/(hmax-hp)^2  cov(hp,m) + Fmax/(hmax-hp) u(m)^2
     * */

    if (verbose) {
        g_print(" u(Su) %g    urel(Su) %g   cov(Su,m) %g \n", uSu, uSu / Su, covSum);
    }

    /*    u(K)   */
    K = (2 * Su - Sl * eps * slbeta) / Su / Sl;

    uK = sqrt(4 * uSl * uSl / Sl / Sl / Sl / Sl + eps * slbeta * eps * slbeta / Su / Su / Su / Su * uSu * uSu + slbeta * slbeta / Su / Su * ueps * ueps + slbeta * eps / Su / Su * depsdm * covSum);

    covKSu = eps * slbeta / Su / Su * uSu * uSu - slbeta / Su * depsdm * covSum;

    /*  K = 2 /Sl - eps*BETA/Su
     *  u(K) = sqrt [ (-2/S1^2 u(S1))^2 + ( eps*BETA/Su^2 u(Su)^2  + ( BETA/Su u(eps))^2 + BETA eps/Su^2 depsdm cov(S,m)  ]
     *  cov(K,Su) = eps*BETA/Su^2 u(Su)^2  -BETA/Su deps/dm cov(Su,m)  */

    if (verbose) {
        g_print(" u(K) %g cov(K,Su) %g    urel(K) %g  \n", uK, covKSu, uK / K);
    }

    /*    u(A)   */
    unc->uA = 2 * Aphc * sqrt(uK * uK / K / K + uF * uF / Fmax / Fmax);

    if (verbose) {
        g_print(" u(A) %g    urel(A) %g  \n", unc->uA, unc->uA / Aphc);
    }

    /*     u(H_IT)   */
    unc->uHit = Hit * sqrt(uF * uF / Fmax / Fmax + 4 * uK * uK / K / K);

    if (verbose) {
        g_print(" u(Hit)  %g   urel(Hit) %g \n", unc->uHit, unc->uHit / Hit);
    }

    /*     u(E_r)   */
    unc->uEr = Er * sqrt(uF * uF / Fmax / Fmax + uK * uK / K / K + uSu * uSu / Su / Su - covKSu / K / Su);

    if (verbose) {
        g_print(" u(Er) %g   urel(Er) %g  \n", unc->uEr, unc->uEr / Er);
    }

    /*     u(E_IT)   */
    //contribution from h,F
    unc->uEit = Eit * Eit / Er / Er / (1 - nu * nu) * unc->uEr; // Eit^2/Er^2 /(1-nu^2) uEr

    //contribution from nu, nui, Ei
    unc->uEitunu =  fabs(2 * nu / (1 - nu * nu) * Eit * unu); // -2nu/(1-nu^2) Eit u(nu)
    unc->uEitunui = fabs(2 * nui / (1 - nu * nu) * Eit * Eit / Ei * unui) * 1e9; // -2nui/(1-nu^2) Eit^2/Ei u(nui)
    unc->uEituEi =  fabs((1 - nui * nui) / (1 - nu * nu) * Eit * Eit / Ei / Ei * uEi) * 1e9; //-(1-nui^2)/(1-nu)^2 Eit^2/Ei^2 u(Ei)

    if (verbose) {
        g_print("u(Eit;h,F) %g    urel(Eit;hF) %g  \n", unc->uEit, unc->uEit / Eit);
        g_print("u(Eit;nu) %g  u(Eit;nui) %g  u(Eit;Ei) %g  urel(Eit;nu) %g  urel(Eit;nui) %g  urel(Eit;Ei) %g\n",
                unc->uEitunu, unc->uEitunui, unc->uEituEi, unc->uEitunu / Eit, unc->uEitunui / Eit, unc->uEituEi / Eit);
    }

    unc->uEit = sqrt(unc->uEit * unc->uEit  + unc->uEitunu * unc->uEitunu + unc->uEitunui * unc->uEitunui + unc->uEituEi * unc->uEituEi);

    if (verbose) {
        g_print("u(Eit;hFnuetc) %g    urel(Eit; hF_nu_etc) %g \n", unc->uEit, unc->uEit / Eit);
    }

    if (verbose) {
        g_print(" \nEnd propagation of uncertainties\n\n");
    }

    g_object_unref(x);
    g_object_unref(y);
    g_free(unc_logfnmload);
    g_free(unc_logfnmunload);
}

gchar *slopes_uncertainty_export_data(const Slopesdata *sldata, const SlopesUncdata *slunc, const FDdata *fddata, enum ExportFormat exportformat)
{
    gint j;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Parameters of two slopes analysis (ODR fit) \n");
    g_string_append_printf(buf, "# u(depth): %g nm\n", slunc->instdata.sigma_h);
    g_string_append_printf(buf, "# u(load): %g nm\n", slunc->instdata.sigma_F);
    g_string_append_printf(buf, "#beta: %g \n", sldata->beta);
    g_string_append_printf(buf, "#nu: %g                 u(nu):  %g  \n", slunc->instdata.nu, slunc->instdata.unu);
    g_string_append_printf(buf, "#nu_indenter: %g       u(nu_indenter):  %g   \n", slunc->instdata.nui, slunc->instdata.unui);
    g_string_append_printf(buf, "#E_indenter: %g        u(E_indenter):   %g GPa\n", slunc->instdata.Ei * 1e-9, slunc->instdata.uEi * 1e-9); //TODO units
    g_string_append_printf(buf, "#Fitted loading range: ");

    if (sldata->Finputload) {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", sldata->loadrange_pct_Fmax[0], sldata->loadrange_pct_Fmax[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] nm \n", sldata->loadrange[0], sldata->loadrange[1]);
    }

    g_string_append_printf(buf, "# corresponds : ");

    if (sldata->Finputload) {
        g_string_append_printf(buf, "[%g , %g ] nm \n", sldata->loadrange[0], sldata->loadrange[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", sldata->loadrange_pct_Fmax[0], sldata->loadrange_pct_Fmax[1]);
    }

    g_string_append_printf(buf, "# %d datapoints used.\n", sldata->nfitloaddata);
    g_string_append_printf(buf, "#Fitted unloading range: ");

    if (sldata->Finputunload) {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", sldata->unloadrange_pct_Fmax[0], sldata->unloadrange_pct_Fmax[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] nm \n", sldata->unloadrange[0], sldata->unloadrange[1]);
    }

    g_string_append_printf(buf, "# corresponds : ");

    if (sldata->Finputunload) {
        g_string_append_printf(buf, "[%g , %g ] nm \n", sldata->unloadrange[0], sldata->unloadrange[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", sldata->unloadrange_pct_Fmax[0], sldata->unloadrange_pct_Fmax[1]);
    }

    g_string_append_printf(buf, "# %d datapoints used.\n", sldata->nfitunloaddata);

    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Results of two slopes analysis (ODR fit) \n");
    g_string_append_printf(buf, "#Ap(hc): %g nm^2\n", sldata->Aphc);
    g_string_append_printf(buf, "#Er: %g GPa\n", sldata->Er);
    g_string_append_printf(buf, "#Hit: %g GPa\n", sldata->Hit);
    g_string_append_printf(buf, "#Eit: %g GPa\n", sldata->Eit);
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Results of the Gaussian propagation of uncertainties \n");
    g_string_append_printf(buf, "#u(n) =  %g  \n", slunc->un);
    g_string_append_printf(buf, "#u(Sload) =  %g mN/nm \n", slunc->uSl);
    g_string_append_printf(buf, "#u(m) =  %g  \n", slunc->um);
    g_string_append_printf(buf, "#u(Sunload) =  %g mN/nm \n", slunc->uSu);
    g_string_append_printf(buf, "#u(Ap(hc)) =  %g nm^2 \n", slunc->uA);
    g_string_append_printf(buf, "#u(Hit) =  %g MPa \n", slunc->uHit);
    g_string_append_printf(buf, "#u(Er) =  %g GPa \n", slunc->uEr);
    g_string_append_printf(buf, "#u(Eit) =  %g GPa \n", slunc->uEit);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of n\n");
    g_string_append_printf(buf, "#u(n;uh) = %g \n", slunc->unuh);
    g_string_append_printf(buf, "#u(n;uF) = %g \n", slunc->unuF);

    g_string_append_printf(buf, "#Contributions to uncertainty of Sload\n");
    g_string_append_printf(buf, "#u(Sl;uh) = %g mN/nm \n", slunc->uSluh);
    g_string_append_printf(buf, "#u(Sl;uF) = %g mN/nm \n", slunc->uSluF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of m\n");
    g_string_append_printf(buf, "#u(m;uh) = %g \n", slunc->umuh);
    g_string_append_printf(buf, "#u(m;uF) = %g \n", slunc->umuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Sunload\n");
    g_string_append_printf(buf, "#u(Su;uh) = %g mN/nm \n", slunc->uSuuh);
    g_string_append_printf(buf, "#u(Su;uF) = %g mN/nm \n", slunc->uSuuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Ap(hc)\n");
    g_string_append_printf(buf, "#u(Ap(hc);uh) = %g nm^2 \n", slunc->uAuh);
    g_string_append_printf(buf, "#u(Ap(hc);uF) = %g nm^2 \n", slunc->uAuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Hit\n");
    g_string_append_printf(buf, "#u(Hit;uh) = %g MPa \n", slunc->uHituh);
    g_string_append_printf(buf, "#u(Hit;uF) = %g MPa \n", slunc->uHituF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Er\n");
    g_string_append_printf(buf, "#u(Er;uh) = %g GPa \n", slunc->uEruh);
    g_string_append_printf(buf, "#u(Er;uF) = %g GPa \n", slunc->uEruF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Eit\n");
    g_string_append_printf(buf, "#u(Eit;uh) = %g GPa \n", slunc->uEituh);
    g_string_append_printf(buf, "#u(Eit;uF) = %g GPa \n", slunc->uEituF);
    g_string_append_printf(buf, "#u(Eit;unu) = %g GPa \n", slunc->uEitunu);
    g_string_append_printf(buf, "#u(Eit;unui) = %g GPa \n", slunc->uEitunui);
    g_string_append_printf(buf, "#u(Eit;uEi) = %g GPa \n", slunc->uEituEi);
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Changes of the Er corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#Er/GPa \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (slunc->Ec[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, slunc->Ec[j + CONTACTMAX]);
        }
    }

    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Changes of the H corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#H/MPa \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (slunc->Hc[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, slunc->Hc[j + CONTACTMAX]);
        }
    }

    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Changes of the Ap corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#A/nm^2 \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (slunc->Ac[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, slunc->Ac[j + CONTACTMAX]);
        }
    }

    return g_string_free(buf, FALSE);
}

void slopes_fit_shift_contact(const Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata, gdouble *Ec, gdouble *Hc, gdouble *Ac, gint j)
{
    Slopesdata slopes;
    FDdata fdc;
    gdouble dh;
    gint istart, iend, ndata;
    GwyDataLine *x, *y;

    Ec[j + CONTACTMAX] = -1;
    Hc[j + CONTACTMAX] = -1;
    Ac[j + CONTACTMAX] = -1;

    if (fddata->i_contact_load + j < 0) {
        return;
    }

    // create shifted FDdata
    dh = create_shifted_FDdata(fddata, &fdc, j);

    // create working copy of Slopesdata
    init_shifted_Slopesdata(sldata, &slopes, dh);

    //check range
    if (sldata->Finputload) {
        if (slopes.loadfrom_pct_Fmax == slopes.loadto_pct_Fmax) {
            return;
        }
    }
    else {
        if (slopes.loadfrom == slopes.loadto) {
            return;
        }
    }

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (sldata->Finputload) {
        range_to_indices(slopes.loadfrom_pct_Fmax * 0.01 * fdc.Fmax, slopes.loadto_pct_Fmax * 0.01 * fdc.Fmax, fdc.Fload, FALSE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(slopes.loadfrom, slopes.loadto, fdc.hload, FALSE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, no fitting can be performed
    ndata = MIN(ndata, fdc.hload->res - istart);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    x = gwy_data_line_part_extract(fdc.hload, istart, ndata);
    y = gwy_data_line_part_extract(fdc.Fload, istart, ndata);

    // perform fit and get results
    slopes_fit_load(x, y, &slopes, &fdc, instdata, NULL);

    // clean up
    g_object_unref(x);
    g_object_unref(y);

    //check range
    if (sldata->Finputunload) {
        if (slopes.unloadfrom_pct_Fmax == slopes.unloadto_pct_Fmax) {
            return;
        }
    }
    else {
        if (slopes.unloadfrom == slopes.unloadto) {
            return;
        }
    }

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (sldata->Finputunload) {
        range_to_indices(slopes.unloadfrom_pct_Fmax * 0.01 * fdc.Fmax, slopes.unloadto_pct_Fmax * 0.01 * fdc.Fmax, fdc.Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(slopes.unloadfrom, slopes.unloadto, fdc.hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, no fitting can be performed
    ndata = MIN(ndata, fdc.hunload->res - istart);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    x = gwy_data_line_part_extract(fdc.hunload, istart, ndata);
    y = gwy_data_line_part_extract(fdc.Funload, istart, ndata);

    // perform fit and get results
    slopes_fit_unload(x, y, &slopes, &fdc, instdata, NULL);

    slopes_combine_fit_results(&slopes, &fdc, instdata);

    if (verbose) {
        g_print("j %d Sl %g Su %g Er %g H %g A %g \n", j, slopes.Sload, slopes.Sunload, slopes.Er, slopes.Hit, slopes.Aphc);
        g_print(" logfnm load %s %s \n", slopes.logfnmload, slopes.logfnmunload);
    }

    // clean up
    g_free(slopes.logfnmload);  //TODO
    g_free(slopes.logfnmunload);// TODO

    g_object_unref(x);
    g_object_unref(y);
    destroy_FDdata(&fdc);

    Ec[j + CONTACTMAX] = slopes.Er;
    Hc[j + CONTACTMAX] = slopes.Hit;
    Ac[j + CONTACTMAX] = slopes.Aphc;
}

void slopes_uncertainty_montecarlo_run_calc(const Slopesdata *slopesdata, SlopesUncdata *slopesunc, SlopesMCdata *slopesmc, FDdata *fddata)
{
    gint istart, iend, ndata;
    gint i, j, Nmc;
    Slopesdata *slopes;
    SlopesMCdata *mc;
    Instdata *instdata;
    GwyDataLine *hl, *Fl;
    GwyDataLine *hlnew, *Flnew;
    GwyDataLine *hlnoise, *Flnoise;
    GwyDataLine *hu, *Fu;
    GwyDataLine *hunew, *Funew;
    GwyDataLine *hunoise, *Funoise;
    gdouble hmax_old, Fmax_old;
    GRand *rng;
    gint nn;
    gdouble **dd;
    gdouble *xhist, *yhist;
    gboolean hist;
    gdouble estimateu[3], estimatel[3];

    instdata = &(slopesunc->instdata);


    //not necessary to check range, couldn't get here in the first place

    //find appropriate loading data, which should be fitted, use same regime as for the main fit
    if (slopesdata->Finputload)
        range_to_indices(slopesdata->loadfrom_pct_Fmax * 0.01 * fddata->Fmax,
                         slopesdata->loadto_pct_Fmax * 0.01 * fddata->Fmax,
                         fddata->Fload, FALSE, &istart, &iend, &ndata);
    else {
        range_to_indices(slopesdata->loadfrom, slopesdata->loadto, fddata->hload, FALSE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, not fitting can be performed
    ndata = MIN(ndata, fddata->hload->res - istart);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    hl = gwy_data_line_part_extract(fddata->hload, istart, ndata);
    Fl = gwy_data_line_part_extract(fddata->Fload, istart, ndata);

    //find appropriate unloadingdata, which should be fitted, use same regime as for the main fit
    if (slopesdata->Finputunload)
        range_to_indices(slopesdata->unloadfrom_pct_Fmax * 0.01 * fddata->Fmax,
                         slopesdata->unloadto_pct_Fmax * 0.01 * fddata->Fmax,
                         fddata->Funload, TRUE, &istart, &iend, &ndata);
    else {
        range_to_indices(slopesdata->unloadfrom, slopesdata->unloadto, fddata->hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, not fitting can be performed
    ndata = MIN(ndata, fddata->hunload->res - istart);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    hu = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
    Fu = gwy_data_line_part_extract(fddata->Funload, istart, ndata);

    //back up hmax, Fmax
    hmax_old = fddata->hmax;
    Fmax_old = fddata->Fmax;

    //use results of previous fit as input estimates
    estimatel[0] = slopesdata->gamma;
    estimatel[1] = slopesdata->h0;
    estimatel[2] = slopesdata->n;
    estimateu[0] = slopesdata->alpha;
    estimateu[1] = slopesdata->hp;
    estimateu[2] = slopesdata->m;

    Nmc = slopesmc->Nmc;
    slopes = slopesmc->slopesmc;

    //prepare data
    for (i = 0; i < Nmc ; i++) {
        init_Slopesdata(slopes);
        //change log filenames
        g_free(slopes->logfnmload);
        g_free(slopes->logfnmunload);
        // send to stdout, which can be redirected to /dev/null
        // otherwise too much IO time
        slopes->logfnmload = NULL;
        slopes->logfnmunload = NULL;
        //slopes->logfnmload = create_log_filename_path(mclogfnmload);
        //slopes->logfnmunload = create_log_filename_path(mclogfnmunload);
        slopes++;
    }

    slopes = slopesmc->slopesmc;

    //copy input parameters
    for (i = 0; i < Nmc ; i++) {
        /* opodr->nu = args->opodrdata.nu; */
        /* opodr->nui = args->opodrdata.nui; */
        /* opodr->Ei = args->opodrdata.Ei; */
        //		opodr->from = args->opodrdata.from;
        //		opodr->to = args->opodrdata.to;
        for (j = 0; j < 3; j++) {
            slopes->loadifixb[j] = slopesdata->loadifixb[j];
        }

        slopes++;
    }

    /*
       for (i=0; i< Nmc ;i++){
       printf ("nu %g \n", args->opodrmc.opodrmc[i].nu);
       }
       */

    // initialize random generator and datafield
    rng = g_rand_new();

    //prepare datafield for noise
    hlnoise = gwy_data_line_new_alike(hl, FALSE);
    Flnoise = gwy_data_line_new_alike(hl, FALSE);
    hunoise = gwy_data_line_new_alike(hu, FALSE);
    Funoise = gwy_data_line_new_alike(hu, FALSE);

    //prepare datafield for varied fields
    hlnew = gwy_data_line_new_alike(hl, FALSE);
    Flnew = gwy_data_line_new_alike(Fl, FALSE);
    hunew = gwy_data_line_new_alike(hu, FALSE);
    Funew = gwy_data_line_new_alike(Fu, FALSE);

    slopes = slopesmc->slopesmc;

    for (i = 0; i < Nmc; i++) {

        //create noise for h-F data
        generate_uncorrelated_noise(hlnoise, slopesunc->instdata.sigma_h, rng);
        generate_uncorrelated_noise(Flnoise, slopesunc->instdata.sigma_F, rng);
        generate_uncorrelated_noise(hunoise, slopesunc->instdata.sigma_h, rng);
        generate_uncorrelated_noise(Funoise, slopesunc->instdata.sigma_F, rng);

        //add noise to original h-F data
        if (!sum_data_lines(hlnew, hl, hlnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        if (!sum_data_lines(Flnew, Fl, Flnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        if (!sum_data_lines(hunew, hu, hunoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        if (!sum_data_lines(Funew, Fu, Funoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        //add random noise to hmax, Fmax
        fddata->hmax = hmax_old + gaussian_random_number(rng) * slopesunc->instdata.sigma_h;
        fddata->Fmax = Fmax_old + gaussian_random_number(rng) * slopesunc->instdata.sigma_F;

        //run fit

        slopes_fit_load(hlnew, Flnew, slopes, fddata, instdata, estimatel);
        slopes_fit_unload(hunew, Funew, slopes, fddata, instdata,  estimateu);
        slopes_combine_fit_results(slopes, fddata, instdata);

        slopes++;
    }

    //reset hmax, Fmax to original value
    fddata->hmax = hmax_old;
    fddata->Fmax = Fmax_old;

    //	printf("run Monte Carlo \n");
    g_object_unref(hl);
    g_object_unref(Fl);
    g_object_unref(hlnew);
    g_object_unref(Flnew);
    g_object_unref(hlnoise);
    g_object_unref(Flnoise);
    g_object_unref(hu);
    g_object_unref(Fu);
    g_object_unref(hunew);
    g_object_unref(Funew);
    g_object_unref(hunoise);
    g_object_unref(Funoise);
    g_rand_free(rng);

    for (i = 0; i < Nmc; i++) {
        //					printf(" %d Er %g Eit %g Hit %g \n", i, args->slopesmc.slopesmc[i].Er, args->slopesmc.slopesmc[i].Eit, args->slopesmc.slopesmc[i].Hit);
    }

    //calculate averages and stdevs
    //hp, Aphc, Hit, Er, Eit, S, m

    mc = slopesmc;

    for (i = 0; i < NPARAM_SLOPES; i++) {
        mc->mcavg[i] = 0;
        mc->mcstd[i] = 0;
    }

    slopes = slopesmc->slopesmc;
    nn = 0;

    for (i = 0; i < Nmc; i++) {
        //check sanity
        if (slopes->m > 1 && slopes->n > 1) {
            mc->mcavg[0] += slopes->Sload;
            mc->mcavg[1] += slopes->Sunload;
            mc->mcavg[2] += slopes->Hit;
            mc->mcavg[3] += slopes->Er;
            mc->mcavg[4] += slopes->Eit;
            mc->mcavg[5] += slopes->Aphc;
            mc->mcavg[6] += slopes->m;
            mc->mcavg[7] += slopes->n;
            mc->mcstd[0] += slopes->Sload * slopes->Sload;
            mc->mcstd[1] += slopes->Sunload * slopes->Sunload;
            mc->mcstd[2] += slopes->Hit * slopes->Hit;
            mc->mcstd[3] += slopes->Er * slopes->Er;
            mc->mcstd[4] += slopes->Eit * slopes->Eit;
            mc->mcstd[5] += slopes->Aphc * slopes->Aphc;
            mc->mcstd[6] += slopes->m * slopes->m;
            mc->mcstd[7] += slopes->n * slopes->n;
            nn++;
        }

        slopes++;
    }

    //use only reasonable data for average and standard deviation
    for (i = 0; i < NPARAM_SLOPES; i++) {
        mc->mcavg[i] /= nn;
        mc->mcstd[i] /= nn;
        mc->mcstd[i] -= mc->mcavg[i] * mc->mcavg[i];
        mc->mcstd[i] = (mc->mcstd[i] < 0) ? 0 : mc->mcstd[i];
        mc->mcstd[i] = sqrt(mc->mcstd[i]);
    }

    // how many datasets were skipped
    slopesmc->skipmc = Nmc - nn;

    /*
       printf("mean  stdev \n");
       for (i = 0; i< NPARAM_SLOPES; i++){
       printf(" %g   %g \n", avg[i], std[i]);
       }
       */

    // Histograms
    //create data structure

    dd = slopesmc->mcdd;
    slopes = slopesmc->slopesmc;

    for (i = 0; i < Nmc; i++) {
        dd[0][i] = slopes->Sload;
        dd[1][i] = slopes->Sunload;
        dd[2][i] = slopes->Hit;
        dd[3][i] = slopes->Er;
        dd[4][i] = slopes->Eit;
        dd[5][i] = slopes->Aphc;
        dd[6][i] = slopes->m;
        dd[7][i] = slopes->n;
        slopes++;
    }

    //prepare histogram structure

    for (i = 0; i < NPARAM_SLOPES; i++) {
        xhist = mc->mcxhist[i];
        yhist = mc->mcyhist[i];
        hist = create_histogram(dd[i], Nmc, xhist, yhist, mc->mcnstat);

        if (!hist) {
            g_printerr("Warning: No histogram could be created.\n");
            mc->mcnstat = -1;
            break;
        }
    }
}

gchar *slopes_mc_export_data(const Slopesdata *slopesdata, const SlopesUncdata *slopesunc, const SlopesMCdata *slopesmc,
                             const FDdata *fddata, enum ExportFormat exportformat)
{
    Slopesdata *sldata;
    gint i, j;
    gchar str[200];
    gboolean is;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Monte Carlo simulation of two slopes analysis (ODR fit)\n");
    g_string_append_printf(buf, "#Input uncertainties \n");
    g_string_append_printf(buf, "#u(h) = %g nm \n", slopesunc->instdata.sigma_h);
    g_string_append_printf(buf, "#u(F) = %g mN \n", slopesunc->instdata.sigma_F);
    g_string_append_printf(buf, "#Number of iterations  = %d  \n", slopesmc->Nmc);
    g_string_append_printf(buf, "\n");

    g_string_append_printf(buf, "#Parameters of two slopes analysis (ODR fit) \n");
    g_string_append_printf(buf, "#beta: %g \n", slopesdata->beta);
    g_string_append_printf(buf, "#nu: %g                 \n", slopesunc->instdata.nu);
    g_string_append_printf(buf, "#nu_indenter: %g        \n", slopesunc->instdata.nui);
    g_string_append_printf(buf, "#E_indenter: %g      GPa\n", slopesunc->instdata.Ei * 1e-9); //TODO units
    g_string_append_printf(buf, "\n");
    g_string_append_printf(buf, "#Fitted loadrange: ");

    if (slopesdata->Finputload) {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", slopesdata->loadrange_pct_Fmax[0], slopesdata->loadrange_pct_Fmax[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] nm \n", slopesdata->loadrange[0], slopesdata->loadrange[1]);
    }

    g_string_append_printf(buf, "# corresponds : ");

    if (slopesdata->Finputload) {
        g_string_append_printf(buf, "[%g , %g ] nm \n", slopesdata->loadrange[0], slopesdata->loadrange[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", slopesdata->loadrange_pct_Fmax[0], slopesdata->loadrange_pct_Fmax[1]);
    }

    g_string_append_printf(buf, "# %d datapoints used.\n", slopesdata->nfitloaddata);
    g_string_append_printf(buf, "\n");
    g_string_append_printf(buf, "#Fitted unloadrange: ");

    if (slopesdata->Finputunload) {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", slopesdata->unloadrange_pct_Fmax[0], slopesdata->unloadrange_pct_Fmax[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] nm \n", slopesdata->unloadrange[0], slopesdata->unloadrange[1]);
    }

    g_string_append_printf(buf, "# corresponds : ");

    if (slopesdata->Finputunload) {
        g_string_append_printf(buf, "[%g , %g ] nm \n", slopesdata->unloadrange[0], slopesdata->unloadrange[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", slopesdata->unloadrange_pct_Fmax[0], slopesdata->unloadrange_pct_Fmax[1]);
    }

    g_string_append_printf(buf, "# %d datapoints used.\n", slopesdata->nfitunloaddata);
    g_string_append_printf(buf, "\n");

    if (slopesmc->skipmc) {
        g_string_append_printf(buf, "#Fit failed in %d cases, the corresponding data are NOT included in the  estimates of mean and uncertainty but ARE included in the histograms.",
                               slopesmc->skipmc);
    }

    g_string_append_printf(buf, "\n");

    if (slopesmc->skipmc) {
        g_string_append_printf(buf, "\n\n");
    }

    g_string_append_printf(buf, "#Results \n");

    for (i = 0; i < NPARAM_SLOPES; i++) {
        is = nparam_to_string(i, str, sizeof(str));
        g_string_append_printf(buf, "#%s = (%g ± %g ) \n", str, slopesmc->mcavg[i], slopesmc->mcstd[i]);
    }

    g_string_append_printf(buf, "\n\n\n");


    // Histograms

    if (slopesmc->mcnstat > 0) {
        g_string_append_printf(buf, "#Histogram\n");

        for (i = 0; i < NPARAM_SLOPES; i++) {
            g_string_append_printf(buf, "#Histogram\n");
            is = nparam_to_string(i, str, sizeof(str));

            if (is) {
                g_string_append_printf(buf, "#%s    count/a.u.  \n", str);

                for (j = 0; j < slopesmc->mcnstat; j++) {
                    g_string_append_printf(buf, "%g %g \n", slopesmc->mcxhist[i][j], slopesmc->mcyhist[i][j]);
                }

                g_string_append_printf(buf, "\n \n");
            }
        }
    }

    g_string_append_printf(buf, "\n\n\n");

    //Full data
    g_string_append_printf(buf, "#Full data \n");
    g_string_append_printf(buf, "#");

    for (j = 0; j < NPARAM_SLOPES; j++) {
        is = nparam_to_string(j, str, sizeof(str)); /* _is_ is not used anywhere */
        g_string_append_printf(buf, "%s ", str);
    }

    g_string_append_printf(buf, "\n");

    sldata = slopesmc->slopesmc;

    for (i = 0; i < slopesmc->Nmc; i++) {
        g_string_append_printf(buf, SLOPE, sldata->Sload);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, SLOPE, sldata->Sunload);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, HARD, sldata->Hit);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, MODUL, sldata->Er);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, MODUL, sldata->Eit);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, AREA, sldata->Aphc);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, NUMBER, sldata->m);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, NUMBER, sldata->n);
        g_string_append_printf(buf, "  \n");
        sldata++;
    }

    return g_string_free(buf, FALSE);
}

void init_SlopesMCdata(SlopesMCdata *slopesmc, gint Nmc)
{
    gint i;

    //number of iterations
    slopesmc->Nmc = Nmc;

    if (verbose) {
        g_print("initialize for %d iterations \n", slopesmc->Nmc);
    }

    //how much had to be skipped
    slopesmc->skipmc = 0;

    //dynamically allocate
    slopesmc->slopesmc = (Slopesdata *)g_malloc(slopesmc->Nmc * sizeof(Slopesdata));
    slopesmc->mcavg = (gdouble *)g_malloc(NPARAM_SLOPES * sizeof(gdouble));
    slopesmc->mcstd = (gdouble *)g_malloc(NPARAM_SLOPES * sizeof(gdouble));
    slopesmc->mcdd = (gdouble **) g_malloc(NPARAM_SLOPES * sizeof(gdouble *));
    slopesmc->mcxhist = (gdouble **) g_malloc(NPARAM_SLOPES * sizeof(gdouble *));
    slopesmc->mcyhist = (gdouble **) g_malloc(NPARAM_SLOPES * sizeof(gdouble *));

    slopesmc->mcnstat = (gint)floor(3.49 * cbrt(slopesmc->Nmc) + 0.5);

    for (i = 0; i < NPARAM_SLOPES; i++) {
        slopesmc->mcdd[i] = (gdouble *) g_malloc(slopesmc->Nmc * sizeof(gdouble));
        slopesmc->mcxhist[i] = (gdouble *) g_malloc(slopesmc->mcnstat * sizeof(gdouble));
        slopesmc->mcyhist[i] = (gdouble *) g_malloc(slopesmc->mcnstat * sizeof(gdouble));
    }
}

void destroy_SlopesMCdata(SlopesMCdata *slopesmc)
{
    gint i;

    //free everything
    g_free(slopesmc->mcavg);
    g_free(slopesmc->mcstd);

    for (i = 0; i < NPARAM_SLOPES; i++) {
        g_free(slopesmc->mcdd[i]);
        g_free(slopesmc->mcxhist[i]);
        g_free(slopesmc->mcyhist[i]);
    }

    g_free(slopesmc->mcdd);
    g_free(slopesmc->mcxhist);
    g_free(slopesmc->mcyhist);

    slopesmc->mcnstat = 0;
    slopesmc->skipmc = 0;

    for (i = 0; i < slopesmc->Nmc; i++) {
        g_free(slopesmc->slopesmc[i].logfnmload);
        g_free(slopesmc->slopesmc[i].logfnmunload);
    }

    g_free(slopesmc->slopesmc);
}

static gboolean nparam_to_string(gint i, gchar *str, gint nstr)
{
    switch (i) {
    case 0:
        g_snprintf(str, nstr, "S_load/(mN/nm)");
        break;

    case 1:
        g_snprintf(str, nstr, "S_unload/(mN/nm)");
        break;

    case 2:
        g_snprintf(str, nstr, "H_IT/MPa");
        break;

    case 3:
        g_snprintf(str, nstr, "E_r/GPa");
        break;

    case 4:
        g_snprintf(str, nstr, "E_IT/GPa");
        break;

    case 5:
        g_snprintf(str, nstr, "A_p(h_c)/nm^2");
        break;

    case 6:
        g_snprintf(str, nstr, "m");
        break;

    case 7:
        g_snprintf(str, nstr, "n");
        break;

    default:
        g_printerr("Unhandled number of results.\n");
        return FALSE;
    }

    return TRUE;
}

gboolean slopes_unc_has_results(const SlopesUncdata *slopesunc)
{
    return (slopesunc->Ec != NULL && slopesunc->Hc != NULL);
}

#endif

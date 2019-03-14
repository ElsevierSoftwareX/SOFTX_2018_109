#ifndef NOFORTRAN

#include "tool_op_odr_unc.h"
#include "tool_op_odr.h"

#include "datatypes.h"
#include "fddata.h"
#include "file-utils.h"
#include "fit-utils.h"
#include "niget-common.h"
#include "mc-utils.h"

#include <math.h>
#include <string.h>

static const gchar *logfnm = "fit.log.op.unc";
//static const gchar *mclogfnm = "fit.log.op.mc";

static void init_shifted_OPODRdata(const OPODRdata *src, OPODRdata *dest, gdouble dh);
static gboolean nparam_to_string(gint i, gchar *str, gint nstr);

void init_OPODRUncdata(OPODRUncdata *opodrunc, Instdata *instdata)
{
    opodrunc->instdata = *instdata;

    /*
      opodrunc->unu = 0.05;
      opodrunc->unui = 0.04;
      opodrunc->uEi = 20e9;
    */
    opodrunc->Nmc = N_MC;

    opodrunc->uhp = 0 ;
    opodrunc->covmhp = 0 ;
    opodrunc->ueps = 0;

    //contributions from uh,uF
    opodrunc->umuh = 0 ;
    opodrunc->umuF = 0 ;
    opodrunc->uhcuh = 0;
    opodrunc->uhcuF = 0;
    opodrunc->uSuh = 0;
    opodrunc->uSuF = 0;
    opodrunc->uAuh = 0;
    opodrunc->uAuF = 0;
    opodrunc->uHituh = 0;
    opodrunc->uHituF = 0;
    opodrunc->uEruh = 0;
    opodrunc->uEruF = 0;
    opodrunc->uEituh = 0;
    opodrunc->uEituF = 0;

    //noise contributions
    opodrunc->uAnoise = 0;
    opodrunc->uHitnoise = 0;
    opodrunc->uErnoise = 0;
    opodrunc->uEitnoise = 0;

    //area coeff. contributions
    opodrunc->uAcoeff = 0;
    opodrunc->uHitcoeff = 0;
    opodrunc->uErcoeff = 0;
    opodrunc->uEitcoeff = 0;

    //radial contributions
    opodrunc->uAradial = 0;
    opodrunc->uHitradial = 0;
    opodrunc->uErradial = 0;
    opodrunc->uEitradial = 0;

    //sample-indenter contributions
    opodrunc->uEitunu = 0;
    opodrunc->uEitunui = 0;
    opodrunc->uEituEi = 0;

    //total contributions
    opodrunc->uhc = 0;
    opodrunc->um = 0 ;
    opodrunc->uS = 0;
    opodrunc->uA = 0;
    opodrunc->uHit = 0;
    opodrunc->uEr = 0;
    opodrunc->uEit = 0;

    // shift of contact
    opodrunc->Ec = NULL;
    opodrunc->Hc = NULL;
    opodrunc->Ac = NULL;
    opodrunc->hc = NULL;
}

static void init_shifted_OPODRdata(const OPODRdata *src, OPODRdata *dest, gdouble dh)
{
    init_OPODRdata(dest);

    // shift fitting range  in h-mode by dh
    // don't shift range for F-mode
    dest->range[0] = src->range[0] - dh;
    dest->range[1] = src->range[1] - dh;
    dest->from = dest->range[0];
    dest->to = dest->range[1];

    dest->from_pct_Fmax = src->from_pct_Fmax;
    dest->to_pct_Fmax = src->to_pct_Fmax;
    dest->range_pct_Fmax[0] = src->range_pct_Fmax[0];
    dest->range_pct_Fmax[1] = src->range_pct_Fmax[1];

    dest->Finput = src->Finput;

    dest->logfnm = create_log_filename_path(logfnm);

    // copy beta
    dest->beta = src->beta;

    dest->radial_corr = src->radial_corr;
    dest->radial_angle = src->radial_angle;
}

void op_odr_propagate_uncertainties_calc(const OPODRdata *opodr, OPODRUncdata *unc, const FDdata *fd, const Area *area)
{
    gdouble uh, uF;
    gdouble umh, umF, uhph, uhpF;
    gdouble uSh, uSF;
    gdouble  depsdm  ;
    gdouble uepsh, uepsF;
    gdouble epsh, epsF, mh, mF, depsdmh, depsdmF;
    gdouble um, uhp, ueps, uS;
    gdouble beta[3], std_beta[3], covariance[3];
    gint ifix[3];
    gint np;
    gint istart, iend, ndata;
    GwyDataLine *x, *y;
    gdouble ux, uy;
    gdouble dAdh;
    gdouble m, S, eps, hc, hp, hmax, Fmax, Aphc, Hit, Er, Eit, nu, nui, Ei;
    gdouble uEi, unu, unui;
    gdouble covmhph, covmhpF, covmhp;
    gdouble covShch, covShcF, covShc;
    gchar *unc_logfnm;
    gdouble SSres, R2, R2adj, chi2;
    gdouble uAp;
    gdouble alpha, acontact;

    //make copies to simplify code
    uh = unc->instdata.sigma_h;
    uF = unc->instdata.sigma_F;
    m = opodr->m;
    eps = opodr->eps;
    S = opodr->S;
    hp = opodr->hp;
    hc = opodr->hc;
    hmax = fd->hmax;
    Fmax = fd->Fmax;
    Aphc = opodr->Aphc;
    Hit = opodr->Hit;
    Eit = opodr->Eit;
    Er = opodr->Er;
    Ei = unc->instdata.Ei;
    nui = unc->instdata.nui;
    nu = unc->instdata.nu;
    uEi = unc->instdata.uEi;
    unui = unc->instdata.unui;
    unu = unc->instdata.unu;

    unc_logfnm =  create_log_filename_path(logfnm);

    if (verbose) {
        g_print("logfnm %s \n", unc_logfnm);
    }

    if (verbose) {
        g_print(" \n\nPropagate uncertainties\n");
    }

    ifix[0] = 1;
    ifix[1] = 1;
    ifix[2] = 1;
    np = 3;

    /*   u(m)    */

    /* prepare odr fit */

    beta[0] = opodr->alpha;
    beta[1] = opodr->hp;
    beta[2] = opodr->m;

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (opodr->Finput) {
        range_to_indices(opodr->range_pct_Fmax[0] * 0.01 * fd->Fmax, opodr->range_pct_Fmax[1] * 0.01 * fd->Fmax, fd->Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(opodr->range[0], opodr->range[1], fd->hunload, TRUE, &istart, &iend, &ndata);
    }

    x = gwy_data_line_part_extract(fd->hunload, istart, ndata);
    y = gwy_data_line_part_extract(fd->Funload, istart, ndata);

    /* fit with only uh */
    ux = uh;
    uy = 0;
    beta[0] = opodr->alpha;
    beta[1] = opodr->hp;
    beta[2] = opodr->m;

    odr_run_fit(x, y, np, beta, ifix, unc_logfnm, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);
    mh = beta[2];
    umh = std_beta[2];
    uhph = std_beta[1];
    covmhph = covariance[0];

    if (verbose)
        g_print(" fit with only uh:\n  alpha %g m %g hp %g \n  u(m) %g u(hp) %g corr(m,hp) %g \n",
                beta[0], beta[1], beta[2], umh, uhph, covmhph / umh / uhph);

    /*   u(eps)    */
    epsh = epsilon(mh);
    depsdmh = der_epsilon(mh);
    uepsh = fabs(depsdmh) * umh;

    if (verbose) {
        g_print("depsdm %g   \n", depsdmh);
        g_print(" u(eps;h) %g  urel(eps;h) %g \n", uepsh,  uepsh / eps);
    }

    /* fit with only uF */
    ux = 0;
    uy = uF;
    odr_run_fit(x, y, np, beta, ifix, unc_logfnm, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);
    mF = beta[2];
    umF = std_beta[2];
    uhpF = std_beta[1];
    covmhpF = covariance[0];

    if (verbose)
        g_print(" fit with only uF:\n  alpha %g m %g hp %g \n  u(m) %g u(hp) %g corr(m,hp) %g \n",
                beta[0], beta[1], beta[2], umF, uhpF, covmhpF / umF / uhpF);

    unc->umuh = umh;
    unc->umuF = umF;

    /*   u(eps)    */
    epsF = epsilon(mF);
    depsdmF = der_epsilon(mF);
    uepsF = fabs(depsdmF) * umF;

    if (verbose) {
        g_print("depsdm %g   \n", depsdmF);
        g_print(" u(eps;h) %g u(eps;F) %g  urel(eps;h) %g urel(eps;F) %g \n", uepsh, uepsF, uepsh / eps, uepsF / eps);
    }

    /*    u(S)   */
    uSh = sqrt(sq(umh) / sq(m) + sq(uh) / sq(hmax - hp) + sq(uhph) / sq(hmax - hp) + 2 * covmhph / m / (hmax - hp)) * S;
    uSF = sqrt(sq(umF) / sq(m) + sq(uF) / sq(Fmax) + sq(uhpF) / sq(hmax - hp) + 2 * covmhpF / m / (hmax - hp)) * S;
    unc->uSuh = uSh;
    unc->uSuF = uSF;
    /* u(S;h) = S *sqrt [ (umh/m)^2 + (uh/(hmax-hp))^2 + (uhph/(hmax-hp))^2  + 2 cov(m,hp;h)/(m*(hmax-hp))]
     * u(S;F) = S *sqrt [ (umF/m)^2 + (uFmax/Fmax)^2   + (uhpF/(hmax-hp))^2  + 2 cov(m,hp;F)/(m*(hmax-hp))]
     */

    if (verbose) {
        g_print(" u(S;h) %g u(S;F) %g  urel(S;h) %g urel(S;F) %g \n", uSh, uSF, uSh / S, uSF / S);
    }

    /*   u(hc)   */
    unc->uhcuh = sqrt(sq(hmax - hp) * sq(epsh / mh - depsdmh) / sq(mh) * sq(umh) + sq(1 - epsh / mh) * sq(uh) +
                      sq(epsh) / sq(mh) * sq(uhph) + 2 * (hmax - hp) * epsh / sq(mh) * (epsh / mh - depsdmh) * covmhph);

    unc->uhcuF = sqrt(sq(hmax - hp) * sq(epsF / mF - depsdmF) / sq(mF) * sq(umF) + sq(epsF) / sq(mF) * sq(uhpF) +
                      2 * (hmax - hp) * epsF / sq(mF) * (epsF / mF - depsdmF) * covmhpF);

    /* u(hc;h) ^2 = [(hmax-hp) (eps/m- deps/dm)/m  u(m;h) ]^2 + [(1-eps/m) u(hmax)]^2 + [eps/m u(hp;h)]^2 +
                    * 2 (hmax-hp)eps/m^2 (eps/m - deps/dm) cov(m,hp;h) */
    /* u(hc;F) ^2 = [(hmax-hp) (eps/m- deps/dm)/m  u(m;F) ]^2 +                         [eps/m u(hp;F)]^2 +
                    * 2 (hmax-hp)eps/m^2 (eps/m - deps/dm) cov(m,hp;F) */

    if (verbose) {
        g_print(" u(hc;h) %g u(hc;F) %g  urel(hc;h) %g urel(hc;F) %g \n", unc->uhcuh, unc->uhcuF, unc->uhcuh / hc, unc->uhcuF / hc);
    }

    /*    u(A)   */
    dAdh =  eval_area_derivative(area, hc);
    unc->uAuh = dAdh * unc->uhcuh;
    unc->uAuF = dAdh * unc->uhcuF;

    if (verbose) {
        g_print(" u(A;h) %g u(A;F) %g  urel(A;h) %g urel(A;F) %g \n", unc->uAuh, unc->uAuF, unc->uAuh / Aphc, unc->uAuF / Aphc);
    }

    /*     u(H_IT)   */
    unc->uHituh = Hit / Aphc * dAdh * unc->uhcuh;
    unc->uHituF = Hit * sqrt(sq(uF) / sq(Fmax) + sq(dAdh) / sq(Aphc) * sq(unc->uhcuh));

    if (verbose) {
        g_print(" u(Hit;h) %g u(Hit;F) %g  urel(Hit;h) %g urel(Hit;F) %g \n", unc->uHituh, unc->uHituF, unc->uHituh / Hit, unc->uHituF / Hit);
    }

    /*     u(E_r)   */

    covShch = - S * (1 - eps) / (hmax - hp) * sq(uh) - S * (hmax - hp) / m * depsdmh * sq(umh) +
              S * eps / (hmax - hp) * sq(uhph) + S * (eps / m - depsdmh) * covmhph;
    covShcF =  - S * (hmax - hp) / m * depsdmF * sq(umF) + S * eps / (hmax - hp) * sq(uhpF) + S * (eps / m - depsdmF) * covmhpF;
    /*cov(S,hc;h) = -S(1-eps)/(hmax-hp)*uh^2 -S(hmax-hp)/m *(deps/dm) *u(m;h)^2 + S*eps/(hmax-hp)*u(hp;h)^2 +S (eps/m-(deps/dm))cov(m,hp;h)
     *cov(S,hc;F) =                          -S(hmax-hp)/m *(deps/dm) *u(m;F)^2 + S*eps/(hmax-hp)*u(hp;F)^2 +S (eps/m-(deps/dm))cov(m,hp;F)
     */

    unc->uEruh = Er * sqrt(sq(uSh) / sq(S) + sq(unc->uAuh) / sq(Aphc) / 4 - 1 / S / Aphc * dAdh * covShch);
    unc->uEruF = Er * sqrt(sq(uSF) / sq(S) + sq(unc->uAuF) / sq(Aphc) / 4 - 1 / S / Aphc * dAdh * covShcF);
    /* u(Er;h) = sqrt[ (Er/S)^2 u(S;h)^2 + (Er/2A)^2 u(A;h)^2  - Er^2/(S Aphc)* (dA/dh) cov(S,hc;h)]
     * u(Er;F) = sqrt[ (Er/S)^2 u(S;F)^2 + (Er/2A)^2 u(A;F)^2  - Er^2/(S Aphc)* (dA/dh) cov(S,hc;F)]
     */

    if (verbose) {
        g_print(" u(Er;h) %g u(Er;F) %g  urel(Er;h) %g urel(Er;F) %g \n", unc->uEruh, unc->uEruF, unc->uEruh / Er, unc->uEruF / Er);
    }

    /*     u(E_IT)   */
    unc->uEituh = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruh; // Eit^2/Er^2 /(1-nu^2) u(Er;h)
    unc->uEituF = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruF; // Eit^2/Er^2 /(1-nu^2) u(Er;F)

    if (verbose) {
        g_print("u(Eit;h) %g  u(Eit;F) %g  urel(Eit;h) %g  urel(Eit;F) %g\n", unc->uEituh, unc->uEituF, unc->uEituh / Eit, unc->uEituF / Eit);
    }

    // total uncertainties
    // should take um as given by odr_wrap, otherwise doublecounting
    /* fit with both uh and uF */
    ux = uh;
    uy = uF;
    odr_run_fit(x, y, np, beta, ifix, unc_logfnm, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);

    if (verbose) {
        g_print(" fit with both ux %g and uy %g \n", ux, uy);
    }

    um = std_beta[2];
    uhp = std_beta[1];
    covmhp = covariance[0];

    unc->um = um;
    unc->uhp = uhp;
    unc->covmhp = covmhp;

    if (verbose)
        g_print(" fit with both uh and uF:\n  alpha %g hp %g m %g \n  u(m) %g u(hp) %g corr(m,hp) %g \n",
                beta[0], beta[1], beta[2], um, uhp, covmhp / um / uhp);

    /*   u(eps)    */
    depsdm = der_epsilon(m);
    ueps = fabs(depsdm) * um;
    unc->ueps = ueps;

    if (verbose) {
        g_print("eps %g    \n", eps);
        g_print("depsdm %g  \n", depsdm);
        g_print(" u(eps) %g    urel(eps) %g  \n", ueps, ueps / eps);

        gdouble e1, e2, e3;

        e2 = m * (1 - 2 * (m - 1) * tgamma(m / 2 / (m - 1)) / tgamma(1.0 / 2.0 / (m - 1)) / sqrt(M_PI));
        e1 = (m - 0.1) * (1 - 2 * ((m - 0.1) - 1) * tgamma((m - 0.1) / 2 / ((m - 0.1) - 1)) / tgamma(1.0 / 2.0 / ((m - 0.1) - 1)) / sqrt(M_PI));
        e3 = (m + 0.1) * (1 - 2 * ((m + 0.1) - 1) * tgamma((m + 0.1) / 2 / ((m + 0.1) - 1)) / tgamma(1.0 / 2.0 / ((m + 0.1) - 1)) / sqrt(M_PI));
        g_print(" eps(m-dm) %g eps (m) %g eps (m+dm) %g \n", e1, e2, e3);
        g_print(" left deps %g right depts %g symm deps %g \n", (e2 - e1) / 0.1, (e3 - e2) / 0.1, (e3 - e1) / 0.2);
    }

    /*    u(S)   */
    uS = sqrt(sq(um) / sq(m) + sq(uh) / sq(hmax - hp) +  sq(uF) / sq(Fmax) + sq(uhp) / sq(hmax - hp) + 2 * covmhp / m / (hmax - hp)) * S;
    unc->uS = uS;
    /* u(S) = S *sqrt [ (umh/m)^2 +  (uhmax/(hmax-hp))^2 + (uFmax/Fmax)^2 +  (uhp/(hmax-hp))^2 + 2 cov(m,hp)/(m*(hmax-hp)) ] */

    if (verbose) {
        g_print(" u(S) %g    urel(S) %g \n", uS, uS / S);
    }

    /*   u(hc)   */
    unc->uhc = sqrt(sq(hmax - hp) * sq(eps / m - depsdm) / sq(m) * sq(um) + sq(1 - eps / m) * sq(uh) +
                    sq(eps) / sq(m) * sq(uhp) + 2 * (hmax - hp) * eps / sq(m) * (eps / m - depsdm) * covmhp);
    /* u(hc) ^2 = [(hmax- hp) (eps/m- deps/dm)/m  u(m) ]^2 + [(1-eps/m) u(hmax)]^2 + [eps/m u(hp)]^2 + 2 (hmax-hp)eps/m^2 (eps/m - deps/dm) cov(m,hp) */

    if (verbose) {
        g_print(" u(hc) %g    urel(hc) %g \n", unc->uhc, unc->uhc / hc);
        /*
          g_print(" %g %g %g %g sum %g \n", (hmax-hc)*(hmax-hc)* (eps/m - depsdm)*(eps/m-depsdm)/m/m *um*um, (1-eps/m)*(1-eps/m)*uh*uh,
          eps*eps/m/m*uhp*uhp, 2*(hmax-hp)*eps/m/m*(eps/m-depsdm)*covmhp,
          ((hmax-hc)*(hmax-hc)* (eps/m - depsdm)*(eps/m-depsdm)/m/m *um*um + (1-eps/m)*(1-eps/m)*uh*uh + eps*eps/m/m*uhp*uhp + 2*(hmax-hp)*eps/m/m*(eps/m-depsdm)*covmhp)
          );
        */
    }

    /*    u(A)   */
    dAdh =  eval_area_derivative(area, hc) ;
    unc->uAnoise = dAdh * unc->uhc;
    uAp = eval_area_uncertainty(area, hc);

    if (uAp < 0) {
        unc->uA = unc->uAnoise;
        unc->uAcoeff = -1;
    }
    else {
        unc->uAcoeff = sqrt(eval_area_uncertainty(area, hc));
        unc->uA = sqrt(sq(unc->uAnoise) + sq(unc->uAcoeff));
    }

    if (verbose) {
        g_print(" u(A; noise) %g  u(A;coeff) %g \n", unc->uAnoise, unc->uAcoeff);
        g_print("relative u(A; noise) %g  u(A;coeff) %g \n", unc->uAnoise / opodr->Aphc, unc->uAcoeff / opodr->Aphc);
        g_print(" u(A) %g    urel(A) %g  \n", unc->uA, unc->uA / Aphc);
    }

    /*     u(H_IT)   */
    unc->uHit = Hit * sqrt(sq(uF) / sq(Fmax) + sq(unc->uA) / sq(opodr->Aphc));
    //contributions from noise and area coefficients
    unc->uHitnoise = Hit * sqrt(sq(uF) / sq(Fmax) + sq(unc->uAnoise) / sq(opodr->Aphc));
    unc->uHitcoeff = Hit * unc->uAcoeff / opodr->Aphc;

    if (verbose) {
        g_print(" u(Hit;h) %g u(Hit;F) %g  urel(Hit;h) %g urel(Hit;F) %g \n", unc->uHituh, unc->uHituF, unc->uHituh / Hit, unc->uHituF / Hit);
        g_print(" u(Hit;noise) %g u(Hit;coeff) %g urel(Hit;noise) %g urel(Hit;coeff) %g \n", unc->uHitnoise, unc->uHitcoeff, unc->uHitnoise / Hit, unc->uHitcoeff / Hit);
    }

    /*     u(E_r)   */
    covShc = - S * (1 - eps) / (hmax - hp) * sq(uh) - S * (hmax - hp) / m * depsdm * sq(um) +
             S * eps / (hmax - hp) * sq(uhp) + S * (eps / m - depsdm) * covmhp;
    /*cov(S,hc;h) = -S(1-eps)/(hmax-hp)*uh^2 -S(hmax-hp)/m *(deps/dm) *u(m)^2 + S*eps/(hmax-hp)*u(hp)^2 +S (eps/m-(deps/dm))cov(m,hp) */

    unc->uEr = Er * sqrt(sq(uS) / sq(S) + sq(unc->uA) / sq(Aphc) / 4 - 1 / S / Aphc * dAdh * covShc);
    /* u(Er) = sqrt[ (Er/S)^2 u(S)^2 + (Er/2A)^2 u(A)^2  - Er^2/(S Aphc)* (dA/dh) cov(S,hc)] */

    //contributions from noise and area coefficients
    unc->uErnoise = Er * sqrt(sq(uS) / sq(S) + sq(unc->uAnoise) / sq(Aphc) / 4 - 1 / S / Aphc * dAdh * covShc);
    unc->uErcoeff = Er * unc->uAcoeff / Aphc / 2;

    if (verbose) {
        g_print(" u(Er) %g   urel(Er) %g  \n", unc->uEr, unc->uEr / Er);
        g_print(" u(Er;noise) %g   u(Er;coeff) %g  urel(Er;noise) %g   urel(Er;coeff) %g   \n", unc->uErnoise, unc->uErcoeff, unc->uErnoise / Er, unc->uErcoeff / Er);
    }

    /*     u(E_IT)   */
    //contribution from h,F
    unc->uEit = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEr; // Eit^2/Er^2 /(1-nu^2) uEr
    //contributions from noise and area coefficients
    unc->uEitnoise = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uErnoise; // Eit^2/Er^2 /(1-nu^2) uEr
    unc->uEitcoeff = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uErcoeff; // Eit^2/Er^2 /(1-nu^2) uEr

    //contribution from nu, nui, Ei
    unc->uEitunu =  fabs(2 * nu / (1 - sq(nu)) * Eit * unu); // -2nu/(1-nu^2) Eit u(nu)
    unc->uEitunui = fabs(2 * nui / (1 - sq(nu)) * sq(Eit) / Ei * unui) * 1e9; // -2nui/(1-nu^2) Eit^2/Ei u(nui)
    unc->uEituEi =  fabs((1 - sq(nui)) / (1 - sq(nu)) * sq(Eit) / sq(Ei) * uEi) * 1e9; //-(1-nui^2)/(1-nu)^2 Eit^2/Ei^2 u(Ei)

    if (verbose) {
        g_print("u(Eit;h,F) %g    urel(Eit;hF) %g  \n", unc->uEit, unc->uEit / Eit);
        g_print(" u(Eit) %g   urel(Eit) %g  \n", unc->uEit, unc->uEit / Eit);
        g_print(" u(Eit;noise) %g   u(Eit;coeff) %g  urel(Eit;noise) %g   urel(Eit;coeff) %g   \n", unc->uEitnoise, unc->uEitcoeff, unc->uEitnoise / Eit, unc->uEitcoeff / Eit);
        g_print("u(Eit;nu) %g  u(Eit;nui) %g  u(Eit;Ei) %g  urel(Eit;nu) %g  urel(Eit;nui) %g  urel(Eit;Ei) %g\n",
                unc->uEitunu, unc->uEitunui, unc->uEituEi, unc->uEitunu / Eit, unc->uEitunui / Eit, unc->uEituEi / Eit);
    }

    unc->uEit = sqrt(sq(unc->uEit) + sq(unc->uEitunu) + sq(unc->uEitunui) + sq(unc->uEituEi));

    if (verbose) {
        g_print("u(Eit;hFnuetc) %g    urel(Eit; hF_nu_etc) %g \n", unc->uEit, unc->uEit / Eit);
    }

    if (opodr->radial_corr) {
        if (verbose) {
            g_print("Start calculation of uncertainty from radial correction angle.\n");
        }

        acontact = sqrt(Fmax / (M_PI * Hit) * 1e-9) * 1e9;
        alpha = M_PI / 2 - atan(hp / acontact);
        alpha *= 180. / M_PI;

        if (verbose) {
            g_print("Lower limit to radial correction alpha = %g degree \n", alpha / M_PI * 180.);
        }

        //The choice of interval is the same and is done at the beginning of this function, the fitting was also already done, we only have to evaluate H,Er,etc. and apply the radial correction

        gdouble Aphc0, Hit0, Er0, Eit0;
        gdouble Aphc2, Hit2, Er2, Eit2;

        Aphc0 = eval_area(area, hc);

        Hit0 = fd->Fmax / Aphc0 * 1e9;
        Er0 = sqrt(M_PI / Aphc0) * S / (2 * opodr->beta) * 1e15;
        Eit0 = (1 - sq(nu)) / (1 / Er0 - (1 - sq(nui)) / Ei);
        Eit0 *= 1e-9;
        Er0 *= 1e-9;

        Aphc2 = Aphc0;
        Hit2 = Hit0;
        Eit2 = Eit0;
        Er2 = Er0;

        radial_correction_op(alpha, &Hit2, &Er2, &Eit2, &Aphc2, &unc->instdata);

        // Aphc,Hit, Eit, Er are with facet angle, Aphc1,Hi1, Er1,Eit1 with arctan(hp/a) angle
        unc->uAradial = (Aphc - Aphc2) / 2 / sqrt(3.);
        unc->uHitradial = (Hit2 - Hit) / 2 / sqrt(3.);
        unc->uEitradial = (Eit2 - Eit) / 2 / sqrt(3.);
        unc->uErradial = (Er2 - Er) / 2 / sqrt(3.);

        unc->uA = sqrt(sq(unc->uA) + sq(unc->uAradial));
        unc->uHit = sqrt(sq(unc->uHit) + sq(unc->uHitradial));
        unc->uEit = sqrt(sq(unc->uEit) + sq(unc->uEitradial));
        unc->uEr = sqrt(sq(unc->uEr) + sq(unc->uErradial));

        if (verbose) {
            g_print("uA radial %g \n", unc->uAradial);
            g_print("uH radial %g \n", unc->uHitradial);
            g_print("uEr radial %g \n", unc->uErradial);
            g_print("uEit radial %g \n", unc->uEitradial);
        }
    }

    if (verbose) {
        g_print(" \nEnd propagation of uncertainties\n\n");
    }

    g_object_unref(x);
    g_object_unref(y);
    g_free(unc_logfnm);
}

gchar *op_odr_uncertainty_export_data(const OPODRdata *opodrdata, const OPODRUncdata *opodrunc,
                                      const FDdata *fddata, const Area *area, enum ExportFormat exportformat)
{
    gint i, j;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Parameters of Oliver Pharr analysis (ODR fit) \n");
    g_string_append_printf(buf, "# u(depth): %g nm\n", opodrunc->instdata.sigma_h);
    g_string_append_printf(buf, "# u(load): %g nm\n", opodrunc->instdata.sigma_F);
    g_string_append_printf(buf, "#beta: %g \n", opodrdata->beta);
    g_string_append_printf(buf, "#nu: %g                 u(nu):  %g  \n", opodrunc->instdata.nu, opodrunc->instdata.unu);
    g_string_append_printf(buf, "#nu_indenter: %g       u(nu_indenter):  %g   \n", opodrunc->instdata.nui, opodrunc->instdata.unui);
    g_string_append_printf(buf, "#E_indenter: %g        u(E_indenter):   %g GPa\n", opodrunc->instdata.Ei * 1e-9, opodrunc->instdata.uEi * 1e-9); //TODO units
    g_string_append_printf(buf, "#Fitted range: ");

    if (opodrdata->Finput) {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", opodrdata->range_pct_Fmax[0], opodrdata->range_pct_Fmax[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] nm \n", opodrdata->range[0], opodrdata->range[1]);
    }

    g_string_append_printf(buf, "# corresponds : ");

    if (opodrdata->Finput) {
        g_string_append_printf(buf, "[%g , %g ] nm \n", opodrdata->range[0], opodrdata->range[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", opodrdata->range_pct_Fmax[0], opodrdata->range_pct_Fmax[1]);
    }

    g_string_append_printf(buf, "# %d datapoints used.\n", opodrdata->nfitdata);
    g_string_append_printf(buf, "#\n");

    if (opodrdata->radial_corr) {
        g_string_append_printf(buf, "#Applied radial correction using α: %.4g°\n", opodrdata->radial_angle);
    }
    else {
        g_string_append_printf(buf, "#Radial correction not applied\n");
    }

    if (opodrunc->uAcoeff != 0) {
        g_string_append_printf(buf, "#Area function polynomial coefficients\n");
        g_string_append_printf(buf, "#p = (");

        for (i = 0; i < area->npoly; i++) {
            g_string_append_printf(buf, " %g ", area->polycoef[i]);
        }

        g_string_append_printf(buf, ")\n");
        g_string_append_printf(buf, "#Covariance matrix of area function polynomial coefficients\n");

        for (i = 0; i < area->npoly; i++) {
            g_string_append_printf(buf, "#");

            for (j = 0; j < area->npoly; j++) {
                g_string_append_printf(buf, " %g ", area->covpolycoef[i][j]);
            }

            g_string_append_printf(buf, "\n");
        }

        g_string_append_printf(buf, "\n");
        g_string_append_printf(buf, "\n");
        g_string_append_printf(buf, "\n");
    }

    g_string_append_printf(buf, "#Results of Oliver Pharr analysis (ODR fit) \n");
    g_string_append_printf(buf, "#hc: %g nm\n", opodrdata->hc);
    g_string_append_printf(buf, "#Ap(hc): %g nm^2\n", opodrdata->Aphc);
    g_string_append_printf(buf, "#Er: %g GPa\n", opodrdata->Er);
    g_string_append_printf(buf, "#Hit: %g MPa\n", opodrdata->Hit);
    g_string_append_printf(buf, "#Eit: %g GPa\n", opodrdata->Eit);
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Results of the Gaussian propagation of uncertainties \n");
    g_string_append_printf(buf, "#u(m) =  %g  \n", opodrunc->um);
    g_string_append_printf(buf, "#u(hp) =  %g nm \n", opodrunc->uhp);
    g_string_append_printf(buf, "#cov(m,hp) =  %g nm \n", opodrunc->covmhp);
    g_string_append_printf(buf, "#u(eps) = %g  \n", opodrunc->ueps);
    g_string_append_printf(buf, "#u(S) =  %g mN/nm \n", opodrunc->uS);
    g_string_append_printf(buf, "#u(hc) =  %g nm \n", opodrunc->uhc);

    g_string_append_printf(buf, "#u(Ap(hc)) =  %g nm^2 \n", opodrunc->uA);
    g_string_append_printf(buf, "#u(Hit) =  %g MPa \n", opodrunc->uHit);
    g_string_append_printf(buf, "#u(Er) =  %g GPa \n", opodrunc->uEr);
    g_string_append_printf(buf, "#u(Eit) =  %g GPa \n", opodrunc->uEit);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of m\n");
    g_string_append_printf(buf, "#u(m;uh) = %g \n", opodrunc->umuh);
    g_string_append_printf(buf, "#u(m;uF) = %g \n", opodrunc->umuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of S\n");
    g_string_append_printf(buf, "#u(S;uh) = %g mN/nm \n", opodrunc->uSuh);
    g_string_append_printf(buf, "#u(S;uF) = %g mN/nm \n", opodrunc->uSuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of hc\n");
    g_string_append_printf(buf, "#u(hc;uh) = %g nm \n", opodrunc->uhcuh);
    g_string_append_printf(buf, "#u(hc;uF) = %g nm \n", opodrunc->uhcuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Ap(hc)\n");
    g_string_append_printf(buf, "#u(Ap(hc);uh) = %g nm^2 \n", opodrunc->uAuh);
    g_string_append_printf(buf, "#u(Ap(hc);uF) = %g nm^2 \n", opodrunc->uAuF);
    g_string_append_printf(buf, "#u(Ap(hc);noise) =  %g nm^2 \n", opodrunc->uAnoise);

    if (opodrunc->uAcoeff == -1) {
        g_string_append_printf(buf, "#covariance of polynomial coefficients is not positive definite\n");
    }
    else {
        g_string_append_printf(buf, "#u(Ap(hc);coeff) =  %g nm^2 \n", opodrunc->uAcoeff);
    }

    g_string_append_printf(buf, "#u(Ap(hc);radial) =  %g nm^2 \n", opodrunc->uAradial);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Hit\n");
    g_string_append_printf(buf, "#u(Hit;uh) = %g MPa \n", opodrunc->uHituh);
    g_string_append_printf(buf, "#u(Hit;uF) = %g MPa \n", opodrunc->uHituF);
    g_string_append_printf(buf, "#u(Hit;noise) =  %g MPa \n", opodrunc->uHitnoise);
    g_string_append_printf(buf, "#u(Hit;coeff) =  %g MPa \n", opodrunc->uHitcoeff);
    g_string_append_printf(buf, "#u(Hit;radial) =  %g MPa \n", opodrunc->uHitradial);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Er\n");
    g_string_append_printf(buf, "#u(Er;uh) = %g GPa \n", opodrunc->uEruh);
    g_string_append_printf(buf, "#u(Er;uF) = %g GPa \n", opodrunc->uEruF);
    g_string_append_printf(buf, "#u(Er;noise) =  %g GPa \n", opodrunc->uErnoise);
    g_string_append_printf(buf, "#u(Er;coeff) =  %g GPa \n", opodrunc->uErcoeff);
    g_string_append_printf(buf, "#u(Er;radial) =  %g GPa \n", opodrunc->uErradial);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Eit\n");
    g_string_append_printf(buf, "#u(Eit;uh) = %g GPa \n", opodrunc->uEituh);
    g_string_append_printf(buf, "#u(Eit;uF) = %g GPa \n", opodrunc->uEituF);
    g_string_append_printf(buf, "#u(Eit;noise) =  %g GPa \n", opodrunc->uEitnoise);
    g_string_append_printf(buf, "#u(Eit;coeff) =  %g GPa \n", opodrunc->uEitcoeff);
    g_string_append_printf(buf, "#u(Eit;radial) =  %g GPa \n", opodrunc->uEitradial);
    g_string_append_printf(buf, "#u(Eit;unu) = %g GPa \n", opodrunc->uEitunu);
    g_string_append_printf(buf, "#u(Eit;unui) = %g GPa \n", opodrunc->uEitunui);
    g_string_append_printf(buf, "#u(Eit;uEi) = %g GPa \n", opodrunc->uEituEi);
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Changes of the Er corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#Er/GPa \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (opodrunc->Ec[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, opodrunc->Ec[j + CONTACTMAX]);
        }
    }

    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#Changes of the H corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#H/MPa \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (opodrunc->Hc[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, opodrunc->Hc[j + CONTACTMAX]);
        }
    }

    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#Changes of the area corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#Aphc/nm^2 \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (opodrunc->Ac[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, opodrunc->Ac[j + CONTACTMAX]);
        }
    }

    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#Changes of the contact depth corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#hc/nm \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (opodrunc->hc[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, opodrunc->hc[j + CONTACTMAX]);
        }
    }

    return g_string_free(buf, FALSE);
}

void op_odr_fit_shift_contact(const OPODRdata *opodrdata, const FDdata *fddata, const Instdata *instdata, const Area *area,
                              gdouble *Ec, gdouble *Hc, gdouble *Ac, gdouble *hc, gint j)
{
    OPODRdata opodr;
    FDdata fdc;
    gdouble dh;
    gint istart, iend, ndata;
    GwyDataLine *x, *y;

    Ec[j + CONTACTMAX] = -1;
    Hc[j + CONTACTMAX] = -1;
    Ac[j + CONTACTMAX] = -1;
    hc[j + CONTACTMAX] = -1;

    if (fddata->i_contact_load + j < 0) {
        return;
    }

    // create shifted FDdata
    dh = create_shifted_FDdata(fddata, &fdc, j);

    // create working copy of OPODRdata
    init_shifted_OPODRdata(opodrdata, &opodr, dh);

    if (verbose) {
        g_print("logfnm %s \n", opodr.logfnm);
    }

    //check range

    if (opodrdata->Finput) {
        if (opodr.from_pct_Fmax == opodr.to_pct_Fmax) {
            return;
        }
    }
    else {
        if (opodr.from == opodr.to) {
            return;
        }
    }

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (opodrdata->Finput) {
        range_to_indices(opodr.from_pct_Fmax * 0.01 * fdc.Fmax, opodr.to_pct_Fmax * 0.01 * fdc.Fmax, fdc.Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(opodr.from, opodr.to, fdc.hunload, TRUE, &istart, &iend, &ndata);
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
    op_odr_fit(x, y, &opodr, &fdc, instdata, area, NULL);

    if (opodr.radial_corr) {
        radial_correction_op(opodr.radial_angle, &(opodr.Hit), &(opodr.Er), &(opodr.Eit), &(opodr.Aphc), instdata);
    }

    if (verbose) {
        g_print("j %d Su %g Er %g H %g A %g \n", j, opodr.S, opodr.Er, opodr.Hit, opodr.Aphc);
    }

    // clean up
    g_free(opodr.logfnm);
    g_object_unref(x);
    g_object_unref(y);
    destroy_FDdata(&fdc);

    Ec[j + CONTACTMAX] = opodr.Er;
    Hc[j + CONTACTMAX] = opodr.Hit;
    Ac[j + CONTACTMAX] = opodr.Aphc;
    hc[j + CONTACTMAX] = opodr.hc;
}

#define NPARAM 9  //hc, Aphc, Hit, Er, Eit, S, m, hp, eps

void init_OPODRMCdata(OPODRMCdata *opodrmc, gint Nmc)
{
    gint i;

    //number of iterations
    opodrmc->Nmc = Nmc;

    if (verbose) {
        g_print("initialize for %d iterations \n", opodrmc->Nmc);
    }

    //how much had to be skipped
    opodrmc->skipmc = 0;

    //dynamically allocate
    opodrmc->opodrmc = (OPODRdata *)g_malloc(opodrmc->Nmc * sizeof(OPODRdata));
    opodrmc->mcavg = (gdouble *)g_malloc(NPARAM * sizeof(gdouble));
    opodrmc->mcstd = (gdouble *)g_malloc(NPARAM * sizeof(gdouble));
    opodrmc->mcdd = (gdouble **) g_malloc(NPARAM * sizeof(gdouble *));
    opodrmc->mcxhist = (gdouble **) g_malloc(NPARAM * sizeof(gdouble *));
    opodrmc->mcyhist = (gdouble **) g_malloc(NPARAM * sizeof(gdouble *));

    opodrmc->mcnstat = (gint)floor(3.49 * cbrt(opodrmc->Nmc) + 0.5);

    for (i = 0; i < NPARAM; i++) {
        opodrmc->mcdd[i] = (gdouble *) g_malloc(opodrmc->Nmc * sizeof(gdouble));
        opodrmc->mcxhist[i] = (gdouble *) g_malloc(opodrmc->mcnstat * sizeof(gdouble));
        opodrmc->mcyhist[i] = (gdouble *) g_malloc(opodrmc->mcnstat * sizeof(gdouble));
    }
}

void destroy_OPODRMCdata(OPODRMCdata *opodrmc)
{
    gint i;

    //free everything
    g_free(opodrmc->mcavg);
    g_free(opodrmc->mcstd);

    for (i = 0; i < NPARAM; i++) {
        g_free(opodrmc->mcdd[i]);
        g_free(opodrmc->mcxhist[i]);
        g_free(opodrmc->mcyhist[i]);
    }

    g_free(opodrmc->mcdd);
    g_free(opodrmc->mcxhist);
    g_free(opodrmc->mcyhist);

    opodrmc->mcnstat = 0;
    opodrmc->skipmc = 0;

    for (i = 0; i < opodrmc->Nmc; i++) {
        g_free(opodrmc->opodrmc[i].logfnm);
    }

    g_free(opodrmc->opodrmc);
}

/* tady kdyz uz, tak deklarovat nejake konstantni pole a volat pak primo cislem z toho */
static gboolean nparam_to_string(gint i, gchar *str, gint nstr)
{
    switch (i) {
    case 0:
        g_snprintf(str, nstr, "h_c/nm");
        break;

    case 1:
        g_snprintf(str, nstr, "A_p(h_c)/nm^2");
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
        g_snprintf(str, nstr, "S/(mN/nm)");
        break;

    case 6:
        g_snprintf(str, nstr, "m");
        break;

    case 7:
        g_snprintf(str, nstr, "h_p/nm");
        break;

    case 8:
        g_snprintf(str, nstr, "eps");
        break;

    default:
        g_printerr("Unhandled number of results.\n");
        return FALSE;
    }

    return TRUE;
}

gchar *op_odr_mc_export_data(const OPODRdata *opodrdata, const OPODRUncdata *opodrunc, const OPODRMCdata *opodrmc,
                             const FDdata *fddata, enum ExportFormat exportformat)
{
    gint i, j;
    OPODRdata *opodr;
    gchar str[200];
    gboolean is;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Monte Carlo simulation of Oliver Pharr analysis (ODR fit)\n");
    g_string_append_printf(buf, "#Input uncertainties \n");
    g_string_append_printf(buf, "#u(h) = %g nm \n", opodrunc->instdata.sigma_h);
    g_string_append_printf(buf, "#u(F) = %g mN \n", opodrunc->instdata.sigma_F);
    g_string_append_printf(buf, "#Number of iterations  = %d  \n", opodrmc->Nmc);
    g_string_append_printf(buf, "\n");

    g_string_append_printf(buf, "#Parameters of Oliver Pharr analysis (ODR fit) \n");
    g_string_append_printf(buf, "#beta: %g \n", opodrdata->beta);
    g_string_append_printf(buf, "#nu: %g                 \n", opodrunc->instdata.nu);
    g_string_append_printf(buf, "#nu_indenter: %g        \n", opodrunc->instdata.nui);
    g_string_append_printf(buf, "#E_indenter: %g      GPa\n", opodrunc->instdata.Ei * 1e-9); //TODO units
    g_string_append_printf(buf, "#Fitted range: ");

    if (opodrdata->Finput) {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", opodrdata->range_pct_Fmax[0], opodrdata->range_pct_Fmax[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] nm \n", opodrdata->range[0], opodrdata->range[1]);
    }

    g_string_append_printf(buf, "# corresponds : ");

    if (opodrdata->Finput) {
        g_string_append_printf(buf, "[%g , %g ] nm \n", opodrdata->range[0], opodrdata->range[1]);
    }
    else {
        g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", opodrdata->range_pct_Fmax[0], opodrdata->range_pct_Fmax[1]);
    }

    g_string_append_printf(buf, "# %d datapoints used.\n", opodrdata->nfitdata);

    if (opodrdata->radial_corr) {
        g_string_append_printf(buf, "#Applied radial correction using α: %.4g°\n", opodrdata->radial_angle);
    }
    else {
        g_string_append_printf(buf, "#Radial correction not applied\n");
    }

    g_string_append_printf(buf, "\n\n");

    g_string_append_printf(buf, "#");

    if (opodrmc->skipmc) {
        g_string_append_printf(buf, "Fit failed in %d cases, the corresponding data area NOT included in the  estimates of mean and uncertainty but ARE included in the histograms.", opodrmc->skipmc);
    }

    g_string_append_printf(buf, "\n\n");

    g_string_append_printf(buf, "#Results \n");

    for (i = 0; i < NPARAM; i++) {
        is = nparam_to_string(i, str, sizeof(str));
        g_string_append_printf(buf, "#%s = (%g ± %g ) \n", str, opodrmc->mcavg[i], opodrmc->mcstd[i]);
    }

    g_string_append_printf(buf, "\n\n\n");

    // Histograms

    if (opodrmc->mcnstat > 0) {
        g_string_append_printf(buf, "#Histogram\n");

        for (i = 0; i < NPARAM; i++) {
            g_string_append_printf(buf, "#Histogram\n");
            is = nparam_to_string(i, str, sizeof(str));

            if (is) {
                g_string_append_printf(buf, "#%s    count/a.u.  \n", str);

                for (j = 0; j < opodrmc->mcnstat; j++) {
                    g_string_append_printf(buf, "%g %g \n", opodrmc->mcxhist[i][j], opodrmc->mcyhist[i][j]);
                }

                g_string_append_printf(buf, "\n \n");
            }
        }
    }

    g_string_append_printf(buf, "\n\n\n");

    //Full data
    g_string_append_printf(buf, "#Full data \n");
    g_string_append_printf(buf, "#");

    for (j = 0; j < NPARAM; j++) {
        is = nparam_to_string(j, str, sizeof(str)); /* _is_ is not used anywhere */
        g_string_append_printf(buf, "%s ", str);
    }

    g_string_append_printf(buf, "\n");

    opodr = opodrmc->opodrmc;

    for (i = 0; i < opodrmc->Nmc; i++) {
        g_string_append_printf(buf, DEPTH, opodr->hc);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, AREA, opodr->Aphc);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, HARD, opodr->Hit);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, MODUL, opodr->Er);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, MODUL, opodr->Eit);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, SLOPE, opodr->S);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, NUMBER, opodr->m);
        g_string_append_printf(buf, "  \n");
        g_string_append_printf(buf, DEPTH, opodr->hp);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, NUMBER, opodr->eps);
        g_string_append_printf(buf, "  \n");
        opodr++;
    }

    return g_string_free(buf, FALSE);
}

void op_odr_uncertainty_montecarlo_run_calc(const OPODRdata *opodrdata, OPODRUncdata *opodrunc, OPODRMCdata *opodrmc,
        FDdata *fddata, const Area *area)
{
    gint istart, iend, ndata;
    gint i, Nmc;
    OPODRdata *opodr;
    OPODRMCdata *mc;
    Instdata *instdata;
    GwyDataLine *h, *F;
    GwyDataLine *hnew, *Fnew;
    GwyDataLine *hnoise, *Fnoise;
    gdouble hmax_old, Fmax_old;
    GRand *rng;
    gint nn;
    gdouble **dd;
    gdouble *xhist, *yhist;
    gboolean hist;
    gdouble estimate[3];

    instdata = &(opodrunc->instdata);

    //not necessary to check range, couldn't get here in the first place

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (opodrdata->Finput)
        range_to_indices(opodrdata->from_pct_Fmax * 0.01 * fddata->Fmax,
                         opodrdata->to_pct_Fmax * 0.01 * fddata->Fmax,
                         fddata->Funload, TRUE, &istart, &iend, &ndata);
    else {
        range_to_indices(opodrdata->from, opodrdata->to, fddata->hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, not fitting can be performed
    ndata = MIN(ndata, fddata->hunload->res - istart);

    if (ndata < 2) {
        return;
    }

    //back up hmax, Fmax
    hmax_old = fddata->hmax;
    Fmax_old = fddata->Fmax;

    //use results of previous fit as input estimates
    estimate[0] = opodrdata->alpha;
    estimate[1] = opodrdata->hp;
    estimate[2] = opodrdata->m;

    //get data to be fitted
    h = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
    F = gwy_data_line_part_extract(fddata->Funload, istart, ndata);

    Nmc = opodrmc->Nmc;
    opodr = opodrmc->opodrmc;

    //prepare data
    for (i = 0; i < Nmc; i++) {
        if (verbose) {
            g_print("i %d \n", i);
        }

        init_OPODRdata(opodr);

        if (verbose) {
            g_print("initialized %s \n", opodr->logfnm);
        }

        g_free(opodr->logfnm);
        // send to stdout, which can be redirected to /dev/null
        // otherwise too much IO time
        opodr->logfnm = NULL;
        //opodr->logfnm = create_log_filename_path(mclogfnm);
        opodr++;
    }

    opodr = opodrmc->opodrmc;

    //copy input parameters
    for (i = 0; i < Nmc ; i++) {
        //    	opodr->radial_corr = args->opodrdata.radial_corr;
        //    	opodr->radial_angle = args->opodrdata.radial_angle;
        /* opodr->nu = args->opodrdata.nu; */
        /* opodr->nui = args->opodrdata.nui; */
        /* opodr->Ei = args->opodrdata.Ei; */
        //		opodr->from = args->opodrdata.from;
        //		opodr->to = args->opodrdata.to;
        opodr++;
    }

    /*
       for (i=0; i< Nmc ;i++){
       printf ("nu %g \n", args->opodrmc.opodrmc[i].nu);
       }
       */

    // initialize random generator and datafield
    rng = g_rand_new();

    //prepare datafield for noise
    hnoise = gwy_data_line_new_alike(h, FALSE);
    Fnoise = gwy_data_line_new_alike(h, FALSE);

    //prepare datafield for varied fields
    hnew = gwy_data_line_new_alike(h, FALSE);
    Fnew = gwy_data_line_new_alike(F, FALSE);

    //prepare auxiliary field
    opodr = opodrmc->opodrmc;

    for (i = 0; i < Nmc; i++) {

        //create noise for h-F data
        generate_uncorrelated_noise(hnoise, instdata->sigma_h, rng);
        generate_uncorrelated_noise(Fnoise, instdata->sigma_F, rng);

        //add noise to original h-F data
        if (!sum_data_lines(hnew, h, hnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        if (!sum_data_lines(Fnew, F, Fnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        //add random noise to hmax, Fmax
        fddata->hmax = hmax_old + gaussian_random_number(rng) * instdata->sigma_h;
        fddata->Fmax = Fmax_old + gaussian_random_number(rng) * instdata->sigma_F;

        //run fit
        op_odr_fit(hnew, Fnew, opodr, fddata, instdata, area, estimate);

        if (opodrdata->radial_corr) {
            radial_correction_op(opodrdata->radial_angle, &(opodr->Hit), &(opodr->Er), &(opodr->Eit), &(opodr->Aphc), instdata);
        }

        opodr++;
    }

    //reset hmax, Fmax to original value
    fddata->hmax = hmax_old;
    fddata->Fmax = Fmax_old;

    if (verbose) {
        g_print("run Monte Carlo \n");
    }

    g_object_unref(h);
    g_object_unref(F);
    g_object_unref(hnew);
    g_object_unref(Fnew);
    g_object_unref(hnoise);
    g_object_unref(Fnoise);
    g_rand_free(rng);

    for (i = 0; i < Nmc; i++) {
        //		printf(" %d Er %g Eit %g Hit %g \n", i, args->opodrmc.opodrmc[i].Er, args->opodrmc.opodrmc[i].Eit, args->opodrmc.opodrmc[i].Hit);
    }

    //calculate averages and stdevs
    //hp, Aphc, Hit, Er, Eit, S, m

    mc = opodrmc;

    for (i = 0; i < NPARAM; i++) {
        mc->mcavg[i] = 0;
        mc->mcstd[i] = 0;
    }

    opodr = opodrmc->opodrmc;
    nn = 0;

    for (i = 0; i < Nmc; i++) {
        //check sanity
        if (opodr->hc > 0) {
            mc->mcavg[0] += opodr->hc;
            mc->mcavg[1] += opodr->Aphc;
            mc->mcavg[2] += opodr->Hit;
            mc->mcavg[3] += opodr->Er;
            mc->mcavg[4] += opodr->Eit;
            mc->mcavg[5] += opodr->S;
            mc->mcavg[6] += opodr->m;
            mc->mcavg[7] += opodr->hp;
            mc->mcavg[8] += opodr->eps;
            mc->mcstd[0] += opodr->hc * opodr->hc;
            mc->mcstd[1] += opodr->Aphc * opodr->Aphc;
            mc->mcstd[2] += opodr->Hit * opodr->Hit;
            mc->mcstd[3] += opodr->Er * opodr->Er;
            mc->mcstd[4] += opodr->Eit * opodr->Eit;
            mc->mcstd[5] += opodr->S * opodr->S;
            mc->mcstd[6] += opodr->m * opodr->m;
            mc->mcstd[7] += opodr->hp * opodr->hp;
            mc->mcstd[8] += opodr->eps * opodr->eps;
            nn++;
        }

        opodr++;
    }

    //use only reasonable data for average and standard deviation
    for (i = 0; i < NPARAM; i++) {
        mc->mcavg[i] /= nn;
        mc->mcstd[i] /= nn;
        mc->mcstd[i] -= mc->mcavg[i] * mc->mcavg[i];
        mc->mcstd[i] = (mc->mcstd[i] < 0) ? 0 : mc->mcstd[i];
        mc->mcstd[i] = sqrt(mc->mcstd[i]);
    }

    // how many datasets were skipped
    opodrmc->skipmc = Nmc - nn;

    /*
       printf("mean  stdev \n");
       for (i = 0; i< NPARAM; i++){
       printf(" %g   %g \n", avg[i], std[i]);
       }
       */

    // Histograms
    //create data structure

    dd = opodrmc->mcdd;
    opodr = opodrmc->opodrmc;

    for (i = 0; i < Nmc; i++) {
        dd[0][i] = opodr->hc;
        dd[1][i] = opodr->Aphc;
        dd[2][i] = opodr->Hit;
        dd[3][i] = opodr->Er;
        dd[4][i] = opodr->Eit;
        dd[5][i] = opodr->S;
        dd[6][i] = opodr->m;
        dd[7][i] = opodr->hp;
        dd[8][i] = opodr->eps;
        opodr++;
    }

    //prepare histogram structure

    for (i = 0; i < NPARAM; i++) {
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

gboolean op_odr_unc_has_results(const OPODRUncdata *opodrunc)
{
    return (opodrunc->Ec != NULL && opodrunc->Hc != NULL && opodrunc->Ac != NULL && opodrunc->hc != NULL);
}

#endif

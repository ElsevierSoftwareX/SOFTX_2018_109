#ifndef NOFORTRAN

#include "tool_hertz_odr_unc.h"
#include "tool_hertz_odr.h"

#include "datatypes.h"
#include "fddata.h"
#include "fit-utils.h"
#include "mc-utils.h"
#include "niget-common.h"

#include <math.h>


#define N_FIX 1.5

static const gchar *logfnm = "fit.log.hz.unc";
//static const gchar *mclogfnm = "fit.log.hz.mc";

#define RADIALCORRECTION FALSE

static void init_shifted_HertzODRdata(const HertzODRdata *src, HertzODRdata *dest, gdouble dh);
static gboolean nparam_to_string(gint i, gchar *str, gint nstr);

void init_HertzODRUncdata(HertzODRUncdata *hzodrunc, Instdata *instdata)
{
    hzodrunc->instdata = *instdata;
    /*
      hzodrunc->unu = 0.05;
      hzodrunc->unui = 0.04;
      hzodrunc->uEi = 20e9;
    */
    hzodrunc->uradius = 10e-9;
    hzodrunc->uEr = 1;
    hzodrunc->uEit = 1;
    hzodrunc->Nmc = N_MC;

    hzodrunc->uEruh = 0;
    hzodrunc->uEruF = 0;
    hzodrunc->uEruradius = 0;
    hzodrunc->uErtotal = 0;

    hzodrunc->uEituh = 0;
    hzodrunc->uEituF = 0;
    hzodrunc->uEituradius = 0;
    hzodrunc->uEitunu = 0;
    hzodrunc->uEitunui = 0;
    hzodrunc->uEituEi = 0;
    hzodrunc->uEittotal = 0;

    hzodrunc->uradiusuh = 0;
    hzodrunc->uradiusuF = 0;
    hzodrunc->uradiusunu = 0;
    hzodrunc->uradiusunui = 0;
    hzodrunc->uradiusuEi = 0;
    hzodrunc->uradiustotal = 0;

    hzodrunc->w = NULL;
    hzodrunc->ERc = NULL;
}

static void init_shifted_HertzODRdata(const HertzODRdata *src, HertzODRdata *dest, gdouble dh)
{
    init_HertzODRdata(dest);

    // shift integration range by dh

    dest->range[0] = src->range[0] - dh;
    dest->range[1] = src->range[1] - dh;
    dest->from = dest->range[0];
    dest->to = dest->range[1];

    dest->mode = src->mode;

    switch (dest->mode) {
    case R_MODE:
        dest->radius = src->radius;
        break;

    case ER_MODE:
        dest->Er = src->Er;
        break;

    case EIT_MODE:
        dest->Eit = src->Eit;
        break;
    }

    dest->ifixb[0] = src->ifixb[0];
    dest->ifixb[1] = src->ifixb[1];
    dest->ifixb[2] = src->ifixb[2];

    dest->radial_corr = src->radial_corr;

    dest->logfnm = create_log_filename_path(logfnm);
}

void hertz_odr_propagate_uncertainties_calc(HertzODRdata *hzodr, HertzODRUncdata *unc, const FDdata *fd)
{
    gdouble R, Er, Eit, nu, nui, Ei;
    gdouble uErunu, uErunui, uEruEi, uEruEit;
    gdouble uEi, unu, unui;
    gdouble uradius, uEr, uEit;
    gdouble dradiusdEr ;
    gchar *unc_logfnm;
    gdouble beta[3], std_beta[3], covariance[3];
    gint istart, iend, ndata, info, nremove;
    GwyDataLine *x, *y;
    gdouble ux, uy;
    gdouble uh, uF;
    gdouble ugammah, ugammaF, ugamma;
    gdouble unh, unF, un;
    gdouble SSres, R2, R2adj, chi2;
    gint np;

    unc_logfnm = create_log_filename_path(logfnm);

    if (verbose) {
        g_print("logfnm %s \n", unc_logfnm);
    }

    if (verbose) {
        g_print(" \n\nPropagate uncertainties\n");
    }

    uh = unc->instdata.sigma_h;
    uF = unc->instdata.sigma_F;
    Ei = unc->instdata.Ei;
    nui = unc->instdata.nui;
    nu = unc->instdata.nu;
    uEi = unc->instdata.uEi;
    unui = unc->instdata.unui;
    unu = unc->instdata.unu;
    uradius = unc->uradius;
    uEr = unc->uEr;
    uEit = unc->uEit;

    /* prepare odr fit */

    beta[0] = hzodr->gamma;
    beta[1] = hzodr->h0;
    beta[2] = hzodr->n;
    np = 3;


    //find appropriate data, which should be fitted
    range_to_indices(hzodr->range[0], hzodr->range[1], fd->hload, FALSE, &istart, &iend, &ndata);

    x = gwy_data_line_part_extract(fd->hload, istart, ndata);
    y = gwy_data_line_part_extract(fd->Fload, istart, ndata);
    //remove negative depths from data
    nremove = filter_negative_data(x, y);
    hzodr->nfitdata -= nremove;

    /* fit with only uh */
    ux = uh;
    uy = 0;

    if (hzodr->radial_corr) {
        info = odr_run_fit(x, y, np, beta, hzodr->ifixb, unc_logfnm, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, HERTZ_RADIAL);
    }
    else {
        info = odr_run_fit(x, y, np, beta, hzodr->ifixb, unc_logfnm, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);
    }

    unh = std_beta[2];
    ugammah = std_beta[0];

    if (verbose)
        g_print(" fit with only uh:\n info %d gamma %g n %g h0 %g \n u(gamma) %g u(n) %g u(h0) %g \n",
                info, beta[0], beta[2], beta[1], ugammah, unh, std_beta[1]);

    /* fit with only uF */
    ux = 0;
    uy = uF;

    if (hzodr->radial_corr) {
        info = odr_run_fit(x, y, np, beta, hzodr->ifixb, unc_logfnm, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, HERTZ_RADIAL);
    }
    else {
        info = odr_run_fit(x, y, np, beta, hzodr->ifixb, unc_logfnm, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);
    }

    unF = std_beta[2];
    ugammaF = std_beta[0];

    if (verbose)
        g_print(" fit with only uF:\n info %d gamma %g n %g h0 %g \n  u(gamma) %g u(n) %g u(h0) %g \n",
                info, beta[0], beta[2], beta[1], ugammaF, unF, std_beta[1]);

    unc->unuh = unh;
    unc->unuF = unF;
    unc->ugammauh = ugammah;
    unc->ugammauF = ugammaF;

    // total uncertainties
    // should take as given by odr_wrap, otherwise doublecounting
    /* fit with both uh and uF */

    ux = uh;
    uy = uF;

    if (hzodr->radial_corr) {
        info = odr_run_fit(x, y, np, beta, hzodr->ifixb, unc_logfnm, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, HERTZ_RADIAL);
    }
    else {
        info = odr_run_fit(x, y, np, beta, hzodr->ifixb, unc_logfnm, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, POWER_3PARAM);
    }

    un = std_beta[2];
    ugamma = std_beta[0];

    unc->un = un;
    unc->ugamma = ugamma;

    if (verbose)
        g_print(" fit with both uh, uF:\n info %d gamma %g n %g h0 %g \n  u(gamma) %g u(n) %g u(h0) %g \n",
                info, beta[0], beta[2], beta[1], ugamma, un, std_beta[1]);

    switch (hzodr->mode) {
    case R_MODE:
        //make copies to simplify code
        // input
        Ei = unc->instdata.Ei;
        nui = unc->instdata.nui;
        nu = unc->instdata.nu;
        R = hzodr->radius;
        //output
        Eit = hzodr->Eit;
        Er = hzodr->Er;

        unc->uEruh = ugammah / hzodr->gamma * Er * uh; //Er/a uah
        unc->uEruF = ugammaF / hzodr->gamma * Er * uF; //Er/a uaF
        unc->uEruradius = 0.5 * Er / R * uradius; // -1/2 *Er/R u(R)

        //uEituh = dEit/ dh uh = dEit / dEr * dEr/dh uh = dEit/dEr * uEruh
        //dEit/dEr = Eit^2/Er^2/ (1-nu^2)

        unc->uEituh = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruh;
        unc->uEituF = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruF;
        unc->uEituradius = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruradius;

        unc->uEitunu =  fabs(2 * nu / (1 - sq(nu)) * Eit * unu); // -2nu/(1-nu^2) Eit u(nu)
        unc->uEitunui = fabs(2 * nui / (1 - sq(nu)) * sq(Eit) / Ei * unui) * 1e9; // -2nui/(1-nu^2) Eit^2/Ei u(nui)
        unc->uEituEi =  fabs((1 - sq(nui)) / (1 - sq(nu)) * sq(Eit) / sq(Ei) * uEi) * 1e9; //-(1-nui^2)/(1-nu)^2 Eit^2/Ei^2 u(Ei)

        // total uncertainties
        unc->uErtotal = sqrt(sq(unc->uEruh) + sq(unc->uEruF) + sq(unc->uEruradius));
        unc->uEittotal = sqrt(sq(unc->uEituh) + sq(unc->uEituF) + sq(unc->uEituradius) + sq(unc->uEitunu) + sq(unc->uEitunui) + sq(unc->uEituEi));
        break;

    case ER_MODE:
        //make copies to simplify code
        // input
        Er = hzodr->Er;
        // output
        R = hzodr->radius;

        unc->uradiusuh = R * 2 * ugammah / hzodr->gamma; // 2 R/a uah
        unc->uradiusuF = R * 2 * ugammaF / hzodr->gamma; // 2 R/a uaF
        unc->uradiusuEr = 2 * R / Er * uEr; // -2 *R/Er u(Er)

        // total uncertainties
        unc->uradiustotal = sqrt(sq(unc->uradiusuh) + sq(unc->uradiusuF) + sq(unc->uradiusuEr));
        break;

    case EIT_MODE :
        //make copies to simplify code
        // input
        Eit = hzodr->Eit;
        //output
        R = hzodr->radius;
        Er = hzodr->Er;

        uErunu =  2 * nu * sq(Er) / Eit * unu; // 2*nu/EIT Er^2 unu
        uErunui =  2 * nui * sq(Er) / Ei * unui * 1e9; // 2*nui/Ei Er^2 unui
        uEruEi = (1 - sq(nui)) / sq(Ei) * sq(Er) * uEi * 1e9; // (1-nui^2)/Ei^2 Er^2 uEi
        uEruEit = (1 - sq(nu)) / sq(Eit) * sq(Er) * uEit ; // (1-nu^2)/EIT^2 Er^2 uEIT

        unc->uradiusuh = R * 2 * ugammah / hzodr->gamma; // 2 R/a uah
        unc->uradiusuF = R * 2 * ugammaF / hzodr->gamma; // 2 R/a uaF
        dradiusdEr = 2 * R / Er; // dR/dEr = -2 *R/Er

        unc->uradiusunu = dradiusdEr * uErunu;
        unc->uradiusunui = dradiusdEr * uErunui;
        unc->uradiusuEi = dradiusdEr * uEruEi;
        unc->uradiusuEit = dradiusdEr * uEruEit;

        // total uncertainties
        unc->uradiustotal = sqrt(sq(unc->uradiusuh) + sq(unc->uradiusuF) + sq(unc->uradiusuEit) + sq(unc->uradiusuEi) + sq(unc->uradiusunu) + sq(unc->uradiusunui));
        break;
    }

    g_object_unref(x);
    g_object_unref(y);
    g_free(unc_logfnm);
}

gchar *hertz_odr_uncertainty_export_data(const HertzODRdata *hzodrdata, const HertzODRUncdata *hzodrunc,
        const FDdata *fddata, enum ExportFormat exportformat)
{
    gint j;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Parameters of HertzODR analysis \n");
    g_string_append_printf(buf, "# u(depth): %g nm\n", hzodrunc->instdata.sigma_h);
    g_string_append_printf(buf, "# u(load): %g nm\n", hzodrunc->instdata.sigma_F);
    g_string_append_printf(buf, "#fitted range:  [%g, %g ] nm \n", hzodrdata->range[0], hzodrdata->range[1]);

    if (hzodrdata->ifixb[2] == 0) {
        g_string_append_printf(buf, "#n: %g fixed \n", N_FIX);
    }

    if (hzodrdata->radial_corr) {
        g_string_append_printf(buf, "Fitted using radial correction\n");
    }

    switch (hzodrdata->mode) {
    case R_MODE:
        g_string_append_printf(buf, "#R_mode \n");
        g_string_append_printf(buf, "#tip radius: %g nm       u(tip radius): %g nm\n", hzodrdata->radius * 1e9, hzodrunc->uradius * 1e9);
        g_string_append_printf(buf, "#nu: %g                 u(nu):  %g  \n", hzodrunc->instdata.nu, hzodrunc->instdata.unu);
        g_string_append_printf(buf, "#nu_indenter: %g       u(nu_indenter):  %g   \n", hzodrunc->instdata.nui, hzodrunc->instdata.unui);
        g_string_append_printf(buf, "#E_indenter: %g        u(E_indenter):   %g GPa\n", hzodrunc->instdata.Ei * 1e-9, hzodrunc->instdata.uEi * 1e-9); //TODO units
        g_string_append_printf(buf, "#\n");

        g_string_append_printf(buf, "#Results of HertzODR analysis \n");
        g_string_append_printf(buf, "#Er: %g GPa\n", hzodrdata->Er);
        g_string_append_printf(buf, "#Eit: %g GPa\n", hzodrdata->Eit);
        break;

    case ER_MODE:
        g_string_append_printf(buf, "#Er_mode \n");
        g_string_append_printf(buf, "#Er: %g GPa            u(Er) : %g GPa\n", hzodrdata->Er, hzodrunc->uEr);
        g_string_append_printf(buf, "#\n");

        g_string_append_printf(buf, "#Results of HertzODR analysis \n");
        g_string_append_printf(buf, "#radius: %g nm\n", hzodrdata->radius * 1e9);
        break;

    case EIT_MODE:
        g_string_append_printf(buf, "#Eit_mode \n");
        g_string_append_printf(buf, "#Eit: %g GPa            u(Eit) : %g GPa\n", hzodrdata->Eit, hzodrunc->uEit);
        g_string_append_printf(buf, "#nu: %g                 u(nu):  %g  \n", hzodrunc->instdata.nu, hzodrunc->instdata.unu);
        g_string_append_printf(buf, "#nu_indenter: %g       u(nu_indenter):  %g   \n", hzodrunc->instdata.nui, hzodrunc->instdata.unui);
        g_string_append_printf(buf, "#E_indenter: %g        u(E_indenter):   %g GPa\n", hzodrunc->instdata.Ei * 1e-9, hzodrunc->instdata.uEi * 1e-9); //TODO units
        g_string_append_printf(buf, "#\n");

        g_string_append_printf(buf, "#Results of HertzODR analysis \n");
        g_string_append_printf(buf, "#radius: %g nm\n", hzodrdata->radius * 1e9);
        break;
    }

    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");

    /* Gaussian propagation */
    g_string_append_printf(buf, "#Results of the Gaussian propagation of uncertainties \n");

    switch (hzodrdata->mode) {
    case R_MODE:
        g_string_append_printf(buf, "#u(Er) =  %g GPa \n", hzodrunc->uErtotal);
        g_string_append_printf(buf, "#u(Eit) =  %g GPa \n", hzodrunc->uEittotal);
        g_string_append_printf(buf, "#\n");
        g_string_append_printf(buf, "#Contributions to uncertainty of Er\n");
        g_string_append_printf(buf, "#u(Er;uh) = %g GPa \n", hzodrunc->uEruh);
        g_string_append_printf(buf, "#u(Er;uF) = %g GPa \n", hzodrunc->uEruF);
        g_string_append_printf(buf, "#u(Er;uradius) = %g GPa \n", hzodrunc->uEruradius);
        g_string_append_printf(buf, "#\n");
        g_string_append_printf(buf, "#Contributions to uncertainty of Eit\n");
        g_string_append_printf(buf, "#u(Eit;uh) = %g GPa \n", hzodrunc->uEituh);
        g_string_append_printf(buf, "#u(Eit;uF) = %g GPa \n", hzodrunc->uEituF);
        g_string_append_printf(buf, "#u(Eit;uradius) = %g GPa \n", hzodrunc->uEituradius);
        g_string_append_printf(buf, "#u(Eit;unu) = %g GPa \n", hzodrunc->uEitunu);
        g_string_append_printf(buf, "#u(Eit;unui) = %g GPa \n", hzodrunc->uEitunui);
        g_string_append_printf(buf, "#u(Eit;uEi) = %g GPa \n", hzodrunc->uEituEi);
        g_string_append_printf(buf, "#\n");
        break;

    case ER_MODE:
        g_string_append_printf(buf, "#u(radius) =  %g nm \n", hzodrunc->uradiustotal * 1e9);
        g_string_append_printf(buf, "#\n");
        g_string_append_printf(buf, "#Contributions to uncertainty of radius\n");
        g_string_append_printf(buf, "#u(radius;uh) = %g nm \n", hzodrunc->uradiusuh * 1e9);
        g_string_append_printf(buf, "#u(radius;uF) = %g nm \n", hzodrunc->uradiusuF * 1e9);
        g_string_append_printf(buf, "#u(uradius;uEr) = %g nm \n", hzodrunc->uradiusuEr * 1e9);
        g_string_append_printf(buf, "#\n");
        break;

    case EIT_MODE:
        g_string_append_printf(buf, "#u(radius) =  %g nm \n", hzodrunc->uradiustotal * 1e9);
        g_string_append_printf(buf, "#\n");
        g_string_append_printf(buf, "#Contributions to uncertainty of radius\n");
        g_string_append_printf(buf, "#u(radius;uh) = %g nm \n", hzodrunc->uradiusuh * 1e9);
        g_string_append_printf(buf, "#u(radius;uF) = %g nm \n", hzodrunc->uradiusuF * 1e9);
        g_string_append_printf(buf, "#u(uradius;uEit) = %g nm \n", hzodrunc->uradiusuEit * 1e9);
        g_string_append_printf(buf, "#u(uradius;unu) = %g nm \n", hzodrunc->uradiusunu * 1e9);
        g_string_append_printf(buf, "#u(uradius;unui) = %g nm \n", hzodrunc->uradiusunui * 1e9);
        g_string_append_printf(buf, "#u(uradius;uEi) = %g nm \n", hzodrunc->uradiusuEi * 1e9);
        break;
    }

    g_string_append_printf(buf, "#\n");

    /* contact point uncertainties */
    switch (hzodrdata->mode) {
    case R_MODE:
        g_string_append_printf(buf, "#Changes of the Er corresponding to changes in the contact point \n");
        g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   ");
        g_string_append_printf(buf, "Er/GPa \n");

        for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
            if (hzodrunc->ERc[j + CONTACTMAX] < 0) {
                g_string_append_printf(buf, "# %d   not defined \n", j);
            }
            else {
                g_string_append_printf(buf, "# %d   %g  \n", j, hzodrunc->ERc[j + CONTACTMAX]);
            }
        }

        break;

    case ER_MODE:
    case EIT_MODE:
        g_string_append_printf(buf, "#Changes of the radius corresponding to changes in the contact point \n");
        g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   ");
        g_string_append_printf(buf, "R/nm \n");

        for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
            if (hzodrunc->ERc[j + CONTACTMAX] < 0) {
                g_string_append_printf(buf, "# %d   not defined \n", j);
            }
            else {
                g_string_append_printf(buf, "# %d   %g  \n", j, hzodrunc->ERc[j + CONTACTMAX] * 1e9);
            }
        }

        break;
    }

    return g_string_free(buf, FALSE);
}

void hertz_odr_fit_shift_contact(const HertzODRdata *hzodrdata, const FDdata *fddata, const Instdata *instdata, gdouble *ERc, gint j)
{
    HertzODRdata hzodr;
    FDdata fdc;
    gdouble dh;
    gint istart, iend, ndata;
    gint nremove;
    GwyDataLine *x, *y;

    ERc[j + CONTACTMAX] = -1;

    if (fddata->i_contact_load + j < 0) {
        return;
    }

    // create shifted FDdata

    if (verbose) {
        g_print("\n\nCreate shifted data for j %d \n", j);
    }

    dh = create_shifted_FDdata(fddata, &fdc, j);
    /*
      fdc.i_contact_load = fddata->i_contact_load+j;

      fdc.hall_orig = gwy_data_line_duplicate(fddata->hall_orig);
      fdc.Fall_orig = gwy_data_line_duplicate(fddata->Fall_orig);


      fdc.hload = gwy_data_line_part_extract(fdc.hall_orig, fdc.i_contact_load, fddata->hload->res-j);
      fdc.Fload = gwy_data_line_part_extract(fdc.Fall_orig, fdc.i_contact_load, fddata->Fload->res-j);
      gwy_data_line_add(fdc.hload, -fdc.hall_orig->data[fdc.i_contact_load]);
      gwy_data_line_add(fdc.Fload, -fdc.Fall_orig->data[fdc.i_contact_load]);

      //difference in shift between different contactpoints
      dh = fdc.hall_orig->data[fdc.i_contact_load]-fddata->hall_orig->data[fddata->i_contact_load];
    */

    if (verbose) {
        g_print("shift dh %g  \n", dh);
    }

    // create working copy of HertzODRdat
    init_shifted_HertzODRdata(hzodrdata, &hzodr, dh);

    //check range

    if (verbose) {
        g_print(" hzodr.from %g to %g \n", hzodr.from, hzodr.to);
    }

    if (hzodr.from == hzodr.to) {
        return;
    }

    //find appropriate data, which should be fitted
    range_to_indices(hzodr.from, hzodr.to, fdc.hload, FALSE, &istart, &iend, &ndata);
    // if there are not enough data in the selected range, not fitting can be performed
    ndata = MIN(ndata, fdc.hload->res - istart);

    if (verbose)
        g_print("interval has %d points, from %d to %d, hrange [%g, %g ] nm, Frange [%g, %g] mN \n", ndata, istart, iend,
                fdc.hload->data[istart], fdc.hload->data[iend], fdc.Fload->data[istart], fdc.Fload->data[iend]);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    x = gwy_data_line_part_extract(fdc.hload, istart, ndata);
    y = gwy_data_line_part_extract(fdc.Fload, istart, ndata);
    nremove = filter_negative_data(x, y);
    hzodr.nfitdata = ndata - nremove;

    hertz_odr_fit(x, y, &hzodr, &fdc, instdata, NULL);

    // clean up
    g_free(hzodr.logfnm);
    g_object_unref(x);
    g_object_unref(y);
    destroy_FDdata(&fdc);

    if (hzodr.mode == R_MODE) {
        ERc[j + CONTACTMAX] = hzodr.Er;
    }
    else {
        ERc[j + CONTACTMAX] = hzodr.radius;
    }

    if (verbose) {
        g_print("\n\n");
    }
}

void init_HertzODRMCdata(HertzODRMCdata *hzodrmc, gint Nmc)
{
    gint i;

    //number of iterations
    hzodrmc->Nmc = Nmc;

    if (verbose) {
        g_print("initialize for %d iterations \n", hzodrmc->Nmc);
    }

    //how much had to be skipped
    hzodrmc->skipmc = 0;

    //dynamically allocate
    hzodrmc->hertzodrmc = (HertzODRdata *)g_malloc(hzodrmc->Nmc * sizeof(HertzODRdata));
    hzodrmc->mcavg = (gdouble *)g_malloc(NPARAM_HZODR * sizeof(gdouble));
    hzodrmc->mcstd = (gdouble *)g_malloc(NPARAM_HZODR * sizeof(gdouble));
    hzodrmc->mcdd = (gdouble **) g_malloc(NPARAM_HZODR * sizeof(gdouble *));
    hzodrmc->mcxhist = (gdouble **) g_malloc(NPARAM_HZODR * sizeof(gdouble *));
    hzodrmc->mcyhist = (gdouble **) g_malloc(NPARAM_HZODR * sizeof(gdouble *));

    hzodrmc->mcnstat = (gint)floor(3.49 * cbrt(hzodrmc->Nmc) + 0.5);

    for (i = 0; i < NPARAM_HZODR; i++) {
        hzodrmc->mcdd[i] = (gdouble *) g_malloc(hzodrmc->Nmc * sizeof(gdouble));
        hzodrmc->mcxhist[i] = (gdouble *) g_malloc(hzodrmc->mcnstat * sizeof(gdouble));
        hzodrmc->mcyhist[i] = (gdouble *) g_malloc(hzodrmc->mcnstat * sizeof(gdouble));
    }
}

void destroy_HertzODRMCdata(HertzODRMCdata *hzodrmc)
{
    gint i;

    //free everything
    g_free(hzodrmc->mcavg);
    g_free(hzodrmc->mcstd);

    for (i = 0; i < NPARAM_HZODR; i++) {
        g_free(hzodrmc->mcdd[i]);
        g_free(hzodrmc->mcxhist[i]);
        g_free(hzodrmc->mcyhist[i]);
    }

    g_free(hzodrmc->mcdd);
    g_free(hzodrmc->mcxhist);
    g_free(hzodrmc->mcyhist);

    hzodrmc->mcnstat = 0;
    hzodrmc->skipmc = 0;

    for (i = 0; i < hzodrmc->Nmc; i++) {
        g_free(hzodrmc->hertzodrmc[i].logfnm);
    }

    g_free(hzodrmc->hertzodrmc);
}

static gboolean nparam_to_string(gint i, gchar *str, gint nstr)
{
    switch (i) {
    case 0:
        g_snprintf(str, nstr, "E_r/GPa");
        break;

    case 1:
        g_snprintf(str, nstr, "E_IT/GPa");
        break;

    case 2:
        g_snprintf(str, nstr, "radius/nm");
        break;

    default:
        g_printerr("Unhandled number of results.\n");
        return FALSE;
    }

    return TRUE;
}

gchar *hertz_odr_mc_export_data(const HertzODRdata *hertzodrdata, const HertzODRUncdata *hertzodrunc, const HertzODRMCdata *hertzodrmc,
                                const FDdata *fddata, enum ExportFormat exportformat)
{
    HertzODRdata *hzodr;
    gint i, j;
    gchar str[200];
    gboolean is;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Monte Carlo simulation of HertzODR analysis\n");
    g_string_append_printf(buf, "#Input uncertainties \n");
    g_string_append_printf(buf, "#u(h) = %g nm \n", hertzodrunc->instdata.sigma_h);
    g_string_append_printf(buf, "#u(F) = %g mN \n", hertzodrunc->instdata.sigma_F);
    g_string_append_printf(buf, "#Number of iterations  = %d  \n", hertzodrmc->Nmc);
    g_string_append_printf(buf, "\n");

    g_string_append_printf(buf, "#Parameters of HertzODR analysis \n");

    switch (hertzodrdata->mode) {
    case R_MODE:
        g_string_append_printf(buf, "#tip radius: %g nm\n", hertzodrdata->radius * 1e9);
        g_string_append_printf(buf, "#nu: %g                 \n", hertzodrunc->instdata.nu);
        g_string_append_printf(buf, "#nu_indenter: %g        \n", hertzodrunc->instdata.nui);
        g_string_append_printf(buf, "#E_indenter: %g      GPa\n", hertzodrunc->instdata.Ei * 1e-9); //TODO units
        g_string_append_printf(buf, "#fitted range:  [%g, %g ] nm \n", hertzodrdata->range[0], hertzodrdata->range[1]);
        break;

    case ER_MODE:
        g_string_append_printf(buf, "#Er: %g      GPa\n", hertzodrdata->Er);  //TODO units
        g_string_append_printf(buf, "#fitted range:  [%g, %g ] nm \n", hertzodrdata->range[0], hertzodrdata->range[1]);
        break;

    case EIT_MODE:
        g_string_append_printf(buf, "#Eit: %g      GPa\n", hertzodrdata->Eit);  //TODO units
        g_string_append_printf(buf, "#nu: %g                 \n", hertzodrunc->instdata.nu);
        g_string_append_printf(buf, "#nu_indenter: %g        \n", hertzodrunc->instdata.nui);
        g_string_append_printf(buf, "#E_indenter: %g      GPa\n", hertzodrunc->instdata.Ei * 1e-9); //TODO units
        g_string_append_printf(buf, "#fitted range:  [%g, %g ] nm \n", hertzodrdata->range[0], hertzodrdata->range[1]);
    }

    if (hertzodrdata->ifixb[2] == 0) {
        g_string_append_printf(buf, "#n: %g fixed \n", N_FIX);
    }

    g_string_append_printf(buf, "#");

    if (hertzodrmc->skipmc) {
        g_string_append_printf(buf, "Fit failed in %d cases, the corresponding data area NOT included in the  estimates of mean and uncertainty but ARE included in the histograms.",
                               hertzodrmc->skipmc);
    }

    g_string_append_printf(buf, "\n");
    g_string_append_printf(buf, "#\n\n");

    g_string_append_printf(buf, "#Results \n");

    switch (hertzodrdata->mode) {
    case R_MODE:
        for (i = 0; i < 2; i++) {
            is = nparam_to_string(i, str, sizeof(str));
            g_string_append_printf(buf, "#%s = (%g ± %g ) \n", str, hertzodrmc->mcavg[i], hertzodrmc->mcstd[i]);
        }

        break;

    case ER_MODE:
    case EIT_MODE:
        is = nparam_to_string(2, str, sizeof(str));
        g_string_append_printf(buf, "#%s = (%g ± %g ) \n", str, hertzodrmc->mcavg[2], hertzodrmc->mcstd[2]);
        break;
    }

    g_string_append_printf(buf, "\n\n\n");

    // Histograms

    if (hertzodrmc->mcnstat > 0) {
        g_string_append_printf(buf, "#Histogram\n");

        switch (hertzodrdata->mode) {
        case R_MODE:
            for (i = 0; i < 2; i++) {
                g_string_append_printf(buf, "#Histogram\n");
                is = nparam_to_string(i, str, sizeof(str));

                if (is) {
                    g_string_append_printf(buf, "#%s    count/a.u.  \n", str);

                    for (j = 0; j < hertzodrmc->mcnstat; j++) {
                        g_string_append_printf(buf, "%g %g \n", hertzodrmc->mcxhist[i][j], hertzodrmc->mcyhist[i][j]);
                    }

                    g_string_append_printf(buf, "\n \n");
                }
            }

            break;

        case ER_MODE:
        case EIT_MODE:
            i = 2;
            g_string_append_printf(buf, "#Histogram\n");
            is = nparam_to_string(i, str, sizeof(str));

            if (is) {
                g_string_append_printf(buf, "#%s    count/a.u.  \n", str);

                for (j = 0; j < hertzodrmc->mcnstat; j++) {
                    g_string_append_printf(buf, "%g %g \n", hertzodrmc->mcxhist[i][j], hertzodrmc->mcyhist[i][j]);
                }

                g_string_append_printf(buf, "\n \n");
            }
        }
    }

    g_string_append_printf(buf, "\n\n\n");

    //Full data
    g_string_append_printf(buf, "#Full data \n");
    g_string_append_printf(buf, "#");

    switch (hertzodrdata->mode) {
    case R_MODE:
        for (j = 0; j < 2; j++) {
            is = nparam_to_string(j, str, sizeof(str));
            g_string_append_printf(buf, "%s ", str);
        }

        break;

    case ER_MODE:
    case EIT_MODE:
        is = nparam_to_string(2, str, sizeof(str));
        g_string_append_printf(buf, "%s ", str);
        break;
    }

    g_string_append_printf(buf, "\n");

    hzodr = hertzodrmc->hertzodrmc;

    if (hertzodrdata->mode == R_MODE) {
        for (i = 0; i < hertzodrmc->Nmc; i++) {
            g_string_append_printf(buf, MODUL, hzodr->Er);
            g_string_append_printf(buf, "  ");
            g_string_append_printf(buf, MODUL, hzodr->Eit);
            g_string_append_printf(buf, "  \n");
            hzodr++;
        }
    }
    else {
        for (i = 0; i < hertzodrmc->Nmc; i++) {
            g_string_append_printf(buf, DEPTH, hzodr->radius * 1e9);
            g_string_append_printf(buf, "  \n");
            hzodr++;
        }
    }

    return g_string_free(buf, FALSE);
}

void hertz_odr_uncertainty_montecarlo_run_calc(const HertzODRdata *hertzodrdata, HertzODRUncdata *hertzodrunc, HertzODRMCdata *hertzodrmc,
        const FDdata *fddata)
{
    gint istart, iend, ndata;
    gint i, j, Nmc;
    HertzODRdata *hzodr;
    HertzODRMCdata *mc;
    Instdata *instdata;
    GwyDataLine *h, *F;
    GwyDataLine *hnew, *Fnew;
    GwyDataLine *hnoise, *Fnoise;
    GRand *rng;
    gint nn, nremove;
    gdouble **dd;
    gdouble *xhist, *yhist;
    gboolean hist;
    gdouble estimate[3];

    instdata = &(hertzodrunc->instdata);

    //not necessary to check range, couldn't get here in the first place

    //find appropriate data, which should be fitted, must use real fitted range
    range_to_indices(hertzodrdata->range[0], hertzodrdata->range[1], fddata->hload, FALSE, &istart, &iend, &ndata);
    // if there are not enough data in the selected range, not fitting can be performed
    ndata = MIN(ndata, fddata->hload->res - istart);

    if (ndata < 2) {
        return;
    }

    //use results of previous fit as input estimates
    estimate[0] = hertzodrdata->gamma;
    estimate[1] = hertzodrdata->h0;
    estimate[2] = hertzodrdata->n;

    //get data to be fitted
    h = gwy_data_line_part_extract(fddata->hload, istart, ndata);
    F = gwy_data_line_part_extract(fddata->Fload, istart, ndata);

    Nmc = hertzodrmc->Nmc;
    hzodr = hertzodrmc->hertzodrmc;

    //prepare data
    for (i = 0; i < Nmc ; i++) {
        init_HertzODRdata(hzodr);
        g_free(hzodr->logfnm);
        // send to stdout, which can be redirected to /dev/null
        // otherwise too much IO time
        hzodr->logfnm = NULL;
        //        hzodr->logfnm = create_log_filename_path(mclogfnm);
        hzodr++;
    }

    if (verbose) {
        g_print("copy params \n");
    }

    hzodr = hertzodrmc->hertzodrmc;

    //copy input parameters
    for (i = 0; i < Nmc ; i++) {
        hzodr->mode = hertzodrdata->mode;

        for (j = 0; j < 3; j++) {
            hzodr->ifixb[j] = hertzodrdata->ifixb[j];
        }

        switch (hzodr->mode) {
        case R_MODE:
            hzodr->radius = hertzodrdata->radius;
            break;

        case ER_MODE:
            hzodr->Er = hertzodrdata->Er;
            break;

        case EIT_MODE:
            hzodr->Eit = hertzodrdata->Eit;
            hzodr->Er = calc_Er(hzodr->Eit, instdata->nu, instdata->Ei * 1e-9, instdata->nui);
            break;
        }

        hzodr++;
    }

    // initialize random generator and datafield
    rng = g_rand_new();

    //prepare datafield for noise
    hnoise = gwy_data_line_new_alike(h, FALSE);
    Fnoise = gwy_data_line_new_alike(h, FALSE);

    //prepare datafield for varied fields
    hnew = gwy_data_line_new_alike(h, FALSE);
    Fnew = gwy_data_line_new_alike(F, FALSE);

    hzodr = hertzodrmc->hertzodrmc;

    for (i = 0; i < Nmc; i++) {

        //create noise for h-F data
        generate_uncorrelated_noise(hnoise, hertzodrunc->instdata.sigma_h, rng);
        generate_uncorrelated_noise(Fnoise, hertzodrunc->instdata.sigma_F, rng);

        //add noise to original h-F data
        if (!sum_data_lines(hnew, h, hnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        if (!sum_data_lines(Fnew, F, Fnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        //remove negative depths from data
        nremove = filter_negative_data(hnew, Fnew);

        //run fit
        hertz_odr_fit(hnew, Fnew, hzodr, fddata, instdata, estimate);

        if (nremove) {
            gwy_data_line_resample(hnew, h->res, GWY_INTERPOLATION_NONE);
            gwy_data_line_resample(Fnew, F->res, GWY_INTERPOLATION_NONE);
        }

        hzodr++;
    }

    if (verbose) {
        printf("run Monte Carlo \n");
    }

    g_object_unref(h);
    g_object_unref(F);
    g_object_unref(hnew);
    g_object_unref(Fnew);
    g_object_unref(hnoise);
    g_object_unref(Fnoise);
    g_rand_free(rng);

    //calculate averages and stdevs
    // Er, Eit

    mc = hertzodrmc;

    for (i = 0; i < NPARAM_HZODR; i++) {
        mc->mcavg[i] = 0;
        mc->mcstd[i] = 0;
    }

    hzodr = hertzodrmc->hertzodrmc;
    nn = 0;

    for (i = 0; i < Nmc; i++) {
        //sanity check
        //        if (hzodr->fitinfo < 5) {
        mc->mcavg[0] += hzodr->Er;
        mc->mcavg[1] += hzodr->Eit;
        mc->mcavg[2] += hzodr->radius * 1e9; //TODO UNITS
        mc->mcstd[0] += hzodr->Er * hzodr->Er;
        mc->mcstd[1] += hzodr->Eit * hzodr->Eit;
        mc->mcstd[2] += hzodr->radius * 1e9 * hzodr->radius * 1e9; //TODO UNITS
        nn++;
        //       }

        hzodr++;
    }

    //use only reasonable data for average and standard deviation
    for (i = 0; i < NPARAM_HZODR; i++) {
        mc->mcavg[i] /= nn;
        mc->mcstd[i] /= nn;
        mc->mcstd[i] -= mc->mcavg[i] * mc->mcavg[i];
        mc->mcstd[i] = (mc->mcstd[i] < 0) ? 0 : mc->mcstd[i];
        mc->mcstd[i] = sqrt(mc->mcstd[i]);
    }

    // how many datasets were skipped
    hertzodrmc->skipmc = Nmc - nn;

    // Histograms
    //create data structure

    dd = hertzodrmc->mcdd;
    hzodr = hertzodrmc->hertzodrmc;

    for (i = 0; i < Nmc; i++) {
        dd[0][i] = hzodr->Er;
        dd[1][i] = hzodr->Eit;
        dd[2][i] = hzodr->radius * 1e9;
        hzodr++;
    }

    //prepare histogram structure

    for (i = 0; i < NPARAM_HZODR; i++) {
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

gboolean hertz_odr_unc_has_results(const HertzODRUncdata *hertzodrunc)
{
    return (hertzodrunc->ERc != NULL);
}

#endif

#ifndef NOFORTRAN

#include "tool_stiffness_unc.h"
#include "tool_stiffness.h"

#include "datatypes.h"
#include "fddata.h"
#include "fit-utils.h"
#include "file-utils.h"
#include "niget-common.h"
#include "mc-utils.h"

#include <math.h>
#include <string.h>


static const gchar *logfnmload = "fit.log.stiff.load.unc";
static const gchar *logfnmunload = "fit.log.stiff.unload.unc";

//static const gchar *mclogfnmload = "fit.log.load.stiff.mc";
//static const gchar *mclogfnmunload = "fit.log.unload.stiff.mc";

static gboolean nparam_to_string(gint i, gchar *str, gint nstr);
static void init_shifted_Stiffnessdata(const Stiffnessdata *src, Stiffnessdata *dest, gdouble dh);

void init_StiffnessUncdata(StiffnessUncdata *stiffnessunc, Instdata *instdata)
{
    stiffnessunc->instdata = *instdata;
    stiffnessunc->Nmc = N_MC;

    stiffnessunc->ukluh = 0;
    stiffnessunc->ukluF = 0;
    stiffnessunc->uqluh = 0;
    stiffnessunc->uqluF = 0;
    stiffnessunc->ukuuh = 0;
    stiffnessunc->ukuuF = 0;
    stiffnessunc->uquuh = 0;
    stiffnessunc->uquuF = 0;

    stiffnessunc->klc = NULL;
    stiffnessunc->qlc = NULL;
    stiffnessunc->kuc = NULL;
    stiffnessunc->quc = NULL;
}

static void init_shifted_Stiffnessdata(const Stiffnessdata *src, Stiffnessdata *dest, gdouble dh)
{
    init_Stiffnessdata(dest);

    // shift fitting range  in h-mode by dh
    // don't shift range for F-mode

    dest->loadrange[0] = src->loadrange[0] - dh;
    dest->loadrange[1] = src->loadrange[1] - dh;
    dest->loadfrom = dest->loadrange[0];
    dest->loadto = dest->loadrange[1];

    dest->unloadrange[0] = src->unloadrange[0] - dh;
    dest->unloadrange[1] = src->unloadrange[1] - dh;
    dest->unloadfrom = dest->unloadrange[0];
    dest->unloadto = dest->unloadrange[1];

    g_free(dest->logfnmload);
    g_free(dest->logfnmunload);
    dest->logfnmload = create_log_filename_path(logfnmload);
    dest->logfnmunload = create_log_filename_path(logfnmunload);

}

void stiffness_propagate_uncertainties_calc(const Stiffnessdata *stdata, StiffnessUncdata *unc, const FDdata *fd)
{
    gint istart, iend, ndata;
    gdouble ukl, uql;
    gdouble uklh, uqlh;
    gdouble uklF, uqlF;
    gdouble uku, uqu;
    gdouble ukuh, uquh;
    gdouble ukuF, uquF;
    gdouble uh, uF;

    gdouble beta[2], std_beta[2], covariance[2];
    /*always fit both parameters */
    gint ifixb[2] = {1, 1};
    gint np = 2;
    GwyDataLine *x, *y;
    gdouble ux, uy;

    gchar *unc_logfnmload, *unc_logfnmunload;
    gdouble SSres, R2, R2adj, chi2;

    //make copies to simplify code
    uh = unc->instdata.sigma_h;
    uF = unc->instdata.sigma_F;

    unc_logfnmload =  create_log_filename_path(logfnmload);
    unc_logfnmunload =  create_log_filename_path(logfnmunload);

    if (verbose) {
        printf(" Start propagation \n");
    }

    if (verbose) {
        g_print(" \n\nNo propagation of uncertainties\n");
    }

    /* prepare odr fit */

    //find appropriate data, which should be fitted, use same regime as for the main fit
    range_to_indices(stdata->loadrange[0], stdata->loadrange[1], fd->hload, FALSE, &istart, &iend, &ndata);

    x = gwy_data_line_part_extract(fd->hload, istart, ndata);
    y = gwy_data_line_part_extract(fd->Fload, istart, ndata);

    /* fit with only uh */
    ux = uh;
    uy = 0;
    beta[0] = stdata->kl;
    beta[1] = stdata->ql;
    odr_run_fit(x, y, np, beta, ifixb, unc_logfnmload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, LINEAR);

    uklh = std_beta[0];
    uqlh = std_beta[1];

    /*
      if (verbose)
      g_print(" fit with only uh:\n  gamma %g n %g h0 %g \n  u(n) %g u(h0) %g corr(n,h0) %g \n",
      beta[0], beta[1], beta[2], unh, uh0h, covnh0h / unh / uh0h);
    */

    /* fit with only uF */
    ux = 0;
    uy = uF;
    odr_run_fit(x, y, np, beta, ifixb, unc_logfnmload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, LINEAR);
    uklF = std_beta[0];
    uqlF = std_beta[1];

    /*
    if (verbose)
    g_print(" fit with only uF:\n  gamma %g n %g h0 %g \n  u(n) %g u(h0) %g corr(n,h0) %g \n",
    beta[0], beta[1], beta[2], unF, uh0F, covnh0F / unF / uh0F);
    */

    unc->ukluh = uklh;
    unc->ukluF = uklF;
    unc->uqluh = uqlh;
    unc->uqluF = uqlF;

    /* fit with both uh and uF */
    ux = uh;
    uy = uF;
    odr_run_fit(x, y, np, beta, ifixb, unc_logfnmload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, LINEAR);
    ukl = std_beta[0];
    uql = std_beta[1];

    /*
      if (verbose)
      g_print(" fit with both uh and uF:\n  gamma %g n %g h0 %g \n  u(n) %g u(h0) %g corr(n,h0) %g \n",
      beta[0], beta[1], beta[2], unF, uh0F, covnh0F / unF / uh0F);
    */

    unc->ukl = ukl;
    unc->uql = uql;

    /* prepare odr fit */

    //find appropriate data, which should be fitted, use same regime as for the main fit
    range_to_indices(stdata->unloadrange[0], stdata->unloadrange[1], fd->hunload, TRUE, &istart, &iend, &ndata);

    x = gwy_data_line_part_extract(fd->hunload, istart, ndata);
    y = gwy_data_line_part_extract(fd->Funload, istart, ndata);

    /* fit with only uh */
    ux = uh;
    uy = 0;
    beta[0] = stdata->ku;
    beta[1] = stdata->qu;
    odr_run_fit(x, y, np, beta, ifixb, unc_logfnmunload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, LINEAR);

    ukuh = std_beta[0];
    uquh = std_beta[1];

    /* if (verbose)
       g_print(" fit with only uh:\n  alpha %g m %g hp %g \n  u(m) %g u(hp) %g corr(m,hp) %g \n",
       beta[0], beta[1], beta[2], umh, uhph, covmhph / umh / uhph);
    */

    /* fit with only uF */
    ux = 0;
    uy = uF;
    odr_run_fit(x, y, np, beta, ifixb, unc_logfnmunload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, LINEAR);

    ukuF = std_beta[0];
    uquF = std_beta[1];

    /*
      if (verbose)
      g_print(" fit with only uF:\n  alpha %g m %g hp %g \n  u(m) %g u(hp) %g corr(m,hp) %g \n",
      beta[0], beta[1], beta[2], umF, uhpF, covmhpF / umF / uhpF);
    */

    unc->ukuuh = ukuh;
    unc->ukuuF = ukuF;
    unc->uquuh = uquh;
    unc->uquuF = uquF;

    /* fit with both */
    ux = uh;
    uy = uF;
    odr_run_fit(x, y, np, beta, ifixb, unc_logfnmunload, ux, uy, std_beta, covariance, &SSres, &R2, &R2adj, &chi2, LINEAR);

    uku = std_beta[0];
    uqu = std_beta[1];

    /*
      if (verbose)
      g_print(" fit with both:\n  alpha %g m %g hp %g \n  u(m) %g u(hp) %g corr(m,hp) %g \n",
      beta[0], beta[1], beta[2], umF, uhpF, covmhpF / umF / uhpF);
    */

    unc->uku = uku;
    unc->uqu = uqu;


    if (verbose) {
        g_print("end propagation %g %g %g %g  \n", unc-> uku, unc->uqu, uh, uF);
    }

    g_object_unref(x);
    g_object_unref(y);
    g_free(unc_logfnmload);
    g_free(unc_logfnmunload);
}

gchar *stiffness_uncertainty_export_data(const Stiffnessdata *stdata, const StiffnessUncdata *stunc, const FDdata *fddata, enum ExportFormat exportformat)
{
    gint j;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Parameters of two stiffness analysis (ODR fit) \n");
    g_string_append_printf(buf, "# u(depth): %g nm\n", stunc->instdata.sigma_h);
    g_string_append_printf(buf, "# u(load): %g nm\n", stunc->instdata.sigma_F);
    g_string_append_printf(buf, "#Fitted loading range: ");
    g_string_append_printf(buf, "[%g , %g ] nm \n", stdata->loadrange[0], stdata->loadrange[1]);
    g_string_append_printf(buf, "# %d datapoints used.\n", stdata->nfitloaddata);
    g_string_append_printf(buf, "#Fitted unloading range: ");
    g_string_append_printf(buf, "[%g , %g ] nm \n", stdata->unloadrange[0], stdata->unloadrange[1]);
    g_string_append_printf(buf, "# %d datapoints used.\n", stdata->nfitunloaddata);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Results of stiffness analysis (ODR fit) \n");
    g_string_append_printf(buf, "#kl: %g mN/nm\n", stdata->kl);
    g_string_append_printf(buf, "#ql: %g nm\n", stdata->ql);
    g_string_append_printf(buf, "#ku: %g mN/nm\n", stdata->ku);
    g_string_append_printf(buf, "#qu: %g nm\n", stdata->qu);
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Results of the Gaussian propagation of uncertainties \n");

    g_string_append_printf(buf, "#u(kl) = %g mN/nm\n", stdata->kl);
    g_string_append_printf(buf, "#u(ql) = %g nm\n", stdata->ql);
    g_string_append_printf(buf, "#u(ku) = %g mN/nm\n", stdata->ku);
    g_string_append_printf(buf, "#u(qu) = %g nm\n", stdata->qu);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of k_load\n");
    g_string_append_printf(buf, "#u(kl;uh) = %g mN/nm \n", stunc->ukluh);
    g_string_append_printf(buf, "#u(kl;uF) = %g mN/nm \n", stunc->ukluF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of q_load\n");
    g_string_append_printf(buf, "#u(ql;uh) = %g nm \n", stunc->uqluh);
    g_string_append_printf(buf, "#u(ql;uF) = %g nm \n", stunc->uqluF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of k_unload\n");
    g_string_append_printf(buf, "#u(ku;uh) = %g mN/nm \n", stunc->ukuuh);
    g_string_append_printf(buf, "#u(ku;uF) = %g mN/nm \n", stunc->ukuuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of q_unload\n");
    g_string_append_printf(buf, "#u(qu;uh) = %g nm \n", stunc->uquuh);
    g_string_append_printf(buf, "#u(qu;uF) = %g nm \n", stunc->uquuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Changes of k_load corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#kl/(mN/nm) \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (stunc->klc[j + CONTACTMAX] > -G_MAXINT) {
            g_string_append_printf(buf, "# %d   %g  \n", j, stunc->klc[j + CONTACTMAX]);
        }
        else {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
    }

    g_string_append_printf(buf, "#Changes of q_load corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#ql/nm \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (stunc->qlc[j + CONTACTMAX] > -G_MAXINT) {
            g_string_append_printf(buf, "# %d   %g  \n", j, stunc->qlc[j + CONTACTMAX]);
        }
        else {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
    }

    g_string_append_printf(buf, "#Changes of k_unload corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#ku/(mN/nm) \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (stunc->kuc[j + CONTACTMAX] > -G_MAXINT) {
            g_string_append_printf(buf, "# %d   %g  \n", j, stunc->kuc[j + CONTACTMAX]);
        }
        else {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
    }

    g_string_append_printf(buf, "#Changes of q_unload corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#qu/nm \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (stunc->quc[j + CONTACTMAX] > -G_MAXINT) {
            g_string_append_printf(buf, "# %d   %g  \n", j, stunc->quc[j + CONTACTMAX]);
        }
        else {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
    }

    g_string_append_printf(buf, "#\n");

    return g_string_free(buf, FALSE);
}

void stiffness_fit_shift_contact(const Stiffnessdata *stdata, const FDdata *fddata, const Instdata *instdata, gdouble *klc, gdouble *qlc, gdouble *kuc, gdouble *quc, gint j)
{
    Stiffnessdata stiffness;
    FDdata fdc;
    gdouble dh;
    gint istart, iend, ndata;
    GwyDataLine *x, *y;

    klc[j + CONTACTMAX] = -G_MAXINT;
    qlc[j + CONTACTMAX] = -G_MAXINT;
    kuc[j + CONTACTMAX] = -G_MAXINT;
    quc[j + CONTACTMAX] = -G_MAXINT;

    if (fddata->i_contact_load + j < 0) {
        return;
    }

    // create shifted FDdata
    dh = create_shifted_FDdata(fddata, &fdc, j);

    // create working copy of Stiffnessdata
    init_shifted_Stiffnessdata(stdata, &stiffness, dh);

    //check range
    if (stiffness.loadfrom == stiffness.loadto) {
        return;
    }

    //find appropriate data, which should be fitted, use same regime as for the main fit
    range_to_indices(stiffness.loadfrom, stiffness.loadto, fdc.hload, FALSE, &istart, &iend, &ndata);

    // if there are not enough data in the selected range, no fitting can be performed
    ndata = MIN(ndata, fdc.hload->res - istart);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    x = gwy_data_line_part_extract(fdc.hload, istart, ndata);
    y = gwy_data_line_part_extract(fdc.Fload, istart, ndata);

    // perform fit and get results
    stiffness_fit_load(x, y, &stiffness, &fdc, instdata, NULL);

    // clean up
    g_object_unref(x);
    g_object_unref(y);

    //check range
    if (stiffness.unloadfrom == stiffness.unloadto) {
        return;
    }

    //find appropriate data, which should be fitted, use same regime as for the main fit
    range_to_indices(stiffness.unloadfrom, stiffness.unloadto, fdc.hunload, TRUE, &istart, &iend, &ndata);

    // if there are not enough data in the selected range, no fitting can be performed
    ndata = MIN(ndata, fdc.hunload->res - istart);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    x = gwy_data_line_part_extract(fdc.hunload, istart, ndata);
    y = gwy_data_line_part_extract(fdc.Funload, istart, ndata);

    // perform fit and get results
    stiffness_fit_unload(x, y, &stiffness, &fdc, instdata, NULL);

    stiffness_combine_fit_results(&stiffness, &fdc, instdata);

    // clean up
    g_free(stiffness.logfnmload);  //TODO
    g_free(stiffness.logfnmunload);// TODO

    g_object_unref(x);
    g_object_unref(y);
    destroy_FDdata(&fdc);

    klc[j + CONTACTMAX] = stiffness.kl;
    qlc[j + CONTACTMAX] = stiffness.ql;
    kuc[j + CONTACTMAX] = stiffness.ku;
    quc[j + CONTACTMAX] = stiffness.qu;
}

void stiffness_uncertainty_montecarlo_run_calc(const Stiffnessdata *stiffnessdata, StiffnessUncdata *stiffnessunc, StiffnessMCdata *stiffnessmc,
        FDdata *fddata)
{
    gint istart, iend, ndata;
    gint i,  Nmc;
    Stiffnessdata *stiffness;
    StiffnessMCdata *mc;
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

    instdata = &(stiffnessunc->instdata);

    //not necessary to check range, couldn't get here in the first place

    //find appropriate loading data, which should be fitted, use same regime as for the main fit
    range_to_indices(stiffnessdata->loadfrom, stiffnessdata->loadto, fddata->hload, FALSE, &istart, &iend, &ndata);

    // if there are not enough data in the selected range, not fitting can be performed
    ndata = MIN(ndata, fddata->hload->res - istart);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    hl = gwy_data_line_part_extract(fddata->hload, istart, ndata);
    Fl = gwy_data_line_part_extract(fddata->Fload, istart, ndata);

    //find appropriate unloadingdata, which should be fitted, use same regime as for the main fit
    range_to_indices(stiffnessdata->unloadfrom, stiffnessdata->unloadto, fddata->hunload, TRUE, &istart, &iend, &ndata);

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
    estimatel[0] = stiffnessdata->kl;
    estimatel[1] = stiffnessdata->ql;
    estimateu[0] = stiffnessdata->ku;
    estimateu[1] = stiffnessdata->qu;

    Nmc = stiffnessmc->Nmc;
    stiffness = stiffnessmc->stiffnessmc;

    //prepare data
    for (i = 0; i < Nmc ; i++) {
        init_Stiffnessdata(stiffness);
        //change log filenames
        g_free(stiffness->logfnmload);
        g_free(stiffness->logfnmunload);
        // send to stdout, which can be redirected to /dev/null
        // otherwise too much IO time
        stiffness->logfnmload = NULL;
        stiffness->logfnmunload = NULL;
        //        stiffness->logfnmload = create_log_filename_path(mclogfnmload);
        //       stiffness->logfnmunload = create_log_filename_path(mclogfnmunload);
        stiffness++;
    }

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

    stiffness = stiffnessmc->stiffnessmc;

    for (i = 0; i < Nmc; i++) {

        //create noise for h-F data
        generate_uncorrelated_noise(hlnoise, instdata->sigma_h, rng);
        generate_uncorrelated_noise(Flnoise, instdata->sigma_F, rng);
        generate_uncorrelated_noise(hunoise, instdata->sigma_h, rng);
        generate_uncorrelated_noise(Funoise, instdata->sigma_F, rng);

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
        fddata->hmax = hmax_old + gaussian_random_number(rng) * instdata->sigma_h;
        fddata->Fmax = Fmax_old + gaussian_random_number(rng) * instdata->sigma_F;

        //run fit

        stiffness_fit_load(hlnew, Flnew, stiffness, fddata, instdata,  estimatel);
        stiffness_fit_unload(hunew, Funew, stiffness, fddata, instdata,  estimateu);
        stiffness_combine_fit_results(stiffness, fddata, instdata);

        stiffness++;
    }

    //reset hmax, Fmax to original value
    fddata->hmax = hmax_old;
    fddata->Fmax = Fmax_old;

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

    mc = stiffnessmc;

    for (i = 0; i < NPARAM_STIFFNESS; i++) {
        mc->mcavg[i] = 0;
        mc->mcstd[i] = 0;
    }

    stiffness = stiffnessmc->stiffnessmc;
    nn = 0;

    for (i = 0; i < Nmc; i++) {
        //check sanity
        mc->mcavg[0] += stiffness->kl;
        mc->mcavg[1] += stiffness->ql;
        mc->mcavg[2] += stiffness->ku;
        mc->mcavg[3] += stiffness->qu;
        mc->mcstd[0] += stiffness->kl * stiffness->kl;
        mc->mcstd[1] += stiffness->ql * stiffness->ql;
        mc->mcstd[2] += stiffness->ku * stiffness->ku;
        mc->mcstd[3] += stiffness->qu * stiffness->qu;
        nn++;

        stiffness++;
    }

    nn = Nmc;

    //use only reasonable data for average and standard deviation
    for (i = 0; i < NPARAM_STIFFNESS; i++) {
        mc->mcavg[i] /= nn;
        mc->mcstd[i] /= nn;
        mc->mcstd[i] -= mc->mcavg[i] * mc->mcavg[i];
        mc->mcstd[i] = (mc->mcstd[i] < 0) ? 0 : mc->mcstd[i];
        mc->mcstd[i] = sqrt(mc->mcstd[i]);
    }

    // how many datasets were skipped
    stiffnessmc->skipmc = Nmc - nn;

    // Histograms
    //create data structure

    dd = stiffnessmc->mcdd;
    stiffness = stiffnessmc->stiffnessmc;

    for (i = 0; i < Nmc; i++) {
        dd[0][i] = stiffness->kl;
        dd[1][i] = stiffness->ql;
        dd[2][i] = stiffness->ku;
        dd[3][i] = stiffness->qu;
        stiffness++;
    }

    //prepare histogram structure

    for (i = 0; i < NPARAM_STIFFNESS; i++) {
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

gchar *stiffness_mc_export_data(const Stiffnessdata *stiffnessdata, const StiffnessUncdata *stiffnessunc, const StiffnessMCdata *stiffnessmc,
                                const FDdata *fddata, enum ExportFormat exportformat)
{
    Stiffnessdata *stdata;
    gint i, j;
    gchar str[200];
    gboolean is;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Monte Carlo simulation of two stiffness analysis (ODR fit)\n");
    g_string_append_printf(buf, "#Input uncertainties \n");
    g_string_append_printf(buf, "#u(h) = %g nm \n", stiffnessunc->instdata.sigma_h);
    g_string_append_printf(buf, "#u(F) = %g mN \n", stiffnessunc->instdata.sigma_F);
    g_string_append_printf(buf, "#Number of iterations  = %d  \n", stiffnessmc->Nmc);
    g_string_append_printf(buf, "\n");

    g_string_append_printf(buf, "#Fitted loadrange: ");
    g_string_append_printf(buf, "[%g , %g ] nm \n", stiffnessdata->loadrange[0], stiffnessdata->loadrange[1]);
    g_string_append_printf(buf, "# %d datapoints used.\n", stiffnessdata->nfitloaddata);
    g_string_append_printf(buf, "\n");
    g_string_append_printf(buf, "#Fitted unloadrange: ");
    g_string_append_printf(buf, "[%g , %g ] nm \n", stiffnessdata->unloadrange[0], stiffnessdata->unloadrange[1]);
    g_string_append_printf(buf, "# %d datapoints used.\n", stiffnessdata->nfitunloaddata);
    g_string_append_printf(buf, "\n");

    if (stiffnessmc->skipmc) {
        g_string_append_printf(buf, "#Fit failed in %d cases, the corresponding data are NOT included in the  estimates of mean and uncertainty but ARE included in the histograms.", stiffnessmc->skipmc);
    }

    g_string_append_printf(buf, "\n");

    if (stiffnessmc->skipmc) {
        g_string_append_printf(buf, "\n\n");
    }

    g_string_append_printf(buf, "#Results \n");

    for (i = 0; i < NPARAM_STIFFNESS; i++) {
        is = nparam_to_string(i, str, sizeof(str));
        g_string_append_printf(buf, "#%s = (%g ± %g ) \n", str, stiffnessmc->mcavg[i], stiffnessmc->mcstd[i]);
    }

    g_string_append_printf(buf, "\n\n\n");

    // Histograms

    if (stiffnessmc->mcnstat > 0) {
        g_string_append_printf(buf, "#Histogram\n");

        for (i = 0; i < NPARAM_STIFFNESS; i++) {
            g_string_append_printf(buf, "#Histogram\n");
            is = nparam_to_string(i, str, sizeof(str));

            if (is) {
                g_string_append_printf(buf, "#%s    count/a.u.  \n", str);

                for (j = 0; j < stiffnessmc->mcnstat; j++) {
                    g_string_append_printf(buf, "%g %g \n", stiffnessmc->mcxhist[i][j], stiffnessmc->mcyhist[i][j]);
                }

                g_string_append_printf(buf, "\n \n");
            }
        }
    }

    g_string_append_printf(buf, "\n\n\n");

    //Full data
    g_string_append_printf(buf, "#Full data \n");
    g_string_append_printf(buf, "#");

    for (j = 0; j < NPARAM_STIFFNESS; j++) {
        is = nparam_to_string(j, str, sizeof(str)); /* _is_ is not used anywhere */
        g_string_append_printf(buf, "%s ", str);
    }

    g_string_append_printf(buf, "\n");

    stdata = stiffnessmc->stiffnessmc;

    for (i = 0; i < stiffnessmc->Nmc; i++) {
        g_string_append_printf(buf, SLOPE, stdata->kl);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, DEPTH, stdata->ql);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, SLOPE, stdata->ku);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, DEPTH, stdata->qu);
        g_string_append_printf(buf, "  ");
        stdata++;
    }

    return g_string_free(buf, FALSE);
}

void init_StiffnessMCdata(StiffnessMCdata *stiffnessmc, gint Nmc)
{
    gint i;

    //number of iterations
    stiffnessmc->Nmc = Nmc;

    if (verbose) {
        g_print("initialize for %d iterations \n", stiffnessmc->Nmc);
    }

    //how much had to be skipped
    stiffnessmc->skipmc = 0;

    //dynamically allocate
    stiffnessmc->stiffnessmc = (Stiffnessdata *)g_malloc(stiffnessmc->Nmc * sizeof(Stiffnessdata));
    stiffnessmc->mcavg = (gdouble *)g_malloc(NPARAM_STIFFNESS * sizeof(gdouble));
    stiffnessmc->mcstd = (gdouble *)g_malloc(NPARAM_STIFFNESS * sizeof(gdouble));
    stiffnessmc->mcdd = (gdouble **) g_malloc(NPARAM_STIFFNESS * sizeof(gdouble *));
    stiffnessmc->mcxhist = (gdouble **) g_malloc(NPARAM_STIFFNESS * sizeof(gdouble *));
    stiffnessmc->mcyhist = (gdouble **) g_malloc(NPARAM_STIFFNESS * sizeof(gdouble *));

    stiffnessmc->mcnstat = (gint)floor(3.49 * cbrt(stiffnessmc->Nmc) + 0.5);

    for (i = 0; i < NPARAM_STIFFNESS; i++) {
        stiffnessmc->mcdd[i] = (gdouble *) g_malloc(stiffnessmc->Nmc * sizeof(gdouble));
        stiffnessmc->mcxhist[i] = (gdouble *) g_malloc(stiffnessmc->mcnstat * sizeof(gdouble));
        stiffnessmc->mcyhist[i] = (gdouble *) g_malloc(stiffnessmc->mcnstat * sizeof(gdouble));
    }
}

void destroy_StiffnessMCdata(StiffnessMCdata *stiffnessmc)
{
    gint i;

    //free everything
    g_free(stiffnessmc->mcavg);
    g_free(stiffnessmc->mcstd);

    for (i = 0; i < NPARAM_STIFFNESS; i++) {
        g_free(stiffnessmc->mcdd[i]);
        g_free(stiffnessmc->mcxhist[i]);
        g_free(stiffnessmc->mcyhist[i]);
    }

    g_free(stiffnessmc->mcdd);
    g_free(stiffnessmc->mcxhist);
    g_free(stiffnessmc->mcyhist);

    stiffnessmc->mcnstat = 0;
    stiffnessmc->skipmc = 0;

    for (i = 0; i < stiffnessmc->Nmc; i++) {
        g_free(stiffnessmc->stiffnessmc[i].logfnmload);
        g_free(stiffnessmc->stiffnessmc[i].logfnmunload);
    }

    g_free(stiffnessmc->stiffnessmc);
}

static gboolean nparam_to_string(gint i, gchar *str, gint nstr)
{
    switch (i) {
    case 0:
        g_snprintf(str, nstr, "k_load/(mN/nm)");
        break;

    case 1:
        g_snprintf(str, nstr, "q_load/(mN/nm)");
        break;

    case 2:
        g_snprintf(str, nstr, "k_unload/(mN/nm)");
        break;

    case 3:
        g_snprintf(str, nstr, "q_unload/(mN/nm)");
        break;

    default:
        g_printerr("Unhandled number of results.\n");
        return FALSE;
    }

    return TRUE;
}

gboolean stiffness_unc_has_results(const StiffnessUncdata *stiffnessunc)
{
    return (stiffnessunc->klc != NULL && stiffnessunc->qlc != NULL && stiffnessunc->kuc != NULL && stiffnessunc->quc != NULL);
}

#endif

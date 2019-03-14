#ifndef NOFORTRAN

#include "tool_op_odr.h"

#include "datatypes.h"
#include "file-utils.h"
#include "fit-utils.h"
#include "niget-common.h"

#include <math.h>
#include <stdlib.h>

/* eliminated in favour of standard parameter passing */
/* #define RESCALE FALSE */

static const gchar *logfnm = "fit.log.op";

void init_OPODRdata(OPODRdata *opodrdata)
{
    opodrdata->beta = BETA_DEFAULT_OPODR;

    opodrdata->hp = 0;
    opodrdata->hc = 0;
    opodrdata->hr = 0;
    opodrdata->eps = 0;
    opodrdata->m = 0;
    opodrdata->S = 0;
    opodrdata->Aphc = 0;

    opodrdata->range[0] = 0;
    opodrdata->range[1] = 0;
    opodrdata->range_pct_Fmax[0] = 0;
    opodrdata->range_pct_Fmax[1] = 0;

    opodrdata->radial_corr = FALSE;
    opodrdata->radial_angle = RADIAL_ANGLE_DEFAULT_OPODR;

    opodrdata->has_fit = FALSE;
    opodrdata->from = 0;
    opodrdata->to = 0;
    opodrdata->to_pct_Fmax = 0;
    opodrdata->from_pct_Fmax = 0;
    opodrdata->xfit = NULL;
    opodrdata->yfit = NULL;

    opodrdata->Er = 0;
    opodrdata->Eit = 0;
    opodrdata->Hit = 0;

    opodrdata->Finput = FALSE;
    opodrdata->nfitdata = 0;

    /* quality of fit stuff */
    opodrdata->SSres = -1;
    opodrdata->R2 = -1;
    opodrdata->R2adj = -1;
    opodrdata->chi2 = -1;

    opodrdata->logfnm = create_log_filename_path(logfnm);
    opodrdata->infolog = 0;
}

gboolean op_odr_has_results(const OPODRdata *opodr)
{
    return opodr->has_fit;
}

void op_odr_remove_all_fit_results(OPODRdata *opodrdata)
{
    //reset all results associated with a fit to zero
    opodrdata->SSres = -1;
    opodrdata->R2 = -1;
    opodrdata->R2adj = -1;
    opodrdata->chi2 = -1;
    opodrdata->hp = 0;
    opodrdata->m = 0;
    opodrdata->eps = 0;
    opodrdata->hc = 0;
    opodrdata->S = 0;
    opodrdata->Aphc = 0;
    opodrdata->range[0] = 0;
    opodrdata->range[1] = 0;
    opodrdata->range_pct_Fmax[0] = 0;
    opodrdata->range_pct_Fmax[1] = 0;
    opodrdata->from = 0;
    opodrdata->to = 0;
    opodrdata->Hit = 0;
    opodrdata->Eit = 0;
    opodrdata->Er = 0;
    opodrdata->infolog = 0;
}

gchar* op_odr_export_data(const OPODRdata *opodrdata, const FDdata *fddata, const Instdata *instdata, const Area *area,
			  enum ExportFormat exportformat)
{
    gint i;
    gchar *areatext;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    //Input parameters
    g_string_append_printf(buf, "#Parameters of OP ODR analysis \n");
    g_string_append_printf(buf, "#beta: %g \n", opodrdata->beta);
    g_string_append_printf(buf, "#sigma_h: %g nm\n", instdata->sigma_h);
    g_string_append_printf(buf, "#sigma_F: %g mN\n", instdata->sigma_F);
    g_string_append_printf(buf, "#delta: %g \n", instdata->delta);
    g_string_append_printf(buf, "#nu: %g \n", instdata->nu);
    g_string_append_printf(buf, "#nui: %g \n", instdata->nui);
    g_string_append_printf(buf, "#Ei: %g GPa\n", instdata->Ei * 1e-9); //TODO units
    // Area
    areatext = write_area_file(area);
    g_string_append_printf(buf, "%s \n", areatext);
    g_free(areatext);

    //warning if extrapolation in area
    if (area->mode == AREA_DATA && opodrdata->hc > area->xmax) {
        g_string_append_printf(buf, "#    WARNING: Extrapolation in area calibration! \n");
    }

    if (opodrdata->radial_corr) {
        g_string_append_printf(buf, "#Applied radial correction using α: %.4g°\n", opodrdata->radial_angle);
    }
    else {
        g_string_append_printf(buf, "#Radial correction not applied\n");
    }

    g_string_append_printf(buf, "\n\n");

    g_string_append_printf(buf, "#Results of OP ODR analysis \n");

    g_string_append_printf(buf, "#hmax: %g nm\n", fddata->hmax);
    g_string_append_printf(buf, "#Fmax: %g mN\n", fddata->Fmax);

    g_string_append_printf(buf, "#hp: %g nm\n", opodrdata->hp);
    g_string_append_printf(buf, "#m: %g \n", opodrdata->m);
    g_string_append_printf(buf, "#α: %g mN.nm^{-%.3f} \n", opodrdata->alpha, opodrdata->m);

    g_string_append_printf(buf, "#fit info: %d \n", opodrdata->fitinfo);
    g_string_append_printf(buf, "#SSres: %g \n", opodrdata->SSres);
    g_string_append_printf(buf, "#R-squared: %g \n", opodrdata->R2);
    g_string_append_printf(buf, "#adjusted R-squared: %g \n", opodrdata->R2adj);
    g_string_append_printf(buf, "#chi-squared: %g \n", opodrdata->chi2);

    g_string_append_printf(buf, "#epsilon: %g \n", opodrdata->eps);
    g_string_append_printf(buf, "#hc: %g nm\n", opodrdata->hc);
    g_string_append_printf(buf, "#S: %g mN/nm\n", opodrdata->S);
    g_string_append_printf(buf, "#Aphc: %g nm^2\n", opodrdata->Aphc);
    g_string_append_printf(buf, "#Hit: %g MPa\n", opodrdata->Hit);
    g_string_append_printf(buf, "#Er: %g GPa\n", opodrdata->Er);
    g_string_append_printf(buf, "#Eit: %g GPa\n", opodrdata->Eit);

    if (opodrdata->has_fit) {
        g_string_append_printf(buf, "\n \n");
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

        g_string_append_printf(buf, "\n# S Fit \n");
        g_string_append_printf(buf, "# x/nm y/mN \n");

        for (i = 0; i < opodrdata->xfit->res; i++) {
            g_string_append_printf(buf, "%g %g \n", opodrdata->xfit->data[i], opodrdata->yfit->data[i]);
        }
    }

    return g_string_free(buf, FALSE);
}

void radial_correction_op(gdouble alpha, gdouble *H, gdouble *Er, gdouble *Eit, gdouble *Aphc, const Instdata *instdata)
{
    /* alpha: _equivalent_ cone angle in degrees */
    /* Berkovich: alpha = 70.3 */

    gdouble K;
    gdouble H0, Hprev, Hnew;
    gdouble Er0, Erprev, Ernew;
    gdouble Eit0, Eitprev, Eitnew;
    gdouble nu, nui, Ei;

    nu = instdata->nu;
    nui = instdata->nui;
    Ei = instdata->Ei * 1e-9; /* all calculations here in GPa */

    gint iter = 0;

    /*convert from degrees to radian */
    alpha *= M_PI / 180;
    K = (1 - 2 * nu) * (1 + nu) / 2 * sin(alpha);

    H0 = *H * 1e-3; /* Hit in MPa, but for radial correction better GPa */
    Eit0 = *Eit;
    Er0 = *Er;

    Hnew = H0;
    Ernew = Er0;
    Eitnew = Eit0;

    do {
        Hprev = Hnew;
        Erprev = Ernew;
        Eitprev = Eitnew;

        Hnew = H0 * sq(1 + K * Hprev / Eitprev);
        Ernew = Er0 * (1 + K * Hprev / Erprev);
        Eitnew = calc_Eit(Ernew, nu, Ei, nui);

        iter++;
    }
    while ((fabs(Hnew / Hprev - 1) > 0.5e-4 || fabs(Ernew / Erprev - 1) > 0.5e-4) && (iter < 1000));

    if (verbose)
        g_print("RADIAL CORRECTION with alpha = %g°:\n  converged after %d iterations \n 	hardness from %g to %g \n \
                reduced modulus from %g to %g \n elastic modulus from %g to %g \n",
                alpha /= M_PI / 180, iter, H0 * 1e3, Hnew * 1e3, Er0, Ernew, Eit0, Eitnew);

    *Eit = Eitnew;
    *H = Hnew * 1e3; /* back to MPa */
    *Er = Ernew;
    *Aphc /= sq(1 + K * Hnew / Eitnew);

    /*
      A0 = (S/Er0)^2
      A = A0 ( 1+K H/E)^2
      Er = S/sqrt(A)	=S/sqrt(A0) /(1+K H/E) = Er0 /(1+K H/E)
      H0 = F/A0
      H = F/A = F/A/ (1 + K H/E)^2 = H0 / (1+K H/E)^2
    */
}

gint op_odr_fit(GwyDataLine *x, GwyDataLine *y, OPODRdata *opodrdata, const FDdata *fddata, const Instdata *instdata, const Area *area, const gdouble estimate[3])
{
    gdouble beta[3], std_beta[3], covariance[3];
    gint info;
    gint ifix[3];
    gint np;
    gdouble ux, uy;
    gint istart, iend, ndata = 0;
    gdouble r1 = 0.05, r2 = 0.1, rstep = 0.05;
    gdouble xmin;

    //fit the data with the function beta[0]*(h-beta[1])^beta[2]
    // beta [0] --- alpha, gamma
    // beta [1] --- hp, h0
    // beta [2] --- m, n
    if (estimate) {
        beta[0] = estimate[0];
        beta[1] = estimate[1];
        beta[2] = estimate[2];
    }
    else {
        beta[2] = 1.5;

        while (ndata == 0 && r2 < 1.0) {
            range_to_indices(r1 * fddata->Fmax, r2 * fddata->Fmax, fddata->Funload, TRUE, &istart, &iend, &ndata);
            r1 += rstep;
            r2 += rstep;
        }

        if (ndata > 0) {
            beta[1] = fddata->hunload->data[iend];
        }
        else {
            beta[1] = fddata->hmax * 0.05;
        }

        /* beta[1] (hp estimate) should be smaller than the smallest x */
        /* beta[1] equal to minimum does not work very well either */
        xmin = gwy_data_line_get_min(x);
        beta[1] = (beta[1] < xmin ? beta[1] : xmin * 0.9);

        if (verbose) {
            g_print(" beta[1] %g xmin %g \n", beta[1], gwy_data_line_get_min(x));
        }

        beta[0] = fddata->Fmax / pow(fddata->hmax - beta[1], beta[2]);
    }

    // TODO keep some parameters fixed 0:fixed 1:unfixed
    ifix[0] = 1;
    ifix[1] = 1;
    ifix[2] = 1;

    np = 3;

    // uncertainties given by sigma_h, sigma_F
    ux = instdata->sigma_h;
    uy = instdata->sigma_F;

    if (verbose) {
        g_print("beta estimate %g %g %g  \n", beta[0], beta[1], beta[2]);
    }

    info = odr_run_fit(x, y, np, beta, ifix, opodrdata->logfnm, ux, uy, std_beta, covariance,
                       &(opodrdata->SSres), &(opodrdata->R2), &(opodrdata->R2adj), &(opodrdata->chi2), POWER_3PARAM);

    if (verbose) {
        g_print("beta result %g %g %g  \n\n", beta[0], beta[1], beta[2]);
    }

    // calculate results
    opodrdata->alpha = beta[0];
    opodrdata->hp = beta[1];
    opodrdata->m = beta[2];

    opodrdata->eps = epsilon(opodrdata->m);
    opodrdata->S = calc_S(fddata->hmax, fddata->Fmax, opodrdata->hp, opodrdata->m);
    opodrdata->hc = calc_hc(fddata->hmax, fddata->Fmax, opodrdata->S, opodrdata->eps);
    opodrdata->Aphc = eval_area(area, opodrdata->hc);

    // TODO units
    // calculate results
    opodrdata->Hit = calc_Hit(fddata->Fmax, opodrdata->Aphc) * 1e9;
    opodrdata->Er = calc_Er_OP(opodrdata->Aphc, opodrdata->S, opodrdata->beta) * 1e15;
    opodrdata->Eit = calc_Eit(opodrdata->Er, instdata->nu, instdata->Ei, instdata->nui);
    opodrdata->Eit *= 1e-9;
    opodrdata->Er *= 1e-9;

    //clean up
    opodrdata->range[0] = opodrdata->from;
    opodrdata->range[1] = opodrdata->to;
    opodrdata->has_fit = TRUE;
    opodrdata->range_pct_Fmax[0] = opodrdata->from_pct_Fmax;
    opodrdata->range_pct_Fmax[1] = opodrdata->to_pct_Fmax;

    opodrdata->fitinfo = info;

    return info;
}

#endif

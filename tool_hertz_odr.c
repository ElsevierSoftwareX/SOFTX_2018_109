#ifndef NOFORTRAN

#include "tool_hertz_odr.h"
#include "tool_hertz_odr_unc.h"

#include "datatypes.h"
#include "fit-utils.h"
#include "niget-common.h"

#include <math.h>

#include <stdlib.h>


#define INTERCEPT FALSE

static const gchar *logfnm = "fit.log.hz";

void init_HertzODRdata(HertzODRdata *hzodrdata)
{
    hzodrdata->mode = R_MODE;
    hzodrdata->radius = 100e-9;
    hzodrdata->Er = 0;
    hzodrdata->Eit = 0;

    hzodrdata->has_fit = FALSE;
    hzodrdata->from = 0;
    hzodrdata->to = 0;
    hzodrdata->xfit = NULL;
    hzodrdata->yfit = NULL;

    hzodrdata->range[0] = 0;
    hzodrdata->range[1] = 0;
    hzodrdata->nfitdata = 0;

    /*initially parameters are not fixed */
    hzodrdata->ifixb[0] = 1;
    hzodrdata->ifixb[1] = 1;
    hzodrdata->ifixb[2] = 1;

    /* quality of fit stuff */
    hzodrdata->SSres = -1;
    hzodrdata->R2 = -1;
    hzodrdata->R2adj = -1;
    hzodrdata->chi2 = -1;

    hzodrdata->h0 = 0;
    hzodrdata->n = 0;
    hzodrdata->gamma = 0;
    hzodrdata->logfnm = create_log_filename_path(logfnm);
    hzodrdata->infolog = 0;

    hzodrdata->radial_corr = FALSE;
}

gboolean hertz_odr_has_results(const HertzODRdata *hzodr)
{
    return hzodr->has_fit;
}

void hertz_odr_remove_all_fit_results(HertzODRdata *hzodrdata)
{
    //reset all results associated with a fit to zero
    /*hzodrdata->range[0]=0;
    hzodrdata->range[1]=0;
    hzodrdata->from=0;
    hzodrdata->to=0;*/
    hzodrdata->SSres = -1;
    hzodrdata->R2 = -1;
    hzodrdata->R2adj = -1;
    hzodrdata->chi2 = -1;
    hzodrdata->gamma = 0;
    hzodrdata->h0 = 0;
    hzodrdata->n = 0;
    hzodrdata->Eit = 0;
    hzodrdata->Er = 0;
    hzodrdata->nfitdata = 0;
}

gchar *hertz_odr_export_data(const HertzODRdata *hzodrdata, const FDdata *fddata, const Instdata *instdata, enum ExportFormat exportformat)
{
    gint i;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    //Input parameters
    g_string_append_printf(buf, "#Parameters of HertzODR analysis \n");
    g_string_append_printf(buf, "#sigma_h: %g \n", instdata->sigma_h);
    g_string_append_printf(buf, "#sigma_F: %g \n", instdata->sigma_F);
    g_string_append_printf(buf, "#delta: %g \n", instdata->delta);

    switch (hzodrdata->mode) {
    case R_MODE:
        g_string_append_printf(buf, "#nu: %g \n", instdata->nu);
        g_string_append_printf(buf, "#nui: %g \n", instdata->nui);
        g_string_append_printf(buf, "#Ei: %g GPa\n", instdata->Ei * 1e-9); //TODO units
        g_string_append_printf(buf, "#tip radius: %g nm\n", hzodrdata->radius * 1e9);
        break;

    case ER_MODE:
        g_string_append_printf(buf, "#Er: %g GPa\n", hzodrdata->Er);
        break;

    case EIT_MODE:
        g_string_append_printf(buf, "#nu: %g \n", instdata->nu);
        g_string_append_printf(buf, "#nui: %g \n", instdata->nui);
        g_string_append_printf(buf, "#Ei: %g GPa\n", instdata->Ei * 1e-9); //TODO units
        g_string_append_printf(buf, "#Eit: %g GPa\n", hzodrdata->Eit);
        break;

    default:
        g_printerr("Unhandled mode type \n");
        break;
    }

    if (hzodrdata->radial_corr) {
        g_string_append_printf(buf, "#Fitted using radial correction ");
        g_string_append_printf(buf, "#b: %g nm^{-1/2} \n", hzodrdata->n);
    }
    else {
        g_string_append_printf(buf, "\n");
    }

    g_string_append_printf(buf, "\n");
    g_string_append_printf(buf, "\n");

    if (hzodrdata->has_fit) {
        g_string_append_printf(buf, "#Results of HertzODR analysis \n");

        if (hzodrdata->radial_corr) {
            g_string_append_printf(buf, "#γ: %g mN.nm^{-%g} \n", hzodrdata->gamma, N_FIX_HZODR);
        }
        else {
            g_string_append_printf(buf, "#γ: %g mN.nm^{-%g} \n", hzodrdata->gamma, hzodrdata->n);
        }

        g_string_append_printf(buf, "#h_0: %g nm \n", hzodrdata->h0);

        if (hzodrdata->ifixb[2]) {
            g_string_append_printf(buf, "#n: %g \n", hzodrdata->n);
        }
        else {
            g_string_append_printf(buf, "#n: %g fixed \n", N_FIX_HZODR);
        }

        g_string_append_printf(buf, "#fit info: %d \n", hzodrdata->fitinfo);

        g_string_append_printf(buf, "#SSres: %g \n", hzodrdata->SSres);
        g_string_append_printf(buf, "#R-squared: %g \n", hzodrdata->R2);
        g_string_append_printf(buf, "#adjusted R-squared: %g \n", hzodrdata->R2adj);
        g_string_append_printf(buf, "#chi-squared: %g \n", hzodrdata->chi2);

        switch (hzodrdata->mode) {
        case R_MODE:
            g_string_append_printf(buf, "#Er: %g GPa\n", hzodrdata->Er);
            g_string_append_printf(buf, "#Eit: %g GPa\n", hzodrdata->Eit);
            break;

        case ER_MODE:
        case EIT_MODE:
            g_string_append_printf(buf, "#tip radius: %g nm\n", hzodrdata->radius * 1e9);
            break;

        default:
            g_printerr("Unhandled mode type \n");
            break;
        }

        g_string_append_printf(buf, " \n");
        g_string_append_printf(buf, "#Fitted range: [%g , %g ] nm \n", hzodrdata->range[0], hzodrdata->range[1]);
        g_string_append_printf(buf, "# %d datapoints used.\n", hzodrdata->nfitdata);

        g_string_append_printf(buf, "\n \n ");
        g_string_append_printf(buf, "\n# S Fit \n");
        g_string_append_printf(buf, "# x/nm y/mN \n");

        for (i = 0; i < hzodrdata->xfit->res; i++) {
            g_string_append_printf(buf, "%g %g \n", hzodrdata->xfit->data[i], hzodrdata->yfit->data[i]);
        }
    }

    return g_string_free(buf, FALSE);
}

void radial_correction_hertz(gint mode, gdouble *Er, gdouble *Eit, gdouble *radius, gdouble to, const Instdata *instdata) /* z mode udelat enum */
{
    double gamma;
    double Ei;

    double Er0, Ernew;
    double Eit0, Eitnew;
    double radius0, radiusnew;

    double nu, nui;

    Er0 = *Er;
    /* Eit0 = *Eit; */
    radius0 = *radius * 1e9;

    nu = instdata->nu;
    nui = instdata->nui;
    Ei = instdata->Ei * 1e-9; /* convert to GPa */

    Eit0 = calc_Eit(Er0, nu, instdata->Ei, nui);

    // TODO iterate
    gamma = 1 + 2 / (3 * M_PI) * Ei / Eit0 * (1 - 2 * nu) * (1 + nu) / (Ei / Eit0 * (1 - sq(nu)) + (1 - sq(nui))) * sqrt(to / radius0);
    /* Hay-Wolff (28); a = sqrt(h*R) */
    printf(" gamma %g \n", gamma);

    switch (mode) {
    case R_MODE:
        Ernew = Er0 / gamma;
        Eitnew = calc_Eit(Ernew, nu, Ei, nui);

        if (verbose) {
            g_print(" Er, Eit before radial correction %g, %g \n", Er0, Eit0);
            g_print(" Er, Eit after radial correction %g, %g \n", Ernew, Eitnew);
        }

        *Er = Ernew;
        *Eit = Eitnew;
        break;

    case EIT_MODE:
        break;

    case ER_MODE:
        radiusnew = radius0 / (sq(gamma));

        if (verbose) {
            g_print(" R before radial correction %g \n", radius0);
            g_print(" R after radial correction %g \n", radiusnew);
        }

        *radius = radiusnew * 1e-9;
        break;
    }
}

gint hertz_odr_fit(GwyDataLine *x, GwyDataLine *y, HertzODRdata *hzodrdata, const FDdata *fddata, const Instdata *instdata, const gdouble estimate[3])
{
    gdouble beta[3], std_beta[3], covariance[3];
    gint info;
    gdouble ux, uy;
    gint np;

    np = 3;

    if (hzodrdata->mode != R_MODE && hzodrdata->radial_corr) {
        g_printerr(" radial correction currently working only for R-mode \n");
    }

    //fit the data with the function beta[0]*(h-beta[1])^beta[2]
    // beta [0] --- alpha, gamma
    // beta [1] --- hp, h0
    // beta [2] --- m, n
    // TODO units
    // calculate results according to mode

    if (estimate) {
        beta[0] = estimate[0];
        beta[1] = estimate[1];
        beta[2] = estimate[2];
    }
    else {
        beta[2] = N_FIX_HZODR;
        beta[0] = fddata->Fmax / pow(fddata->hmax, beta[2]);
        beta[1] = 0.9 * gwy_data_line_get_min(x);

        if (verbose) {
            g_print(" xmin %g  0.9*xmin %g beta[1] %g \n", gwy_data_line_get_min(x), 0.9 * gwy_data_line_get_min(x), beta[1]);
        }

        beta[1] = (beta[1] > 0 ? beta[1] : 1e-10);
    }

    if (!hzodrdata->ifixb[2]) {
        beta[2] = N_FIX_HZODR;
    }

    // uncertainties given by sigma_h, sigma_F
    ux = instdata->sigma_h;
    uy = instdata->sigma_F;

    if (hzodrdata->radial_corr) {

        if (verbose) {
            g_print("use radial correction  with b\n");
        }

        beta[2] = 2. / 3. / M_PI * (1 - 2 * instdata->nu) / (1 - instdata->nu) / sqrt(hzodrdata->radius * 1e9); /* must be in nanometers*/

        if (verbose) {
            g_print(" %g \n", beta[2]);
            g_print(" nu %g  radius %g \n", instdata->nu, hzodrdata->radius);
        }

        hzodrdata->ifixb[2] = 0;

        if (verbose) {
            g_print("ifixb %d %d %d \n", hzodrdata->ifixb[0], hzodrdata->ifixb[1], hzodrdata->ifixb[2]);
        }

        info = odr_run_fit(x, y, np, beta, hzodrdata->ifixb, hzodrdata->logfnm, ux, uy, std_beta, covariance,
                           &(hzodrdata->SSres), &(hzodrdata->R2), &(hzodrdata->R2adj), &(hzodrdata->chi2), HERTZ_RADIAL);

        if (verbose) {
            g_print(" info %d beta %g %g %g \n", info, beta[0], beta[1], beta[2]);
        }
    }
    else {
        info = odr_run_fit(x, y, np, beta, hzodrdata->ifixb, hzodrdata->logfnm, ux, uy, std_beta, covariance,
                           &(hzodrdata->SSres), &(hzodrdata->R2), &(hzodrdata->R2adj), &(hzodrdata->chi2), POWER_3PARAM);
    }

    // calculate results
    hzodrdata->gamma = beta[0];
    hzodrdata->h0 = beta[1];
    hzodrdata->n = beta[2];

    //	printf(" hzodrdata->mode %d\n", hzodrdata->mode);
    switch (hzodrdata->mode) {
    case R_MODE:
        hzodrdata->Er = calc_Er_Hertz(hzodrdata->gamma, hzodrdata->radius * 1e9) * 1e15;
        hzodrdata->Eit = calc_Eit(hzodrdata->Er, instdata->nu, instdata->Ei, instdata->nui) * 1e-9;
        hzodrdata->Er *= 1e-9;
        break;

    case EIT_MODE:
        hzodrdata->Er = calc_Er(hzodrdata->Eit, instdata->nu, instdata->Ei * 1e-9, instdata->nui);

    /* RS - missing break? */

    case ER_MODE:
        hzodrdata->radius =  9.0 / 16.0 * hzodrdata->gamma * hzodrdata->gamma / hzodrdata->Er / hzodrdata->Er * 1e3;
        break;
    }

    //clean up
    hzodrdata->range[0] = hzodrdata->from;
    hzodrdata->range[1] = hzodrdata->to;
    hzodrdata->has_fit = TRUE;
    hzodrdata->fitinfo = info;

    return info;
}

#endif

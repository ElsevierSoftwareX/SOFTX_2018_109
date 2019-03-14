#ifndef NOFORTRAN

#include "tool_twoslopes.h"
#include "tool_twoslopes_unc.h"

#include "datatypes.h"
#include "file-utils.h"
#include "fit-utils.h"
#include "niget-common.h"

#include <math.h>
#include <stdlib.h>

static const gchar *logfnmload = "fit.log.slopes.load";
static const gchar *logfnmunload = "fit.log.slopes.unload";

void init_Slopesdata(Slopesdata *sldata)
{
    sldata->beta = BETA_DEFAULT_SLOPES;
    sldata->gamma = 0;
    sldata->h0 = 0;
    sldata->n = 0;

    sldata->alpha = 0;
    sldata->hp = 0;
    sldata->m = 0;

    sldata->eps = 0;
    sldata->Sload = 0;
    sldata->Sunload = 0;

    sldata->Aphc = 0;

    sldata->has_fit_load = FALSE;
    sldata->loadto = 0;
    sldata->loadfrom = 0;
    sldata->loadto_pct_Fmax = 0;
    sldata->loadfrom_pct_Fmax = 0;
    sldata->loadrange[0] = 0;
    sldata->loadrange[1] = 0;
    sldata->loadrange_pct_Fmax[0] = 0;
    sldata->loadrange_pct_Fmax[1] = 0;
    sldata->loadifixb[0] = 1;
    sldata->loadifixb[1] = 1;
    sldata->loadifixb[2] = 1;

    sldata->has_fit_unload = FALSE;
    sldata->unloadto = 0;
    sldata->unloadfrom = 0;
    sldata->unloadto_pct_Fmax = 0;
    sldata->unloadfrom_pct_Fmax = 0;
    sldata->unloadrange[0] = 0;
    sldata->unloadrange[1] = 0;
    sldata->unloadrange_pct_Fmax[0] = 0;
    sldata->unloadrange_pct_Fmax[1] = 0;
    sldata->unloadifixb[0] = 1;
    sldata->unloadifixb[1] = 1;
    sldata->unloadifixb[2] = 1;

    sldata->xloadfit = NULL;
    sldata->yloadfit = NULL;
    sldata->xunloadfit = NULL;
    sldata->yunloadfit = NULL;

    sldata->Er = 0;
    sldata->Eit = 0;
    sldata->Hit = 0;

    sldata->Finputload = FALSE;
    sldata->nfitloaddata = 0;
    sldata->Finputunload = FALSE;
    sldata->nfitunloaddata = 0;

    /* quality of fit stuff */
    sldata->loadSSres = -1;
    sldata->loadR2 = -1;
    sldata->loadR2adj = -1;
    sldata->loadchi2 = -1;
    sldata->unloadSSres = -1;
    sldata->unloadR2 = -1;
    sldata->unloadR2adj = -1;
    sldata->unloadchi2 = -1;

    // TODO build logfilename, modify also in unc and mc
    sldata->logfnmload = create_log_filename_path(logfnmload);
    sldata->infologload = 0;

    sldata->logfnmunload = create_log_filename_path(logfnmunload);
    sldata->infologunload = 0;
}

gboolean slopes_has_results(const Slopesdata *slopes)
{
    return (slopes->has_fit_load && slopes->has_fit_unload);
}

void slopes_remove_all_fit_results_unload(Slopesdata *sldata)
{
    //reset all results associated with the powerlaw fit of the unloading curve to zero
    sldata->alpha = 0;
    sldata->hp = 0;
    sldata->m = 0;
    sldata->unloadSSres = -1;
    sldata->unloadR2 = -1;
    sldata->unloadR2adj = -1;
    sldata->unloadchi2 = -1;

    sldata->eps = 0;
    sldata->Sunload = 0;

    sldata->Aphc = 0;
    sldata->Er = 0;
    sldata->Eit = 0;
    sldata->Hit = 0;

    /*
    sldata->unloadto=0;
    sldata->unloadfrom=0;
    sldata->unloadrange[0]=0;
    sldata->unloadrange[1]=0;
    sldata->unloadrange_pct_Fmax[0]=0;
    sldata->unloadrange_pct_Fmax[1]=0;
    */

    sldata->Finputunload = FALSE;
    sldata->has_fit_unload = FALSE;
    sldata->nfitunloaddata = 0;

    sldata->infologunload = 0;
}

void slopes_remove_all_fit_results_load(Slopesdata *sldata)
{
    //reset all results associated with the powerlaw fit of the loading curve to zero
    sldata->gamma = 0;
    sldata->h0 = 0;
    sldata->n = 0;
    sldata->loadSSres = -1;
    sldata->loadR2 = -1;
    sldata->loadR2adj = -1;
    sldata->loadchi2 = -1;

    sldata->Sload = 0;

    sldata->Aphc = 0;
    sldata->Er = 0;
    sldata->Eit = 0;
    sldata->Hit = 0;

    /*
    sldata->loadto=0;
    sldata->loadfrom=0;
    sldata->loadrange[0]=0;
    sldata->loadrange[1]=0;
    sldata->loadrange_pct_Fmax[0]=0;
    sldata->loadrange_pct_Fmax[1]=0;
    */

    sldata->Finputload = FALSE;
    sldata->has_fit_load = FALSE;
    sldata->nfitloaddata = 0;

    sldata->infologload = 0;
}

gchar* slopes_export_data(const Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata, enum ExportFormat exportformat)
{
    gint i;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    //Input parameters
    g_string_append_printf(buf, "#Parameters of Slopes  analysis \n");
    g_string_append_printf(buf, "#beta: %g \n", sldata->beta);
    g_string_append_printf(buf, "#sigma_h: %g \n", instdata->sigma_h);
    g_string_append_printf(buf, "#sigma_F: %g \n", instdata->sigma_F);
    g_string_append_printf(buf, "#delta: %g \n", instdata->delta);
    g_string_append_printf(buf, "#nu: %g \n", instdata->nu);
    g_string_append_printf(buf, "#nui: %g \n", instdata->nui);
    g_string_append_printf(buf, "#Ei: %g GPa\n", instdata->Ei * 1e-9); //TODO units

    g_string_append_printf(buf, "\n\n");
    g_string_append_printf(buf, "#Results of Slopes  analysis \n");
    g_string_append_printf(buf, "#hmax: %g nm\n", fddata->hmax);
    g_string_append_printf(buf, "#Fmax: %g mN\n", fddata->Fmax);

    g_string_append_printf(buf, "#Fit of loading curve     gamma (h-h0)^n\n");
    g_string_append_printf(buf, "#gamma: %g mN.nm^{-%g}\n", sldata->gamma, sldata->n);
    g_string_append_printf(buf, "#h0: %g nm\n", sldata->h0);

    if (sldata->loadifixb[2]) {
        g_string_append_printf(buf, "#n: %g \n", sldata->n);
    }
    else {
        g_string_append_printf(buf, "#n: %g fixed \n", N_FIX_SLOPES);
    }

    g_string_append_printf(buf, "#fit info: %d \n", sldata->fitinfoload);
    g_string_append_printf(buf, "#SSres: %g \n", sldata->loadSSres);
    g_string_append_printf(buf, "#R-squared: %g \n", sldata->loadR2);
    g_string_append_printf(buf, "#adjusted R-squared: %g \n", sldata->loadR2adj);
    g_string_append_printf(buf, "#chi-squared: %g \n", sldata->loadchi2);

    g_string_append_printf(buf, "#Fit of unloading curve     alpha (h-hp)^m\n");
    g_string_append_printf(buf, "#alpha: %g mN.nm^{-%g}\n", sldata->alpha, sldata->m);
    g_string_append_printf(buf, "#hp: %g nm\n", sldata->hp);
    g_string_append_printf(buf, "#m: %g nm\n", sldata->m);

    g_string_append_printf(buf, "#fit info: %d \n", sldata->fitinfounload);
    g_string_append_printf(buf, "#SSres: %g \n", sldata->unloadSSres);
    g_string_append_printf(buf, "#R-squared: %g \n", sldata->unloadR2);
    g_string_append_printf(buf, "#adjusted R-squared: %g \n", sldata->unloadR2adj);
    g_string_append_printf(buf, "#chi-squared: %g \n", sldata->unloadchi2);

    g_string_append_printf(buf, "#epsilon: %g \n", sldata->eps);
    g_string_append_printf(buf, "#Sload: %g mN/nm\n", sldata->Sload);
    g_string_append_printf(buf, "#Sunload: %g mN/nm\n", sldata->Sunload);
    g_string_append_printf(buf, "#Aphc: %g nm^2\n", sldata->Aphc);
    g_string_append_printf(buf, "#Hit: %g MPa\n", sldata->Hit);
    g_string_append_printf(buf, "#Er: %g GPa\n", sldata->Er);
    g_string_append_printf(buf, "#Eit: %g GPa\n", sldata->Eit);

    if (sldata->has_fit_load) {
        g_string_append_printf(buf, "\n \n");
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

        g_string_append_printf(buf, "\n# Loading Fit \n");
        g_string_append_printf(buf, "# x/nm y/mN \n");

        for (i = 0; i < sldata->xloadfit->res; i++) {
            g_string_append_printf(buf, "%g %g \n", sldata->xloadfit->data[i], sldata->yloadfit->data[i]);
        }
    }

    if (sldata->has_fit_unload) {
        g_string_append_printf(buf, "\n \n");
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

        g_string_append_printf(buf, "\n# Unloading Fit \n");
        g_string_append_printf(buf, "# x/nm y/mN \n");

        for (i = 0; i < sldata->xunloadfit->res; i++) {
            g_string_append_printf(buf, "%g %g \n", sldata->xunloadfit->data[i], sldata->yunloadfit->data[i]);
        }
    }

    return g_string_free(buf, FALSE);
}

void slopes_combine_fit_results(Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata)
{
    gdouble C0 = 24.5, Fmax ;
    gdouble K;
    // TODO 24.5

    Fmax = fddata->Fmax;

    K = (2 * sldata->Sunload - sldata->Sload * sldata->eps * sldata->beta) / sldata->Sunload / sldata->Sload;
    //printf(" K %g K(beta=1) %g  \n", K, (2*sldata->Sunload - sldata->Sload*sldata->eps*BETA)/sldata->Sunload/sldata->Sload);

    //TODO Units
    sldata->Aphc = C0 * Fmax * Fmax * K * K;
    sldata->Hit = 1 / C0 / Fmax / K / K * 1e9;

    sldata->Er = sqrt(M_PI / C0) / 2 / sldata->beta / Fmax / K * sldata->Sunload * 1e6;

    sldata->Eit = (1 - instdata->nu * instdata->nu) / (1 / sldata->Er - (1 - instdata->nui * instdata->nui) / (instdata->Ei * 1e-9));

    /*printf("BETA %g \n", BETA);
    printf("Hit %g \n", sldata->Hit);
    printf("Er %g \n", sldata->Er);
    printf("Eit %g \n", sldata->Eit);
    */
}

gint slopes_fit_load(GwyDataLine *x, GwyDataLine *y, Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata, const gdouble estimate[3])
{
    gdouble beta[3], std_beta[3], covariance[3];
    gint info;
    gdouble ux, uy;
    gint np;

    np = 3;


    //fit the data with the function beta[0]*(h-beta[1])^beta[2]
    if (estimate) {
        beta[0] = estimate[0];
        beta[1] = estimate[1];
        beta[2] = estimate[2];
    }
    else {
        beta[2] = N_FIX_SLOPES;
        beta[0] = fddata->Fmax / pow(fddata->hmax, beta[2]);
        beta[1] = 1;
    }

    if (!sldata->loadifixb[2]) {
        beta[2] = N_FIX_SLOPES;
    }

    // uncertainties given by sigma_h, sigma_F
    ux = instdata->sigma_h;
    uy = instdata->sigma_F;


    info = odr_run_fit(x, y, np, beta, sldata->loadifixb, sldata->logfnmload, ux, uy, std_beta, covariance,
                       &(sldata->loadSSres), &(sldata->loadR2), &(sldata->loadR2adj), &(sldata->loadchi2), POWER_3PARAM);
    //	printf("load after beta %g %g %g  \n\n", beta[0], beta[1], beta[2]);

    // calculate results
    sldata->gamma = beta[0];
    sldata->h0 = beta[1];
    sldata->n = beta[2];

    sldata->Sload = sldata->n * fddata->Fmax / (fddata->hmax - sldata->h0);

    if (sldata->has_fit_unload) {
        slopes_combine_fit_results(sldata, fddata, instdata);
    }

    //clean up
    sldata->loadrange[0] = sldata->loadfrom;
    sldata->loadrange[1] = sldata->loadto;
    sldata->has_fit_load = TRUE;
    sldata->loadrange_pct_Fmax[0] = sldata->loadfrom_pct_Fmax;
    sldata->loadrange_pct_Fmax[1] = sldata->loadto_pct_Fmax;

    sldata->fitinfoload = info;

    return info;
}

gint slopes_fit_unload(GwyDataLine *x, GwyDataLine *y, Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata, const gdouble estimate[3])
{
    gdouble beta[3], std_beta[3], covariance[3];
    gint info;
    gdouble ux, uy;
    gint istart, iend, ndata = 0;
    gdouble r1 = 0.05, r2 = 0.1, rstep = 0.05;
    gint np;
    gdouble xmin;

    np = 3;

    //fit the data with the function beta[0]*(h-beta[1])^beta[2]
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

	beta[0] = fddata->Fmax / pow(fddata->hmax - beta[1], beta[2]);
    }

    // uncertainties given by sigma_h, sigma_F
    ux = instdata->sigma_h;
    uy = instdata->sigma_F;

    info = odr_run_fit(x, y, np, beta, sldata->unloadifixb, sldata->logfnmunload, ux, uy, std_beta, covariance,
                       &(sldata->unloadSSres), &(sldata->unloadR2), &(sldata->unloadR2adj), &(sldata->unloadchi2), POWER_3PARAM);
    //	printf("unload after beta %g %g %g  \n\n", beta[0], beta[1], beta[2]);

    // calculate results
    sldata->alpha = beta[0];
    sldata->hp = beta[1];
    sldata->m = beta[2];

    if (sldata->m <= 1) {
        sldata->eps = 0.75;
    }
    else {
        sldata->eps = epsilon(sldata->m);
    }

    sldata->Sunload = sldata->m * fddata->Fmax / (fddata->hmax - sldata->hp);

    if (sldata->has_fit_load) {
        slopes_combine_fit_results(sldata, fddata, instdata);
    }

    //clean up
    sldata->unloadrange[0] = sldata->unloadfrom;
    sldata->unloadrange[1] = sldata->unloadto;
    sldata->has_fit_unload = TRUE;
    sldata->unloadrange_pct_Fmax[0] = sldata->unloadfrom_pct_Fmax;
    sldata->unloadrange_pct_Fmax[1] = sldata->unloadto_pct_Fmax;

    sldata->fitinfounload = info;

    return info;
}


#endif

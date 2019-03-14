#ifndef NOFORTRAN

#include "tool_stiffness.h"

#include "datatypes.h"
#include "file-utils.h"
#include "fit-utils.h"
#include "niget-common.h"

#include <math.h>
#include <stdlib.h>


static const gchar *logfnmload = "fit.log.stiff.load";
static const gchar *logfnmunload = "fit.log.stiff.unload";


void init_Stiffnessdata(Stiffnessdata *stdata)
{
    stdata->kl = 0;
    stdata->ql = 0;
    stdata->ku = 0;
    stdata->qu = 0;

    stdata->has_fit_load = FALSE;
    stdata->loadto = 0;
    stdata->loadfrom = 0;
    stdata->loadrange[0] = 0;
    stdata->loadrange[1] = 0;

    stdata->has_fit_unload = FALSE;
    stdata->unloadto = 0;
    stdata->unloadfrom = 0;
    stdata->unloadrange[0] = 0;
    stdata->unloadrange[1] = 0;

    stdata->xloadfit = NULL;
    stdata->yloadfit = NULL;
    stdata->xunloadfit = NULL;
    stdata->yunloadfit = NULL;

    stdata->nfitloaddata = 0;
    stdata->nfitunloaddata = 0;

    /* quality of fit stuff */
    stdata->loadSSres = -1;
    stdata->loadR2 = -1;
    stdata->loadR2adj = -1;
    stdata->loadchi2 = -1;
    stdata->unloadSSres = -1;
    stdata->unloadR2 = -1;
    stdata->unloadR2adj = -1;
    stdata->unloadchi2 = -1;

    // TODO build logfilename, modify also in unc and mc
    stdata->logfnmload = create_log_filename_path(logfnmload);
    stdata->infologload = 0;

    stdata->logfnmunload = create_log_filename_path(logfnmunload);
    stdata->infologunload = 0;
}

gboolean stiffness_has_results(const Stiffnessdata *stiffness)
{
    return (stiffness->has_fit_load && stiffness->has_fit_unload);
}

void stiffness_remove_all_fit_results_unload(Stiffnessdata *stdata)
{
    //reset all results associated with the straight line fit of the unloading curve to zero
    stdata->ku = 0;
    stdata->qu = 0;

    stdata->unloadSSres = -1;
    stdata->unloadR2 = -1;
    stdata->unloadR2adj = -1;
    stdata->unloadchi2 = -1;

    stdata->has_fit_unload = FALSE;
    stdata->nfitunloaddata = 0;

    stdata->infologunload = 0;
}

void stiffness_remove_all_fit_results_load(Stiffnessdata *stdata)
{
    //reset all results associated with the straight line fit of the loading curve to zero
    stdata->kl = 0;
    stdata->ql = 0;

    stdata->has_fit_load = FALSE;
    stdata->nfitloaddata = 0;

    stdata->infologload = 0;
}

gchar *stiffness_export_data(const Stiffnessdata *stdata, const FDdata *fddata, const Instdata *instdata, enum ExportFormat exportformat)
{
    gint i;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    //Input parameters
    g_string_append_printf(buf, "#Parameters of Stiffness  analysis \n");
    g_string_append_printf(buf, "#sigma_h: %g \n", instdata->sigma_h);
    g_string_append_printf(buf, "#sigma_F: %g \n", instdata->sigma_F);
    g_string_append_printf(buf, "#delta: %g \n", instdata->delta);

    g_string_append_printf(buf, "\n\n");
    g_string_append_printf(buf, "#Results of Stiffness  analysis \n");
    g_string_append_printf(buf, "#hmax: %g nm\n", fddata->hmax);
    g_string_append_printf(buf, "#Fmax: %g mN\n", fddata->Fmax);

    g_string_append_printf(buf, "#Fit of loading curve     kl h + ql\n");
    g_string_append_printf(buf, "#kl: %g mN.nm^{-1}\n", stdata->kl);
    g_string_append_printf(buf, "#ql: %g nm\n", stdata->ql);

    g_string_append_printf(buf, "#fit info: %d \n", stdata->fitinfoload);
    g_string_append_printf(buf, "#SSres: %g \n", stdata->loadSSres);
    g_string_append_printf(buf, "#R-squared: %g \n", stdata->loadR2);
    g_string_append_printf(buf, "#adjusted R-squared: %g \n", stdata->loadR2adj);
    g_string_append_printf(buf, "#chi-squared: %g \n", stdata->loadchi2);

    g_string_append_printf(buf, "#Fit of unloading curve     ku h + qu\n");
    g_string_append_printf(buf, "#ku: %g mN.nm^{-1}\n", stdata->ku);
    g_string_append_printf(buf, "#qu: %g nm\n", stdata->qu);

    g_string_append_printf(buf, "#fit info: %d \n", stdata->fitinfounload);
    g_string_append_printf(buf, "#SSres: %g \n", stdata->unloadSSres);
    g_string_append_printf(buf, "#R-squared: %g \n", stdata->unloadR2);
    g_string_append_printf(buf, "#adjusted R-squared: %g \n", stdata->unloadR2adj);
    g_string_append_printf(buf, "#chi-squared: %g \n", stdata->unloadchi2);

    if (stdata->has_fit_load) {
        g_string_append_printf(buf, "\n \n");
        g_string_append_printf(buf, "#Fitted loading range: ");

        g_string_append_printf(buf, "[%g , %g ] nm \n", stdata->loadrange[0], stdata->loadrange[1]);
        g_string_append_printf(buf, "# %d datapoints used.\n", stdata->nfitloaddata);

        g_string_append_printf(buf, "\n# Loading Fit \n");
        g_string_append_printf(buf, "# x/nm y/mN \n");

        for (i = 0; i < stdata->xloadfit->res; i++) {
            g_string_append_printf(buf, "%g %g \n", stdata->xloadfit->data[i], stdata->yloadfit->data[i]);
        }
    }

    if (stdata->has_fit_unload) {
        g_string_append_printf(buf, "\n \n");
        g_string_append_printf(buf, "#Fitted unloading range: ");

        g_string_append_printf(buf, "[%g , %g ] nm \n", stdata->unloadrange[0], stdata->unloadrange[1]);
        g_string_append_printf(buf, "# %d datapoints used.\n", stdata->nfitunloaddata);

        g_string_append_printf(buf, "\n# Unloading Fit \n");
        g_string_append_printf(buf, "# x/nm y/mN \n");

        for (i = 0; i < stdata->xunloadfit->res; i++) {
            g_string_append_printf(buf, "%g %g \n", stdata->xunloadfit->data[i], stdata->yunloadfit->data[i]);
        }
    }

    return g_string_free(buf, FALSE);
}

void stiffness_combine_fit_results(Stiffnessdata *stdata, const FDdata *fddata, const Instdata *instdata)
{
}

gint stiffness_fit_load(GwyDataLine *x, GwyDataLine *y, Stiffnessdata *stdata, const FDdata *fddata, const Instdata *instdata, const gdouble estimate[2])
{
    gdouble beta[2], std_beta[2], covariance[2];
    gint info;
    gdouble ux, uy;
    // always fit both parameters
    const gint np = 2;
    const gint ifixb[2] = {1, 1};

    //fit the data with the function beta[0]*h+beta[1]
    if (estimate) {
        beta[0] = estimate[0];
        beta[1] = estimate[1];
    }
    else {
        beta[0] = 1;
        beta[1] = 1;
    }

    // uncertainties given by sigma_h, sigma_F
    ux = instdata->sigma_h;
    uy = instdata->sigma_F;

    info = odr_run_fit(x, y, np, beta, ifixb, stdata->logfnmload, ux, uy, std_beta, covariance,
                       &(stdata->loadSSres), &(stdata->loadR2), &(stdata->loadR2adj), &(stdata->loadchi2), LINEAR);

    // calculate results
    stdata->kl = beta[0];
    stdata->ql = beta[1];

    if (stdata->has_fit_unload) {
        stiffness_combine_fit_results(stdata, fddata, instdata);
    }

    //clean up
    stdata->loadrange[0] = stdata->loadfrom;
    stdata->loadrange[1] = stdata->loadto;
    stdata->has_fit_load = TRUE;
    stdata->fitinfoload = info;

    return info;
}

gint stiffness_fit_unload(GwyDataLine *x, GwyDataLine *y, Stiffnessdata *stdata, const FDdata *fddata, const Instdata *instdata, const gdouble estimate[2])
{
    gdouble beta[2], std_beta[2], covariance[2];
    gint info;
    gdouble ux, uy;
    /*always fit both parameters */
    gint ifixb[2] = {1, 1};
    gint np = 2;

    //fit the data with the function beta[0]*(h-beta[1])^beta[2]
    if (estimate) {
        beta[0] = estimate[0];
        beta[1] = estimate[1];
    }
    else {
        beta[0] = 1;
        beta[1] = 1;
    }

    // uncertainties given by sigma_h, sigma_F
    ux = instdata->sigma_h;
    uy = instdata->sigma_F;

    info = odr_run_fit(x, y, np, beta, ifixb, stdata->logfnmunload, ux, uy, std_beta, covariance,
                       &(stdata->unloadSSres), &(stdata->unloadR2), &(stdata->unloadR2adj), &(stdata->unloadchi2), LINEAR);

    // calculate results
    stdata->ku = beta[0];
    stdata->qu = beta[1];

    if (stdata->has_fit_load) {
        stiffness_combine_fit_results(stdata, fddata, instdata);
    }

    //clean up
    stdata->unloadrange[0] = stdata->unloadfrom;
    stdata->unloadrange[1] = stdata->unloadto;
    stdata->has_fit_unload = TRUE;
    stdata->fitinfounload = info;

    return info;
}

#endif

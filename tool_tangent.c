#include "tool_tangent.h"

#include "datatypes.h"
#include "file-utils.h"
#include "fit-utils.h"
#include "niget-common.h"

#include <math.h>
#include <stdlib.h>

void init_Tangentdata(Tangentdata *tgdata)
{
    tgdata->beta = BETA_DEFAULT_TG;
    tgdata->hc = 0;
    tgdata->S = 0;
    tgdata->Aphc = 0;
    tgdata->eps = EPS;
    tgdata->range[0] = 0;
    tgdata->range[1] = 0;
    tgdata->range_pct_Fmax[0] = 0;
    tgdata->range_pct_Fmax[1] = 0;

    tgdata->has_fit = FALSE;
    tgdata->from = 0;
    tgdata->to = 0;
    tgdata->from_pct_Fmax = 0;
    tgdata->to_pct_Fmax = 0;
    tgdata->xfit = NULL;
    tgdata->yfit = NULL;

    tgdata->Er = 0;
    tgdata->Eit = 0;
    tgdata->Hit = 0;

    tgdata->Finput = FALSE;
    tgdata->nfitdata = 0;

    tgdata->reg.SSres = -1;
    tgdata->reg.R2 = -1;
    tgdata->reg.R2adj = -1;
    tgdata->reg.chi2 = -1;
}

gboolean tangent_has_results(Tangentdata *tgdata)
{
    return tgdata->has_fit;
}

void tangent_remove_all_fit_results(Tangentdata *tgdata)
{
    //reset all results associated with a fit to zero
    tgdata->hc = 0;
    tgdata->S = 0;
    tgdata->Aphc = 0;
    tgdata->range[0] = 0;
    tgdata->range[1] = 0;
    tgdata->from = 0;
    tgdata->to = 0;
    tgdata->range_pct_Fmax[0] = 0;
    tgdata->range_pct_Fmax[1] = 0;
    tgdata->from_pct_Fmax = 0;
    tgdata->to_pct_Fmax = 0;
    tgdata->Hit = 0;
    tgdata->Eit = 0;
    tgdata->Er = 0;
    tgdata->reg.SSres = -1;
    tgdata->reg.R2 = -1;
    tgdata->reg.R2adj = -1;
    tgdata->reg.chi2 = -1;
}

gchar* tangent_export_data(Tangentdata *tgdata, FDdata *fddata, Instdata *instdata, Area *area, enum ExportFormat exportformat)
{
    gint i;
    gchar *areatext;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    //Input parameters
    g_string_append_printf(buf, "#Parameters of tangent analysis \n");
    g_string_append_printf(buf, "#beta: %g \n", tgdata->beta);
    g_string_append_printf(buf, "#nu: %g \n", instdata->nu);
    g_string_append_printf(buf, "#nui: %g \n", instdata->nui);
    g_string_append_printf(buf, "#Ei: %g GPa\n", instdata->Ei * 1e-9); //TODO units
    // Area
    areatext = write_area_file(area);
    g_string_append_printf(buf, "%s \n", areatext);
    g_free(areatext);

    //warning if extrapolation in area
    if (area->mode == AREA_DATA && tgdata->hc > area->xmax) {
        g_string_append_printf(buf, "#    WARNING: Extrapolation in area calibration! \n");
    }

    g_string_append_printf(buf, "\n\n");

    g_string_append_printf(buf, "#Results of tangent analysis \n");
    g_string_append_printf(buf, "#hmax: %g nm\n", fddata->hmax);
    g_string_append_printf(buf, "#Fmax: %g mN\n", fddata->Fmax);
    g_string_append_printf(buf, "#S: %g mN/nm\n", tgdata->S);
    g_string_append_printf(buf, "#Residual sum of squares: %g \n", tgdata->reg.SSres);
    g_string_append_printf(buf, "#R-squared: %g \n", tgdata->reg.R2);
    g_string_append_printf(buf, "#adjusted R-squared: %g \n", tgdata->reg.R2adj);
    g_string_append_printf(buf, "#chi-squared: %g \n", tgdata->reg.chi2);
    g_string_append_printf(buf, "#hc: %g nm\n", tgdata->hc);
    g_string_append_printf(buf, "#Aphc: %g nm^2\n", tgdata->Aphc);
    g_string_append_printf(buf, "#Hit: %g MPa\n", tgdata->Hit);
    g_string_append_printf(buf, "#Er: %g GPa\n", tgdata->Er);
    g_string_append_printf(buf, "#Eit: %g GPa\n", tgdata->Eit);

    if (tgdata->has_fit) {
        g_string_append_printf(buf, "\n \n");
        g_string_append_printf(buf, "#Fitted range: ");

        if (tgdata->Finput) {
            g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", tgdata->range_pct_Fmax[0], tgdata->range_pct_Fmax[1]);
        }
        else {
            g_string_append_printf(buf, "[%g , %g ] nm \n", tgdata->range[0], tgdata->range[1]);
        }

        g_string_append_printf(buf, "# corresponds : ");

        if (tgdata->Finput) {
            g_string_append_printf(buf, "[%g , %g ] nm \n", tgdata->range[0], tgdata->range[1]);
        }
        else {
            g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", tgdata->range_pct_Fmax[0], tgdata->range_pct_Fmax[1]);
        }

        g_string_append_printf(buf, "# %d datapoints used.\n", tgdata->nfitdata);

        g_string_append_printf(buf, "\n# S Fit \n");
        g_string_append_printf(buf, "# x/nm y/mN \n");

        for (i = 0; i < tgdata->xfit->res; i++)
            if (tgdata->xfit->data[i] >= tgdata->range[0] && tgdata->xfit->data[i] <= tgdata->range[1]) {
                g_string_append_printf(buf, "%g %g \n", tgdata->xfit->data[i], tgdata->yfit->data[i]);
            }

        g_string_append_printf(buf, "\n\n");

        g_string_append_printf(buf, "\n# tangent \n");
        g_string_append_printf(buf, "# x/nm y/mN \n");

        for (i = 0; i < tgdata->xfit->res; i++) {
            g_string_append_printf(buf, "%g %g \n", tgdata->xfit->data[i], tgdata->yfit->data[i]);
        }
    }

    return g_string_free(buf, FALSE);
}

void tangent_fit(GwyDataLine *x, GwyDataLine *y, Tangentdata *tgdata, FDdata *fddata, Instdata *instdata, Area *area, gdouble ux, gdouble uy)
{
    //fit the data with a straight line
    //total_least_squares_line_fit(x, y, &(tgdata->reg));
    if (FIT_TLS) {
        total_least_squares_line_fit__(x, y, &(tgdata->reg), ux, uy, CORRECT, RESCALE, DELT);
    }
    else {
        ordinary_least_squares_line_fit__(x, y, &(tgdata->reg), ux, uy, CORRECT, RESCALE);
    }

    // TODO units
    // calculate results
    tgdata->S = tgdata->reg.slope;
    tgdata->hc = fddata->hmax - tgdata->eps * fddata->Fmax / tgdata->S;
    tgdata->Aphc = eval_area(area, tgdata->hc);

    tgdata->Hit = fddata->Fmax / tgdata->Aphc * 1e9;
    tgdata->Er = sqrt(M_PI / tgdata->Aphc) * tgdata->S / (2 * tgdata->beta) * 1e15;
    tgdata->Eit = (1 - instdata->nu * instdata->nu) / (1 / tgdata->Er - (1 - instdata->nui * instdata->nui) / instdata->Ei);
    tgdata->Eit *= 1e-9;
    tgdata->Er *= 1e-9;

    //clean up
    tgdata->range[0] = tgdata->from;
    tgdata->range[1] = tgdata->to;
    tgdata->range_pct_Fmax[0] = tgdata->from_pct_Fmax;
    tgdata->range_pct_Fmax[1] = tgdata->to_pct_Fmax;
    tgdata->has_fit = TRUE;
}


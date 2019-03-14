#include "tool_hertz.h"
#include "tool_hertz_unc.h"

#include "datatypes.h"
#include "fit-utils.h"
#include "niget-common.h"

#include <math.h>
#include <stdlib.h>


#define INTERCEPT FALSE

void init_Hertzdata(Hertzdata *hzdata)
{
    hzdata->mode = R_MODE;
    hzdata->radius = 100e-9;
    hzdata->has_fit = FALSE;
    hzdata->from = 0;
    hzdata->to = 0;
    hzdata->xfit = NULL;
    hzdata->yfit = NULL;

    hzdata->range[0] = 0;
    hzdata->range[1] = 0;
    hzdata->Er = 0;
    hzdata->Eit = 0;
    hzdata->nfitdata = 0;

    hzdata->reg.SSres = -1;
    hzdata->reg.R2 = -1;
    hzdata->reg.R2adj = -1;
}

gboolean hertz_has_results(const Hertzdata *hzdata)
{
    return hzdata->has_fit;
}

void hertz_remove_all_fit_results(Hertzdata *hzdata)
{
    //reset all results associated with a fit to zero
    hzdata->range[0] = 0;
    hzdata->range[1] = 0;
    hzdata->from = 0;
    hzdata->to = 0;
    /* TODO je to nutne */
    hzdata->reg.slope = 0;
    hzdata->reg.SSres = -1;
    hzdata->reg.R2 = -1;
    hzdata->reg.R2adj = -1;
    hzdata->reg.chi2 = -1;
    hzdata->Eit = 0;
    hzdata->Er = 0;
    hzdata->nfitdata = 0;
}

gchar* hertz_export_data(const Hertzdata *hzdata, const FDdata *fddata, const Instdata *instdata, enum ExportFormat exportformat)
{
    gint i;
    GString *buf;
    
    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    //Input parameters
    g_string_append_printf(buf, "#Parameters of Hertz analysis \n");

    switch (hzdata->mode) {
    case R_MODE:
        g_string_append_printf(buf, "#nu: %g \n", instdata->nu);
        g_string_append_printf(buf, "#nui: %g \n", instdata->nui);
        g_string_append_printf(buf, "#Ei: %g GPa\n", instdata->Ei * 1e-9); //TODO units
        g_string_append_printf(buf, "#tip radius: %g nm\n", hzdata->radius * 1e9);
        break;

    case ER_MODE:
        g_string_append_printf(buf, "#Er: %g GPa\n", hzdata->Er);
        break;

    case EIT_MODE:
        g_string_append_printf(buf, "#nu: %g \n", instdata->nu);
        g_string_append_printf(buf, "#nui: %g \n", instdata->nui);
        g_string_append_printf(buf, "#Ei: %g GPa\n", instdata->Ei * 1e-9); //TODO units
        g_string_append_printf(buf, "#Eit: %g GPa\n", hzdata->Eit);
        break;

    default:
        g_printerr("Unhandled mode type \n");
        break;
    }

    g_string_append_printf(buf, "\n");

    if (hzdata->has_fit) {
        g_string_append_printf(buf, "#Results of Hertz analysis \n");
        g_string_append_printf(buf, "#a: %g mN.nm^{-3/2} \n", hzdata->reg.slope);

        g_string_append_printf(buf, "#SSres: %g \n", hzdata->reg.SSres);
        g_string_append_printf(buf, "#R-squared: %g \n", hzdata->reg.R2);
        g_string_append_printf(buf, "#adjusted R-squared: %g \n", hzdata->reg.R2adj);
        g_string_append_printf(buf, "#chi-squared: %g \n", hzdata->reg.chi2);

        switch (hzdata->mode) {
        case R_MODE:
            g_string_append_printf(buf, "#Er: %g GPa\n", hzdata->Er);
            g_string_append_printf(buf, "#Eit: %g GPa\n", hzdata->Eit);
            break;

        case ER_MODE:
        case EIT_MODE:
            g_string_append_printf(buf, "#tip radius: %g nm\n", hzdata->radius * 1e9);
            break;

        default:
            g_printerr("Unhandled mode type \n");
            break;
        }

        g_string_append_printf(buf, " \n");
        g_string_append_printf(buf, "#Fitted range: [%g , %g ] nm \n", hzdata->range[0], hzdata->range[1]);
        g_string_append_printf(buf, "# %d datapoints used.\n", hzdata->nfitdata);

        g_string_append_printf(buf, "\n \n ");
        g_string_append_printf(buf, "\n# S Fit \n");
        g_string_append_printf(buf, "# x/nm y/mN \n");

        for (i = 0; i < hzdata->xfit->res; i++) {
            g_string_append_printf(buf, "%g %g \n", hzdata->xfit->data[i], hzdata->yfit->data[i]);
        }
    }

    return g_string_free(buf, FALSE);
}

void hertz_fit(GwyDataLine *x, GwyDataLine *y, Hertzdata *hzdata, const Instdata *instdata)
{
    if (INTERCEPT) {
        //fit the data with  F=a*x^(3/2)+b
        least_squares_32power_fit_two(x, y, &(hzdata->reg)); //TODO
    }
    else {
        //fit the data with  F=a*x^(3/2)
        filter_negative_data(x, y);
        least_squares_32power_fit(x, y, instdata->sigma_h, instdata->sigma_F, &(hzdata->reg)); //TODO
    }

    // TODO units
    // calculate results according to mode

    //printf(" hzdata->mode %d\n", hzdata->mode);
    switch (hzdata->mode) {
    case R_MODE:
        hzdata->Er = calc_Er_Hertz(hzdata->reg.slope, hzdata->radius * 1e9) * 1e15;
        hzdata->Eit = calc_Eit(hzdata->Er, instdata->nu, instdata->Ei, instdata->nui) * 1e-9;
        hzdata->Er *= 1e-9;
        break;

    case EIT_MODE:
        hzdata->Er = calc_Er(hzdata->Eit, instdata->nu, instdata->Ei *1e-9, instdata->nui);
	/* RS - missing break? */
	
    case ER_MODE:
        hzdata->radius =  calc_R_Hertz(hzdata->reg.slope, hzdata->Er) * 1e3;
        break;
    }

    //	printf("radius %g slope %g  Er %g \n", hzdata->radius,hzdata->reg.slope, hzdata->Er);

    //clean up
    hzdata->range[0] = hzdata->from;
    hzdata->range[1] = hzdata->to;
    hzdata->has_fit = TRUE;
}

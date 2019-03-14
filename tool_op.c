#include "tool_op.h"

#include "datatypes.h"
#include "file-utils.h"
#include "fit-utils.h"
#include "niget-common.h"

#include <math.h>
#include <stdlib.h>

void init_OPdata(OPdata *opdata)
{
    opdata->beta = BETA_DEFAULT_OP;
    opdata->hp = 0;
    opdata->hc = 0;
    opdata->hr = 0;
    opdata->eps = 0;
    opdata->m = 0;
    opdata->S = 0;
    opdata->Aphc = 0;

    opdata->has_fit_hp = FALSE;
    opdata->hpto = 0;
    opdata->hpfrom = 0;
    opdata->hpto_pct_Fmax = 0;
    opdata->hpfrom_pct_Fmax = 0;
    opdata->hprange[0] = 0;
    opdata->hprange[1] = 0;
    opdata->hprange_pct_Fmax[0] = 0;
    opdata->hprange_pct_Fmax[1] = 0;
    opdata->has_fit_S = FALSE;
    opdata->Srange[0] = 0;
    opdata->Srange[1] = 0;
    opdata->Srange_pct_Fmax[0] = 0;
    opdata->Srange_pct_Fmax[1] = 0;
    opdata->Sto = 0;
    opdata->Sfrom = 0;
    opdata->Sto_pct_Fmax = 0;
    opdata->Sfrom_pct_Fmax = 0;
    opdata->xhpfit = NULL;
    opdata->yhpfit = NULL;
    opdata->xSfit = NULL;
    opdata->ySfit = NULL;

    opdata->Er = 0;
    opdata->Eit = 0;
    opdata->Hit = 0;

    opdata->Finputhp = FALSE;
    opdata->nfithpdata = 0;
    opdata->FinputS = FALSE;
    opdata->nfitSdata = 0;

    opdata->reghp.SSres = -1;
    opdata->reghp.R2 = -1;
    opdata->reghp.R2adj = -1;
    opdata->reghp.chi2 = -1;

    opdata->regS.SSres = -1;
    opdata->regS.R2 = -1;
    opdata->regS.R2adj = -1;
    opdata->regS.chi2 = -1;
}


gboolean op_has_results(const OPdata *opdata)
{
    return opdata->has_fit_hp;
}

void op_remove_all_fit_results_S(OPdata *opdata)
{
    //reset all results associated with the powerlaw (S) fit to zero
    opdata->hc = 0;
    opdata->m = 0;
    opdata->eps = 0;
    opdata->S = 0;
    opdata->Aphc = 0;
    opdata->Srange[0] = 0;
    opdata->Srange[1] = 0;
    opdata->Srange_pct_Fmax[0] = 0;
    opdata->Srange_pct_Fmax[1] = 0;
    opdata->Sfrom = 0;
    opdata->Sto = 0;
    opdata->Hit = 0;
    opdata->Eit = 0;
    opdata->Er = 0;
    opdata->regS.SSres = -1;
    opdata->regS.R2 = -1;
    opdata->regS.R2adj = -1;
    opdata->regS.chi2 = -1;
}

void op_remove_all_fit_results_hp(OPdata *opdata)
{
    //reset all results associated with the linear (hp) fit to zero, results associated with the subsequent powerlaw (S) fit must be also invalidated
    opdata->hp = 0;
    opdata->hprange[0] = 0;
    opdata->hprange[1] = 0;
    opdata->hprange_pct_Fmax[0] = 0;
    opdata->hprange_pct_Fmax[1] = 0;
    opdata->hpfrom = 0;
    opdata->hpto = 0;
    opdata->reghp.SSres = -1;
    opdata->reghp.R2 = -1;
    opdata->reghp.R2adj = -1;
    opdata->reghp.chi2 = -1;
    op_remove_all_fit_results_S(opdata);
}

/* exportformat presently ignored */
gchar* op_export_data(const OPdata *opdata, const FDdata *fddata, const Instdata *instdata, const Area *area,
		      enum ExportFormat exportformat)
{
    gchar *areatext;
    gint i;
    GString *buf;

    buf = g_string_new(NULL);

    // after removing the Args dependence, add !NULL checks for the provided pointers

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    //Input parameters
    g_string_append_printf(buf, "#Parameters of Oliver Pharr  analysis \n");
    g_string_append_printf(buf, "#beta: %g \n", opdata->beta);
    g_string_append_printf(buf, "#nu: %g \n", instdata->nu);
    g_string_append_printf(buf, "#nui: %g \n", instdata->nui);
    g_string_append_printf(buf, "#Ei: %g GPa\n", instdata->Ei * 1e-9); //TODO units
    // Area
    areatext = write_area_file(area);
    g_string_append_printf(buf, "%s \n", areatext);
    g_free(areatext);

    //warning if extrapolation in area
    if (area->mode == AREA_DATA && opdata->hc > area->xmax) {
        g_string_append_printf(buf, "#    WARNING: Extrapolation in area calibration! \n");
    }

    g_string_append_printf(buf, "\n\n");
    g_string_append_printf(buf, "#Results of OP  analysis \n");
    g_string_append_printf(buf, "#hmax: %g nm\n", fddata->hmax);
    g_string_append_printf(buf, "#Fmax: %g mN\n", fddata->Fmax);
    g_string_append_printf(buf, "#hp: %g nm\n", opdata->hp);
    g_string_append_printf(buf, "# hp fit quality:\n");
    g_string_append_printf(buf, "#Residual sum of squares: %g \n", opdata->reghp.SSres);
    g_string_append_printf(buf, "#R-squared: %g \n", opdata->reghp.R2);
    g_string_append_printf(buf, "#adjusted R-squared: %g \n", opdata->reghp.R2adj);
    g_string_append_printf(buf, "#chi-squared: %g \n", opdata->reghp.chi2);
    g_string_append_printf(buf, "#m: %g nm\n", opdata->m);
    g_string_append_printf(buf, "# S fit quality:\n");
    g_string_append_printf(buf, "#Residual sum of squares: %g \n", opdata->regS.SSres);
    g_string_append_printf(buf, "#R-squared: %g \n", opdata->regS.R2);
    g_string_append_printf(buf, "#adjusted R-squared: %g \n", opdata->regS.R2adj);
    g_string_append_printf(buf, "#chi-squared: %g \n", opdata->regS.chi2);

    g_string_append_printf(buf, "#epsilon: %g \n", opdata->eps);
    g_string_append_printf(buf, "#hc: %g nm\n", opdata->hc);
    g_string_append_printf(buf, "#S: %g mN/nm\n", opdata->S);
    g_string_append_printf(buf, "#Aphc: %g nm^2\n", opdata->Aphc);
    g_string_append_printf(buf, "#Hit: %g MPa\n", opdata->Hit);
    g_string_append_printf(buf, "#Er: %g GPa\n", opdata->Er);
    g_string_append_printf(buf, "#Eit: %g GPa\n", opdata->Eit);

    if (opdata->has_fit_hp) {
        g_string_append_printf(buf, "\n \n");
        g_string_append_printf(buf, "#Fitted hprange: ");

        if (opdata->Finputhp) {
            g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", opdata->hprange_pct_Fmax[0], opdata->hprange_pct_Fmax[1]);
        }
        else {
            g_string_append_printf(buf, "[%g , %g ] nm \n", opdata->hprange[0], opdata->hprange[1]);
        }

        g_string_append_printf(buf, "# corresponds : ");

        if (opdata->Finputhp) {
            g_string_append_printf(buf, "[%g , %g ] nm \n", opdata->hprange[0], opdata->hprange[1]);
        }
        else {
            g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", opdata->hprange_pct_Fmax[0], opdata->hprange_pct_Fmax[1]);
        }

        g_string_append_printf(buf, "# %d datapoints used.\n", opdata->nfithpdata);

        g_string_append_printf(buf, "\n# hp Fit \n");
        g_string_append_printf(buf, "# x/nm y/mN \n");

        for (i = 0; i < opdata->xhpfit->res; i++) {
            g_string_append_printf(buf, "%g %g \n", opdata->xhpfit->data[i], opdata->yhpfit->data[i]);
        }
    }

    if (opdata->has_fit_S) {
        g_string_append_printf(buf, "\n \n");
        g_string_append_printf(buf, "#Fitted S range: ");

        if (opdata->FinputS) {
            g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", opdata->Srange_pct_Fmax[0], opdata->Srange_pct_Fmax[1]);
        }
        else {
            g_string_append_printf(buf, "[%g , %g ] nm \n", opdata->Srange[0], opdata->Srange[1]);
        }

        g_string_append_printf(buf, "# corresponds : ");

        if (opdata->FinputS) {
            g_string_append_printf(buf, "[%g , %g ] nm \n", opdata->Srange[0], opdata->Srange[1]);
        }
        else {
            g_string_append_printf(buf, "[%g , %g ] %% Fmax \n", opdata->Srange_pct_Fmax[0], opdata->Srange_pct_Fmax[1]);
        }

        g_string_append_printf(buf, "# %d datapoints used.\n", opdata->nfitSdata);

        g_string_append_printf(buf, "\n# S Fit \n");
        g_string_append_printf(buf, "# x/nm y/mN \n");

        for (i = 0; i < opdata->xSfit->res; i++) {
            g_string_append_printf(buf, "%g %g \n", opdata->xSfit->data[i], opdata->ySfit->data[i]);
        }
    }

    return g_string_free(buf, FALSE);
}

void op_fit_hp(GwyDataLine *x, GwyDataLine *y, OPdata *opdata, const FDdata *fddata)
{
    //fit the data with a line, equal weights for x and y => delta=1
    total_least_squares_line_fit_delta(x, y, &(opdata->reghp), 1);

    opdata->hp = -opdata->reghp.intercept / opdata->reghp.slope;

    //clean up
    opdata->hprange[0] = opdata->hpfrom;
    opdata->hprange[1] = opdata->hpto;
    opdata->hprange_pct_Fmax[0] = opdata->hpfrom_pct_Fmax;
    opdata->hprange_pct_Fmax[1] = opdata->hpto_pct_Fmax;
    opdata->has_fit_hp = TRUE;
}

void op_fit_S_lr(GwyDataLine *x, GwyDataLine *y, OPdata *opdata, const FDdata *fddata, const Instdata *instdata, const Area *area)
{
    // fit the log-log data with a straight line, equal weights x and y => delta=1
    total_least_squares_line_fit_delta(x, y, &(opdata->regS), 1);

    opdata->m = opdata->regS.slope;

    // calculate results
    opdata->eps = epsilon(opdata->m);
    opdata->S = calc_S(fddata->hmax, fddata->Fmax, opdata->hp, opdata->m);
    opdata->hc = calc_hc(fddata->hmax, fddata->Fmax, opdata->S, opdata->eps);
    opdata->Aphc = eval_area(area, opdata->hc);

    // TODO units
    // calculate results
    opdata->Hit = calc_Hit(fddata->Fmax, opdata->Aphc) * 1e9;
    opdata->Er = calc_Er_OP(opdata->Aphc, opdata->S, opdata->beta) * 1e15;
    opdata->Eit = calc_Eit(opdata->Er, instdata->nu, instdata->Ei, instdata->nui);
    opdata->Eit *= 1e-9;
    opdata->Er *= 1e-9;

    //clean up
    opdata->Srange[0] = opdata->Sfrom;
    opdata->Srange[1] = opdata->Sto;
    opdata->Srange_pct_Fmax[0] = opdata->Sfrom_pct_Fmax;
    opdata->Srange_pct_Fmax[1] = opdata->Sto_pct_Fmax;
    opdata->has_fit_S = TRUE;
}

void op_fit_S_nl(GwyDataLine *x, GwyDataLine *y, OPdata *opdata, const FDdata *fddata, const Instdata *instdata, const Area *area)
{
    gdouble beta[3];

    // make a nonlinear fit of the data with the function beta[0]*(h-beta[1])^beta[2]

    beta[0] = fddata->Fmax / pow(fddata->hmax - opdata->hp, 1.5);
    beta[1] = opdata->hp;
    beta[2] = 1.5;

    nonlinear_least_squares_power_fit(x, y, beta);

    opdata->alpha = beta[0];
    opdata->hp = beta[1];
    opdata->m = beta[2];
    
    //calculate results
    opdata->eps = epsilon(opdata->m);
    opdata->S = calc_S(fddata->hmax, fddata->Fmax, opdata->hp, opdata->m);
    opdata->hc = calc_hc(fddata->hmax, fddata->Fmax, opdata->S, opdata->eps);
    opdata->Aphc = eval_area(area, opdata->hc);

    // TODO units
    // calculate results
    opdata->Hit = calc_Hit(fddata->Fmax, opdata->Aphc) * 1e9;
    opdata->Er = calc_Er_OP(opdata->Aphc, opdata->S, opdata->beta) * 1e15;
    opdata->Eit = calc_Eit(opdata->Er, instdata->nu, instdata->Ei, instdata->nui);
    opdata->Eit *= 1e-9;
    opdata->Er *= 1e-9;

    //clean up
    opdata->Srange[0] = opdata->Sfrom;
    opdata->Srange[1] = opdata->Sto;
    opdata->Srange_pct_Fmax[0] = opdata->Sfrom_pct_Fmax;
    opdata->Srange_pct_Fmax[1] = opdata->Sto_pct_Fmax;
    opdata->has_fit_S = TRUE;
}

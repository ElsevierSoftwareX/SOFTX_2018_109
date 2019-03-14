#include "tool_tangent_unc.h"
#include "tool_tangent.h"

#include "datatypes.h"
#include "fddata.h"
#include "fit-utils.h"
#include "niget-common.h"
#include "mc-utils.h"

#include <math.h>


#define DEBUG FALSE
#define DEBUG2 FALSE //save individual data (h,F) which were fitted, TAKES TIME!!!!

static gboolean nparam_to_string(gint i, gchar *str, gint nstr);

void init_TangentUncdata(TangentUncdata *tgunc)
{
    tgunc->uh = 1;
    tgunc->uF = 0.001;
    /*
      tgunc->unu = 0.05;
      tgunc->unui = 0.04;
      tgunc->uEi = 20e9;
    */
    tgunc->Nmc = N_MC;

    tgunc->uSuh = 0;
    tgunc->uSuF = 0;
    tgunc->uS = 0;

    tgunc->uhcuh = 0;
    tgunc->uhcuF = 0;
    tgunc->uhc = 0;

    tgunc->uAuh = 0;
    tgunc->uAuF = 0;
    tgunc->uA = 0;

    tgunc->uHituh = 0;
    tgunc->uHituF = 0;
    tgunc->uHit = 0;

    tgunc->uEruh = 0;
    tgunc->uEruF = 0;
    tgunc->uEr = 0;

    tgunc->uEituh = 0;
    tgunc->uEituF = 0;
    tgunc->uEitunu = 0;
    tgunc->uEitunui = 0;
    tgunc->uEituEi = 0;
    tgunc->uEit = 0;

    tgunc->w = NULL;
    tgunc->Sc = NULL;
    tgunc->Ec = NULL;
    tgunc->Hc = NULL;
}

void tangent_propagate_uncertainties_calc(Tangentdata *tg, TangentUncdata *unc, FDdata *fd, Instdata *inst, Area *area)
{
    gdouble uh, uF;
    gdouble dAdh;
    gdouble S, eps, hc,  hmax, Fmax, Aphc, Hit, Er, Eit, nu, nui, Ei;
    gint istart, iend, ndata;
    GwyDataLine *x, *y;

    //make copies to simplify code
    uh = unc->uh;
    uF = unc->uF;
    eps = tg->eps;
    S = tg->S;
    hc = tg->hc;
    hmax = fd->hmax;
    Fmax = fd->Fmax;
    Aphc = tg->Aphc;
    Hit = tg->Hit;
    Eit = tg->Eit;
    Er = tg->Er;
    Ei = inst->Ei;
    nui = inst->nui;
    nu = inst->nu;

    if (verbose) {
        g_print(" \n\nPropagate uncertainties\n");
    }

    /* constant uncertainty for all data points */
    if (tg->Finput) {
        range_to_indices(tg->range_pct_Fmax[0] * 0.01 * fd->Fmax, tg->range_pct_Fmax[1] * 0.01 * fd->Fmax, fd->Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(tg->range[0], tg->range[1], fd->hunload, TRUE, &istart, &iend, &ndata);
    }

    x = gwy_data_line_part_extract(fd->hunload, istart, ndata);
    y = gwy_data_line_part_extract(fd->Funload, istart, ndata);

    if (FIT_TLS) {
        total_least_squares_line_fit_const(x, y, &(tg->reg), uh, uF);
    }
    else {
        ordinary_least_squares_line_fit_unc_const(x, y, &(tg->reg), uh, uF);
    }

    g_object_unref(x);
    g_object_unref(y);

    unc->uSuh = tg->reg.unc_slope_x;
    unc->uSuF = tg->reg.unc_slope_y;

    if (verbose) {
        g_print("u(S;h) %g u(S;F) %g  urel(S;h) %g urel(S;F) %g  \n", unc->uSuh, unc->uSuF, unc->uSuh / S, unc->uSuF / S);
    }

    /*    u(hc)   */
    unc->uhcuh = sqrt(sq(uh) + sq(hmax - hc) / sq(S) * sq(unc->uSuh));
    unc->uhcuF = sqrt(sq(eps) / sq(S) * sq(uF) + sq(hmax - hc) / sq(S) * sq(unc->uSuF));
    /* u(hc;h) = sqrt[uh^2+ ((hmax-hc)/S^2)^2 uSh^2]
     * u(hc;F) = sqrt[(eps/S) uF^2 + ((hmax-hc)/S^2)^2 uSF^2 )
     */

    if (verbose) {
        g_print("u(hc;h) %g u(hc;F) %g  urel(hc;h) %g urel(hc;F) %g  \n", unc->uhcuh, unc->uhcuF, unc->uhcuh / hc, unc->uhcuF / hc);
    }

    /*    u(A)   */
    dAdh =  eval_area_derivative(area, hc) ;
    unc->uAuh = dAdh * unc->uhcuh;
    unc->uAuF = dAdh * unc->uhcuF;
    /* u(Ap;h) = dA/dh * u(hc;h)
     * u(Ap;F) = dA/dh * u(hc;F)
     */

    if (verbose) {
        g_print("u(A;h) %g u(A;F) %g urel(A;h) %g urel(A;F) %g \n", unc->uAuh, unc->uAuF, unc->uAuh / Aphc, unc->uAuF / Aphc);
    }

    /*     u(H_IT)   */
    unc->uHituh = Hit / Aphc * dAdh * unc->uhcuh;
    unc->uHituF = Hit * sqrt(sq(uF) / sq(Fmax) + sq(dAdh) / sq(Aphc) * sq(unc->uhcuF));
    /* u(Hit;h) = Fmax/Aphc^2 * dA/dh *u(hc;h)
     * u(Hit;F) = sqrt[ (1/Aphc)^2 u(F)^2 + (Fmax/Aphc^2 * dA/dh)^2 *u(hc;F)^2] = Hit sqrt[ u(F)^2/Fmax^2 + (dA/dh  / Aphc)^2 u(hc;F)^2 ]
     */

    if (verbose) {
        g_print("u(H;h) %g  u(H;F) %g  urel(H;h) %g  urel(H;F) %g\n", unc->uHituh, unc->uHituF, unc->uHituh / Hit, unc->uHituF / Hit);
    }

    /*     u(E_r)   */
    unc->uEruh = Er * sqrt(sq(unc->uSuh) / sq(S) + sq(unc->uAuh) / (sq(Aphc) * 4) - dAdh / Aphc * eps * Fmax / S * sq(unc->uSuh) / sq(S));
    unc->uEruF = Er * sqrt(sq(unc->uSuF) / sq(S) + sq(unc->uAuF) / (sq(Aphc) * 4) - dAdh / Aphc * eps * Fmax / S * sq(unc->uSuF) / sq(S));
    /* u(Er;h) = sqrt( (Er/S)^2 u(S;h)^2 + (Er/2A)^2 u(A;h)^2 - Er^2/(S Aphc) *dAdh*eps*Fmax/S^2 *uSh^2 )
     * u(Er;F) = sqrt( (Er/S)^2 u(S;F)^2 + (Er/2A)^2 u(A;F)^2 - Er^2/(S Aphc) *dAdh*eps*Fmax/S^2 *uSh^2 )
     */

    if (verbose) {
        g_print("u(Er;h)^2 terms %g %g %g \n",
                sq(unc->uSuh) / sq(S), sq(unc->uAuh) / (sq(Aphc) * 4), -dAdh / Aphc * eps * Fmax / S * sq(unc->uSuh) / sq(S));
        g_print("u(Er;F)^2 terms %g %g %g \n",
                sq(unc->uSuF) / sq(S), sq(unc->uAuF) / (sq(Aphc) * 4), -dAdh / Aphc * eps * Fmax / S * sq(unc->uSuF) / sq(S));
        g_print("u(Er;h) %g  u(Er;F) %g  urel(Er;h) %g  urel(Er;F) %g\n", unc->uEruh, unc->uEruF, unc->uEruh / Er, unc->uEruF / Er);
    }

    /*     u(E_IT)   */
    unc->uEituh = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruh;
    unc->uEituF = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruF;
    unc->uEitunu = fabs(2 * nu / (1 - sq(nu)) * Eit * inst->unu);
    unc->uEitunui = fabs(2 * nui / (1 - sq(nu)) * sq(Eit) / Ei * inst->unui) * 1e9;
    unc->uEituEi = fabs((1 - sq(nui)) / (1 - sq(nu)) * sq(Eit) / sq(Ei) * inst->uEi) * 1e9;
    /* u(Eit;h) =  Eit^2/Er^2 /(1-nu^2) u(Er;h)
     * u(Eit;F) =  Eit^2/Er^2 /(1-nu^2) u(Er;F)
     * u(Eit;nu) = -2nu/(1-nu^2) Eit u(nu)
     * u(Eit;nui) =  -2nui/(1-nu^2) Eit^2/Ei u(nui)
     * u(Eit;Ei) = -(1-nui^2)/(1-nu)^2 Eit^2/Ei^2 u(Ei)
     */

    if (verbose) {
        g_print("u(Eit;h) %g  u(Eit;F) %g  urel(Eit;h) %g  urel(Eit;F) %g\n", unc->uEituh, unc->uEituF, unc->uEituh / Eit, unc->uEituF / Eit);
        g_print("u(Eit;nu) %g  u(Eit;nui) %g  u(Eit;Ei) %g  urel(Eit;nu) %g  urel(Eit;nui) %g  urel(Eit;Ei) %g\n",
                unc->uEitunu, unc->uEitunui, unc->uEituEi, unc->uEitunu / Eit, unc->uEitunui / Eit, unc->uEituEi / Eit);
    }

    if (verbose) {
        g_print(" \nEnd propagation of uncertainties\n\n");
    }

    // total uncertainties
    unc->uS = sqrt(sq(unc->uSuh) + sq(unc->uSuF));
    unc->uhc = sqrt(sq(unc->uhcuh) + sq(unc->uhcuF));
    unc->uA = sqrt(sq(unc->uAuh) + sq(unc->uAuF));
    unc->uHit = sqrt(sq(unc->uHituh) + sq(unc->uHituF));
    unc->uEr = sqrt(sq(unc->uEruh) + sq(unc->uEruF));
    unc->uEit = sqrt(sq(unc->uEituh) + sq(unc->uEituF) + sq(unc->uEitunu) + sq(unc->uEitunui) + sq(unc->uEituEi));
}

gchar* tangent_uncertainty_export_data(Tangentdata *tgdata, TangentUncdata *tgunc, FDdata *fddata, Instdata *instdata, enum ExportFormat exportformat)
{
    gint j;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Parameters of Tangent analysis \n");
    g_string_append_printf(buf, "# u(depth): %g nm\n", tgunc->uh);
    g_string_append_printf(buf, "# u(load): %g nm\n", tgunc->uF);
    g_string_append_printf(buf, "#beta: %g \n", tgdata->beta);
    g_string_append_printf(buf, "#nu: %g                 u(nu):  %g  \n", instdata->nu, instdata->unu);
    g_string_append_printf(buf, "#nu_indenter: %g       u(nu_indenter):  %g   \n", instdata->nui, instdata->unui);
    g_string_append_printf(buf, "#E_indenter: %g        u(E_indenter):   %g GPa\n", instdata->Ei * 1e-9, instdata->uEi * 1e-9); //TODO units
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

    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#Results of tangent analysis \n");
    g_string_append_printf(buf, "#S = %g mN/nm\n", tgdata->S);
    g_string_append_printf(buf, "#hc = %g nm\n", tgdata->hc);
    g_string_append_printf(buf, "#Ap(hc) = %g nm^2\n", tgdata->Aphc);
    g_string_append_printf(buf, "#Er = %g GPa\n", tgdata->Er);
    g_string_append_printf(buf, "#Hit = %g GPa\n", tgdata->Hit);
    g_string_append_printf(buf, "#Eit = %g GPa\n", tgdata->Eit);
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Results of the Gaussian propagation of uncertainties \n");
    g_string_append_printf(buf, "#u(S) =  %g mN/nm \n", tgunc->uS);
    g_string_append_printf(buf, "#u(hc) =  %g nm \n", tgunc->uhc);
    g_string_append_printf(buf, "#u(Ap(hc)) =  %g nm^2 \n", tgunc->uA);
    g_string_append_printf(buf, "#u(Hit) =  %g MPa \n", tgunc->uHit);
    g_string_append_printf(buf, "#u(Er) =  %g GPa \n", tgunc->uEr);
    g_string_append_printf(buf, "#u(Eit) =  %g GPa \n", tgunc->uEit);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of S\n");
    g_string_append_printf(buf, "#u(S;uh) = %g mN/nm \n", tgunc->uSuh);
    g_string_append_printf(buf, "#u(S;uF) = %g mN/nm \n", tgunc->uSuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of hc\n");
    g_string_append_printf(buf, "#u(hc;uh) = %g nm \n", tgunc->uhcuh);
    g_string_append_printf(buf, "#u(hc;uF) = %g nm \n", tgunc->uhcuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Ap(hc)\n");
    g_string_append_printf(buf, "#u(Ap(hc);uh) = %g nm^2 \n", tgunc->uAuh);
    g_string_append_printf(buf, "#u(Ap(hc);uF) = %g nm^2 \n", tgunc->uAuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Hit\n");
    g_string_append_printf(buf, "#u(Hit;uh) = %g MPa \n", tgunc->uHituh);
    g_string_append_printf(buf, "#u(Hit;uF) = %g MPa \n", tgunc->uHituF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Er\n");
    g_string_append_printf(buf, "#u(Er;uh) = %g GPa \n", tgunc->uEruh);
    g_string_append_printf(buf, "#u(Er;uF) = %g GPa \n", tgunc->uEruF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Eit\n");
    g_string_append_printf(buf, "#u(Eit;uh) = %g GPa \n", tgunc->uEituh);
    g_string_append_printf(buf, "#u(Eit;uF) = %g GPa \n", tgunc->uEituF);
    g_string_append_printf(buf, "#u(Eit;unu) = %g GPa \n", tgunc->uEitunu);
    g_string_append_printf(buf, "#u(Eit;unui) = %g GPa \n", tgunc->uEitunui);
    g_string_append_printf(buf, "#u(Eit;uEi) = %g GPa \n", tgunc->uEituEi);
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Changes of the S corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#S/(mN/nm) \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (tgunc->Sc[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "#d    not defined \n");
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, tgunc->Sc[j + CONTACTMAX]);
        }
    }

    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#Changes of the Er corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#Er/GPa \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (tgunc->Ec[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "#d    not defined \n");
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, tgunc->Ec[j + CONTACTMAX]);
        }
    }

    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#Changes of the H corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#H/MPa \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (tgunc->Hc[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, tgunc->Hc[j + CONTACTMAX]);
        }
    }

    return g_string_free(buf, FALSE);
}

void tangent_fit_shift_contact(Tangentdata *tgdata, FDdata *fddata, Instdata *inst, Area *area, gdouble *Sc, gdouble *Ec, gdouble *Hc, gint j)
{
    Tangentdata tg;
    FDdata fdc;
    gdouble dh, dF;
    gint istart, iend, ndata;
    GwyDataLine *x, *y;

    Sc[j + CONTACTMAX] = -1;
    Ec[j + CONTACTMAX] = -1;
    Hc[j + CONTACTMAX] = -1;

    // create shifted FDdata

    if (verbose) {
        g_print("\n\nCreate shifted data for j %d \n", j);
    }

    if (fddata->i_contact_load + j < 0) {
        return;
    }

    fdc.i_contact_load = fddata->i_contact_load + j;

    fdc.horig = gwy_data_line_duplicate(fddata->horig);
    fdc.Forig = gwy_data_line_duplicate(fddata->Forig);

    fdc.hunload = gwy_data_line_duplicate(fddata->hunload);
    fdc.Funload = gwy_data_line_duplicate(fddata->Funload);

    /*
       gwy_data_line_add(fdc.hunload, -fdc.hall_orig->data[fdc.i_contact_load]);
       gwy_data_line_add(fdc.Funload, -fdc.Fall_orig->data[fdc.i_contact_load]);

       fdc.hmax = fddata.hmax- fdc.hall_orig->data[fdc.i_contact_load];
       fdc.Fmax = fddata.Fmax- fdc.Fall_orig->data[fdc.i_contact_load];
       */

    //difference in shift between different contactpoints
    dh = fdc.horig->data[fdc.i_contact_load] - fddata->horig->data[fddata->i_contact_load];
    dF = fdc.Forig->data[fdc.i_contact_load] - fddata->Forig->data[fddata->i_contact_load];

    if (verbose) {
        g_print("shift dh %g dF %g \n", dh, dF);
    }

    gwy_data_line_add(fdc.hunload, -dh);
    gwy_data_line_add(fdc.Funload, -dF);

    fdc.hmax = fddata->hmax - dh;
    fdc.Fmax = fddata->Fmax - dF;

    // create working copy of Tangentdata
    // shift fitting range
    tg = *tgdata;
    tg.range[0] += -dh;
    tg.range[1] += -dh;
    tg.from = tg.range[0];
    tg.to = tg.range[1];

    //check range
    if (tgdata->Finput) {
        if (verbose) {
            g_print(" from %g to %g \n", tg.from_pct_Fmax, tg.to_pct_Fmax);
        }

        if (tg.from_pct_Fmax == tg.to_pct_Fmax) {
            return;
        }
    }
    else {
        if (verbose) {
            g_print(" from %g to %g \n", tg.from, tg.to);
        }

        if (tg.from == tg.to) {
            return;
        }
    }

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (tgdata->Finput) {
        range_to_indices(tg.range_pct_Fmax[0] * 0.01 * fdc.Fmax, tg.range_pct_Fmax[1] * 0.01 * fdc.Fmax, fdc.Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(tg.range[0], tg.range[1], fdc.hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, not fitting can be performed
    ndata = MIN(ndata, fdc.hunload->res - istart);

    if (verbose)
        g_print("interval has %d points, from %d to %d, hrange [%g, %g ] nm, Frange [%g, %g] mN \n", ndata, istart, iend,
                fdc.hunload->data[iend], fdc.hunload->data[istart], fdc.Funload->data[iend], fdc.Funload->data[istart]);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    x = gwy_data_line_part_extract(fdc.hunload, istart, ndata);
    y = gwy_data_line_part_extract(fdc.Funload, istart, ndata);

    // perform fit and get results
    tangent_fit(x, y, &tg, &fdc, inst, area, 0, 0);

    // clean up
    g_object_unref(x);
    g_object_unref(y);
    g_object_unref(fdc.horig);
    g_object_unref(fdc.Forig);
    g_object_unref(fdc.hunload);
    g_object_unref(fdc.Funload);

    Sc[j + CONTACTMAX] = tg.S;
    Ec[j + CONTACTMAX] = tg.Er;
    Hc[j + CONTACTMAX] = tg.Hit;

    if (verbose) {
        g_print("\n\n");
    }
}

void init_TangentMCdata(TangentMCdata *tgmc, gint Nmc)
{
    gint i;

    //number of iterations
    tgmc->Nmc = Nmc;

    //how much had to be skipped
    tgmc->skipmc = 0;

    //dynamically allocate
    tgmc->tgmc = (Tangentdata *)g_malloc(tgmc->Nmc * sizeof(Tangentdata));
    tgmc->mcavg = (gdouble *)g_malloc(NPARAM_TG * sizeof(gdouble));
    tgmc->mcstd = (gdouble *)g_malloc(NPARAM_TG * sizeof(gdouble));
    tgmc->mcdd = (gdouble **) g_malloc(NPARAM_TG * sizeof(gdouble *));
    tgmc->mcxhist = (gdouble **) g_malloc(NPARAM_TG * sizeof(gdouble *));
    tgmc->mcyhist = (gdouble **) g_malloc(NPARAM_TG * sizeof(gdouble *));

    tgmc->mcnstat = (gint)floor(3.49 * cbrt(tgmc->Nmc) + 0.5);

    for (i = 0; i < NPARAM_TG; i++) {
        tgmc->mcdd[i] = (gdouble *) g_malloc(tgmc->Nmc * sizeof(gdouble));
        tgmc->mcxhist[i] = (gdouble *) g_malloc(tgmc->mcnstat * sizeof(gdouble));
        tgmc->mcyhist[i] = (gdouble *) g_malloc(tgmc->mcnstat * sizeof(gdouble));
    }
}

void destroy_TangentMCdata(TangentMCdata *tgmc)
{
    gint i;

    //free everything
    g_free(tgmc->mcavg);
    g_free(tgmc->mcstd);

    for (i = 0; i < NPARAM_TG; i++) {
        g_free(tgmc->mcdd[i]);
        g_free(tgmc->mcxhist[i]);
        g_free(tgmc->mcyhist[i]);
    }

    g_free(tgmc->mcdd);
    g_free(tgmc->mcxhist);
    g_free(tgmc->mcyhist);

    tgmc->mcnstat = 0;
    tgmc->skipmc = 0;
    g_free(tgmc->tgmc);
}

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

    default:
        g_printerr("Unhandled number of results.\n");
        return FALSE;
    }

    return TRUE;
}

gchar* tangent_mc_export_data(Tangentdata *tgdata, TangentUncdata *tgunc, TangentMCdata *tgmc, FDdata *fddata, Instdata *instdata, enum ExportFormat exportformat)
{
    Tangentdata *tg;
    gint i, j;
    gchar str[200];
    gboolean is;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Monte Carlo simulation of tangent analysis\n");
    g_string_append_printf(buf, "#Input uncertainties \n");
    g_string_append_printf(buf, "#u(h) = %g nm \n", tgunc->uh);
    g_string_append_printf(buf, "#u(F) = %g mN \n", tgunc->uF);
    g_string_append_printf(buf, "#Number of iterations  = %d  \n", tgmc->Nmc);
    g_string_append_printf(buf, "\n");

    g_string_append_printf(buf, "#Parameters of Tangent analysis \n");
    g_string_append_printf(buf, "#beta: %g \n", tgdata->beta);
    g_string_append_printf(buf, "#nu: %g                 \n", instdata->nu);
    g_string_append_printf(buf, "#nu_indenter: %g        \n", instdata->nui);
    g_string_append_printf(buf, "#E_indenter: %g      GPa\n", instdata->Ei * 1e-9); //TODO units
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

    g_string_append_printf(buf, "\n\n");

    g_string_append_printf(buf, "#Results \n");

    for (i = 0; i < NPARAM_TG; i++) {
        is = nparam_to_string(i, str, sizeof(str));
        g_string_append_printf(buf, "#%s = (%g ± %g ) \n", str, tgmc->mcavg[i], tgmc->mcstd[i]);
    }

    g_string_append_printf(buf, "\n\n\n");

    // Histograms

    if (tgmc->mcnstat > 0) {
        g_string_append_printf(buf, "#Histogram\n");

        for (i = 0; i < NPARAM_TG; i++) {
            g_string_append_printf(buf, "#Histogram\n");
            is = nparam_to_string(i, str, sizeof(str));

            if (is) {
                g_string_append_printf(buf, "#%s    count/a.u.  \n", str);

                for (j = 0; j < tgmc->mcnstat; j++) {
                    g_string_append_printf(buf, "%g %g \n", tgmc->mcxhist[i][j], tgmc->mcyhist[i][j]);
                }

                g_string_append_printf(buf, "\n \n");
            }
        }
    }

    g_string_append_printf(buf, "\n\n\n");

    //Full data
    g_string_append_printf(buf, "#Full data \n");
    g_string_append_printf(buf, "#");

    for (j = 0; j < NPARAM_TG; j++) {
        is = nparam_to_string(j, str, sizeof(str)); /* _is_ is not used anywhere */
        g_string_append_printf(buf, "%s ", str);
    }

    g_string_append_printf(buf, "\n");

    tg = tgmc->tgmc;

    for (i = 0; i < tgmc->Nmc; i++) {
        g_string_append_printf(buf, DEPTH, tg->hc);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, AREA, tg->Aphc);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, HARD, tg->Hit);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, MODUL, tg->Er);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, MODUL, tg->Eit);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, SLOPE, tg->S);
        g_string_append_printf(buf, "  \n");
        tg++;
    }

    return g_string_free(buf, FALSE);
}

void tangent_uncertainty_montecarlo_run_calc(Tangentdata *tgdata, TangentUncdata *tgunc, TangentMCdata *tgmc, FDdata *fddata, Instdata *instdata, Area *area)
{
    gint istart, iend, ndata;
    gint i, Nmc, j;
    Tangentdata *tg;
    TangentMCdata *mc;
    GwyDataLine *h, *F;
    GwyDataLine *hnew, *Fnew;
    GwyDataLine *hnoise, *Fnoise;
    gdouble hmax_old, Fmax_old;
    GRand *rng;
    gint nn;
    gdouble **dd;
    gdouble *xhist, *yhist;
    gboolean hist;
    FILE *aux;
    FILE *aux2;
    gchar auxstr[1000];

    //not necessary to check range, couldn't get here in the first place

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (tgdata->Finput)
        range_to_indices(tgdata->from_pct_Fmax * 0.01 * fddata->Fmax, tgdata->to_pct_Fmax * 0.01 * fddata->Fmax,
                         fddata->Funload, TRUE, &istart, &iend, &ndata);
    else {
        range_to_indices(tgdata->from, tgdata->to, fddata->hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, not fitting can be performed
    ndata = MIN(ndata, fddata->hunload->res - istart);

    if (ndata < 2) {
        return;
    }

    //back up hmax, Fmax
    hmax_old = fddata->hmax;
    Fmax_old = fddata->Fmax;

    //get data to be fitted
    h = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
    F = gwy_data_line_part_extract(fddata->Funload, istart, ndata);

    Nmc = tgmc->Nmc;
    tg = tgmc->tgmc;

    //prepare data
    for (i = 0; i < Nmc ; i++) {
        init_Tangentdata(tg);
        tg++;
    }

    tg = tgmc->tgmc;

    //copy input parameters
    for (i = 0; i < Nmc ; i++) {
        /* tg->nu = args->tgdata.nu; */
        /* tg->nui = args->tgdata.nui; */
        /* tg->Ei = args->tgdata.Ei; */
        //		tg->from = args->tgdata.from;
        //		tg->to = args->tgdata.to;
        tg++;
    }

    /*
       for (i=0; i< Nmc ;i++){
       printf ("nu %g \n", args->tgmc.tgmc[i].nu);
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

    tg = tgmc->tgmc;

    if (SAVEMC) {
        aux = fopen("debug_tangent_mc.dat", "w");
    }

    for (i = 0; i < Nmc; i++) {

        //create noise for h-F data
        generate_uncorrelated_noise(hnoise, tgunc->uh, rng);
        generate_uncorrelated_noise(Fnoise, tgunc->uF, rng);

        if (SAVEMC) {
            fprintf(aux, "#MC %d \n", i);
            fprintf(aux, "#noise\n");

            for (j = 0 ; j < hnoise->res; j++) {
                fprintf(aux, "%d %g %g \n", j, hnoise->data[j], Fnoise->data[j]);
            }

            fprintf(aux, "\n\n");
        }


        //add noise to original h-F data
        if (!sum_data_lines(hnew, h, hnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        if (!sum_data_lines(Fnew, F, Fnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        if (SAVEMC) {
            fprintf(aux, "#MC %d \n", i);
            fprintf(aux, "#new data\n");

            for (j = 0 ; j < hnew->res; j++) {
                fprintf(aux, "%d %g %g \n", j, hnew->data[j], Fnew->data[j]);
            }

            fprintf(aux, "\n\n");
        }

        if (SAVEMC2) {
            g_snprintf(auxstr, sizeof(auxstr), "debug%d.txt", i);
            aux2 = fopen(auxstr, "w");

            for (j = 0 ; j < hnew->res; j++) {
                fprintf(aux2, "%g %g \n", hnew->data[j], Fnew->data[j]);
            }

            fclose(aux2);
        }

        //add random noise to hmax, Fmax
        fddata->hmax = hmax_old + gaussian_random_number(rng) * tgunc->uh;
        fddata->Fmax = Fmax_old + gaussian_random_number(rng) * tgunc->uF;

        //run fit
        tangent_fit(hnew, Fnew, tg, fddata, instdata, area, tgunc->uh, tgunc->uF);
        tg++;
    }

    if (SAVEMC) {
        fclose(aux);
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
        //		printf(" %d Er %g Eit %g Hit %g \n", i, args->tgmc.tgmc[i].Er, args->tgmc.tgmc[i].Eit, args->tgmc.tgmc[i].Hit);
    }

    //calculate averages and stdevs
    //hp, Aphc, Hit, Er, Eit, S

    mc = tgmc;

    for (i = 0; i < NPARAM_TG; i++) {
        mc->mcavg[i] = 0;
        mc->mcstd[i] = 0;
    }

    tg = tgmc->tgmc;
    nn = 0;

    for (i = 0; i < Nmc; i++) {
        //check sanity
        if (tg->hc > 0) {
            mc->mcavg[0] += tg->hc;
            mc->mcavg[1] += tg->Aphc;
            mc->mcavg[2] += tg->Hit;
            mc->mcavg[3] += tg->Er;
            mc->mcavg[4] += tg->Eit;
            mc->mcavg[5] += tg->S;
            mc->mcstd[0] += tg->hc * tg->hc;
            mc->mcstd[1] += tg->Aphc * tg->Aphc;
            mc->mcstd[2] += tg->Hit * tg->Hit;
            mc->mcstd[3] += tg->Er * tg->Er;
            mc->mcstd[4] += tg->Eit * tg->Eit;
            mc->mcstd[5] += tg->S * tg->S;
            nn++;
        }

        tg++;
    }

    //use only reasonable data for average and standard deviation
    for (i = 0; i < NPARAM_TG; i++) {
        mc->mcavg[i] /= nn;
        mc->mcstd[i] /= nn;
        mc->mcstd[i] -= mc->mcavg[i] * mc->mcavg[i];
        mc->mcstd[i] = (mc->mcstd[i] < 0) ? 0 : mc->mcstd[i];
        mc->mcstd[i] = sqrt(mc->mcstd[i]);
    }

    // how many datasets were skipped
    tgmc->skipmc = Nmc - nn;

    /*
       printf("mean  stdev \n");
       for (i = 0; i< NPARAM_TG; i++){
       printf(" %g   %g \n", avg[i], std[i]);
       }
       */

    // Histograms
    //create data structure

    dd = tgmc->mcdd;
    tg = tgmc->tgmc;

    for (i = 0; i < Nmc; i++) {
        dd[0][i] = tg->hc;
        dd[1][i] = tg->Aphc;
        dd[2][i] = tg->Hit;
        dd[3][i] = tg->Er;
        dd[4][i] = tg->Eit;
        dd[5][i] = tg->S;
        tg++;
    }

    //prepare histogram structure

    for (i = 0; i < NPARAM_TG; i++) {
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

gboolean tangent_unc_has_results(TangentUncdata *tgunc)
{
    return (tgunc->Ec != NULL && tgunc->Hc != NULL);
}

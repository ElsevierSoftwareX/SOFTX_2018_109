#include "tool_hertz_unc.h"
#include "tool_hertz.h"

#include "datatypes.h"
#include "fddata.h"
#include "mc-utils.h"
#include "niget-common.h"

#include <math.h>

static void init_shifted_Hertzdata(const Hertzdata *src, Hertzdata *dest, gdouble dh);
static gboolean nparam_to_string(gint i, gchar *str, gint nstr);

void init_HertzUncdata(HertzUncdata *hzunc)
{
    hzunc->uh = 1;
    hzunc->uF = 0.001;
    /*
      hzunc->unu = 0.05;
      hzunc->unui = 0.04;
      hzunc->uEi = 20e9;
    */
    hzunc->uradius = 10e-9;
    hzunc->uEr = 1;
    hzunc->uEit = 1;
    hzunc->Nmc = N_MC;

    hzunc->uEruh = 0;
    hzunc->uEruF = 0;
    hzunc->uEruradius = 0;
    hzunc->uErtotal = 0;

    hzunc->uEituh = 0;
    hzunc->uEituF = 0;
    hzunc->uEituradius = 0;
    hzunc->uEitunu = 0;
    hzunc->uEitunui = 0;
    hzunc->uEituEi = 0;
    hzunc->uEittotal = 0;

    hzunc->uradiusuh = 0;
    hzunc->uradiusuF = 0;
    hzunc->uradiusunu = 0;
    hzunc->uradiusunui = 0;
    hzunc->uradiusuEi = 0;
    hzunc->uradiustotal = 0;

    hzunc->w = NULL;
    hzunc->ERc = NULL;
}

static void init_shifted_Hertzdata(const Hertzdata *src, Hertzdata *dest, gdouble dh)
{
    init_Hertzdata(dest);

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
}

void hertz_propagate_uncertainties_calc(const Hertzdata *hz, HertzUncdata *unc, const Instdata *inst)
{
    gdouble R, Er, Eit, nu, nui, Ei;
    gdouble uErunu, uErunui, uEruEi, uEruEit;
    gdouble dradiusdEr ;

    switch (hz->mode) {
    case R_MODE:
        //make copies to simplify code
        // input
        Ei = inst->Ei;
        nui = inst->nui;
        nu = inst->nu;
        R = hz->radius;
        //output
        Eit = hz->Eit;
        Er = hz->Er;

        unc->uEruh = hz->reg.unc_slope_x / hz->reg.slope * Er * unc->uh; //Er/a uah = Er/a ua_x*uh
        unc->uEruF = hz->reg.unc_slope_y / hz->reg.slope * Er * unc->uF; //Er/a uaF = Er/a ua_y*uF
        unc->uEruradius = 0.5 * Er / R * unc->uradius; // -1/2 *Er/R u(R)

        //uEituh = dEit/ dh uh = dEit / dEr * dEr/dh uh = dEit/dEr * uEruh
        //dEit/dEr = Eit^2/Er^2/ (1-nu^2)

        unc->uEituh = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruh;
        unc->uEituF = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruF;
        unc->uEituradius = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruradius;

        unc->uEitunu =  fabs(2 * nu / (1 - sq(nu)) * Eit * inst->unu); // -2nu/(1-nu^2) Eit u(nu)
        unc->uEitunui = fabs(2 * nui / (1 - sq(nu)) * sq(Eit) / Ei * inst->unui) * 1e9; // -2nui/(1-nu^2) Eit^2/Ei u(nui)
        unc->uEituEi =  fabs((1 - sq(nui)) / (1 - sq(nu)) * sq(Eit) / sq(Ei) * inst->uEi) * 1e9; //-(1-nui^2)/(1-nu)^2 Eit^2/Ei^2 u(Ei)

        // total uncertainties
        unc->uErtotal = sqrt(sq(unc->uEruh) + sq(unc->uEruF) + sq(unc->uEruradius));
        unc->uEittotal = sqrt(sq(unc->uEituh) + sq(unc->uEituF) + sq(unc->uEituradius) + sq(unc->uEitunu) + sq(unc->uEitunui) + sq(unc->uEituEi));
        break;

    case ER_MODE:
        //make copies to simplify code
        // input
        Er = hz->Er;
        // output
        R = hz->radius;

        if (verbose) {
            g_print("unc_slope %g %g \n", hz->reg.unc_slope_x, hz->reg.unc_slope_y);
        }

        unc->uradiusuh = R * 2 * hz->reg.unc_slope_x * unc->uh / hz->reg.slope; // 2 R/a uah = 2R ua_x*uh/a
        unc->uradiusuF = R * 2 * hz->reg.unc_slope_y * unc->uF / hz->reg.slope; // 2 R/a uaF = 2R ua_y*uF/a
        unc->uradiusuEr = 2 * R / Er * unc->uEr; // -2 *R/Er u(Er)

        // total uncertainties
        unc->uradiustotal = sqrt(unc->uradiusuh * unc->uradiusuh + unc->uradiusuF * unc->uradiusuF + unc->uradiusuEr * unc->uradiusuEr);

        break;

    case EIT_MODE :
        //make copies to simplify code
        // input
        Ei = inst->Ei;
        nui = inst->nui;
        nu = inst->nu;
        Eit = hz->Eit;
        //output
        R = hz->radius;
        Er = hz->Er;

        uErunu =  2 * nu * sq(Er) / Eit * inst->unu; // 2*nu/EIT Er^2 unu
        uErunui =  2 * nui * sq(Er) / Ei * inst->unui * 1e9; // 2*nui/Ei Er^2 unui
        uEruEi = (1 - sq(nui)) / sq(Ei) * sq(Er) * inst->uEi * 1e9; // (1-nui^2)/Ei^2 Er^2 uEi
        uEruEit = (1 - sq(nu)) / sq(Eit) * sq(Er) * unc->uEit ; // (1-nu^2)/EIT^2 Er^2 uEIT

        unc->uradiusuh = R * 2 * hz->reg.unc_slope_x * unc->uh / hz->reg.slope; // 2 R/a uah = 2R ua_x*uh/a
        unc->uradiusuF = R * 2 * hz->reg.unc_slope_y * unc->uF / hz->reg.slope; // 2 R/a uaF = 2R ua_y*uF/a
        dradiusdEr = 2 * R / Er; // dR/dEr = -2 *R/Er

        unc->uradiusunu = dradiusdEr * uErunu;
        unc->uradiusunui = dradiusdEr * uErunui;
        unc->uradiusuEi = dradiusdEr * uEruEi;
        unc->uradiusuEit = dradiusdEr * uEruEit;

        // total uncertainties
        unc->uradiustotal = sqrt(sq(unc->uradiusuh) + sq(unc->uradiusuF) + sq(unc->uradiusuEit) + sq(unc->uradiusuEi) + sq(unc->uradiusunu) + sq(unc->uradiusunui));

        break;
    }
}

gchar* hertz_uncertainty_export_data(const Hertzdata *hzdata, const HertzUncdata *hzunc,
				     const FDdata *fddata, const Instdata *instdata, enum ExportFormat exportformat)
{
    gint j;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Parameters of Hertz analysis \n");
    g_string_append_printf(buf, "# u(depth): %g nm\n", hzunc->uh);
    g_string_append_printf(buf, "# u(load): %g nm\n", hzunc->uF);
    g_string_append_printf(buf, "#fitted range:  [%g, %g ] nm \n", hzdata->range[0], hzdata->range[1]);

    switch (hzdata->mode) {
    case R_MODE:
        g_string_append_printf(buf, "#R_mode \n");
        g_string_append_printf(buf, "#tip radius: %g nm       u(tip radius): %g nm\n", hzdata->radius * 1e9, hzunc->uradius * 1e9);
        g_string_append_printf(buf, "#nu: %g                 u(nu):  %g  \n", instdata->nu, instdata->unu);
        g_string_append_printf(buf, "#nu_indenter: %g       u(nu_indenter):  %g   \n", instdata->nui, instdata->unui);
        g_string_append_printf(buf, "#E_indenter: %g        u(E_indenter):   %g GPa\n", instdata->Ei * 1e-9, instdata->uEi * 1e-9); //TODO units
        g_string_append_printf(buf, "#\n");

        g_string_append_printf(buf, "#Results of Hertz analysis \n");
        g_string_append_printf(buf, "#Er: %g GPa\n", hzdata->Er);
        g_string_append_printf(buf, "#Eit: %g GPa\n", hzdata->Eit);
        break;

    case ER_MODE:
        g_string_append_printf(buf, "#Er_mode \n");
        g_string_append_printf(buf, "#Er: %g GPa            u(Er) : %g GPa\n", hzdata->Er, hzunc->uEr);
        g_string_append_printf(buf, "#\n");

        g_string_append_printf(buf, "#Results of Hertz analysis \n");
        g_string_append_printf(buf, "#radius: %g nm\n", hzdata->radius * 1e9);
        break;

    case EIT_MODE:
        g_string_append_printf(buf, "#Eit_mode \n");
        g_string_append_printf(buf, "#Eit: %g GPa            u(Eit) : %g GPa\n", hzdata->Eit, hzunc->uEit);
        g_string_append_printf(buf, "#nu: %g                 u(nu):  %g  \n", instdata->nu, instdata->unu);
        g_string_append_printf(buf, "#nu_indenter: %g       u(nu_indenter):  %g   \n", instdata->nui, instdata->unui);
        g_string_append_printf(buf, "#E_indenter: %g        u(E_indenter):   %g GPa\n", instdata->Ei * 1e-9, instdata->uEi * 1e-9); //TODO units
        g_string_append_printf(buf, "#\n");

        g_string_append_printf(buf, "#Results of Hertz analysis \n");
        g_string_append_printf(buf, "#radius: %g nm\n", hzdata->radius * 1e9);
        break;
    }

    g_string_append_printf(buf, "#\n");

    /* Gaussian propagation */
    g_string_append_printf(buf, "#Results of the Gaussian propagation of uncertainties \n");

    switch (hzdata->mode) {
    case R_MODE:
        g_string_append_printf(buf, "#u(Er) =  %g GPa \n", hzunc->uErtotal);
        g_string_append_printf(buf, "#u(Eit) =  %g GPa \n", hzunc->uEittotal);
        g_string_append_printf(buf, "#\n");
        g_string_append_printf(buf, "#Contributions to uncertainty of Er\n");
        g_string_append_printf(buf, "#u(Er;uh) = %g GPa \n", hzunc->uEruh);
        g_string_append_printf(buf, "#u(Er;uF) = %g GPa \n", hzunc->uEruF);
        g_string_append_printf(buf, "#u(Er;uradius) = %g GPa \n", hzunc->uEruradius);
        g_string_append_printf(buf, "#\n");
        g_string_append_printf(buf, "#Contributions to uncertainty of Eit\n");
        g_string_append_printf(buf, "#u(Eit;uh) = %g GPa \n", hzunc->uEituh);
        g_string_append_printf(buf, "#u(Eit;uF) = %g GPa \n", hzunc->uEituF);
        g_string_append_printf(buf, "#u(Eit;uradius) = %g GPa \n", hzunc->uEituradius);
        g_string_append_printf(buf, "#u(Eit;unu) = %g GPa \n", hzunc->uEitunu);
        g_string_append_printf(buf, "#u(Eit;unui) = %g GPa \n", hzunc->uEitunui);
        g_string_append_printf(buf, "#u(Eit;uEi) = %g GPa \n", hzunc->uEituEi);
        g_string_append_printf(buf, "#\n");
        break;

    case ER_MODE:
        g_string_append_printf(buf, "#u(radius) =  %g nm \n", hzunc->uradiustotal * 1e9);
        g_string_append_printf(buf, "#\n");
        g_string_append_printf(buf, "#Contributions to uncertainty of radius\n");
        g_string_append_printf(buf, "#u(radius;uh) = %g nm \n", hzunc->uradiusuh * 1e9);
        g_string_append_printf(buf, "#u(radius;uF) = %g nm \n", hzunc->uradiusuF * 1e9);
        g_string_append_printf(buf, "#u(uradius;uEr) = %g nm \n", hzunc->uradiusuEr * 1e9);
        g_string_append_printf(buf, "#\n");
        break;

    case EIT_MODE:
        g_string_append_printf(buf, "#u(radius) =  %g nm \n", hzunc->uradiustotal * 1e9);
        g_string_append_printf(buf, "#\n");
        g_string_append_printf(buf, "#Contributions to uncertainty of radius\n");
        g_string_append_printf(buf, "#u(radius;uh) = %g nm \n", hzunc->uradiusuh * 1e9);
        g_string_append_printf(buf, "#u(radius;uF) = %g nm \n", hzunc->uradiusuF * 1e9);
        g_string_append_printf(buf, "#u(uradius;uEit) = %g nm \n", hzunc->uradiusuEit * 1e9);
        g_string_append_printf(buf, "#u(uradius;unu) = %g nm \n", hzunc->uradiusunu * 1e9);
        g_string_append_printf(buf, "#u(uradius;unui) = %g nm \n", hzunc->uradiusunui * 1e9);
        g_string_append_printf(buf, "#u(uradius;uEi) = %g nm \n", hzunc->uradiusuEi * 1e9);
        break;
    }

    g_string_append_printf(buf, "#\n");

    /* contact point uncertainties */
    switch (hzdata->mode) {
    case R_MODE:
        g_string_append_printf(buf, "#Changes of the Er corresponding to changes in the contact point \n");
        g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   ");
        g_string_append_printf(buf, "Er/GPa \n");

        for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
            if (hzunc->ERc[j + CONTACTMAX] < 0) {
                g_string_append_printf(buf, "# %d   not defined \n", j);
            }
            else {
                g_string_append_printf(buf, "# %d   %g  \n", j, hzunc->ERc[j + CONTACTMAX]);
            }
        }

        break;

    case ER_MODE:
    case EIT_MODE:
        g_string_append_printf(buf, "#Changes of the radius corresponding to changes in the contact point \n");
        g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   ");
        g_string_append_printf(buf, "R/nm \n");

        for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
            if (hzunc->ERc[j + CONTACTMAX] < 0) {
                g_string_append_printf(buf, "# %d   not defined \n", j);
            }
            else {
                g_string_append_printf(buf, "# %d   %g  \n", j, hzunc->ERc[j + CONTACTMAX] * 1e9);
            }
        }

        break;
    }

    return g_string_free(buf, FALSE);
}

void hertz_fit_shift_contact(const Hertzdata *hzdata, const FDdata *fddata, const Instdata *instdata, gdouble *ERc, gint j)
{
    Hertzdata hz;
    FDdata fdc;
    gdouble dh;
    gint istart, iend, ndata;
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

    // create working copy of Hertzdat
    init_shifted_Hertzdata(hzdata, &hz, dh);

    //check range

    if (verbose) {
        g_print(" hz.from %g to %g \n", hz.from, hz.to);
    }

    if (hz.from == hz.to) {
        return;
    }

    //find appropriate data, which should be fitted
    range_to_indices(hz.from, hz.to, fdc.hload, FALSE, &istart, &iend, &ndata);
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

    if (verbose) {
        g_print(" hz.mode %d \n", hz.mode);
    }

    hertz_fit(x, y, &hz, instdata);

    // clean up
    g_object_unref(x);
    g_object_unref(y);
    destroy_FDdata(&fdc);

    if (hz.mode == R_MODE) {
        ERc[j + CONTACTMAX] = hz.Er;
    }
    else {
        ERc[j + CONTACTMAX] = hz.radius;
    }

    if (verbose) {
        g_print("\n\n");
    }
}

void init_HertzMCdata(HertzMCdata *hzmc, gint Nmc)
{
    gint i;

    //number of iterations
    hzmc->Nmc = Nmc;

    if (verbose) {
        g_print("initialize for %d iterations \n", hzmc->Nmc);
    }

    //how much had to be skipped
    hzmc->skipmc = 0;

    //dynamically allocate
    hzmc->hertzmc = (Hertzdata *)g_malloc(hzmc->Nmc * sizeof(Hertzdata));
    hzmc->mcavg = (gdouble *)g_malloc(NPARAM_HZ * sizeof(gdouble));
    hzmc->mcstd = (gdouble *)g_malloc(NPARAM_HZ * sizeof(gdouble));
    hzmc->mcdd = (gdouble **) g_malloc(NPARAM_HZ * sizeof(gdouble *));
    hzmc->mcxhist = (gdouble **) g_malloc(NPARAM_HZ * sizeof(gdouble *));
    hzmc->mcyhist = (gdouble **) g_malloc(NPARAM_HZ * sizeof(gdouble *));

    hzmc->mcnstat = (gint)floor(3.49 * cbrt(hzmc->Nmc) + 0.5);

    for (i = 0; i < NPARAM_HZ; i++) {
        hzmc->mcdd[i] = (gdouble *) g_malloc(hzmc->Nmc * sizeof(gdouble));
        hzmc->mcxhist[i] = (gdouble *) g_malloc(hzmc->mcnstat * sizeof(gdouble));
        hzmc->mcyhist[i] = (gdouble *) g_malloc(hzmc->mcnstat * sizeof(gdouble));
    }
}

void destroy_HertzMCdata(HertzMCdata *hzmc)
{
    gint i;

    //free everything
    g_free(hzmc->mcavg);
    g_free(hzmc->mcstd);

    for (i = 0; i < NPARAM_HZ; i++) {
        g_free(hzmc->mcdd[i]);
        g_free(hzmc->mcxhist[i]);
        g_free(hzmc->mcyhist[i]);
    }

    g_free(hzmc->mcdd);
    g_free(hzmc->mcxhist);
    g_free(hzmc->mcyhist);

    hzmc->mcnstat = 0;
    hzmc->skipmc = 0;
    g_free(hzmc->hertzmc);
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

gchar* hertz_mc_export_data(const Hertzdata *hertzdata, const HertzUncdata *hertzunc, const HertzMCdata *hertzmc,
			    const FDdata *fddata, const Instdata *instdata, enum ExportFormat exportformat)
{
    Hertzdata *hz;
    gint i, j;
    gchar str[200];
    gboolean is; /* RS - this is not used for anything, can be removed? */
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Monte Carlo simulation of Hertz analysis\n");
    g_string_append_printf(buf, "#Input uncertainties \n");
    g_string_append_printf(buf, "#u(h) = %g nm \n", hertzunc->uh);
    g_string_append_printf(buf, "#u(F) = %g mN \n", hertzunc->uF);
    g_string_append_printf(buf, "#Number of iterations  = %d  \n", hertzmc->Nmc);
    g_string_append_printf(buf, "\n");

    g_string_append_printf(buf, "#Parameters of Hertz analysis \n");

    switch (hertzdata->mode) {
    case R_MODE:
        g_string_append_printf(buf, "#tip radius: %g nm\n", hertzdata->radius * 1e9);
        g_string_append_printf(buf, "#nu: %g                 \n", instdata->nu);
        g_string_append_printf(buf, "#nu_indenter: %g        \n", instdata->nui);
        g_string_append_printf(buf, "#E_indenter: %g      GPa\n", instdata->Ei * 1e-9); //TODO units
        g_string_append_printf(buf, "#fitted range:  [%g, %g ] nm \n", hertzdata->range[0], hertzdata->range[1]);
        break;

    case ER_MODE:

        g_string_append_printf(buf, "#Er: %g      GPa\n", hertzdata->Er);  //TODO units
        g_string_append_printf(buf, "#fitted range:  [%g, %g ] nm \n", hertzdata->range[0], hertzdata->range[1]);
        break;

    case EIT_MODE:
        g_string_append_printf(buf, "#Eit: %g      GPa\n", hertzdata->Eit);  //TODO units
        g_string_append_printf(buf, "#nu: %g                 \n", instdata->nu);
        g_string_append_printf(buf, "#nu_indenter: %g        \n", instdata->nui);
        g_string_append_printf(buf, "#E_indenter: %g      GPa\n", instdata->Ei * 1e-9); //TODO units
        g_string_append_printf(buf, "#fitted range:  [%g, %g ] nm \n", hertzdata->range[0], hertzdata->range[1]);
    }

    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n\n");

    g_string_append_printf(buf, "#Results \n");

    switch (hertzdata->mode) {
    case R_MODE:
        for (i = 0; i < 2; i++) {
            is = nparam_to_string(i, str, sizeof(str));
            g_string_append_printf(buf, "#%s = (%g ± %g ) \n", str, hertzmc->mcavg[i], hertzmc->mcstd[i]);
        }

        break;

    case ER_MODE:
    case EIT_MODE:
        is = nparam_to_string(2, str, sizeof(str));
        g_string_append_printf(buf, "#%s = (%g ± %g ) \n", str, hertzmc->mcavg[2], hertzmc->mcstd[2]);
        break;
    }

    g_string_append_printf(buf, "\n\n\n");

    // Histograms

    if (hertzmc->mcnstat > 0) {
        g_string_append_printf(buf, "#Histogram\n");

        switch (hertzdata->mode) {
        case R_MODE:
            for (i = 0; i < 2; i++) {
                g_string_append_printf(buf, "#Histogram\n");
                is = nparam_to_string(i, str, sizeof(str));

                if (is) {
                    g_string_append_printf(buf, "#%s    count/a.u.  \n", str);

                    for (j = 0; j < hertzmc->mcnstat; j++) {
                        g_string_append_printf(buf, "%g %g \n", hertzmc->mcxhist[i][j], hertzmc->mcyhist[i][j]);
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

                for (j = 0; j < hertzmc->mcnstat; j++) {
                    g_string_append_printf(buf, "%g %g \n", hertzmc->mcxhist[i][j], hertzmc->mcyhist[i][j]);
                }

                g_string_append_printf(buf, "\n \n");
            }
        }
    }

    g_string_append_printf(buf, "\n\n\n");

    //Full data
    g_string_append_printf(buf, "#Full data \n");
    g_string_append_printf(buf, "#");

    switch (hertzdata->mode) {
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

    hz = hertzmc->hertzmc;

    if (hertzdata->mode == R_MODE) {
        for (i = 0; i < hertzmc->Nmc; i++) {
            g_string_append_printf(buf, MODUL, hz->Er);
            g_string_append_printf(buf, "  ");
            g_string_append_printf(buf, MODUL, hz->Eit);
            g_string_append_printf(buf, "  \n");
            hz++;
        }
    }
    else {
        for (i = 0; i < hertzmc->Nmc; i++) {
            g_string_append_printf(buf, DEPTH, hz->radius * 1e9);
            g_string_append_printf(buf, "  \n");
            hz++;
        }
    }

    return g_string_free(buf, FALSE);
}

void hertz_uncertainty_montecarlo_run_calc(const Hertzdata *hertzdata, const HertzUncdata *hertzunc, HertzMCdata *hertzmc,
					   const FDdata *fddata, const Instdata *instdata)
{
    gint istart, iend, ndata;
    gint i, Nmc;
    Hertzdata *hz;
    HertzMCdata *mc;
    GwyDataLine *h, *F;
    GwyDataLine *hnew, *Fnew;
    GwyDataLine *hnoise, *Fnoise;
    GRand *rng;
    gint nn;
    gdouble **dd;
    gdouble *xhist, *yhist;
    gboolean hist;

    //not necessary to check range, couldn't get here in the first place

    //find appropriate data, which should be fitted, must use real fitted range
    range_to_indices(hertzdata->range[0], hertzdata->range[1], fddata->hload, FALSE, &istart, &iend, &ndata);
    // if there are not enough data in the selected range, not fitting can be performed
    ndata = MIN(ndata, fddata->hload->res - istart);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    h = gwy_data_line_part_extract(fddata->hload, istart, ndata);
    F = gwy_data_line_part_extract(fddata->Fload, istart, ndata);

    Nmc = hertzmc->Nmc;
    hz = hertzmc->hertzmc;

    //prepare data
    for (i = 0; i < Nmc ; i++) {
        init_Hertzdata(hz);
        hz++;
    }

    if (verbose) {
        g_print("copy params \n");
    }

    hz = hertzmc->hertzmc;

    //copy input parameters
    for (i = 0; i < Nmc ; i++) {
        hz->mode = hertzdata->mode;

        switch (hz->mode) {
        case R_MODE:
            hz->radius = hertzdata->radius;
            break;

        case ER_MODE:
            hz->Er = hertzdata->Er;
            break;

        case EIT_MODE:
            hz->Eit = hertzdata->Eit;
            hz->Er = calc_Er(hz->Eit, instdata->nu, instdata->Ei * 1e-9, instdata->nui);
            break;
        }

        hz++;
    }

    // initialize random generator and datafield
    rng = g_rand_new();

    //prepare datafield for noise
    hnoise = gwy_data_line_new_alike(h, FALSE);
    Fnoise = gwy_data_line_new_alike(h, FALSE);

    //prepare datafield for varied fields
    hnew = gwy_data_line_new_alike(h, FALSE);
    Fnew = gwy_data_line_new_alike(F, FALSE);

    hz = hertzmc->hertzmc;

    for (i = 0; i < Nmc; i++) {

        //create noise for h-F data
        generate_uncorrelated_noise(hnoise, hertzunc->uh, rng);
        generate_uncorrelated_noise(Fnoise, hertzunc->uF, rng);

        //add noise to original h-F data
        if (!sum_data_lines(hnew, h, hnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        if (!sum_data_lines(Fnew, F, Fnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        //run fit
        hertz_fit(hnew, Fnew, hz, instdata);
        /*
        		printf("-------------------\n");
        		gint j;
        		for (j=0; j< hnew->res;j++)
        			printf("%g %g \n", hnew->data[j], Fnew->data[j]);
        */
        hz++;
    }

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

    //calculate averages and stdevs
    // Er, Eit

    mc = hertzmc;

    for (i = 0; i < NPARAM_HZ; i++) {
        mc->mcavg[i] = 0;
        mc->mcstd[i] = 0;
    }

    hz = hertzmc->hertzmc;
    nn = 0;

    for (i = 0; i < Nmc; i++) {
        //sanity check
        {
            mc->mcavg[0] += hz->Er;
            mc->mcavg[1] += hz->Eit;
            mc->mcavg[2] += hz->radius * 1e9; //TODO UNITS
            mc->mcstd[0] += hz->Er * hz->Er;
            mc->mcstd[1] += hz->Eit * hz->Eit;
            mc->mcstd[2] += hz->radius * 1e9 * hz->radius * 1e9; //TODO UNITS
            nn++;
        }
        hz++;
    }

    //use only reasonable data for average and standard deviation
    for (i = 0; i < NPARAM_HZ; i++) {
        mc->mcavg[i] /= nn;
        mc->mcstd[i] /= nn;
        mc->mcstd[i] -= mc->mcavg[i] * mc->mcavg[i];
        mc->mcstd[i] = (mc->mcstd[i] < 0) ? 0 : mc->mcstd[i];
        mc->mcstd[i] = sqrt(mc->mcstd[i]);
    }

    // how many datasets were skipped
    hertzmc->skipmc = Nmc - nn;

    // Histograms
    //create data structure

    dd = hertzmc->mcdd;
    hz = hertzmc->hertzmc;

    for (i = 0; i < Nmc; i++) {
        dd[0][i] = hz->Er;
        dd[1][i] = hz->Eit;
        dd[2][i] = hz->radius * 1e9;
        hz++;
    }

    //prepare histogram structure

    for (i = 0; i < NPARAM_HZ; i++) {
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

gboolean hertz_unc_has_results(const HertzUncdata *hzunc)
{
    return (hzunc->ERc != NULL);
}

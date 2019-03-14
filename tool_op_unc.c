#include "tool_op_unc.h"
#include "tool_op.h"

#include "datatypes.h"
#include "fddata.h"
#include "fit-utils.h"
#include "niget-common.h"
#include "mc-utils.h"

#include <math.h>

#include <libprocess/gwyprocess.h>

static void init_shifted_OPdata(const OPdata *src, OPdata *dest, gdouble dh);
static gboolean nparam_to_string(gint i, gchar *str, gint nstr);

void init_OPUncdata(OPUncdata *opunc)
{
    opunc->uh = 1;
    opunc->uF = 0.001;
    /* 
       opunc->unu = 0.05;
       opunc->unui = 0.04;
       opunc->uEi = 20e9;
    */
    opunc->Nmc = N_MC;

    opunc->uhcuh = 0;
    opunc->uhcuF = 0;
    opunc->uhc = 0;

    opunc->uAuh = 0;
    opunc->uAuF = 0;
    opunc->uA = 0;

    opunc->uHituh = 0;
    opunc->uHituF = 0;
    opunc->uHit = 0;

    opunc->uEruh = 0;
    opunc->uEruF = 0;
    opunc->uEr = 0;

    opunc->uEituh = 0;
    opunc->uEituF = 0;
    opunc->uEitunu = 0;
    opunc->uEitunui = 0;
    opunc->uEituEi = 0;
    opunc->uEit = 0;

    opunc->Ec = NULL;
    opunc->Hc = NULL;
}

static void init_shifted_OPdata(const OPdata *src, OPdata *dest, gdouble dh)
{
    init_OPdata(dest);

    // shift fitting range  in h-mode by dh
    // don't shift range for F-mode
    // both for hp fit and for S fit
    dest->hprange[0] = src->hprange[0] - dh;
    dest->hprange[1] = src->hprange[1] - dh;
    dest->hpfrom = dest->hprange[0];
    dest->hpto = dest->hprange[1];

    dest->hprange_pct_Fmax[0] = src->hprange_pct_Fmax[0];
    dest->hprange_pct_Fmax[1] = src->hprange_pct_Fmax[1];
    dest->hpfrom_pct_Fmax = src->hpfrom_pct_Fmax;
    dest->hpto_pct_Fmax = src->hpto_pct_Fmax;

    dest->Finputhp = src->Finputhp;

    dest->Srange[0] = src->Srange[0] - dh;
    dest->Srange[1] = src->Srange[1] - dh;
    dest->Sfrom = dest->Srange[0];
    dest->Sto = dest->Srange[1];

    dest->Srange_pct_Fmax[0] = src->Srange_pct_Fmax[0];
    dest->Srange_pct_Fmax[1] = src->Srange_pct_Fmax[1];
    dest->Sfrom_pct_Fmax = src->Sfrom_pct_Fmax;
    dest->Sto_pct_Fmax = src->Sto_pct_Fmax;

    dest->FinputS = src->FinputS;

    // copy beta
    dest->beta = src->beta;
}

/* do we need to modify OPdata LinRegs? */
void op_propagate_uncertainties_calc(const OPdata *op, OPUncdata *unc, const FDdata *fd, const Instdata *inst, const Area *area)
{
    // TODO
    gdouble uh, uF;
    gdouble uSh, uSF ;
    gdouble uhph, uhpF;
    gdouble uah, uaF, ubh, ubF;
    gdouble umh, umF,  depsdm, uepsh, uepsF;
    gint istart, iend, ndata;
    gint i;
    gdouble covabh, covabF;
    gdouble covShch, covShcF;
    gdouble covmhph, covmhpF;
    gdouble dAdh;
    gdouble m, S, eps, hc, hp, hmax, Fmax, Aphc, Hit, Er, Eit, nu, nui, Ei;
    GwyDataLine *x, *y, *ux, *uy;

    //make copies to simplify code
    uh = unc->uh;
    uF = unc->uF;
    m = op->m;
    eps = op->eps;
    S = op->S;
    hp = op->hp;
    hc = op->hc;
    hmax = fd->hmax;
    Fmax = fd->Fmax;
    Aphc = op->Aphc;
    Hit = op->Hit;
    Eit = op->Eit;
    Er = op->Er;
    Ei = inst->Ei;
    nui = inst->nui;
    nu = inst->nu;

    if (verbose) {
        g_print(" \n\nPropagate uncertainties\n");
    }

    /* constant uncertainty for all data points */
    if (op->Finputhp) {
        range_to_indices(op->hprange_pct_Fmax[0] * 0.01 * fd->Fmax, op->hprange_pct_Fmax[1] * 0.01 * fd->Fmax, fd->Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(op->hprange[0], op->hprange[1], fd->hunload, TRUE, &istart, &iend, &ndata);
    }

    x = gwy_data_line_part_extract(fd->hunload, istart, ndata);
    y = gwy_data_line_part_extract(fd->Funload, istart, ndata);

    total_least_squares_line_fit_const(x, y, &(op->reghp), uh, uF);
    g_object_unref(x);
    g_object_unref(y);

    //for hp fit a*x+b
    uah = op->reghp.unc_slope_x;
    uaF = op->reghp.unc_slope_y;
    ubh = op->reghp.unc_intercept_x;
    ubF = op->reghp.unc_intercept_y;
    covabh = op->reghp.cov_slope_intercept_x;
    covabF = op->reghp.cov_slope_intercept_y;

    if (verbose)
        g_print(" a %g b %g hp %g  uah %g uaF %g ubh %g ubF %g \n",
                op->reghp.slope, op->reghp.intercept, -op->reghp.intercept / op->reghp.slope, uah, uaF, ubh, ubF);

    /*   u(hp)   */
    uhph = hp * sqrt((sq(uah) / sq(op->reghp.slope)) + (sq(ubh) / sq(op->reghp.intercept)) - 2 * covabh);
    uhpF = hp * sqrt((sq(uaF) / sq(op->reghp.slope)) + (sq(ubF) / sq(op->reghp.intercept)) - 2 * covabF);

    /* hp = -b/a
     * uhp^2 = ( -ub/a)^2 + (ua*b/a^2) ^2 + 2 cov(a,b) 1/a^2 = hp^2 * ( (ub/b)^2 + ( ua/a)^2 ) -2 corr(a,b) ua ub hp/(ab)
     */

    if (verbose) {
        g_print("u(hp;h) %g u(hp;F) %g      urel(hp;h) %g urel(hp;F) %g \n", uhph, uhpF, uhph / hp, uhpF / hp);
    }

    /*   u(m)    */

    if (op->FinputS) {
        range_to_indices(op->Srange_pct_Fmax[0] * 0.01 * fd->Fmax, op->Srange_pct_Fmax[1] * 0.01 * fd->Fmax, fd->Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(op->Srange[0], op->Srange[1], fd->hunload, TRUE, &istart, &iend, &ndata);
    }

    x = gwy_data_line_new(ndata, 1.0, TRUE);
    y = gwy_data_line_new(ndata, 1.0, TRUE);
    ux = gwy_data_line_new(ndata, 1.0, TRUE);
    uy = gwy_data_line_new(ndata, 1.0, TRUE);

    for (i = 0; i < ndata; i++) {
        x->data[i] = log(fd->hunload->data[istart + i] - hp);
        y->data[i] = log(fd->Funload->data[istart + i]);
        ux->data[i] = sqrt(sq(uh) + sq(uhph)) / (fd->hunload->data[istart + i] - hp); //TODO chybi u(hp;F)
        uy->data[i] = uF / fd->Funload->data[istart + i];
    }

    total_least_squares_line_fit_nonconst(x, y, &(op->regS), ux, uy);
    umh = op->regS.unc_slope_x;
    umF = op->regS.unc_slope_y;

    covmhpF = 0;
    covmhph = 0;

    for (i = 0; i < x->res; i ++) {
        covmhph -= op->regS.auxtoberemoved_axi[i] / exp(x->data[i]);
    }

    if (verbose) {
        g_print(" covmhph %g \n", covmhph);
    }

    g_object_unref(x);
    g_object_unref(y);
    g_object_unref(ux);
    g_object_unref(uy);

    if (verbose) {
        g_print(" new  slope %g intercept %g slope_x %g slope_y %g intercept_x %g intercept_y %g cov_x %g cov_y %g  \n",
                op->regS.slope, op->regS.intercept,
                op->regS.unc_slope_x, op->regS.unc_slope_y,
                op->regS.unc_intercept_x, op->regS.unc_intercept_y,
                op->regS.cov_slope_intercept_x, op->regS.cov_slope_intercept_y);
        g_print(" m %g umh %g umF %g \n", m, umh, umF);
    }

    /*    u(eps)   */
    depsdm = der_epsilon(m);
    uepsh = depsdm * umh;
    uepsF = depsdm * umF;

    if (verbose) {
        g_print("depsdm %g  \n", depsdm);
        g_print(" u(eps;h) %g u(eps;F) %g  urel(eps;h) %g urel(eps;F) %g \n", uepsh, uepsF, uepsh / eps, uepsF / eps);
    }

    /*    u(S)   */
    uSh = sqrt(sq(umh) / sq(m) + sq(uh) / sq(hmax - hp) + sq(uhph) / sq(hmax - hp) + 2 * covmhph / m / (hmax - hp)) * S;
    uSF = sqrt(sq(umF) / sq(m) + sq(uF) / sq(Fmax) + sq(uhpF) / sq(hmax - hp) + 2 * covmhpF / m / (hmax - hp)) * S;
    unc->uSuh = uSh;
    unc->uSuF = uSF;
    /* u(S;h) = S *sqrt [ (umh/m)^2 + (uh/(hmax-hp))^2 + (uhph/(hmax-hp))^2  + 2 cov(m,hp;h)/(m*(hmax-hp))]
     * u(S;F) = S *sqrt [ (umF/m)^2 + (uFmax/Fmax)^2   + (uhpF/(hmax-hp))^2  + 2 cov(m,hp;F)/(m*(hmax-hp))]
     */

    if (verbose) {
        g_print(" u(S;h) %g u(S;F) %g  urel(S;h) %g urel(S;F) %g \n", uSh, uSF, uSh / S, uSF / S);
    }

    /*   u(hc)   */
    unc->uhcuh = sqrt(sq(hmax - hp) * sq(eps / m - depsdm) / sq(m) * sq(umh) + sq(1 - eps / m) * sq(uh) + sq(eps) / sq(m) * sq(uhph));
    unc->uhcuF = sqrt(sq(hmax - hp) * sq(eps / m - depsdm) / sq(m) * sq(umF) + sq(eps) / sq(m) * sq(uhpF));

    /* u(hc;h) ^2 = [(hmax-hp) (eps/m- deps/dm)/m  u(m;h) ]^2 + [(1-eps/m) u(hmax)]^2 + [eps/m u(hp;h)]^2
     u(hc;F) ^2 = [(hmax-hp) (eps/m- deps/dm)/m  u(m;F) ]^2 +                         [eps/m u(hp;F)]^2
     */

    if (verbose) {
        g_print(" u(hc;h) %g u(hc;F) %g  urel(hc;h) %g urel(hc;F) %g \n", unc->uhcuh, unc->uhcuF, unc->uhcuh / hc, unc->uhcuF / hc);
    }

    /*    u(A)   */
    dAdh =  eval_area_derivative(area, hc) ;
    unc->uAuh = dAdh * unc->uhcuh;
    unc->uAuF = dAdh * unc->uhcuF;

    if (verbose) {
        g_print(" u(A;h) %g u(A;F) %g  urel(A;h) %g urel(A;F) %g \n", unc->uAuh, unc->uAuF, unc->uAuh / Aphc, unc->uAuF / Aphc);
    }

    /*     u(H_IT)   */
    unc->uHituh = Hit / Aphc * dAdh * unc->uhcuh;
    unc->uHituF = Hit * sqrt(sq(uF) / sq(Fmax) + sq(dAdh) / sq(Aphc) * sq(unc->uhcuh));

    if (verbose) {
        g_print(" u(Hit;h) %g u(Hit;F) %g  urel(Hit;h) %g urel(Hit;F) %g \n", unc->uHituh, unc->uHituF, unc->uHituh / Hit, unc->uHituF / Hit);
    }

    /*     u(E_r)   */

    covShch = - S * (1 - eps) / (hmax - hp) * sq(uh) - S * (hmax - hp) / m * depsdm * sq(umh) +	S * eps / (hmax - hp) * sq(uhph);
    covShcF =  - S * (hmax - hp) / m * depsdm * sq(umF) + S * eps / (hmax - hp) * sq(uhpF);

    /*cov(S,hc;h) = -S(1-eps)/(hmax-hp)*uh^2 -S(hmax-hp)/m *(deps/dm) *u(m;h)^2 + S*eps/(hmax-hp)*u(hp;h)^2
     *cov(S,hc;F) =                          -S(hmax-hp)/m *(deps/dm) *u(m;F)^2 + S*eps/(hmax-hp)*u(hp;F)^2
     */

    unc->uEruh = Er * sqrt(sq(uSh) / sq(S) + sq(unc->uAuh) / sq(Aphc) / 4 - 1 / S / Aphc * dAdh * covShch);
    unc->uEruF = Er * sqrt(sq(uSF) / sq(S) + sq(unc->uAuF) / sq(Aphc) / 4 - 1 / S / Aphc * dAdh * covShcF);
    /* u(Er;h) = sqrt[ (Er/S)^2 u(S;h)^2 + (Er/2A)^2 u(A;h)^2  - Er^2/(S Aphc)* (dA/dh) cov(S,hc;h) ]
     * u(Er;F) = sqrt[ (Er/S)^2 u(S;F)^2 + (Er/2A)^2 u(A;F)^2  - Er^2/(S Aphc)* (dA/dh) cov(S,hc;F) ]
     */

    if (verbose) {
        g_print(" u(Er;h) %g u(Er;F) %g  urel(Er;h) %g urel(Er;F) %g \n", unc->uEruh, unc->uEruF, unc->uEruh / Er, unc->uEruF / Er);
    }

    /*     u(E_IT)   */
    unc->uEituh = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruh; // Eit^2/Er^2 /(1-nu^2) u(Er;h)
    unc->uEituF = sq(Eit) / sq(Er) / (1 - sq(nu)) * unc->uEruF; // Eit^2/Er^2 /(1-nu^2) u(Er;F)

    unc->uEitunu =  fabs(2 * nu / (1 - sq(nu)) * Eit * inst->unu); // -2nu/(1-nu^2) Eit u(nu)
    unc->uEitunui = fabs(2 * nui / (1 - sq(nu)) * sq(Eit) / Ei * inst->unui) * 1e9; // -2nui/(1-nu^2) Eit^2/Ei u(nui)
    unc->uEituEi =  fabs((1 - sq(nui)) / (1 - sq(nu)) * sq(Eit) / sq(Ei) * inst->uEi) * 1e9; //-(1-nui^2)/(1-nu)^2 Eit^2/Ei^2 u(Ei)

    if (verbose) {
        g_print("u(Eit;h) %g  u(Eit;F) %g  urel(Eit;h) %g  urel(Eit;F) %g\n",
                unc->uEituh, unc->uEituF, unc->uEituh / Eit, unc->uEituF / Eit);
        g_print("u(Eit;nu) %g  u(Eit;nui) %g  u(Eit;Ei) %g  urel(Eit;nu) %g  urel(Eit;nui) %g  urel(Eit;Ei) %g\n",
                unc->uEitunu, unc->uEitunui, unc->uEituEi, unc->uEitunu / Eit, unc->uEitunui / Eit, unc->uEituEi / Eit);
    }

    if (verbose) {
        g_print(" \nEnd propagation of uncertainties\n\n");
    }
}

gchar* op_uncertainty_export_data(const OPdata *opdata, const OPUncdata *opunc, const FDdata *fddata, const Instdata *instdata,
				 enum ExportFormat exportformat)
{
    gint j;
    GString *buf;
    
    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Parameters of Oliver Pharr analysis  \n");
    g_string_append_printf(buf, "# u(depth): %g nm\n", opunc->uh);
    g_string_append_printf(buf, "# u(load): %g nm\n", opunc->uF);
    g_string_append_printf(buf, "#beta: %g \n", opdata->beta);
    g_string_append_printf(buf, "#nu: %g                 u(nu):  %g  \n", instdata->nu, instdata->unu);
    g_string_append_printf(buf, "#nu_indenter: %g       u(nu_indenter):  %g   \n", instdata->nui, instdata->unui);
    g_string_append_printf(buf, "#E_indenter: %g        u(E_indenter):   %g GPa\n", instdata->Ei * 1e-9, instdata->uEi * 1e-9); //TODO units
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
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Results of Oliver Pharr analysis \n");
    g_string_append_printf(buf, "#S = %g mN/nm\n", opdata->S);
    g_string_append_printf(buf, "#hc = %g nm\n", opdata->hc);
    g_string_append_printf(buf, "#Ap(hc) = %g nm^2\n", opdata->Aphc);
    g_string_append_printf(buf, "#Er = %g GPa\n", opdata->Er);
    g_string_append_printf(buf, "#Hit = %g GPa\n", opdata->Hit);
    g_string_append_printf(buf, "#Eit = %g GPa\n", opdata->Eit);
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Results of the Gaussian propagation of uncertainties \n");
    g_string_append_printf(buf, "#u(S) =  %g nm \n", opunc->uS);
    g_string_append_printf(buf, "#u(hc) =  %g nm \n", opunc->uhc);
    g_string_append_printf(buf, "#u(Ap(hc)) =  %g nm^2 \n", opunc->uA);
    g_string_append_printf(buf, "#u(Hit) =  %g MPa \n", opunc->uHit);
    g_string_append_printf(buf, "#u(Er) =  %g GPa \n", opunc->uEr);
    g_string_append_printf(buf, "#u(Eit) =  %g GPa \n", opunc->uEit);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of S\n");
    g_string_append_printf(buf, "#u(S;uh) = %g mN/nm \n", opunc->uSuh);
    g_string_append_printf(buf, "#u(S;uF) = %g mN/nm \n", opunc->uSuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of hc\n");
    g_string_append_printf(buf, "#u(hc;uh) = %g nm \n", opunc->uhcuh);
    g_string_append_printf(buf, "#u(hc;uF) = %g nm \n", opunc->uhcuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Ap(hc)\n");
    g_string_append_printf(buf, "#u(Ap(hc);uh) = %g nm^2 \n", opunc->uAuh);
    g_string_append_printf(buf, "#u(Ap(hc);uF) = %g nm^2 \n", opunc->uAuF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Hit\n");
    g_string_append_printf(buf, "#u(Hit;uh) = %g MPa \n", opunc->uHituh);
    g_string_append_printf(buf, "#u(Hit;uF) = %g MPa \n", opunc->uHituF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Er\n");
    g_string_append_printf(buf, "#u(Er;uh) = %g GPa \n", opunc->uEruh);
    g_string_append_printf(buf, "#u(Er;uF) = %g GPa \n", opunc->uEruF);
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Contributions to uncertainty of Eit\n");
    g_string_append_printf(buf, "#u(Eit;uh) = %g GPa \n", opunc->uEituh);
    g_string_append_printf(buf, "#u(Eit;uF) = %g GPa \n", opunc->uEituF);
    g_string_append_printf(buf, "#u(Eit;unu) = %g GPa \n", opunc->uEitunu);
    g_string_append_printf(buf, "#u(Eit;unui) = %g GPa \n", opunc->uEitunui);
    g_string_append_printf(buf, "#u(Eit;uEi) = %g GPa \n", opunc->uEituEi);
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#\n");

    g_string_append_printf(buf, "#Changes of the S corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#S/(mN/nm) \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (opunc->Sc[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, opunc->Sc[j + CONTACTMAX]);
        }
    }

    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#Changes of the Er corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#Er/GPa \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (opunc->Ec[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, opunc->Ec[j + CONTACTMAX]);
        }
    }

    g_string_append_printf(buf, "#\n");
    g_string_append_printf(buf, "#Changes of the H corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   \n");
    g_string_append_printf(buf, "#H/MPa \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        if (opunc->Hc[j + CONTACTMAX] < 0) {
            g_string_append_printf(buf, "# %d   not defined \n", j);
        }
        else {
            g_string_append_printf(buf, "# %d   %g  \n", j, opunc->Hc[j + CONTACTMAX]);
        }
    }

    return g_string_free(buf, FALSE);
}

void op_fit_shift_contact(const OPdata *opdata, const FDdata *fddata, const Instdata *instdata, const Area *area,
			  gdouble *Sc, gdouble *Ec, gdouble *Hc, gint j)
{
    OPdata op;
    FDdata fdc;
    gdouble dh;
    gint istart, iend, ndata;
    GwyDataLine *x, *y;
    gint i;

    Sc[j + CONTACTMAX] = -1;
    Ec[j + CONTACTMAX] = -1;
    Hc[j + CONTACTMAX] = -1;

    if (fddata->i_contact_load + j < 0) {
        return;
    }

    // create shifted FDdata
    dh = create_shifted_FDdata(fddata, &fdc, j);

    if (verbose) {
        g_print("\n\nCreate shifted data for j %d \n", j);
    }

    // create working copy of OPdata
    init_shifted_OPdata(opdata, &op, dh);
    // shift fitting range both for hp fit and for S fit
    /*
    op = *opdata;
    op.hprange[0] += -dh;
    op.hprange[1] += -dh;
    op.hpfrom = op.hprange[0];
    op.hpto = op.hprange[1];
    op.Srange[0] += -dh;
    op.Srange[1] += -dh;
    op.Sfrom = op.Srange[0];
    op.Sto = op.Srange[1];
    */

    /* hp fit */

    /* check range */
    if (opdata->Finputhp) {
        if (verbose) {
            g_print(" hpfrom %g hpto %g \n", op.hpfrom_pct_Fmax, op.hpto_pct_Fmax);
        }

        if (op.hpfrom_pct_Fmax == op.hpto_pct_Fmax) {
            return;
        }
    }
    else {
        if (verbose) {
            g_print(" hpfrom %g hpto %g \n", op.hpfrom, op.hpto);
        }

        if (op.hpfrom == op.hpto) {
            return;
        }
    }

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (op.Finputhp)
        range_to_indices(op.hprange_pct_Fmax[0] * 0.01 * fdc.Fmax, op.hprange_pct_Fmax[1] * 0.01 * fdc.Fmax,
                         fdc.Funload, TRUE, &istart, &iend, &ndata);
    else {
        range_to_indices(op.hprange[0], op.hprange[1], fdc.hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, no fitting can be performed
    ndata = MIN(ndata, fdc.hunload->res - istart);

    if (verbose)
        g_print("interval has %d points, from %d to %d, hrange [%g, %g ] nm, Frange [%g, %g] mN \n",
                ndata, istart, iend,
                fdc.hunload->data[iend], fdc.hunload->data[istart], fdc.Funload->data[iend], fdc.Funload->data[istart]);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    x = gwy_data_line_part_extract(fdc.hunload, istart, ndata);
    y = gwy_data_line_part_extract(fdc.Funload, istart, ndata);

    // perform fit and get results
    op_fit_hp(x, y, &op, &fdc);

    // clean up
    g_object_unref(x);
    g_object_unref(y);

    /*  S fit */

    /*check range */
    if (opdata->FinputS) {
        if (verbose) {
            g_print(" Sfrom %g Sto %g \n", op.Sfrom_pct_Fmax, op.Sto_pct_Fmax);
        }

        if (op.Sfrom_pct_Fmax == op.Sto_pct_Fmax) {
            return;
        }
    }
    else {
        if (verbose) {
            g_print(" Sfrom %g Sto %g \n", op.Sfrom, op.Sto);
        }

        if (op.Sfrom == op.Sto) {
            return;
        }
    }

    //find appropriate data, which should be fitted, use same regime as for the main fit
    if (opdata->FinputS)
        range_to_indices(op.Srange_pct_Fmax[0] * 0.01 * fdc.Fmax, op.Srange_pct_Fmax[1] * 0.01 * fdc.Fmax,
                         fdc.Funload, TRUE, &istart, &iend, &ndata);
    else {
        range_to_indices(op.Srange[0], op.Srange[1], fdc.hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, no fitting can be performed
    ndata = MIN(ndata, fdc.hunload->res - istart);

    if (verbose)
        g_print("interval has %d points, from %d to %d, hrange [%g, %g ] nm, Frange [%g, %g] mN \n",
                ndata, istart, iend,
                fdc.hunload->data[iend], fdc.hunload->data[istart], fdc.Funload->data[iend], fdc.Funload->data[istart]);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    if (FIT_S_LR) {
        //linear regression
        x = gwy_data_line_new(ndata, 1.0, TRUE);
        y = gwy_data_line_new(ndata, 1.0, TRUE);

        for (i = 0; i < ndata; i++) {
            x->data[i] = log(fdc.hunload->data[istart + i] - op.hp);
            y->data[i] = log(fdc.Funload->data[istart + i]);
        }
    }
    else {
        //nonlinear regression
        x = gwy_data_line_part_extract(fdc.hunload, istart, ndata);
        y = gwy_data_line_part_extract(fdc.Funload, istart, ndata);
    }

    // perform fit and get results
    if (FIT_S_LR) {
        //linear regression
        op_fit_S_lr(x, y, &op, &fdc, instdata, area);
    }
    else {
        //nonlinear regression
        op_fit_S_nl(x, y, &op, &fdc, instdata, area);
    }

    // clean up
    g_object_unref(x);
    g_object_unref(y);
    destroy_FDdata(&fdc);

    Sc[j + CONTACTMAX] = op.S;
    Ec[j + CONTACTMAX] = op.Er;
    Hc[j + CONTACTMAX] = op.Hit;

    if (verbose) {
        g_print("\n\n");
    }
}

void op_uncertainty_montecarlo_run_calc(const OPdata *opdata, const OPUncdata *opunc, OPMCdata *opmc, FDdata *fddata,
					const Instdata *instdata, const Area *area)
{
    gint istart, iend, ndata;
    gint i, j, Nmc;
    OPdata *op;
    OPMCdata *mc;
    GwyDataLine *h, *F;
    GwyDataLine *hnew, *Fnew;
    GwyDataLine *hnoise, *Fnoise;
    gdouble hmax_old, Fmax_old;
    GRand *rng;
    gint nn;
    gdouble **dd;
    gdouble *xhist, *yhist;
    gboolean hist;

    // initialize data
    Nmc = opmc->Nmc;
    op = opmc->opmc;

    //prepare data
    for (i = 0; i < Nmc ; i++) {
        init_OPdata(op);
        op++;
    }

    op = opmc->opmc;

    //copy input parameters
    for (i = 0; i < Nmc ; i++) {
        /* OPdata no longer contains these; use Instdata instead */
        /* op->nu = args->instdata.nu; */
        /* op->nui = args->instdata.nui; */
        /* op->Ei = args->instdata.Ei; */
        //		op->from = args->opdata.from;
        //		op->to = args->opdata.to;
        op++;
    }

    /*
      if (verbose) {
      for (i=0; i< Nmc ;i++)
      g_print("nu %g \n", args->opmc.opmc[i].nu);
      }
    */

    // initialize random generator and datafield
    rng = g_rand_new();


    /* hp fit */

    //not necessary to check range, couldn't get here in the first place

    //find appropriate data, which should be fitted, must use real fitted range, use same regime as for the main fit
    if (opdata->Finputhp)
        range_to_indices(opdata->hprange_pct_Fmax[0] * 0.01 * fddata->Fmax,
                         opdata->hprange_pct_Fmax[1] * 0.01 * fddata->Fmax,
                         fddata->Funload, TRUE, &istart, &iend, &ndata);
    else {
        range_to_indices(opdata->hprange[0], opdata->hprange[1], fddata->hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, not fitting can be performed
    ndata = MIN(ndata, fddata->hunload->res - istart);

    if (ndata < 2) {
        g_rand_free(rng);
        return;
    }

    //back up hmax, Fmax
    hmax_old = fddata->hmax;
    Fmax_old = fddata->Fmax;

    //get data to be fitted
    h = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
    F = gwy_data_line_part_extract(fddata->Funload, istart, ndata);

    //prepare datafield for noise
    hnoise = gwy_data_line_new_alike(h, FALSE);
    Fnoise = gwy_data_line_new_alike(h, FALSE);

    //prepare datafield for varied fields
    hnew = gwy_data_line_new_alike(h, FALSE);
    Fnew = gwy_data_line_new_alike(F, FALSE);

    op = opmc->opmc;

    for (i = 0; i < Nmc; i++) {

        //create noise for h-F data
        generate_uncorrelated_noise(hnoise, opunc->uh, rng);
        generate_uncorrelated_noise(Fnoise, opunc->uF, rng);

        //add noise to original h-F data
        if (!sum_data_lines(hnew, h, hnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        if (!sum_data_lines(Fnew, F, Fnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        //add random noise to hmax, Fmax
        fddata->hmax = hmax_old + gaussian_random_number(rng) * opunc->uh;
        fddata->Fmax = Fmax_old + gaussian_random_number(rng) * opunc->uF;

        //run fit
        // if (verbose) {
        //	g_print("i %d \n", i);
        //	g_print("-----\n");
        //	for (j =0;j < hnew->res; j++)
        //		g_print("%g %g \n", hnew->data[j], Fnew->data[j]);
        // }

        op_fit_hp(hnew, Fnew, op, fddata);
        // if (verbose)
        // g_print(" hp %g \n", op->hp);
        op++;
    }

    // clean up
    g_object_unref(h);
    g_object_unref(F);
    g_object_unref(hnew);
    g_object_unref(Fnew);
    g_object_unref(hnoise);
    g_object_unref(Fnoise);

    /*            S fit */

    //not necessary to check range, couldn't get here in the first place

    //find appropriate data, which should be fitted, must use real fitted range, use same regime as for the main fit
    if (opdata->FinputS)
        range_to_indices(opdata->Srange_pct_Fmax[0] * 0.01 * fddata->Fmax,
                         opdata->Srange_pct_Fmax[1] * 0.01 * fddata->Fmax,
                         fddata->Funload, TRUE, &istart, &iend, &ndata);
    else {
        range_to_indices(opdata->Srange[0], opdata->Srange[1], fddata->hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, not fitting can be performed
    ndata = MIN(ndata, fddata->hunload->res - istart);

    if (ndata < 2) {
        return;
    }

    //get data to be fitted
    h = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
    F = gwy_data_line_part_extract(fddata->Funload, istart, ndata);

    //prepare datafield for noise
    hnoise = gwy_data_line_new_alike(h, FALSE);
    Fnoise = gwy_data_line_new_alike(h, FALSE);

    //prepare datafield for varied fields
    hnew = gwy_data_line_new_alike(h, FALSE);
    Fnew = gwy_data_line_new_alike(F, FALSE);

    op = opmc->opmc;

    for (i = 0; i < Nmc; i++) {

        //create noise for h-F data
        generate_uncorrelated_noise(hnoise, opunc->uh, rng);
        generate_uncorrelated_noise(Fnoise, opunc->uF, rng);

        //add noise to original h-F data
        if (!sum_data_lines(hnew, h, hnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        if (!sum_data_lines(Fnew, F, Fnoise)) {
            g_printerr("Should not get here.\n");
            return;
        }

        //add random noise to hmax, Fmax
        fddata->hmax = hmax_old + gaussian_random_number(rng) * opunc->uh;
        fddata->Fmax = Fmax_old + gaussian_random_number(rng) * opunc->uF;

        //if necessary calculate logarithm
        if (FIT_S_LR) {
            for (j = 0 ; j < hnew->res; j++) {
                hnew->data[j] = log(hnew->data[j] - op->hp);
                Fnew->data[j] = log(Fnew->data[j]);
            }
        }

        //		printf("i %d \n", i);
        //		printf("-----\n");
        //		printf("hnew, Fnew \n");
        //		for (j =0;j < hnew->res; j++)
        //			printf("%g %g \n", hnew->data[j], Fnew->data[j]);


        //run fit
        op_fit_hp(hnew, Fnew, op, fddata);

        if (FIT_S_LR) {	//linear regression
            op_fit_S_lr(hnew, Fnew, op, fddata, instdata, area);
        }
        else {	//nonlinear regression
            op_fit_S_nl(hnew, Fnew, op, fddata, instdata, area);
        }

        //		printf("H %g E %g \n", op->Hit, op->Eit);
        op++;
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
        //		printf(" %d Er %g Eit %g Hit %g \n", i, args->opmc.opmc[i].Er, args->opmc.opmc[i].Eit, args->opmc.opmc[i].Hit);
    }

    //calculate averages and stdevs
    //hp, Aphc, Hit, Er, Eit

    mc = opmc;

    for (i = 0; i < NPARAM_OP; i++) {
        mc->mcavg[i] = 0;
        mc->mcstd[i] = 0;
    }

    op = opmc->opmc;
    nn = 0;

    for (i = 0; i < Nmc; i++) {
        //check sanity
        if (op->hc > 0) {
            mc->mcavg[0] += op->hc;
            mc->mcavg[1] += op->Aphc;
            mc->mcavg[2] += op->Hit;
            mc->mcavg[3] += op->Er;
            mc->mcavg[4] += op->Eit;
            mc->mcavg[5] += op->S;
            mc->mcstd[0] += op->hc * op->hc;
            mc->mcstd[1] += op->Aphc * op->Aphc;
            mc->mcstd[2] += op->Hit * op->Hit;
            mc->mcstd[3] += op->Er * op->Er;
            mc->mcstd[4] += op->Eit * op->Eit;
            mc->mcstd[5] += op->S * op->S;
            nn++;
        }

        op++;
    }

    //use only reasonable data for average and standard deviation
    for (i = 0; i < NPARAM_OP; i++) {
        mc->mcavg[i] /= nn;
        mc->mcstd[i] /= nn;
        mc->mcstd[i] -= mc->mcavg[i] * mc->mcavg[i];
        mc->mcstd[i] = (mc->mcstd[i] < 0) ? 0 : mc->mcstd[i];
        mc->mcstd[i] = sqrt(mc->mcstd[i]);
    }

    // how many datasets were skipped
    opmc->skipmc = Nmc - nn;

    /*
       printf("mean  stdev \n");
       for (i = 0; i< NPARAM_OP; i++){
       printf(" %g   %g \n", avg[i], std[i]);
       }
       */

    // Histograms
    //create data structure

    dd = opmc->mcdd;
    op = opmc->opmc;

    for (i = 0; i < Nmc; i++) {
        dd[0][i] = op->hc;
        dd[1][i] = op->Aphc;
        dd[2][i] = op->Hit;
        dd[3][i] = op->Er;
        dd[4][i] = op->Eit;
        dd[5][i] = op->S;
        op++;
    }

    //prepare histogram structure
    for (i = 0; i < NPARAM_OP; i++) {
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

gchar* op_mc_export_data(const OPdata *opdata, const OPUncdata *opunc, const OPMCdata *opmc, const FDdata *fddata, const Instdata *instdata,
			 enum ExportFormat exportformat)
{
    OPdata *op;
    gint i, j;
    gchar str[200];
    gboolean is;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Monte Carlo simulation of Oliver Pharr analysis\n");
    g_string_append_printf(buf, "#Input uncertainties \n");
    g_string_append_printf(buf, "#u(h) = %g nm \n", opunc->uh);
    g_string_append_printf(buf, "#u(F) = %g mN \n", opunc->uF);
    g_string_append_printf(buf, "#Number of iterations  = %d  \n", opmc->Nmc);
    g_string_append_printf(buf, "\n");

    g_string_append_printf(buf, "#Parameters of Oliver Pharr analysis  \n");
    g_string_append_printf(buf, "#beta: %g \n", opdata->beta);
    g_string_append_printf(buf, "#nu: %g                 \n", instdata->nu);
    g_string_append_printf(buf, "#nu_indenter: %g        \n", instdata->nui);
    g_string_append_printf(buf, "#E_indenter: %g      GPa\n", instdata->Ei * 1e-9); //TODO units
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

    g_string_append_printf(buf, "\n\n");

    g_string_append_printf(buf, "#Results \n");

    for (i = 0; i < NPARAM_OP; i++) {
        is = nparam_to_string(i, str, sizeof(str));
        g_string_append_printf(buf, "#%s = (%g ± %g ) \n", str, opmc->mcavg[i], opmc->mcstd[i]);
    }

    g_string_append_printf(buf, "\n\n\n");

    // Histograms

    if (opmc->mcnstat > 0) {
        g_string_append_printf(buf, "#Histogram\n");

        for (i = 0; i < NPARAM_OP; i++) {
            g_string_append_printf(buf, "#Histogram\n");
            is = nparam_to_string(i, str, sizeof(str));

            if (is) {
                g_string_append_printf(buf, "#%s    count/a.u.  \n", str);

                for (j = 0; j < opmc->mcnstat; j++) {
                    g_string_append_printf(buf, "%g %g \n", opmc->mcxhist[i][j], opmc->mcyhist[i][j]);
                }

                g_string_append_printf(buf, "\n \n");
            }
        }
    }

    g_string_append_printf(buf, "\n\n\n");

    //Full data
    g_string_append_printf(buf, "#Full data \n");
    g_string_append_printf(buf, "#");

    for (j = 0; j < NPARAM_OP; j++) {
        is = nparam_to_string(j, str, sizeof(str)); /* _is_ is not used anywhere */
        g_string_append_printf(buf, "%s ", str);
    }

    g_string_append_printf(buf, "\n");

    op = opmc->opmc;

    for (i = 0; i < opmc->Nmc; i++) {
        g_string_append_printf(buf, DEPTH, op->hc);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, AREA, op->Aphc);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, HARD, op->Hit);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, MODUL, op->Er);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, MODUL, op->Eit);
        g_string_append_printf(buf, "  ");
        g_string_append_printf(buf, SLOPE, op->S);
        g_string_append_printf(buf, "  \n");
        op++;
    }

    return g_string_free(buf, FALSE);
}

void init_OPMCdata(OPMCdata *opmc, gint Nmc)
{
    gint i;

    //number of iterations
    opmc->Nmc = Nmc;

    if (verbose) {
        g_print("initialize for %d iterations \n", opmc->Nmc);
    }

    //how much had to be skipped
    opmc->skipmc = 0;

    //dynamically allocate
    opmc->opmc = (OPdata *)g_malloc(opmc->Nmc * sizeof(OPdata));
    opmc->mcavg = (gdouble *)g_malloc(NPARAM_OP * sizeof(gdouble));
    opmc->mcstd = (gdouble *)g_malloc(NPARAM_OP * sizeof(gdouble));
    opmc->mcdd = (gdouble **) g_malloc(NPARAM_OP * sizeof(gdouble *));
    opmc->mcxhist = (gdouble **) g_malloc(NPARAM_OP * sizeof(gdouble *));
    opmc->mcyhist = (gdouble **) g_malloc(NPARAM_OP * sizeof(gdouble *));

    opmc->mcnstat = (gint)floor(3.49 * cbrt(opmc->Nmc) + 0.5);

    for (i = 0; i < NPARAM_OP; i++) {
        opmc->mcdd[i] = (gdouble *) g_malloc(opmc->Nmc * sizeof(gdouble));
        opmc->mcxhist[i] = (gdouble *) g_malloc(opmc->mcnstat * sizeof(gdouble));
        opmc->mcyhist[i] = (gdouble *) g_malloc(opmc->mcnstat * sizeof(gdouble));
    }
}

void destroy_OPMCdata(OPMCdata *opmc)
{
    gint i;

    //free everything
    g_free(opmc->mcavg);
    g_free(opmc->mcstd);

    for (i = 0; i < NPARAM_OP; i++) {
        g_free(opmc->mcdd[i]);
        g_free(opmc->mcxhist[i]);
        g_free(opmc->mcyhist[i]);
    }

    g_free(opmc->mcdd);
    g_free(opmc->mcxhist);
    g_free(opmc->mcyhist);

    opmc->mcnstat = 0;
    opmc->skipmc = 0;
    g_free(opmc->opmc);
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

gboolean op_unc_has_results(const OPUncdata *opuncdata)
{
    return (opuncdata->Ec != NULL && opuncdata->Hc != NULL);
}

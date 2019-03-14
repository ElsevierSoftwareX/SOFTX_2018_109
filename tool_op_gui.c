#include "tool_op_gui.h"
#include "tool_op_unc_gui.h"
#include "tool_op.h"

#include "controls.h"
#include "datatypes.h"
#include "fit-utils.h"
#include "gui-utils.h"
#include "niget-common.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gtk/gtk.h>
#include <libgwydgets/gwydgets.h>

#include <app/file.h>

static void op_remove_all_fit_labels_S(OPControls *opcontrols);
static void op_remove_all_fit_labels_hp(OPControls *opcontrols);
static void op_remove_all_fit_curves(OPControls *opcontrols);

static void op_save_data(Args *args);

static void update_models_op(OPControls *opcontrols);
static void op_zoom_in(Data *data);
static void op_zoom_out(Data *data);


void op_remove_fit_results_labels(Data *data)
{
    gwy_selection_clear(data->controls->opcontrols->selection);
    op_remove_all_fit_results_hp(&(data->args->opdata));
    op_remove_all_fit_labels_hp(data->controls->opcontrols);
    op_remove_all_fit_curves(data->controls->opcontrols);

    if (data->args->opunc.Ec != NULL || data->args->opunc.Hc != NULL) {
        op_unc_close_window(data);
    }
}

static void op_remove_all_fit_labels_S(OPControls *opcontrols)
{
    /* remove all labels that result from a powerlaw (S) fit */

    entry_clear(opcontrols->from);
    entry_clear(opcontrols->to);
    entry_clear(opcontrols->from_pct_Fmax);
    entry_clear(opcontrols->to_pct_Fmax);

    label_clear(opcontrols->SR);
    label_clear(opcontrols->Schi2);
    label_clear(opcontrols->SSSres);
    label_clear(opcontrols->SR2);
    label_clear(opcontrols->SR2adj);
    label_clear(opcontrols->S);
    label_clear(opcontrols->hc);

    label_clear(opcontrols->m);
    label_clear(opcontrols->eps);
    label_clear(opcontrols->Aphc);
    label_clear(opcontrols->extrapol);

    label_clear(opcontrols->Hit);
    label_clear(opcontrols->Er);
    label_clear(opcontrols->Eit);
    label_clear(opcontrols->Srange);
    label_clear(opcontrols->Srange_pct_Fmax);
}

static void op_remove_all_fit_labels_hp(OPControls *opcontrols)
{
    /* remove all labels that result from a linear (hp) fit */
    /* labels associated with the results of a powerlaw (S) fit must be removed as well */

    entry_clear(opcontrols->from);
    entry_clear(opcontrols->to);
    entry_clear(opcontrols->from_pct_Fmax);
    entry_clear(opcontrols->to_pct_Fmax);

    label_clear(opcontrols->hp);
    label_clear(opcontrols->hpR);
    label_clear(opcontrols->hprange);
    label_clear(opcontrols->hprange_pct_Fmax);
    label_clear(opcontrols->hpchi2);
    label_clear(opcontrols->hpSSres);
    label_clear(opcontrols->hpR2);
    label_clear(opcontrols->hpR2adj);

    op_remove_all_fit_labels_S(opcontrols);
}

static void op_remove_all_fit_curves(OPControls *opcontrols)
{
    gwy_graph_model_remove_curve_by_description(opcontrols->graph_model, "S fit");
    gwy_graph_model_remove_curve_by_description(opcontrols->graph_model, "hp fit");
}

static void recalculate_Eit(OPControls *opcontrols, OPdata *opdata, Instdata *instdata)
{
    gdouble Er;
    
    Er = opdata->Er * 1e9;

    // recalc for changed  nu, nui or Ei
    if (Er == 0) {
        return;
    }

    opdata->Eit = calc_Eit(Er, instdata->nu, instdata->Ei, instdata->nui) * 1e-9;
    label_set_gdouble_format(opcontrols->Eit, opdata->Eit, MODUL);
}

void op_recalculate(Data *data)
{
    recalculate_Eit(data->controls->opcontrols, &(data->args->opdata), &(data->args->instdata));
    op_propagate_uncertainties(data->controls->opunccontrols, &(data->args->opdata), &(data->args->opunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

void op_redraw(Data *data)
{
    // update loading curve and Fmax and hmax labels
    if (data->args->fddata.hunload) {
        gwy_graph_curve_model_set_data(data->controls->opcontrols->cmodelunload,
                                       data->args->fddata.hunload->data, data->args->fddata.Funload->data,
                                       data->args->fddata.hunload->res);

        label_set_gdouble_format(data->controls->opcontrols->hmax, data->args->fddata.hmax, DEPTH);
        label_set_gdouble_format(data->controls->opcontrols->Fmax, data->args->fddata.Fmax, FORCE);

        if (gwy_graph_model_get_n_curves(data->controls->opcontrols->graph_model) >= 2) {
            gwy_graph_curve_model_set_data(data->controls->opcontrols->cmodelhpfit,
                                           data->args->opdata.xhpfit->data, data->args->opdata.yhpfit->data,  data->args->opdata.xhpfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(data->controls->opcontrols->cmodelhpfit);
#endif
	    
            if (gwy_graph_model_get_n_curves(data->controls->opcontrols->graph_model) == 3) {
                gwy_graph_curve_model_set_data(data->controls->opcontrols->cmodelSfit,
                                               data->args->opdata.xSfit->data, data->args->opdata.ySfit->data,  data->args->opdata.xSfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
                gwy_graph_curve_model_enforce_order(data->controls->opcontrols->cmodelSfit);
#endif
            }
        }
    }
}

static void op_recalculate_beta(OPControls *opcontrols, OPdata *opdata, const Instdata *instdata)
{
    if (!opdata->has_fit_S) {
        return;
    }

    /* duplicate code with tool_op.c */
    opdata->Er = calc_Er_OP(opdata->Aphc, opdata->S, opdata->beta) * 1e15;
    opdata->Eit = calc_Eit(opdata->Er, instdata->nu, instdata->Ei, instdata->nui);
    opdata->Eit *= 1e-9;
    opdata->Er *= 1e-9;

    label_set_gdouble_format(opcontrols->Eit, opdata->Eit, MODUL);
    label_set_gdouble_format(opcontrols->Er, opdata->Er, MODUL);
}

static void op_beta_changed(Data *data)
{
    OPControls *fcontrols;

    fcontrols = data->controls->opcontrols;

    if (beta_changed(fcontrols->beta, &(data->args->opdata.beta))) {
        op_recalculate_beta(data->controls->opcontrols, &(data->args->opdata), &(data->args->instdata));
    }
}

static void op_reset_beta(Data *data)
{
    data->args->opdata.beta =  BETA_DEFAULT_OP;
    entry_set_gdouble_format(data->controls->opcontrols->beta, data->args->opdata.beta,  "%.03f");
    op_recalculate_beta(data->controls->opcontrols, &(data->args->opdata), &(data->args->instdata));
}

static void op_range_from_changed(Data *data)
{
    OPControls *fcontrols;
    OPdata *fdata;

    fcontrols = data->controls->opcontrols;
    fdata = &(data->args->opdata);

    range_from_changed(fcontrols->from, fcontrols->to, fcontrols->selection, &(fdata->hpfrom), &(fdata->Sfrom));
}

static void op_range_to_changed(Data *data)
{
    OPControls *tcontrols;
    OPdata *tdata;

    tcontrols = data->controls->opcontrols;
    tdata = &(data->args->opdata);

    range_to_changed(tcontrols->from, tcontrols->to, tcontrols->selection, &(tdata->hpto), &(tdata->Sto));
}

static gboolean op_range_from_pct_Fmax_changed(Data *data)
{
    OPControls *rcontrols;
    OPdata *rdata;

    rcontrols = data->controls->opcontrols;
    rdata = &(data->args->opdata);

    return range_from_pct_Fmax_changed(rcontrols->from_pct_Fmax, rcontrols->to_pct_Fmax, rcontrols->selection, data->args->fddata.Fmax,
				       data->args->fddata.hunload, data->args->fddata.Funload,
				       &(rdata->hpfrom_pct_Fmax), &(rdata->Sfrom_pct_Fmax), &(rdata->hpto_pct_Fmax), &(rdata->Sto_pct_Fmax),
				       &(rdata->Finputhp), &(rdata->FinputS));
}

static gboolean op_range_to_pct_Fmax_changed(Data *data)
{
    OPControls *rcontrols;
    OPdata *rdata;

    rcontrols = data->controls->opcontrols;
    rdata = &(data->args->opdata);

    return range_to_pct_Fmax_changed(rcontrols->from_pct_Fmax, rcontrols->to_pct_Fmax, rcontrols->selection, data->args->fddata.Fmax,
				     data->args->fddata.hunload, data->args->fddata.Funload,
                                     &(rdata->hpfrom_pct_Fmax), &(rdata->Sfrom_pct_Fmax), &(rdata->hpto_pct_Fmax), &(rdata->Sto_pct_Fmax),
				     &(rdata->Finputhp), &(rdata->FinputS));
}

static void op_graph_selected(GwySelection *selection, gint i, Data *data)
{
    gdouble range[2];
    gint istart, iend, ndata;
    gdouble Fmax;

    g_return_if_fail(i <= 0);

    //take data from selection and put them in the "from" and "to" values for both hp and S (at this moment, we don't know yet which fit will be chosen), create labels only once
    //OPControls has entries from and to, which are filled here
    // and labels hp,Srange, where the values of the ranges that were actually used in the corresponding fits  are shown
    //OPdata has hpfrom,hpto, Sfrom,Sto which are the data from the entries, i.e. from the selection
    //OPdata; also has hprange, Srange where the ranges that were actually used in the fit

    gwy_selection_get_object(selection, 0, range);

    data->args->opdata.Sfrom = MIN(range[0], range[1]);
    data->args->opdata.Sto = MAX(range[0], range[1]);

    data->args->opdata.hpfrom = MIN(range[0], range[1]);
    data->args->opdata.hpto = MAX(range[0], range[1]);

    entry_set_gdouble_format(data->controls->opcontrols->from, data->args->opdata.Sfrom, DEPTH);
    entry_set_gdouble_format(data->controls->opcontrols->to, data->args->opdata.Sto, DEPTH);

    Fmax = data->args->fddata.Fmax;
    range_to_indices(data->args->opdata.Sfrom, data->args->opdata.Sto, data->args->fddata.hunload, TRUE, &istart, &iend, &ndata);
    data->args->opdata.hpfrom_pct_Fmax = data->args->fddata.Funload->data[iend] / Fmax * 100;
    data->args->opdata.hpto_pct_Fmax = data->args->fddata.Funload->data[istart] / Fmax * 100;
    data->args->opdata.Sfrom_pct_Fmax = data->args->fddata.Funload->data[iend] / Fmax * 100;
    data->args->opdata.Sto_pct_Fmax = data->args->fddata.Funload->data[istart] / Fmax * 100;

    entry_set_gdouble_format(data->controls->opcontrols->from_pct_Fmax, data->args->opdata.Sfrom_pct_Fmax, DEPTH);
    entry_set_gdouble_format(data->controls->opcontrols->to_pct_Fmax, data->args->opdata.Sto_pct_Fmax, DEPTH);

    data->args->opdata.Finputhp = FALSE;
    data->args->opdata.FinputS = FALSE;
}

static void op_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (op_has_results(&(args->opdata))) {
	buffer = op_export_data(&(args->opdata), &(args->fddata), &(args->instdata), &(args->area), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

void op_fit_and_update_hp(Data *data)
{
    LinReg reg;
    gchar str[300];
    gint nc;
    gint istart, iend, ndata;
    gboolean fit = TRUE;
    GwyDataLine *x, *y;

    OPdata *opdata;
    OPControls *opcontrols;
    FDdata *fddata;

    //TEMPORARY update GtkEntries for ranges
    /*
    if (data->args->opdata.Finputhp){
    	range_from_pct_Fmax_changed(data);
    	range_to_pct_Fmax_changed(data);
    }
    else{
    	range_from_changed(data);
    	range_to_changed(data);
    }
    */

    //create pointers for easier reading
    opdata = &(data->args->opdata);
    opcontrols = data->controls->opcontrols;
    fddata = &(data->args->fddata);

    //if there was a previous fit, invalidate previous fitted curve
    if (opdata->has_fit_hp) {
        g_object_unref(opdata->xhpfit);
        g_object_unref(opdata->yhpfit);
    }

    //check range
    if (opdata->hpfrom == opdata->hpto) {
        fit = FALSE;
    }

    //find appropriate data, which should be fitted
    //default regime is governed by depth, can be switched to load regime by Finput (if the pct_Fmax ranges were used)
    if (opdata->Finputhp) {
        range_to_indices(opdata->hpfrom_pct_Fmax * 0.01 * fddata->Fmax, opdata->hpto_pct_Fmax * 0.01 * fddata->Fmax,
                         fddata->Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(opdata->hpfrom, opdata->hpto, fddata->hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, no fitting can be performed
    ndata = MIN(ndata, fddata->hunload->res - istart);

    if (ndata < 2) {
        fit = FALSE;
    }

    opdata->nfithpdata = ndata;

    if (fit) {
        //get data to be fitted
        x = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
        y = gwy_data_line_part_extract(fddata->Funload, istart, ndata);

        //fit the data and evaluate
        op_fit_hp(x, y, opdata, &(data->args->fddata));
        reg = opdata->reghp;


        // remove all previous S fit results, with the exception of the range which is common for both hp fit and S fit
        if (opdata->has_fit_S) {
            g_object_unref(opdata->xSfit);
            g_object_unref(opdata->ySfit);
        }

        op_remove_all_fit_labels_S(opcontrols);
        op_remove_all_fit_results_S(opdata);
        //update range labels
        entry_set_gdouble_format(opcontrols->from, data->args->opdata.hpfrom, DEPTH);
        entry_set_gdouble_format(opcontrols->to, data->args->opdata.hpto, DEPTH);
        entry_set_gdouble_format(opcontrols->from_pct_Fmax, data->args->opdata.hpfrom_pct_Fmax, DEPTH);
        entry_set_gdouble_format(opcontrols->to_pct_Fmax, data->args->opdata.hpto_pct_Fmax, DEPTH);

        //update labels
        label_set_gdouble_format(opcontrols->hp, opdata->hp, DEPTH);
        label_set_gdouble_format(opcontrols->hpR, opdata->reghp.corr, NUMBER);
        label_set_gdouble_format(opcontrols->hpSSres, opdata->reghp.SSres, "%.4g");
        label_set_gdouble_format(opcontrols->hpR2, opdata->reghp.R2, "%.4g");
        label_set_gdouble_format(opcontrols->hpR2adj, opdata->reghp.R2adj, "%.4g");
        label_set_gdouble_format(opcontrols->hpchi2, opdata->reghp.chi2, "%.4g");

        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", opdata->hpfrom, opdata->hpto);
        gtk_label_set_text(GTK_LABEL(opcontrols->hprange), str);
        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", opdata->hpfrom_pct_Fmax, opdata->hpto_pct_Fmax);
        gtk_label_set_text(GTK_LABEL(opcontrols->hprange_pct_Fmax), str);

        //calculate curve of fit
        opdata->xhpfit = gwy_data_line_duplicate(x);
        opdata->yhpfit = gwy_data_line_new_alike(x, FALSE);
        gwy_data_line_copy(opdata->xhpfit, opdata->yhpfit);
        gwy_data_line_multiply(opdata->yhpfit, reg.slope);
        gwy_data_line_add(opdata->yhpfit, reg.intercept);

        // clean up
        g_object_unref(x);
        g_object_unref(y);

        //update fitted curve, if necessary create it
        nc = gwy_graph_model_get_n_curves(opcontrols->graph_model);

        if (nc >= 2) {
            gwy_graph_curve_model_set_data(opcontrols->cmodelhpfit, opdata->xhpfit->data, opdata->yhpfit->data, opdata->xhpfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(opcontrols->cmodelhpfit);
#endif
        }
        else {
            gwy_graph_curve_model_set_data(opcontrols->cmodelhpfit, opdata->xhpfit->data, opdata->yhpfit->data, opdata->xhpfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(opcontrols->cmodelhpfit);
#endif
            gwy_graph_model_add_curve(opcontrols->graph_model, opcontrols->cmodelhpfit);
            g_object_set(opcontrols->cmodelhpfit,
                         "mode", GWY_GRAPH_CURVE_LINE,
                         "description", "hp fit",
                         "color", CURVE_COLOR_FIT,
                         "line-style", GDK_LINE_SOLID,
                         "line-width", 3,
                         NULL);
        }

        opdata->has_fit_S = FALSE;
        gtk_widget_set_sensitive(opcontrols->button_fit_S, TRUE);
    }
    else {
        nc = gwy_graph_model_get_n_curves(opcontrols->graph_model);
        op_remove_all_fit_labels_hp(opcontrols);
        op_remove_all_fit_results_hp(opdata);
        opdata->has_fit_hp = FALSE;
        gtk_widget_set_sensitive(opcontrols->button_fit_S, FALSE);
        opdata->has_fit_S = FALSE;

        if (nc >= 2) {
            if (nc == 3) {
                gwy_graph_model_remove_curve_by_description(opcontrols->graph_model, "S fit");
            }

            gwy_graph_model_remove_curve_by_description(opcontrols->graph_model, "hp fit");
        }
    }

    if (data->args->opunc.Ec != NULL || data->args->opunc.Hc != NULL) {
        op_unc_close_window(data);
    }

    gtk_widget_set_sensitive(opcontrols->button_unc, FALSE);
}

void op_fit_and_update_S(Data *data)
{
    LinReg reg;
    gchar str[300];
    gint nc, i;
    gint istart, iend, ndata;
    gboolean fit = TRUE;
    GwyDataLine *x, *y;

    OPdata *opdata;
    OPControls *opcontrols;
    FDdata *fddata;
    GdkColor color;

    gdk_color_parse("red", &color);

    // TEMPORARY update GtkEntries for ranges
    /*
    if (data->args->opdata.FinputS){
    	range_from_pct_Fmax_changed(data);
    	range_to_pct_Fmax_changed(data);
    }
    else{
    	range_from_changed(data);
    	range_to_changed(data);
    }
    */

    //create pointers for easier reading
    opdata = &(data->args->opdata);
    opcontrols = data->controls->opcontrols;
    fddata = &(data->args->fddata);

    //if there was a previous fit, free the old data and create new data
    if (opdata->has_fit_S) {
        g_object_unref(opdata->xSfit);
        g_object_unref(opdata->ySfit);
    }

    // check range
    if (opdata->Sfrom == opdata->Sto) {
        fit = FALSE;
    }

    //find appropriate data, which should be fitted
    //default regime is governed by depth, can be switched to load regime by Finput (if the pct_Fmax ranges were used)
    if (opdata->FinputS) {
        range_to_indices(opdata->Sfrom_pct_Fmax * 0.01 * fddata->Fmax, opdata->Sto_pct_Fmax * 0.01 * fddata->Fmax,
                         fddata->Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(opdata->Sfrom, opdata->Sto, fddata->hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range,no fitting can be performed
    ndata = MIN(ndata, fddata->hunload->res - istart);

    if (ndata < 2) {
        fit = FALSE;
    }

    opdata->nfitSdata = ndata;

    if (fit) {
        //get data to be fitted
        if (FIT_S_LR) {
            //linear regression
            x = gwy_data_line_new(ndata, 1.0, TRUE);
            y = gwy_data_line_new(ndata, 1.0, TRUE);

            for (i = 0; i < ndata; i++) {
                x->data[i] = log(fddata->hunload->data[istart + i] - opdata->hp);
                y->data[i] = log(fddata->Funload->data[istart + i]);
            }

            op_fit_S_lr(x, y, opdata, &(data->args->fddata), &(data->args->instdata), &(data->args->area));
            reg = opdata->regS;
        }
        else {
            //nonlinear regression
            x = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
            y = gwy_data_line_part_extract(fddata->Funload, istart, ndata);
            op_fit_S_nl(x, y, opdata, &(data->args->fddata), &(data->args->instdata), &(data->args->area));
        }

        //update labels
        label_set_gdouble_format(opcontrols->hp, opdata->hp, DEPTH);

        label_set_gdouble_format(opcontrols->SSSres, opdata->regS.SSres, "%.4g");
        label_set_gdouble_format(opcontrols->SR2, opdata->regS.R2, "%.4g");
        label_set_gdouble_format(opcontrols->SR2adj, opdata->regS.R2adj, "%.4g");
        label_set_gdouble_format(opcontrols->Schi2, opdata->regS.chi2, "%.4g");

        label_set_gdouble_format(opcontrols->m, opdata->m, NUMBER);
        label_set_gdouble_format(opcontrols->eps, opdata->eps, NUMBER);
        label_set_gdouble_format(opcontrols->S, opdata->S, SLOPE);
        label_set_gdouble_format(opcontrols->hc, opdata->hc, DEPTH);
        label_set_gdouble_format(opcontrols->Aphc, opdata->Aphc, AREA);

        if (data->args->area.mode == AREA_DATA && opdata->hc > data->args->area.xmax) {
            gtk_label_set_text(GTK_LABEL(opcontrols->extrapol), "Warning: Extrapolation in area calibration!");
            gtk_widget_modify_fg(opcontrols->extrapol, GTK_STATE_NORMAL, &color);
        }
        else {
            label_clear(opcontrols->extrapol);
        }

        /*
        label_set_gdouble_format(opcontrols->Schi2, opdata->regS.chi, NUMBER);
        label_set_gdouble_format(opcontrols->SR, opdata->regS.corr, NUMBER);
        */

        label_set_gdouble_format(opcontrols->Hit, opdata->Hit, HARD);
        label_set_gdouble_format(opcontrols->Eit, opdata->Eit, MODUL);
        label_set_gdouble_format(opcontrols->Er, opdata->Er, MODUL);

        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", opdata->Sfrom, opdata->Sto);
        gtk_label_set_text(GTK_LABEL(opcontrols->Srange), str);
        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", opdata->Sfrom_pct_Fmax, opdata->Sto_pct_Fmax);
        gtk_label_set_text(GTK_LABEL(opcontrols->Srange_pct_Fmax), str);

        //calculate curve of fit

        if (FIT_S_LR) {
            //linear regression
            opdata->xSfit = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
            opdata->ySfit = gwy_data_line_new_alike(y, FALSE);

            for (i = 0; i < opdata->xSfit->res; i++) {
                opdata->ySfit->data[i] = exp(reg.intercept) * pow(opdata->xSfit->data[i] - opdata->hp, opdata->m);
            }
        }
        else {
            //, nonlinear regression
            opdata->xSfit = gwy_data_line_duplicate(x);
            opdata->ySfit = gwy_data_line_new_alike(y, FALSE);

            for (i = 0; i < opdata->xSfit->res; i++) {
                opdata->ySfit->data[i] =  opdata->alpha * pow(opdata->xSfit->data[i] - opdata->hp, opdata->m);
            }
        }

        // clean up
        g_object_unref(x);
        g_object_unref(y);

        //update fitted curve, if necessary create it
        nc = gwy_graph_model_get_n_curves(opcontrols->graph_model);

        if (nc == 3) {
            gwy_graph_curve_model_set_data(opcontrols->cmodelSfit, opdata->xSfit->data, opdata->ySfit->data, opdata->xSfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(opcontrols->cmodelSfit);
#endif
        }
        else {
            gwy_graph_curve_model_set_data(opcontrols->cmodelSfit, opdata->xSfit->data, opdata->ySfit->data, opdata->xSfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(opcontrols->cmodelSfit);
#endif
            gwy_graph_model_add_curve(opcontrols->graph_model, opcontrols->cmodelSfit);
            g_object_set(opcontrols->cmodelSfit,
                         "mode", GWY_GRAPH_CURVE_LINE,
                         "description", "S fit",
                         "color", CURVE_COLOR_FIT2,
                         "line-style", GDK_LINE_SOLID,
                         "line-width", 3,
                         NULL);
        }
    }
    else {
        nc = gwy_graph_model_get_n_curves(opcontrols->graph_model);

        op_remove_all_fit_labels_S(opcontrols);
        op_remove_all_fit_results_S(opdata);
        opdata->has_fit_S = FALSE;

        if (nc == 3) {
            gwy_graph_model_remove_curve_by_description(opcontrols->graph_model, "S fit");
        }
    }

    // update previous uncertainties if existing
    if (opdata->has_fit_S) {
        if (data->args->opunc.Ec != NULL || data->args->opunc.Hc != NULL) {
            op_unc_close_window(data);
        }

        gtk_widget_set_sensitive(opcontrols->button_unc, TRUE);
    }
    else {
        gtk_widget_set_sensitive(opcontrols->button_unc, FALSE);
    }
}

GtkWidget *tool_op_create(Data *data)
{
    GwyGraphArea *area;
    Args *args;
    OPControls *opcontrols;
    gint i;

    args = data->args;
    opcontrols = data->controls->opcontrols;

    // GraphArea, three curves can be displayed: loading F-d curve, the linear (hp) fit and the powerlaw (S) fit

    opcontrols->graph_model = gwy_graph_model_new();
    g_object_set(opcontrols->graph_model, "title", "Unload", NULL);
    gwy_graph_model_set_axis_label(opcontrols->graph_model, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(opcontrols->graph_model, GTK_POS_LEFT, "F / mN");

    opcontrols->zoomstack = NULL;
    opcontrols->zoomstack = g_slist_prepend(opcontrols->zoomstack, opcontrols->graph_model);

    opcontrols->hbox_graph = gtk_hbox_new(FALSE, CTRLS_SPACING);

    opcontrols->graph = gwy_graph_new(opcontrols->graph_model);
    /* gtk_widget_set_size_request(GTK_WIDGET(opcontrols->graph), 600, 400); */

    gwy_graph_enable_user_input(GWY_GRAPH(opcontrols->graph), FALSE);
    gwy_graph_set_status(GWY_GRAPH(opcontrols->graph), GWY_GRAPH_STATUS_XSEL);


    area = GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(opcontrols->graph)));
    opcontrols->selection = gwy_graph_area_get_selection(area, GWY_GRAPH_STATUS_XSEL);
    gwy_selection_set_max_objects(opcontrols->selection, 1);
    g_signal_connect(opcontrols->selection, "changed", G_CALLBACK(op_graph_selected), data);

    opcontrols->cmodelunload = gwy_graph_curve_model_new();
    opcontrols->cmodelhpfit = gwy_graph_curve_model_new();
    opcontrols->cmodelSfit = gwy_graph_curve_model_new();

    g_object_set(opcontrols->cmodelunload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Unloading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_UNLOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);
    gwy_graph_model_add_curve(opcontrols->graph_model, opcontrols->cmodelunload);


    // Table of input data and results

    opcontrols->ctrls = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* Info section */

    opcontrols->frame_info = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(opcontrols->frame_info), GTK_WIDGET(gwy_label_new_header("Info")));
    opcontrols->table_info = gtk_table_new(2, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(opcontrols->table_info), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(opcontrols->table_info), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(opcontrols->table_info), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(opcontrols->frame_info), opcontrols->table_info);

    i = 0;

    /* h_max */

    gtk_table_attach(GTK_TABLE(opcontrols->table_info), label_new_with_markup_left("h<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    opcontrols->hmax = gtk_label_new(NULL);
    label_set_gdouble_format(opcontrols->hmax, args->fddata.hmax, DEPTH);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->hmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_info), opcontrols->hmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(opcontrols->table_info), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /* F_max */

    gtk_table_attach(GTK_TABLE(opcontrols->table_info), label_new_with_markup_left("F<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    opcontrols->Fmax = gtk_label_new(NULL);
    label_set_gdouble_format(opcontrols->Fmax, args->fddata.Fmax, FORCE);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_info), opcontrols->Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_info), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    /* Parameters section */

    opcontrols->frame_parameters = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(opcontrols->frame_parameters), GTK_WIDGET(gwy_label_new_header("Parameters")));
    opcontrols->table_parameters = gtk_table_new(5, 4, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(opcontrols->table_parameters), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(opcontrols->table_parameters), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(opcontrols->table_parameters), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(opcontrols->frame_parameters), opcontrols->table_parameters);
    i = 0;

    gtk_table_attach(GTK_TABLE(opcontrols->table_parameters), label_new_with_markup_left("Fit range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->from = gtk_entry_new();
    g_object_set_data(G_OBJECT(opcontrols->from), "id", (gpointer)"from");
    gtk_entry_set_width_chars(GTK_ENTRY(opcontrols->from), 8);
    gtk_table_attach(GTK_TABLE(opcontrols->table_parameters), opcontrols->from, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(opcontrols->from, TRUE);
    g_signal_connect_swapped(opcontrols->from, "activate", G_CALLBACK(op_range_from_changed), data);
    /*	g_signal_connect_swapped(opcontrols->from, "focus-out-event", G_CALLBACK(range_from_changed),data); */

    opcontrols->to = gtk_entry_new();
    g_object_set_data(G_OBJECT(opcontrols->to), "id", (gpointer)"to");
    gtk_entry_set_width_chars(GTK_ENTRY(opcontrols->to), 8);
    gtk_table_attach(GTK_TABLE(opcontrols->table_parameters), opcontrols->to, 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(opcontrols->to, TRUE);
    g_signal_connect_swapped(opcontrols->to, "activate", G_CALLBACK(op_range_to_changed), data);
    /*	g_signal_connect_swapped(opcontrols->to, "focus-out-event", G_CALLBACK(range_to_changed),data); */

    gtk_table_attach(GTK_TABLE(opcontrols->table_parameters), label_new_with_markup_left("nm"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    opcontrols->from_pct_Fmax = gtk_entry_new();
    g_object_set_data(G_OBJECT(opcontrols->from_pct_Fmax), "id", (gpointer)"from_pct_Fmax");
    gtk_entry_set_width_chars(GTK_ENTRY(opcontrols->from_pct_Fmax), 8);
    gtk_table_attach(GTK_TABLE(opcontrols->table_parameters), opcontrols->from_pct_Fmax, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(opcontrols->from_pct_Fmax, TRUE);
    g_signal_connect_swapped(opcontrols->from_pct_Fmax, "activate", G_CALLBACK(op_range_from_pct_Fmax_changed), data);
    /*	g_signal_connect_swapped(opcontrols->from_pct_Fmax, "focus-out-event", G_CALLBACK(range_from_pct_Fmax_changed),data);  */


    opcontrols->to_pct_Fmax = gtk_entry_new();
    g_object_set_data(G_OBJECT(opcontrols->to_pct_Fmax), "id", (gpointer)"to_pct_Fmax");
    gtk_entry_set_width_chars(GTK_ENTRY(opcontrols->to_pct_Fmax), 8);
    gtk_table_attach(GTK_TABLE(opcontrols->table_parameters), opcontrols->to_pct_Fmax, 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(opcontrols->to_pct_Fmax, TRUE);
    g_signal_connect_swapped(opcontrols->to_pct_Fmax, "activate", G_CALLBACK(op_range_to_pct_Fmax_changed), data);
    /*g_signal_connect_swapped(opcontrols->to_pct_Fmax, "focus-out-event", G_CALLBACK(range_to_pct_Fmax_changed),data);  */

    gtk_table_attach(GTK_TABLE(opcontrols->table_parameters), label_new_with_markup_left("% F<sub>max</sub>"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_parameters), label_new_with_markup_left("β"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->beta = gtk_entry_new();
    //	g_object_set_data(G_OBJECT(opcontrols->beta), "id", (gpointer)"beta");
    gtk_entry_set_width_chars(GTK_ENTRY(opcontrols->beta), 8);
    entry_set_gdouble_format(opcontrols->beta, data->args->opdata.beta, "%.03f");
    gwy_widget_set_activate_on_unfocus(opcontrols->beta, TRUE);
    gtk_table_attach(GTK_TABLE(opcontrols->table_parameters), opcontrols->beta, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(opcontrols->beta, "activate", G_CALLBACK(op_beta_changed), data);

    opcontrols->reset_beta = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(opcontrols->reset_beta), gtk_image_new_from_stock(GTK_STOCK_REVERT_TO_SAVED, GTK_ICON_SIZE_BUTTON));
    gtk_table_attach(GTK_TABLE(opcontrols->table_parameters), opcontrols->reset_beta, 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(opcontrols->reset_beta, "clicked", G_CALLBACK(op_reset_beta), data);


    /* Fit buttons */

    opcontrols->hbox_fitbuttons = gtk_hbox_new(TRUE, 0);

    opcontrols->button_fit_hp = gtk_button_new_with_mnemonic(gwy_sgettext("_Fit hp"));
    //opcontrols->button_fit_hp = gtk_button_new_with_mnemonic (NULL);
    //gtk_container_add (GTK_CONTAINER (opcontrols->button_fit_hp), label_new_with_markup_with_mnemonic ("_Fit h<sub>p</sub>"));
    gtk_button_set_image(GTK_BUTTON(opcontrols->button_fit_hp), gtk_image_new_from_stock(GTK_STOCK_EXECUTE, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(opcontrols->hbox_fitbuttons), opcontrols->button_fit_hp, TRUE, TRUE, 0);
    g_signal_connect_swapped(opcontrols->button_fit_hp, "clicked", G_CALLBACK(op_fit_and_update_hp), data);

    opcontrols->button_fit_S = gtk_button_new_with_mnemonic(gwy_sgettext("verb|F_it S"));
    gtk_button_set_image(GTK_BUTTON(opcontrols->button_fit_S), gtk_image_new_from_stock(GTK_STOCK_EXECUTE, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(opcontrols->hbox_fitbuttons), opcontrols->button_fit_S, TRUE, TRUE, 0);
    gtk_widget_set_sensitive(opcontrols->button_fit_S, FALSE);
    g_signal_connect_swapped(opcontrols->button_fit_S, "clicked", G_CALLBACK(op_fit_and_update_S), data);


    /* Results section */

    opcontrols->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(opcontrols->frame_results), GTK_WIDGET(gwy_label_new_header("Results")));

    opcontrols->scrolledwindow_results = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(opcontrols->scrolledwindow_results), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(opcontrols->frame_results), opcontrols->scrolledwindow_results);

    opcontrols->table_results = gtk_table_new(15, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(opcontrols->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(opcontrols->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(opcontrols->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(opcontrols->scrolledwindow_results), opcontrols->table_results);
    i = 0;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("h<sub>p</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->hp = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->hp), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->hp, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;


    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("SS<sub>res</sub> hp fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->hpSSres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->hpSSres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->hpSSres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("R<sup>2</sup> hp fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->hpR2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->hpR2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->hpR2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("R<sup>2</sup><sub>adj</sub> hp fit "), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->hpR2adj = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->hpR2adj), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->hpR2adj, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("χ<sup>2</sup> hp fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->hpchi2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->hpchi2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->hpchi2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    //	gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("R h<sub>p</sub> fit"), 0, 1, i, i+1, GTK_FILL, 0, 0, 0);
    opcontrols->hpR = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->hpR), 0.0, 0.5);
    //	gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->hpR, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("m"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->m = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->m), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->m, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("SS<sub>res</sub> S fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->SSSres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->SSSres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->SSSres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("R<sup>2</sup> S fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->SR2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->SR2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->SR2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("R<sup>2</sup><sub>adj</sub> S fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->SR2adj = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->SR2adj), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->SR2adj, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("χ<sup>2</sup> S fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->Schi2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->Schi2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->Schi2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("ε"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->eps = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->eps), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->eps, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("h<sub>c</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->hc = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->hc), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->hc, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("S"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->S = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->S), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->S, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("mN/nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("A<sub>p</sub>(h<sub>c</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->Aphc = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->Aphc), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->Aphc, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("nm<sup>2</sup>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    opcontrols->extrapol = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->extrapol, 0, 3, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    //	gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("R S fit"), 0, 1, i, i+1, GTK_FILL, 0, 0, 0);
    opcontrols->SR = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->SR), 0.0, 0.5);
    //	gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->SR, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("H<sub>IT</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->Hit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->Hit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->Hit, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("MPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("E<sub>r</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->Er = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->Er), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->Er, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("E<sub>IT</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->Eit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->Eit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->Eit, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("Fitted h<sub>p</sub> range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->hprange = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->hprange), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->hprange, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    opcontrols->hprange_pct_Fmax = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->hprange_pct_Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->hprange_pct_Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("% F<sub>max</sub>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("Fitted S range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opcontrols->Srange = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->Srange), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->Srange, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    opcontrols->Srange_pct_Fmax = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opcontrols->Srange_pct_Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), opcontrols->Srange_pct_Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opcontrols->table_results), label_new_with_markup_left("% F<sub>max</sub>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    // Unc, Save Buttons
    opcontrols->button_unc = gtk_button_new_with_mnemonic("_Uncertainties");
    gtk_button_set_image(GTK_BUTTON(opcontrols->button_unc), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(opcontrols->button_unc, "clicked", G_CALLBACK(op_uncertainty), data);
    gtk_widget_set_sensitive(opcontrols->button_unc, FALSE);

    opcontrols->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    //gtk_table_attach(GTK_TABLE(opcontrols->ctrls), opcontrols->button_save, 0, 4, i, i+1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(opcontrols->button_save, "clicked", G_CALLBACK(op_save_data), args);

    gtk_box_pack_start(GTK_BOX(opcontrols->ctrls), opcontrols->frame_info, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(opcontrols->ctrls), opcontrols->frame_parameters, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(opcontrols->ctrls), opcontrols->hbox_fitbuttons, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(opcontrols->ctrls), opcontrols->frame_results, TRUE, TRUE, 0);
    // gtk_box_pack_start (GTK_BOX (opcontrols->ctrls), opcontrols->button_unc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(opcontrols->ctrls), opcontrols->button_save, FALSE, FALSE, 0);

    /* Zooming */

    opcontrols->vbox_zoom_buttons = gtk_vbox_new(TRUE, 0);
    opcontrols->vbox_zoom_outer = gtk_vbox_new(FALSE, 0);

    //opcontrols->button_zoom_in = gtk_button_new_with_label ("Zoom in");
    opcontrols->button_zoom_in = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(opcontrols->button_zoom_in), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(opcontrols->vbox_zoom_buttons), opcontrols->button_zoom_in, FALSE, FALSE, 0);
    g_signal_connect_swapped(opcontrols->button_zoom_in, "clicked", G_CALLBACK(op_zoom_in), data);

    //opcontrols->button_zoom_out = gtk_button_new_with_label ("Zoom out");
    opcontrols->button_zoom_out = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(opcontrols->button_zoom_out), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(opcontrols->vbox_zoom_buttons), opcontrols->button_zoom_out, FALSE, FALSE, 0);
    g_signal_connect_swapped(opcontrols->button_zoom_out, "clicked", G_CALLBACK(op_zoom_out), data);

    //opcontrols->button_zoom_restore = gtk_button_new_with_label ("Restore");
    opcontrols->button_zoom_restore = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(opcontrols->button_zoom_restore), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(opcontrols->vbox_zoom_buttons), opcontrols->button_zoom_restore, FALSE, FALSE, 0);
    g_signal_connect_swapped(opcontrols->button_zoom_restore, "clicked", G_CALLBACK(op_zoom_restore), data);

    gtk_box_pack_start(GTK_BOX(opcontrols->vbox_zoom_outer), opcontrols->vbox_zoom_buttons, TRUE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(opcontrols->hbox_graph), opcontrols->graph, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(opcontrols->hbox_graph), opcontrols->vbox_zoom_outer, FALSE, FALSE, 0);


    //zoom and unzoom buttons
    /* opcontrols->button_zoom_in = gtk_button_new(); */
    /* gtk_button_set_image (GTK_BUTTON (opcontrols->button_zoom_in), gtk_image_new_from_stock (GWY_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON)); */
    /* g_signal_connect_swapped(opcontrols->button_zoom_in, "clicked", G_CALLBACK (zoom_in), data); */

    /* opcontrols->button_zoom_out = gtk_button_new(); */
    /* gtk_button_set_image (GTK_BUTTON (opcontrols->button_zoom_out), gtk_image_new_from_stock (GWY_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON)); */
    /* g_signal_connect_swapped(opcontrols->button_zoom_out, "clicked", G_CALLBACK (zoom_out), data); */


    /* vbox = gtk_vbox_new(FALSE,FALSE); */
    /* gtk_box_pack_start (GTK_BOX(vbox), opcontrols->graph, TRUE, TRUE, 0); */
    /* gtk_box_pack_start (GTK_BOX(vbox), opcontrols->button_zoom_in, FALSE, FALSE, 0); */
    /* gtk_box_pack_start (GTK_BOX(vbox), opcontrols->button_zoom_out, FALSE, FALSE, 0); */


    //final

    opcontrols->op_gui = gtk_hbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(opcontrols->op_gui), opcontrols->ctrls, FALSE, FALSE, 0);   /* expand, fill, padding */
    gtk_box_pack_start(GTK_BOX(opcontrols->op_gui), opcontrols->hbox_graph, TRUE, TRUE, 0);

    gtk_widget_show(opcontrols->op_gui);

    return opcontrols->op_gui;
}

static void update_models_op(OPControls *opcontrols)
{
    GwyGraphCurveModel *cmodel;

    opcontrols->graph_model = gwy_graph_get_model(GWY_GRAPH(opcontrols->graph));
    cmodel = gwy_graph_model_get_curve_by_description(opcontrols->graph_model, "Unloading");

    if (cmodel) {
        opcontrols->cmodelunload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(opcontrols->graph_model, "hp fit");

    if (cmodel) {
        opcontrols->cmodelhpfit = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(opcontrols->graph_model, "S fit");

    if (cmodel) {
        opcontrols->cmodelSfit = cmodel;
    }
}

static void op_zoom_in(Data *data)
{
    OPControls *opcontrols;

    opcontrols = data->controls->opcontrols;

    rs_zoom_in(GWY_GRAPH(opcontrols->graph), &(opcontrols->zoomstack), FALSE);   /* do not rescale y */
    update_models_op(opcontrols);
}

static void op_zoom_out(Data *data)
{
    OPControls *opcontrols;

    opcontrols = data->controls->opcontrols;

    rs_zoom_out(GWY_GRAPH(opcontrols->graph), &(opcontrols->zoomstack));
    update_models_op(opcontrols);
}

void op_zoom_restore(Data *data)
{
    OPControls *opcontrols;

    opcontrols = data->controls->opcontrols;

    rs_zoom_restore(GWY_GRAPH(opcontrols->graph), &(opcontrols->zoomstack));
    update_models_op(opcontrols);
}

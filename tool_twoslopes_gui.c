#ifndef NOFORTRAN

#include "tool_twoslopes_gui.h"
#include "tool_twoslopes.h"
#include "tool_twoslopes_unc_gui.h"

#include "controls.h"
#include "datatypes.h"
#include "file-utils.h"
#include "fit-utils.h"
#include "gui-utils.h"
#include "niget-common.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>
#include <app/file.h>


static void slopes_remove_all_fit_labels_load(SlopesControls *slcontrols);
static void slopes_remove_all_fit_labels_unload(SlopesControls *slcontrols);
static void slopes_remove_all_fit_curves(SlopesControls *slcontrols);

static void recalculate_Eit(SlopesControls *slopescontrols, Slopesdata *slopesdata, const Instdata *instdata);
static void slopes_recalculate_beta(SlopesControls *slcontrols, Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata);

static void update_models_slopes(SlopesControls *slcontrols);
static void slopes_zoom_in(Data *data);
static void slopes_zoom_out(Data *data);

static void slopes_show_log_load(Data *data);
static void slopes_show_log_unload(Data *data);


void slopes_remove_fit_results_labels(Data *data)
{
    gwy_selection_clear(data->controls->slopescontrols->selection);
    slopes_remove_all_fit_results_load(&(data->args->slopesdata));
    slopes_remove_all_fit_results_unload(&(data->args->slopesdata));
    slopes_remove_all_fit_labels_load(data->controls->slopescontrols);
    slopes_remove_all_fit_labels_unload(data->controls->slopescontrols);
    slopes_remove_all_fit_curves(data->controls->slopescontrols);
}

static void slopes_remove_all_fit_labels_unload(SlopesControls *slcontrols)
{
    /* remove all labels that result from the powerlaw fit of the unloading curve */

    /*
    entry_clear(slcontrols->from);
    entry_clear(slcontrols->to);
    entry_clear(slcontrols->from_pct_Fmax);
    entry_clear(slcontrols->to_pct_Fmax);
    */

    label_clear(slcontrols->fitinfounload);
    label_clear(slcontrols->hp);
    label_clear(slcontrols->Sunload);
    label_clear(slcontrols->m);
    label_clear(slcontrols->eps);
    label_clear(slcontrols->Aphc);

    label_clear(slcontrols->Hit);
    label_clear(slcontrols->Er);
    label_clear(slcontrols->Eit);
    label_clear(slcontrols->unloadrange);
    label_clear(slcontrols->unloadrange_pct_Fmax);

    label_clear(slcontrols->unloadSSres);
    label_clear(slcontrols->unloadR2);
    label_clear(slcontrols->unloadR2adj);
    label_clear(slcontrols->unloadchi2);
}

static void slopes_remove_all_fit_labels_load(SlopesControls *slcontrols)
{
    /* remove all labels that result from the powerlaw fit of the loading curve */
    /*
    entry_clear(slcontrols->from);
    entry_clear(slcontrols->to);
    entry_clear(slcontrols->from_pct_Fmax);
    entry_clear(slcontrols->to_pct_Fmax);
    */

    label_clear(slcontrols->fitinfoload);
    label_clear(slcontrols->h0);
    label_clear(slcontrols->Sload);
    label_clear(slcontrols->n);
    label_clear(slcontrols->Aphc);

    label_clear(slcontrols->Hit);
    label_clear(slcontrols->Er);
    label_clear(slcontrols->Eit);
    label_clear(slcontrols->loadrange);
    label_clear(slcontrols->loadrange_pct_Fmax);

    label_clear(slcontrols->loadSSres);
    label_clear(slcontrols->loadR2);
    label_clear(slcontrols->loadR2adj);
    label_clear(slcontrols->loadchi2);
}

static void slopes_remove_all_fit_curves(SlopesControls *slcontrols)
{
    gwy_graph_model_remove_curve_by_description(slcontrols->graph_model, "Loading fit");
    gwy_graph_model_remove_curve_by_description(slcontrols->graph_model, "Unloading fit");
}

static void recalculate_Eit(SlopesControls *slopescontrols, Slopesdata *slopesdata, const Instdata *instdata)
{
    gdouble Er;

    Er = slopesdata->Er * 1e9;

    // recalc for changed  nu, nui or Ei
    if (Er == 0) {
        return;
    }

    slopesdata->Eit = calc_Eit(Er, instdata->nu, instdata->Ei, instdata->nui) * 1e-9;
    label_set_gdouble_format(slopescontrols->Eit, slopesdata->Eit, MODUL);
}

static void fixgamma_changed(GtkWidget *button, Data *data)
{
    /*fix gamma 1-fitted, 0-fixed */
    Slopesdata *sldata;
    SlopesControls *slctrls;

    sldata = &(data->args->slopesdata);
    slctrls = data->controls->slopescontrols;

    sldata->loadifixb[2] = 1 - gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));

    slopes_remove_all_fit_labels_load(slctrls);
    slopes_remove_all_fit_results_load(sldata);
    slopes_remove_all_fit_curves(slctrls);

    if (verbose) {
        if (sldata->loadifixb[2]) {
            g_print("Gamma fitted \n");
        }
        else {
            g_print("Gamma fixed \n");
        }
    }
}

void slopes_redraw(Data *data)
{
    GwyGraphCurveModel *cmodel;

    // update loading curve and Fmax and hmax labels
    if (data->args->fddata.hunload && data->args->fddata.hload && data->args->fddata.hhold) {

        gwy_graph_curve_model_set_data(data->controls->slopescontrols->cmodelload,
                                       data->args->fddata.hload->data, data->args->fddata.Fload->data,
                                       data->args->fddata.hload->res);

        gwy_graph_curve_model_set_data(data->controls->slopescontrols->cmodelhold,
                                       data->args->fddata.hhold->data, data->args->fddata.Fhold->data,
                                       data->args->fddata.hhold->res);

        gwy_graph_curve_model_set_data(data->controls->slopescontrols->cmodelunload,
                                       data->args->fddata.hunload->data, data->args->fddata.Funload->data,
                                       data->args->fddata.hunload->res);

        label_set_gdouble_format(data->controls->slopescontrols->hmax, data->args->fddata.hmax, DEPTH);
        label_set_gdouble_format(data->controls->slopescontrols->Fmax, data->args->fddata.Fmax, FORCE);

        cmodel = gwy_graph_model_get_curve_by_description(data->controls->slopescontrols->graph_model, "Loading fit");

        if (cmodel) {
            gwy_graph_curve_model_set_data(data->controls->slopescontrols->cmodelloadfit,
                                           data->args->slopesdata.xloadfit->data, data->args->slopesdata.yloadfit->data,  data->args->slopesdata.xloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(data->controls->slopescontrols->cmodelloadfit);
#endif
        }

        cmodel = gwy_graph_model_get_curve_by_description(data->controls->slopescontrols->graph_model, "Unloading fit");

        if (cmodel) {
            gwy_graph_curve_model_set_data(data->controls->slopescontrols->cmodelunloadfit,
                                           data->args->slopesdata.xunloadfit->data, data->args->slopesdata.yunloadfit->data,  data->args->slopesdata.xunloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(data->controls->slopescontrols->cmodelunloadfit);
#endif
        }
    }
}

static void slopes_recalculate_beta(SlopesControls *slcontrols, Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata)
{
    if (!sldata->has_fit_load || ! sldata->has_fit_unload) {
        return;
    }

    slopes_combine_fit_results(sldata, fddata, instdata);

    label_set_gdouble_format(slcontrols->Hit, sldata->Hit, HARD);
    label_set_gdouble_format(slcontrols->Eit, sldata->Eit, MODUL);
    label_set_gdouble_format(slcontrols->Er, sldata->Er, MODUL);
}

static void slopes_beta_changed(Data *data)
{
    SlopesControls *fcontrols;

    fcontrols = data->controls->slopescontrols;

    if (beta_changed(fcontrols->beta, &(data->args->slopesdata.beta))) {
        slopes_recalculate_beta(data->controls->slopescontrols, &(data->args->slopesdata), &(data->args->fddata), &(data->args->instdata));
    }
}

static void slopes_reset_beta(Data *data)
{
    data->args->slopesdata.beta =  BETA_DEFAULT_SLOPES;
    entry_set_gdouble_format(data->controls->slopescontrols->beta, data->args->slopesdata.beta, "%.03f");
    slopes_recalculate_beta(data->controls->slopescontrols, &(data->args->slopesdata), &(data->args->fddata), &(data->args->instdata));
}

static void slopes_range_from_changed(Data *data)
{
    SlopesControls *fcontrols;
    Slopesdata *fdata;

    fcontrols = data->controls->slopescontrols;
    fdata = &(data->args->slopesdata);

    range_from_changed(fcontrols->from, fcontrols->to, fcontrols->selection, &(fdata->loadfrom), &(fdata->unloadfrom));
}

static void slopes_range_to_changed(Data *data)
{
    SlopesControls *fcontrols;
    Slopesdata *fdata;

    fcontrols = data->controls->slopescontrols;
    fdata = &(data->args->slopesdata);

    range_to_changed(fcontrols->from, fcontrols->to, fcontrols->selection, &(fdata->loadto), &(fdata->unloadto));
}

static gboolean slopes_range_from_pct_Fmax_changed(Data *data)
{
    SlopesControls *rcontrols;
    Slopesdata *rdata;

    rcontrols = data->controls->slopescontrols;
    rdata = &(data->args->slopesdata);

    return range_from_pct_Fmax_changed(rcontrols->from_pct_Fmax, rcontrols->to_pct_Fmax, rcontrols->selection, data->args->fddata.Fmax,
                                       data->args->fddata.hunload, data->args->fddata.Funload,
                                       &(rdata->loadfrom_pct_Fmax), &(rdata->unloadfrom_pct_Fmax), &(rdata->loadto_pct_Fmax), &(rdata->unloadto_pct_Fmax),
                                       &(rdata->Finputload), &(rdata->Finputunload));
}

static gboolean slopes_range_to_pct_Fmax_changed(Data *data)
{
    SlopesControls *rcontrols;
    Slopesdata *rdata;

    rcontrols = data->controls->slopescontrols;
    rdata = &(data->args->slopesdata);

    return range_to_pct_Fmax_changed(rcontrols->from_pct_Fmax, rcontrols->to_pct_Fmax, rcontrols->selection, data->args->fddata.Fmax,
                                     data->args->fddata.hunload, data->args->fddata.Funload,
                                     &(rdata->loadfrom_pct_Fmax), &(rdata->unloadfrom_pct_Fmax), &(rdata->loadto_pct_Fmax), &(rdata->unloadto_pct_Fmax),
                                     &(rdata->Finputload), &(rdata->Finputunload));
}

static void slopes_graph_selected(GwySelection *selection, gint i, Data *data)
{
    gdouble range[2];
    gint istart, iend, ndata;
    gdouble Fmax;

    g_return_if_fail(i <= 0);

    //take data from selection and put them in the "from" and "to" values for both hp and S (at this moment, we don't know yet which fit will be chosen), create labels only once
    //SlopesControls has entries from and to, which are filled here
    // and labels hp,unloadrange, where the values of the ranges that were actually used in the corresponding fits  are shown
    //Slopesdata has loadfrom,loadto, unloadfrom,unloadto which are the data from the entries, i.e. from the selection
    //Slopesdata; also has loadrange, unloadrange where the ranges that were actually used in the fit

    gwy_selection_get_object(selection, 0, range);

    data->args->slopesdata.unloadfrom = MIN(range[0], range[1]);
    data->args->slopesdata.unloadto = MAX(range[0], range[1]);

    data->args->slopesdata.loadfrom = MIN(range[0], range[1]);
    data->args->slopesdata.loadto = MAX(range[0], range[1]);

    entry_set_gdouble_format(data->controls->slopescontrols->from, data->args->slopesdata.unloadfrom, DEPTH);
    entry_set_gdouble_format(data->controls->slopescontrols->to, data->args->slopesdata.unloadto, DEPTH);

    Fmax = data->args->fddata.Fmax;
    range_to_indices(data->args->slopesdata.unloadfrom, data->args->slopesdata.unloadto, data->args->fddata.hunload, TRUE, &istart, &iend, &ndata);
    data->args->slopesdata.loadfrom_pct_Fmax = data->args->fddata.Funload->data[iend] / Fmax * 100;
    data->args->slopesdata.loadto_pct_Fmax = data->args->fddata.Funload->data[istart] / Fmax * 100;
    data->args->slopesdata.unloadfrom_pct_Fmax = data->args->fddata.Funload->data[iend] / Fmax * 100;
    data->args->slopesdata.unloadto_pct_Fmax = data->args->fddata.Funload->data[istart] / Fmax * 100;

    entry_set_gdouble_format(data->controls->slopescontrols->from_pct_Fmax, data->args->slopesdata.unloadfrom_pct_Fmax, DEPTH);
    entry_set_gdouble_format(data->controls->slopescontrols->to_pct_Fmax, data->args->slopesdata.unloadto_pct_Fmax, DEPTH);

    data->args->slopesdata.Finputload = FALSE;
    data->args->slopesdata.Finputunload = FALSE;
}

void slopes_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (slopes_has_results(&(args->slopesdata))) {
	buffer = slopes_export_data(&(args->slopesdata), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

void slopes_fit_and_update_load(Data *data)
{
    gchar str[300];
    gint i;
    gint istart, iend, ndata;
    gboolean fit = TRUE;
    GwyDataLine *x, *y;
    gint fitinfo;
    gdouble range[2];

    Slopesdata *sldata;
    SlopesControls *slcontrols;
    FDdata *fddata;
    Instdata *instdata;
    GdkColor color;
    GwyGraphCurveModel *cmodel;

    gdk_color_parse("red", &color);

    //create pointers for easier reading
    sldata = &(data->args->slopesdata);
    fddata = &(data->args->fddata);
    instdata = &(data->args->instdata);
    slcontrols = data->controls->slopescontrols;

    //if there was a previous fit, invalidate previous fitted curve
    if (sldata->has_fit_load) {
        g_object_unref(sldata->xloadfit);
        g_object_unref(sldata->yloadfit);
    }

    if (verbose) {
        g_print("from %g to %g \n", sldata->loadfrom, sldata->loadto);
    }

    //check range
    if (sldata->loadfrom == sldata->loadto) {
        fit = FALSE;
    }

    //find appropriate data, which should be fitted
    //default regime is governed by depth, can be switched to load regime by Finput (if the pct_Fmax ranges were used)
    if (sldata->Finputload) {
        //		printf("Finputload  from %g to %g %%Fmax, \n", sldata->loadfrom_pct_Fmax, sldata->loadto_pct_Fmax);
        double oldfromF, oldtoF;
        oldfromF = sldata->loadfrom_pct_Fmax;
        oldtoF = sldata->loadto_pct_Fmax;
        //		printf(" corresponding to %g %g \n", sldata->loadfrom_pct_Fmax*0.01*fddata->Fmax, sldata->loadto_pct_Fmax*0.01*fddata->Fmax);
        range_to_indices(sldata->loadfrom_pct_Fmax * 0.01 * fddata->Fmax, sldata->loadto_pct_Fmax * 0.01 * fddata->Fmax,
                         fddata->Fload, FALSE, &istart, &iend, &ndata);
        //		printf(" istart %d iend %d \n", istart, iend);
        //	printf(" hstart %g hend %g \n", fddata->hload->data[istart], fddata->hload->data[iend]);
        //printf(" Fstart %g Fend %g \n", fddata->Fload->data[istart], fddata->Fload->data[iend]);
        /* force to use recalculated loads suing loading curve (not default unloading) */
        sldata->loadfrom = fddata->hload->data[istart];
        sldata->loadto = fddata->hload->data[iend];

        /* force to use recalculated loads suing loading curve (not default unloading) */
        range[0] = sldata->loadfrom;
        range[1] = sldata->loadto;
        gwy_selection_set_object(data->controls->slopescontrols->selection, 0, range);
        sldata->loadfrom_pct_Fmax = oldfromF;
        sldata->loadto_pct_Fmax = oldtoF;
    }
    else {
        range_to_indices(sldata->loadfrom, sldata->loadto, fddata->hload, FALSE, &istart, &iend, &ndata);
        /* force to use recalculated loads suing loading curve (not default unloading) */
        sldata->loadfrom_pct_Fmax = fddata->Fload->data[istart] / fddata->Fmax * 100;
        sldata->loadto_pct_Fmax = fddata->Fload->data[iend] / fddata->Fmax * 100;
    }

    // if there are not enough data in the selected range, no fitting can be performed
    ndata = MIN(ndata, fddata->hload->res - istart);

    if (ndata < 2) {
        fit = FALSE;
    }

    sldata->nfitloaddata = ndata;

    if (fit) {
        //get data to be fitted
        x = gwy_data_line_part_extract(fddata->hload, istart, ndata);
        y = gwy_data_line_part_extract(fddata->Fload, istart, ndata);

        //fit the data and evaluate
        fitinfo = slopes_fit_load(x, y, sldata, fddata, instdata, NULL);

        if (sldata->has_fit_unload) {
            slopes_combine_fit_results(sldata, fddata, instdata);
        }

        if (fitinfo >= 1e6) {
            sldata->infologload = fitinfo / 1000000;
            fitinfo = fitinfo % 1000000;
        }

        if (fitinfo >= 5) {
            gtk_label_set_text(GTK_LABEL(slcontrols->fitinfoload), "Warning Load fit: Questionable results \n or fatal errors in fitting!");
            gtk_widget_modify_fg(slcontrols->fitinfoload, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(slcontrols->fitinfoload), TRUE);
        }
        else if (fitinfo == 4) {
            gtk_label_set_text(GTK_LABEL(slcontrols->fitinfoload), "Warning Load fit: Iteration limit reached!");
            gtk_widget_modify_fg(slcontrols->fitinfoload, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(slcontrols->fitinfoload), TRUE);
        }
        else {
            label_clear(slcontrols->fitinfoload);
        }

        //update range labels
        entry_set_gdouble_format(slcontrols->from, sldata->loadfrom, DEPTH);
        entry_set_gdouble_format(slcontrols->to, sldata->loadto, DEPTH);

        entry_set_gdouble_format(slcontrols->from_pct_Fmax, sldata->loadfrom_pct_Fmax, DEPTH);
        entry_set_gdouble_format(slcontrols->to_pct_Fmax, sldata->loadto_pct_Fmax, DEPTH);

        //update labels
        label_set_gdouble_format(slcontrols->h0, sldata->h0, DEPTH);
        label_set_gdouble_format(slcontrols->n, sldata->n, NUMBER);
        label_set_gdouble_format(slcontrols->Sload, sldata->Sload, SLOPE);
        label_set_gdouble_format(slcontrols->Aphc, sldata->Aphc, AREA);

        label_set_gdouble_format(slcontrols->loadSSres, sldata->loadSSres, "%.4g");
        label_set_gdouble_format(slcontrols->loadR2, sldata->loadR2, "%.4g");
        label_set_gdouble_format(slcontrols->loadR2adj, sldata->loadR2adj, "%.4g");
        label_set_gdouble_format(slcontrols->loadchi2, sldata->loadchi2, "%.4g");

        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", sldata->loadfrom, sldata->loadto);
        gtk_label_set_text(GTK_LABEL(slcontrols->loadrange), str);
        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", sldata->loadfrom_pct_Fmax, sldata->loadto_pct_Fmax);
        gtk_label_set_text(GTK_LABEL(slcontrols->loadrange_pct_Fmax), str);

        label_set_gdouble_format(slcontrols->Hit, sldata->Hit, HARD);
        label_set_gdouble_format(slcontrols->Eit, sldata->Eit, MODUL);
        label_set_gdouble_format(slcontrols->Er, sldata->Er, MODUL);

        //calculate curve of fit
        sldata->xloadfit = gwy_data_line_duplicate(x);
        sldata->yloadfit = gwy_data_line_new_alike(y, FALSE);

        for (i = 0; i < sldata->xloadfit->res; i++) {
            sldata->yloadfit->data[i] =  sldata->gamma * pow(sldata->xloadfit->data[i] - sldata->h0, sldata->n);
        }

        // clean up
        g_object_unref(x);
        g_object_unref(y);

        //update fitted curve, if necessary create it
        cmodel = gwy_graph_model_get_curve_by_description(slcontrols->graph_model, "Loading fit");

        if (cmodel) {
            gwy_graph_curve_model_set_data(slcontrols->cmodelloadfit, sldata->xloadfit->data, sldata->yloadfit->data, sldata->xloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(slcontrols->cmodelloadfit);
#endif
        }
        else {
            gwy_graph_curve_model_set_data(slcontrols->cmodelloadfit, sldata->xloadfit->data, sldata->yloadfit->data, sldata->xloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(slcontrols->cmodelloadfit);
#endif
            gwy_graph_model_add_curve(slcontrols->graph_model, slcontrols->cmodelloadfit);
            g_object_set(slcontrols->cmodelloadfit,
                         "mode", GWY_GRAPH_CURVE_LINE,
                         "description", "Loading fit",
                         "color", CURVE_COLOR_FIT,
                         "line-style", GDK_LINE_SOLID,
                         "line-width", 3,
                         NULL);
        }
    }
    else {
        slopes_remove_all_fit_labels_load(slcontrols);
        slopes_remove_all_fit_results_load(sldata);
        sldata->has_fit_load = FALSE;
        gwy_graph_model_remove_curve_by_description(slcontrols->graph_model, "Loading fit");
    }

    // log of fitting procedure is available
    if (sldata->has_fit_load) {
        gtk_widget_set_sensitive(slcontrols->button_log_load, TRUE);
    }
    else {
        gtk_widget_set_sensitive(slcontrols->button_log_load, FALSE);
    }

    if (sldata->has_fit_load &&  sldata->has_fit_unload) {
        if (data->args->slopesunc.Ec != NULL || data->args->slopesunc.Hc != NULL || data->args->slopesunc.Ac != NULL) {
            slopes_unc_close_window(data);
        }

        gtk_widget_set_sensitive(slcontrols->button_unc, TRUE);
    }
    else {
        gtk_widget_set_sensitive(slcontrols->button_unc, FALSE);
    }
}

void slopes_fit_and_update_unload(Data *data)
{
    gchar str[300];
    gint  i;
    gint istart, iend, ndata;
    gint fitinfo;
    gboolean fit = TRUE;
    GwyDataLine *x, *y;

    Slopesdata *sldata;
    SlopesControls *slcontrols;
    FDdata *fddata;
    Instdata *instdata;
    GdkColor color;
    GwyGraphCurveModel *cmodel;

    gdk_color_parse("red", &color);

    //create pointers for easier reading
    sldata = &(data->args->slopesdata);
    fddata = &(data->args->fddata);
    instdata = &(data->args->instdata);
    slcontrols = data->controls->slopescontrols;

    //if there was a previous fit, free the old data and create new data
    if (sldata->has_fit_unload) {
        g_object_unref(sldata->xunloadfit);
        g_object_unref(sldata->yunloadfit);
    }

    if (verbose) {
        g_print("from %g to %g \n", sldata->unloadfrom, sldata->unloadto);
    }

    // check range
    if (sldata->unloadfrom == sldata->unloadto) {
        fit = FALSE;
    }

    //find appropriate data, which should be fitted
    //default regime is governed by depth, can be switched to load regime by Finput (if the pct_Fmax ranges were used)
    if (sldata->Finputunload) {
        range_to_indices(sldata->unloadfrom_pct_Fmax * 0.01 * fddata->Fmax, sldata->unloadto_pct_Fmax * 0.01 * fddata->Fmax,
                         fddata->Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(sldata->unloadfrom, sldata->unloadto, fddata->hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range,no fitting can be performed
    ndata = MIN(ndata, fddata->hunload->res - istart);

    if (ndata < 2) {
        fit = FALSE;
    }

    sldata->nfitunloaddata = ndata;

    if (fit) {
        //get data to be fitted
        x = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
        y = gwy_data_line_part_extract(fddata->Funload, istart, ndata);

        //fit the data and evaluate
        fitinfo = slopes_fit_unload(x, y, sldata, fddata, instdata, NULL);

        if (sldata->has_fit_load) {
            slopes_combine_fit_results(sldata, fddata, instdata);
        }

        if (fitinfo >= 1e6) {
            sldata->infologunload = fitinfo / 1000000;
            fitinfo = fitinfo % 1000000;
        }

        if (fitinfo >= 5) {
            gtk_label_set_text(GTK_LABEL(slcontrols->fitinfounload), "Warning Unload fit: Questionable results \n or fatal errors in fitting!");
            gtk_widget_modify_fg(slcontrols->fitinfounload, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(slcontrols->fitinfounload), TRUE);
        }
        else if (fitinfo == 4) {
            gtk_label_set_text(GTK_LABEL(slcontrols->fitinfounload), "Warning Unload fit: Iteration limit reached!");
            gtk_widget_modify_fg(slcontrols->fitinfounload, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(slcontrols->fitinfounload), TRUE);
        }
        else if (sldata->m > 2) {
            gtk_label_set_text(GTK_LABEL(slcontrols->fitinfounload), "Warning Unload fit: Unphysical fit results!");
            gtk_widget_modify_fg(slcontrols->fitinfounload, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(slcontrols->fitinfounload), TRUE);
        }
        else {
            label_clear(slcontrols->fitinfounload);
        }

        //update range labels
        entry_set_gdouble_format(slcontrols->from, sldata->unloadfrom, DEPTH);
        entry_set_gdouble_format(slcontrols->to, sldata->unloadto, DEPTH);

        entry_set_gdouble_format(slcontrols->from_pct_Fmax, sldata->unloadfrom_pct_Fmax, DEPTH);
        entry_set_gdouble_format(slcontrols->to_pct_Fmax, sldata->unloadto_pct_Fmax, DEPTH);

        //update labels
        label_set_gdouble_format(slcontrols->hp, sldata->hp, DEPTH);
        label_set_gdouble_format(slcontrols->m, sldata->m, NUMBER);
        label_set_gdouble_format(slcontrols->eps, sldata->eps, NUMBER);
        label_set_gdouble_format(slcontrols->Sunload, sldata->Sunload, SLOPE);
        label_set_gdouble_format(slcontrols->Aphc, sldata->Aphc, AREA);

        label_set_gdouble_format(slcontrols->unloadSSres, sldata->unloadSSres, "%.4g");
        label_set_gdouble_format(slcontrols->unloadR2, sldata->unloadR2, "%.4g");
        label_set_gdouble_format(slcontrols->unloadR2adj, sldata->unloadR2adj, "%.4g");
        label_set_gdouble_format(slcontrols->unloadchi2, sldata->unloadchi2, "%.4g");

        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", sldata->unloadfrom, sldata->unloadto);
        gtk_label_set_text(GTK_LABEL(slcontrols->unloadrange), str);
        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", sldata->unloadfrom_pct_Fmax, sldata->unloadto_pct_Fmax);
        gtk_label_set_text(GTK_LABEL(slcontrols->unloadrange_pct_Fmax), str);

        label_set_gdouble_format(slcontrols->Hit, sldata->Hit, HARD);
        label_set_gdouble_format(slcontrols->Eit, sldata->Eit, MODUL);
        label_set_gdouble_format(slcontrols->Er, sldata->Er, MODUL);

        //calculate curve of fit
        sldata->xunloadfit = gwy_data_line_duplicate(x);
        sldata->yunloadfit = gwy_data_line_new_alike(y, FALSE);

        for (i = 0; i < sldata->xunloadfit->res; i++) {
            sldata->yunloadfit->data[i] =  sldata->alpha * pow(sldata->xunloadfit->data[i] - sldata->hp, sldata->m);
        }

        // clean up
        g_object_unref(x);
        g_object_unref(y);

        //update fitted curve, if necessary create it
        cmodel = gwy_graph_model_get_curve_by_description(slcontrols->graph_model, "Unloading fit");

        if (cmodel) {
            gwy_graph_curve_model_set_data(slcontrols->cmodelunloadfit, sldata->xunloadfit->data, sldata->yunloadfit->data, sldata->xunloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(slcontrols->cmodelunloadfit);
#endif
	}
        else {
            gwy_graph_curve_model_set_data(slcontrols->cmodelunloadfit, sldata->xunloadfit->data, sldata->yunloadfit->data, sldata->xunloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(slcontrols->cmodelunloadfit);
#endif
            gwy_graph_model_add_curve(slcontrols->graph_model, slcontrols->cmodelunloadfit);
            g_object_set(slcontrols->cmodelunloadfit,
                         "mode", GWY_GRAPH_CURVE_LINE,
                         "description", "Unloading fit",
                         "color", CURVE_COLOR_FIT2,
                         "line-style", GDK_LINE_SOLID,
                         "line-width", 3,
                         NULL);
        }
    }
    else {
        slopes_remove_all_fit_labels_unload(slcontrols);
        slopes_remove_all_fit_results_unload(sldata);
        sldata->has_fit_unload = FALSE;
        gwy_graph_model_remove_curve_by_description(slcontrols->graph_model, "Unloading fit");
    }

    // log of fitting procedure is available
    if (sldata->has_fit_unload) {
        gtk_widget_set_sensitive(slcontrols->button_log_unload, TRUE);
    }
    else {
        gtk_widget_set_sensitive(slcontrols->button_log_unload, FALSE);
    }

    if (sldata->has_fit_load &&  sldata->has_fit_unload) {
        if (data->args->slopesunc.Ec != NULL || data->args->slopesunc.Hc != NULL || data->args->slopesunc.Ac != NULL) {
            slopes_unc_close_window(data);
        }

        gtk_widget_set_sensitive(slcontrols->button_unc, TRUE);
    }
    else {
        gtk_widget_set_sensitive(slcontrols->button_unc, FALSE);
    }
}

static void slopes_show_log_load(Data *data)
{
    gint infolog;

    infolog = data->args->slopesdata.infologload;

    /* check if there are logs */
    switch (infolog) {
    case 0:
	show_log(data->controls->slopescontrols->load_log_error_window, data->args->slopesdata.logfnmload, ".err", "Fit error log loading");
	show_log(data->controls->slopescontrols->load_log_report_window, data->args->slopesdata.logfnmload, ".rpt", "Fit report log loading");
	break;
	
    case 1:
	show_log(data->controls->slopescontrols->load_log_error_window,  data->args->slopesdata.logfnmload, ".err", "Fit error log loading");
	show_warning("Fit report log could not be created.");
	break;

    case 2:
	show_log(data->controls->slopescontrols->load_log_report_window, data->args->slopesdata.logfnmload, ".rpt", "Fit report log loading");
	show_warning("Fit error log could not be created.");
	break;
	
    case 3:
	show_warning("Fit error log loading could not be created.");
	show_warning("Fit report log loading could not be created.");
	break;

    default:
	g_printerr("Should not get here.\n");
	break;
    }
}

static void slopes_show_log_unload(Data *data)
{
    gint infolog;

    infolog = data->args->slopesdata.infologunload;

    /* check if there are logs */
    switch (infolog) {
    case 0:
	show_log(data->controls->slopescontrols->unload_log_error_window, data->args->slopesdata.logfnmunload, ".err", "Fit error log unloading");
	show_log(data->controls->slopescontrols->unload_log_report_window, data->args->slopesdata.logfnmunload, ".rpt", "Fit report log unloading");
	break;
	
    case 1:
	show_log(data->controls->slopescontrols->unload_log_error_window,  data->args->slopesdata.logfnmunload, ".err", "Fit error log unloading");
	show_warning("Fit report log could not be created.");
	break;

    case 2:
	show_log(data->controls->slopescontrols->unload_log_report_window, data->args->slopesdata.logfnmunload, ".rpt", "Fit report log unloading");
	show_warning("Fit error log could not be created.");
	break;
	
    case 3:
	show_warning("Fit error log unloading could not be created.");
	show_warning("Fit report log unloading could not be created.");
	break;

    default:
	g_printerr("Should not get here.\n");
	break;
    }
}

GtkWidget *tool_slopes_create(Data *data)
{
    gchar str[300];
    GwyGraphArea *area;
    Args *args;
    SlopesControls *slcontrols;
    gint i;

    args = data->args;
    slcontrols = data->controls->slopescontrols;

    // GraphArea, three curves can be displayed: loading F-d curve, the linear (hp) fit and the powerlaw (S) fit

    slcontrols->graph_model = gwy_graph_model_new();
    g_object_set(slcontrols->graph_model, "title", "Unload", NULL);
    gwy_graph_model_set_axis_label(slcontrols->graph_model, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(slcontrols->graph_model, GTK_POS_LEFT, "F / mN");

    slcontrols->zoomstack = NULL;
    slcontrols->zoomstack = g_slist_prepend(slcontrols->zoomstack, slcontrols->graph_model);

    slcontrols->hbox_graph = gtk_hbox_new(FALSE, CTRLS_SPACING);

    slcontrols->graph = gwy_graph_new(slcontrols->graph_model);
    /* gtk_widget_set_size_request(GTK_WIDGET(slcontrols->graph), 600, 400); */

    gwy_graph_enable_user_input(GWY_GRAPH(slcontrols->graph), FALSE);
    gwy_graph_set_status(GWY_GRAPH(slcontrols->graph), GWY_GRAPH_STATUS_XSEL);

    area = GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(slcontrols->graph)));
    slcontrols->selection = gwy_graph_area_get_selection(area, GWY_GRAPH_STATUS_XSEL);
    gwy_selection_set_max_objects(slcontrols->selection, 1);
    g_signal_connect(slcontrols->selection, "changed", G_CALLBACK(slopes_graph_selected), data);

    slcontrols->cmodelload = gwy_graph_curve_model_new();
    slcontrols->cmodelhold = gwy_graph_curve_model_new();
    slcontrols->cmodelunload = gwy_graph_curve_model_new();
    slcontrols->cmodelloadfit = gwy_graph_curve_model_new();
    slcontrols->cmodelunloadfit = gwy_graph_curve_model_new();

    g_object_set(slcontrols->cmodelload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Loading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_LOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(slcontrols->cmodelhold,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Hold",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_HOLD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(slcontrols->cmodelunload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Unloading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_UNLOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    gwy_graph_model_add_curve(slcontrols->graph_model, slcontrols->cmodelload);
    gwy_graph_model_add_curve(slcontrols->graph_model, slcontrols->cmodelhold);
    gwy_graph_model_add_curve(slcontrols->graph_model, slcontrols->cmodelunload);


    // Table of input data and results

    slcontrols->ctrls = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* Info section */

    slcontrols->frame_info = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(slcontrols->frame_info), GTK_WIDGET(gwy_label_new_header("Info")));
    slcontrols->table_info = gtk_table_new(2, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(slcontrols->table_info), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(slcontrols->table_info), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(slcontrols->table_info), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(slcontrols->frame_info), slcontrols->table_info);

    i = 0;

    /* h_max */

    gtk_table_attach(GTK_TABLE(slcontrols->table_info), label_new_with_markup_left("h<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    slcontrols->hmax = gtk_label_new(NULL);
    label_set_gdouble_format(slcontrols->hmax, args->fddata.hmax, DEPTH);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->hmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_info), slcontrols->hmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(slcontrols->table_info), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /* F_max */

    gtk_table_attach(GTK_TABLE(slcontrols->table_info), label_new_with_markup_left("F<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    slcontrols->Fmax = gtk_label_new(NULL);
    label_set_gdouble_format(slcontrols->Fmax, args->fddata.Fmax, FORCE);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_info), slcontrols->Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_info), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    /* Parameters section */

    slcontrols->frame_parameters = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(slcontrols->frame_parameters), GTK_WIDGET(gwy_label_new_header("Parameters")));
    slcontrols->table_parameters = gtk_table_new(5, 4, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(slcontrols->table_parameters), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(slcontrols->table_parameters), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(slcontrols->table_parameters), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(slcontrols->frame_parameters), slcontrols->table_parameters);
    i = 0;

    gtk_table_attach(GTK_TABLE(slcontrols->table_parameters), label_new_with_markup_left("Fit range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->from = gtk_entry_new();
    g_object_set_data(G_OBJECT(slcontrols->from), "id", (gpointer)"from");
    gtk_entry_set_width_chars(GTK_ENTRY(slcontrols->from), 8);
    gtk_table_attach(GTK_TABLE(slcontrols->table_parameters), slcontrols->from, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(slcontrols->from, TRUE);
    g_signal_connect_swapped(slcontrols->from, "activate", G_CALLBACK(slopes_range_from_changed), data);
    /*	g_signal_connect_swapped(slcontrols->from, "focus-out-event", G_CALLBACK(range_from_changed),data); */

    slcontrols->to = gtk_entry_new();
    g_object_set_data(G_OBJECT(slcontrols->to), "id", (gpointer)"to");
    gtk_entry_set_width_chars(GTK_ENTRY(slcontrols->to), 8);
    gtk_table_attach(GTK_TABLE(slcontrols->table_parameters), slcontrols->to, 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(slcontrols->to, TRUE);
    g_signal_connect_swapped(slcontrols->to, "activate", G_CALLBACK(slopes_range_to_changed), data);
    /*	g_signal_connect_swapped(slcontrols->to, "focus-out-event", G_CALLBACK(range_to_changed),data); */

    gtk_table_attach(GTK_TABLE(slcontrols->table_parameters), label_new_with_markup_left("nm"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    slcontrols->from_pct_Fmax = gtk_entry_new();
    g_object_set_data(G_OBJECT(slcontrols->from_pct_Fmax), "id", (gpointer)"from_pct_Fmax");
    gtk_entry_set_width_chars(GTK_ENTRY(slcontrols->from_pct_Fmax), 8);
    gtk_table_attach(GTK_TABLE(slcontrols->table_parameters), slcontrols->from_pct_Fmax, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(slcontrols->from_pct_Fmax, TRUE);
    g_signal_connect_swapped(slcontrols->from_pct_Fmax, "activate", G_CALLBACK(slopes_range_from_pct_Fmax_changed), data);
    /*	g_signal_connect_swapped(slcontrols->from_pct_Fmax, "focus-out-event", G_CALLBACK(range_from_pct_Fmax_changed),data);  */

    slcontrols->to_pct_Fmax = gtk_entry_new();
    g_object_set_data(G_OBJECT(slcontrols->to_pct_Fmax), "id", (gpointer)"to_pct_Fmax");
    gtk_entry_set_width_chars(GTK_ENTRY(slcontrols->to_pct_Fmax), 8);
    gtk_table_attach(GTK_TABLE(slcontrols->table_parameters), slcontrols->to_pct_Fmax, 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(slcontrols->to_pct_Fmax, TRUE);
    g_signal_connect_swapped(slcontrols->to_pct_Fmax, "activate", G_CALLBACK(slopes_range_to_pct_Fmax_changed), data);
    /*g_signal_connect_swapped(slcontrols->to_pct_Fmax, "focus-out-event", G_CALLBACK(range_to_pct_Fmax_changed),data);  */

    gtk_table_attach(GTK_TABLE(slcontrols->table_parameters), label_new_with_markup_left("% F<sub>max</sub>"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(slcontrols->table_parameters), label_new_with_markup_left("β"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->beta = gtk_entry_new();
    //	g_object_set_data(G_OBJECT(slcontrols->beta), "id", (gpointer)"beta");
    gtk_entry_set_width_chars(GTK_ENTRY(slcontrols->beta), 8);
    entry_set_gdouble_format(slcontrols->beta, data->args->slopesdata.beta, "%.03f");
    gwy_widget_set_activate_on_unfocus(slcontrols->beta, TRUE);
    gtk_table_attach(GTK_TABLE(slcontrols->table_parameters), slcontrols->beta, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(slcontrols->beta, "activate", G_CALLBACK(slopes_beta_changed), data);

    slcontrols->reset_beta = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(slcontrols->reset_beta), gtk_image_new_from_stock(GTK_STOCK_REVERT_TO_SAVED, GTK_ICON_SIZE_BUTTON));
    gtk_table_attach(GTK_TABLE(slcontrols->table_parameters), slcontrols->reset_beta, 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(slcontrols->reset_beta, "clicked", G_CALLBACK(slopes_reset_beta), data);
    i++;

    g_snprintf(str, sizeof(str), "Fix exponent at %g", N_FIX_SLOPES);
    slcontrols->fixgamma = gtk_check_button_new_with_label(str);
    //slcontrols->fixgamma = gtk_check_button_new_with_label("Fix exponent at 1.5");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(slcontrols->fixgamma), FALSE);
    gtk_table_attach(GTK_TABLE(slcontrols->table_parameters), slcontrols->fixgamma, 0, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(slcontrols->fixgamma, "toggled", G_CALLBACK(fixgamma_changed), data);


    /* Fit buttons */

    slcontrols->hbox_fitbuttons = gtk_hbox_new(TRUE, 0);

    slcontrols->button_fit_load = gtk_button_new_with_mnemonic(gwy_sgettext("_Fit loading"));
    gtk_button_set_image(GTK_BUTTON(slcontrols->button_fit_load), gtk_image_new_from_stock(GTK_STOCK_EXECUTE, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(slcontrols->hbox_fitbuttons), slcontrols->button_fit_load, TRUE, TRUE, 0);
    g_signal_connect_swapped(slcontrols->button_fit_load, "clicked", G_CALLBACK(slopes_fit_and_update_load), data);

    slcontrols->button_fit_unload = gtk_button_new_with_mnemonic(gwy_sgettext("verb|F_it unloading"));
    gtk_button_set_image(GTK_BUTTON(slcontrols->button_fit_unload), gtk_image_new_from_stock(GTK_STOCK_EXECUTE, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(slcontrols->hbox_fitbuttons), slcontrols->button_fit_unload, TRUE, TRUE, 0);
    g_signal_connect_swapped(slcontrols->button_fit_unload, "clicked", G_CALLBACK(slopes_fit_and_update_unload), data);


    /* Results section */

    slcontrols->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(slcontrols->frame_results), GTK_WIDGET(gwy_label_new_header("Results")));

    slcontrols->scrolledwindow_results = gtk_scrolled_window_new(NULL, NULL);
    /* TODO: resolve table column widths so that horizontal scroll policy can become NEVER */
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(slcontrols->scrolledwindow_results), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(slcontrols->frame_results), slcontrols->scrolledwindow_results);

    slcontrols->table_results = gtk_table_new(20, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(slcontrols->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(slcontrols->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(slcontrols->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(slcontrols->scrolledwindow_results), slcontrols->table_results);
    i = 0;

    slcontrols->fitinfoload = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->fitinfoload, 0, 3, i, i + 2, GTK_FILL, GTK_FILL, 0, 0);
    i += 2;


    /*    h0    */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("h<sub>0</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->h0 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->h0), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->h0, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /*    n    */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("n"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->n = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->n), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->n, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* SSres */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("SS<sub>res</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->loadSSres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->loadSSres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->loadSSres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* R2 */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("R<sup>2</sup>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->loadR2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->loadR2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->loadR2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* R2adj */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("R<sup>2</sup><sub>adj</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->loadR2adj = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->loadR2adj), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->loadR2adj, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* chi2 */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("χ<sup>2</sup> fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->loadchi2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->loadchi2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->loadchi2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;


    /* unload */

    slcontrols->fitinfounload = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->fitinfounload, 0, 3, i, i + 2, GTK_FILL, GTK_FILL, 0, 0);
    i += 2;

    /*    hp    */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("h<sub>p</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->hp = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->hp), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->hp, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /*    m    */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("m"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->m = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->m), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->m, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* SSres */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("SS<sub>res</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->unloadSSres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->unloadSSres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->unloadSSres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* R2 */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("R<sup>2</sup>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->unloadR2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->unloadR2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->unloadR2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* R2adj */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("R<sup>2</sup><sub>adj</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->unloadR2adj = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->unloadR2adj), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->unloadR2adj, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* chi2 */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("χ<sup>2</sup> fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->unloadchi2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->unloadchi2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->unloadchi2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /*    eps    */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("ε"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->eps = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->eps), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->eps, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /*    Sload    */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("S<sub>load</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->Sload = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->Sload), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->Sload, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("mN/nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /*    Sunload    */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("S<sub>unload</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->Sunload = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->Sunload), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->Sunload, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("mN/nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /*    Aphc    */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("A<sub>p</sub>(h<sub>c</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->Aphc = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->Aphc), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->Aphc, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("nm<sup>2</sup>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /*    HIT    */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("H<sub>IT</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->Hit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->Hit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->Hit, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("MPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /*    Er    */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("E<sub>r</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->Er = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->Er), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->Er, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /*    EIT    */
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("E<sub>IT</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->Eit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->Eit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->Eit, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("Fitted loading range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->loadrange = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->loadrange), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->loadrange, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    slcontrols->loadrange_pct_Fmax = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->loadrange_pct_Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->loadrange_pct_Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("% F<sub>max</sub>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("Fitted unloading range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    slcontrols->unloadrange = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->unloadrange), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->unloadrange, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    slcontrols->unloadrange_pct_Fmax = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(slcontrols->unloadrange_pct_Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), slcontrols->unloadrange_pct_Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(slcontrols->table_results), label_new_with_markup_left("% F<sub>max</sub>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);


    /* Log buttons */
    slcontrols->hbox_logbuttons = gtk_hbox_new(TRUE, 0);

    slcontrols->button_log_load = gtk_button_new_with_mnemonic("Show _log load");
    gtk_button_set_image(GTK_BUTTON(slcontrols->button_log_load), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(slcontrols->button_log_load, "clicked", G_CALLBACK(slopes_show_log_load), data);
    gtk_widget_set_sensitive(slcontrols->button_log_load, FALSE);
    slcontrols->load_log_error_window = NULL;
    slcontrols->load_log_report_window = NULL;
    gtk_box_pack_start(GTK_BOX(slcontrols->hbox_logbuttons), slcontrols->button_log_load, TRUE, TRUE, 0);

    slcontrols->button_log_unload = gtk_button_new_with_mnemonic("Show lo_g unload");
    gtk_button_set_image(GTK_BUTTON(slcontrols->button_log_unload), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(slcontrols->button_log_unload, "clicked", G_CALLBACK(slopes_show_log_unload), data);
    gtk_widget_set_sensitive(slcontrols->button_log_unload, FALSE);
    slcontrols->unload_log_error_window = NULL;
    slcontrols->unload_log_report_window = NULL;
    gtk_box_pack_start(GTK_BOX(slcontrols->hbox_logbuttons), slcontrols->button_log_unload, TRUE, TRUE, 0);


    // Uncertainties Button
    slcontrols->button_unc = gtk_button_new_with_mnemonic("_Uncertainties");
    gtk_button_set_image(GTK_BUTTON(slcontrols->button_unc), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(slcontrols->button_unc, "clicked", G_CALLBACK(slopes_uncertainty), data);
    gtk_widget_set_sensitive(slcontrols->button_unc, FALSE);


    // Save Button

    slcontrols->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(slcontrols->button_save, "clicked", G_CALLBACK(slopes_save_data), args);

    gtk_box_pack_start(GTK_BOX(slcontrols->ctrls), slcontrols->frame_info, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(slcontrols->ctrls), slcontrols->frame_parameters, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(slcontrols->ctrls), slcontrols->hbox_fitbuttons, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(slcontrols->ctrls), slcontrols->frame_results, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(slcontrols->ctrls), slcontrols->hbox_logbuttons, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(slcontrols->ctrls), slcontrols->button_unc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(slcontrols->ctrls), slcontrols->button_save, FALSE, FALSE, 0);

    /* Zooming */

    slcontrols->vbox_zoom_buttons = gtk_vbox_new(TRUE, 0);
    slcontrols->vbox_zoom_outer = gtk_vbox_new(FALSE, 0);

    slcontrols->button_zoom_in = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(slcontrols->button_zoom_in), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(slcontrols->vbox_zoom_buttons), slcontrols->button_zoom_in, FALSE, FALSE, 0);
    g_signal_connect_swapped(slcontrols->button_zoom_in, "clicked", G_CALLBACK(slopes_zoom_in), data);

    slcontrols->button_zoom_out = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(slcontrols->button_zoom_out), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(slcontrols->vbox_zoom_buttons), slcontrols->button_zoom_out, FALSE, FALSE, 0);
    g_signal_connect_swapped(slcontrols->button_zoom_out, "clicked", G_CALLBACK(slopes_zoom_out), data);

    slcontrols->button_zoom_restore = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(slcontrols->button_zoom_restore), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(slcontrols->vbox_zoom_buttons), slcontrols->button_zoom_restore, FALSE, FALSE, 0);
    g_signal_connect_swapped(slcontrols->button_zoom_restore, "clicked", G_CALLBACK(slopes_zoom_restore), data);

    gtk_box_pack_start(GTK_BOX(slcontrols->vbox_zoom_outer), slcontrols->vbox_zoom_buttons, TRUE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(slcontrols->hbox_graph), slcontrols->graph, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(slcontrols->hbox_graph), slcontrols->vbox_zoom_outer, FALSE, FALSE, 0);

    //final

    slcontrols->slopes_gui = gtk_hbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(slcontrols->slopes_gui), slcontrols->ctrls, FALSE, FALSE, 0);   /* expand, fill, padding */
    gtk_box_pack_start(GTK_BOX(slcontrols->slopes_gui), slcontrols->hbox_graph, TRUE, TRUE, 0);

    gtk_widget_show(slcontrols->slopes_gui);

    return slcontrols->slopes_gui;
}

void slopes_recalculate(Data *data)
{
    recalculate_Eit(data->controls->slopescontrols, &(data->args->slopesdata), &(data->args->instdata));
}

static void update_models_slopes(SlopesControls *slcontrols)
{
    GwyGraphCurveModel *cmodel;

    slcontrols->graph_model = gwy_graph_get_model(GWY_GRAPH(slcontrols->graph));

    cmodel = gwy_graph_model_get_curve_by_description(slcontrols->graph_model, "Loading");

    if (cmodel) {
        slcontrols->cmodelload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(slcontrols->graph_model, "Hold");

    if (cmodel) {
        slcontrols->cmodelhold = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(slcontrols->graph_model, "Unloading");

    if (cmodel) {
        slcontrols->cmodelunload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(slcontrols->graph_model, "Unloading fit");

    if (cmodel) {
        slcontrols->cmodelloadfit = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(slcontrols->graph_model, "Loading fit");

    if (cmodel) {
        slcontrols->cmodelunloadfit = cmodel;
    }
}

static void slopes_zoom_in(Data *data)
{
    SlopesControls *slcontrols;

    slcontrols = data->controls->slopescontrols;

    rs_zoom_in(GWY_GRAPH(slcontrols->graph), &(slcontrols->zoomstack), FALSE);   /* do not rescale y */
    update_models_slopes(slcontrols);
}

static void slopes_zoom_out(Data *data)
{
    SlopesControls *slcontrols;

    slcontrols = data->controls->slopescontrols;

    rs_zoom_out(GWY_GRAPH(slcontrols->graph), &(slcontrols->zoomstack));
    update_models_slopes(slcontrols);
}

void slopes_zoom_restore(Data *data)
{
    SlopesControls *slcontrols;

    slcontrols = data->controls->slopescontrols;

    rs_zoom_restore(GWY_GRAPH(slcontrols->graph), &(slcontrols->zoomstack));
    update_models_slopes(slcontrols);
}

#endif

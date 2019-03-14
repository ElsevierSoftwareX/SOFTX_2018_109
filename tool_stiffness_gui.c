#ifndef NOFORTRAN

#include "tool_stiffness_gui.h"
#include "tool_stiffness.h"
#include "tool_stiffness_unc_gui.h"

#include "controls.h"
#include "datatypes.h"
#include "file-utils.h"
#include "fit-utils.h"
#include "gui-utils.h"
#include "niget-common.h"

#include <stdlib.h>
#include <string.h>

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>
#include <app/file.h>

static void stiffness_remove_all_fit_labels_load(StiffnessControls *stcontrols);
static void stiffness_remove_all_fit_labels_unload(StiffnessControls *stcontrols);
static void stiffness_remove_all_fit_curves(StiffnessControls *stcontrols);

static void update_models_stiffness(StiffnessControls *stcontrols);
static void stiffness_zoom_in(Data *data);
static void stiffness_zoom_out(Data *data);

static void stiffness_show_log_load(Data *data);
static void stiffness_show_log_unload(Data *data);


void stiffness_remove_fit_results_labels(Data *data)
{
    gwy_selection_clear(data->controls->stiffnesscontrols->selection);
    stiffness_remove_all_fit_results_load(&(data->args->stiffnessdata));
    stiffness_remove_all_fit_results_unload(&(data->args->stiffnessdata));
    stiffness_remove_all_fit_labels_load(data->controls->stiffnesscontrols);
    stiffness_remove_all_fit_labels_unload(data->controls->stiffnesscontrols);
    stiffness_remove_all_fit_curves(data->controls->stiffnesscontrols);
}

static void stiffness_remove_all_fit_labels_unload(StiffnessControls *stcontrols)
{
    /* remove all labels that result from the straight line fit of the unloading curve */

    label_clear(stcontrols->fitinfounload);
    label_clear(stcontrols->ku);
    label_clear(stcontrols->qu);

    label_clear(stcontrols->unloadrange);

    label_clear(stcontrols->unloadSSres);
    label_clear(stcontrols->unloadR2);
    label_clear(stcontrols->unloadR2adj);
    label_clear(stcontrols->unloadchi2);
}

static void stiffness_remove_all_fit_labels_load(StiffnessControls *stcontrols)
{
    /* remove all labels that result from the straight line fit of the loading curve */

    label_clear(stcontrols->fitinfoload);
    label_clear(stcontrols->kl);
    label_clear(stcontrols->ql);

    label_clear(stcontrols->loadrange);

    label_clear(stcontrols->loadSSres);
    label_clear(stcontrols->loadR2);
    label_clear(stcontrols->loadR2adj);
    label_clear(stcontrols->loadchi2);
}

static void stiffness_remove_all_fit_curves(StiffnessControls *stcontrols)
{
    gwy_graph_model_remove_curve_by_description(stcontrols->graph_model, "Loading fit");
    gwy_graph_model_remove_curve_by_description(stcontrols->graph_model, "Unloading fit");
}

void stiffness_redraw(Data *data)
{
    GwyGraphCurveModel *cmodel;

    // update loading curve and Fmax and hmax labels
    if (data->args->fddata.hunload && data->args->fddata.hload && data->args->fddata.hhold) {

        gwy_graph_curve_model_set_data(data->controls->stiffnesscontrols->cmodelload,
                                       data->args->fddata.hload->data, data->args->fddata.Fload->data,
                                       data->args->fddata.hload->res);

        gwy_graph_curve_model_set_data(data->controls->stiffnesscontrols->cmodelhold,
                                       data->args->fddata.hhold->data, data->args->fddata.Fhold->data,
                                       data->args->fddata.hhold->res);

        gwy_graph_curve_model_set_data(data->controls->stiffnesscontrols->cmodelunload,
                                       data->args->fddata.hunload->data, data->args->fddata.Funload->data,
                                       data->args->fddata.hunload->res);

        label_set_gdouble_format(data->controls->stiffnesscontrols->hmax, data->args->fddata.hmax, DEPTH);
        label_set_gdouble_format(data->controls->stiffnesscontrols->Fmax, data->args->fddata.Fmax, FORCE);

        cmodel = gwy_graph_model_get_curve_by_description(data->controls->stiffnesscontrols->graph_model, "Loading fit");

        if (cmodel) {
            gwy_graph_curve_model_set_data(data->controls->stiffnesscontrols->cmodelloadfit,
                                           data->args->stiffnessdata.xloadfit->data, data->args->stiffnessdata.yloadfit->data,  data->args->stiffnessdata.xloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(data->controls->stiffnesscontrols->cmodelloadfit);
#endif
        }

        cmodel = gwy_graph_model_get_curve_by_description(data->controls->stiffnesscontrols->graph_model, "Unloading fit");

        if (cmodel) {
            gwy_graph_curve_model_set_data(data->controls->stiffnesscontrols->cmodelunloadfit,
                                           data->args->stiffnessdata.xunloadfit->data, data->args->stiffnessdata.yunloadfit->data,  data->args->stiffnessdata.xunloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(data->controls->stiffnesscontrols->cmodelunloadfit);
#endif
        }
    }
}

static void stiffness_range_from_changed(Data *data)
{
    StiffnessControls *fcontrols;
    Stiffnessdata *fdata;

    fcontrols = data->controls->stiffnesscontrols;
    fdata = &(data->args->stiffnessdata);

    range_from_changed(fcontrols->from, fcontrols->to, fcontrols->selection, &(fdata->loadfrom), &(fdata->unloadfrom));
}

static void stiffness_range_to_changed(Data *data)
{
    StiffnessControls *fcontrols;
    Stiffnessdata *fdata;

    fcontrols = data->controls->stiffnesscontrols;
    fdata = &(data->args->stiffnessdata);

    range_to_changed(fcontrols->from, fcontrols->to, fcontrols->selection, &(fdata->loadto), &(fdata->unloadto));
}

static void stiffness_graph_selected(GwySelection *selection, gint i, Data *data)
{
    gdouble range[2];

    g_return_if_fail(i <= 0);

    //take data from selection and put them in the "from" and "to" values for both hp and S (at this moment, we don't know yet which fit will be chosen), create labels only once
    //StiffnessControls has entries from and to, which are filled here
    // and labels hp,unloadrange, where the values of the ranges that were actually used in the corresponding fits  are shown
    //Stiffnessdata has loadfrom,loadto, unloadfrom,unloadto which are the data from the entries, i.e. from the selection
    //Stiffnessdata; also has loadrange, unloadrange where the ranges that were actually used in the fit

    gwy_selection_get_object(selection, 0, range);

    data->args->stiffnessdata.unloadfrom = MIN(range[0], range[1]);
    data->args->stiffnessdata.unloadto = MAX(range[0], range[1]);

    data->args->stiffnessdata.loadfrom = MIN(range[0], range[1]);
    data->args->stiffnessdata.loadto = MAX(range[0], range[1]);

    entry_set_gdouble_format(data->controls->stiffnesscontrols->from, data->args->stiffnessdata.unloadfrom, DEPTH);
    entry_set_gdouble_format(data->controls->stiffnesscontrols->to, data->args->stiffnessdata.unloadto, DEPTH);

}

void stiffness_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (stiffness_has_results(&(args->stiffnessdata))) {
	buffer = stiffness_export_data(&(args->stiffnessdata), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

void stiffness_fit_and_update_load(Data *data)
{
    gchar str[300];
    gint i;
    gint istart, iend, ndata;
    gboolean fit = TRUE;
    GwyDataLine *x, *y;
    gint fitinfo;

    Stiffnessdata *stdata;
    StiffnessControls *stcontrols;
    FDdata *fddata;
    Instdata *instdata;
    GdkColor color;
    GwyGraphCurveModel *cmodel;

    gdk_color_parse("red", &color);

    //create pointers for easier reading
    stdata = &(data->args->stiffnessdata);
    fddata = &(data->args->fddata);
    instdata = &(data->args->instdata);
    stcontrols = data->controls->stiffnesscontrols;

    //if there was a previous fit, invalidate previous fitted curve
    if (stdata->has_fit_load) {
        g_object_unref(stdata->xloadfit);
        g_object_unref(stdata->yloadfit);
    }

    if (verbose) {
        g_print("from %g to %g \n", stdata->loadfrom, stdata->loadto);
    }

    //check range
    if (stdata->loadfrom == stdata->loadto) {
        fit = FALSE;
    }

    //find appropriate data, which should be fitted

    range_to_indices(stdata->loadfrom, stdata->loadto, fddata->hload, FALSE, &istart, &iend, &ndata);

    // if there are not enough data in the selected range, no fitting can be performed
    ndata = MIN(ndata, fddata->hload->res - istart);

    if (ndata < 2) {
        fit = FALSE;
    }

    stdata->nfitloaddata = ndata;

    if (fit) {
        //get data to be fitted
        x = gwy_data_line_part_extract(fddata->hload, istart, ndata);
        y = gwy_data_line_part_extract(fddata->Fload, istart, ndata);

        //fit the data and evaluate
        fitinfo = stiffness_fit_load(x, y, stdata, fddata, instdata, NULL);

        if (stdata->has_fit_unload) {
            stiffness_combine_fit_results(stdata, fddata, instdata);
        }

        if (fitinfo >= 1e6) {
            stdata->infologload = fitinfo / 1000000;
            fitinfo = fitinfo % 1000000;
        }

        if (fitinfo >= 5) {
            gtk_label_set_text(GTK_LABEL(stcontrols->fitinfoload), "Warning Load fit: Questionable results \n or fatal errors in fitting!");
            gtk_widget_modify_fg(stcontrols->fitinfoload, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(stcontrols->fitinfoload), TRUE);
        }
        else if (fitinfo == 4) {
            gtk_label_set_text(GTK_LABEL(stcontrols->fitinfoload), "Warning Load fit: Iteration limit reached!");
            gtk_widget_modify_fg(stcontrols->fitinfoload, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(stcontrols->fitinfoload), TRUE);
        }
        else {
            label_clear(stcontrols->fitinfoload);
        }

        //update range labels
        entry_set_gdouble_format(stcontrols->from, stdata->loadfrom, DEPTH);
        entry_set_gdouble_format(stcontrols->to, stdata->loadto, DEPTH);

        //update labels
        label_set_gdouble_format(stcontrols->kl, stdata->kl, SLOPE);
        label_set_gdouble_format(stcontrols->ql, stdata->ql, DEPTH);

        label_set_gdouble_format(stcontrols->loadSSres, stdata->loadSSres, "%.4g");
        label_set_gdouble_format(stcontrols->loadR2, stdata->loadR2, "%.4g");
        label_set_gdouble_format(stcontrols->loadR2adj, stdata->loadR2adj, "%.4g");
        label_set_gdouble_format(stcontrols->loadchi2, stdata->loadchi2, "%.4g");

        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", stdata->loadfrom, stdata->loadto);
        gtk_label_set_text(GTK_LABEL(stcontrols->loadrange), str);

        //calculate curve of fit
        stdata->xloadfit = gwy_data_line_duplicate(x);
        stdata->yloadfit = gwy_data_line_new_alike(y, FALSE);

        for (i = 0; i < stdata->xloadfit->res; i++) {
            stdata->yloadfit->data[i] =  stdata->kl * stdata-> xloadfit->data[i] + stdata->ql;
        }

        // clean up
        g_object_unref(x);
        g_object_unref(y);

        //update fitted curve, if necessary create it
        cmodel = gwy_graph_model_get_curve_by_description(stcontrols->graph_model, "Loading fit");

        if (cmodel) {
            gwy_graph_curve_model_set_data(stcontrols->cmodelloadfit, stdata->xloadfit->data, stdata->yloadfit->data, stdata->xloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(stcontrols->cmodelloadfit);
#endif
        }
        else {
            gwy_graph_curve_model_set_data(stcontrols->cmodelloadfit, stdata->xloadfit->data, stdata->yloadfit->data, stdata->xloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(stcontrols->cmodelloadfit);
#endif
            gwy_graph_model_add_curve(stcontrols->graph_model, stcontrols->cmodelloadfit);
            g_object_set(stcontrols->cmodelloadfit,
                         "mode", GWY_GRAPH_CURVE_LINE,
                         "description", "Loading fit",
                         "color", CURVE_COLOR_FIT,
                         "line-style", GDK_LINE_SOLID,
                         "line-width", 3,
                         NULL);
        }
    }
    else {
        stiffness_remove_all_fit_labels_load(stcontrols);
        stiffness_remove_all_fit_results_load(stdata);
        stdata->has_fit_load = FALSE;
        gwy_graph_model_remove_curve_by_description(stcontrols->graph_model, "Loading fit");
    }

    // log of fitting procedure is available
    if (stdata->has_fit_load) {
        gtk_widget_set_sensitive(stcontrols->button_log_load, TRUE);
    }
    else {
        gtk_widget_set_sensitive(stcontrols->button_log_load, FALSE);
    }

    if (stdata->has_fit_load &&  stdata->has_fit_unload) {
        if (data->args->stiffnessunc.klc != NULL || data->args->stiffnessunc.qlc != NULL || data->args->stiffnessunc.kuc != NULL || data->args->stiffnessunc.quc != NULL) {
            stiffness_unc_close_window(data);
        }

        gtk_widget_set_sensitive(stcontrols->button_unc, TRUE);
    }
    else {
        gtk_widget_set_sensitive(stcontrols->button_unc, FALSE);
    }
}

void stiffness_fit_and_update_unload(Data *data)
{
    gchar str[300];
    gint  i;
    gint istart, iend, ndata;
    gint fitinfo;
    gboolean fit = TRUE;
    GwyDataLine *x, *y;

    Stiffnessdata *stdata;
    StiffnessControls *stcontrols;
    FDdata *fddata;
    Instdata *instdata;
    GdkColor color;
    GwyGraphCurveModel *cmodel;

    gdk_color_parse("red", &color);

    //create pointers for easier reading
    stdata = &(data->args->stiffnessdata);
    fddata = &(data->args->fddata);
    instdata = &(data->args->instdata);
    stcontrols = data->controls->stiffnesscontrols;

    //if there was a previous fit, free the old data and create new data
    if (stdata->has_fit_unload) {
        g_object_unref(stdata->xunloadfit);
        g_object_unref(stdata->yunloadfit);
    }

    if (verbose) {
        g_print("from %g to %g \n", stdata->unloadfrom, stdata->unloadto);
    }

    // check range
    if (stdata->unloadfrom == stdata->unloadto) {
        fit = FALSE;
    }

    //find appropriate data, which should be fitted
    range_to_indices(stdata->unloadfrom, stdata->unloadto, fddata->hunload, TRUE, &istart, &iend, &ndata);

    // if there are not enough data in the selected range,no fitting can be performed
    ndata = MIN(ndata, fddata->hunload->res - istart);

    if (ndata < 2) {
        fit = FALSE;
    }

    stdata->nfitunloaddata = ndata;

    if (fit) {
        //get data to be fitted
        x = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
        y = gwy_data_line_part_extract(fddata->Funload, istart, ndata);

        //fit the data and evaluate
        fitinfo = stiffness_fit_unload(x, y, stdata, fddata, instdata, NULL);

        if (stdata->has_fit_load) {
            stiffness_combine_fit_results(stdata, fddata, instdata);
        }

        if (fitinfo >= 1e6) {
            stdata->infologunload = fitinfo / 1000000;
            fitinfo = fitinfo % 1000000;
        }

        if (fitinfo >= 5) {
            gtk_label_set_text(GTK_LABEL(stcontrols->fitinfounload), "Warning Unload fit: Questionable results \n or fatal errors in fitting!");
            gtk_widget_modify_fg(stcontrols->fitinfounload, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(stcontrols->fitinfounload), TRUE);
        }
        else if (fitinfo == 4) {
            gtk_label_set_text(GTK_LABEL(stcontrols->fitinfounload), "Warning Unload fit: Iteration limit reached!");
            gtk_widget_modify_fg(stcontrols->fitinfounload, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(stcontrols->fitinfounload), TRUE);
        }
        else {
            label_clear(stcontrols->fitinfounload);
        }

        //update range labels
        entry_set_gdouble_format(stcontrols->from, stdata->unloadfrom, DEPTH);
        entry_set_gdouble_format(stcontrols->to, stdata->unloadto, DEPTH);

        //update labels
        label_set_gdouble_format(stcontrols->ku, stdata->ku, SLOPE);
        label_set_gdouble_format(stcontrols->qu, stdata->qu, DEPTH);

        label_set_gdouble_format(stcontrols->unloadSSres, stdata->unloadSSres, "%.4g");
        label_set_gdouble_format(stcontrols->unloadR2, stdata->unloadR2, "%.4g");
        label_set_gdouble_format(stcontrols->unloadR2adj, stdata->unloadR2adj, "%.4g");
        label_set_gdouble_format(stcontrols->unloadchi2, stdata->unloadchi2, "%.4g");

        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", stdata->unloadfrom, stdata->unloadto);
        gtk_label_set_text(GTK_LABEL(stcontrols->unloadrange), str);

        //calculate curve of fit
        stdata->xunloadfit = gwy_data_line_duplicate(x);
        stdata->yunloadfit = gwy_data_line_new_alike(y, FALSE);

        for (i = 0; i < stdata->xunloadfit->res; i++) {
            stdata->yunloadfit->data[i] =  stdata->ku * stdata->xunloadfit->data[i] + stdata->qu;
        }

        // clean up
        g_object_unref(x);
        g_object_unref(y);

        //update fitted curve, if necessary create it
        cmodel = gwy_graph_model_get_curve_by_description(stcontrols->graph_model, "Unloading fit");

        if (cmodel) {
            gwy_graph_curve_model_set_data(stcontrols->cmodelunloadfit, stdata->xunloadfit->data, stdata->yunloadfit->data, stdata->xunloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(stcontrols->cmodelunloadfit);
#endif
        }
        else {
            gwy_graph_curve_model_set_data(stcontrols->cmodelunloadfit, stdata->xunloadfit->data, stdata->yunloadfit->data, stdata->xunloadfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(stcontrols->cmodelunloadfit);
#endif
            gwy_graph_model_add_curve(stcontrols->graph_model, stcontrols->cmodelunloadfit);
            g_object_set(stcontrols->cmodelunloadfit,
                         "mode", GWY_GRAPH_CURVE_LINE,
                         "description", "Unloading fit",
                         "color", CURVE_COLOR_FIT2,
                         "line-style", GDK_LINE_SOLID,
                         "line-width", 3,
                         NULL);
        }
    }
    else {
        stiffness_remove_all_fit_labels_unload(stcontrols);
        stiffness_remove_all_fit_results_unload(stdata);
        stdata->has_fit_unload = FALSE;
        gwy_graph_model_remove_curve_by_description(stcontrols->graph_model, "Unloading fit");
    }

    // log of fitting procedure is available
    if (stdata->has_fit_unload) {
        gtk_widget_set_sensitive(stcontrols->button_log_unload, TRUE);
    }
    else {
        gtk_widget_set_sensitive(stcontrols->button_log_unload, FALSE);
    }

    if (stdata->has_fit_load &&  stdata->has_fit_unload) {
        if (data->args->stiffnessunc.kuc != NULL || data->args->stiffnessunc.quc != NULL || data->args->stiffnessunc.klc != NULL || data->args->stiffnessunc.qlc) {
            stiffness_unc_close_window(data);
        }

        gtk_widget_set_sensitive(stcontrols->button_unc, TRUE);
    }
    else {
        gtk_widget_set_sensitive(stcontrols->button_unc, FALSE);
    }
}

static void stiffness_show_log_load(Data *data)
{
    gint infolog;

    infolog = data->args->stiffnessdata.infologload;

    /* check if there are logs */
    switch (infolog) {
    case 0:
	show_log(data->controls->stiffnesscontrols->load_log_error_window, data->args->stiffnessdata.logfnmload, ".err", "Fit error log loading");
	show_log(data->controls->stiffnesscontrols->load_log_report_window, data->args->stiffnessdata.logfnmload, ".rpt", "Fit report log loading");
	break;
	
    case 1:
	show_log(data->controls->stiffnesscontrols->load_log_error_window,  data->args->stiffnessdata.logfnmload, ".err", "Fit error log loading");
	show_warning("Fit report log could not be created.");
	break;

    case 2:
	show_log(data->controls->stiffnesscontrols->load_log_report_window, data->args->stiffnessdata.logfnmload, ".rpt", "Fit report log loading");
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

static void stiffness_show_log_unload(Data *data)
{
    gint infolog;

    infolog = data->args->stiffnessdata.infologunload;

    /* check if there are logs */
    switch (infolog) {
    case 0:
	show_log(data->controls->stiffnesscontrols->unload_log_error_window, data->args->stiffnessdata.logfnmunload, ".err", "Fit error log unloading");
	show_log(data->controls->stiffnesscontrols->unload_log_report_window, data->args->stiffnessdata.logfnmunload, ".rpt", "Fit report log unloading");
	break;
	
    case 1:
	show_log(data->controls->stiffnesscontrols->unload_log_error_window,  data->args->stiffnessdata.logfnmunload, ".err", "Fit error log unloading");
	show_warning("Fit report log could not be created.");
	break;

    case 2:
	show_log(data->controls->stiffnesscontrols->unload_log_report_window, data->args->stiffnessdata.logfnmunload, ".rpt", "Fit report log unloading");
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

GtkWidget *tool_stiffness_create(Data *data)
{
    GwyGraphArea *area;
    Args *args;
    StiffnessControls *stcontrols;
    gint i;

    args = data->args;
    stcontrols = data->controls->stiffnesscontrols;


    stcontrols->graph_model = gwy_graph_model_new();
    g_object_set(stcontrols->graph_model, "title", "Unload", NULL);
    gwy_graph_model_set_axis_label(stcontrols->graph_model, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(stcontrols->graph_model, GTK_POS_LEFT, "F / mN");

    stcontrols->zoomstack = NULL;
    stcontrols->zoomstack = g_slist_prepend(stcontrols->zoomstack, stcontrols->graph_model);

    stcontrols->hbox_graph = gtk_hbox_new(FALSE, CTRLS_SPACING);

    stcontrols->graph = gwy_graph_new(stcontrols->graph_model);

    gwy_graph_enable_user_input(GWY_GRAPH(stcontrols->graph), FALSE);
    gwy_graph_set_status(GWY_GRAPH(stcontrols->graph), GWY_GRAPH_STATUS_XSEL);

    area = GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(stcontrols->graph)));
    stcontrols->selection = gwy_graph_area_get_selection(area, GWY_GRAPH_STATUS_XSEL);
    gwy_selection_set_max_objects(stcontrols->selection, 1);
    g_signal_connect(stcontrols->selection, "changed", G_CALLBACK(stiffness_graph_selected), data);

    stcontrols->cmodelload = gwy_graph_curve_model_new();
    stcontrols->cmodelhold = gwy_graph_curve_model_new();
    stcontrols->cmodelunload = gwy_graph_curve_model_new();
    stcontrols->cmodelloadfit = gwy_graph_curve_model_new();
    stcontrols->cmodelunloadfit = gwy_graph_curve_model_new();

    g_object_set(stcontrols->cmodelload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Loading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_LOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(stcontrols->cmodelhold,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Hold",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_HOLD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(stcontrols->cmodelunload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Unloading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_UNLOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    gwy_graph_model_add_curve(stcontrols->graph_model, stcontrols->cmodelload);
    gwy_graph_model_add_curve(stcontrols->graph_model, stcontrols->cmodelhold);
    gwy_graph_model_add_curve(stcontrols->graph_model, stcontrols->cmodelunload);


    // Table of input data and results

    stcontrols->ctrls = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* Info section */

    stcontrols->frame_info = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(stcontrols->frame_info), GTK_WIDGET(gwy_label_new_header("Info")));
    stcontrols->table_info = gtk_table_new(2, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(stcontrols->table_info), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(stcontrols->table_info), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(stcontrols->table_info), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(stcontrols->frame_info), stcontrols->table_info);

    i = 0;

    /* h_max */

    gtk_table_attach(GTK_TABLE(stcontrols->table_info), label_new_with_markup_left("h<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    stcontrols->hmax = gtk_label_new(NULL);
    label_set_gdouble_format(stcontrols->hmax, args->fddata.hmax, DEPTH);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->hmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_info), stcontrols->hmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(stcontrols->table_info), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /* F_max */

    gtk_table_attach(GTK_TABLE(stcontrols->table_info), label_new_with_markup_left("F<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    stcontrols->Fmax = gtk_label_new(NULL);
    label_set_gdouble_format(stcontrols->Fmax, args->fddata.Fmax, FORCE);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_info), stcontrols->Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(stcontrols->table_info), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    /* Parameters section */

    stcontrols->frame_parameters = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(stcontrols->frame_parameters), GTK_WIDGET(gwy_label_new_header("Parameters")));
    stcontrols->table_parameters = gtk_table_new(5, 4, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(stcontrols->table_parameters), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(stcontrols->table_parameters), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(stcontrols->table_parameters), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(stcontrols->frame_parameters), stcontrols->table_parameters);
    i = 0;

    gtk_table_attach(GTK_TABLE(stcontrols->table_parameters), label_new_with_markup_left("Fit range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->from = gtk_entry_new();
    g_object_set_data(G_OBJECT(stcontrols->from), "id", (gpointer)"from");
    gtk_entry_set_width_chars(GTK_ENTRY(stcontrols->from), 8);
    gtk_table_attach(GTK_TABLE(stcontrols->table_parameters), stcontrols->from, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(stcontrols->from, TRUE);
    g_signal_connect_swapped(stcontrols->from, "activate", G_CALLBACK(stiffness_range_from_changed), data);

    stcontrols->to = gtk_entry_new();
    g_object_set_data(G_OBJECT(stcontrols->to), "id", (gpointer)"to");
    gtk_entry_set_width_chars(GTK_ENTRY(stcontrols->to), 8);
    gtk_table_attach(GTK_TABLE(stcontrols->table_parameters), stcontrols->to, 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(stcontrols->to, TRUE);
    g_signal_connect_swapped(stcontrols->to, "activate", G_CALLBACK(stiffness_range_to_changed), data);

    gtk_table_attach(GTK_TABLE(stcontrols->table_parameters), label_new_with_markup_left("nm"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;


    /* Fit buttons */

    stcontrols->hbox_fitbuttons = gtk_hbox_new(TRUE, 0);

    stcontrols->button_fit_load = gtk_button_new_with_mnemonic(gwy_sgettext("_Fit loading"));
    gtk_button_set_image(GTK_BUTTON(stcontrols->button_fit_load), gtk_image_new_from_stock(GTK_STOCK_EXECUTE, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(stcontrols->hbox_fitbuttons), stcontrols->button_fit_load, TRUE, TRUE, 0);
    g_signal_connect_swapped(stcontrols->button_fit_load, "clicked", G_CALLBACK(stiffness_fit_and_update_load), data);

    stcontrols->button_fit_unload = gtk_button_new_with_mnemonic(gwy_sgettext("verb|F_it unloading"));
    gtk_button_set_image(GTK_BUTTON(stcontrols->button_fit_unload), gtk_image_new_from_stock(GTK_STOCK_EXECUTE, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(stcontrols->hbox_fitbuttons), stcontrols->button_fit_unload, TRUE, TRUE, 0);
    g_signal_connect_swapped(stcontrols->button_fit_unload, "clicked", G_CALLBACK(stiffness_fit_and_update_unload), data);


    /* Results section */

    stcontrols->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(stcontrols->frame_results), GTK_WIDGET(gwy_label_new_header("Results")));

    stcontrols->scrolledwindow_results = gtk_scrolled_window_new(NULL, NULL);
    /* TODO: resolve table column widths so that horizontal scroll policy can become NEVER */
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(stcontrols->scrolledwindow_results), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(stcontrols->frame_results), stcontrols->scrolledwindow_results);

    stcontrols->table_results = gtk_table_new(20, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(stcontrols->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(stcontrols->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(stcontrols->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(stcontrols->scrolledwindow_results), stcontrols->table_results);
    i = 0;

    stcontrols->fitinfoload = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->fitinfoload, 0, 3, i, i + 2, GTK_FILL, GTK_FILL, 0, 0);
    i += 2;


    /*    kl    */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("k<sub>load</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->kl = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->kl), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->kl, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("N/nm<sup>-1</sup>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /*    ql    */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("q<sub>load</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->ql = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->ql), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->ql, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("N/nm<sup>-1</sup>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;


    /* SSres */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("SS<sub>res</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->loadSSres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->loadSSres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->loadSSres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* R2 */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("R<sup>2</sup>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->loadR2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->loadR2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->loadR2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* R2adj */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("R<sup>2</sup><sub>adj</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->loadR2adj = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->loadR2adj), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->loadR2adj, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* chi2 */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("χ<sup>2</sup> fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->loadchi2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->loadchi2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->loadchi2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;


    /* unload */

    stcontrols->fitinfounload = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->fitinfounload, 0, 3, i, i + 2, GTK_FILL, GTK_FILL, 0, 0);
    i += 2;

    /*    ku    */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("k<sub>unload</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->ku = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->ku), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->ku, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("N/nm<sup>-1</sup>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /*    qu    */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("q<sub>unload</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->qu = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->qu), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->qu, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("N/nm<sup>-1</sup>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;


    /* SSres */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("SS<sub>res</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->unloadSSres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->unloadSSres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->unloadSSres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* R2 */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("R<sup>2</sup>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->unloadR2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->unloadR2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->unloadR2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* R2adj */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("R<sup>2</sup><sub>adj</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->unloadR2adj = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->unloadR2adj), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->unloadR2adj, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* chi2 */
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("χ<sup>2</sup> fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->unloadchi2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->unloadchi2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->unloadchi2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;


    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("Fitted loading range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->loadrange = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->loadrange), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->loadrange, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;


    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("Fitted unloading range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    stcontrols->unloadrange = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(stcontrols->unloadrange), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), stcontrols->unloadrange, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(stcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;


    /* Log buttons */
    stcontrols->hbox_logbuttons = gtk_hbox_new(TRUE, 0);

    stcontrols->button_log_load = gtk_button_new_with_mnemonic("Show _log load");
    gtk_button_set_image(GTK_BUTTON(stcontrols->button_log_load), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(stcontrols->button_log_load, "clicked", G_CALLBACK(stiffness_show_log_load), data);
    gtk_widget_set_sensitive(stcontrols->button_log_load, FALSE);
    stcontrols->load_log_error_window = NULL;
    stcontrols->load_log_report_window = NULL;
    gtk_box_pack_start(GTK_BOX(stcontrols->hbox_logbuttons), stcontrols->button_log_load, TRUE, TRUE, 0);

    stcontrols->button_log_unload = gtk_button_new_with_mnemonic("Show lo_g unload");
    gtk_button_set_image(GTK_BUTTON(stcontrols->button_log_unload), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(stcontrols->button_log_unload, "clicked", G_CALLBACK(stiffness_show_log_unload), data);
    gtk_widget_set_sensitive(stcontrols->button_log_unload, FALSE);
    stcontrols->unload_log_error_window = NULL;
    stcontrols->unload_log_report_window = NULL;
    gtk_box_pack_start(GTK_BOX(stcontrols->hbox_logbuttons), stcontrols->button_log_unload, TRUE, TRUE, 0);


    // Uncertainties Button
    stcontrols->button_unc = gtk_button_new_with_mnemonic("_Uncertainties");
    gtk_button_set_image(GTK_BUTTON(stcontrols->button_unc), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(stcontrols->button_unc, "clicked", G_CALLBACK(stiffness_uncertainty), data);
    gtk_widget_set_sensitive(stcontrols->button_unc, FALSE);


    // Save Button

    stcontrols->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(stcontrols->button_save, "clicked", G_CALLBACK(stiffness_save_data), args);

    gtk_box_pack_start(GTK_BOX(stcontrols->ctrls), stcontrols->frame_info, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(stcontrols->ctrls), stcontrols->frame_parameters, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(stcontrols->ctrls), stcontrols->hbox_fitbuttons, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(stcontrols->ctrls), stcontrols->frame_results, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(stcontrols->ctrls), stcontrols->hbox_logbuttons, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(stcontrols->ctrls), stcontrols->button_unc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(stcontrols->ctrls), stcontrols->button_save, FALSE, FALSE, 0);

    /* Zooming */

    stcontrols->vbox_zoom_buttons = gtk_vbox_new(TRUE, 0);
    stcontrols->vbox_zoom_outer = gtk_vbox_new(FALSE, 0);

    stcontrols->button_zoom_in = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(stcontrols->button_zoom_in), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(stcontrols->vbox_zoom_buttons), stcontrols->button_zoom_in, FALSE, FALSE, 0);
    g_signal_connect_swapped(stcontrols->button_zoom_in, "clicked", G_CALLBACK(stiffness_zoom_in), data);

    stcontrols->button_zoom_out = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(stcontrols->button_zoom_out), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(stcontrols->vbox_zoom_buttons), stcontrols->button_zoom_out, FALSE, FALSE, 0);
    g_signal_connect_swapped(stcontrols->button_zoom_out, "clicked", G_CALLBACK(stiffness_zoom_out), data);

    stcontrols->button_zoom_restore = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(stcontrols->button_zoom_restore), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(stcontrols->vbox_zoom_buttons), stcontrols->button_zoom_restore, FALSE, FALSE, 0);
    g_signal_connect_swapped(stcontrols->button_zoom_restore, "clicked", G_CALLBACK(stiffness_zoom_restore), data);

    gtk_box_pack_start(GTK_BOX(stcontrols->vbox_zoom_outer), stcontrols->vbox_zoom_buttons, TRUE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(stcontrols->hbox_graph), stcontrols->graph, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(stcontrols->hbox_graph), stcontrols->vbox_zoom_outer, FALSE, FALSE, 0);

    //final

    stcontrols->stiffness_gui = gtk_hbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(stcontrols->stiffness_gui), stcontrols->ctrls, FALSE, FALSE, 0);   /* expand, fill, padding */
    gtk_box_pack_start(GTK_BOX(stcontrols->stiffness_gui), stcontrols->hbox_graph, TRUE, TRUE, 0);

    gtk_widget_show(stcontrols->stiffness_gui);

    return stcontrols->stiffness_gui;
}

static void update_models_stiffness(StiffnessControls *stcontrols)
{
    GwyGraphCurveModel *cmodel;

    stcontrols->graph_model = gwy_graph_get_model(GWY_GRAPH(stcontrols->graph));

    cmodel = gwy_graph_model_get_curve_by_description(stcontrols->graph_model, "Loading");

    if (cmodel) {
        stcontrols->cmodelload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(stcontrols->graph_model, "Hold");

    if (cmodel) {
        stcontrols->cmodelhold = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(stcontrols->graph_model, "Unloading");

    if (cmodel) {
        stcontrols->cmodelunload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(stcontrols->graph_model, "Unloading fit");

    if (cmodel) {
        stcontrols->cmodelloadfit = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(stcontrols->graph_model, "Loading fit");

    if (cmodel) {
        stcontrols->cmodelunloadfit = cmodel;
    }
}

static void stiffness_zoom_in(Data *data)
{
    StiffnessControls *stcontrols;

    stcontrols = data->controls->stiffnesscontrols;

    rs_zoom_in(GWY_GRAPH(stcontrols->graph), &(stcontrols->zoomstack), FALSE);   /* do not rescale y */
    update_models_stiffness(stcontrols);
}

static void stiffness_zoom_out(Data *data)
{
    StiffnessControls *stcontrols;

    stcontrols = data->controls->stiffnesscontrols;

    rs_zoom_out(GWY_GRAPH(stcontrols->graph), &(stcontrols->zoomstack));
    update_models_stiffness(stcontrols);
}

void stiffness_zoom_restore(Data *data)
{
    StiffnessControls *stcontrols;

    stcontrols = data->controls->stiffnesscontrols;

    rs_zoom_restore(GWY_GRAPH(stcontrols->graph), &(stcontrols->zoomstack));
    update_models_stiffness(stcontrols);
}

#endif

#ifndef NOFORTRAN

#include "tool_op_odr_gui.h"
#include "tool_op_odr.h"
#include "tool_op_odr_unc_gui.h"

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

#include <app/file.h>
#include <libgwydgets/gwydgets.h>
#include <libprocess/gwyprocess.h>

static void op_odr_show_log(Data *data);

static void radial_toggled(GtkToggleButton *button, Data *data);
static void op_odr_remove_all_fit_labels(OPODRControls *opodrcontrols);
static void op_odr_remove_all_fit_curves(OPODRControls *opodrcontrols);

static void op_odr_range_from_changed(Data *data);
static void op_odr_range_to_changed(Data *data);
static gboolean op_odr_range_from_pct_Fmax_changed(Data *data);
static gboolean op_odr_range_to_pct_Fmax_changed(Data *data);

static void update_models_op_odr(OPODRControls *opodrcontrols);
static void op_odr_zoom_in(Data *data);
static void op_odr_zoom_out(Data *data);

void op_odr_fit_and_update(Data *data);


void op_odr_remove_fit_results_labels(Data *data)
{
    gwy_selection_clear(data->controls->opodrcontrols->selection);
    op_odr_remove_all_fit_results(&(data->args->opodrdata));
    op_odr_remove_all_fit_labels(data->controls->opodrcontrols);
    op_odr_remove_all_fit_curves(data->controls->opodrcontrols);

    if (data->args->opodrunc.Ec != NULL || data->args->opodrunc.Hc != NULL) {
        op_odr_unc_close_window(data);
    }
}

static void op_odr_remove_all_fit_labels(OPODRControls *opodrcontrols)
{
    /* remove all labels that result from a fit */
    entry_clear(opodrcontrols->from);
    entry_clear(opodrcontrols->to);
    entry_clear(opodrcontrols->from_pct_Fmax);
    entry_clear(opodrcontrols->to_pct_Fmax);

    /* R currently not available */
    /* label_clear(opodrcontrols->R); */
    label_clear(opodrcontrols->SSres);
    label_clear(opodrcontrols->R2);
    label_clear(opodrcontrols->R2adj);
    label_clear(opodrcontrols->chi2);
    label_clear(opodrcontrols->fitinfo);
    label_clear(opodrcontrols->S);
    label_clear(opodrcontrols->hp);
    label_clear(opodrcontrols->m);
    label_clear(opodrcontrols->alpha);
    label_clear(opodrcontrols->eps);
    label_clear(opodrcontrols->hc);
    label_clear(opodrcontrols->Aphc);
    label_clear(opodrcontrols->extrapol);

    label_clear(opodrcontrols->Hit);
    label_clear(opodrcontrols->Er);
    label_clear(opodrcontrols->Eit);
    label_clear(opodrcontrols->range);
    label_clear(opodrcontrols->range_pct_Fmax);
}

static void op_odr_remove_all_fit_curves(OPODRControls *opodrcontrols)
{
    gwy_graph_model_remove_curve_by_description(opodrcontrols->graph_model, "Fit");
}

static void recalculate_Eit(OPODRControls *opodrcontrols, OPODRdata *opodrdata, const Instdata *instdata)
{
    gdouble Er;

    Er = opodrdata->Er * 1e9;

    // recalc for changed  nu, nui or Ei
    if (Er == 0) {
        return;
    }

    opodrdata->Eit = calc_Eit(Er, instdata->nu, instdata->Ei, instdata->nui) * 1e-9;
    label_set_gdouble_format(opodrcontrols->Eit, opodrdata->Eit, MODUL);
}

void op_odr_redraw(Data *data)
{
    // update loading curve and Fmax and hmax labels
    if (data->args->fddata.hunload) {
        gwy_graph_curve_model_set_data(data->controls->opodrcontrols->cmodelunload,
                                       data->args->fddata.hunload->data, data->args->fddata.Funload->data,
                                       data->args->fddata.hunload->res);

        label_set_gdouble_format(data->controls->opodrcontrols->hmax, data->args->fddata.hmax, DEPTH);
        label_set_gdouble_format(data->controls->opodrcontrols->Fmax, data->args->fddata.Fmax, FORCE);

        if (gwy_graph_model_get_n_curves(data->controls->opodrcontrols->graph_model) == 2) {
            gwy_graph_curve_model_set_data(data->controls->opodrcontrols->cmodelSfit,
                                           data->args->opodrdata.xfit->data, data->args->opodrdata.yfit->data,  data->args->opodrdata.xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(data->controls->opodrcontrols->cmodelSfit);
#endif
        }
    }
}

static void opodr_recalculate_beta(OPODRControls *opodrcontrols, OPODRdata *opodrdata, const Instdata *instdata)
{
    if (!opodrdata->has_fit) {
        return;
    }

    opodrdata->Er = calc_Er_OP(opodrdata->Aphc, opodrdata->S, opodrdata->beta) * 1e15;
    opodrdata->Eit = calc_Eit(opodrdata->Er, instdata->nu, instdata->Ei, instdata->nui);
    opodrdata->Eit *= 1e-9;
    opodrdata->Er *= 1e-9;

    label_set_gdouble_format(opodrcontrols->Eit, opodrdata->Eit, MODUL);
    label_set_gdouble_format(opodrcontrols->Er, opodrdata->Er, MODUL);
}

static void opodr_beta_changed(Data *data)
{
    OPODRControls *fcontrols;

    fcontrols = data->controls->opodrcontrols;

    if (beta_changed(fcontrols->beta, &(data->args->opodrdata.beta))) {
        opodr_recalculate_beta(data->controls->opodrcontrols, &(data->args->opodrdata), &(data->args->instdata));
    }
}

static void opodr_reset_defaults(Data *data)
{
    OPODRControls *opodrcontrols;
    OPODRdata *opodrdata;

    opodrcontrols = data->controls->opodrcontrols;
    opodrdata = &(data->args->opodrdata);
    
    opodrdata->beta = BETA_DEFAULT_OPODR;
    opodrdata->from_pct_Fmax = RANGE_FROM_DEFAULT_OPODR;
    opodrdata->to_pct_Fmax = RANGE_TO_DEFAULT_OPODR;
    entry_set_gdouble_format(opodrcontrols->beta, opodrdata->beta, "%.03f");
    entry_set_gdouble_format(opodrcontrols->from_pct_Fmax, opodrdata->from_pct_Fmax, "%.02f");
    entry_set_gdouble_format(opodrcontrols->to_pct_Fmax, opodrdata->to_pct_Fmax, "%.02f");
    op_odr_range_from_pct_Fmax_changed(data);
    op_odr_range_to_pct_Fmax_changed(data);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opodrcontrols->toggle_radial), TRUE);
    entry_set_gdouble_format(opodrcontrols->radial_angle, RADIAL_ANGLE_DEFAULT_OPODR, "%.1f");
}

static void op_odr_range_from_changed(Data *data)
{
    OPODRControls *fcontrols;
    OPODRdata *fdata;

    fcontrols = data->controls->opodrcontrols;
    fdata = &(data->args->opodrdata);

    range_from_changed(fcontrols->from, fcontrols->to, fcontrols->selection, &(fdata->from), &(fdata->from));
}

static void op_odr_range_to_changed(Data *data)
{
    OPODRControls *tcontrols;
    OPODRdata *tdata;

    tcontrols = data->controls->opodrcontrols;
    tdata = &(data->args->opodrdata);

    range_to_changed(tcontrols->from, tcontrols->to, tcontrols->selection, &(tdata->to), &(tdata->to));
}

static gboolean op_odr_range_from_pct_Fmax_changed(Data *data)
{
    OPODRControls *rcontrols;
    OPODRdata *rdata;

    rcontrols = data->controls->opodrcontrols;
    rdata = &(data->args->opodrdata);

    return range_from_pct_Fmax_changed(rcontrols->from_pct_Fmax, rcontrols->to_pct_Fmax, rcontrols->selection, data->args->fddata.Fmax,
                                       data->args->fddata.hunload, data->args->fddata.Funload,
                                       &(rdata->from_pct_Fmax), &(rdata->from_pct_Fmax), &(rdata->to_pct_Fmax), &(rdata->to_pct_Fmax),
                                       &(rdata->Finput), &(rdata->Finput));
}

static gboolean op_odr_range_to_pct_Fmax_changed(Data *data)
{
    OPODRControls *rcontrols;
    OPODRdata *rdata;

    rcontrols = data->controls->opodrcontrols;
    rdata = &(data->args->opodrdata);

    return range_to_pct_Fmax_changed(rcontrols->from_pct_Fmax, rcontrols->to_pct_Fmax, rcontrols->selection, data->args->fddata.Fmax,
                                     data->args->fddata.hunload, data->args->fddata.Funload,
                                     &(rdata->from_pct_Fmax), &(rdata->from_pct_Fmax), &(rdata->to_pct_Fmax), &(rdata->to_pct_Fmax),
                                     &(rdata->Finput), &(rdata->Finput));
}

static void op_odr_graph_selected(GwySelection *selection, gint i, Data *data)
{
    gdouble range[2];
    gint istart, iend, ndata;
    gdouble Fmax;

    g_return_if_fail(i <= 0);

    //take data from selection and put them in the "from" and "to" values
    gwy_selection_get_object(selection, 0, range);

    data->args->opodrdata.from = MIN(range[0], range[1]);
    data->args->opodrdata.to = MAX(range[0], range[1]);

    //update labels
    entry_set_gdouble_format(data->controls->opodrcontrols->from, data->args->opodrdata.from, DEPTH);
    entry_set_gdouble_format(data->controls->opodrcontrols->to, data->args->opodrdata.to, DEPTH);

    Fmax = data->args->fddata.Fmax;
    range_to_indices(data->args->opodrdata.from, data->args->opodrdata.to, data->args->fddata.hunload, TRUE, &istart, &iend, &ndata);
    data->args->opodrdata.from_pct_Fmax = data->args->fddata.Funload->data[iend] / Fmax * 100;
    data->args->opodrdata.to_pct_Fmax = data->args->fddata.Funload->data[istart] / Fmax * 100;

    entry_set_gdouble_format(data->controls->opodrcontrols->from_pct_Fmax, data->args->opodrdata.from_pct_Fmax, DEPTH);
    entry_set_gdouble_format(data->controls->opodrcontrols->to_pct_Fmax, data->args->opodrdata.to_pct_Fmax, DEPTH);

    data->args->opodrdata.Finput = FALSE;
}

void op_odr_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (op_odr_has_results(&(args->opodrdata))) {
        buffer = op_odr_export_data(&(args->opodrdata), &(args->fddata), &(args->instdata), &(args->area), EXPORT_FORMAT_PLAIN);
        buffer_to_file_dialog(buffer);
        g_free(buffer);
    }
}

void op_odr_fit_and_update(Data *data)
{
    gchar str[300];
    gint nc, i;
    gint istart, iend, ndata;
    gint fitinfo;
    gboolean fit = TRUE;
    GwyDataLine *x, *y;

    OPODRdata *opodrdata;
    FDdata *fddata;
    OPODRControls *opodrcontrols;
    GdkColor color;

    gdk_color_parse("red", &color);

    // TEMPORARY update GtkEntries for ranges
    /*
    if (data->args->opodrdata.Finput){
    	range_from_pct_Fmax_changed(data);
    	range_to_pct_Fmax_changed(data);
    }
    else{
    	range_from_changed(data);
    	range_to_changed(data);
    }
    */

    //create pointers for easier reading
    opodrdata = &(data->args->opodrdata);
    opodrcontrols = data->controls->opodrcontrols;
    fddata = &(data->args->fddata);

    //if there was a previous fit, invalidate previous fitted curve
    if (opodrdata->has_fit) {
        g_object_unref(opodrdata->xfit);
        g_object_unref(opodrdata->yfit);
    }

    // check range
    if (opodrdata->from == opodrdata->to) {
        fit = FALSE;
    }

    opodrdata->radial_corr = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opodrcontrols->toggle_radial));
    opodrdata->radial_angle = strtod(gtk_entry_get_text(GTK_ENTRY(opodrcontrols->radial_angle)), NULL);

    //find appropriate data, which should be fitted
    //default regime is governed by depth, can be switched to load regime by Finput (if the pct_Fmax ranges were used)
    if (opodrdata->Finput) {
        range_to_indices(opodrdata->from_pct_Fmax * 0.01 * fddata->Fmax, opodrdata->to_pct_Fmax * 0.01 * fddata->Fmax,
                         fddata->Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(opodrdata->from, opodrdata->to, fddata->hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range,no fitting can be performed
    ndata = MIN(ndata, fddata->hunload->res - istart);

    if (ndata < 2) {
        fit = FALSE;
    }

    opodrdata->nfitdata = ndata;

    if (fit) {
        //get data to be fitted
        x = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
        y = gwy_data_line_part_extract(fddata->Funload, istart, ndata);

        //fit the data and evaluate
        fitinfo = op_odr_fit(x, y, opodrdata, &(data->args->fddata), &(data->args->instdata), &(data->args->area), NULL);

        if (opodrdata->radial_corr) {
            radial_correction_op(opodrdata->radial_angle, &(opodrdata->Hit), &(opodrdata->Er), &(opodrdata->Eit), &(opodrdata->Aphc), &(data->args->instdata));
        }

        if (verbose) {
            g_print("fitinfo pre %d infolog %d \n", fitinfo, opodrdata->infolog);
        }

        /* fatal errors end with fitinfo >= 5, infolog first digit, fitinfo last digit  */
        if (fitinfo >= 1e6) {
            opodrdata->infolog = fitinfo / 1000000;
            fitinfo = fitinfo % 1000000;
        }

        if (verbose) {
            g_print("fitinfo post %d infolog %d \n", fitinfo, opodrdata->infolog);
        }

        if (fitinfo >= 5) {
            gtk_label_set_text(GTK_LABEL(opodrcontrols->fitinfo), "Warning: Questionable results \n or fatal errors in fitting!");
            gtk_widget_modify_fg(opodrcontrols->fitinfo, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(opodrcontrols->fitinfo), TRUE);
        }
        else if (fitinfo == 4) {
            gtk_label_set_text(GTK_LABEL(opodrcontrols->fitinfo), "Warning: Iteration limit reached!");
            gtk_widget_modify_fg(opodrcontrols->fitinfo, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(opodrcontrols->fitinfo), TRUE);
        }
        else if (opodrdata->m > 2) {
            gtk_label_set_text(GTK_LABEL(opodrcontrols->fitinfo), "Warning: Unphysical fit results!");
            gtk_widget_modify_fg(opodrcontrols->fitinfo, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(opodrcontrols->fitinfo), TRUE);
        }
        else {
            label_clear(opodrcontrols->fitinfo);
        }

        //update range labels
        entry_set_gdouble_format(opodrcontrols->from, opodrdata->from, DEPTH);
        entry_set_gdouble_format(opodrcontrols->to, opodrdata->to, DEPTH);
        entry_set_gdouble_format(opodrcontrols->from_pct_Fmax, opodrdata->from_pct_Fmax, DEPTH);
        entry_set_gdouble_format(opodrcontrols->to_pct_Fmax, opodrdata->to_pct_Fmax, DEPTH);

        //update labels
        label_set_gdouble_format(opodrcontrols->hp, opodrdata->hp, DEPTH);
        label_set_gdouble_format(opodrcontrols->m, opodrdata->m, NUMBER);
        label_set_gdouble_format(opodrcontrols->alpha, opodrdata->alpha, "%.4g");

        label_set_gdouble_format(opodrcontrols->SSres, opodrdata->SSres, "%.4g");
        label_set_gdouble_format(opodrcontrols->R2, opodrdata->R2, "%.4g");
        label_set_gdouble_format(opodrcontrols->R2adj, opodrdata->R2adj, "%.4g");
        label_set_gdouble_format(opodrcontrols->chi2, opodrdata->chi2, "%.4g");

        label_set_gdouble_format(opodrcontrols->eps, opodrdata->eps, NUMBER);
        label_set_gdouble_format(opodrcontrols->S, opodrdata->S, SLOPE);
        label_set_gdouble_format(opodrcontrols->hc, opodrdata->hc, DEPTH);
        label_set_gdouble_format(opodrcontrols->Aphc, opodrdata->Aphc, AREA);

        if (data->args->area.mode == AREA_DATA && opodrdata->hc > data->args->area.xmax) {
            gtk_label_set_text(GTK_LABEL(opodrcontrols->extrapol), "Warning: Extrapolation in area calibration!");
            gtk_widget_modify_fg(opodrcontrols->extrapol, GTK_STATE_NORMAL, &color);
        }
        else {
            label_clear(opodrcontrols->extrapol);
        }

        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", opodrdata->from, opodrdata->to);
        gtk_label_set_text(GTK_LABEL(opodrcontrols->range), str);
        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", opodrdata->from_pct_Fmax, opodrdata->to_pct_Fmax);
        gtk_label_set_text(GTK_LABEL(opodrcontrols->range_pct_Fmax), str);

        label_set_gdouble_format(opodrcontrols->Hit, opodrdata->Hit, HARD);
        label_set_gdouble_format(opodrcontrols->Eit, opodrdata->Eit, MODUL);
        label_set_gdouble_format(opodrcontrols->Er, opodrdata->Er, MODUL);

        if (opodrdata->radial_corr) {
            label_set_gdouble_format(opodrcontrols->radial_info, opodrdata->radial_angle, "Applied radial correction using α = %.4g°");
        }
        else {
            gtk_label_set_text(GTK_LABEL(opodrcontrols->radial_info), "Radial correction not applied");
        }

        //calculate curve of fit
        opodrdata->xfit = gwy_data_line_duplicate(x);
        opodrdata->yfit = gwy_data_line_new_alike(y, FALSE);

        for (i = 0; i < opodrdata->xfit->res; i++) {
            opodrdata->yfit->data[i] =  opodrdata->alpha * pow(opodrdata->xfit->data[i] - opodrdata->hp, opodrdata->m);
        }

        // clean up
        g_object_unref(x);
        g_object_unref(y);

        //update fitted curve, if necessary create it
        nc = gwy_graph_model_get_n_curves(opodrcontrols->graph_model);

        if (nc == 2) {
            gwy_graph_curve_model_set_data(opodrcontrols->cmodelSfit, opodrdata->xfit->data, opodrdata->yfit->data, opodrdata->xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(opodrcontrols->cmodelSfit);
#endif
        }
        else {
            gwy_graph_curve_model_set_data(opodrcontrols->cmodelSfit, opodrdata->xfit->data, opodrdata->yfit->data, opodrdata->xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(opodrcontrols->cmodelSfit);
#endif
            gwy_graph_model_add_curve(opodrcontrols->graph_model, opodrcontrols->cmodelSfit);
            g_object_set(opodrcontrols->cmodelSfit,
                         "mode", GWY_GRAPH_CURVE_LINE,
                         "description", "Fit",
                         "color", CURVE_COLOR_FIT,
                         "line-style", GDK_LINE_SOLID,
                         "line-width", 3,
                         NULL);
        }
    }
    else {
        nc = gwy_graph_model_get_n_curves(opodrcontrols->graph_model);
        op_odr_remove_all_fit_labels(opodrcontrols);
        op_odr_remove_all_fit_results(opodrdata);
        opodrdata->has_fit = FALSE;

        if (nc == 2) {
            gwy_graph_model_remove_curve_by_description(opodrcontrols->graph_model, "Fit");
        }
    }

    // log of fitting procedure is available
    if (opodrdata->has_fit) {
        gtk_widget_set_sensitive(opodrcontrols->button_log, TRUE);
    }
    else {
        gtk_widget_set_sensitive(opodrcontrols->button_log, FALSE);
    }

    // update previous uncertainties if existing, set button
    if (opodrdata->has_fit) {
        if (data->args->opodrunc.Ec != NULL || data->args->opodrunc.Hc != NULL || data->args->opodrunc.Ac != NULL || data->args->opodrunc.hc != NULL) {
            op_odr_unc_close_window(data);
        }

        gtk_widget_set_sensitive(opodrcontrols->button_unc, TRUE);
    }
    else {
        gtk_widget_set_sensitive(opodrcontrols->button_unc, FALSE);
    }
}

static void op_odr_show_log(Data *data)
{
    gint infolog;

    infolog = data->args->opodrdata.infolog;

    /* check if there are logs */
    switch (infolog) {
    case 0:
        show_log(data->controls->opodrcontrols->log_error_window, data->args->opodrdata.logfnm, ".err", "Fit error log");
        show_log(data->controls->opodrcontrols->log_report_window, data->args->opodrdata.logfnm, ".rpt", "Fit report log");
        break;

    case 1:
        show_log(data->controls->opodrcontrols->log_error_window,  data->args->opodrdata.logfnm, ".err", "Fit error log");
        show_warning("Fit report log could not be created.");
        break;

    case 2:
        show_log(data->controls->opodrcontrols->log_report_window, data->args->opodrdata.logfnm, ".rpt", "Fit report log");
        show_warning("Fit error log could not be created.");
        break;

    case 3:
        show_warning("Fit error log could not be created.");
        show_warning("Fit report log could not be created.");
        break;

    default:
        g_printerr("Should not get here.\n");
        break;
    }
}

GtkWidget *tool_op_odr_create(Data *data)
{
    GwyGraphArea *area;
    Args *args;
    OPODRControls *opodrcontrols;
    gint i;
    //	gchar str[300];

    args = data->args;
    opodrcontrols = data->controls->opodrcontrols;

    // GraphArea, two curves can be displayed: loading F-d curve and the resulting fit
    opodrcontrols->graph_model = gwy_graph_model_new();
    g_object_set(opodrcontrols->graph_model, "title", "Unload", NULL);
    gwy_graph_model_set_axis_label(opodrcontrols->graph_model, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(opodrcontrols->graph_model, GTK_POS_LEFT, "F / mN");

    opodrcontrols->zoomstack = NULL;
    opodrcontrols->zoomstack = g_slist_prepend(opodrcontrols->zoomstack, opodrcontrols->graph_model);

    opodrcontrols->hbox_graph = gtk_hbox_new(FALSE, CTRLS_SPACING);

    opodrcontrols->graph = gwy_graph_new(opodrcontrols->graph_model);
    /* gtk_widget_set_size_request(GTK_WIDGET(opodrcontrols->graph), 600, 400); */

    gwy_graph_enable_user_input(GWY_GRAPH(opodrcontrols->graph), FALSE);
    gwy_graph_set_status(GWY_GRAPH(opodrcontrols->graph), GWY_GRAPH_STATUS_XSEL);

    area = GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(opodrcontrols->graph)));
    opodrcontrols->selection = gwy_graph_area_get_selection(area, GWY_GRAPH_STATUS_XSEL);
    gwy_selection_set_max_objects(opodrcontrols->selection, 1);
    g_signal_connect(opodrcontrols->selection, "changed", G_CALLBACK(op_odr_graph_selected), data);

    opodrcontrols->cmodelunload = gwy_graph_curve_model_new();
    opodrcontrols->cmodelSfit = gwy_graph_curve_model_new();

    g_object_set(opodrcontrols->cmodelunload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Unloading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_UNLOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);
    gwy_graph_model_add_curve(opodrcontrols->graph_model, opodrcontrols->cmodelunload);

    // Table of input data and results
    //opodrcontrols->ctrls = gtk_table_new(13,4,TRUE);
    opodrcontrols->ctrls = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* Info section */

    opodrcontrols->frame_info = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(opodrcontrols->frame_info), GTK_WIDGET(gwy_label_new_header("Info")));
    opodrcontrols->table_info = gtk_table_new(2, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(opodrcontrols->table_info), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(opodrcontrols->table_info), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(opodrcontrols->table_info), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(opodrcontrols->frame_info), opodrcontrols->table_info);

    i = 0;

    /* h_max */

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_info), label_new_with_markup_left("h<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    opodrcontrols->hmax = gtk_label_new(NULL);
    label_set_gdouble_format(opodrcontrols->hmax, args->fddata.hmax, DEPTH);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->hmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_info), opodrcontrols->hmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_info), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /* F_max */

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_info), label_new_with_markup_left("F<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    opodrcontrols->Fmax = gtk_label_new(NULL);
    label_set_gdouble_format(opodrcontrols->Fmax, args->fddata.Fmax, FORCE);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_info), opodrcontrols->Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_info), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    /* Parameters section */

    opodrcontrols->frame_parameters = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(opodrcontrols->frame_parameters), GTK_WIDGET(gwy_label_new_header("Parameters")));
    opodrcontrols->table_parameters = gtk_table_new(4, 4, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(opodrcontrols->table_parameters), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(opodrcontrols->table_parameters), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(opodrcontrols->table_parameters), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(opodrcontrols->frame_parameters), opodrcontrols->table_parameters);
    i = 0;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_parameters), label_new_with_markup_left("Fit range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->from = gtk_entry_new();
    g_object_set_data(G_OBJECT(opodrcontrols->from), "id", (gpointer)"from");
    gtk_entry_set_width_chars(GTK_ENTRY(opodrcontrols->from), 8);
    gtk_entry_set_text(GTK_ENTRY(opodrcontrols->from), "");
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_parameters), opodrcontrols->from, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(opodrcontrols->from, TRUE);
    g_signal_connect_swapped(opodrcontrols->from, "activate", G_CALLBACK(op_odr_range_from_changed), data);
    /*	g_signal_connect_swapped(opodrcontrols->from, "focus-out-event", G_CALLBACK(range_from_changed),data); */

    opodrcontrols->to = gtk_entry_new();
    g_object_set_data(G_OBJECT(opodrcontrols->to), "id", (gpointer)"to");
    gtk_entry_set_width_chars(GTK_ENTRY(opodrcontrols->to), 8);
    gtk_entry_set_text(GTK_ENTRY(opodrcontrols->to), "");
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_parameters), opodrcontrols->to, 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(opodrcontrols->to, TRUE);
    g_signal_connect_swapped(opodrcontrols->to, "activate", G_CALLBACK(op_odr_range_to_changed), data);
    /*	g_signal_connect_swapped(opodrcontrols->to, "focus-out-event", G_CALLBACK(range_to_changed),data); */

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_parameters), label_new_with_markup_left("nm"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    opodrcontrols->from_pct_Fmax = gtk_entry_new();
    g_object_set_data(G_OBJECT(opodrcontrols->from_pct_Fmax), "id", (gpointer)"from_pct_Fmax");
    gtk_entry_set_width_chars(GTK_ENTRY(opodrcontrols->from_pct_Fmax), 8);
    gtk_entry_set_text(GTK_ENTRY(opodrcontrols->from_pct_Fmax), "");
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_parameters), opodrcontrols->from_pct_Fmax, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(opodrcontrols->from_pct_Fmax, TRUE);
    g_signal_connect_swapped(opodrcontrols->from_pct_Fmax, "activate", G_CALLBACK(op_odr_range_from_pct_Fmax_changed), data);
    /*	g_signal_connect_swapped(opodrcontrols->from_pct_Fmax, "focus-out-event", G_CALLBACK(range_from_pct_Fmax_changed),data); */

    opodrcontrols->to_pct_Fmax = gtk_entry_new();
    g_object_set_data(G_OBJECT(opodrcontrols->to_pct_Fmax), "id", (gpointer)"to_pct_Fmax");
    gtk_entry_set_width_chars(GTK_ENTRY(opodrcontrols->to_pct_Fmax), 8);
    gtk_entry_set_text(GTK_ENTRY(opodrcontrols->to_pct_Fmax), "");
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_parameters), opodrcontrols->to_pct_Fmax, 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(opodrcontrols->to_pct_Fmax, TRUE);
    g_signal_connect_swapped(opodrcontrols->to_pct_Fmax, "activate", G_CALLBACK(op_odr_range_to_pct_Fmax_changed), data);
    /*	g_signal_connect_swapped(opodrcontrols->to_pct_Fmax, "focus-out-event", G_CALLBACK(range_to_pct_Fmax_changed),data); */

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_parameters), label_new_with_markup_left("% F<sub>max</sub>"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /* beta */
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_parameters), label_new_with_markup_left("β"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->beta = gtk_entry_new();
    //	g_object_set_data(G_OBJECT(opodrcontrols->beta), "id", (gpointer)"beta");
    gtk_entry_set_width_chars(GTK_ENTRY(opodrcontrols->beta), 8);
    entry_set_gdouble_format(opodrcontrols->beta, data->args->opodrdata.beta, "%.03f");
    gwy_widget_set_activate_on_unfocus(opodrcontrols->beta, TRUE);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_parameters), opodrcontrols->beta, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(opodrcontrols->beta, "activate", G_CALLBACK(opodr_beta_changed), data);
    i++;

    opodrcontrols->hbox_radial = gtk_hbox_new(FALSE, 0); /* not using CTRLS_SPACING here */
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_parameters), opodrcontrols->hbox_radial, 0, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    opodrcontrols->toggle_radial = gtk_check_button_new_with_label("Radial correction using α = ");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opodrcontrols->toggle_radial), FALSE);
    g_signal_connect(opodrcontrols->toggle_radial, "toggled", G_CALLBACK(radial_toggled), data);
    gtk_box_pack_start(GTK_BOX(opodrcontrols->hbox_radial), opodrcontrols->toggle_radial, FALSE, FALSE, 0);

    opodrcontrols->radial_angle = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(opodrcontrols->radial_angle), 5);
    gtk_entry_set_max_length(GTK_ENTRY(opodrcontrols->radial_angle), 5);
    entry_set_gdouble_format(opodrcontrols->radial_angle, data->args->opodrdata.radial_angle, "%.4g");
    gtk_widget_set_sensitive(opodrcontrols->radial_angle, FALSE);
    gtk_box_pack_start(GTK_BOX(opodrcontrols->hbox_radial), opodrcontrols->radial_angle, FALSE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(opodrcontrols->hbox_radial), gtk_label_new("°"), FALSE, FALSE, 0);
    i++;


    opodrcontrols->reset_beta = gtk_button_new_with_label("Set defaults");
    gtk_button_set_image(GTK_BUTTON(opodrcontrols->reset_beta), gtk_image_new_from_stock(GTK_STOCK_REVERT_TO_SAVED, GTK_ICON_SIZE_BUTTON));
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_parameters), opodrcontrols->reset_beta, 0, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(opodrcontrols->reset_beta, "clicked", G_CALLBACK(opodr_reset_defaults), data);

    /* Fit button */

    opodrcontrols->button_fit = gtk_button_new_with_mnemonic(gwy_sgettext("verb|_Fit"));
    gtk_button_set_image(GTK_BUTTON(opodrcontrols->button_fit), gtk_image_new_from_stock(GTK_STOCK_EXECUTE, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(opodrcontrols->button_fit, "clicked", G_CALLBACK(op_odr_fit_and_update), data);


    /* Results section */

    opodrcontrols->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(opodrcontrols->frame_results), GTK_WIDGET(gwy_label_new_header("Results")));

    opodrcontrols->scrolledwindow_results = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(opodrcontrols->scrolledwindow_results), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(opodrcontrols->frame_results), opodrcontrols->scrolledwindow_results);

    opodrcontrols->table_results = gtk_table_new(13, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(opodrcontrols->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(opodrcontrols->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(opodrcontrols->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(opodrcontrols->scrolledwindow_results), opodrcontrols->table_results);
    i = 0;

    opodrcontrols->fitinfo = gtk_label_new(NULL);
    gtk_label_set_line_wrap(GTK_LABEL(opodrcontrols->fitinfo), TRUE);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->fitinfo, 0, 3, i, i + 2, GTK_FILL, GTK_FILL, 0, 0);
    i += 2;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("h<sub>p</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->hp = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->hp), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->hp, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("m"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->m = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->m), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->m, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("α"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->alpha = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->alpha), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->alpha, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    //	g_snprintf(str, sizeof(str), "mN.nm<sup>-%.2g</sup>", data->args->opodrdata.m);
    //	gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left(str), 2, 3, i, i+1,  GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("SS<sub>res</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->SSres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->SSres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->SSres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("R<sup>2</sup>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->R2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->R2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->R2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("R<sup>2</sup><sub>adj</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->R2adj = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->R2adj), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->R2adj, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("χ<sup>2</sup> fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->chi2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->chi2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->chi2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("ε"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->eps = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->eps), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->eps, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("h<sub>c</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->hc = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->hc), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->hc, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("S"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->S = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->S), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->S, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("mN/nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("A<sub>p</sub>(h<sub>c</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->Aphc = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->Aphc), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->Aphc, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("nm<sup>2</sup>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    opodrcontrols->extrapol = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->extrapol, 0, 3, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /* R not available yet */

    /* gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("R fit"), 0, 1, i, i+1, GTK_FILL, 0, 0, 0);  */
    /* opodrcontrols->R = gtk_label_new(NULL);  */
    /* gtk_misc_set_alignment(GTK_MISC(opodrcontrols->R), 0.0, 0.5);  */
    /* gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->R, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);  */
    /* i++; */

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("H<sub>IT</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->Hit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->Hit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->Hit, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("MPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("E<sub>r</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->Er = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->Er), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->Er, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("E<sub>IT</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->Eit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->Eit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->Eit, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("Fitted range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    opodrcontrols->range = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->range), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->range, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    opodrcontrols->range_pct_Fmax = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->range_pct_Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->range_pct_Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), label_new_with_markup_left("% F<sub>max</sub>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    opodrcontrols->radial_info = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(opodrcontrols->radial_info), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(opodrcontrols->table_results), opodrcontrols->radial_info, 0, 3, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    /* Log, uncertainties, save buttons */

    opodrcontrols->button_log = gtk_button_new_with_mnemonic("Show _log");
    gtk_button_set_image(GTK_BUTTON(opodrcontrols->button_log), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(opodrcontrols->button_log, "clicked", G_CALLBACK(op_odr_show_log), data);
    gtk_widget_set_sensitive(opodrcontrols->button_log, FALSE);
    opodrcontrols->log_error_window = NULL;
    opodrcontrols->log_report_window = NULL;

    opodrcontrols->button_unc = gtk_button_new_with_mnemonic("_Uncertainties");
    gtk_button_set_image(GTK_BUTTON(opodrcontrols->button_unc), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(opodrcontrols->button_unc, "clicked", G_CALLBACK(op_odr_uncertainty), data);
    gtk_widget_set_sensitive(opodrcontrols->button_unc, FALSE);

    opodrcontrols->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(opodrcontrols->button_save, "clicked", G_CALLBACK(op_odr_save_data), args);


    gtk_box_pack_start(GTK_BOX(opodrcontrols->ctrls), opodrcontrols->frame_info, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(opodrcontrols->ctrls), opodrcontrols->frame_parameters, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(opodrcontrols->ctrls), opodrcontrols->button_fit, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(opodrcontrols->ctrls), opodrcontrols->frame_results, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(opodrcontrols->ctrls), opodrcontrols->button_log, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(opodrcontrols->ctrls), opodrcontrols->button_unc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(opodrcontrols->ctrls), opodrcontrols->button_save, FALSE, FALSE, 0);

    /* Zooming */

    opodrcontrols->vbox_zoom_buttons = gtk_vbox_new(TRUE, 0);
    opodrcontrols->vbox_zoom_outer = gtk_vbox_new(FALSE, 0);

    opodrcontrols->button_zoom_in = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(opodrcontrols->button_zoom_in), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(opodrcontrols->vbox_zoom_buttons), opodrcontrols->button_zoom_in, TRUE, TRUE, 0);
    g_signal_connect_swapped(opodrcontrols->button_zoom_in, "clicked", G_CALLBACK(op_odr_zoom_in), data);

    opodrcontrols->button_zoom_out = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(opodrcontrols->button_zoom_out), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(opodrcontrols->vbox_zoom_buttons), opodrcontrols->button_zoom_out, TRUE, TRUE, 0);
    g_signal_connect_swapped(opodrcontrols->button_zoom_out, "clicked", G_CALLBACK(op_odr_zoom_out), data);

    opodrcontrols->button_zoom_restore = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(opodrcontrols->button_zoom_restore), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(opodrcontrols->vbox_zoom_buttons), opodrcontrols->button_zoom_restore, TRUE, TRUE, 0);
    g_signal_connect_swapped(opodrcontrols->button_zoom_restore, "clicked", G_CALLBACK(op_odr_zoom_restore), data);

    gtk_box_pack_start(GTK_BOX(opodrcontrols->vbox_zoom_outer), opodrcontrols->vbox_zoom_buttons, TRUE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(opodrcontrols->hbox_graph), opodrcontrols->graph, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(opodrcontrols->hbox_graph), opodrcontrols->vbox_zoom_outer, FALSE, FALSE, 0);


    //zoom and unzoom buttons
    /* opodrcontrols->button_zoom_in = gtk_button_new(); */
    /* gtk_button_set_image (GTK_BUTTON (opodrcontrols->button_zoom_in), gtk_image_new_from_stock (GWY_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON)); */
    /* g_signal_connect_swapped(opodrcontrols->button_zoom_in, "clicked", G_CALLBACK (zoom_in), data); */

    /* opodrcontrols->button_zoom_out = gtk_button_new(); */
    /* gtk_button_set_image (GTK_BUTTON (opodrcontrols->button_zoom_out), gtk_image_new_from_stock (GWY_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON)); */
    /* g_signal_connect_swapped(opodrcontrols->button_zoom_out, "clicked", G_CALLBACK (zoom_out), data); */

    /* vbox = gtk_vbox_new(FALSE,FALSE); */
    /* gtk_box_pack_start (GTK_BOX(vbox), opodrcontrols->graph, TRUE, TRUE, 0); */
    /* gtk_box_pack_start (GTK_BOX(vbox), opodrcontrols->button_zoom_in, FALSE, FALSE, 0); */
    /* gtk_box_pack_start (GTK_BOX(vbox), opodrcontrols->button_zoom_out, FALSE, FALSE, 0); */


    //final assembly

    opodrcontrols->opodr_gui = gtk_hbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(opodrcontrols->opodr_gui), opodrcontrols->ctrls, FALSE, FALSE, 0);
    //	gtk_box_pack_start (GTK_BOX (opodrcontrols->opodr_gui), opodrcontrols->graph, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(opodrcontrols->opodr_gui), opodrcontrols->hbox_graph, TRUE, TRUE, 0);

    gtk_widget_show(opodrcontrols->opodr_gui);

    return opodrcontrols->opodr_gui;
}

static void radial_toggled(GtkToggleButton *button, Data *data)
{
    if (gtk_toggle_button_get_active(button)) {
        gtk_widget_set_sensitive(data->controls->opodrcontrols->radial_angle, TRUE);
    }
    else {
        gtk_widget_set_sensitive(data->controls->opodrcontrols->radial_angle, FALSE);
    }
}

void op_odr_recalculate(Data *data)
{
    recalculate_Eit(data->controls->opodrcontrols, &(data->args->opodrdata), &(data->args->instdata));
    op_odr_propagate_uncertainties(data->controls->opodrunccontrols, &(data->args->opodrdata), &(data->args->opodrunc),
                                   &(data->args->fddata), &(data->args->area));
}

static void update_models_op_odr(OPODRControls *opodrcontrols)
{
    GwyGraphCurveModel *cmodel;

    opodrcontrols->graph_model = gwy_graph_get_model(GWY_GRAPH(opodrcontrols->graph));
    cmodel = gwy_graph_model_get_curve_by_description(opodrcontrols->graph_model, "Unloading");

    if (cmodel) {
        opodrcontrols->cmodelunload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(opodrcontrols->graph_model, "Fit");

    if (cmodel) {
        opodrcontrols->cmodelSfit = cmodel;
    }
}

static void op_odr_zoom_in(Data *data)
{
    OPODRControls *opodrcontrols;

    opodrcontrols = data->controls->opodrcontrols;

    rs_zoom_in(GWY_GRAPH(opodrcontrols->graph), &(opodrcontrols->zoomstack), FALSE);   /* do not rescale y */
    update_models_op_odr(opodrcontrols);
}

static void op_odr_zoom_out(Data *data)
{
    OPODRControls *opodrcontrols;

    opodrcontrols = data->controls->opodrcontrols;

    rs_zoom_out(GWY_GRAPH(opodrcontrols->graph), &(opodrcontrols->zoomstack));
    update_models_op_odr(opodrcontrols);
}

void op_odr_zoom_restore(Data *data)
{
    OPODRControls *opodrcontrols;

    opodrcontrols = data->controls->opodrcontrols;

    rs_zoom_restore(GWY_GRAPH(opodrcontrols->graph), &(opodrcontrols->zoomstack));
    update_models_op_odr(opodrcontrols);
}

#endif

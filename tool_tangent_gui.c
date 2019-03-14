#include "tool_tangent_gui.h"
#include "tool_tangent.h"
#include "tool_tangent_unc_gui.h"

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


static void tangent_remove_all_fit_labels(TangentControls *tgcontrols);
static void tangent_remove_all_fit_curves(TangentControls *tgcontrols);

static void update_models_tg(TangentControls *tgcontrols);
static void tg_zoom_in(Data *data);
static void tg_zoom_out(Data *data);

void tangent_remove_fit_results_labels(Data *data)
{
    gwy_selection_clear(data->controls->tgcontrols->selection);
    tangent_remove_all_fit_results(&(data->args->tgdata));
    tangent_remove_all_fit_labels(data->controls->tgcontrols);
    tangent_remove_all_fit_curves(data->controls->tgcontrols);

    if (data->args->tgunc.Ec != NULL || data->args->tgunc.Hc != NULL) {
        tangent_unc_close_window(data);
    }
}

static void tangent_remove_all_fit_labels(TangentControls *tgcontrols)
{
    /* remove all labels that result from a fit */
    entry_clear(tgcontrols->from);
    entry_clear(tgcontrols->to);
    entry_clear(tgcontrols->from_pct_Fmax);
    entry_clear(tgcontrols->to_pct_Fmax);

    label_clear(tgcontrols->chi2);
    label_clear(tgcontrols->SSres);
    label_clear(tgcontrols->R2);
    label_clear(tgcontrols->R2adj);
    label_clear(tgcontrols->R);
    label_clear(tgcontrols->S);
    label_clear(tgcontrols->hc);
    label_clear(tgcontrols->Aphc);
    label_clear(tgcontrols->extrapol);

    label_clear(tgcontrols->Hit);
    label_clear(tgcontrols->Er);
    label_clear(tgcontrols->Eit);
    label_clear(tgcontrols->range);
    label_clear(tgcontrols->range_pct_Fmax);
}

static void tangent_remove_all_fit_curves(TangentControls *tgcontrols)
{
    gwy_graph_model_remove_curve_by_description(tgcontrols->graph_model, "Fit");
}

static void recalculate_Eit(TangentControls *tgcontrols, Tangentdata *tgdata, const Instdata *instdata)
{
    gdouble Er;

    Er = tgdata->Er * 1e9;

    // recalc for changed  nu, nui or Ei
    if (Er == 0) {
        return;
    }

    tgdata->Eit = calc_Eit(Er, instdata->nu, instdata->Ei, instdata->nui) * 1e-9;
    label_set_gdouble_format(tgcontrols->Eit, tgdata->Eit, MODUL);
}

void tangent_redraw(Data *data)
{
    // update loading curve and Fmax and hmax labels
    if (data->args->fddata.hunload) {
        gwy_graph_curve_model_set_data(data->controls->tgcontrols->cmodelunload,
                                       data->args->fddata.hunload->data, data->args->fddata.Funload->data,
                                       data->args->fddata.hunload->res);

        label_set_gdouble_format(data->controls->tgcontrols->hmax, data->args->fddata.hmax, DEPTH);
        label_set_gdouble_format(data->controls->tgcontrols->Fmax, data->args->fddata.Fmax, FORCE);

        if (gwy_graph_model_get_n_curves(data->controls->tgcontrols->graph_model) == 2) {
            gwy_graph_curve_model_set_data(data->controls->tgcontrols->cmodelSfit,
                                           data->args->tgdata.xfit->data, data->args->tgdata.yfit->data,  data->args->tgdata.xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(data->controls->tgcontrols->cmodelSfit);
#endif
        }
    }
}

static void tangent_recalculate_beta(TangentControls *tgcontrols, Tangentdata *tgdata, const Instdata *instdata)
{
    if (!tgdata->has_fit) {
        return;
    }

    tgdata->Er = calc_Er_OP(tgdata->Aphc, tgdata->S, tgdata->beta) * 1e15;
    tgdata->Eit = calc_Eit(tgdata->Er, instdata->nu, instdata->Ei, instdata->nui);
    tgdata->Eit *= 1e-9;
    tgdata->Er *= 1e-9;

    label_set_gdouble_format(tgcontrols->Eit, tgdata->Eit, MODUL);
    label_set_gdouble_format(tgcontrols->Er, tgdata->Er, MODUL);
}

static void tangent_beta_changed(Data *data)
{
    TangentControls *fcontrols;

    fcontrols = data->controls->tgcontrols;

    if (beta_changed(fcontrols->beta, &(data->args->tgdata.beta))) {
        tangent_recalculate_beta(data->controls->tgcontrols, &(data->args->tgdata), &(data->args->instdata));
    }
}

static void tangent_reset_beta(Data *data)
{
    data->args->tgdata.beta =  BETA_DEFAULT_TG;
    entry_set_gdouble_format(data->controls->tgcontrols->beta, data->args->tgdata.beta, "%.03f");
    tangent_recalculate_beta(data->controls->tgcontrols, &(data->args->tgdata), &(data->args->instdata));
}

static void tangent_range_from_changed(Data *data)
{
    TangentControls *fcontrols;
    Tangentdata *fdata;

    fcontrols = data->controls->tgcontrols;
    fdata = &(data->args->tgdata);

    range_from_changed(fcontrols->from, fcontrols->to, fcontrols->selection, &(fdata->from), &(fdata->from));
}

static void tangent_range_to_changed(Data *data)
{
    TangentControls *tcontrols;
    Tangentdata *tdata;

    tcontrols = data->controls->tgcontrols;
    tdata = &(data->args->tgdata);

    range_to_changed(tcontrols->from, tcontrols->to, tcontrols->selection, &(tdata->to), &(tdata->to));
}

static gboolean tangent_range_from_pct_Fmax_changed(Data *data)
{
    TangentControls *rcontrols;
    Tangentdata *rdata;

    rcontrols = data->controls->tgcontrols;
    rdata = &(data->args->tgdata);

    return range_from_pct_Fmax_changed(rcontrols->from_pct_Fmax, rcontrols->to_pct_Fmax, rcontrols->selection, data->args->fddata.Fmax,
				       data->args->fddata.hunload, data->args->fddata.Funload,
                                       &(rdata->from_pct_Fmax), &(rdata->from_pct_Fmax), &(rdata->to_pct_Fmax), &(rdata->to_pct_Fmax),
				       &(rdata->Finput), &(rdata->Finput));
}

static gboolean tangent_range_to_pct_Fmax_changed(Data *data)
{
    TangentControls *rcontrols;
    Tangentdata *rdata;

    rcontrols = data->controls->tgcontrols;
    rdata = &(data->args->tgdata);

    return range_to_pct_Fmax_changed(rcontrols->from_pct_Fmax, rcontrols->to_pct_Fmax, rcontrols->selection, data->args->fddata.Fmax,
				     data->args->fddata.hunload, data->args->fddata.Funload,
                                     &(rdata->from_pct_Fmax), &(rdata->from_pct_Fmax), &(rdata->to_pct_Fmax), &(rdata->to_pct_Fmax),
				     &(rdata->Finput), &(rdata->Finput));
}

static void tangent_graph_selected(GwySelection *selection, gint i, Data *data)
{
    gdouble range[2];
    gint istart, iend, ndata;
    gdouble Fmax;

    g_return_if_fail(i <= 0);

    //take data from selection and put them in the "from" and "to" values
    gwy_selection_get_object(selection, 0, range);

    data->args->tgdata.from = MIN(range[0], range[1]);
    data->args->tgdata.to = MAX(range[0], range[1]);

    //update labels
    entry_set_gdouble_format(data->controls->tgcontrols->from, data->args->tgdata.from, DEPTH);
    entry_set_gdouble_format(data->controls->tgcontrols->to, data->args->tgdata.to, DEPTH);

    Fmax = data->args->fddata.Fmax;
    range_to_indices(data->args->tgdata.from, data->args->tgdata.to, data->args->fddata.hunload, TRUE, &istart, &iend, &ndata);
    data->args->tgdata.from_pct_Fmax = data->args->fddata.Funload->data[iend] / Fmax * 100;
    data->args->tgdata.to_pct_Fmax = data->args->fddata.Funload->data[istart] / Fmax * 100;

    entry_set_gdouble_format(data->controls->tgcontrols->from_pct_Fmax, data->args->tgdata.from_pct_Fmax, DEPTH);
    entry_set_gdouble_format(data->controls->tgcontrols->to_pct_Fmax, data->args->tgdata.to_pct_Fmax, DEPTH);

    data->args->tgdata.Finput = FALSE;
}

static void tangent_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (tangent_has_results(&(args->tgdata))) {
	buffer = tangent_export_data(&(args->tgdata), &(args->fddata), &(args->instdata), &(args->area), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

void tangent_fit_and_update(Data *data)
{
    LinReg reg;
    gchar str[300];
    gint nc;
    gint istart, iend, ndata;
    gdouble h1;
    gboolean fit = TRUE;
    GwyDataLine *x, *y;

    Tangentdata *tgdata;
    FDdata *fddata;
    TangentControls *tgcontrols;
    GdkColor color;

    gdk_color_parse("red", &color);

    // TEMPORARY update GtkEntries for ranges
    /*
    if (data->args->tgdata.Finput){
    	range_from_pct_Fmax_changed(data);
    	range_to_pct_Fmax_changed(data);
    }
    else{
    	range_from_changed(data);
    	range_to_changed(data);
    }
    */

    //create pointers for easier reading
    tgdata = &(data->args->tgdata);
    tgcontrols = data->controls->tgcontrols;
    fddata = &(data->args->fddata);

    //if there was a previous fit, invalidate previous fitted curve
    if (tgdata->has_fit) {
        g_object_unref(tgdata->xfit);
        g_object_unref(tgdata->yfit);
    }

    //check range
    if (tgdata->from == tgdata->to) {
        fit = FALSE;
    }

    //find appropriate data, which should be fitted
    //default regime is governed by depth, can be switched to load regime by Finput (if the pct_Fmax ranges were used)
    if (tgdata->Finput) {
        range_to_indices(tgdata->from_pct_Fmax * 0.01 * fddata->Fmax, tgdata->to_pct_Fmax * 0.01 * fddata->Fmax,
                         fddata->Funload, TRUE, &istart, &iend, &ndata);
    }
    else {
        range_to_indices(tgdata->from, tgdata->to, fddata->hunload, TRUE, &istart, &iend, &ndata);
    }

    // if there are not enough data in the selected range, no fitting can be performed
    ndata = MIN(ndata, fddata->hunload->res - istart);

    if (ndata < 2) {
        fit = FALSE;
    }

    tgdata->nfitdata = ndata;

    if (fit) {
        //get data to be fitted
        x = gwy_data_line_part_extract(fddata->hunload, istart, ndata);
        y = gwy_data_line_part_extract(fddata->Funload, istart, ndata);
        //
        //fit the data and evaluate
        tangent_fit(x, y, tgdata, &(data->args->fddata), &(data->args->instdata), &(data->args->area), 0, 0);
        reg = tgdata->reg;

        //update range labels
        entry_set_gdouble_format(tgcontrols->from, tgdata->from, DEPTH);
        entry_set_gdouble_format(tgcontrols->to, tgdata->to, DEPTH);
        entry_set_gdouble_format(tgcontrols->from_pct_Fmax, tgdata->from_pct_Fmax, DEPTH);
        entry_set_gdouble_format(tgcontrols->to_pct_Fmax, tgdata->to_pct_Fmax, DEPTH);

        //update labels
        /*
          label_set_gdouble_format(tgcontrols->chi2, reg.chi, NUMBER);
          label_set_gdouble_format(tgcontrols->R, reg.R, NUMBER);
        */

        label_set_gdouble_format(tgcontrols->SSres, tgdata->reg.SSres, "%.4g");
        label_set_gdouble_format(tgcontrols->R2, tgdata->reg.R2, "%.4g");
        label_set_gdouble_format(tgcontrols->R2adj, tgdata->reg.R2adj, "%.4g");
        label_set_gdouble_format(tgcontrols->chi2, tgdata->reg.chi2, "%.4g");
        label_set_gdouble_format(tgcontrols->S, tgdata->S, SLOPE);
        label_set_gdouble_format(tgcontrols->hc, tgdata->hc, DEPTH);
        label_set_gdouble_format(tgcontrols->Aphc, tgdata->Aphc, AREA);

        if (data->args->area.mode == AREA_DATA && tgdata->hc > data->args->area.xmax) {
            gtk_label_set_text(GTK_LABEL(tgcontrols->extrapol), "Warning: Extrapolation in area calibration!");
            gtk_widget_modify_fg(tgcontrols->extrapol, GTK_STATE_NORMAL, &color);
        }
        else {
            label_clear(tgcontrols->extrapol);
        }

        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", tgdata->from, tgdata->to);
        gtk_label_set_text(GTK_LABEL(tgcontrols->range), str);
        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", tgdata->from_pct_Fmax, tgdata->to_pct_Fmax);
        gtk_label_set_text(GTK_LABEL(tgcontrols->range_pct_Fmax), str);

        label_set_gdouble_format(tgcontrols->Hit, tgdata->Hit, HARD);
        label_set_gdouble_format(tgcontrols->Eit, tgdata->Eit, MODUL);
        label_set_gdouble_format(tgcontrols->Er, tgdata->Er, MODUL);

        // create fit curve data
        // the straight line corresponding to the fit will be shown going up to the intersection with the y-axis, find appropriate data
        h1 = -tgdata->reg.intercept / tgdata->reg.slope;

        while (iend < fddata->hunload->res && fddata->hunload->data[iend] > h1) {
            iend++;
        }

        ndata = iend + 2;
        ndata = MIN(ndata, fddata->hunload->res);
        tgdata->xfit = gwy_data_line_part_extract(fddata->hunload, 0, ndata);
        tgdata->yfit = gwy_data_line_duplicate(tgdata->xfit);
        gwy_data_line_multiply(tgdata->yfit, reg.slope);
        gwy_data_line_add(tgdata->yfit, reg.intercept);

        // clean up
        g_object_unref(x);
        g_object_unref(y);

        // update fitted curve, if necessary create it.
        nc = gwy_graph_model_get_n_curves(tgcontrols->graph_model);

        if (nc == 2) {
            gwy_graph_curve_model_set_data(tgcontrols->cmodelSfit, tgdata->xfit->data, tgdata->yfit->data, tgdata->xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(tgcontrols->cmodelSfit);
#endif
        }
        else {
            gwy_graph_curve_model_set_data(tgcontrols->cmodelSfit, tgdata->xfit->data, tgdata->yfit->data, tgdata->xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(tgcontrols->cmodelSfit);
#endif
            gwy_graph_model_add_curve(tgcontrols->graph_model, tgcontrols->cmodelSfit);
            g_object_set(tgcontrols->cmodelSfit,
                         "mode", GWY_GRAPH_CURVE_LINE,
                         "description", "Fit",
                         "color", CURVE_COLOR_FIT,
                         "line-style", GDK_LINE_SOLID,
                         "line-width", 3,
                         NULL);
        }
    }
    else {
        nc = gwy_graph_model_get_n_curves(tgcontrols->graph_model);
        tangent_remove_all_fit_labels(tgcontrols);
        tangent_remove_all_fit_results(tgdata);
        tgdata->has_fit = FALSE;

        if (nc == 2) {
            gwy_graph_model_remove_curve_by_description(tgcontrols->graph_model, "Fit");
        }
    }

    // update previous uncertainties if existing, set button
    if (tgdata->has_fit) {
        if (data->args->tgunc.Ec != NULL || data->args->tgunc.Hc != NULL) {
            tangent_unc_close_window(data);
        }

        gtk_widget_set_sensitive(tgcontrols->button_unc, TRUE);
    }
    else {
        gtk_widget_set_sensitive(tgcontrols->button_unc, FALSE);
    }
}

GtkWidget *tool_tangent_create(Data *data)
{
    GwyGraphArea *area;
    Args *args;
    TangentControls *tgcontrols;
    gint i;

    args = data->args;
    tgcontrols = data->controls->tgcontrols;

    // GraphArea, two curves can be displayed: loading F-d curve and the resulting fit (in this case a straight line going all the way to the y-axis, NOT the fitted data)
    tgcontrols->graph_model = gwy_graph_model_new();
    g_object_set(tgcontrols->graph_model, "title", "Unload", NULL);
    gwy_graph_model_set_axis_label(tgcontrols->graph_model, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(tgcontrols->graph_model, GTK_POS_LEFT, "F / mN");

    tgcontrols->zoomstack = NULL;
    tgcontrols->zoomstack = g_slist_prepend(tgcontrols->zoomstack, tgcontrols->graph_model);

    tgcontrols->hbox_graph = gtk_hbox_new(FALSE, CTRLS_SPACING);

    tgcontrols->graph = gwy_graph_new(tgcontrols->graph_model);
    /* gtk_widget_set_size_request(GTK_WIDGET(tgcontrols->graph), 600, 400); */

    gwy_graph_enable_user_input(GWY_GRAPH(tgcontrols->graph), FALSE);
    gwy_graph_set_status(GWY_GRAPH(tgcontrols->graph), GWY_GRAPH_STATUS_XSEL);

    area = GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(tgcontrols->graph)));
    tgcontrols->selection = gwy_graph_area_get_selection(area, GWY_GRAPH_STATUS_XSEL);
    gwy_selection_set_max_objects(tgcontrols->selection, 1);
    g_signal_connect(tgcontrols->selection, "changed", G_CALLBACK(tangent_graph_selected), data);

    tgcontrols->cmodelunload = gwy_graph_curve_model_new();
    tgcontrols->cmodelSfit = gwy_graph_curve_model_new();

    g_object_set(tgcontrols->cmodelunload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Unloading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_UNLOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);
    
    gwy_graph_model_add_curve(tgcontrols->graph_model, tgcontrols->cmodelunload);


    /* Ctrls */

    tgcontrols->ctrls = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* Info section */

    tgcontrols->frame_info = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(tgcontrols->frame_info), GTK_WIDGET(gwy_label_new_header("Info")));
    tgcontrols->table_info = gtk_table_new(2, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(tgcontrols->table_info), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(tgcontrols->table_info), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(tgcontrols->table_info), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(tgcontrols->frame_info), tgcontrols->table_info);

    i = 0;

    /* h_max */

    gtk_table_attach(GTK_TABLE(tgcontrols->table_info), label_new_with_markup_left("h<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    tgcontrols->hmax = gtk_label_new(NULL);
    label_set_gdouble_format(tgcontrols->hmax, args->fddata.hmax, DEPTH);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->hmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_info), tgcontrols->hmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(tgcontrols->table_info), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /* F_max */

    gtk_table_attach(GTK_TABLE(tgcontrols->table_info), label_new_with_markup_left("F<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    tgcontrols->Fmax = gtk_label_new(NULL);
    label_set_gdouble_format(tgcontrols->Fmax, args->fddata.Fmax, FORCE);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_info), tgcontrols->Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(tgcontrols->table_info), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    /* Parameters section */

    tgcontrols->frame_parameters = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(tgcontrols->frame_parameters), GTK_WIDGET(gwy_label_new_header("Parameters")));
    tgcontrols->table_parameters = gtk_table_new(5, 4, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(tgcontrols->table_parameters), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(tgcontrols->table_parameters), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(tgcontrols->table_parameters), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(tgcontrols->frame_parameters), tgcontrols->table_parameters);

    i = 0;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_parameters), label_new_with_markup_left("Fit range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->from = gtk_entry_new();
    g_object_set_data(G_OBJECT(tgcontrols->from), "id", (gpointer)"from");
    gtk_entry_set_width_chars(GTK_ENTRY(tgcontrols->from), 8);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_parameters), tgcontrols->from, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(tgcontrols->from, TRUE);
    g_signal_connect_swapped(tgcontrols->from, "activate", G_CALLBACK(tangent_range_from_changed), data);
    /*	g_signal_connect_swapped(tgcontrols->from, "focus-out-event", G_CALLBACK(range_from_changed),data); */

    tgcontrols->to = gtk_entry_new();
    g_object_set_data(G_OBJECT(tgcontrols->to), "id", (gpointer)"to");
    gtk_entry_set_width_chars(GTK_ENTRY(tgcontrols->to), 8);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_parameters), tgcontrols->to, 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_parameters), label_new_with_markup_left("nm"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(tgcontrols->to, TRUE);
    g_signal_connect_swapped(tgcontrols->to, "activate", G_CALLBACK(tangent_range_to_changed), data);
    /*	g_signal_connect_swapped(tgcontrols->to, "focus-out-event", G_CALLBACK(range_to_changed),data); */
    i++;

    tgcontrols->from_pct_Fmax = gtk_entry_new();
    g_object_set_data(G_OBJECT(tgcontrols->from_pct_Fmax), "id", (gpointer)"from_pct_Fmax");
    gtk_entry_set_width_chars(GTK_ENTRY(tgcontrols->from_pct_Fmax), 8);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_parameters), tgcontrols->from_pct_Fmax, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(tgcontrols->from_pct_Fmax, TRUE);
    g_signal_connect_swapped(tgcontrols->from_pct_Fmax, "activate", G_CALLBACK(tangent_range_from_pct_Fmax_changed), data);
    /*	g_signal_connect_swapped(tgcontrols->from_pct_Fmax, "focus-out-event", G_CALLBACK(range_from_pct_Fmax_changed),data); */

    tgcontrols->to_pct_Fmax = gtk_entry_new();
    g_object_set_data(G_OBJECT(tgcontrols->to_pct_Fmax), "id", (gpointer)"to_pct_Fmax");
    gtk_entry_set_width_chars(GTK_ENTRY(tgcontrols->to_pct_Fmax), 8);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_parameters), tgcontrols->to_pct_Fmax, 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(tgcontrols->to_pct_Fmax, TRUE);
    g_signal_connect_swapped(tgcontrols->to_pct_Fmax, "activate", G_CALLBACK(tangent_range_to_pct_Fmax_changed), data);
    /*	g_signal_connect_swapped(tgcontrols->to_pct_Fmax, "focus-out-event", G_CALLBACK(range_to_pct_Fmax_changed),data); */

    gtk_table_attach(GTK_TABLE(tgcontrols->table_parameters), label_new_with_markup_left("% F<sub>max</sub>"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_parameters), label_new_with_markup_left("β"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->beta = gtk_entry_new();
    //	g_object_set_data(G_OBJECT(tgcontrols->beta), "id", (gpointer)"beta");
    gtk_entry_set_width_chars(GTK_ENTRY(tgcontrols->beta), 8);
    entry_set_gdouble_format(tgcontrols->beta, data->args->tgdata.beta, "%.03f");
    gwy_widget_set_activate_on_unfocus(tgcontrols->beta, TRUE);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_parameters), tgcontrols->beta, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(tgcontrols->beta, "activate", G_CALLBACK(tangent_beta_changed), data);

    tgcontrols->reset_beta = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(tgcontrols->reset_beta), gtk_image_new_from_stock(GTK_STOCK_REVERT_TO_SAVED, GTK_ICON_SIZE_BUTTON));
    gtk_table_attach(GTK_TABLE(tgcontrols->table_parameters), tgcontrols->reset_beta, 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(tgcontrols->reset_beta, "clicked", G_CALLBACK(tangent_reset_beta), data);

    /* Fit button */
    tgcontrols->button_fit = gtk_button_new_with_mnemonic(gwy_sgettext("verb|_Fit"));
    gtk_button_set_image(GTK_BUTTON(tgcontrols->button_fit), gtk_image_new_from_stock(GTK_STOCK_EXECUTE, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(tgcontrols->button_fit, "clicked", G_CALLBACK(tangent_fit_and_update), data);

    /* Results section */

    tgcontrols->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(tgcontrols->frame_results), GTK_WIDGET(gwy_label_new_header("Results")));

    tgcontrols->scrolledwindow_results = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(tgcontrols->scrolledwindow_results), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(tgcontrols->frame_results), tgcontrols->scrolledwindow_results);

    tgcontrols->table_results = gtk_table_new(12, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(tgcontrols->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(tgcontrols->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(tgcontrols->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(tgcontrols->scrolledwindow_results), tgcontrols->table_results);

    i = 0;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("h<sub>c</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->hc = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->hc), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->hc, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("S"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->S = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->S), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->S, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("mN/nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("A<sub>p</sub>(h<sub>c</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->Aphc = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->Aphc), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->Aphc, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("nm<sup>2</sup>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    tgcontrols->extrapol = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->extrapol, 0, 3, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("SS<sub>res</sub>  fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->SSres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->SSres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->SSres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("R<sup>2</sup>  fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->R2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->R2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->R2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("R<sup>2</sup><sub>adj</sub>  fit "), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->R2adj = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->R2adj), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->R2adj, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("χ<sup>2</sup>  fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->chi2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->chi2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->chi2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    //	gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("R fit"), 0, 1, i, i+1, GTK_FILL, 0, 0, 0);
    tgcontrols->R = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->R), 0.0, 0.5);
    //	gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->R, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("H<sub>IT</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->Hit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->Hit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->Hit, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("MPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("E<sub>r</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->Er = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->Er), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->Er, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("E<sub>IT</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->Eit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->Eit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->Eit, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("Fitted range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    tgcontrols->range = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->range), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->range, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    tgcontrols->range_pct_Fmax = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(tgcontrols->range_pct_Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), tgcontrols->range_pct_Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(tgcontrols->table_results), label_new_with_markup_left("% F<sub>max</sub>"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);


    tgcontrols->button_unc = gtk_button_new_with_mnemonic("_Uncertainties");
    gtk_button_set_image(GTK_BUTTON(tgcontrols->button_unc), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(tgcontrols->button_unc, "clicked", G_CALLBACK(tangent_uncertainty), data);
    gtk_widget_set_sensitive(tgcontrols->button_unc, FALSE);

    tgcontrols->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(tgcontrols->button_save, "clicked", G_CALLBACK(tangent_save_data), args);

    //final  assembly

    gtk_box_pack_start(GTK_BOX(tgcontrols->ctrls), tgcontrols->frame_info, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(tgcontrols->ctrls), tgcontrols->frame_parameters, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(tgcontrols->ctrls), tgcontrols->button_fit, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(tgcontrols->ctrls), tgcontrols->frame_results, TRUE, TRUE, 0);
    //	gtk_box_pack_start (GTK_BOX (tgcontrols->ctrls), tgcontrols->button_unc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(tgcontrols->ctrls), tgcontrols->button_save, FALSE, FALSE, 0);

    /* Zooming */

    tgcontrols->vbox_zoom_buttons = gtk_vbox_new(TRUE, 0);
    tgcontrols->vbox_zoom_outer = gtk_vbox_new(FALSE, 0);

    tgcontrols->button_zoom_in = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(tgcontrols->button_zoom_in), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(tgcontrols->vbox_zoom_buttons), tgcontrols->button_zoom_in, TRUE, TRUE, 0);
    g_signal_connect_swapped(tgcontrols->button_zoom_in, "clicked", G_CALLBACK(tg_zoom_in), data);

    tgcontrols->button_zoom_out = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(tgcontrols->button_zoom_out), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(tgcontrols->vbox_zoom_buttons), tgcontrols->button_zoom_out, TRUE, TRUE, 0);
    g_signal_connect_swapped(tgcontrols->button_zoom_out, "clicked", G_CALLBACK(tg_zoom_out), data);

    tgcontrols->button_zoom_restore = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(tgcontrols->button_zoom_restore), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(tgcontrols->vbox_zoom_buttons), tgcontrols->button_zoom_restore, TRUE, TRUE, 0);
    g_signal_connect_swapped(tgcontrols->button_zoom_restore, "clicked", G_CALLBACK(tg_zoom_restore), data);

    gtk_box_pack_start(GTK_BOX(tgcontrols->vbox_zoom_outer), tgcontrols->vbox_zoom_buttons, TRUE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(tgcontrols->hbox_graph), tgcontrols->graph, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(tgcontrols->hbox_graph), tgcontrols->vbox_zoom_outer, FALSE, FALSE, 0);

    //zoom and unzoom buttons
    /* tgcontrols->button_zoom_in = gtk_button_new(); */
    /* gtk_button_set_image (GTK_BUTTON (tgcontrols->button_zoom_in), gtk_image_new_from_stock (GWY_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON)); */
    /* g_signal_connect_swapped(tgcontrols->button_zoom_in, "clicked", G_CALLBACK (zoom_in), data); */

    /* tgcontrols->button_zoom_out = gtk_button_new(); */
    /* gtk_button_set_image (GTK_BUTTON (tgcontrols->button_zoom_out), gtk_image_new_from_stock (GWY_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON)); */
    /* g_signal_connect_swapped(tgcontrols->button_zoom_out, "clicked", G_CALLBACK (zoom_out), data); */

    /* vbox = gtk_vbox_new(FALSE,FALSE); */
    /* gtk_box_pack_start (GTK_BOX(vbox), tgcontrols->graph, TRUE, TRUE, 0); */
    /* gtk_box_pack_start (GTK_BOX(vbox), tgcontrols->button_zoom_in, FALSE, FALSE, 0); */
    /* gtk_box_pack_start (GTK_BOX(vbox), tgcontrols->button_zoom_out, FALSE, FALSE, 0); */

    tgcontrols->tangent_gui = gtk_hbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(tgcontrols->tangent_gui), tgcontrols->ctrls, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(tgcontrols->tangent_gui), tgcontrols->hbox_graph, TRUE, TRUE, 0);

    gtk_widget_show(tgcontrols->tangent_gui);

    return tgcontrols->tangent_gui;
}

void tangent_recalculate(Data *data)
{
    recalculate_Eit(data->controls->tgcontrols, &(data->args->tgdata), &(data->args->instdata));
    tangent_propagate_uncertainties(data->controls->tgunccontrols, &(data->args->tgdata), &(data->args->tgunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

static void update_models_tg(TangentControls *tgcontrols)
{
    GwyGraphCurveModel *cmodel;

    tgcontrols->graph_model = gwy_graph_get_model(GWY_GRAPH(tgcontrols->graph));
    cmodel = gwy_graph_model_get_curve_by_description(tgcontrols->graph_model, "Unloading");

    if (cmodel) {
        tgcontrols->cmodelunload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(tgcontrols->graph_model, "Fit");

    if (cmodel) {
        tgcontrols->cmodelSfit = cmodel;
    }
}

static void tg_zoom_in(Data *data)
{
    TangentControls *tgcontrols;

    tgcontrols = data->controls->tgcontrols;

    rs_zoom_in(GWY_GRAPH(tgcontrols->graph), &(tgcontrols->zoomstack), FALSE);   /* do not rescale y */
    update_models_tg(tgcontrols);
}

static void tg_zoom_out(Data *data)
{
    TangentControls *tgcontrols;

    tgcontrols = data->controls->tgcontrols;

    rs_zoom_out(GWY_GRAPH(tgcontrols->graph), &(tgcontrols->zoomstack));
    update_models_tg(tgcontrols);
}

void tg_zoom_restore(Data *data)
{
    TangentControls *tgcontrols;

    tgcontrols = data->controls->tgcontrols;

    rs_zoom_restore(GWY_GRAPH(tgcontrols->graph), &(tgcontrols->zoomstack));
    update_models_tg(tgcontrols);
}

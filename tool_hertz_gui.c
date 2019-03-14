#include "tool_hertz_gui.h"
#include "tool_hertz.h"
#include "tool_hertz_unc_gui.h"

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

static void hertz_remove_all_fit_labels(HertzControls *hzcontrols);
static void hertz_remove_all_fit_curves(HertzControls *hzcontrols);

static void update_models_hz(HertzControls *hzcontrols);
static void hz_zoom_in(Data *data);
static void hz_zoom_out(Data *data);

void hertz_remove_fit_results_labels(Data *data)
{
    gwy_selection_clear(data->controls->hertzcontrols->selection);
    hertz_remove_all_fit_results(&(data->args->hertzdata));
    hertz_remove_all_fit_labels(data->controls->hertzcontrols);
    hertz_remove_all_fit_curves(data->controls->hertzcontrols);

    if (data->args->hertzunc.ERc != NULL) {
        hertz_unc_close_window(data);
    }
}

static void hertz_remove_all_fit_labels(HertzControls *hzcontrols)
{
    /* remove all labels that result from a fit */
    entry_clear(hzcontrols->from);
    entry_clear(hzcontrols->to);

    label_clear(hzcontrols->gamma);
    label_clear(hzcontrols->chi2);
    label_clear(hzcontrols->SSres);
    label_clear(hzcontrols->R2);
    label_clear(hzcontrols->R2adj);
    label_clear(hzcontrols->Erres);
    label_clear(hzcontrols->Eitres);
    label_clear(hzcontrols->range);
}

static void hertz_remove_all_fit_curves(HertzControls *hzcontrols)
{
    gwy_graph_model_remove_curve_by_description(hzcontrols->graph_model, "Fit");
}

static void recalculate_Eit(HertzControls *hertzcontrols, Hertzdata *hertzdata, const Instdata *instdata)
{
    if (hertzdata->has_fit) {
        hertzdata->Eit = calc_Eit(hertzdata->Er * 1e9, instdata->nu, instdata->Ei, instdata->nui) * 1e-9;
        label_set_gdouble_format(hertzcontrols->Eitres, hertzdata->Eit, MODUL);
    }
}

static void recalculate_Er(HertzControls *hertzcontrols, Hertzdata *hertzdata)
{
    hertzdata->Er = calc_Er_Hertz(hertzdata->reg.slope, hertzdata->radius * 1e9) * 1e6;
    label_set_gdouble_format(hertzcontrols->Erres, hertzdata->Er, MODUL);
}

static void recalculate_radius(HertzControls *hertzcontrols, Hertzdata *hertzdata)
{
    hertzdata->radius = calc_R_Hertz(hertzdata->reg.slope, hertzdata->Er) * 1e3;
    label_set_gdouble_format(hertzcontrols->radiusres, hertzdata->radius * 1e9, DEPTH);
}

static void radius_changed(Data *data)
{
    data->args->hertzdata.radius = gtk_adjustment_get_value(GTK_ADJUSTMENT(gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(data->controls->hertzcontrols->radiusinp)))) * 1e-9;

    if (data->args->hertzdata.has_fit) {
        recalculate_Er(data->controls->hertzcontrols, &(data->args->hertzdata));
        recalculate_Eit(data->controls->hertzcontrols, &(data->args->hertzdata), &(data->args->instdata));
    }

    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}

static void Er_changed(Data *data)
{
    data->args->hertzdata.Er = gtk_adjustment_get_value(GTK_ADJUSTMENT(gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(data->controls->hertzcontrols->Erinp))));

    if (data->args->hertzdata.has_fit) {
        recalculate_radius(data->controls->hertzcontrols, &(data->args->hertzdata));
    }

    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}

static void Eit_changed(Data *data)
{
    gdouble nu, nui, Ei;

    data->args->hertzdata.Eit = gtk_adjustment_get_value(GTK_ADJUSTMENT(gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(data->controls->hertzcontrols->Eitinp))));

    nu = data->args->instdata.nu;
    nui = data->args->instdata.nui;
    Ei = data->args->instdata.Ei;

    data->args->hertzdata.Er = calc_Er(data->args->hertzdata.Eit, nu, Ei * 1e-9, nui);

    if (data->args->hertzdata.has_fit) {
        recalculate_radius(data->controls->hertzcontrols, &(data->args->hertzdata));
    }

    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}

static void hertz_radio_callback(GtkWidget *widget, Data *data)
{
    HertzControls *hzcontrols;

    hzcontrols = data->controls->hertzcontrols;

    if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget))) {
        gtk_label_set_text(GTK_LABEL(data->controls->hertzcontrols->radiusres), NULL);
        gtk_label_set_text(GTK_LABEL(data->controls->hertzcontrols->Erres), NULL);
        gtk_label_set_text(GTK_LABEL(data->controls->hertzcontrols->Eitres), NULL);
        return;
    }

    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(hzcontrols->radio_R))) {
        data->args->hertzdata.mode = R_MODE;
        gtk_widget_set_sensitive(hzcontrols->radiusinp, TRUE);
        gtk_widget_set_sensitive(hzcontrols->Erinp, FALSE);
        gtk_widget_set_sensitive(hzcontrols->Eitinp, FALSE);
        radius_changed(data);
    }
    else if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(hzcontrols->radio_Er))) {
        data->args->hertzdata.mode = ER_MODE;
        gtk_widget_set_sensitive(hzcontrols->radiusinp, FALSE);
        gtk_widget_set_sensitive(hzcontrols->Erinp, TRUE);
        gtk_widget_set_sensitive(hzcontrols->Eitinp, FALSE);
        Er_changed(data);
    }
    else if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(hzcontrols->radio_Eit))) {
        data->args->hertzdata.mode = EIT_MODE;
        gtk_widget_set_sensitive(hzcontrols->radiusinp, FALSE);
        gtk_widget_set_sensitive(hzcontrols->Erinp, FALSE);
        gtk_widget_set_sensitive(hzcontrols->Eitinp, TRUE);
        Eit_changed(data);
    }
    else {
        g_printerr("You should never get here! Something weird is going on \n");
    }
}

void hertz_redraw(Data *data)
{
    // update loading curve and Fmax and hmax labels
    if (data->args->fddata.hload) {
        gwy_graph_curve_model_set_data(data->controls->hertzcontrols->cmodelload,
                                       data->args->fddata.hload->data, data->args->fddata.Fload->data,
                                       data->args->fddata.hload->res);

        label_set_gdouble_format(data->controls->hertzcontrols->hmax, data->args->fddata.hmax, DEPTH);
        label_set_gdouble_format(data->controls->hertzcontrols->Fmax, data->args->fddata.Fmax, FORCE);

        if (gwy_graph_model_get_n_curves(data->controls->hertzcontrols->graph_model) == 2) {
            gwy_graph_curve_model_set_data(data->controls->hertzcontrols->cmodelfit,
                                           data->args->hertzdata.xfit->data, data->args->hertzdata.yfit->data,  data->args->hertzdata.xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(data->controls->hertzcontrols->cmodelfit);
#endif
        }
    }
}

static void hz_range_from_changed(Data *data)
{
    HertzControls *fcontrols;
    Hertzdata *fdata;

    fcontrols = data->controls->hertzcontrols;
    fdata = &(data->args->hertzdata);

    range_from_changed(fcontrols->from, fcontrols->to, fcontrols->selection, &(fdata->from), &(fdata->from));
}

static void hz_range_to_changed(Data *data)
{
    HertzControls *tcontrols;
    Hertzdata *tdata;

    tcontrols = data->controls->hertzcontrols;
    tdata = &(data->args->hertzdata);

    range_to_changed(tcontrols->from, tcontrols->to, tcontrols->selection, &(tdata->to), &(tdata->to));
}

static void hertz_graph_selected(GwySelection *selection, gint i, Data *data)
{
    gdouble range[2];

    g_return_if_fail(i <= 0);

    //take data from selection and put them in the "from" and "to" values
    gwy_selection_get_object(selection, 0, range);

    data->args->hertzdata.from = MIN(range[0], range[1]);
    data->args->hertzdata.to = MAX(range[0], range[1]);

    //update labels
    entry_set_gdouble_format(data->controls->hertzcontrols->from, data->args->hertzdata.from, DEPTH);
    entry_set_gdouble_format(data->controls->hertzcontrols->to, data->args->hertzdata.to, DEPTH);
}

static void hertz_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (hertz_has_results(&(args->hertzdata))) {
	buffer = hertz_export_data(&(args->hertzdata), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

void hertz_fit_and_update(Data *data)
{
    LinReg reg;
    gchar str[30];
    gint nc, i;
    gint istart, iend, ndata;
    gboolean fit = TRUE;
    GwyDataLine *x, *y;

    Hertzdata *hzdata;
    FDdata *fddata;
    Instdata *instdata;
    HertzControls *hzcontrols;
    GdkColor color;

    gdk_color_parse("red", &color);

    //create pointers for easier reading
    hzdata = &(data->args->hertzdata);
    hzcontrols = data->controls->hertzcontrols;
    fddata = &(data->args->fddata);
    instdata = &(data->args->instdata);

    //if there was a previous fit, invalidate previous fitted curve
    if (hzdata->has_fit) {
        g_object_unref(hzdata->xfit);
        g_object_unref(hzdata->yfit);
    }

    // check range
    if (hzdata->from == hzdata->to) {
        fit = FALSE;
    }

    //find appropriate data, which should be fitted
    range_to_indices(hzdata->from, hzdata->to, fddata->hload, FALSE, &istart, &iend, &ndata);
    // if there are not enough data in the selected range, return fail
    ndata = MIN(ndata, fddata->hload->res - istart);

    if (ndata < 2) {
        fit = FALSE;
    }

    hzdata->nfitdata = ndata;

    if (fit) {
        //get data to be fitted
        x = gwy_data_line_part_extract(fddata->hload, istart, ndata);
        y = gwy_data_line_part_extract(fddata->Fload, istart, ndata);
        //fit the data and evaluate
        hertz_fit(x, y, hzdata, instdata);
        reg = hzdata->reg;

        //update labels
        entry_set_gdouble_format(hzcontrols->from, hzdata->from, DEPTH);
        entry_set_gdouble_format(hzcontrols->to, hzdata->to, DEPTH);
        label_set_gdouble_format(hzcontrols->chi2, reg.chi, NUMBER);

        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", hzdata->from, hzdata->to);
        gtk_label_set_text(GTK_LABEL(hzcontrols->range), str);

        label_set_gdouble_format(hzcontrols->gamma, hzdata->reg.slope, "%.4g");
        label_set_gdouble_format(hzcontrols->SSres, hzdata->reg.SSres, "%.4g");
        label_set_gdouble_format(hzcontrols->R2, hzdata->reg.R2, "%.4g");
        label_set_gdouble_format(hzcontrols->R2adj, hzdata->reg.R2adj, "%.4g");
        label_set_gdouble_format(hzcontrols->chi2, hzdata->reg.chi2, "%.4g");

        if (hzdata->mode == R_MODE) {
            label_set_gdouble_format(hzcontrols->Eitres, hzdata->Eit, MODUL);
            label_set_gdouble_format(hzcontrols->Erres, hzdata->Er, MODUL);
        }
        else {
            label_set_gdouble_format(hzcontrols->radiusres, hzdata->radius * 1e9, DEPTH);
        }

        // create fit curve data
        hzdata->xfit = gwy_data_line_duplicate(x);
        hzdata->yfit = gwy_data_line_part_extract(fddata->hload, 0, ndata);

        for (i = 0; i < hzdata->xfit->res; i++) {
            hzdata->yfit->data[i] =  reg.slope * pow(hzdata->xfit->data[i], 1.5) + reg.intercept;
        }

        //clean up
        g_object_unref(x);
        g_object_unref(y);

        // update fitted curve, if necessary create it.
        nc = gwy_graph_model_get_n_curves(hzcontrols->graph_model);

        if (nc == 2) {
            gwy_graph_curve_model_set_data(hzcontrols->cmodelfit, hzdata->xfit->data, hzdata->yfit->data, hzdata->xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(hzcontrols->cmodelfit);
#endif
        }
        else {
            gwy_graph_curve_model_set_data(hzcontrols->cmodelfit, hzdata->xfit->data, hzdata->yfit->data, hzdata->xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(hzcontrols->cmodelfit);
#endif
            gwy_graph_model_add_curve(hzcontrols->graph_model, hzcontrols->cmodelfit);
            g_object_set(hzcontrols->cmodelfit,
                         "mode", GWY_GRAPH_CURVE_LINE,
                         "description", "Fit",
                         "color", CURVE_COLOR_FIT,
                         "line-style", GDK_LINE_SOLID,
                         "line-width", 3,
                         NULL);
        }
    }
    else {
        nc = gwy_graph_model_get_n_curves(hzcontrols->graph_model);
        hertz_remove_all_fit_labels(hzcontrols);
        hertz_remove_all_fit_results(hzdata);
        hzdata->has_fit = FALSE;

        if (nc == 2) {
            gwy_graph_model_remove_curve_by_description(hzcontrols->graph_model, "Fit");
        }
    }

    // update previous uncertainties if existing, set button
    if (hzdata->has_fit) {
        hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
        gtk_widget_set_sensitive(hzcontrols->button_unc, TRUE);
    }
    else {
        gtk_widget_set_sensitive(hzcontrols->button_unc, FALSE);
    }
}

GtkWidget *tool_hertz_create(Data *data)
{
    GtkWidget *label;
    GtkObject *adj;
    GwyGraphArea *area;
    Args *args;
    HertzControls *hzcontrols;
    gint i;

    args = data->args;
    hzcontrols = data->controls->hertzcontrols;

    // Buttons

    hzcontrols->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(hzcontrols->button_save, "clicked", G_CALLBACK(hertz_save_data), args);

    hzcontrols->button_unc = gtk_button_new_with_mnemonic("_Uncertainties");
    gtk_button_set_image(GTK_BUTTON(hzcontrols->button_unc), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(hzcontrols->button_unc, "clicked", G_CALLBACK(hertz_uncertainty), data);
    gtk_widget_set_sensitive(hzcontrols->button_unc, FALSE);

    // GraphArea, two curves will be displayed: loading curve (F-d) and the resulting fit
    hzcontrols->graph_model = gwy_graph_model_new();
    g_object_set(hzcontrols->graph_model, "title", "Load", NULL);
    gwy_graph_model_set_axis_label(hzcontrols->graph_model, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(hzcontrols->graph_model, GTK_POS_LEFT, "F / mN");

    hzcontrols->zoomstack = NULL;
    hzcontrols->zoomstack = g_slist_prepend(hzcontrols->zoomstack, hzcontrols->graph_model);

    hzcontrols->hbox_graph = gtk_hbox_new(FALSE, CTRLS_SPACING);

    hzcontrols->graph = gwy_graph_new(hzcontrols->graph_model);
    /* gtk_widget_set_size_request(GTK_WIDGET(hzcontrols->graph), 600, 400); */

    gwy_graph_enable_user_input(GWY_GRAPH(hzcontrols->graph), FALSE);
    gwy_graph_set_status(GWY_GRAPH(hzcontrols->graph), GWY_GRAPH_STATUS_XSEL);

    area = GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(hzcontrols->graph)));
    hzcontrols->selection = gwy_graph_area_get_selection(area, GWY_GRAPH_STATUS_XSEL);
    gwy_selection_set_max_objects(hzcontrols->selection, 1);
    g_signal_connect(hzcontrols->selection, "changed", G_CALLBACK(hertz_graph_selected), data);

    hzcontrols->cmodelload = gwy_graph_curve_model_new();
    hzcontrols->cmodelfit = gwy_graph_curve_model_new();

    g_object_set(hzcontrols->cmodelload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Loading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_LOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);
    gwy_graph_model_add_curve(hzcontrols->graph_model, hzcontrols->cmodelload);


    /* Ctrls */

    hzcontrols->ctrls = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* Info section */

    hzcontrols->frame_info = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(hzcontrols->frame_info), GTK_WIDGET(gwy_label_new_header("Info")));
    hzcontrols->table_info = gtk_table_new(2, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(hzcontrols->table_info), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(hzcontrols->table_info), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(hzcontrols->table_info), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(hzcontrols->frame_info), hzcontrols->table_info);

    i = 0;

    /* h_max */
    gtk_table_attach(GTK_TABLE(hzcontrols->table_info), label_new_with_markup_left("h<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    hzcontrols->hmax = gtk_label_new(NULL);
    label_set_gdouble_format(hzcontrols->hmax, args->fddata.hmax, DEPTH);
    gtk_misc_set_alignment(GTK_MISC(hzcontrols->hmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_info), hzcontrols->hmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_info), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /* F_max */
    gtk_table_attach(GTK_TABLE(hzcontrols->table_info), label_new_with_markup_left("F<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    hzcontrols->Fmax = gtk_label_new(NULL);
    label_set_gdouble_format(hzcontrols->Fmax, args->fddata.Fmax, FORCE);
    gtk_misc_set_alignment(GTK_MISC(hzcontrols->Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_info), hzcontrols->Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_info), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    /* Parameters section */

    hzcontrols->frame_parameters = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(hzcontrols->frame_parameters), GTK_WIDGET(gwy_label_new_header("Parameters")));
    hzcontrols->table_parameters = gtk_table_new(4, 4, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(hzcontrols->table_parameters), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(hzcontrols->table_parameters), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(hzcontrols->table_parameters), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(hzcontrols->frame_parameters), hzcontrols->table_parameters);

    i = 0;

    hzcontrols->radio_R   = gtk_radio_button_new_with_label(NULL, "Radius");
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), hzcontrols->radio_R, 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(G_OBJECT(hzcontrols->radio_R), "toggled", G_CALLBACK(hertz_radio_callback), data);
    /*hzcontrols->radiusinp = gtk_adjustment_new(args->hertzdata.radius*1e9, 1e-9, 1e6, 1, 10,0);
    gtk_widget_set_sensitive(hzcontrols->radiusinp, TRUE);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(hzcontrols->radiusinp), 0.1, 1);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), spin, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    */
    adj = gtk_adjustment_new(args->hertzdata.radius * 1e9, 1e-9, 1e6, 1, 10, 0);
    hzcontrols->radiusinp = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 0.1, 1);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(hzcontrols->radiusinp), TRUE);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), hzcontrols->radiusinp, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_widget_set_sensitive(hzcontrols->radiusinp, TRUE);
    data->args->hertzdata.mode = R_MODE;
    g_signal_connect_swapped(hzcontrols->radiusinp, "value-changed", G_CALLBACK(radius_changed), data);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    hzcontrols->radio_Er  = gtk_radio_button_new_from_widget(GTK_RADIO_BUTTON(hzcontrols->radio_R));
    label = label_new_with_markup("E<sub>r</sub>");
    gtk_container_add(GTK_CONTAINER(hzcontrols->radio_Er), label);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), hzcontrols->radio_Er, 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(G_OBJECT(hzcontrols->radio_Er), "toggled", G_CALLBACK(hertz_radio_callback), data);
    /*
    hzcontrols->Erinp = gtk_adjustment_new(args->hertzdata.Er*1e-9, 1e-9, 1e6, 1, 10,0);
    gtk_widget_set_sensitive(hzcontrols->Erinp, FALSE);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(hzcontrols->Erinp), 0.1, 1);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), spin, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(hzcontrols->Erinp, "value-changed", G_CALLBACK(Er_changed), data);
    */
    adj = gtk_adjustment_new(10, 1e-3, 1e6, 1, 10, 0);
    hzcontrols->Erinp = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 0.001, 3);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(hzcontrols->Erinp), TRUE);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), hzcontrols->Erinp, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_widget_set_sensitive(hzcontrols->Erinp, FALSE);
    g_signal_connect_swapped(hzcontrols->Erinp, "value-changed", G_CALLBACK(Er_changed), data);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    hzcontrols->radio_Eit  = gtk_radio_button_new_from_widget(GTK_RADIO_BUTTON(hzcontrols->radio_R));
    label = label_new_with_markup("E<sub>IT</sub>");
    gtk_container_add(GTK_CONTAINER(hzcontrols->radio_Eit), label);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), hzcontrols->radio_Eit, 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(G_OBJECT(hzcontrols->radio_Eit), "toggled", G_CALLBACK(hertz_radio_callback), data);
    /*
    hzcontrols->Eitinp = gtk_adjustment_new(args->hertzdata.Eit*1e-9, 1e-9, 1e6, 1, 10,0);
    gtk_widget_set_sensitive(GTK_WIDGET(hzcontrols->Eitinp), FALSE);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(hzcontrols->Eitinp), 0.1, 1);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), spin, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(hzcontrols->Eitinp, "value-changed", G_CALLBACK(Eit_changed), data);
    */
    adj = gtk_adjustment_new(10, 1e-3, 1e6, 1, 10, 0);
    hzcontrols->Eitinp = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 0.001, 3);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(hzcontrols->Eitinp), TRUE);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), hzcontrols->Eitinp, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_widget_set_sensitive(hzcontrols->Eitinp, FALSE);
    g_signal_connect_swapped(hzcontrols->Eitinp, "value-changed", G_CALLBACK(Eit_changed), data);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), label_new_with_markup_left("Fit range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzcontrols->from = gtk_entry_new();
    g_object_set_data(G_OBJECT(hzcontrols->from), "id", (gpointer)"from");
    gtk_entry_set_width_chars(GTK_ENTRY(hzcontrols->from), 8);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), hzcontrols->from, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(hzcontrols->from, TRUE);
    g_signal_connect_swapped(hzcontrols->from, "activate", G_CALLBACK(hz_range_from_changed), data);
    /*	g_signal_connect_swapped(hzcontrols->from, "focus-out-event", G_CALLBACK(range_from_changed),data); */

    hzcontrols->to = gtk_entry_new();
    g_object_set_data(G_OBJECT(hzcontrols->to), "id", (gpointer)"to");
    gtk_entry_set_width_chars(GTK_ENTRY(hzcontrols->to), 8);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), hzcontrols->to, 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(hzcontrols->to, TRUE);
    g_signal_connect_swapped(hzcontrols->to, "activate", G_CALLBACK(hz_range_to_changed), data);
    /*	g_signal_connect_swapped(hzcontrols->to, "focus-out-event", G_CALLBACK(range_to_changed),data);*/

    gtk_table_attach(GTK_TABLE(hzcontrols->table_parameters), label_new_with_markup_left("nm"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    /* Fit button */
    hzcontrols->button_fit = gtk_button_new_with_mnemonic(gwy_sgettext("verb|_Fit"));
    gtk_button_set_image(GTK_BUTTON(hzcontrols->button_fit), gtk_image_new_from_stock(GTK_STOCK_EXECUTE, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(hzcontrols->button_fit, "clicked", G_CALLBACK(hertz_fit_and_update), data);

    /* Results section */

    hzcontrols->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(hzcontrols->frame_results), GTK_WIDGET(gwy_label_new_header("Results")));

    hzcontrols->scrolledwindow_results = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(hzcontrols->scrolledwindow_results), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(hzcontrols->frame_results), hzcontrols->scrolledwindow_results);

    hzcontrols->table_results = gtk_table_new(5, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(hzcontrols->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(hzcontrols->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(hzcontrols->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(hzcontrols->scrolledwindow_results), hzcontrols->table_results);

    i = 0;

    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("γ"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzcontrols->gamma = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzcontrols->gamma), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), hzcontrols->gamma, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /*
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("Tip radius"), 0, 1, i, i+1, GTK_FILL, 0, 0, 0);
    hzcontrols->radius = gtk_adjustment_new(args->hertzdata.radius*1e9, 1e-9, 1e6, 1, 10,0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(hzcontrols->radius), 0.1, 1);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), spin, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(hzcontrols->radius, "value-changed", G_CALLBACK(radius_changed), data);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i+1, GTK_FILL, 0, 0, 0);
    i++;
    */

    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("SS<sub>res</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzcontrols->SSres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzcontrols->SSres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), hzcontrols->SSres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("R<sup>2</sup>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzcontrols->R2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzcontrols->R2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), hzcontrols->R2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("R<sup>2</sup><sub>adj</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzcontrols->R2adj = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzcontrols->R2adj), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), hzcontrols->R2adj, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("χ<sup>2</sup> fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzcontrols->chi2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzcontrols->chi2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), hzcontrols->chi2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("Tip radius"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzcontrols->radiusres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzcontrols->radiusres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), hzcontrols->radiusres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("E<sub>r</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzcontrols->Erres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzcontrols->Erres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), hzcontrols->Erres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("E<sub>IT</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzcontrols->Eitres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzcontrols->Eitres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), hzcontrols->Eitres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("Fitted range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzcontrols->range = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzcontrols->range), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), hzcontrols->range, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);


    //final assembly

    gtk_box_pack_start(GTK_BOX(hzcontrols->ctrls), hzcontrols->frame_info, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hzcontrols->ctrls), hzcontrols->frame_parameters, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hzcontrols->ctrls), hzcontrols->button_fit, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hzcontrols->ctrls), hzcontrols->frame_results, TRUE, TRUE, 0);
    //	gtk_box_pack_start (GTK_BOX (hzcontrols->ctrls), hzcontrols->button_unc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hzcontrols->ctrls), hzcontrols->button_save, FALSE, FALSE, 0);


    /* Zooming */

    hzcontrols->vbox_zoom_buttons = gtk_vbox_new(TRUE, 0);
    hzcontrols->vbox_zoom_outer = gtk_vbox_new(FALSE, 0);

    hzcontrols->button_zoom_in = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(hzcontrols->button_zoom_in), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(hzcontrols->vbox_zoom_buttons), hzcontrols->button_zoom_in, TRUE, TRUE, 0);
    g_signal_connect_swapped(hzcontrols->button_zoom_in, "clicked", G_CALLBACK(hz_zoom_in), data);

    hzcontrols->button_zoom_out = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(hzcontrols->button_zoom_out), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(hzcontrols->vbox_zoom_buttons), hzcontrols->button_zoom_out, TRUE, TRUE, 0);
    g_signal_connect_swapped(hzcontrols->button_zoom_out, "clicked", G_CALLBACK(hz_zoom_out), data);

    hzcontrols->button_zoom_restore = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(hzcontrols->button_zoom_restore), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(hzcontrols->vbox_zoom_buttons), hzcontrols->button_zoom_restore, TRUE, TRUE, 0);
    g_signal_connect_swapped(hzcontrols->button_zoom_restore, "clicked", G_CALLBACK(hz_zoom_restore), data);

    gtk_box_pack_start(GTK_BOX(hzcontrols->vbox_zoom_outer), hzcontrols->vbox_zoom_buttons, TRUE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(hzcontrols->hbox_graph), hzcontrols->graph, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(hzcontrols->hbox_graph), hzcontrols->vbox_zoom_outer, FALSE, FALSE, 0);

    //zoom and unzoom buttons
    /* hzcontrols->button_zoom_in = gtk_button_new(); */
    /* gtk_button_set_image (GTK_BUTTON (hzcontrols->button_zoom_in), gtk_image_new_from_stock (GWY_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON)); */
    /* g_signal_connect_swapped(hzcontrols->button_zoom_in, "clicked", G_CALLBACK (zoom_in), data); */

    /* hzcontrols->button_zoom_out = gtk_button_new(); */
    /* gtk_button_set_image (GTK_BUTTON (hzcontrols->button_zoom_out), gtk_image_new_from_stock (GWY_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON)); */
    /* g_signal_connect_swapped(hzcontrols->button_zoom_out, "clicked", G_CALLBACK (zoom_out), data); */

    /* vbox = gtk_vbox_new(FALSE,FALSE); */
    /* gtk_box_pack_start (GTK_BOX(vbox), hzcontrols->graph, TRUE, TRUE, 0); */
    /* gtk_box_pack_start (GTK_BOX(vbox), hzcontrols->button_zoom_in, FALSE, FALSE, 0); */
    /* gtk_box_pack_start (GTK_BOX(vbox), hzcontrols->button_zoom_out, FALSE, FALSE, 0); */

    hzcontrols->hertz_gui = gtk_hbox_new(FALSE, CTRLS_SPACING);

    gtk_box_pack_start(GTK_BOX(hzcontrols->hertz_gui), hzcontrols->ctrls, FALSE, FALSE, 0);   /* expand, fill, padding */
    gtk_box_pack_start(GTK_BOX(hzcontrols->hertz_gui), hzcontrols->hbox_graph, TRUE, TRUE, 0);

    gtk_widget_show(hzcontrols->hertz_gui);

    return hzcontrols->hertz_gui;
}

void hertz_recalculate(Data *data)
{
    if (data->args->hertzdata.mode == R_MODE) {
        recalculate_Eit(data->controls->hertzcontrols, &(data->args->hertzdata), &(data->args->instdata));
    }
    else if (data->args->hertzdata.mode == EIT_MODE) {
	data->args->hertzdata.Er = calc_Er(data->args->hertzdata.Eit, data->args->instdata.nu, data->args->instdata.Ei * 1e-9, data->args->instdata.nui);
        recalculate_radius(data->controls->hertzcontrols, &(data->args->hertzdata));
    }

    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}

static void update_models_hz(HertzControls *hzcontrols)
{
    GwyGraphCurveModel *cmodel;

    hzcontrols->graph_model = gwy_graph_get_model(GWY_GRAPH(hzcontrols->graph));
    cmodel = gwy_graph_model_get_curve_by_description(hzcontrols->graph_model, "Loading");

    if (cmodel) {
        hzcontrols->cmodelload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(hzcontrols->graph_model, "Fit");

    if (cmodel) {
        hzcontrols->cmodelfit = cmodel;
    }
}

static void hz_zoom_in(Data *data)
{
    HertzControls *hzcontrols;

    hzcontrols = data->controls->hertzcontrols;

    rs_zoom_in(GWY_GRAPH(hzcontrols->graph), &(hzcontrols->zoomstack), FALSE);   /* do not rescale y */
    update_models_hz(hzcontrols);
}

static void hz_zoom_out(Data *data)
{
    HertzControls *hzcontrols;

    hzcontrols = data->controls->hertzcontrols;

    rs_zoom_out(GWY_GRAPH(hzcontrols->graph), &(hzcontrols->zoomstack));
    update_models_hz(hzcontrols);
}

void hz_zoom_restore(Data *data)
{
    HertzControls *hzcontrols;

    hzcontrols = data->controls->hertzcontrols;

    rs_zoom_restore(GWY_GRAPH(hzcontrols->graph), &(hzcontrols->zoomstack));
    update_models_hz(hzcontrols);
}

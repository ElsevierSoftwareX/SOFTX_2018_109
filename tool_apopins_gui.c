#include "tool_apopins_gui.h"
#include "tool_apopins.h"

#include "controls.h"
#include "datatypes.h"
#include "gui-utils.h"
#include "niget-common.h"

#include <math.h>

#include <gtk/gtk.h>

#include <app/file.h>
#include <libprocess/gwyprocess.h>


#define NDER 3
#define NLOAD 2


static void apopin_remove_labels(ApopinControls *apopcontrols);
static void apopin_remove_curves(ApopinControls *apopcontrols);

static void update_models_apop_h(ApopinControls *apopcontrols);
static void apop_zoom_in_linked_h(Data *data);

static void update_models_apop_der(ApopinControls *apopcontrols);
static void apop_zoom_in_linked_der(Data *data);

static void apop_zoom_out_linked(Data *data);

static void apopin_calc_and_update(Data *data);

static void hpop_changed(Data *data)
{
    data->args->apopdata.hpop = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->apopcontrols->hpop));
    apopin_calc_and_update(data);
}

static void wpop_changed(Data *data)
{
    data->args->apopdata.wpop = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->apopcontrols->wpop));
    apopin_calc_and_update(data);
}

static void thresh_changed(Data *data)
{
    data->args->apopdata.thresh = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->apopcontrols->thresh));
    gwy_data_line_resample(data->args->apopdata.dth, data->args->fddata.hload->res, GWY_INTERPOLATION_NONE);
    gwy_data_line_clear(data->args->apopdata.dth);
    gwy_data_line_add(data->args->apopdata.dth, data->args->apopdata.thresh);

    apopin_calc_and_update(data);
}

static void thresh2_changed(Data *data)
{
    data->args->apopdata.thresh2 = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->apopcontrols->thresh2));
    gwy_data_line_resample(data->args->apopdata.dth2, data->args->fddata.hload->res, GWY_INTERPOLATION_NONE);
    gwy_data_line_clear(data->args->apopdata.dth2);
    gwy_data_line_add(data->args->apopdata.dth2, data->args->apopdata.thresh2);
    apopin_calc_and_update(data);
}

static void nmove_changed(Data *data)
{
    data->args->apopdata.nmove = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->apopcontrols->nmove));
    apopin_calc_and_update(data);
}

void apopin_remove_fit_results_labels(Data *data)
{
    apopin_remove_results(&(data->args->apopdata));
    apopin_remove_labels(data->controls->apopcontrols);
    apopin_remove_curves(data->controls->apopcontrols);
}

static void apopin_remove_labels(ApopinControls *apopcontrols)
{
    /* remove table of popins */
    remove_table_content(apopcontrols->res);

    /* remove number of popins */
    label_clear(apopcontrols->npopin);
}

static void apopin_remove_curves(ApopinControls *apopcontrols)
{
    gwy_graph_model_remove_all_curves(apopcontrols->graph_model);
    gwy_graph_model_remove_all_curves(apopcontrols->graph_model_der);
    gwy_graph_model_add_curve(apopcontrols->graph_model, apopcontrols->cmodelload);
}

static void reset_default(Data *data)
{
    data->args->apopdata.nmove = 1;
    data->args->apopdata.thresh = THRESH_DEFAULT;
    data->args->apopdata.thresh2 = THRESH2_DEFAULT;
    data->args->apopdata.wpop = WPOP_DEFAULT;
    data->args->apopdata.hpop = HPOP_DEFAULT;

    gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->apopcontrols->nmove), data->args->apopdata.nmove);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->apopcontrols->thresh), data->args->apopdata.thresh);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->apopcontrols->thresh2), data->args->apopdata.thresh2);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->apopcontrols->wpop), data->args->apopdata.wpop);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->apopcontrols->hpop), data->args->apopdata.hpop);

    apopin_remove_results(&(data->args->apopdata));
    apopin_remove_labels(data->controls->apopcontrols);
    apopin_remove_curves(data->controls->apopcontrols);
    apopin_redraw(data);
}

void apopin_redraw(Data *data)
{
    gint nc;
    Apopindata *apopdata;
    ApopinControls *apopctrls;

    apopdata = &(data->args->apopdata);
    apopctrls = data->controls->apopcontrols;

    apopin_create_aux_datalines(apopdata, &data->args->fddata);

    gwy_graph_curve_model_set_data(apopctrls->cmodelload,
                                   apopdata->t->data, data->args->fddata.hload->data,
                                   data->args->fddata.hload->res);

    nc = gwy_graph_model_get_n_curves(apopctrls->graph_model_der);

    if (nc == NDER) {
        gwy_graph_curve_model_set_data(apopctrls->cmodelder,
                                       apopdata->t->data, apopdata->dh->data,
                                       apopdata->dh->res);
        gwy_graph_curve_model_set_data(apopctrls->cmodelthresh,
                                       apopdata->t->data, apopdata->dth->data,
                                       apopdata->dth->res);
        gwy_graph_curve_model_set_data(apopctrls->cmodelthresh2,
                                       apopdata->t->data, apopdata->dth2->data,
                                       apopdata->dth2->res);
    }
}

void apopin_calc_and_update(Data *data)
{
    gint i, j;

    gchar str[30];
    GtkWidget *label;
    gint nc;
    GwyDataLine *hpopin, *tpopin, *dhpopin;

    Apopindata *apopdata;
    ApopinControls *apopctrls;
    GwyGraphCurveModel **cmodelpopin;
    GwyGraphCurveModel **cmodeldpopin;

    //local copies for simpler reading
    apopdata = &(data->args->apopdata);
    apopctrls = data->controls->apopcontrols;

    apopin_calc(apopdata, &data->args->fddata);

    apop_zoom_restore_linked(data);

    // if no sliding average data displayed, add them. update them in any case.
    nc = gwy_graph_model_get_n_curves(apopctrls->graph_model);

    if (nc == 1) {
        gwy_graph_model_add_curve(apopctrls->graph_model, apopctrls->cmodelhavg);
    }

    gwy_graph_curve_model_set_data(apopctrls->cmodelhavg, apopdata->t->data, apopdata->havg->data, apopdata->havg->res);

    // if popin data displayed, remove them
    nc = gwy_graph_model_get_n_curves(apopctrls->graph_model);

    for (i = 0; i < nc - NLOAD; i++) {
        g_snprintf(str, sizeof(str), "Pop-in #%d", i + 1);
        gwy_graph_model_remove_curve_by_description(apopctrls->graph_model, str);
    }

    nc = gwy_graph_model_get_n_curves(apopctrls->graph_model_der);

    for (i = 0; i < nc - NDER; i++) {
        g_snprintf(str, sizeof(str), "Pop-in #%d", i + 1);
        gwy_graph_model_remove_curve_by_description(apopctrls->graph_model_der, str);
    }

    // if no derivative data displayed, add curve; update them in anycase.
    if (gwy_graph_model_get_n_curves(apopctrls->graph_model_der) == 0) {
        gwy_graph_model_add_curve(apopctrls->graph_model_der, apopctrls->cmodelder);
        gwy_graph_model_add_curve(apopctrls->graph_model_der, apopctrls->cmodelthresh);
        gwy_graph_model_add_curve(apopctrls->graph_model_der, apopctrls->cmodelthresh2);
    }

    gwy_graph_curve_model_set_data(apopctrls->cmodelder, apopdata->t->data, apopdata->dh->data, apopdata->t->res);
    gwy_graph_curve_model_set_data(apopctrls->cmodelthresh, apopdata->t->data, apopdata->dth->data, apopdata->t->res);
    gwy_graph_curve_model_set_data(apopctrls->cmodelthresh2, apopdata->t->data, apopdata->dth2->data, apopdata->t->res);

    //remove all previous results, that were displayed in the table
    remove_table_content(apopctrls->res);
    gtk_table_resize(GTK_TABLE(apopctrls->res), 1, 4);

    // update npopins
    label_set_gint_format(apopctrls->npopin, apopdata->npopin, "%d");

    //create new results table
    i = 0;

    gtk_table_resize(GTK_TABLE(apopctrls->res), 1 + apopdata->npopin, 4);

    label = gtk_label_new("Load / mN");
    gtk_table_attach(GTK_TABLE(apopctrls->res), label, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    label = gtk_label_new("Depth / nm");
    gtk_table_attach(GTK_TABLE(apopctrls->res), label, 2, 3, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    label = gtk_label_new("Depth difference / nm");
    gtk_table_attach(GTK_TABLE(apopctrls->res), label, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    i++;

    for (j = 0; j < apopdata->npopin; j++) {
        g_snprintf(str, sizeof(str), "#%d:", j + 1);
        label = gtk_label_new(str);
        gtk_table_attach(GTK_TABLE(apopctrls->res), label, 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        g_snprintf(str, sizeof(str), FORCE, apopdata->Fpopin[j]);
        label = gtk_label_new(str);
        gtk_table_attach(GTK_TABLE(apopctrls->res), label, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        g_snprintf(str, sizeof(str), DEPTH, apopdata->hpopin[j]);
        label = gtk_label_new(str);
        gtk_table_attach(GTK_TABLE(apopctrls->res), label, 2, 3, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        g_snprintf(str, sizeof(str), DEPTH, apopdata->dhpopin[j]);
        label = gtk_label_new(str);
        gtk_table_attach(GTK_TABLE(apopctrls->res), label, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        i++;
    }

    gtk_widget_show_all(apopctrls->res);

    //create curve for every pop-in and display
    cmodelpopin = (GwyGraphCurveModel **)g_malloc(apopdata->npopin * sizeof(GwyGraphCurveModel *));

    for (i = 0; i < apopdata->npopin; i++) {
        cmodelpopin[i] = gwy_graph_curve_model_new();
        gwy_graph_model_add_curve(apopctrls->graph_model, cmodelpopin[i]);
        g_snprintf(str, sizeof(str), "Pop-in #%d", i + 1);
        g_object_set(cmodelpopin[i],
                     "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                     "description", str,
                     "point-size", 0,
                     "point-type", GWY_GRAPH_POINT_CIRCLE,
                     "color", gwy_graph_get_preset_color((guint)(i + fmax(NLOAD, NDER) + 1)),
                     "line-style", GDK_LINE_SOLID,
                     "line-width", 2,
                     NULL);
        tpopin = gwy_data_line_part_extract(apopdata->t, apopdata->ileft[i], apopdata->iright[i] - apopdata->ileft[i] + 1);
        hpopin = gwy_data_line_part_extract(data->args->fddata.hload, apopdata->ileft[i], apopdata->iright[i] - apopdata->ileft[i] + 1);
        gwy_graph_curve_model_set_data(cmodelpopin[i], tpopin->data, hpopin->data, tpopin->res);
    }

    g_free(cmodelpopin);

    cmodeldpopin = (GwyGraphCurveModel **)g_malloc(apopdata->npopin * sizeof(GwyGraphCurveModel *));

    for (i = 0; i < apopdata->npopin; i++) {
        cmodeldpopin[i] = gwy_graph_curve_model_new();
        gwy_graph_model_add_curve(apopctrls->graph_model_der, cmodeldpopin[i]);
        g_snprintf(str, sizeof(str), "Pop-in #%d", i + 1);
        g_object_set(cmodeldpopin[i],
                     "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                     "description", str,
                     "point-size", 0,
                     "point-type", GWY_GRAPH_POINT_CIRCLE,
                     "color", gwy_graph_get_preset_color((guint)(i + fmax(NLOAD, NDER) + 1)),
                     "line-style", GDK_LINE_SOLID,
                     "line-width", 2,
                     NULL);
        tpopin = gwy_data_line_part_extract(apopdata->t, apopdata->ileft[i], apopdata->iright[i] - apopdata->ileft[i] + 1);
        dhpopin = gwy_data_line_part_extract(apopdata->dh, apopdata->ileft[i], apopdata->iright[i] - apopdata->ileft[i] + 1);
        gwy_graph_curve_model_set_data(cmodeldpopin[i], tpopin->data, dhpopin->data, tpopin->res);
    }

    g_free(cmodeldpopin);
}

static void apopin_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (apopin_has_results(&(args->apopdata))) {
	buffer = apopin_export_data(&(args->apopdata), &(args->fddata), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

// create gui
GtkWidget *tool_apopins_create(Data *data)
{
    GtkWidget *spin;
    Args *args;
    ApopinControls *apopcontrols;
    gint i;

    args = data->args;
    apopcontrols = data->controls->apopcontrols;

    // GraphAreas: one for the dependence depth-(pseudo)time (== index), one for the derivative of the smoothed depth
    apopcontrols->graph_model = gwy_graph_model_new();
    g_object_set(apopcontrols->graph_model, "title", "Depth", NULL);
    gwy_graph_model_set_axis_label(apopcontrols->graph_model, GTK_POS_BOTTOM, "i");
    gwy_graph_model_set_axis_label(apopcontrols->graph_model, GTK_POS_LEFT, "h / nm");

    apopcontrols->zoomstack_h = NULL;
    apopcontrols->zoomstack_h = g_slist_prepend(apopcontrols->zoomstack_h, apopcontrols->graph_model);

    apopcontrols->graph_h = gwy_graph_new(apopcontrols->graph_model);
    gtk_widget_set_size_request(GTK_WIDGET(apopcontrols->graph_h), 600, 250);

    gwy_graph_enable_user_input(GWY_GRAPH(apopcontrols->graph_h), FALSE);
    gwy_graph_set_status(GWY_GRAPH(apopcontrols->graph_h), GWY_GRAPH_STATUS_XSEL);

    apopcontrols->cmodelload = gwy_graph_curve_model_new();
    apopcontrols->cmodelhavg = gwy_graph_curve_model_new();

    g_object_set(apopcontrols->cmodelload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Depth",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_LOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);
    gwy_graph_model_add_curve(apopcontrols->graph_model, apopcontrols->cmodelload);

    g_object_set(apopcontrols->cmodelhavg,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Smoothed depth",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_SMOOTH_LOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    apopcontrols->graph_model_der = gwy_graph_model_new();
    g_object_set(apopcontrols->graph_model_der, "title", "Derivative", NULL);
    gwy_graph_model_set_axis_label(apopcontrols->graph_model_der, GTK_POS_BOTTOM, "i");
    gwy_graph_model_set_axis_label(apopcontrols->graph_model_der, GTK_POS_LEFT, "dh/di / nm");

    apopcontrols->zoomstack_der = NULL;
    apopcontrols->zoomstack_der = g_slist_prepend(apopcontrols->zoomstack_der, apopcontrols->graph_model_der);

    apopcontrols->graph_der = gwy_graph_new(apopcontrols->graph_model_der);
    gtk_widget_set_size_request(GTK_WIDGET(apopcontrols->graph_der), 600, 150);

    gwy_graph_enable_user_input(GWY_GRAPH(apopcontrols->graph_der), FALSE);
    gwy_graph_set_status(GWY_GRAPH(apopcontrols->graph_der), GWY_GRAPH_STATUS_XSEL);

    apopcontrols->cmodelder = gwy_graph_curve_model_new();
    apopcontrols->cmodelthresh = gwy_graph_curve_model_new();
    apopcontrols->cmodelthresh2 = gwy_graph_curve_model_new();

    g_object_set(apopcontrols->cmodelder,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Derivative",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_DER,
                 "line-style", GDK_LINE_SOLID,
                 NULL);

    g_object_set(apopcontrols->cmodelthresh,
                 "mode", GWY_GRAPH_CURVE_LINE,
                 "description", "Threshold",
                 "color", CURVE_COLOR_TH,
                 "line-style", GDK_LINE_SOLID,
                 NULL);

    g_object_set(apopcontrols->cmodelthresh2,
                 "mode", GWY_GRAPH_CURVE_LINE,
                 "description", "Threshold minimum",
                 "color", CURVE_COLOR_TH2,
                 "line-style", GDK_LINE_SOLID,
                 NULL);

    apopcontrols->ctrls = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* Parameters section */

    apopcontrols->frame_parameters = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(apopcontrols->frame_parameters), GTK_WIDGET(gwy_label_new_header("Parameters")));
    apopcontrols->table_parameters = gtk_table_new(5, 2, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(apopcontrols->table_parameters), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(apopcontrols->table_parameters), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(apopcontrols->table_parameters), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(apopcontrols->frame_parameters), apopcontrols->table_parameters);

    i = 0;

    gtk_table_attach(GTK_TABLE(apopcontrols->table_parameters), label_new_with_markup_left("Moving average window width"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    apopcontrols->nmove = gtk_adjustment_new(args->apopdata.nmove, 1, 101, 2, 10, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(apopcontrols->nmove), 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_spin_button_set_snap_to_ticks(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(apopcontrols->table_parameters), spin, 1, 2, i, i + 1,  GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(apopcontrols->nmove, "value-changed", G_CALLBACK(nmove_changed), data);
    i++;

    gtk_table_attach(GTK_TABLE(apopcontrols->table_parameters), label_new_with_markup_left("Derivative threshhold (identification)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    apopcontrols->thresh = gtk_adjustment_new(args->apopdata.thresh, 0, 10000, 0.001, 0.1, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(apopcontrols->thresh), 0.01, 3);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_spin_button_set_snap_to_ticks(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(apopcontrols->table_parameters), spin, 1, 2, i, i + 1,  GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(apopcontrols->thresh, "value-changed", G_CALLBACK(thresh_changed), data);
    i++;

    gtk_table_attach(GTK_TABLE(apopcontrols->table_parameters), label_new_with_markup_left("Derivative threshhold (width)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    apopcontrols->thresh2 = gtk_adjustment_new(args->apopdata.thresh2, 0, 10000, 0.001, 0.1, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(apopcontrols->thresh2), 0.01, 3);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_spin_button_set_snap_to_ticks(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(apopcontrols->table_parameters), spin, 1, 2, i, i + 1,  GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(apopcontrols->thresh2, "value-changed", G_CALLBACK(thresh2_changed), data);
    i++;

    gtk_table_attach(GTK_TABLE(apopcontrols->table_parameters), label_new_with_markup_left("Minimum pop-in width"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    apopcontrols->wpop = gtk_adjustment_new(args->apopdata.wpop, 0, 100, 1, 10, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(apopcontrols->wpop), 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_spin_button_set_snap_to_ticks(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(apopcontrols->table_parameters), spin, 1, 2, i, i + 1,  GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(apopcontrols->wpop, "value-changed", G_CALLBACK(wpop_changed), data);
    i++;

    gtk_table_attach(GTK_TABLE(apopcontrols->table_parameters), label_new_with_markup_left("Minimum pop-in height"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    apopcontrols->hpop = gtk_adjustment_new(args->apopdata.hpop, 0, 100, 0.001, 0.1, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(apopcontrols->hpop), 0.01, 3);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_spin_button_set_snap_to_ticks(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(apopcontrols->table_parameters), spin, 1, 2, i, i + 1,  GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(apopcontrols->hpop, "value-changed", G_CALLBACK(hpop_changed), data);
    i++;

    /* Reset button */
    apopcontrols->button_reset =  gtk_button_new_with_label("Reset to defaults");
    gtk_table_attach(GTK_TABLE(apopcontrols->table_parameters), apopcontrols->button_reset, 0, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(apopcontrols->button_reset, "clicked", G_CALLBACK(reset_default), data);

    /* Run button */
    apopcontrols->button_run = gtk_button_new_with_mnemonic(gwy_sgettext("verb|_Run"));
    gtk_button_set_image(GTK_BUTTON(apopcontrols->button_run), gtk_image_new_from_stock(GTK_STOCK_MEDIA_PLAY, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(apopcontrols->button_run, "clicked", G_CALLBACK(apopin_calc_and_update), data);
    //	g_signal_connect_swapped (apopcontrols->button_run, "clicked", G_CALLBACK (apopin_run), data);

    /* Results section */

    apopcontrols->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(apopcontrols->frame_results), GTK_WIDGET(gwy_label_new_header("Results")));

    apopcontrols->scrolledwindow_results = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(apopcontrols->scrolledwindow_results), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(apopcontrols->frame_results), apopcontrols->scrolledwindow_results);

    apopcontrols->table_results = gtk_table_new(2, 2, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(apopcontrols->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(apopcontrols->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(apopcontrols->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(apopcontrols->scrolledwindow_results), apopcontrols->table_results);

    gtk_table_attach(GTK_TABLE(apopcontrols->table_results), label_new_with_markup_left("Found pop-ins:"), 0, 1, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    apopcontrols->npopin = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(apopcontrols->npopin), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(apopcontrols->table_results), data->controls->apopcontrols->npopin, 1, 2, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    apopcontrols->res = gtk_table_new(1, 4, FALSE);
    gtk_table_attach(GTK_TABLE(apopcontrols->table_results), apopcontrols->res, 0, 2, 1, 2, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    /* Save button */
    apopcontrols->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(apopcontrols->button_save, "clicked", G_CALLBACK(apopin_save_data), args);

    /* Zooming */

    apopcontrols->vbox_zoom_buttons_h = gtk_vbox_new(TRUE, 0);
    apopcontrols->vbox_zoom_outer_h = gtk_vbox_new(FALSE, 0);

    apopcontrols->button_zoom_in_h = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(apopcontrols->button_zoom_in_h), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(apopcontrols->vbox_zoom_buttons_h), apopcontrols->button_zoom_in_h, TRUE, TRUE, 0);
    g_signal_connect_swapped(apopcontrols->button_zoom_in_h, "clicked", G_CALLBACK(apop_zoom_in_linked_h), data);

    apopcontrols->button_zoom_out_h = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(apopcontrols->button_zoom_out_h), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(apopcontrols->vbox_zoom_buttons_h), apopcontrols->button_zoom_out_h, TRUE, TRUE, 0);
    g_signal_connect_swapped(apopcontrols->button_zoom_out_h, "clicked", G_CALLBACK(apop_zoom_out_linked), data);

    apopcontrols->button_zoom_restore_h = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(apopcontrols->button_zoom_restore_h), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(apopcontrols->vbox_zoom_buttons_h), apopcontrols->button_zoom_restore_h, TRUE, TRUE, 0);
    g_signal_connect_swapped(apopcontrols->button_zoom_restore_h, "clicked", G_CALLBACK(apop_zoom_restore_linked), data);

    apopcontrols->vbox_zoom_buttons_der = gtk_vbox_new(TRUE, 0);
    apopcontrols->vbox_zoom_outer_der = gtk_vbox_new(FALSE, 0);

    apopcontrols->button_zoom_in_der = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(apopcontrols->button_zoom_in_der), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(apopcontrols->vbox_zoom_buttons_der), apopcontrols->button_zoom_in_der, TRUE, TRUE, 0);
    g_signal_connect_swapped(apopcontrols->button_zoom_in_der, "clicked", G_CALLBACK(apop_zoom_in_linked_der), data);

    apopcontrols->button_zoom_out_der = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(apopcontrols->button_zoom_out_der), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(apopcontrols->vbox_zoom_buttons_der), apopcontrols->button_zoom_out_der, TRUE, TRUE, 0);
    g_signal_connect_swapped(apopcontrols->button_zoom_out_der, "clicked", G_CALLBACK(apop_zoom_out_linked), data);

    apopcontrols->button_zoom_restore_der = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(apopcontrols->button_zoom_restore_der), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(apopcontrols->vbox_zoom_buttons_der), apopcontrols->button_zoom_restore_der, TRUE, TRUE, 0);
    g_signal_connect_swapped(apopcontrols->button_zoom_restore_der, "clicked", G_CALLBACK(apop_zoom_restore_linked), data);

    //final assembly
    apopcontrols->vbox_graphs = gtk_vbox_new(FALSE, CTRLS_SPACING);
    apopcontrols->hbox_graph_h = gtk_hbox_new(FALSE, CTRLS_SPACING);
    apopcontrols->hbox_graph_der = gtk_hbox_new(FALSE, CTRLS_SPACING);

    gtk_box_pack_start(GTK_BOX(apopcontrols->hbox_graph_h), apopcontrols->graph_h, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(apopcontrols->hbox_graph_h), apopcontrols->vbox_zoom_outer_h, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(apopcontrols->hbox_graph_der), apopcontrols->graph_der, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(apopcontrols->hbox_graph_der), apopcontrols->vbox_zoom_outer_der, FALSE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(apopcontrols->vbox_zoom_outer_h), apopcontrols->vbox_zoom_buttons_h, TRUE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(apopcontrols->vbox_zoom_outer_der), apopcontrols->vbox_zoom_buttons_der, TRUE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(apopcontrols->ctrls), apopcontrols->frame_parameters, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(apopcontrols->ctrls), apopcontrols->button_run, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(apopcontrols->ctrls), apopcontrols->frame_results, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(apopcontrols->ctrls), apopcontrols->button_save, FALSE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(apopcontrols->vbox_graphs), apopcontrols->hbox_graph_h, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(apopcontrols->vbox_graphs), apopcontrols->hbox_graph_der, TRUE, TRUE, 0);

    apopcontrols->apop_gui = gtk_hbox_new(FALSE, CTRLS_SPACING);

    gtk_box_pack_start(GTK_BOX(apopcontrols->apop_gui), apopcontrols->ctrls, FALSE, FALSE, 0);   /* expand, fill, padding */
    gtk_box_pack_start(GTK_BOX(apopcontrols->apop_gui), apopcontrols->vbox_graphs, TRUE, TRUE, 0);

    gtk_widget_show(apopcontrols->apop_gui);

    return apopcontrols->apop_gui;
}

static void update_models_apop_h(ApopinControls *apopcontrols)
{
    GwyGraphCurveModel *cmodel;

    apopcontrols->graph_model = gwy_graph_get_model(GWY_GRAPH(apopcontrols->graph_h));
    cmodel = gwy_graph_model_get_curve_by_description(apopcontrols->graph_model, "Depth");

    if (cmodel) {
        apopcontrols->cmodelload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(apopcontrols->graph_model, "Smoothed depth");

    if (cmodel) {
        apopcontrols->cmodelhavg = cmodel;
    }
}

static void apop_zoom_in_linked_h(Data *data)
{
    ApopinControls *apopcontrols;
    gdouble range[2];
    GwySelection *selection;

    apopcontrols = data->controls->apopcontrols;

    selection = gwy_graph_area_get_selection(GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(apopcontrols->graph_h))), GWY_GRAPH_STATUS_PLAIN);

    if (gwy_selection_get_object(selection, 0, range)) {
        selection = gwy_graph_area_get_selection(GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(apopcontrols->graph_der))), GWY_GRAPH_STATUS_PLAIN);
        gwy_selection_set_object(selection, 0, range);

        rs_zoom_in(GWY_GRAPH(apopcontrols->graph_h), &(apopcontrols->zoomstack_h), FALSE);   /* do not rescale y */
        rs_zoom_in(GWY_GRAPH(apopcontrols->graph_der), &(apopcontrols->zoomstack_der), FALSE);   /* do not rescale y */
        update_models_apop_h(apopcontrols);
        update_models_apop_der(apopcontrols);
    }
}

static void update_models_apop_der(ApopinControls *apopcontrols)
{
    GwyGraphCurveModel *cmodel;

    apopcontrols->graph_model_der = gwy_graph_get_model(GWY_GRAPH(apopcontrols->graph_der));
    cmodel = gwy_graph_model_get_curve_by_description(apopcontrols->graph_model_der, "Derivative");

    if (cmodel) {
        apopcontrols->cmodelder = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(apopcontrols->graph_model_der, "Threshold");

    if (cmodel) {
        apopcontrols->cmodelthresh = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(apopcontrols->graph_model_der, "Threshold minimum");

    if (cmodel) {
        apopcontrols->cmodelthresh2 = cmodel;
    }
}

static void apop_zoom_in_linked_der(Data *data)
{
    ApopinControls *apopcontrols;
    gdouble range[2];
    GwySelection *selection;

    apopcontrols = data->controls->apopcontrols;

    selection = gwy_graph_area_get_selection(GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(apopcontrols->graph_der))), GWY_GRAPH_STATUS_PLAIN);

    if (gwy_selection_get_object(selection, 0, range)) {
        selection = gwy_graph_area_get_selection(GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(apopcontrols->graph_h))), GWY_GRAPH_STATUS_PLAIN);
        gwy_selection_set_object(selection, 0, range);

        rs_zoom_in(GWY_GRAPH(apopcontrols->graph_der), &(apopcontrols->zoomstack_der), FALSE);   /* do not rescale y */
        rs_zoom_in(GWY_GRAPH(apopcontrols->graph_h), &(apopcontrols->zoomstack_h), FALSE);   /* do not rescale y */
        update_models_apop_der(apopcontrols);
        update_models_apop_h(apopcontrols);
    }
}

static void apop_zoom_out_linked(Data *data)
{
    ApopinControls *apopcontrols;

    apopcontrols = data->controls->apopcontrols;

    rs_zoom_out(GWY_GRAPH(apopcontrols->graph_h), &(apopcontrols->zoomstack_h));
    rs_zoom_out(GWY_GRAPH(apopcontrols->graph_der), &(apopcontrols->zoomstack_der));
    update_models_apop_h(apopcontrols);
    update_models_apop_der(apopcontrols);
}

void apop_zoom_restore_linked(Data *data)
{
    ApopinControls *apopcontrols;

    apopcontrols = data->controls->apopcontrols;

    rs_zoom_restore(GWY_GRAPH(apopcontrols->graph_h), &(apopcontrols->zoomstack_h));
    rs_zoom_restore(GWY_GRAPH(apopcontrols->graph_der), &(apopcontrols->zoomstack_der));
    update_models_apop_h(apopcontrols);
    update_models_apop_der(apopcontrols);
}

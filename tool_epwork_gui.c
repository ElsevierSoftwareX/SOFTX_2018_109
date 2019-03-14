#include "tool_epwork_gui.h"
#include "tool_epwork.h"
#include "tool_epwork_unc_gui.h"

#include "controls.h"
#include "datatypes.h"
#include "file-utils.h"
#include "gui-utils.h"
#include "niget-common.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>
#include <app/file.h>

static void work_remove_labels(WorkControls *workcontrols);

static void update_models_work(WorkControls *workcontrols);
static void work_zoom_in(Data *data);
static void work_zoom_out(Data *data);

static void work_run(Data *data)
{
    FDdata *fddata;
    Workdata *workdata;

    // local copies
    fddata = &(data->args->fddata);
    workdata = &(data->args->workdata);

    work_calc(workdata, fddata);

    work_zoom_restore(data);

    // update smoothed loading curve
    gwy_graph_curve_model_set_data(data->controls->workcontrols->cmodelloadavg, workdata->hloadavg->data, workdata->Floadavg->data, workdata->hloadavg->res);

    // update smoothed holding curve
    gwy_graph_curve_model_set_data(data->controls->workcontrols->cmodelholdavg, workdata->hholdavg->data, workdata->Fholdavg->data, workdata->hholdavg->res);

    // update smoothed unloading curve
    gwy_graph_curve_model_set_data(data->controls->workcontrols->cmodelunloadavg, workdata->hunloadavg->data, workdata->Funloadavg->data, workdata->hunloadavg->res);

    //update labels
    label_set_gdouble_format(data->controls->workcontrols->workelastic, workdata->workelastic, EPWORK);
    label_set_gdouble_format(data->controls->workcontrols->workplastic, workdata->workplastic, EPWORK);
    label_set_gdouble_format(data->controls->workcontrols->eta, workdata->eta, "%.1f");
}

static void nmove_changed(Data *data)
{
    data->args->workdata.nmove = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->workcontrols->nmove));
    work_run(data);
}

void work_remove_fit_results_labels(Data *data)
{
    work_remove_results(&(data->args->workdata));
    work_remove_labels(data->controls->workcontrols);

    if (data->args->workunc.We != NULL || data->args->workunc.Wp != NULL) {
        work_unc_close_window(data);
    }
}

static void work_remove_labels(WorkControls *workcontrols)
{
    label_clear(workcontrols->workplastic);
    label_clear(workcontrols->workelastic);
    label_clear(workcontrols->eta);
}

static void work_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (work_has_results(&(args->workdata))) {
	buffer = work_export_data(&(args->workdata), &(args->fddata), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

void work_redraw(Data *data)
{
    //update loading, hold and unloading curves
    if (data->args->fddata.hload) {
        work_calc(&(data->args->workdata), &(data->args->fddata));
        gwy_graph_curve_model_set_data(data->controls->workcontrols->cmodelloadavg,
                                       data->args->workdata.hloadavg->data, data->args->workdata.Floadavg->data,
                                       data->args->workdata.hloadavg->res);
        gwy_graph_curve_model_set_data(data->controls->workcontrols->cmodelholdavg,
                                       data->args->workdata.hholdavg->data, data->args->workdata.Fholdavg->data,
                                       data->args->workdata.hholdavg->res);
        gwy_graph_curve_model_set_data(data->controls->workcontrols->cmodelunloadavg,
                                       data->args->workdata.hunloadavg->data, data->args->workdata.Funloadavg->data,
                                       data->args->workdata.hunloadavg->res);
    }
}

GtkWidget *tool_epwork_create(Data *data)
{
    GtkWidget *spin;
    Args *args;
    WorkControls *workcontrols;
    gint i;

    args = data->args;
    workcontrols = data->controls->workcontrols;

    // GraphArea: loading, hold and unloading curves, and their smoothed counterparts
    //
    workcontrols->graph_model = gwy_graph_model_new();
    g_object_set(workcontrols->graph_model, "title", "Load-Hold-Unload", NULL);
    gwy_graph_model_set_axis_label(workcontrols->graph_model, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(workcontrols->graph_model, GTK_POS_LEFT, "F / mN");

    workcontrols->zoomstack = NULL;
    workcontrols->zoomstack = g_slist_prepend(workcontrols->zoomstack, workcontrols->graph_model);

    workcontrols->hbox_graph = gtk_hbox_new(FALSE, CTRLS_SPACING);

    workcontrols->graph = gwy_graph_new(workcontrols->graph_model);
    gtk_widget_set_size_request(GTK_WIDGET(workcontrols->graph), 600, 400);

    gwy_graph_enable_user_input(GWY_GRAPH(workcontrols->graph), FALSE);
    gwy_graph_set_status(GWY_GRAPH(workcontrols->graph), GWY_GRAPH_STATUS_XSEL);

    workcontrols->cmodelload = gwy_graph_curve_model_new();
    workcontrols->cmodelloadavg = gwy_graph_curve_model_new();
    workcontrols->cmodelhold = gwy_graph_curve_model_new();
    workcontrols->cmodelholdavg = gwy_graph_curve_model_new();
    workcontrols->cmodelunload = gwy_graph_curve_model_new();
    workcontrols->cmodelunloadavg = gwy_graph_curve_model_new();

    g_object_set(workcontrols->cmodelload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Loading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_LOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(workcontrols->cmodelloadavg,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Smoothed loading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color",  CURVE_COLOR_LOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(workcontrols->cmodelhold,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Hold",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_HOLD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(workcontrols->cmodelholdavg,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Smoothed hold",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_HOLD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(workcontrols->cmodelunload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Unloading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_UNLOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(workcontrols->cmodelunloadavg,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Smoothed unloading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_UNLOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    /* only the smoothed curves are shown */
    gwy_graph_model_add_curve(workcontrols->graph_model, workcontrols->cmodelloadavg);
    gwy_graph_model_add_curve(workcontrols->graph_model, workcontrols->cmodelholdavg);
    gwy_graph_model_add_curve(workcontrols->graph_model, workcontrols->cmodelunloadavg);

    // Table of input data and results
    workcontrols->ctrls = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* Parameters section */

    workcontrols->frame_parameters = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(workcontrols->frame_parameters), GTK_WIDGET(gwy_label_new_header("Parameters")));
    workcontrols->table_parameters = gtk_table_new(1, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(workcontrols->table_parameters), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(workcontrols->table_parameters), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(workcontrols->table_parameters), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(workcontrols->frame_parameters), workcontrols->table_parameters);
    i = 0;

    gtk_table_attach(GTK_TABLE(workcontrols->table_parameters), label_new_with_markup_left("Moving average window width"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    workcontrols->nmove = gtk_adjustment_new(args->workdata.nmove, 1, 101, 2, 10, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(workcontrols->nmove), 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_spin_button_set_snap_to_ticks(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(workcontrols->table_parameters), spin, 1, 2, i, i + 1,  GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(workcontrols->nmove, "value-changed", G_CALLBACK(nmove_changed), data);

    /* Run button */
    workcontrols->button_run = gtk_button_new_with_mnemonic(gwy_sgettext("verb|_Run"));
    gtk_button_set_image(GTK_BUTTON(workcontrols->button_run), gtk_image_new_from_stock(GTK_STOCK_MEDIA_PLAY, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(workcontrols->button_run, "clicked", G_CALLBACK(work_run), data);

    /* Results section */

    workcontrols->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(workcontrols->frame_results), GTK_WIDGET(gwy_label_new_header("Results")));
    workcontrols->table_results = gtk_table_new(3, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(workcontrols->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(workcontrols->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(workcontrols->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(workcontrols->frame_results), workcontrols->table_results);
    i = 0;

    gtk_table_attach(GTK_TABLE(workcontrols->table_results), label_new_with_markup_left("Plastic work"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    workcontrols->workplastic = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(workcontrols->workplastic), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(workcontrols->table_results), workcontrols->workplastic, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(workcontrols->table_results), label_new_with_markup_left("pJ"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(workcontrols->table_results), label_new_with_markup_left("Elastic work"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    workcontrols->workelastic = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(workcontrols->workelastic), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(workcontrols->table_results), workcontrols->workelastic, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(workcontrols->table_results), label_new_with_markup_left("pJ"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(workcontrols->table_results), label_new_with_markup_left("η"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    workcontrols->eta = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(workcontrols->eta), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(workcontrols->table_results), workcontrols->eta, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(workcontrols->table_results), label_new_with_markup_left("%"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    /* Unc button */
    /*
    workcontrols->button_unc = gtk_button_new_with_mnemonic("_Uncertainties");
    gtk_button_set_image (GTK_BUTTON (workcontrols->button_unc), gtk_image_new_from_stock (GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped (workcontrols->button_unc, "clicked", G_CALLBACK (work_uncertainty), data);
    */

    /* Save button */
    workcontrols->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(workcontrols->button_save, "clicked", G_CALLBACK(work_save_data), args);

    gtk_box_pack_start(GTK_BOX(workcontrols->ctrls), workcontrols->frame_parameters, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(workcontrols->ctrls), workcontrols->button_run, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(workcontrols->ctrls), workcontrols->frame_results, FALSE, FALSE, 0);
    /*gtk_box_pack_start (GTK_BOX (workcontrols->ctrls), workcontrols->button_unc, FALSE, FALSE, 0); */
    gtk_box_pack_start(GTK_BOX(workcontrols->ctrls), workcontrols->button_save, FALSE, FALSE, 0);

    /* Zooming */

    workcontrols->vbox_zoom_buttons = gtk_vbox_new(TRUE, 0);
    workcontrols->vbox_zoom_outer = gtk_vbox_new(FALSE, 0);

    workcontrols->button_zoom_in = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(workcontrols->button_zoom_in), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(workcontrols->vbox_zoom_buttons), workcontrols->button_zoom_in, TRUE, TRUE, 0);
    g_signal_connect_swapped(workcontrols->button_zoom_in, "clicked", G_CALLBACK(work_zoom_in), data);

    workcontrols->button_zoom_out = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(workcontrols->button_zoom_out), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(workcontrols->vbox_zoom_buttons), workcontrols->button_zoom_out, TRUE, TRUE, 0);
    g_signal_connect_swapped(workcontrols->button_zoom_out, "clicked", G_CALLBACK(work_zoom_out), data);

    workcontrols->button_zoom_restore = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(workcontrols->button_zoom_restore), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(workcontrols->vbox_zoom_buttons), workcontrols->button_zoom_restore, TRUE, TRUE, 0);
    g_signal_connect_swapped(workcontrols->button_zoom_restore, "clicked", G_CALLBACK(work_zoom_restore), data);

    gtk_box_pack_start(GTK_BOX(workcontrols->vbox_zoom_outer), workcontrols->vbox_zoom_buttons, TRUE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(workcontrols->hbox_graph), workcontrols->graph, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(workcontrols->hbox_graph), workcontrols->vbox_zoom_outer, FALSE, FALSE, 0);


    //final

    workcontrols->work_gui = gtk_hbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(workcontrols->work_gui), workcontrols->ctrls, FALSE, FALSE, 0);   /* expand, fill, padding */
    gtk_box_pack_start(GTK_BOX(workcontrols->work_gui), workcontrols->hbox_graph, TRUE, TRUE, 0);

    gtk_widget_show(workcontrols->work_gui);

    return workcontrols->work_gui;
}

static void update_models_work(WorkControls *workcontrols)
{
    GwyGraphCurveModel *cmodel;

    workcontrols->graph_model = gwy_graph_get_model(GWY_GRAPH(workcontrols->graph));
    cmodel = gwy_graph_model_get_curve_by_description(workcontrols->graph_model, "Loading");

    if (cmodel) {
        workcontrols->cmodelload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(workcontrols->graph_model, "Smoothed loading");

    if (cmodel) {
        workcontrols->cmodelloadavg = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(workcontrols->graph_model, "Hold");

    if (cmodel) {
        workcontrols->cmodelhold = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(workcontrols->graph_model, "Smoothed hold");

    if (cmodel) {
        workcontrols->cmodelholdavg = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(workcontrols->graph_model, "Unloading");

    if (cmodel) {
        workcontrols->cmodelunload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(workcontrols->graph_model, "Smoothed unloading");

    if (cmodel) {
        workcontrols->cmodelunloadavg = cmodel;
    }
}

static void work_zoom_in(Data *data)
{
    WorkControls *workcontrols;

    workcontrols = data->controls->workcontrols;

    rs_zoom_in(GWY_GRAPH(workcontrols->graph), &(workcontrols->zoomstack), FALSE);   /* do not rescale y */
    update_models_work(workcontrols);
}

static void work_zoom_out(Data *data)
{
    WorkControls *workcontrols;

    workcontrols = data->controls->workcontrols;

    rs_zoom_out(GWY_GRAPH(workcontrols->graph), &(workcontrols->zoomstack));
    update_models_work(workcontrols);

}

void work_zoom_restore(Data *data)
{
    WorkControls *workcontrols;

    workcontrols = data->controls->workcontrols;

    rs_zoom_restore(GWY_GRAPH(workcontrols->graph), &(workcontrols->zoomstack));
    update_models_work(workcontrols);
}

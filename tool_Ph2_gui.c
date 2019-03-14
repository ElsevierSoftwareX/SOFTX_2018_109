#include "tool_Ph2_gui.h"
#include "tool_Ph2.h"

#include "controls.h"
#include "datatypes.h"
#include "gui-utils.h"
#include "niget-common.h"

#include <math.h>

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>
#include <app/file.h>
#include <libprocess/gwyprocess.h>


static void Ph2_remove_curves(Ph2Controls *ph2controls);

static void update_models_ph2_Fh(Ph2Controls *ph2controls);
static void ph2_zoom_in_linked_Fh(Data *data);

static void update_models_ph2_H(Ph2Controls *ph2controls);
static void ph2_zoom_in_linked_H(Data *data);

static void ph2_zoom_out_linked(Data *data);


void Ph2_remove_fit_results_labels(Data *data)
{
    Ph2_remove_curves(data->controls->ph2controls);
}

static void Ph2_remove_curves(Ph2Controls *ph2controls)
{
    if (gwy_graph_model_get_n_curves(ph2controls->graph_modelFh) > 1) { // je potreba testovat pocet?
        gwy_graph_model_remove_curve_by_description(ph2controls->graph_modelFh, "Smoothed loading");
    }

    if (gwy_graph_model_get_n_curves(ph2controls->graph_modelPh2) == 2) { // dtto
        gwy_graph_model_remove_curve_by_description(ph2controls->graph_modelPh2, "P/h<sup>2</sup>");
        gwy_graph_model_remove_curve_by_description(ph2controls->graph_modelPh2, "dP/dh<sup>2</sup>");
    }
}

static void Ph2_calc_and_update(Data *data)
{
    GwyDataLine *Ph2, *dPdh2, *havg;
    Ph2data *ph2data;
    FDdata *fddata;
    gint n;
    GdkColor color;

    gdk_color_parse("red", &color);

    // local copies
    fddata = &(data->args->fddata);
    ph2data = &(data->args->ph2data);
    n = fddata->hload->res;

    Ph2_calc(fddata, ph2data);

    ph2_zoom_restore_linked(data);

    // update smoothed loading curve, if necessary add it to graph
    if (gwy_graph_model_get_n_curves(data->controls->ph2controls->graph_modelFh) == 1) {
        gwy_graph_model_add_curve(data->controls->ph2controls->graph_modelFh, data->controls->ph2controls->cmodelavg);
    }

    gwy_graph_curve_model_set_data(data->controls->ph2controls->cmodelavg, ph2data->havg->data, ph2data->Favg->data, n);

    havg = copy_data_remove_ind(ph2data->havg, ph2data->ind_remove);
    Ph2 = copy_data_remove_ind(ph2data->Ph2,  ph2data->ind_remove);
    dPdh2 = copy_data_remove_ind(ph2data->dPdh2, ph2data->ind_remove);

    // update curves, if necessary add them to graph
    if (!gwy_graph_model_get_n_curves(data->controls->ph2controls->graph_modelPh2)) {
        gwy_graph_model_add_curve(data->controls->ph2controls->graph_modelPh2, data->controls->ph2controls->cmodelPh2);
        gwy_graph_model_add_curve(data->controls->ph2controls->graph_modelPh2, data->controls->ph2controls->cmodeldPdh2);
    }

    gwy_graph_curve_model_set_data(data->controls->ph2controls->cmodelPh2, havg->data, Ph2->data, havg->res);
    gwy_graph_curve_model_set_data(data->controls->ph2controls->cmodeldPdh2, havg->data, dPdh2->data, havg->res);

    g_object_unref(havg);
    g_object_unref(Ph2);
    g_object_unref(dPdh2);
}

static void nmove_changed(Data *data)
{
    data->args->ph2data.nmove = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->ph2controls->nmove));
    Ph2_calc_and_update(data);
}

static void Ph2_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (Ph2_has_results(&(args->ph2data))) {
	buffer = Ph2_export_data(&(args->ph2data), &(args->fddata), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

void Ph2_redraw(Data *data)
{
    GwyDataLine *havg, *dPdh2, *Ph2;

    //update loading and smoothed loading curve, update P/h^2 and dP/dh^2 curves (if available)
    if (data->args->fddata.hload) {
        gwy_graph_curve_model_set_data(data->controls->ph2controls->cmodelload,
                                       data->args->fddata.hload->data, data->args->fddata.Fload->data,
                                       data->args->fddata.hload->res);

        if (gwy_graph_model_get_n_curves(data->controls->ph2controls->graph_modelFh) == 2) {
            gwy_graph_curve_model_set_data(data->controls->ph2controls->cmodelavg,
                                           data->args->ph2data.havg->data, data->args->ph2data.Favg->data,
                                           data->args->ph2data.havg->res);
        }


        if (gwy_graph_model_get_n_curves(data->controls->ph2controls->graph_modelPh2)) {
            havg = copy_data_remove_ind(data->args->ph2data.havg,  data->args->ph2data.ind_remove);
            Ph2 = copy_data_remove_ind(data->args->ph2data.Ph2,  data->args->ph2data.ind_remove);
            dPdh2 = copy_data_remove_ind(data->args->ph2data.dPdh2, data->args->ph2data.ind_remove);

            gwy_graph_curve_model_set_data(data->controls->ph2controls->cmodelPh2,
                                           havg->data, Ph2->data, havg->res);
            gwy_graph_curve_model_set_data(data->controls->ph2controls->cmodeldPdh2,
                                           havg->data, dPdh2->data, havg->res);

            g_object_unref(havg);
            g_object_unref(Ph2);
            g_object_unref(dPdh2);
        }
    }
}

GtkWidget *tool_Ph2_create(Data *data)
{
    GtkWidget *spin;
    Args *args;
    Ph2Controls *ph2controls;
    gint i;

    args = data->args;
    ph2controls = data->controls->ph2controls;

    // GraphArea, two graph areas: one for the loading curve and the smoother loading curve, one for P/h^2 and dP/dh^2
    ph2controls->graph_modelPh2 = gwy_graph_model_new();
    g_object_set(ph2controls->graph_modelPh2, "title", "P-h<sup>2</sup>", NULL);
    gwy_graph_model_set_axis_label(ph2controls->graph_modelPh2, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(ph2controls->graph_modelPh2, GTK_POS_LEFT, "P/h<sup>2</sup> / MPa");

    ph2controls->zoomstack_H = NULL;
    ph2controls->zoomstack_H = g_slist_prepend(ph2controls->zoomstack_H, ph2controls->graph_modelPh2);

    ph2controls->graph_modelFh = gwy_graph_model_new();
    g_object_set(ph2controls->graph_modelFh, "title", "Fh", NULL);
    gwy_graph_model_set_axis_label(ph2controls->graph_modelFh, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(ph2controls->graph_modelFh, GTK_POS_LEFT, "F / mN");

    ph2controls->zoomstack_Fh = NULL;
    ph2controls->zoomstack_Fh = g_slist_prepend(ph2controls->zoomstack_Fh, ph2controls->graph_modelFh);

    ph2controls->graph_H = gwy_graph_new(ph2controls->graph_modelPh2);
    gtk_widget_set_size_request(GTK_WIDGET(ph2controls->graph_H), 600, 200);
    ph2controls->graph_Fh = gwy_graph_new(ph2controls->graph_modelFh);
    gtk_widget_set_size_request(GTK_WIDGET(ph2controls->graph_Fh), 600, 200);

    gwy_graph_enable_user_input(GWY_GRAPH(ph2controls->graph_H), FALSE);
    gwy_graph_set_status(GWY_GRAPH(ph2controls->graph_H), GWY_GRAPH_STATUS_XSEL);
    gwy_graph_enable_user_input(GWY_GRAPH(ph2controls->graph_Fh), FALSE);
    gwy_graph_set_status(GWY_GRAPH(ph2controls->graph_Fh), GWY_GRAPH_STATUS_XSEL);

    ph2controls->cmodelload = gwy_graph_curve_model_new();
    ph2controls->cmodelavg = gwy_graph_curve_model_new();
    ph2controls->cmodelPh2 = gwy_graph_curve_model_new();
    ph2controls->cmodeldPdh2 = gwy_graph_curve_model_new();

    g_object_set(ph2controls->cmodelload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Loading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_LOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(ph2controls->cmodelavg,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Smoothed loading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_SMOOTH_LOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(ph2controls->cmodelPh2,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "P/h<sup>2</sup>",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_PH2,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    g_object_set(ph2controls->cmodeldPdh2,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "dP/dh<sup>2</sup>",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_DIFF_PH2,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);

    gwy_graph_model_add_curve(ph2controls->graph_modelFh, ph2controls->cmodelload);

    ph2controls->ctrls = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* Parameters section */

    ph2controls->frame_parameters = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ph2controls->frame_parameters), GTK_WIDGET(gwy_label_new_header("Parameters")));
    ph2controls->table_parameters = gtk_table_new(1, 2, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(ph2controls->table_parameters), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ph2controls->table_parameters), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ph2controls->table_parameters), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ph2controls->frame_parameters), ph2controls->table_parameters);

    i = 0;

    /* Mov avg win width */

    gtk_table_attach(GTK_TABLE(ph2controls->table_parameters), label_new_with_markup_left("Moving average window width"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    ph2controls->nmove = gtk_adjustment_new(args->ph2data.nmove, 1, 101, 2, 10, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ph2controls->nmove), 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_spin_button_set_snap_to_ticks(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ph2controls->table_parameters), spin, 1, 2, i, i + 1,  GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ph2controls->nmove, "value-changed", G_CALLBACK(nmove_changed), data);

    /* Run button */
    ph2controls->button_run = gtk_button_new_with_mnemonic(gwy_sgettext("verb|_Run"));
    gtk_button_set_image(GTK_BUTTON(ph2controls->button_run), gtk_image_new_from_stock(GTK_STOCK_MEDIA_PLAY, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ph2controls->button_run, "clicked", G_CALLBACK(Ph2_calc_and_update), data);

    /* Save button */

    ph2controls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ph2controls->button_save, "clicked", G_CALLBACK(Ph2_save_data), args);

    /* Zooming */

    ph2controls->vbox_zoom_buttons_Fh = gtk_vbox_new(TRUE, 0);
    ph2controls->vbox_zoom_outer_Fh = gtk_vbox_new(FALSE, 0);

    ph2controls->button_zoom_in_Fh = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ph2controls->button_zoom_in_Fh), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(ph2controls->vbox_zoom_buttons_Fh), ph2controls->button_zoom_in_Fh, TRUE, TRUE, 0);
    g_signal_connect_swapped(ph2controls->button_zoom_in_Fh, "clicked", G_CALLBACK(ph2_zoom_in_linked_Fh), data);

    ph2controls->button_zoom_out_Fh = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ph2controls->button_zoom_out_Fh), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(ph2controls->vbox_zoom_buttons_Fh), ph2controls->button_zoom_out_Fh, TRUE, TRUE, 0);
    g_signal_connect_swapped(ph2controls->button_zoom_out_Fh, "clicked", G_CALLBACK(ph2_zoom_out_linked), data);

    ph2controls->button_zoom_restore_Fh = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ph2controls->button_zoom_restore_Fh), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(ph2controls->vbox_zoom_buttons_Fh), ph2controls->button_zoom_restore_Fh, TRUE, TRUE, 0);
    g_signal_connect_swapped(ph2controls->button_zoom_restore_Fh, "clicked", G_CALLBACK(ph2_zoom_restore_linked), data);

    ph2controls->vbox_zoom_buttons_H = gtk_vbox_new(TRUE, 0);
    ph2controls->vbox_zoom_outer_H = gtk_vbox_new(FALSE, 0);

    ph2controls->button_zoom_in_H = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ph2controls->button_zoom_in_H), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(ph2controls->vbox_zoom_buttons_H), ph2controls->button_zoom_in_H, TRUE, TRUE, 0);
    g_signal_connect_swapped(ph2controls->button_zoom_in_H, "clicked", G_CALLBACK(ph2_zoom_in_linked_H), data);

    ph2controls->button_zoom_out_H = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ph2controls->button_zoom_out_H), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(ph2controls->vbox_zoom_buttons_H), ph2controls->button_zoom_out_H, TRUE, TRUE, 0);
    g_signal_connect_swapped(ph2controls->button_zoom_out_H, "clicked", G_CALLBACK(ph2_zoom_out_linked), data);

    ph2controls->button_zoom_restore_H = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ph2controls->button_zoom_restore_H), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(ph2controls->vbox_zoom_buttons_H), ph2controls->button_zoom_restore_H, TRUE, TRUE, 0);
    g_signal_connect_swapped(ph2controls->button_zoom_restore_H, "clicked", G_CALLBACK(ph2_zoom_restore_linked), data);


    //final

    gtk_box_pack_start(GTK_BOX(ph2controls->ctrls), ph2controls->frame_parameters, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ph2controls->ctrls), ph2controls->button_run, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ph2controls->ctrls), ph2controls->button_save, FALSE, FALSE, 0);

    ph2controls->vbox_graphs = gtk_vbox_new(FALSE, CTRLS_SPACING);
    ph2controls->hbox_graph_Fh = gtk_hbox_new(FALSE, CTRLS_SPACING);
    ph2controls->hbox_graph_H = gtk_hbox_new(FALSE, CTRLS_SPACING);

    gtk_box_pack_start(GTK_BOX(ph2controls->vbox_zoom_outer_Fh), ph2controls->vbox_zoom_buttons_Fh, TRUE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ph2controls->vbox_zoom_outer_H), ph2controls->vbox_zoom_buttons_H, TRUE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(ph2controls->hbox_graph_Fh), ph2controls->graph_Fh, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(ph2controls->hbox_graph_Fh), ph2controls->vbox_zoom_outer_Fh, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ph2controls->hbox_graph_H), ph2controls->graph_H, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(ph2controls->hbox_graph_H), ph2controls->vbox_zoom_outer_H, FALSE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(ph2controls->vbox_graphs), ph2controls->hbox_graph_Fh, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(ph2controls->vbox_graphs), ph2controls->hbox_graph_H, TRUE, TRUE, 0);

    ph2controls->ph2_gui = gtk_hbox_new(FALSE, CTRLS_SPACING);

    gtk_box_pack_start(GTK_BOX(ph2controls->ph2_gui), ph2controls->ctrls, FALSE, FALSE, 0);   /* expand, fill, padding */
    gtk_box_pack_start(GTK_BOX(ph2controls->ph2_gui), ph2controls->vbox_graphs, TRUE, TRUE, 0);

    gtk_widget_show(ph2controls->ph2_gui);

    return ph2controls->ph2_gui;
}

static void update_models_ph2_Fh(Ph2Controls *ph2controls)
{
    GwyGraphCurveModel *cmodel;

    ph2controls->graph_modelFh = gwy_graph_get_model(GWY_GRAPH(ph2controls->graph_Fh));
    cmodel = gwy_graph_model_get_curve_by_description(ph2controls->graph_modelFh, "Loading");

    if (cmodel) {
        ph2controls->cmodelload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(ph2controls->graph_modelFh, "Smoothed loading");

    if (cmodel) {
        ph2controls->cmodelavg = cmodel;
    }
}

static void ph2_zoom_in_linked_Fh(Data *data)
{
    Ph2Controls *ph2controls;
    gdouble range[2];
    GwySelection *selection;

    ph2controls = data->controls->ph2controls;

    selection = gwy_graph_area_get_selection(GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(ph2controls->graph_Fh))), GWY_GRAPH_STATUS_PLAIN);

    if (gwy_selection_get_object(selection, 0, range)) {
        selection = gwy_graph_area_get_selection(GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(ph2controls->graph_H))), GWY_GRAPH_STATUS_PLAIN);
        gwy_selection_set_object(selection, 0, range);

        rs_zoom_in(GWY_GRAPH(ph2controls->graph_Fh), &(ph2controls->zoomstack_Fh), FALSE);   /* do not rescale y */
        rs_zoom_in(GWY_GRAPH(ph2controls->graph_H), &(ph2controls->zoomstack_H), FALSE);   /* do not rescale y */
        update_models_ph2_Fh(ph2controls);
        update_models_ph2_H(ph2controls);
    }
}

static void update_models_ph2_H(Ph2Controls *ph2controls)
{
    GwyGraphCurveModel *cmodel;

    ph2controls->graph_modelPh2 = gwy_graph_get_model(GWY_GRAPH(ph2controls->graph_H));
    cmodel = gwy_graph_model_get_curve_by_description(ph2controls->graph_modelPh2, "P/h<sup>2</sup>");

    if (cmodel) {
        ph2controls->cmodelPh2 = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(ph2controls->graph_modelPh2, "dP/dh<sup>2</sup>");

    if (cmodel) {
        ph2controls->cmodeldPdh2 = cmodel;
    }
}

static void ph2_zoom_in_linked_H(Data *data)
{
    Ph2Controls *ph2controls;
    gdouble range[2];
    GwySelection *selection;

    ph2controls = data->controls->ph2controls;

    selection = gwy_graph_area_get_selection(GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(ph2controls->graph_H))), GWY_GRAPH_STATUS_PLAIN);

    if (gwy_selection_get_object(selection, 0, range)) {
        selection = gwy_graph_area_get_selection(GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(ph2controls->graph_Fh))), GWY_GRAPH_STATUS_PLAIN);
        gwy_selection_set_object(selection, 0, range);

        rs_zoom_in(GWY_GRAPH(ph2controls->graph_H), &(ph2controls->zoomstack_H), FALSE);   /* do not rescale y */
        rs_zoom_in(GWY_GRAPH(ph2controls->graph_Fh), &(ph2controls->zoomstack_Fh), FALSE);   /* do not rescale y */
        update_models_ph2_H(ph2controls);
        update_models_ph2_Fh(ph2controls);
    }
}

static void ph2_zoom_out_linked(Data *data)
{
    Ph2Controls *ph2controls;

    ph2controls = data->controls->ph2controls;

    rs_zoom_out(GWY_GRAPH(ph2controls->graph_Fh), &(ph2controls->zoomstack_Fh));
    rs_zoom_out(GWY_GRAPH(ph2controls->graph_H), &(ph2controls->zoomstack_H));
    update_models_ph2_Fh(ph2controls);
    update_models_ph2_H(ph2controls);
}

void ph2_zoom_restore_linked(Data *data)
{
    Ph2Controls *ph2controls;

    ph2controls = data->controls->ph2controls;

    rs_zoom_restore(GWY_GRAPH(ph2controls->graph_Fh), &(ph2controls->zoomstack_Fh));
    rs_zoom_restore(GWY_GRAPH(ph2controls->graph_H), &(ph2controls->zoomstack_H));
    update_models_ph2_Fh(ph2controls);
    update_models_ph2_H(ph2controls);
}

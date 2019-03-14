#include "tool_contact.h"

#include "datatypes.h"
#include "controls.h"

#include "fddata.h"
#include "instrument.h"
#include "gui-utils.h"
#include "niget-common.h"

#include "tool_op_gui.h"
#include "tool_op_odr_gui.h"
#include "tool_tangent_gui.h"
#include "tool_hertz_gui.h"
#include "tool_twoslopes_gui.h"
#include "dialog_area_calibration.h"

#include <math.h>

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>
#include <libprocess/gwyprocess.h>
#include <app/file.h>
#include <app/gwyapp.h>

const gchar *point_sel_status[4] = {"Not set", "Auto set", "Man set", "From file"};

static void plot_fd_curves(ContactControls *contactcontrols, FDdata *fddata);

static void radio_switch_to_fd_toggled(GtkWidget *widget, Data *data);
static void radio_switch_to_fd_time_toggled(GtkWidget *widget, Data *data);

static void update_gui_after_split(Data *data);
static gboolean one_point_selected(GwySelection *selection, GwyDataLine *x, gint *pindex);
static void button_contact_point_clicked(GtkWidget *widget, Data *data);
static void button_load_hold_clicked(GtkWidget *widget, Data *data);
static void button_hold_unload_clicked(GtkWidget *widget, Data *data);
static void button_end_unload_clicked(GtkWidget *widget, Data *data);
static void button_split_auto_clicked(GtkWidget *widget, Data *data);
static void button_save_data_clicked(GtkWidget *widget, Data *data);

static void save_preprocessed_data(gchar *filename, FDdata *fddata);

static void nu_changed(Data *data);
static void nui_changed(Data *data);
static void Ei_changed(Data *data);
static void reset_default(Data *data);
static void sigma_changed(Data *data);
//static void reset_delta(Data *data);

static void update_point_markline(GwyGraphModel *gmodel, GwyGraphCurveModel *cmodel,
                                  gdouble xmin, gdouble xmax, gdouble ymin, gdouble ymax, gint i, gchar *description);

static void contact_zoom_in_fd(Data *data);
static void contact_zoom_out_fd(Data *data);
/* static void contact_zoom_restore_fd(Data *data); */

static void contact_zoom_in_fd_time(Data *data);
static void contact_zoom_out_fd_time(Data *data);
/* static void contact_zoom_restore_fd_time(Data *data); */

gboolean lock_update = FALSE;


GtkWidget *tool_contact_create(Data *data)
{
    ContactControls *contactcontrols;

    GtkWidget *spin;

    gint    row;
    gchar   str[256];
    gchar   *str2;

    /* for better legibility */
    contactcontrols = data->controls->contactcontrols;

    /* Controls */
    contactcontrols->vbox_ctrls = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* Info section */

    contactcontrols->frame_info = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(contactcontrols->frame_info), GTK_WIDGET(gwy_label_new_header("Info")));
    contactcontrols->table_info = gtk_table_new(2, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(contactcontrols->table_info), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(contactcontrols->table_info), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(contactcontrols->table_info), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(contactcontrols->frame_info), contactcontrols->table_info);
    row = 0;

    /* h_max */
    gtk_table_attach(GTK_TABLE(contactcontrols->table_info), label_new_with_markup_left("h<sub>max</sub>"), 0, 1, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    contactcontrols->hmax = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_info), contactcontrols->hmax, 1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_misc_set_alignment(GTK_MISC(contactcontrols->hmax), 0.0, 0.5);

    gtk_table_attach(GTK_TABLE(contactcontrols->table_info), label_new_with_markup_left("nm"), 2, 3, row, row + 1, GTK_FILL, 0, 0, 0);

    row++;

    /* F_max */
    gtk_table_attach(GTK_TABLE(contactcontrols->table_info), label_new_with_markup_left("F<sub>max</sub>"), 0, 1, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    contactcontrols->Fmax = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_info), contactcontrols->Fmax, 1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_misc_set_alignment(GTK_MISC(contactcontrols->Fmax), 0.0, 0.5);

    gtk_table_attach(GTK_TABLE(contactcontrols->table_info), label_new_with_markup_left("mN"), 2, 3, row, row + 1, GTK_FILL, 0, 0, 0);

    /* switch graph mode */

    contactcontrols->hbox_switch = gtk_hbox_new(FALSE, CTRLS_SPACING);

    gtk_box_pack_start(GTK_BOX(contactcontrols->hbox_switch), gtk_label_new("Graph data as:"), FALSE, FALSE, 0);

    contactcontrols->radio_fd = gtk_radio_button_new_with_label(NULL, "F-d");
    g_signal_connect(contactcontrols->radio_fd, "toggled", G_CALLBACK(radio_switch_to_fd_toggled), data);
    gtk_box_pack_start(GTK_BOX(contactcontrols->hbox_switch), contactcontrols->radio_fd, TRUE, FALSE, 0);

    contactcontrols->radio_fd_time = gtk_radio_button_new_with_label_from_widget
                                     (GTK_RADIO_BUTTON(contactcontrols->radio_fd), "F(t), d(t)");
    g_signal_connect(contactcontrols->radio_fd_time, "toggled", G_CALLBACK(radio_switch_to_fd_time_toggled), data);
    gtk_box_pack_start(GTK_BOX(contactcontrols->hbox_switch), contactcontrols->radio_fd_time, TRUE, FALSE, 0);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(contactcontrols->radio_fd), TRUE);
    contactcontrols->fd_plot_mode = PLOT_FD;

    /* Contact & hold section */

    contactcontrols->frame_contact_hold = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(contactcontrols->frame_contact_hold), GTK_WIDGET(gwy_label_new_header("Contact-load-hold-unload")));
    contactcontrols->table_contact_hold = gtk_table_new(4, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(contactcontrols->table_contact_hold), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(contactcontrols->table_contact_hold), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(contactcontrols->table_contact_hold), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(contactcontrols->frame_contact_hold), contactcontrols->table_contact_hold);

    row = 0;

    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), label_new_left("Contact point"), 0, 1, row, row + 1, GTK_FILL, GTK_FILL, 0, 0);

    contactcontrols->status_contact_load = gtk_label_new(point_sel_status[POINT_STATUS_NOT_SET]);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), contactcontrols->status_contact_load, 1, 2, row, row + 1, GTK_FILL, GTK_FILL, 0, 0);

    contactcontrols->button_contact_point = gtk_button_new_with_mnemonic("Set _contact point");
    // gtk_button_set_image (GTK_BUTTON (contactcontrols->button_contact_point), gtk_image_new_from_stock (GWY_STOCK_FIX_ZERO, GTK_ICON_SIZE_BUTTON));
    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), contactcontrols->button_contact_point, 2, 3, row, row + 1, GTK_FILL | GTK_EXPAND, GTK_FILL | GTK_EXPAND, 0, 0);
    gtk_widget_set_sensitive(contactcontrols->button_contact_point, FALSE);
    g_signal_connect(contactcontrols->button_contact_point, "clicked", G_CALLBACK(button_contact_point_clicked), data);

    row++;

    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), label_new_left("Hold begin"), 0, 1, row, row + 1, GTK_FILL, GTK_FILL, 0, 0);

    contactcontrols->status_load_hold = gtk_label_new(point_sel_status[POINT_STATUS_NOT_SET]);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), contactcontrols->status_load_hold, 1, 2, row, row + 1, GTK_FILL, GTK_FILL, 0, 0);

    contactcontrols->button_load_hold = gtk_button_new_with_mnemonic("Set _hold begin");
    // gtk_button_set_image (GTK_BUTTON (contactcontrols->button_load_hold), gtk_image_new_from_stock (GWY_STOCK_FIX_ZERO, GTK_ICON_SIZE_BUTTON));
    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), contactcontrols->button_load_hold, 2, 3, row, row + 1, GTK_FILL | GTK_EXPAND, GTK_FILL | GTK_EXPAND, 0, 0);
    gtk_widget_set_sensitive(contactcontrols->button_load_hold, FALSE);
    g_signal_connect(contactcontrols->button_load_hold, "clicked", G_CALLBACK(button_load_hold_clicked), data);

    row++;

    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), label_new_left("Unload begin"), 0, 1, row, row + 1, GTK_FILL, GTK_FILL, 0, 0);

    contactcontrols->status_hold_unload = gtk_label_new(point_sel_status[POINT_STATUS_NOT_SET]);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), contactcontrols->status_hold_unload, 1, 2, row, row + 1, GTK_FILL, GTK_FILL, 0, 0);

    contactcontrols->button_hold_unload = gtk_button_new_with_mnemonic("Set _unload begin");
    // gtk_button_set_image (GTK_BUTTON (contactcontrols->button_hold_unload), gtk_image_new_from_stock (GWY_STOCK_FIX_ZERO, GTK_ICON_SIZE_BUTTON));
    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), contactcontrols->button_hold_unload, 2, 3, row, row + 1, GTK_FILL | GTK_EXPAND, GTK_FILL | GTK_EXPAND, 0, 0);
    gtk_widget_set_sensitive(contactcontrols->button_hold_unload, FALSE);
    g_signal_connect(contactcontrols->button_hold_unload, "clicked", G_CALLBACK(button_hold_unload_clicked), data);

    row++;

    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), label_new_left("Unload end"), 0, 1, row, row + 1, GTK_FILL, GTK_FILL, 0, 0);

    contactcontrols->status_end_unload = gtk_label_new(point_sel_status[POINT_STATUS_NOT_SET]);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), contactcontrols->status_end_unload, 1, 2, row, row + 1, GTK_FILL, GTK_FILL, 0, 0);

    contactcontrols->button_end_unload = gtk_button_new_with_mnemonic("Set unload _end");
    // gtk_button_set_image (GTK_BUTTON (contactcontrols->button_hold_unload), gtk_image_new_from_stock (GWY_STOCK_FIX_ZERO, GTK_ICON_SIZE_BUTTON));
    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), contactcontrols->button_end_unload, 2, 3, row, row + 1, GTK_FILL | GTK_EXPAND, GTK_FILL | GTK_EXPAND, 0, 0);
    gtk_widget_set_sensitive(contactcontrols->button_end_unload, FALSE);
    g_signal_connect(contactcontrols->button_end_unload, "clicked", G_CALLBACK(button_end_unload_clicked), data);

    row++;

    contactcontrols->button_split_auto = gtk_button_new_with_mnemonic("Reset to automatic split");
    gtk_table_attach(GTK_TABLE(contactcontrols->table_contact_hold), contactcontrols->button_split_auto, 0, 3, row, row + 1, GTK_FILL | GTK_EXPAND, GTK_FILL | GTK_EXPAND, 0, 0);
    g_signal_connect(contactcontrols->button_split_auto, "clicked", G_CALLBACK(button_split_auto_clicked), data);


    /* Save data button */

    contactcontrols->button_save_data = gtk_button_new_with_mnemonic("Save preprocessed data");
    gtk_button_set_image(GTK_BUTTON(contactcontrols->button_save_data), gtk_image_new_from_stock(GTK_STOCK_SAVE, GTK_ICON_SIZE_BUTTON));
    gtk_widget_set_sensitive(contactcontrols->button_save_data, FALSE);
    g_signal_connect(contactcontrols->button_save_data, "clicked", G_CALLBACK(button_save_data_clicked), data);

    /* Area section */

    contactcontrols->frame_area = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(contactcontrols->frame_area), GTK_WIDGET(gwy_label_new_header("Area")));
    contactcontrols->table_area = gtk_table_new(2, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(contactcontrols->table_area), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(contactcontrols->table_area), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(contactcontrols->table_area), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(contactcontrols->frame_area), contactcontrols->table_area);
    row = 0;

    /* area function button */
    contactcontrols->area_function_button = gtk_button_new_with_mnemonic("Area _function");
    gtk_button_set_image(GTK_BUTTON(contactcontrols->area_function_button), gtk_image_new_from_stock(GWY_STOCK_DATA_MEASURE, GTK_ICON_SIZE_BUTTON));
    gtk_table_attach(GTK_TABLE(contactcontrols->table_area), contactcontrols->area_function_button, 0, 3, row, row + 1, GTK_EXPAND | GTK_FILL, GTK_FILL, 0, 0);
    g_signal_connect(contactcontrols->area_function_button, "clicked", G_CALLBACK(area_calibration), data);

    row++;

    /* area function label */
    gtk_table_attach(GTK_TABLE(contactcontrols->table_area), label_new_with_markup_left("Area function"), 0, 1, row, row + 1, GTK_FILL, 0, 0, 0);

    row++;

    str2 = write_area_label(&(data->args->area));
    contactcontrols->area_function_label = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(contactcontrols->area_function_label), str2);
    g_free(str2);
    gtk_label_set_line_wrap(GTK_LABEL(contactcontrols->area_function_label), TRUE);
    gtk_label_set_line_wrap_mode(GTK_LABEL(contactcontrols->area_function_label), PANGO_WRAP_WORD_CHAR);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_area), contactcontrols->area_function_label, 0, 3, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);


    /* Sample & indenter section */

    contactcontrols->frame_sample_indenter = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(contactcontrols->frame_sample_indenter), GTK_WIDGET(gwy_label_new_header("Sample &amp; indenter")));
    contactcontrols->table_sample_indenter = gtk_table_new(4, 4, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(contactcontrols->table_sample_indenter), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(contactcontrols->table_sample_indenter), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(contactcontrols->table_sample_indenter), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(contactcontrols->frame_sample_indenter), contactcontrols->table_sample_indenter);
    row = 0;

    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), label_new_with_markup_left("ν"), 0, 1, row, row + 1, 0, 0, 0, 0);
    contactcontrols->nu = gtk_adjustment_new(data->args->instdata.nu, -1.0, 0.5, 0.01, 0.1, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(contactcontrols->nu), 0.01, 2);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), spin, 1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(contactcontrols->nu, "value-changed", G_CALLBACK(nu_changed), data);
    row++;

    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), label_new_with_markup_left("E<sub>i</sub>"), 0, 1, row, row + 1, 0, 0, 0, 0);
    contactcontrols->Ei = gtk_adjustment_new(data->args->instdata.Ei * 1e-9, 0, 9e4, 10, 100, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(contactcontrols->Ei), 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), spin, 1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(contactcontrols->Ei, "value-changed", G_CALLBACK(Ei_changed), data);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), label_new_with_markup_left("GPa"), 2, 3, row, row + 1, 0, 0, 0, 0);
    row++;

    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), label_new_with_markup_left("ν<sub>i</sub>"), 0, 1, row, row + 1, 0, 0, 0, 0);
    contactcontrols->nui = gtk_adjustment_new(data->args->instdata.nui, -1.0, 0.5, 0.01, 0.1, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(contactcontrols->nui), 0.01, 2);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), spin, 1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(contactcontrols->nui, "value-changed", G_CALLBACK(nui_changed), data);
    row++;

    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), label_new_with_markup_left("σ<sub>h</sub>"), 0, 1, row, row + 1, 0, 0, 0, 0);
    contactcontrols->sigma_h = gtk_adjustment_new(data->args->instdata.sigma_h, 0.0, 1e3, 0.01, 0.1, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(contactcontrols->sigma_h), 0.001, 3);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), spin, 1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(contactcontrols->sigma_h, "value-changed", G_CALLBACK(sigma_changed), data);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), label_new_with_markup_left("nm"), 2, 3, row, row + 1, 0, 0, 0, 0);
    row++;

    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), label_new_with_markup_left("σ<sub>F</sub>"), 0, 1, row, row + 1, 0, 0, 0, 0);
    contactcontrols->sigma_F = gtk_adjustment_new(data->args->instdata.sigma_F, 0.0, 1e3, 0.01, 0.1, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(contactcontrols->sigma_F), 0.000001, 6);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), spin, 1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(contactcontrols->sigma_F, "value-changed", G_CALLBACK(sigma_changed), data);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), label_new_with_markup_left("mN"), 2, 3, row, row + 1, 0, 0, 0, 0);
    row++;

    //gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), label_new_with_markup_left("δ = σ<sub>F</sub><sup>2</sup>/σ<sub>h</sub><sup>2</sup>"), 0, 2, row, row+1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), label_new_with_markup_left("δ"), 0, 1, row, row + 1, 0, 0, 0, 0);
    contactcontrols->delta = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(contactcontrols->delta), 0.0, 0.5);

    g_snprintf(str, sizeof(str), DELTA, data->args->instdata.delta);
    gtk_label_set_markup(GTK_LABEL(data->controls->contactcontrols->delta), str);

    //label_set_gdouble_format(contactcontrols->delta, data->args->instdata.delta, DELTA);
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), contactcontrols->delta, 1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    row++;

    contactcontrols->button_reset =  gtk_button_new_with_label("Reset to defaults");
    gtk_button_set_image(GTK_BUTTON(contactcontrols->button_reset), gtk_image_new_from_stock(GTK_STOCK_REVERT_TO_SAVED, GTK_ICON_SIZE_BUTTON));
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), contactcontrols->button_reset, 0, 4, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(contactcontrols->button_reset, "clicked", G_CALLBACK(reset_default), data);

    /*
    contactcontrols->reset_delta = gtk_button_new();
    gtk_button_set_image (GTK_BUTTON(contactcontrols->reset_delta), gtk_image_new_from_stock (GTK_STOCK_REVERT_TO_SAVED, GTK_ICON_SIZE_BUTTON));
    gtk_table_attach(GTK_TABLE(contactcontrols->table_sample_indenter), contactcontrols->reset_delta, 3,4,i,i+1, GTK_FILL, 0,0,0);
    g_signal_connect_swapped (contactcontrols->reset_delta, "clicked", G_CALLBACK (reset_delta), data);
    i++;
    */

    /* ======== */

    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_ctrls), contactcontrols->frame_info, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_ctrls), contactcontrols->hbox_switch, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_ctrls), contactcontrols->frame_contact_hold, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_ctrls), contactcontrols->button_save_data, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_ctrls), contactcontrols->frame_area, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_ctrls), contactcontrols->frame_sample_indenter, FALSE, FALSE, 0);

    /* ======== */

    contactcontrols->hbox_graph = gtk_hbox_new(FALSE, CTRLS_SPACING);

    /* Graph model FD */

    contactcontrols->graph_model_fd = gwy_graph_model_new();
    g_object_set(contactcontrols->graph_model_fd, "title", "Load-Hold-Unload", NULL);
    gwy_graph_model_set_axis_label(contactcontrols->graph_model_fd, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(contactcontrols->graph_model_fd, GTK_POS_LEFT, "F / mN");

    contactcontrols->zoomstack_fd = NULL;
    contactcontrols->zoomstack_fd = g_slist_prepend(contactcontrols->zoomstack_fd, contactcontrols->graph_model_fd);

    contactcontrols->cmodelload = gwy_graph_curve_model_new();
    contactcontrols->cmodelhold = gwy_graph_curve_model_new();
    contactcontrols->cmodelunload = gwy_graph_curve_model_new();

    g_object_set(contactcontrols->cmodelload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_LOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 "description", "Loading",
                 NULL);
    g_object_set(contactcontrols->cmodelhold,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_HOLD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 "description", "Hold",
                 NULL);
    g_object_set(contactcontrols->cmodelunload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_UNLOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 "description", "Unloading",
                 NULL);

    gwy_graph_model_add_curve(contactcontrols->graph_model_fd, contactcontrols->cmodelload);
    gwy_graph_model_add_curve(contactcontrols->graph_model_fd, contactcontrols->cmodelhold);
    gwy_graph_model_add_curve(contactcontrols->graph_model_fd, contactcontrols->cmodelunload);

    /* Graphs */

    /* Graph F-d */

    contactcontrols->graph_fd = gwy_graph_new(contactcontrols->graph_model_fd); /* start with fd curves */
    gtk_widget_set_size_request(contactcontrols->graph_fd, 600, 400);

    gwy_graph_enable_user_input(GWY_GRAPH(contactcontrols->graph_fd), FALSE);
    gwy_graph_set_status(GWY_GRAPH(contactcontrols->graph_fd), GWY_GRAPH_STATUS_XSEL); /* use selection as interval for selecting contact point */

    contactcontrols->selection_fd =
        gwy_graph_area_get_selection(GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(contactcontrols->graph_fd))), GWY_GRAPH_STATUS_XSEL);
    gwy_selection_set_max_objects(contactcontrols->selection_fd, 1);

    /* Zooming F-d */

    contactcontrols->vbox_zoom_buttons_fd = gtk_vbox_new(TRUE, 0);
    contactcontrols->vbox_zoom_outer_fd = gtk_vbox_new(FALSE, 0);

    contactcontrols->button_zoom_in_fd = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(contactcontrols->button_zoom_in_fd), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_zoom_buttons_fd), contactcontrols->button_zoom_in_fd, TRUE, TRUE, 0);
    g_signal_connect_swapped(contactcontrols->button_zoom_in_fd, "clicked", G_CALLBACK(contact_zoom_in_fd), data);

    contactcontrols->button_zoom_out_fd = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(contactcontrols->button_zoom_out_fd), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_zoom_buttons_fd), contactcontrols->button_zoom_out_fd, TRUE, TRUE, 0);
    g_signal_connect_swapped(contactcontrols->button_zoom_out_fd, "clicked", G_CALLBACK(contact_zoom_out_fd), data);

    contactcontrols->button_zoom_restore_fd = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(contactcontrols->button_zoom_restore_fd), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_zoom_buttons_fd), contactcontrols->button_zoom_restore_fd, TRUE, TRUE, 0);
    g_signal_connect_swapped(contactcontrols->button_zoom_restore_fd, "clicked", G_CALLBACK(contact_zoom_restore_fd), data);

    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_zoom_outer_fd), contactcontrols->vbox_zoom_buttons_fd, TRUE, FALSE, 0);

    /* Graph model F,d in time */

    contactcontrols->graph_model_fd_time = gwy_graph_model_new();
    g_object_set(contactcontrols->graph_model_fd_time, "title", "Force and depth in time", NULL);
    gwy_graph_model_set_axis_label(contactcontrols->graph_model_fd_time, GTK_POS_BOTTOM, "Step");
    /* gwy_graph_model_set_axis_label (contactcontrols->graph_model_fd_time, GTK_POS_LEFT, "F / mN"); */
    /* gwy_graph_model_set_axis_label (contactcontrols->graph_model_fd_time, GTK_POS_RIGHT, "h / nm"); */

    contactcontrols->zoomstack_fd_time = NULL;
    contactcontrols->zoomstack_fd_time = g_slist_prepend(contactcontrols->zoomstack_fd_time, contactcontrols->graph_model_fd_time);

    contactcontrols->cmodelFtime = gwy_graph_curve_model_new();
    contactcontrols->cmodelhtime = gwy_graph_curve_model_new();
    contactcontrols->cmodelFc = gwy_graph_curve_model_new();

    contactcontrols->cmodel_contact_load = gwy_graph_curve_model_new();
    contactcontrols->cmodel_load_hold = gwy_graph_curve_model_new();
    contactcontrols->cmodel_hold_unload = gwy_graph_curve_model_new();
    contactcontrols->cmodel_end_unload = gwy_graph_curve_model_new();

    g_object_set(contactcontrols->cmodelFtime,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "point-size", 0,
                 "color", CURVE_COLOR_F,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 "description", "F",
                 NULL);

    g_object_set(contactcontrols->cmodelhtime,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "point-size", 0,
                 "color", CURVE_COLOR_H,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 "description", "h",
                 NULL);

    g_object_set(contactcontrols->cmodelFc,
                 "mode", GWY_GRAPH_CURVE_LINE,
                 "color", CURVE_COLOR_FC,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 1,
                 "description", "Fc",
                 NULL);

    g_object_set(contactcontrols->cmodel_contact_load,
                 "mode", GWY_GRAPH_CURVE_LINE,
                 "color", CURVE_COLOR_LOAD,
                 "line-style", GDK_LINE_ON_OFF_DASH,
                 "line-width", 2,
                 "description", "Contact point",
                 NULL);

    g_object_set(contactcontrols->cmodel_load_hold,
                 "mode", GWY_GRAPH_CURVE_LINE,
                 "color", CURVE_COLOR_HOLD,
                 "line-style", GDK_LINE_ON_OFF_DASH,
                 "line-width", 2,
                 "description", "Loading/hold",
                 NULL);

    g_object_set(contactcontrols->cmodel_hold_unload,
                 "mode", GWY_GRAPH_CURVE_LINE,
                 "color", CURVE_COLOR_UNLOAD,
                 "line-style", GDK_LINE_ON_OFF_DASH,
                 "line-width", 2,
                 "description", "Hold/unloading",
                 NULL);

    g_object_set(contactcontrols->cmodel_end_unload,
                 "mode", GWY_GRAPH_CURVE_LINE,
                 "color", CURVE_COLOR_END_UNLOAD,
                 "line-style", GDK_LINE_ON_OFF_DASH,
                 "line-width", 2,
                 "description", "Unloading end",
                 NULL);

    gwy_graph_model_add_curve(contactcontrols->graph_model_fd_time, contactcontrols->cmodelFtime);
    gwy_graph_model_add_curve(contactcontrols->graph_model_fd_time, contactcontrols->cmodelhtime);
    gwy_graph_model_add_curve(contactcontrols->graph_model_fd_time, contactcontrols->cmodelFc);

    gwy_graph_model_add_curve(contactcontrols->graph_model_fd_time, contactcontrols->cmodel_contact_load);
    gwy_graph_model_add_curve(contactcontrols->graph_model_fd_time, contactcontrols->cmodel_load_hold);
    gwy_graph_model_add_curve(contactcontrols->graph_model_fd_time, contactcontrols->cmodel_hold_unload);
    gwy_graph_model_add_curve(contactcontrols->graph_model_fd_time, contactcontrols->cmodel_end_unload);

    /* Graph F-d-time */

    contactcontrols->graph_fd_time = gwy_graph_new(contactcontrols->graph_model_fd_time); /* start with fd curves */
    gtk_widget_set_size_request(contactcontrols->graph_fd_time, 600, 400);

    gwy_graph_set_axis_visible(GWY_GRAPH(contactcontrols->graph_fd_time), GTK_POS_LEFT, FALSE);
    gwy_graph_enable_user_input(GWY_GRAPH(contactcontrols->graph_fd_time), FALSE);
    gwy_graph_set_status(GWY_GRAPH(contactcontrols->graph_fd_time), GWY_GRAPH_STATUS_XSEL);  /* use selection as interval for selecting contact point */

    contactcontrols->selection_fd_time =
        gwy_graph_area_get_selection(GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(contactcontrols->graph_fd_time))), GWY_GRAPH_STATUS_XSEL);
    gwy_selection_set_max_objects(contactcontrols->selection_fd_time, 1);

    /* Zooming F-d-time */

    contactcontrols->vbox_zoom_buttons_fd_time = gtk_vbox_new(TRUE, 0);
    contactcontrols->vbox_zoom_outer_fd_time = gtk_vbox_new(FALSE, 0);

    contactcontrols->button_zoom_in_fd_time = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(contactcontrols->button_zoom_in_fd_time), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_zoom_buttons_fd_time), contactcontrols->button_zoom_in_fd_time, TRUE, TRUE, 0);
    g_signal_connect_swapped(contactcontrols->button_zoom_in_fd_time, "clicked", G_CALLBACK(contact_zoom_in_fd_time), data);

    contactcontrols->button_zoom_out_fd_time = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(contactcontrols->button_zoom_out_fd_time), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_zoom_buttons_fd_time), contactcontrols->button_zoom_out_fd_time, TRUE, TRUE, 0);
    g_signal_connect_swapped(contactcontrols->button_zoom_out_fd_time, "clicked", G_CALLBACK(contact_zoom_out_fd_time), data);

    contactcontrols->button_zoom_restore_fd_time = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(contactcontrols->button_zoom_restore_fd_time), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_zoom_buttons_fd_time), contactcontrols->button_zoom_restore_fd_time, TRUE, TRUE, 0);
    g_signal_connect_swapped(contactcontrols->button_zoom_restore_fd_time, "clicked", G_CALLBACK(contact_zoom_restore_fd_time), data);

    gtk_box_pack_start(GTK_BOX(contactcontrols->vbox_zoom_outer_fd_time), contactcontrols->vbox_zoom_buttons_fd_time, TRUE, FALSE, 0);

    /* Graph hbox */

    gtk_box_pack_start(GTK_BOX(contactcontrols->hbox_graph), contactcontrols->graph_fd, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(contactcontrols->hbox_graph), contactcontrols->vbox_zoom_outer_fd, FALSE, FALSE, 0);

    /* ================================ */

    /* Main */

    contactcontrols->contact_gui = gtk_hbox_new(FALSE, 10);  /* homogeneous, spacing */
    gtk_box_pack_start(GTK_BOX(contactcontrols->contact_gui), contactcontrols->vbox_ctrls, FALSE, FALSE, 0);   /* expand, fill, padding */
    gtk_box_pack_start(GTK_BOX(contactcontrols->contact_gui), contactcontrols->hbox_graph, TRUE, TRUE, 0);
    //contactcontrols->tool_label = gtk_label_new ("Contact point & area");

    return contactcontrols->contact_gui;
}

void contact_switch_to_fd(Data *data)
{
    ContactControls *contactcontrols;

    contactcontrols = data->controls->contactcontrols;

    if (contactcontrols->fd_plot_mode == PLOT_FD_TIME) {

        contactcontrols->fd_plot_mode = PLOT_FD;

        /* inc. ref. count to protect widget from being destroyed */
        g_object_ref(contactcontrols->graph_fd_time);
        g_object_ref(contactcontrols->vbox_zoom_outer_fd_time);

        gtk_container_remove(GTK_CONTAINER(contactcontrols->hbox_graph), contactcontrols->graph_fd_time);
        gtk_container_remove(GTK_CONTAINER(contactcontrols->hbox_graph), contactcontrols->vbox_zoom_outer_fd_time);

        gtk_box_pack_start(GTK_BOX(contactcontrols->hbox_graph), contactcontrols->graph_fd, TRUE, TRUE, 0);
        gtk_box_pack_start(GTK_BOX(contactcontrols->hbox_graph), contactcontrols->vbox_zoom_outer_fd, FALSE, FALSE, 0);

        gtk_widget_show(contactcontrols->graph_fd);
        gtk_widget_show_all(contactcontrols->vbox_zoom_outer_fd);

        gtk_widget_set_sensitive(contactcontrols->button_contact_point, FALSE);
        gtk_widget_set_sensitive(contactcontrols->button_load_hold, FALSE);
        gtk_widget_set_sensitive(contactcontrols->button_hold_unload, FALSE);
        gtk_widget_set_sensitive(contactcontrols->button_end_unload, FALSE);
    }
}

static void radio_switch_to_fd_toggled(GtkWidget *widget, Data *data)
{
    contact_switch_to_fd(data);
}

void contact_switch_to_fd_time(Data *data)
{
    ContactControls *contactcontrols;
    FDdata *fddata;

    contactcontrols = data->controls->contactcontrols;
    fddata = &(data->args->fddata);

    if (contactcontrols->fd_plot_mode == PLOT_FD) {

        contactcontrols->fd_plot_mode = PLOT_FD_TIME;

        update_point_marklines(contactcontrols, fddata);

        /* inc. ref. counts to protect widgets from being destroyed */
        g_object_ref(contactcontrols->graph_fd);
        g_object_ref(contactcontrols->vbox_zoom_outer_fd);

        gtk_container_remove(GTK_CONTAINER(contactcontrols->hbox_graph), contactcontrols->graph_fd);
        gtk_container_remove(GTK_CONTAINER(contactcontrols->hbox_graph), contactcontrols->vbox_zoom_outer_fd);

        gtk_box_pack_start(GTK_BOX(contactcontrols->hbox_graph), contactcontrols->graph_fd_time, TRUE, TRUE, 0);
        gtk_box_pack_start(GTK_BOX(contactcontrols->hbox_graph), contactcontrols->vbox_zoom_outer_fd_time, FALSE, FALSE, 0);

        gtk_widget_show(contactcontrols->graph_fd_time);
        gtk_widget_show_all(contactcontrols->vbox_zoom_outer_fd_time);

        gtk_widget_set_sensitive(contactcontrols->button_contact_point, TRUE);
        gtk_widget_set_sensitive(contactcontrols->button_load_hold, (fddata->status_contact_load != POINT_STATUS_NOT_SET));
        gtk_widget_set_sensitive(contactcontrols->button_hold_unload, (fddata->status_load_hold != POINT_STATUS_NOT_SET));
        gtk_widget_set_sensitive(contactcontrols->button_end_unload, (fddata->status_hold_unload != POINT_STATUS_NOT_SET));
    }
}

static void radio_switch_to_fd_time_toggled(GtkWidget *widget, Data *data)
{
    contact_switch_to_fd_time(data);
}

static gboolean one_point_selected(GwySelection *selection, GwyDataLine *x, gint *pindex)
{
    gint istart, iend, ndata, nselections;
    gdouble from, to;
    gdouble range[2];

    nselections = gwy_selection_get_data(selection, NULL);

    if (nselections > 0) { /* that means 1 */

        gwy_selection_get_object(selection, 0, range);
        from = MIN(range[0], range[1]);
        to = MAX(range[0], range[1]);

        range_to_indices_dataline(from, to, x, FALSE, &istart, &iend, &ndata);

        if (ndata == 0) {
            show_error("No point selected. Please select a wider interval that contains just one point.");
            return FALSE;
        }
        else if (ndata >= 2) {
            show_error("Too many points selected. Please select a shorter interval that contains just one point.");
            return FALSE;
        }
        else { /* ndata = 1 */
            *pindex = istart;
            return TRUE;
        }
    }
    else {
        return FALSE;
    }
}

enum PointType {
    POINT_CONTACT_LOAD = 0,
    POINT_LOAD_HOLD,
    POINT_HOLD_UNLOAD,
    POINT_END_UNLOAD
};

static void update_gui_after_split(Data *data)
{
    FDdata *fddata;
    ContactControls *contactcontrols;

    fddata = &(data->args->fddata);
    contactcontrols = data->controls->contactcontrols;

    update_max_labels(contactcontrols, fddata);

    gwy_graph_model_remove_all_curves(contactcontrols->graph_model_fd);
    plot_fd_curves(contactcontrols, fddata);
    update_point_marklines(contactcontrols, fddata);

    data->functions->remove_all_tool_fit_results_labels(data, FALSE);
    data->functions->redraw_all_tools(data);
}

static void button_contact_point_clicked(GtkWidget *widget, Data *data)
{
    FDdata *fddata;
    ContactControls *contactcontrols;

    gboolean splitok;
    gint point;

    fddata = &(data->args->fddata);
    contactcontrols = data->controls->contactcontrols;

    if (one_point_selected(contactcontrols->selection_fd_time, fddata->time, &point)) {

        fddata->status_contact_load = POINT_STATUS_SET_MAN;
        gtk_label_set_text(GTK_LABEL(contactcontrols->status_contact_load), point_sel_status[POINT_STATUS_SET_MAN]);

        if (point != fddata->i_contact_load) { /* selected different point from previous */

            fddata->i_contact_load = point;
            fddata->Fc = fddata->Forig->data[point];

            splitok = split_hF_data_auto(fddata);

            if (splitok) {
                gtk_label_set_text(GTK_LABEL(contactcontrols->status_load_hold), point_sel_status[POINT_STATUS_SET_AUTO]);
                gtk_label_set_text(GTK_LABEL(contactcontrols->status_hold_unload), point_sel_status[POINT_STATUS_SET_AUTO]);
                gtk_label_set_text(GTK_LABEL(contactcontrols->status_end_unload), point_sel_status[POINT_STATUS_SET_AUTO]);
            }
            else {
                show_warning("Could not split data into loading and unloading curves. Check your data.");
            }

            update_gui_after_split(data);

            gtk_widget_set_sensitive(contactcontrols->button_load_hold, TRUE);
            gtk_widget_set_sensitive(contactcontrols->button_hold_unload, TRUE);
            gtk_widget_set_sensitive(contactcontrols->button_end_unload, TRUE);
        } /* selected different point from previous */
    } /* one point selected */
}

static void button_load_hold_clicked(GtkWidget *widget, Data *data)
{
    FDdata *fddata;
    ContactControls *contactcontrols;
    gint point;

    fddata = &(data->args->fddata);
    contactcontrols = data->controls->contactcontrols;

    if (one_point_selected(contactcontrols->selection_fd_time, fddata->time, &point)) {
        if (point >= fddata->i_contact_load) {

            fddata->status_load_hold = POINT_STATUS_SET_MAN;
            gtk_label_set_text(GTK_LABEL(contactcontrols->status_load_hold), point_sel_status[POINT_STATUS_SET_MAN]);

            if (point != fddata->i_load_hold) { /* selected different point from previous */

                fddata->i_load_hold = point;

                if (point > fddata->i_hold_unload) {
                    fddata->i_hold_unload = point;
                    fddata->status_hold_unload = POINT_STATUS_SET_MAN;
                    gtk_label_set_text(GTK_LABEL(contactcontrols->status_hold_unload), point_sel_status[POINT_STATUS_SET_MAN]);
                }

                if (point > fddata->i_end_unload) {
                    fddata->i_end_unload = point;
                    fddata->status_end_unload = POINT_STATUS_SET_MAN;
                    gtk_label_set_text(GTK_LABEL(contactcontrols->status_end_unload), point_sel_status[POINT_STATUS_SET_MAN]);
                }

                split_hF_data_spec(fddata, fddata->i_contact_load, fddata->i_load_hold, fddata->i_hold_unload, fddata->i_end_unload);
                update_gui_after_split(data);
            } /* selected different point from previous */
        }
        else {
            show_warning("The hold-begin point cannot preceed the contact point.");
        }
    } /* one point selected */
}

static void button_hold_unload_clicked(GtkWidget *widget, Data *data)
{
    FDdata *fddata;
    ContactControls *contactcontrols;
    gint point;

    fddata = &(data->args->fddata);
    contactcontrols = data->controls->contactcontrols;

    if (one_point_selected(contactcontrols->selection_fd_time, fddata->time, &point)) {
        if (point >= fddata->i_contact_load) {

            fddata->status_hold_unload = POINT_STATUS_SET_MAN;
            gtk_label_set_text(GTK_LABEL(contactcontrols->status_hold_unload), point_sel_status[POINT_STATUS_SET_MAN]);

            if (point != fddata->i_hold_unload) { /* selected different point from previous */

                fddata->i_hold_unload = point;

                if (point < fddata->i_load_hold) {
                    fddata->i_load_hold = point;
                    fddata->status_load_hold = POINT_STATUS_SET_MAN;
                    gtk_label_set_text(GTK_LABEL(contactcontrols->status_load_hold), point_sel_status[POINT_STATUS_SET_MAN]);
                }

                if (point > fddata->i_end_unload) {
                    fddata->i_end_unload = point;
                    fddata->status_end_unload = POINT_STATUS_SET_MAN;
                    gtk_label_set_text(GTK_LABEL(contactcontrols->status_end_unload), point_sel_status[POINT_STATUS_SET_MAN]);
                }

                split_hF_data_spec(fddata, fddata->i_contact_load, fddata->i_load_hold, fddata->i_hold_unload, fddata->i_end_unload);
                update_gui_after_split(data);
            } /* selected different point from previous */
        }
        else {
            show_warning("The hold-end point cannot preceed the contact point.");
        }
    } /* one point selected */
}

static void button_end_unload_clicked(GtkWidget *widget, Data *data)
{
    FDdata *fddata;
    ContactControls *contactcontrols;
    gint point;

    fddata = &(data->args->fddata);
    contactcontrols = data->controls->contactcontrols;

    if (one_point_selected(contactcontrols->selection_fd_time, fddata->time, &point)) {
        if (point >= fddata->i_contact_load) {

            fddata->status_end_unload = POINT_STATUS_SET_MAN;
            gtk_label_set_text(GTK_LABEL(contactcontrols->status_end_unload), point_sel_status[POINT_STATUS_SET_MAN]);

            if (point != fddata->i_end_unload) { /* selected different point from previous */

                fddata->i_end_unload = point;

                if (point < fddata->i_hold_unload) {
                    fddata->i_hold_unload = point;
                    fddata->status_hold_unload = POINT_STATUS_SET_MAN;
                    gtk_label_set_text(GTK_LABEL(contactcontrols->status_hold_unload), point_sel_status[POINT_STATUS_SET_MAN]);
                }

                if (point < fddata->i_load_hold) {
                    fddata->i_load_hold = point;
                    fddata->status_load_hold = POINT_STATUS_SET_MAN;
                    gtk_label_set_text(GTK_LABEL(contactcontrols->status_load_hold), point_sel_status[POINT_STATUS_SET_MAN]);
                }

                split_hF_data_spec(fddata, fddata->i_contact_load, fddata->i_load_hold, fddata->i_hold_unload, fddata->i_end_unload);
                update_gui_after_split(data);
            } /* selected different point from previous */
        }
        else {
            show_warning("The unload-end point cannot preceed the contact point.");
        }
    } /* one point selected */
}

static void button_split_auto_clicked(GtkWidget *widget, Data *data)
{
    ContactControls *contactcontrols;
    gboolean splitok;

    contactcontrols = data->controls->contactcontrols;

    splitok = split_hF_data_auto(&(data->args->fddata));

    if (!splitok) {
        show_warning("Could not split data into loading and unloading curves. Check your data.");
    }

    update_fd_plot(data->controls->contactcontrols, &(data->args->fddata));
    update_fd_time_plot(data->controls->contactcontrols, &(data->args->fddata));

    gtk_label_set_text(GTK_LABEL(contactcontrols->status_load_hold), point_sel_status[POINT_STATUS_SET_AUTO]);
    gtk_label_set_text(GTK_LABEL(contactcontrols->status_hold_unload), point_sel_status[POINT_STATUS_SET_AUTO]);
    gtk_label_set_text(GTK_LABEL(contactcontrols->status_end_unload), point_sel_status[POINT_STATUS_SET_AUTO]);

    update_gui_after_split(data);
}

static void button_save_data_clicked(GtkWidget *widget, Data *data)
{
    GtkWidget *dialog;
    gchar *filename;
    gchar *basename;
    gchar *ext;
    gchar *premise_filename = NULL;

    basename = g_path_get_basename(data->args->fddata.filename);
    ext = g_utf8_strrchr(basename, g_utf8_strlen(basename, -1), '.');

    if (NULL != ext) {
        *ext = '\0';
    }

    premise_filename = g_strconcat(basename, ".nig", NULL);
    g_free(basename);

    dialog = gtk_file_chooser_dialog_new("Save preprocessed data", NULL,
                                         GTK_FILE_CHOOSER_ACTION_SAVE,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT, NULL);
    gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), gwy_app_get_current_directory());
    gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog), TRUE);
    gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(dialog), premise_filename);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
        save_preprocessed_data(filename, &(data->args->fddata));
        g_free(filename);

        if (NULL != premise_filename) {
            g_free(premise_filename);
        }
    }

    gtk_widget_destroy(dialog);
}

static void save_preprocessed_data(gchar *filename, FDdata *fddata)
{
    FILE *file;
    gint i;

    file = fopen(filename, "w");

    if (!file) {
        show_error("Could not save file.");
        return;
    }

    fprintf(file, DATAFORMAT_HEADER_FORMAT_NEWLINE,
            fddata->i_contact_load, fddata->i_load_hold, fddata->i_hold_unload, fddata->i_end_unload);

    for (i = 0; i < fddata->ndata; i++) {
        fprintf(file, "%g\t%g\n", fddata->horig->data[i], fddata->Forig->data[i]);
    }

    fclose(file);
}

static void nu_changed(Data *data)
{
    data->args->instdata.nu = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->contactcontrols->nu));

    op_recalculate(data);
    tangent_recalculate(data);
    hertz_recalculate(data);
#ifndef NOFORTRAN
    op_odr_recalculate(data);
    slopes_recalculate(data);
#endif
}

static void nui_changed(Data *data)
{
    data->args->instdata.nui = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->contactcontrols->nui));

    op_recalculate(data);
    tangent_recalculate(data);
    hertz_recalculate(data);
#ifndef NOFORTRAN
    op_odr_recalculate(data);
    slopes_recalculate(data);
#endif
}

static void Ei_changed(Data *data)
{
    data->args->instdata.Ei = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->contactcontrols->Ei)) * 1e9;

    op_recalculate(data);
    tangent_recalculate(data);
    hertz_recalculate(data);
#ifndef NOFORTRAN
    op_odr_recalculate(data);
    slopes_recalculate(data);
#endif
}

static void reset_default(Data *data)
{
    lock_update = TRUE;

    init_Instdata(&(data->args->instdata));

    gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->contactcontrols->nu), data->args->instdata.nu);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->contactcontrols->nui), data->args->instdata.nui);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->contactcontrols->Ei), data->args->instdata.Ei * 1e-9);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->contactcontrols->sigma_h), data->args->instdata.sigma_h);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->contactcontrols->sigma_F), data->args->instdata.sigma_F);

    op_recalculate(data);
    tangent_recalculate(data);
    hertz_recalculate(data);
#ifndef NOFORTRAN
    op_odr_fit_and_update(data);
    slopes_fit_and_update_load(data);
    slopes_fit_and_update_unload(data);
#endif

    lock_update = FALSE;
}

static void sigma_changed(Data *data)
{
    if (lock_update) {
        return;
    }

    data->args->instdata.sigma_h = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->contactcontrols->sigma_h));
    data->args->instdata.sigma_F = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->contactcontrols->sigma_F));
    data->args->instdata.delta = data->args->instdata.sigma_F * data->args->instdata.sigma_F / data->args->instdata.sigma_h / data->args->instdata.sigma_h;

    //label_set_gdouble_format(data->controls->contactcontrols->delta, data->args->instdata.delta, DELTA);

    gchar str[256];
    g_snprintf(str, sizeof(str), DELTA, data->args->instdata.delta);
    gtk_label_set_markup(GTK_LABEL(data->controls->contactcontrols->delta), str);

    /* recalculate fit and results */
#ifndef NOFORTRAN
    op_odr_fit_and_update(data);
    slopes_fit_and_update_load(data);
    slopes_fit_and_update_unload(data);
#endif
}

/*
static void
reset_delta(Data *data)
{
	data->args->instdata.sigma_h =  SIGMA_H_DEFAULT;
	data->args->instdata.sigma_F =  SIGMA_F_DEFAULT;
	data->args->instdata.delta = SIGMA_F_DEFAULT*SIGMA_F_DEFAULT/SIGMA_H_DEFAULT/SIGMA_H_DEFAULT;

	gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->contactcontrols->sigma_h), data->args->instdata.sigma_h);
	gtk_adjustment_set_value(GTK_ADJUSTMENT(data->controls->contactcontrols->sigma_F), data->args->instdata.sigma_F);

	label_set_gdouble_format(data->controls->contactcontrols->delta), data->args->instdata.delta, DELTA);

	// recalculate fit and results
	op_odr_fit_and_update(data);
	slopes_fit_and_update(data);

}
*/

void update_fd_plot(ContactControls *contactcontrols, FDdata *fddata)
{
    if (GWY_IS_DATA_LINE(fddata->hload) && GWY_IS_DATA_LINE(fddata->Fload))
        gwy_graph_curve_model_set_data(contactcontrols->cmodelload,
                                       fddata->hload->data, fddata->Fload->data,
                                       fddata->hload->res);

    if (GWY_IS_DATA_LINE(fddata->hhold) && GWY_IS_DATA_LINE(fddata->Fhold))
        gwy_graph_curve_model_set_data(contactcontrols->cmodelhold,
                                       fddata->hhold->data, fddata->Fhold->data,
                                       fddata->hhold->res);

    if (GWY_IS_DATA_LINE(fddata->hunload) && GWY_IS_DATA_LINE(fddata->Funload))
        gwy_graph_curve_model_set_data(contactcontrols->cmodelunload,
                                       fddata->hunload->data, fddata->Funload->data,
                                       fddata->hunload->res);
}

void update_fd_time_plot(ContactControls *contactcontrols, FDdata *fddata)
{
    gint timeres;
    gdouble xmin, xmax;
    gdouble Ftimexdata[2], Fcconst[2];

    timeres = fddata->ndata;

    if (GWY_IS_DATA_LINE(fddata->time)) {
        if (GWY_IS_DATA_LINE(fddata->htime))
            gwy_graph_curve_model_set_data(contactcontrols->cmodelhtime,
                                           fddata->time->data, fddata->htime->data,
                                           timeres);

        if (GWY_IS_DATA_LINE(fddata->Fload))
            gwy_graph_curve_model_set_data(contactcontrols->cmodelFtime,
                                           fddata->time->data, fddata->Ftime->data,
                                           timeres);
    }

    if (gwy_graph_model_get_x_range(contactcontrols->graph_model_fd_time, &(xmin), &(xmax))) {
        Ftimexdata[0] = (gint)(xmin);
        Ftimexdata[1] = (gint)(xmax);
        Fcconst[0] = fddata->Fc;
        Fcconst[1] = fddata->Fc;
        gwy_graph_curve_model_set_data(contactcontrols->cmodelFc, Ftimexdata, Fcconst, 2);
    }
}

static void plot_fd_curves(ContactControls *contactcontrols, FDdata *fddata)
{
    gwy_graph_model_add_curve(contactcontrols->graph_model_fd, contactcontrols->cmodelload);
    gwy_graph_model_add_curve(contactcontrols->graph_model_fd, contactcontrols->cmodelhold);
    gwy_graph_model_add_curve(contactcontrols->graph_model_fd, contactcontrols->cmodelunload);

    update_fd_plot(contactcontrols, fddata);
}

void update_max_labels(ContactControls *contactcontrols, const FDdata *fddata)
{
    label_set_gdouble_format(contactcontrols->hmax, fddata->hmax, DEPTH);
    label_set_gdouble_format(contactcontrols->Fmax, fddata->Fmax, FORCE);
}

void update_models_fd(ContactControls *contactcontrols)
{
    GwyGraphCurveModel *cmodel;

    contactcontrols->graph_model_fd = gwy_graph_get_model(GWY_GRAPH(contactcontrols->graph_fd));
    cmodel = gwy_graph_model_get_curve_by_description(contactcontrols->graph_model_fd, "Loading");

    if (cmodel) {
        contactcontrols->cmodelload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(contactcontrols->graph_model_fd, "Hold");

    if (cmodel) {
        contactcontrols->cmodelhold = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(contactcontrols->graph_model_fd, "Unloading");

    if (cmodel) {
        contactcontrols->cmodelunload = cmodel;
    }
}

void update_models_fd_time(ContactControls *contactcontrols)
{
    GwyGraphCurveModel *cmodel;

    contactcontrols->graph_model_fd_time = gwy_graph_get_model(GWY_GRAPH(contactcontrols->graph_fd_time));

    cmodel = gwy_graph_model_get_curve_by_description(contactcontrols->graph_model_fd_time, "h");

    if (cmodel) {
        contactcontrols->cmodelhtime = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(contactcontrols->graph_model_fd_time, "F");

    if (cmodel) {
        contactcontrols->cmodelFtime = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(contactcontrols->graph_model_fd_time, "Fc");

    if (cmodel) {
        contactcontrols->cmodelFc = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(contactcontrols->graph_model_fd_time, "Contact point");

    if (cmodel) {
        contactcontrols->cmodel_contact_load = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(contactcontrols->graph_model_fd_time, "Loading/hold");

    if (cmodel) {
        contactcontrols->cmodel_load_hold = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(contactcontrols->graph_model_fd_time, "Hold/unloading");

    if (cmodel) {
        contactcontrols->cmodel_hold_unload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(contactcontrols->graph_model_fd_time, "Unloading end");

    if (cmodel) {
        contactcontrols->cmodel_end_unload = cmodel;
    }
}

static void update_point_markline(GwyGraphModel *gmodel, GwyGraphCurveModel *cmodel,
                                  gdouble xmin, gdouble xmax, gdouble ymin, gdouble ymax, gint i, gchar *description)
{
    gdouble d;
    gdouble markline_x[2], markline_y[2];

    if ((xmin <= i) && (i <= xmax)) {
        d = (xmax - xmin) * 0.0005; /* little cheat to have the line displayed even when on the very edge of the graph */
        markline_x[0] = i - d;
        markline_x[1] = i + d;
        markline_y[0] = ymin;
        markline_y[1] = ymax;

        if (!gwy_graph_model_get_curve_by_description(gmodel, description)) {
            gwy_graph_model_add_curve(gmodel, cmodel);
        }

        gwy_graph_curve_model_set_data(cmodel, markline_x, markline_y, 2);
    }
    else {
        gwy_graph_model_remove_curve_by_description(gmodel, description);
    }
}

void update_point_marklines(ContactControls *contactcontrols, FDdata *fddata)
{
    gint i;
    gdouble Fi, xmin, xmax, ymin, ymax, yminh, yminF, ymaxh, ymaxF;
    gdouble Fcconst[2], Ftimexdata[2];

    if (gwy_graph_model_get_x_range(contactcontrols->graph_model_fd_time, &(xmin), &(xmax))) {

        if (gwy_graph_curve_model_get_ndata(contactcontrols->cmodelFtime) > 0) {
            Ftimexdata[0] = xmin;
            Ftimexdata[1] = xmax;

            if (((gint)xmin <= fddata->i_contact_load) && (fddata->i_contact_load <= (gint)xmax)) { // kontaktni bod je v grafu
                i = fddata->i_contact_load - (gint)xmin;
                Fi = gwy_graph_curve_model_get_ydata(contactcontrols->cmodelFtime)[i];
                Fcconst[0] = Fi;
                Fcconst[1] = Fi;
                gwy_graph_curve_model_set_data(contactcontrols->cmodelFc, Ftimexdata, Fcconst, 2);
                gwy_graph_model_remove_curve_by_description(contactcontrols->graph_model_fd_time, "Fc");
                gwy_graph_model_add_curve(contactcontrols->graph_model_fd_time, contactcontrols->cmodelFc);
            }
        }

        if (gwy_graph_curve_model_get_y_range(contactcontrols->cmodelhtime, &yminh, &ymaxh) &&
                gwy_graph_curve_model_get_y_range(contactcontrols->cmodelFtime, &yminF, &ymaxF)) {

            ymin = MIN(yminh, yminF);
            ymax = MAX(ymaxh, ymaxF);

            update_point_markline(contactcontrols->graph_model_fd_time, contactcontrols->cmodel_contact_load,
                                  xmin, xmax, ymin, ymax, fddata->i_contact_load, "Contact point");

            update_point_markline(contactcontrols->graph_model_fd_time, contactcontrols->cmodel_load_hold,
                                  xmin, xmax, ymin, ymax, fddata->i_load_hold, "Loading/hold");

            update_point_markline(contactcontrols->graph_model_fd_time, contactcontrols->cmodel_hold_unload,
                                  xmin, xmax, ymin, ymax, fddata->i_hold_unload, "Hold/unloading");

            update_point_markline(contactcontrols->graph_model_fd_time, contactcontrols->cmodel_end_unload,
                                  xmin, xmax, ymin, ymax, fddata->i_end_unload, "Unloading end");
        }
    }
}

static void contact_zoom_in_fd(Data *data)
{
    ContactControls *contactcontrols;

    contactcontrols = data->controls->contactcontrols;

    rs_zoom_in(GWY_GRAPH(contactcontrols->graph_fd), &(contactcontrols->zoomstack_fd), FALSE);   /* do not rescale y */
    update_models_fd(contactcontrols);
}

static void contact_zoom_in_fd_time(Data *data)
{
    ContactControls *contactcontrols;

    contactcontrols = data->controls->contactcontrols;

    rs_zoom_in(GWY_GRAPH(contactcontrols->graph_fd_time), &(contactcontrols->zoomstack_fd_time), TRUE);   /* do rescale y */
    update_models_fd_time(contactcontrols);
    update_point_marklines(contactcontrols, &(data->args->fddata));
}

static void contact_zoom_out_fd(Data *data)
{
    ContactControls *contactcontrols;

    contactcontrols = data->controls->contactcontrols;

    rs_zoom_out(GWY_GRAPH(contactcontrols->graph_fd), &(contactcontrols->zoomstack_fd));
    update_models_fd(contactcontrols);
}

static void contact_zoom_out_fd_time(Data *data)
{
    ContactControls *contactcontrols;

    contactcontrols = data->controls->contactcontrols;

    rs_zoom_out(GWY_GRAPH(contactcontrols->graph_fd_time), &(contactcontrols->zoomstack_fd_time));
    update_models_fd_time(contactcontrols);
    update_point_marklines(contactcontrols, &(data->args->fddata));
}

void contact_zoom_restore_fd(Data *data)
{
    ContactControls *contactcontrols;

    contactcontrols = data->controls->contactcontrols;

    rs_zoom_restore(GWY_GRAPH(contactcontrols->graph_fd), &(contactcontrols->zoomstack_fd));
    update_models_fd(contactcontrols);
}

void contact_zoom_restore_fd_time(Data *data)
{
    ContactControls *contactcontrols;

    contactcontrols = data->controls->contactcontrols;

    rs_zoom_restore(GWY_GRAPH(contactcontrols->graph_fd_time), &(contactcontrols->zoomstack_fd_time));
    update_models_fd_time(contactcontrols);
    update_point_marklines(contactcontrols, &(data->args->fddata));
}

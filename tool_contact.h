#ifndef TOOL_CONTACT_H
#define TOOL_CONTACT_H

#include "controls.h"
#include "datatypes.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

enum FDPlotMode {
    PLOT_FD = 0,
    PLOT_FD_TIME
};

struct ContactControls_ {
    GtkWidget *contact_gui;
    //GtkWidget *tool_label;

    /* Controls (left) */
    GtkWidget *vbox_ctrls;

    /* info */
    GtkWidget *frame_info;
    GtkWidget *table_info;
    GtkWidget *hmax, *Fmax;

    /* graph mode */
    GtkWidget *hbox_switch;
    GtkWidget *radio_fd, *radio_fd_time;
    enum FDPlotMode fd_plot_mode;

    /* NEW contact, load, hold, unload */
    GtkWidget *frame_contact_hold;
    GtkWidget *table_contact_hold;
    GtkWidget *status_contact_load, *status_load_hold, *status_hold_unload,
              *status_end_unload;
    GtkWidget *button_contact_point, *button_load_hold,	*button_hold_unload,
              *button_end_unload;
    GtkWidget *button_split_auto;

    /* save data */
    GtkWidget *button_save_data;

    /* area */
    GtkWidget *frame_area;
    GtkWidget *table_area;
    GtkWidget *area_function_button;
    GtkWidget *area_function_label;

    /* sample & indenter */
    GtkWidget *frame_sample_indenter;
    GtkWidget *table_sample_indenter;
    GtkObject *nu, *Ei, *nui;
    GtkWidget *button_reset;
    GtkObject *sigma_h, *sigma_F;
    GtkWidget  *delta, *reset_delta;

    /* Graph (right) */
    GtkWidget *hbox_graph;

    GtkWidget *graph_fd;
    GtkWidget *vbox_zoom_buttons_fd, *vbox_zoom_outer_fd;
    GtkWidget *button_zoom_in_fd, *button_zoom_out_fd, *button_zoom_restore_fd;

    GtkWidget *graph_fd_time;
    GtkWidget *vbox_zoom_buttons_fd_time, *vbox_zoom_outer_fd_time;
    GtkWidget *button_zoom_in_fd_time, *button_zoom_out_fd_time, *button_zoom_restore_fd_time;

    GwySelection *selection_fd, *selection_fd_time;

    /* stacks of graph models for zooming */
    GSList *zoomstack_fd, *zoomstack_fd_time;

    /* these should point to structures used by currently displayed graph */
    GwyGraphModel *graph_model_fd, *graph_model_fd_time;
    GwyGraphCurveModel *cmodelload, *cmodelhold, *cmodelunload;
    GwyGraphCurveModel *cmodelFtime, *cmodelhtime, *cmodelFc;
    GwyGraphCurveModel *cmodel_contact_load, *cmodel_load_hold, *cmodel_hold_unload, *cmodel_end_unload;
};


GtkWidget *tool_contact_create(Data *data);

void update_max_labels(ContactControls *contactcontrols, const FDdata *fddata);
void update_fd_plot(ContactControls *contactcontrols, FDdata *fddata);
void update_fd_time_plot(ContactControls *contactcontrols, FDdata *fddata);
void update_point_marklines(ContactControls *contactcontrols, FDdata *fddata);
void update_dummy_marklines(ContactControls *contactcontrols, FDdata *fddata);

void update_models_fd(ContactControls *contactcontrols);
void update_models_fd_time(ContactControls *contactcontrols);

void contact_switch_to_fd(Data *data);
void contact_switch_to_fd_time(Data *data);

void contact_zoom_restore_fd(Data *data);
void contact_zoom_restore_fd_time(Data *data);

#endif

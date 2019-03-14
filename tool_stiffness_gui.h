#ifndef NOFORTRAN

#ifndef TOOL_STIFFNESS_GUI_H
#define TOOL_STIFFNESS_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

struct StiffnessControls_ {
    GtkWidget *stiffness_gui;

    /* Ctrls (left) */
    GtkWidget *ctrls;

    /* info */
    GtkWidget *frame_info;
    GtkWidget *table_info;
    GtkWidget *hmax, *Fmax;

    /* parameters */
    GtkWidget *frame_parameters;
    GtkWidget *table_parameters;

    GtkWidget *from, *to;
    GtkWidget *beta, *reset_beta;

    /* fit buttons */
    GtkWidget *hbox_fitbuttons;
    GtkWidget *button_fit_load, *button_fit_unload;

    /* results */
    GtkWidget *frame_results;
    GtkWidget *scrolledwindow_results;
    GtkWidget *table_results;

    GtkWidget *fitinfoload, *fitinfounload;
    GtkWidget *kl, *ql, *ku, *qu;
    GtkWidget *loadrange;
    GtkWidget *unloadrange;
    GtkWidget *loadSSres, *loadR2, *loadR2adj, *loadchi2; /*quality of fit */
    GtkWidget *unloadSSres, *unloadR2, *unloadR2adj, *unloadchi2; /*quality of fit */

    /* unc, save */
    GtkWidget *button_save, *button_unc;

    /* Log */
    GtkWidget *hbox_logbuttons;
    GtkWidget *button_log_load, *button_log_unload;
    GtkWidget *load_log_error_window, *unload_log_error_window;
    GtkWidget *load_log_report_window, *unload_log_report_window;

    /* Graph (right) */
    GtkWidget *hbox_graph;

    GtkWidget *graph;
    GtkWidget *vbox_zoom_buttons, *vbox_zoom_outer;
    GtkWidget *button_zoom_in, *button_zoom_out, *button_zoom_restore;

    GwySelection *selection;

    GSList *zoomstack;

    /* these should point to structures used by currently displayed graph */
    GwyGraphModel *graph_model;
    GwyGraphCurveModel *cmodelload, *cmodelhold, *cmodelunload, *cmodelloadfit, *cmodelunloadfit;
};


GtkWidget *tool_stiffness_create(Data *data);
void stiffness_recalculate(Data *data);
void stiffness_redraw(Data *data);
void stiffness_remove_fit_results_labels(Data *data);
void stiffness_zoom_restore(Data *data);
void stiffness_fit_and_update_load(Data *data);
void stiffness_fit_and_update_unload(Data *data);

#endif

#endif

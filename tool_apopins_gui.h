#ifndef TOOL_APOPINS_GUI_H
#define TOOL_APOPINS_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

struct ApopinControls_ {
    GtkWidget *apop_gui;

    /* Controls (left) */
    GtkWidget *ctrls;

    /* parameters */
    GtkWidget *frame_parameters;
    GtkWidget *table_parameters;

    GtkObject *nmove, *thresh, *thresh2, *wpop, *hpop;
    GtkWidget *button_reset;

    /* run button */
    GtkWidget *button_run;

    /* results */
    GtkWidget *frame_results;
    GtkWidget *scrolledwindow_results;
    GtkWidget *table_results;

    GtkWidget *npopin;
    GtkWidget *res;

    /* save button */
    GtkWidget *button_save;

    /* Graphs (right) */
    GtkWidget *vbox_graphs;

    GtkWidget *hbox_graph_h;
    GtkWidget *graph_h;
    GtkWidget *vbox_zoom_buttons_h, *vbox_zoom_outer_h;
    GtkWidget *button_zoom_in_h, *button_zoom_out_h, *button_zoom_restore_h;

    GtkWidget *hbox_graph_der;
    GtkWidget *graph_der;
    GtkWidget *vbox_zoom_buttons_der, *vbox_zoom_outer_der;
    GtkWidget *button_zoom_in_der, *button_zoom_out_der, *button_zoom_restore_der;

    GSList *zoomstack_h, *zoomstack_der;

    GwyGraphModel *graph_model, *graph_model_der;
    GwyGraphCurveModel *cmodelload, *cmodelhavg, *cmodelder, *cmodelthresh, *cmodelthresh2;
};


GtkWidget *tool_apopins_create(Data *data);
void apopin_redraw(Data *data);
void apopin_remove_fit_results_labels(Data *data);
void apop_zoom_restore_linked(Data *data);

#endif

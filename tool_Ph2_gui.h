#ifndef TOOL_PH2_GUI_H
#define TOOL_PH2_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

struct Ph2Controls_ {
    GtkWidget *ph2_gui;

    /* Controls (left) */
    GtkWidget *ctrls;

    /* parameters */
    GtkWidget *frame_parameters;
    GtkWidget *table_parameters;

    GtkObject *nmove;

    /* run, save */
    GtkWidget *button_save, *button_run;

    /* Graphs (right) */
    GtkWidget *vbox_graphs;

    GtkWidget *hbox_graph_Fh;
    GtkWidget *graph_Fh;
    GtkWidget *vbox_zoom_buttons_Fh, *vbox_zoom_outer_Fh;
    GtkWidget *button_zoom_in_Fh, *button_zoom_out_Fh, *button_zoom_restore_Fh;

    GtkWidget *hbox_graph_H;
    GtkWidget *graph_H;
    GtkWidget *vbox_zoom_buttons_H, *vbox_zoom_outer_H;
    GtkWidget *button_zoom_in_H, *button_zoom_out_H, *button_zoom_restore_H;

    GSList *zoomstack_Fh, *zoomstack_H;

    GwyGraphModel *graph_modelFh, *graph_modelPh2;
    GwyGraphCurveModel *cmodelload, *cmodelavg, *cmodelPh2, *cmodeldPdh2;
};


GtkWidget *tool_Ph2_create(Data *data);
void Ph2_redraw(Data *data);
void Ph2_remove_fit_results_labels(Data *data);
void ph2_zoom_restore_linked(Data *data);

#endif

#ifndef TOOL_EPWORK_GUI_H
#define TOOL_EPWORK_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

struct WorkControls_ {
    GtkWidget *work_gui;

    /* Controls (left) */
    GtkWidget *ctrls;

    /* parameters */
    GtkWidget *frame_parameters;
    GtkWidget *table_parameters;

    GtkObject *nmove;

    /* run button */
    GtkWidget *button_run;

    /* results */
    GtkWidget *frame_results;
    GtkWidget *table_results;

    GtkWidget *workplastic, *workelastic;
    GtkWidget *eta;

    /* unc, save */
    GtkWidget *button_save,  *button_unc;

    /* Graph (right) */
    GtkWidget *hbox_graph;

    GtkWidget *graph;
    GtkWidget *vbox_zoom_buttons, *vbox_zoom_outer;
    GtkWidget *button_zoom_in, *button_zoom_out, *button_zoom_restore;

    GSList *zoomstack;

    /* these should point to structures used by currently displayed graph */
    GwyGraphModel *graph_model;
    GwyGraphCurveModel *cmodelload, *cmodelloadavg;
    GwyGraphCurveModel *cmodelhold, *cmodelholdavg;
    GwyGraphCurveModel *cmodelunload, *cmodelunloadavg;
};


GtkWidget *tool_epwork_create(Data *data);
void work_redraw(Data *data);
void work_remove_fit_results_labels(Data *data);
void work_zoom_restore(Data *data);

#endif

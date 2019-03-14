#ifndef TOOL_HERTZ_GUI_H
#define TOOL_HERTZ_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

struct HertzControls_ {
    GtkWidget *hertz_gui;

    /* Controls (left) */
    GtkWidget *ctrls;

    /* info */
    GtkWidget *frame_info;
    GtkWidget *table_info;

    GtkWidget *hmax, *Fmax;

    /* parameters */
    GtkWidget *frame_parameters;
    GtkWidget *table_parameters;

    GtkWidget *from, *to;

    /*radio buttons to switch between and R and E mode */
    GtkWidget *radio_R, *radio_Er, *radio_Eit;

    GtkWidget *radiusinp, *Eitinp, *Erinp;	/*input according to radio buttons */

    /* fit button */
    GtkWidget *button_fit;

    /* results */
    GtkWidget *frame_results;
    GtkWidget *scrolledwindow_results;
    GtkWidget *table_results;

    GtkWidget *gamma;
    GtkWidget *SSres, *R2, *R2adj, *chi2; /*quality of fit */
    GtkWidget *Erres, *Eitres, *radiusres; /*display results */
    GtkWidget *range;

    /* unc, save */
    GtkWidget *button_unc, *button_save;

    /* Graph (right) */
    GtkWidget *hbox_graph;

    GtkWidget *graph;
    GtkWidget *vbox_zoom_buttons, *vbox_zoom_outer;
    GtkWidget *button_zoom_in, *button_zoom_out, *button_zoom_restore;

    GwySelection *selection;

    GSList *zoomstack;

    GwyGraphModel *graph_model;
    GwyGraphCurveModel *cmodelload, *cmodelfit;
};


GtkWidget *tool_hertz_create(Data *data);
void hertz_recalculate(Data *data);
void hertz_redraw(Data *data);
void hertz_remove_fit_results_labels(Data *data);
void hz_zoom_restore(Data *data);

#endif

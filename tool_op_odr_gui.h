#ifndef NOFORTRAN

#ifndef TOOL_OP_ODR_GUI_H
#define TOOL_OP_ODR_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

struct OPODRControls_ {
    GtkWidget *opodr_gui;

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
    GtkWidget *from_pct_Fmax, *to_pct_Fmax;
    GtkWidget *beta, *reset_beta;
    GtkWidget *hbox_radial, *toggle_radial, *radial_angle;

    /* fit button */
    GtkWidget *button_fit;

    /* results */
    GtkWidget *frame_results;
    GtkWidget *scrolledwindow_results;
    GtkWidget *table_results;

    GtkWidget *fitinfo;
    GtkWidget *hp, *m, *alpha, *eps, *hc, *S, *Aphc;
    GtkWidget *SSres, *R2, *R2adj, *chi2; /*quality of fit */
    GtkWidget *extrapol;
    GtkWidget *Hit, *Er, *Eit;
    GtkWidget *range, *range_pct_Fmax;
    GtkWidget *radial_info;

    /* log, unc, save */
    GtkWidget *button_log, *button_unc, *button_save;

    /* log windows */
    GtkWidget *log_error_window, *log_report_window;

    /* Graph (right) */
    GtkWidget *hbox_graph;

    GtkWidget *graph;
    GtkWidget *vbox_zoom_buttons, *vbox_zoom_outer;
    GtkWidget *button_zoom_in, *button_zoom_out, *button_zoom_restore;

    GwySelection *selection;

    GSList *zoomstack;

    /* these should point to structures used by currently displayed graph */
    GwyGraphModel *graph_model;
    GwyGraphCurveModel *cmodelunload, *cmodelSfit;
};

GtkWidget *tool_op_odr_create(Data *data);
void op_odr_recalculate(Data *data);
void op_odr_redraw(Data *data);
void op_odr_remove_fit_results_labels(Data *data);
void op_odr_zoom_restore(Data *data);
void op_odr_fit_and_update(Data *data);

#endif

#endif

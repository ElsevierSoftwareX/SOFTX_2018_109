#ifndef NOFORTRAN

#ifndef TOOL_SLOPES_GUI_H
#define TOOL_SLOPES_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

struct SlopesControls_ {
    GtkWidget *slopes_gui;

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
    GtkWidget *from_pct_Fmax, *to_pct_Fmax;
    GtkWidget *beta, *reset_beta;
    GtkWidget *fixgamma; /*whether or not gamma should be fixed */

    /* fit buttons */
    GtkWidget *hbox_fitbuttons;
    GtkWidget *button_fit_load, *button_fit_unload;

    /* results */
    GtkWidget *frame_results;
    GtkWidget *scrolledwindow_results;
    GtkWidget *table_results;

    GtkWidget *fitinfoload, *fitinfounload;
    GtkWidget *hp, *m, *eps, *h0, *n, *Sload, *Sunload;
    GtkWidget *Hit, *Er, *Eit, *Aphc;
    GtkWidget *loadrange, *loadrange_pct_Fmax;
    GtkWidget *unloadrange, *unloadrange_pct_Fmax;
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


GtkWidget *tool_slopes_create(Data *data);
void slopes_recalculate(Data *data);
void slopes_redraw(Data *data);
void slopes_remove_fit_results_labels(Data *data);
void slopes_zoom_restore(Data *data);
void slopes_fit_and_update_load(Data *data);
void slopes_fit_and_update_unload(Data *data);

#endif

#endif

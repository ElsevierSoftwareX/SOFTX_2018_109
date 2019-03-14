#ifndef TOOL_TANGENT_GUI_H
#define TOOL_TANGENT_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

struct TangentControls_ {
    GtkWidget *tangent_gui;

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

    /* fit button */
    GtkWidget *button_fit;

    /* results */
    GtkWidget *frame_results;
    GtkWidget *scrolledwindow_results;
    GtkWidget *table_results;

    GtkWidget *hc, *S, *Aphc;
    GtkWidget *SSres, *R2, *R2adj, *chi2; /*quality of fit */
    GtkWidget *extrapol;
    GtkWidget  *R; /* created but not attached */
    GtkWidget *Hit, *Er, *Eit;
    GtkWidget *range, *range_pct_Fmax;

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
    GwyGraphCurveModel *cmodelunload, *cmodelSfit;
};

GtkWidget *tool_tangent_create(Data *data);
void tangent_recalculate(Data *data);
void tangent_redraw(Data *data);
void tangent_remove_fit_results_labels(Data *data);
void tg_zoom_restore(Data *data);

#define FIT_TLS TRUE
#define CORRECT FALSE
#define RESCALE FALSE
#define DELT  TRUE

#endif

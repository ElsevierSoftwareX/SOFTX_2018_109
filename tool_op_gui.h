#ifndef TOOL_OP_GUI_H
#define TOOL_OP_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

struct OPControls_ {
    GtkWidget *op_gui;

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

    /* fit buttons */
    GtkWidget *hbox_fitbuttons;
    GtkWidget *button_fit_hp, *button_fit_S;

    /* results */
    GtkWidget *frame_results;
    GtkWidget *scrolledwindow_results;
    GtkWidget *table_results;

    GtkWidget *hp, *m, *eps, *hc, *S, *Aphc;
    GtkWidget *hpSSres, *hpR2, *hpR2adj, *hpchi2; /*quality of fit */
    GtkWidget *SSSres, *SR2, *SR2adj, *Schi2; /*quality of fit */
    GtkWidget *hpR; /* created but not attached */
    GtkWidget *SR; /* created but not attached */
    GtkWidget *extrapol;
    GtkWidget *Hit, *Er, *Eit;
    GtkWidget *hprange, *hprange_pct_Fmax;
    GtkWidget *Srange, *Srange_pct_Fmax;

    /* unc, save */
    GtkWidget *button_save, *button_unc;

    /* Graph (right) */
    GtkWidget *hbox_graph;

    GtkWidget *graph;
    GtkWidget *vbox_zoom_buttons, *vbox_zoom_outer;
    GtkWidget *button_zoom_in, *button_zoom_out, *button_zoom_restore;

    GwySelection *selection;

    GSList *zoomstack;

    /* these should point to structures used by currently displayed graph */
    GwyGraphModel *graph_model;
    GwyGraphCurveModel *cmodelunload, *cmodelhpfit, *cmodelSfit;
};

#define FIT_S_LR  TRUE // TRUE - Linear REgreeession, FALSE - Nonlinear Regression

GtkWidget* tool_op_create(Data *data);
void op_recalculate(Data *data);
void op_redraw(Data *data);
void op_remove_fit_results_labels(Data *data);
void op_zoom_restore(Data *data);

#endif

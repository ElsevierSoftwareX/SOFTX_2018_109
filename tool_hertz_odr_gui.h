#ifndef NOFORTRAN

#ifndef TOOL_HERTZ_ODR_GUI_H
#define TOOL_HERTZ_ODR_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

struct HertzODRControls_ {
    GtkWidget *hertz_odr_gui;

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

    GtkWidget *fixgamma; /*whether or not gamma should be fixed */
    GtkWidget *toggle_radial;

    /* fit button */
    GtkWidget *button_fit;

    /* results */
    GtkWidget *frame_results;
    GtkWidget *scrolledwindow_results;
    GtkWidget *table_results;

    GtkWidget *fitinfo;
    GtkWidget *gamma, *n, *h0;
    GtkWidget *SSres, *R2, *R2adj, *chi2; /*quality of fit */
    GtkWidget *Erres, *Eitres, *radiusres; /*display results */
    GtkWidget *range;

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

    GwyGraphModel *graph_model;
    GwyGraphCurveModel *cmodelload, *cmodelfit;
};


GtkWidget *tool_hertz_odr_create(Data *data);
void hertz_odr_recalculate(Data *data);
void hertz_odr_redraw(Data *data);
void hertz_odr_remove_fit_results_labels(Data *data);
void hz_odr_zoom_restore(Data *data);

#endif

#endif

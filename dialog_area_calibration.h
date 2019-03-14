#ifndef DIALOG_AREA_CALIBRATION_H
#define DIALOG_AREA_CALIBRATION_H

#include "controls.h"

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

struct AreaControls_ {
    GtkWidget *dialog;
    GtkWidget *hbox;

    GtkWidget *area_table;

    /* OP */
    GtkWidget *opp_button;
    GtkWidget **area_polynom;

    /* Raw data */
    GtkWidget *rd_button;
    GtkWidget *area_entry;
    GtkWidget *area_file_button;

    /* Polynom data */
    GtkWidget *coeff_button;
    GtkWidget *coeff_entry;
    GtkWidget *coeff_file_button;

    /* Graph */
    GtkWidget *area_graph;
    GwyGraphModel *gmodel;
    GwyGraphCurveModel *cmodel;
};


void area_calibration(GtkWidget *widget, Data *data);

#endif

#ifndef NOFORTRAN

#ifndef TOOL_STIFFNESS_UNC_GUI_H
#define TOOL_STIFFNESS_UNC_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

struct StiffnessUncControls_ {
    GtkObject *uh, *uF;

    GtkWidget *ukl, *uql, *uku, *uqu;

    GtkWidget *window;

    GtkObject *Nmc;
    GtkWidget *spin_Nmc;
    GtkWidget *mc_warn;

    GtkWidget *vbox;
    GtkWidget *frame_unc_input, *frame_unc_contact, *frame_mc;
    GtkWidget *table_unc_input, *table_unc_contact, *table_mc;
    GtkWidget *button_save, *button_mc;

};

struct StiffnessMCControls_ {
    GtkWidget *window;
    gboolean has_window;
    GtkWidget *vbox;
    GtkWidget *frame_input, *frame_results;
    GtkWidget *table_input, *table_results;
    GtkWidget *button_save;
    GtkWidget **button_hist;
    GtkWidget **hist_window;
    //	gboolean *has_hist_window;
};


void stiffness_uncertainty(Data *data);
void stiffness_propagate_uncertainties(StiffnessUncControls *ctrls, const Stiffnessdata *stiffnessdata, StiffnessUncdata *unc, const FDdata *fddata);
void stiffness_unc_close_window(Data *data);

void stiffness_mc_close_window(Data *data);
void stiffness_uncertainty_montecarlo(Data *data);

#endif

#endif

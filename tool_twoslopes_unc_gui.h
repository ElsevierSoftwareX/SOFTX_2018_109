#ifndef NOFORTRAN

#ifndef TOOL_SLOPES_UNC_GUI_H
#define TOOL_SLOPES_UNC_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

struct SlopesUncControls_ {
    GtkObject *uh, *uF;
    GtkObject *unu, *unui, *uEi;

    GtkWidget *uEitunu, *uEitunui, *uEituEi;

    GtkWidget *uSu, *uSl, *un, *um, *uA, *uHit, *uEr, *uEit;

    GtkWidget *window;

    GtkObject *Nmc;
    GtkWidget *spin_Nmc;
    GtkWidget *mc_warn;

    GtkWidget *vbox;
    GtkWidget *frame_unc_input, *frame_unc_contact, *frame_mc;
    GtkWidget *table_unc_input, *table_unc_contact, *table_mc;
    GtkWidget *button_save, *button_mc;
};

struct SlopesMCControls_ {
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


void slopes_uncertainty(Data *data);
void slopes_propagate_uncertainties(SlopesUncControls *ctrls, const Slopesdata *slopesdata, SlopesUncdata *unc, const FDdata *fddata);
void slopes_unc_close_window(Data *data);

void slopes_mc_close_window(Data *data);
void slopes_uncertainty_montecarlo(Data *data);

#endif

#endif

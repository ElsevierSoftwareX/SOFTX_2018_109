#ifndef TOOL_TANGENT_UNC_GUI_H
#define TOOL_TANGENT_UNC_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

struct TangentUncControls_ {
    GtkObject *uh, *uF;
    GtkObject *unu, *unui, *uEi;

    GtkWidget *uEitunu, *uEitunui, *uEituEi;

    GtkWidget *uS, *uhc, *uA, *uHit, *uEr, *uEit;

    GtkObject *Nmc;
    GtkWidget *mc_warn;

    GtkWidget *window;
    GtkWidget *vbox;
    GtkWidget *frame_unc_input, *frame_unc_contact, *frame_mc;
    GtkWidget *table_unc_input, *table_unc_contact, *table_mc;
    GtkWidget *button_save, *button_mc;
};

struct TangentMCControls_ {
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


void tangent_uncertainty(Data *data);
void tangent_propagate_uncertainties(TangentUncControls *ctrls, Tangentdata *tgdata, TangentUncdata *unc, FDdata *fddata, Instdata *instdata, Area *area);
void tangent_unc_close_window(Data *data);

void tangent_mc_close_window(Data *data);
void tangent_uncertainty_montecarlo(Data *data);

#endif

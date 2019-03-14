#ifndef TOOL_OP_UNC_GUI_H
#define TOOL_OP_UNC_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

struct OPUncControls_ {
    GtkObject *uh, *uF;
    GtkObject *unu, *unui, *uEi;

    GtkWidget *uEitunu, *uEitunui, *uEituEi;

    GtkWidget *uS, *uhc, *uA, *uHit, *uEr, *uEit;

    GtkWidget *window;

    GtkObject *Nmc;
    GtkWidget *mc_warn;

    GtkWidget *vbox;
    GtkWidget *frame_unc_input, *frame_unc_contact, *frame_mc;
    GtkWidget *table_unc_input, *table_unc_contact, *table_mc;
    GtkWidget *button_save, *button_mc;
};

struct OPMCControls_ {
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

void op_uncertainty(Data *data);
void op_propagate_uncertainties(OPUncControls *ctrls, const OPdata *opdata, OPUncdata *opunc,
				const FDdata *fddata, const Instdata *instdata, const Area *area);
void op_unc_close_window(Data *data);

void op_mc_close_window(Data *data);
void op_uncertainty_montecarlo(Data *data);

#endif

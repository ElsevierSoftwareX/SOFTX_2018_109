#ifndef NOFORTRAN

#ifndef TOOL_OPODR_UNC_GUI_H
#define TOOL_OPODR_UNC_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

struct OPODRUncControls_ {
    GtkObject *uh, *uF;
    GtkObject *unu, *unui, *uEi;

    GtkWidget *uEitunu, *uEitunui, *uEituEi;

    GtkWidget *uS, *um, *uhc, *uA, *uHit, *uEr, *uEit;

    GtkWidget *window;

    GtkObject *Nmc;
    GtkWidget *spin_Nmc;
    GtkWidget *mc_warn;

    GtkWidget *vbox;
    GtkWidget *frame_unc_input, *frame_unc_contact, *frame_mc;
    GtkWidget *table_unc_input, *table_unc_contact, *table_mc;
    GtkWidget *button_save, *button_mc;
};

struct OPODRMCControls_ {
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

struct OPODRUnc2Controls_ {
    GtkWidget *window;
    GtkWidget *notebook;

    /* Instrument (unc. prop.) tab */
    GtkWidget *vbox_prop;
    GtkWidget *table_unc_input;

    GtkObject *uh, *uF, *unu, *unui, *uEi;

    GtkWidget *uEitunu, *uEitunui, *uEituEi;
    GtkWidget *um, *uS, *uhc, *uA, *uHit, *uEr, *uEit;

    GtkWidget *button_save;

    /* MC tab */
    GtkWidget *vbox_mc;

    GtkWidget *frame_input;
    GtkWidget *table_input;

    GtkObject *uh_mc, *uF_mc, *N_mc;
    GtkWidget *spin_N_mc;
    GtkWidget *mc_warn;

    GtkWidget *button_run;

    GtkWidget *frame_results;
    GtkWidget *table_results;

    GtkWidget *button_save_mc;
    GtkWidget **button_hist;
    GtkWidget **hist_window;
    //gboolean *has_hist_window;

    /* Contact Point tab */
    GtkWidget *vbox_cp;
    GtkObject *shift_range;
    GtkWidget *table_unc_contact;
};


void op_odr_uncertainty(Data *data);
void op_odr_uncertainty2(Data *data);
void op_odr_propagate_uncertainties(OPODRUncControls *ctrls, const OPODRdata *opodrdata, OPODRUncdata *unc,
                                    const FDdata *fddata, const Area *area);
void op_odr_unc_close_window(Data *data);

void op_odr_mc_close_window(Data *data);

#endif

#endif

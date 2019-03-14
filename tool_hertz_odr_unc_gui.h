#ifndef TOOL_HERTZ_ODR_UNC_GUI_H
#define TOOL_HERTZ_ODR_UNC_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

struct HertzODRUncControls_ {
    GtkObject *uh, *uF;
    GtkObject *unu, *unui, *uEi;
    GtkObject *uradius;
    GtkObject *uEr, *uEit;

    GtkWidget *uEitunu, *uEitunui, *uEituEi;
    GtkWidget *uEruradius, *uEituradius;

    GtkWidget *uErtotal, *uEittotal, *uradiustotal;
    GtkWidget *uradiusuh, *uradiusuF;
    GtkWidget *uradiusunu, *uradiusunui, *uradiusuEi;
    GtkWidget *uradiusuEr, *uradiusuEit;

    GtkObject *Nmc;
    GtkWidget *mc_warn;

    GtkWidget *window;
    GtkWidget *vbox;
    GtkWidget *frame_unc_input, *frame_unc_contact, *frame_mc;
    GtkWidget *table_unc_input, *table_unc_contact, *table_mc;
    GtkWidget *button_save, *button_mc;
};

struct HertzODRMCControls_ {
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


void hertz_odr_uncertainty(Data *data);
void hertz_odr_propagate_uncertainties(HertzODRUncControls *ctrls, HertzODRdata *hzodr, HertzODRUncdata *unc, const FDdata *fddata);
void hertz_odr_unc_close_window(Data *data);

#endif

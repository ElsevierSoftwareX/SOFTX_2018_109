#ifndef TOOL_HERTZ_UNC_GUI_H
#define TOOL_HERTZ_UNC_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

struct HertzUncControls_ {
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

struct HertzMCControls_ {
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


void hertz_uncertainty(Data *data);
void hertz_propagate_uncertainties(HertzUncControls *ctrls, const Hertzdata *hz, HertzUncdata *unc, const FDdata *fddata, const Instdata *instdata);
void hertz_unc_close_window(Data *data);

void hertz_mc_close_window(Data *data);
void hertz_uncertainty_montecarlo(Data *data);

#endif

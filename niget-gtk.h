#ifndef NIGET_H
#define NIGET_H

#include <gtk/gtk.h>

struct MainControls_ {
    GtkWidget *window_main;
    GtkWidget *vbox_main;

    GtkWidget *notebook_tools;
    GtkWidget *hbox_buttons;

    GtkWidget *tab_contact, *tab_op, *tab_op_odr, *tab_tangent,
              *tab_hertz, *tab_hertz_odr, *tab_stiffness,
              *tab_slopes, *tab_ph2, *tab_apopins, *tab_epwork ;

    GtkWidget *tab_contact_label, *tab_op_label, *tab_op_odr_label, *tab_tangent_label,
              *tab_hertz_label, *tab_hertz_odr_label, *tab_stiffness_label,
              *tab_slopes_label, *tab_ph2_label, *tab_apopins_label, *tab_epwork_label;

    GtkWidget *button_open, *button_save, *button_help, *button_about, *button_quit;
    GtkWidget *dialog_about;
};

#endif

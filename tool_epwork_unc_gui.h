#ifndef TOOL_EPWORK_UNC_GUI_H
#define TOOL_EPWORK_UNC_GUI_H

#include "controls.h"

#include <gtk/gtk.h>

struct WorkUncControls_ {
    GtkWidget *window;
    GtkWidget *table_contact;
};

void work_uncertainty(Data *data);
void work_propagate_uncertainties(Data *data);
void work_unc_close_window(Data *data);

#endif

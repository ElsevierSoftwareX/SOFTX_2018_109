#include "tool_epwork_unc_gui.h"
#include "tool_epwork_unc.h"

#include "controls.h"
#include "datatypes.h"
#include "file-utils.h"
#include "gui-utils.h"
#include "niget-common.h"

#include <gtk/gtk.h>

#include <app/file.h>

static void work_uncertainty_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (work_unc_has_results(&(args->workunc))) {
	buffer = work_uncertainty_export_data(&(args->workunc), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

void work_unc_close_window(Data *data)
{
    g_free(data->args->workunc.We);
    g_free(data->args->workunc.Wp);
    data->args->workunc.We = NULL;
    data->args->workunc.Wp = NULL;
    gtk_widget_destroy(data->controls->workunccontrols->window);
}

void work_uncertainty(Data *data)
{
    GtkWidget  *vbox, *label;
    GtkWidget *button_save;
    gint i, j;
    WorkUncControls *ctrls;
    WorkUncdata *unc;
    gchar str[50];

    unc = &(data->args->workunc);
    ctrls =	data->controls->workunccontrols;

    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Uncertainty in elastic/plastic work");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(work_unc_close_window), data);

    ctrls->table_contact = gtk_table_new(20, 5, FALSE);
    i = 0;
    gtk_table_attach(GTK_TABLE(ctrls->table_contact), gwy_label_new_header("Uncertainties due to choice of contact point"), 0, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;
    gtk_table_attach(GTK_TABLE(ctrls->table_contact), gwy_label_new_header("shift contact point"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_contact), gwy_label_new_header("W<sub>el</sub>/GPa"), 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_contact), gwy_label_new_header("W<sub>pl</sub>/GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    //	gtk_table_attach(GTK_TABLE(ctrls->table_contact), gwy_label_new_header("weight"), 2, 3, i,i+1, GTK_FILL, 0, 0, 0);
    i++;

    unc->We = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
    unc->Wp = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));

    for (j = -CONTACTMAX; j <= CONTACTMAX ; j++) {
        label = gtk_label_new(NULL);
        g_snprintf(str, sizeof(str), "%d\n", j);
        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_table_attach(GTK_TABLE(ctrls->table_contact), label, 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        work_run_shift_contact(&(data->args->workdata), &(data->args->fddata), unc->We, unc->Wp, j);

        label = gtk_label_new(NULL);

        if (unc->We[j + CONTACTMAX] > 0) {
            g_snprintf(str, sizeof(str), EPWORK, unc->We[j + CONTACTMAX]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_contact), label, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        label = gtk_label_new(NULL);

        if (unc->Wp[j + CONTACTMAX] > 0) {
            g_snprintf(str, sizeof(str), EPWORK, unc->Wp[j + CONTACTMAX]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_contact), label, 2, 3, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        i++;
    }

    button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(button_save, "clicked", G_CALLBACK(work_uncertainty_save_data), data->args);

    vbox = gtk_vbox_new(TRUE, 4);
    gtk_box_pack_start(GTK_BOX(vbox), ctrls->table_contact, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), button_save,  TRUE, FALSE, 0);
    gtk_container_add(GTK_CONTAINER(ctrls->window), vbox);
    gtk_widget_show_all(ctrls->window);
}

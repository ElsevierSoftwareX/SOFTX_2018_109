#include "tool_hertz_unc_gui.h"
#include "tool_hertz_unc.h"
#include "tool_hertz_gui.h" /* only because of touching hertzcontrols */

#include "controls.h"
#include "datatypes.h"
#include "gui-utils.h"
#include "niget-common.h"

#include "niget-gtk.h" /* controls definitions */

#include <math.h>

#include <gtk/gtk.h>

#define MC_WHEN_WARN 1e6

void hertz_propagate_uncertainties(HertzUncControls *ctrls, const Hertzdata *hz, HertzUncdata *unc,
				   const FDdata *fddata, const Instdata *instdata)
{
    gchar str[30];
    gint i, j;
    GtkWidget *label;

    if (!hz->has_fit) {
        return;
    }

    if (unc->ERc == NULL) {
        return;
    }

    // calculate the uncertainties
    hertz_propagate_uncertainties_calc(hz, unc, instdata);

    /* show */
    switch (hz->mode) {
    case R_MODE:
        // radius
        label_set_gdouble_format(ctrls->uEruradius, unc->uEruradius, MODUL);
        label_set_gdouble_format(ctrls->uEituradius, unc->uEituradius, MODUL);

        //other contributions
        label_set_gdouble_format(ctrls->uEitunu, unc->uEitunu, MODUL);
        label_set_gdouble_format(ctrls->uEitunui, unc->uEitunui, MODUL);
        label_set_gdouble_format(ctrls->uEituEi, unc->uEituEi, MODUL);

        // total uncertainties
        label_set_gdouble_format(ctrls->uErtotal, unc->uErtotal, MODUL);
        label_set_gdouble_format(ctrls->uEittotal, unc->uEittotal, MODUL);

        break;

    case ER_MODE:
        /* show */
        label_set_gdouble_format(ctrls->uradiusuh, unc->uradiusuh * 1e9, DEPTH);
        label_set_gdouble_format(ctrls->uradiusuF, unc->uradiusuF * 1e9, DEPTH);
        label_set_gdouble_format(ctrls->uradiusuEr, unc->uradiusuEr * 1e9, DEPTH);

        // total uncertainties
        label_set_gdouble_format(ctrls->uradiustotal, unc->uradiustotal * 1e9, DEPTH);

        break;

    case EIT_MODE :
        /* show */
        label_set_gdouble_format(ctrls->uradiusuh, unc->uradiusuh * 1e9, DEPTH);
        label_set_gdouble_format(ctrls->uradiusuF, unc->uradiusuF * 1e9, DEPTH);
        label_set_gdouble_format(ctrls->uradiusuEi, unc->uradiusuEi * 1e9, DEPTH);
        label_set_gdouble_format(ctrls->uradiusuEit, unc->uradiusuEit * 1e9, DEPTH);
        label_set_gdouble_format(ctrls->uradiusunui, unc->uradiusunui * 1e9, DEPTH);
        label_set_gdouble_format(ctrls->uradiusunu, unc->uradiusunu * 1e9, DEPTH);

        // total uncertainties
        label_set_gdouble_format(ctrls->uradiustotal, unc->uradiustotal * 1e9, DEPTH);

        break;
    }

    /* contact point uncertainties */

    //destroy the labels of the old table_unc_contact
    remove_table_content(ctrls->table_unc_contact);

    //create new table
    i = 0;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), gwy_label_new_header("Uncertainties due to choice of contact point"), 0, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), gwy_label_new_header("shift contact point"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    if (hz->mode == R_MODE)
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                         gwy_label_new_header("E<sub>r</sub>/GPa"), 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    else
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                         gwy_label_new_header("radius/nm"), 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    i++;

    for (j = -CONTACTMAX; j <= CONTACTMAX ; j++) {
        label = gtk_label_new(NULL);
        label_set_gint_format(label, j, "%d\n");
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        hertz_fit_shift_contact(hz, fddata, instdata, unc->ERc, j);
        label = gtk_label_new(NULL);

        if (unc->ERc[j + CONTACTMAX] > 0)
            if (hz->mode == R_MODE) {
                g_snprintf(str, sizeof(str), MODUL, unc->ERc[j + CONTACTMAX]);
            }
            else {
                g_snprintf(str, sizeof(str), DEPTH, unc->ERc[j + CONTACTMAX] * 1e9);
            }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        label = gtk_label_new(NULL);
        g_snprintf(str, sizeof(str), "1/%d", 2 * CONTACTMAX + 1); /* co s tim pak? */
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 2, 3, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        i++;

    }

    gtk_widget_show_all(ctrls->table_unc_contact);
}

static void mc_iterations_changed(Data *data)
{
    data->args->hertzunc.Nmc = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->hertzunccontrols->Nmc));

    if (data->args->hertzunc.Nmc > MC_WHEN_WARN) {
        label_set_gint_format(data->controls->hertzunccontrols->mc_warn, data->args->hertzunc.Nmc,
                              "Warning: Monte Carlo calculation for %d iterations may take very long!");
    }
    else {
        label_clear(data->controls->hertzunccontrols->mc_warn);
    }
}

static void unc_depth_changed(Data *data)
{
    data->args->hertzunc.uh = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->hertzunccontrols->uh));
    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}

static void unc_load_changed(Data *data)
{
    data->args->hertzunc.uF = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->hertzunccontrols->uF));
    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}

static void unc_radius_changed(Data *data)
{
    data->args->hertzunc.uradius = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->hertzunccontrols->uradius)) * 1e-9;
    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}

static void unc_nu_changed(Data *data)
{
    data->args->instdata.unu = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->hertzunccontrols->unu));
    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}

static void unc_nui_changed(Data *data)
{
    data->args->instdata.unui = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->hertzunccontrols->unui));
    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}

static void unc_Ei_changed(Data *data)
{
    data->args->instdata.uEi = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->hertzunccontrols->uEi)) * 1e9;
    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}

static void unc_Er_changed(Data *data)
{
    data->args->hertzunc.uEr = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->hertzunccontrols->uEr));
    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}


static void unc_Eit_changed(Data *data)
{
    data->args->hertzunc.uEit = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->hertzunccontrols->uEit));
    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));
}

static void hertz_uncertainty_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (hertz_unc_has_results(&(args->hertzunc))) {
	buffer = hertz_uncertainty_export_data(&(args->hertzdata), &(args->hertzunc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

void hertz_unc_close_window(Data *data)
{
    g_free(data->args->hertzunc.ERc);
    data->args->hertzunc.ERc = NULL;
    hertz_mc_close_window(data);

    gtk_widget_set_sensitive(GTK_WIDGET(data->controls->hertzcontrols->radio_Er), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(data->controls->hertzcontrols->radio_Eit), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(data->controls->hertzcontrols->radio_R), TRUE);
    gtk_widget_destroy(data->controls->hertzunccontrols->window);
}

void hertz_uncertainty(Data *data)
{
    GtkWidget  *spin;
    gint i;
    HertzUncControls *ctrls;
    HertzUncdata *unc;
    Instdata *inst;

    if (data->args->hertzunc.ERc != NULL) {
        gtk_window_present(GTK_WINDOW(data->controls->hertzunccontrols->window));
        return;
    }

    /* deactivate radiobutton for mode, so that it cannot be changed while the uncertainty dialog is open
     * adjustment of value remains active */

    switch (data->args->hertzdata.mode) {
    case R_MODE:
        gtk_widget_set_sensitive(GTK_WIDGET(data->controls->hertzcontrols->radio_Er), FALSE);
        gtk_widget_set_sensitive(GTK_WIDGET(data->controls->hertzcontrols->radio_Eit), FALSE);
        break;

    case ER_MODE:
        gtk_widget_set_sensitive(GTK_WIDGET(data->controls->hertzcontrols->radio_R), FALSE);
        gtk_widget_set_sensitive(GTK_WIDGET(data->controls->hertzcontrols->radio_Eit), FALSE);
        break;

    case EIT_MODE:
        gtk_widget_set_sensitive(GTK_WIDGET(data->controls->hertzcontrols->radio_R), FALSE);
        gtk_widget_set_sensitive(GTK_WIDGET(data->controls->hertzcontrols->radio_Er), FALSE);
        break;
    }

    unc = &(data->args->hertzunc);
    ctrls =	data->controls->hertzunccontrols;
    inst = &(data->args->instdata);

    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Uncertainty in Hertz analysis");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(hertz_unc_close_window), data);
    gtk_window_set_transient_for(GTK_WINDOW(ctrls->window), GTK_WINDOW(data->controls->maincontrols->window_main));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(ctrls->window), TRUE);

    /* Uncertainties due to uncertainties in input values */
    ctrls->frame_unc_input = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_unc_input),
                               GTK_WIDGET(gwy_label_new_header("Uncertainties due to uncertainties in input values")));

    ctrls->table_unc_input = gtk_table_new(6, 8, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_unc_input), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_unc_input), ctrls->table_unc_input);

    i = 0;

    /* create separately for individual modes */

    switch (data->args->hertzdata.mode) {
    case R_MODE:
        /*             RMODE           */
        /*    labels   */
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(E<sub>r</sub>)/GPa"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(E<sub>it</sub>)/GPa"), 4, 5, i, i + 1, GTK_FILL, 0, 0, 0);

        i++;

        /* u(h) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(h)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uh = gtk_adjustment_new(unc->uh, 0, 1e6, 0.01, 1, 0);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uh), 0.001, 3);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->uh, "value-changed", G_CALLBACK(unc_depth_changed), data);

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

        i++;

        /* u(F) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uF = gtk_adjustment_new(unc->uF, 0, 1e6, 0.0001, 0.01, 0.);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uF), 0.0001, 4);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->uF, "value-changed", G_CALLBACK(unc_load_changed), data);

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

        i++;

        /* u(radius) */
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(radius)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uradius = gtk_adjustment_new(unc->uradius * 1e9, 0, 1e6, 1, 10, 0);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uradius), 0.1, 1);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->uradius, "value-changed", G_CALLBACK(unc_radius_changed), data);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uEruradius = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uEruradius), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEruradius, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        ctrls->uEituradius = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uEituradius), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEituradius, 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /* u(nu) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(ν)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->unu = gtk_adjustment_new(inst->unu, 0., 0.5, 0.01, 0.1, 0.);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->unu), 0.01, 2);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->unu, "value-changed", G_CALLBACK(unc_nu_changed), data);

        ctrls->uEitunu = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uEitunu), 0.0, 0.5);

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEitunu, 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /* u(nu_i) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(ν<sub>i</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->unui = gtk_adjustment_new(inst->unui, 0, 0.5, 0.01, 0.1, 0.);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->unui), 0.01, 2);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->unui, "value-changed", G_CALLBACK(unc_nui_changed), data);

        ctrls->uEitunui = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uEitunui), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEitunui, 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /* u(E_i) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(E<sub>i</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uEi = gtk_adjustment_new(inst->uEi * 1e-9, 0, 1e5, 1, 10, 0.);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uEi), 1, 0);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->uEi, "value-changed", G_CALLBACK(unc_Ei_changed), data);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uEituEi = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uEituEi), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEituEi, 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /* total */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("Total"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uErtotal = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uErtotal), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uErtotal, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        ctrls->uEittotal = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uEittotal), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEittotal, 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        break;

    case ER_MODE:
        /*             ER_MODE           */
        /*    labels   */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(radius)/nm"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

        i++;

        /* u(h) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(h)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uh = gtk_adjustment_new(unc->uh, 0, 1e6, 0.01, 1, 0);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uh), 0.001, 3);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->uh, "value-changed", G_CALLBACK(unc_depth_changed), data);

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uradiusuh = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uradiusuh), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uradiusuh, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /* u(F) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uF = gtk_adjustment_new(unc->uF, 0, 1e6, 0.0001, 0.01, 0.);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uF), 0.0001, 4);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->uF, "value-changed", G_CALLBACK(unc_load_changed), data);

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uradiusuF = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uradiusuF), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uradiusuF, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /*  u(Er) */
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(E<sub>r</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uEr = gtk_adjustment_new(unc->uEr, 0, 1e6, 1, 10, 0);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uEr), 0.1, 1);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->uEr, "value-changed", G_CALLBACK(unc_Er_changed), data);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uradiusuEr = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uradiusuEr), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uradiusuEr, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /* total */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("Total"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uradiustotal = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uradiustotal), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uradiustotal, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        break;

    case EIT_MODE:
        /*             EIT_MODE           */
        /*    labels   */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(radius)/nm"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

        i++;

        /* u(h) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(h)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uh = gtk_adjustment_new(unc->uh, 0, 1e6, 0.01, 1, 0);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uh), 0.001, 3);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->uh, "value-changed", G_CALLBACK(unc_depth_changed), data);

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uradiusuh = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uradiusuh), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uradiusuh, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /* u(F) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uF = gtk_adjustment_new(unc->uF, 0, 1e6, 0.0001, 0.01, 0.);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uF), 0.0001, 4);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->uF, "value-changed", G_CALLBACK(unc_load_changed), data);

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uradiusuF = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uradiusuF), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uradiusuF, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /*  u(Eit) */
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(E<sub>IT</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uEit = gtk_adjustment_new(unc->uEit, 0, 1e6, 1, 10, 0);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uEit), 0.1, 1);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->uEit, "value-changed", G_CALLBACK(unc_Eit_changed), data);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uradiusuEit = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uradiusuEit), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uradiusuEit, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /* u(nu) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(ν)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->unu = gtk_adjustment_new(inst->unu, 0., 0.5, 0.01, 0.1, 0.);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->unu), 0.01, 2);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->unu, "value-changed", G_CALLBACK(unc_nu_changed), data);

        ctrls->uradiusunu = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uradiusunu), 0.0, 0.5);

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uradiusunu, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /* u(nu_i) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(ν<sub>i</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->unui = gtk_adjustment_new(inst->unui, 0, 0.5, 0.01, 0.1, 0.);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->unui), 0.01, 2);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->unui, "value-changed", G_CALLBACK(unc_nui_changed), data);

        ctrls->uradiusunui = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uradiusunui), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uradiusunui, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /* u(E_i) */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(E<sub>i</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uEi = gtk_adjustment_new(inst->uEi * 1e-9, 0, 1e5, 1, 10, 0.);
        spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uEi), 1, 0);
        gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(ctrls->uEi, "value-changed", G_CALLBACK(unc_Ei_changed), data);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uradiusuEi = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uradiusuEi), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uradiusuEi, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        i++;

        /* total */

        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("Total"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

        ctrls->uradiustotal = gtk_label_new(NULL);
        gtk_misc_set_alignment(GTK_MISC(ctrls->uradiustotal), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uradiustotal, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        break;
    }

    /* Uncertainties due to contact point */

    ctrls->frame_unc_contact = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_unc_contact),
                               GTK_WIDGET(gwy_label_new_header("Uncertainties due to choice of contact point")));

    ctrls->table_unc_contact = gtk_table_new(2 * CONTACTMAX + 1 + 1, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_unc_contact), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_unc_contact), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_unc_contact), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_unc_contact), ctrls->table_unc_contact);

    unc->ERc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));

    /* table is filled in propagate_uncertainties, because it is updated at the change of R */
    hertz_propagate_uncertainties(data->controls->hertzunccontrols, &(data->args->hertzdata), &(data->args->hertzunc), &(data->args->fddata), &(data->args->instdata));

    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(hertz_uncertainty_save_data), data->args);

    /* Monte Carlo calculation */

    ctrls->frame_mc = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_mc),
                               GTK_WIDGET(gwy_label_new_header("Monte Carlo calculation of uncertainties")));

    ctrls->table_mc = gtk_table_new(6, 7, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_mc), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_mc), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_mc), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_mc), ctrls->table_mc);

    i = 0;
    gtk_table_attach(GTK_TABLE(ctrls->table_mc), gwy_label_new_header("Number of iterations"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->Nmc = gtk_adjustment_new(unc->Nmc, 1, 1e9, 10, 100, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->Nmc), 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_mc), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->Nmc, "value-changed", G_CALLBACK(mc_iterations_changed), data);

    i++;

    ctrls->mc_warn = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(ctrls->table_mc), ctrls->mc_warn, 0, 3, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    if (unc->Nmc > MC_WHEN_WARN) {
        label_set_gint_format(ctrls->mc_warn, unc->Nmc, "Warning: Monte Carlo calculation for %d iterations may take very long!");
    }

    //no Monte Carlo window without the uncertainty dialog
    data->controls->hertzmccontrols->has_window = FALSE;
    ctrls->button_mc = gtk_button_new_with_mnemonic("_Monte Carlo");
    gtk_button_set_image(GTK_BUTTON(ctrls->button_mc), gtk_image_new_from_stock(GWY_STOCK_POLYNOM, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_mc, "clicked", G_CALLBACK(hertz_uncertainty_montecarlo), data);

    ctrls->vbox = gtk_vbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_unc_input, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_unc_contact, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_save,  FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_mc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_mc, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->vbox);
    gtk_widget_show_all(ctrls->window);
}

static void hertz_mc_hist_Er(Data *data)
{
    show_histogram("E_r", "GPa", data->args->hertzmc.mcxhist[0], data->args->hertzmc.mcyhist[0], data->args->hertzmc.mcnstat, data->controls->hertzmccontrols->window);
}

static void hertz_mc_hist_Eit(Data *data)
{
    show_histogram("E_IT", "GPa", data->args->hertzmc.mcxhist[1], data->args->hertzmc.mcyhist[1], data->args->hertzmc.mcnstat, data->controls->hertzmccontrols->window);
}

static void hertz_mc_hist_radius(Data *data)
{
    show_histogram("R", "nm", data->args->hertzmc.mcxhist[2], data->args->hertzmc.mcyhist[2], data->args->hertzmc.mcnstat, data->controls->hertzmccontrols->window);
}

static void hertz_mc_save_data(Args *args)
{
    gchar *buffer = NULL;

    /* add check for mc results */
    // if (hertz_mc_has_results(&(args->hertzmc))) {
    buffer = hertz_mc_export_data(&(args->hertzdata), &(args->hertzunc), &(args->hertzmc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
	// }
}

void hertz_mc_close_window(Data *data)
{

    gint j;

    if (data->controls->hertzmccontrols->has_window) {
        gtk_widget_destroy(data->controls->hertzmccontrols->window);

        destroy_HertzMCdata(&(data->args->hertzmc));
        data->controls->hertzmccontrols->has_window = FALSE;

        for (j = 0; j < NPARAM_HZ; j++) {
            if (GTK_IS_WIDGET(data->controls->hertzmccontrols->hist_window[j])) {
                gtk_widget_destroy(data->controls->hertzmccontrols->hist_window[j]);
            }
        }

        g_free(data->controls->hertzmccontrols->hist_window);
        g_free(data->controls->hertzmccontrols->button_hist);
    }
}

void hertz_uncertainty_montecarlo(Data *data)
{
    gint i, j;
    GtkWidget *label;
    gchar str[200];
    GdkColor color;

    HertzMCControls *ctrls;
    HertzMCdata *hzmc;

    gdk_color_parse("red", &color);

    //if there was a previous window, close it and any histograms, free data
    hertz_mc_close_window(data);

    hzmc = &(data->args->hertzmc);

    // initialize for number of iterations which was chosen in the Uncertainties dialog
    init_HertzMCdata(hzmc, data->args->hertzunc.Nmc);

    // run calculation
    hertz_uncertainty_montecarlo_run_calc(&(data->args->hertzdata), &(data->args->hertzunc), &(data->args->hertzmc),
					  &(data->args->fddata), &(data->args->instdata));

    //create window and display results
    ctrls = data->controls->hertzmccontrols;

    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Monte Carlo calculation of uncertainties in Hertz analysis");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(hertz_mc_close_window), data);
    gtk_window_set_transient_for(GTK_WINDOW(ctrls->window), GTK_WINDOW(data->controls->maincontrols->window_main));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(ctrls->window), TRUE);

    ctrls->frame_input = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_input),
                               GTK_WIDGET(gwy_label_new_header("Input for Monte Carlo simulation")));

    ctrls->table_input = gtk_table_new(3, 4, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_input), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_input), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_input), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_input), ctrls->table_input);

    i = 0;

    /*          u(h)             */

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("u(h)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    g_snprintf(str, sizeof(str), "%g", data->args->hertzunc.uh);
    gtk_label_set_text(GTK_LABEL(label), str);
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /*          u(F)             */

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    g_snprintf(str, sizeof(str), "%g", data->args->hertzunc.uF);
    gtk_label_set_text(GTK_LABEL(label), str);
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /*          Number of iteration */

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("Number of iterations "), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    g_snprintf(str, sizeof(str), "%d", hzmc->Nmc);
    gtk_label_set_text(GTK_LABEL(label), str);
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    // Results
    ctrls->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_results),
                               GTK_WIDGET(gwy_label_new_header("Results of Monte Carlo simulation")));

    ctrls->table_results = gtk_table_new(NPARAM_HZ + 1, 4, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_results), ctrls->table_results);

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("mean"), 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("stdev"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    ctrls->button_hist = (GtkWidget **)g_malloc(NPARAM_HZ * sizeof(GtkWidget *));
    ctrls->hist_window = (GtkWidget **)g_malloc(NPARAM_HZ * sizeof(GtkWidget *));

    for (j = 0; j < NPARAM_HZ; j++) {
        ctrls->hist_window[j] = NULL;
    }

    if (data->args->hertzdata.mode == R_MODE) {
        /*                       E_r           */

        gtk_table_attach(GTK_TABLE(ctrls->table_results),
                         gwy_label_new_header("E<sub>r</sub>/GPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

        label = gtk_label_new(NULL);
        g_snprintf(str, sizeof(str), MODUL, hzmc->mcavg[i - 1]);
        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        label = gtk_label_new(NULL);
        g_snprintf(str, sizeof(str), MODUL, hzmc->mcstd[i - 1]);
        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        ctrls->button_hist[i - 1] = gtk_button_new();
        gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
        g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(hertz_mc_hist_Er), data);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
        i++;

        /*                       E_IT          */

        gtk_table_attach(GTK_TABLE(ctrls->table_results),
                         gwy_label_new_header("E<sub>it</sub>/GPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

        label = gtk_label_new(NULL);
        g_snprintf(str, sizeof(str), MODUL, hzmc->mcavg[i - 1]);
        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        label = gtk_label_new(NULL);
        g_snprintf(str, sizeof(str), MODUL, hzmc->mcstd[i - 1]);
        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        ctrls->button_hist[i - 1] = gtk_button_new();
        gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
        g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(hertz_mc_hist_Eit), data);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
        i++;
    }
    else {
        /*                       radius                */

        i = 3;
        gtk_table_attach(GTK_TABLE(ctrls->table_results),
                         gwy_label_new_header("radius/nm"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

        label = gtk_label_new(NULL);
        g_snprintf(str, sizeof(str), DEPTH, hzmc->mcavg[i - 1]);
        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        label = gtk_label_new(NULL);
        g_snprintf(str, sizeof(str), DEPTH, hzmc->mcstd[i - 1]);
        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        ctrls->button_hist[i - 1] = gtk_button_new();
        gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
        g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(hertz_mc_hist_radius), data);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
        i++;
    }

    if (verbose) {
        g_print(" Nmc %d \n", hzmc->Nmc);
    }

    if (100.*data->args->hertzmc.skipmc / hzmc->Nmc > MAXBAD) {
        label = gtk_label_new(NULL);
        g_snprintf(str, sizeof(str), "Warning: Approx %.1f %% of the results had to be removed!", 100.*data->args->hertzmc.skipmc / hzmc->Nmc);
        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_widget_modify_fg(label, GTK_STATE_NORMAL, &color);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 0, 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    }

    // save data
    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(hertz_mc_save_data), data->args);

    ctrls->vbox = gtk_vbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_input, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_results, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_save, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->vbox);
    gtk_widget_show_all(ctrls->window);
    ctrls->has_window = TRUE;
}

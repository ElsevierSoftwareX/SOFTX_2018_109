#include "tool_op_unc_gui.h"
#include "tool_op_unc.h"

#include "controls.h"
#include "datatypes.h"
#include "gui-utils.h"
#include "niget-common.h"

#include "niget-gtk.h" /* controls definitions */

#include <math.h>

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

#define MC_WHEN_WARN 1e5

static void op_mc_hist_hc(Data *data);
static void op_mc_hist_Aphc(Data *data);
static void op_mc_hist_Hit(Data *data);
static void op_mc_hist_Er(Data *data);
static void op_mc_hist_Eit(Data *data);
static void op_mc_hist_S(Data *data);

static void op_uncertainty_save_data(Args *args);
static void op_mc_save_data(Args *args);


void op_propagate_uncertainties(OPUncControls *ctrls, const OPdata *opdata, OPUncdata *opunc,
				const FDdata *fddata, const Instdata *instdata, const Area *area)
{
    if (!opdata->has_fit_S) {
        return;
    }

    if (opunc->Sc == NULL || opunc->Ec == NULL || opunc->Hc == NULL) {
        return;
    }

    // calculate the uncertainties
    op_propagate_uncertainties_calc(opdata, opunc, fddata, instdata, area);

    /* show */
    label_set_gdouble_format(ctrls->uEitunu, opunc->uEitunu, MODUL);
    label_set_gdouble_format(ctrls->uEitunui, opunc->uEitunui, MODUL);
    label_set_gdouble_format(ctrls->uEituEi, opunc->uEituEi, MODUL);

    // total uncertainties
    opunc->uS = sqrt(sq(opunc->uSuh) + sq(opunc->uSuF));
    label_set_gdouble_format(ctrls->uS, opunc->uS, SLOPE);

    opunc->uhc = sqrt(sq(opunc->uhcuh) + sq(opunc->uhcuF));
    label_set_gdouble_format(ctrls->uhc, opunc->uhc, DEPTH);

    opunc->uA = sqrt(sq(opunc->uAuh) + sq(opunc->uAuF));
    label_set_gdouble_format(ctrls->uA, opunc->uA, AREA);

    opunc->uHit = sqrt(sq(opunc->uHituh) + sq(opunc->uHituF));
    label_set_gdouble_format(ctrls->uHit, opunc->uHit, HARD);

    opunc->uEr = sqrt(sq(opunc->uEruh) + sq(opunc->uEruF));
    label_set_gdouble_format(ctrls->uEr, opunc->uEr, MODUL);

    opunc->uEit = sqrt(sq(opunc->uEituh) + sq(opunc->uEituF) + sq(opunc->uEitunu) + sq(opunc->uEitunui) + sq(opunc->uEituEi));
    label_set_gdouble_format(ctrls->uEit, opunc->uEit, MODUL);

    //table_contact doesn't change with any of the adjustables nu, nui or Ei, uh, uF, unu, uniu, uEi
    //no need to update
}

static void mc_iterations_changed(Data *data)
{
    data->args->opunc.Nmc = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opunccontrols->Nmc));

    if (data->args->opunc.Nmc > MC_WHEN_WARN) {
        label_set_gint_format(data->controls->opunccontrols->mc_warn, data->args->opunc.Nmc,
                              "Warning: Monte Carlo calculation for %d iterations may take very long!");
    }
    else {
        label_clear(data->controls->opunccontrols->mc_warn);
    }
}

static void unc_depth_changed(Data *data)
{
    data->args->opunc.uh = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opunccontrols->uh));
    op_propagate_uncertainties(data->controls->opunccontrols, &(data->args->opdata), &(data->args->opunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

static void unc_load_changed(Data *data)
{
    data->args->opunc.uF = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opunccontrols->uF));
    op_propagate_uncertainties(data->controls->opunccontrols, &(data->args->opdata), &(data->args->opunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

static void unc_nu_changed(Data *data)
{
    data->args->instdata.unu = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opunccontrols->unu));
    op_propagate_uncertainties(data->controls->opunccontrols, &(data->args->opdata), &(data->args->opunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

static void unc_nui_changed(Data *data)
{
    data->args->instdata.unui = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opunccontrols->unui));
    op_propagate_uncertainties(data->controls->opunccontrols, &(data->args->opdata), &(data->args->opunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

static void unc_Ei_changed(Data *data)
{
    data->args->instdata.uEi = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opunccontrols->uEi)) * 1e9;
    op_propagate_uncertainties(data->controls->opunccontrols, &(data->args->opdata), &(data->args->opunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

static void op_uncertainty_save_data(Args *args)
{
    gchar *buffer = NULL;
    
    if (op_unc_has_results(&(args->opunc))) {
	buffer = op_uncertainty_export_data(&(args->opdata), &(args->opunc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

void op_uncertainty(Data *data)
{
    GtkWidget *label, *spin;
    gint i, j;
    OPUncControls *ctrls;
    OPUncdata *unc;
    Instdata *inst;
    gchar str[50];

    if ((data->args->opunc.Sc != NULL) && (data->args->opunc.Ec != NULL) && (data->args->opunc.Hc != NULL)) {
        gtk_window_present(GTK_WINDOW(data->controls->opunccontrols->window));
        return;
    }

    unc = &(data->args->opunc);
    ctrls = data->controls->opunccontrols;
    inst = &(data->args->instdata);

    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Uncertainty in Oliver Pharr analysis");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(op_unc_close_window), data);
    gtk_window_set_transient_for(GTK_WINDOW(ctrls->window), GTK_WINDOW(data->controls->maincontrols->window_main));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(ctrls->window), TRUE);

    /* Uncertainties due to uncertainties in input values */

    ctrls->frame_unc_input = gtk_frame_new(NULL);
    //	gtk_frame_set_label_widget (GTK_FRAME (ctrls->frame_unc_input),
    //				    GTK_WIDGET (gwy_label_new_header("Uncertainties due to uncertainties in input values")));
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_unc_input),
                               GTK_WIDGET(gwy_label_new_header("Uncertainties in input values")));

    ctrls->table_unc_input = gtk_table_new(6, 10, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_unc_input), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_unc_input), ctrls->table_unc_input);

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(S/(mN/nm)"), 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(h<sub>c</sub>)/nm"), 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(A<sub>p</sub>)/nm<sup>2</sup>"), 5, 6, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(H<sub>it</sub>)/MPa"), 6, 7, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(E<sub>r</sub>)/GPa"), 7, 8, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(E<sub>it</sub>)/GPa"), 8, 9, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    i++;

    /* u(h) */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(h)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->uh = gtk_adjustment_new(unc->uh, 0, 1e6, 0.01, 1, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uh), 0.001, 3);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->uh, "value-changed", G_CALLBACK(unc_depth_changed), data);

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /* u(F) */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->uF = gtk_adjustment_new(unc->uF, 0, 1e6, 0.0001, 0.01, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uF), 0.0001, 4);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->uF, "value-changed", G_CALLBACK(unc_load_changed), data);

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /* u(nu) */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(ν)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->unu = gtk_adjustment_new(inst->unu, 0, 0.5, 0.01, 0.1, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->unu), 0.01, 2);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->unu, "value-changed", G_CALLBACK(unc_nu_changed), data);

    ctrls->uEitunu = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEitunu), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEitunu, 8, 9, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

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
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEitunui, 8, 9, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    i++;

    /* u(E_i) */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(E<sub>i</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->uEi = gtk_adjustment_new(inst->uEi * 1e-9, 0, 1e5, 1, 10, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uEi), 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->uEi, "value-changed", G_CALLBACK(unc_Ei_changed), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->uEituEi = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEituEi), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEituEi, 8, 9, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    i++;

    /* total */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("Total"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->uS = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uS), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uS, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uhc = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uhc), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uhc, 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uA = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uA), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uA, 5, 6, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uHit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uHit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uHit, 6, 7, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uEr = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEr), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEr, 7, 8, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uEit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEit, 8, 9, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    /* Uncertainties due to contact point */

    ctrls->frame_unc_contact = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_unc_contact),
                               GTK_WIDGET(gwy_label_new_header("Uncertainties due to choice of contact point")));

    ctrls->table_unc_contact = gtk_table_new(2 * CONTACTMAX + 1 + 1, 4, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_unc_contact), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_unc_contact), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_unc_contact), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_unc_contact), ctrls->table_unc_contact);

    /* table headers */

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("Contact point shift"), 0, 1, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("S/(mN/nm)"), 1, 2, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("E<sub>r</sub>/GPa"), 2, 3, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("H<sub>it</sub>/GPa"), 3, 4, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    // gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
    // gwy_label_new_header("weight"), 2, 3, i,i+1,  GTK_EXPAND | GTK_FILL, GTK_FILL, 0, 0);

    i++;

    // Ec is Er, Hc is Hit, S slope powerlaw
    unc->Sc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
    unc->Ec = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
    unc->Hc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));

    /* values */

    for (j = -CONTACTMAX; j <= CONTACTMAX; j++)	{

        /* vertical alignments: j: 1.0 (bottom), others: 0.0 (top); looks quite nice, still probably incorrect */

        label = gtk_label_new(NULL);
        label_set_gint_format(label, j, "%d\n");
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 1.0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 0, 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        op_fit_shift_contact(&(data->args->opdata), &(data->args->fddata), &(data->args->instdata), &(data->args->area), unc->Sc, unc->Ec, unc->Hc, j);

        //fill in S
        label = gtk_label_new(NULL);

        if (unc->Sc[j + CONTACTMAX] > 0) {
            g_snprintf(str, sizeof(str), SLOPE, unc->Sc[j + CONTACTMAX]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);


        //fill in Er
        label = gtk_label_new(NULL);

        if (unc->Ec[j + CONTACTMAX] > 0) {
            g_snprintf(str, sizeof(str), MODUL, unc->Ec[j + CONTACTMAX]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        /*
        label = gtk_label_new(NULL);
        g_snprintf(str, sizeof(str),"1/%d", 2*CONTACTMAX+1);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 2, 3, i, i+1, GTK_FILL, GTK_FILL, 0, 0);
        */

        //fill in Hit
        label = gtk_label_new(NULL);

        if (unc->Hc[j + CONTACTMAX] > 0) {
            g_snprintf(str, sizeof(str), HARD, unc->Hc[j + CONTACTMAX]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        i++;
    }

    op_propagate_uncertainties(data->controls->opunccontrols, &(data->args->opdata), &(data->args->opunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));

    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(op_uncertainty_save_data), data->args);

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

    // no Monte Carlo window without the uncertainty dialog
    data->controls->opmccontrols->has_window = FALSE;
    ctrls->button_mc = gtk_button_new_with_mnemonic("_Monte Carlo");
    gtk_button_set_image(GTK_BUTTON(ctrls->button_mc), gtk_image_new_from_stock(GWY_STOCK_POLYNOM, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_mc, "clicked", G_CALLBACK(op_uncertainty_montecarlo), data);

    ctrls->vbox = gtk_vbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_unc_input, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_unc_contact, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_save, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_mc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_mc, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->vbox);
    gtk_widget_show_all(ctrls->window);
}

void op_unc_close_window(Data *data)
{
    g_free(data->args->opunc.Sc);
    g_free(data->args->opunc.Ec);
    g_free(data->args->opunc.Hc);
    data->args->opunc.Sc = NULL;
    data->args->opunc.Ec = NULL;
    data->args->opunc.Hc = NULL;
    op_mc_close_window(data);
    gtk_widget_destroy(data->controls->opunccontrols->window);
}

void op_uncertainty_montecarlo(Data *data)
{
    gint i, j;
    GtkWidget *label;
    GdkColor color;

    OPMCControls *ctrls;
    OPMCdata *opmc;

    gdk_color_parse("red", &color);

    //if there was a previous window, close it and any histograms, free data
    op_mc_close_window(data);

    opmc = &(data->args->opmc);

    // initialize for number of iterations which was chosen in the Uncertainties dialog
    init_OPMCdata(opmc, data->args->opunc.Nmc);

    // run calculation
    op_uncertainty_montecarlo_run_calc(&(data->args->opdata), &(data->args->opunc), &(data->args->opmc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));


    //create window and display results
    ctrls = data->controls->opmccontrols;

    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);

    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Monte Carlo calculation of uncertainties in Oliver Pharr analysis");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(op_mc_close_window), data);
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
    label_set_gdouble_format(label, data->args->opunc.uh, "%g");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /*          u(F)             */

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opunc.uF, "%g");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /*          Number of iteration */

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("Number of iterations "), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gint_format(label, data->args->opunc.Nmc, "%d");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    // Results
    ctrls->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_results),
                               GTK_WIDGET(gwy_label_new_header("Results of Monte Carlo simulation")));

    ctrls->table_results = gtk_table_new(NPARAM_OP + 1, 4, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_results), ctrls->table_results);

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("mean"), 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("stdev"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    ctrls->button_hist = (GtkWidget **)g_malloc(NPARAM_OP * sizeof(GtkWidget *));
    ctrls->hist_window = (GtkWidget **)g_malloc(NPARAM_OP * sizeof(GtkWidget *));

    for (j = 0; j < NPARAM_OP; j++) {
        ctrls->hist_window[j] = NULL;
    }

    /*                       h_c           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("h<sub>c</sub>/nm"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcavg[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcstd[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_mc_hist_hc), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       A_p(h_c)           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("A<sub>p</sub>/nm<sup>2</sup>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcavg[i - 1], AREA);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcstd[i - 1], AREA);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_mc_hist_Aphc), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       H_IT           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("H<sub>it</sub>/MPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcavg[i - 1], HARD);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcstd[i - 1], HARD);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_mc_hist_Hit), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       E_r           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("E<sub>r</sub>/GPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcavg[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcstd[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_mc_hist_Er), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       E_IT          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("E<sub>it</sub>/GPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcavg[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcstd[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_mc_hist_Eit), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       S           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("S/(mN/nm)"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcavg[i - 1], "%g ");
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opmc.mcstd[i - 1], SLOPE);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_mc_hist_S), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    if (100.*data->args->opmc.skipmc / opmc->Nmc > MAXBAD) {
        label = gtk_label_new(NULL);
        label_set_gdouble_format(label, 100.*data->args->opmc.skipmc / opmc->Nmc, "Warning: Approx %.1f %% of the results had to be removed!");
        gtk_widget_modify_fg(label, GTK_STATE_NORMAL, &color);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 0, 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    }

    // save data
    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(op_mc_save_data), data->args);

    ctrls->vbox = gtk_vbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_input, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_results, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_save, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->vbox);
    gtk_widget_show_all(ctrls->window);
    ctrls->has_window = TRUE;
}

void op_mc_close_window(Data *data)
{
    gint j;

    if (data->controls->opmccontrols->has_window) {
        gtk_widget_destroy(data->controls->opmccontrols->window);

        destroy_OPMCdata(&(data->args->opmc));
        data->controls->opmccontrols->has_window = FALSE;

        for (j = 0; j < NPARAM_OP; j++) {
            if (GTK_IS_WIDGET(data->controls->opmccontrols->hist_window[j])) {
                gtk_widget_destroy(data->controls->opmccontrols->hist_window[j]);
            }
        }

        g_free(data->controls->opmccontrols->hist_window);
        g_free(data->controls->opmccontrols->button_hist);
    }
}

static void op_mc_save_data(Args *args)
{
    gchar *buffer = NULL;

    /* ADD CHECK FOR EXISTING RESULTS */
    buffer = op_mc_export_data(&(args->opdata), &(args->opunc), &(args->opmc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN);
    buffer_to_file_dialog(buffer);
    g_free(buffer);
}

static void op_mc_hist_hc(Data *data)
{
    show_histogram("h_c", "nm", data->args->opmc.mcxhist[0], data->args->opmc.mcyhist[0], data->args->opmc.mcnstat, data->controls->opmccontrols->window);
}

static void op_mc_hist_Aphc(Data *data)
{
    show_histogram("A_p(h_c)", "nm^2", data->args->opmc.mcxhist[1], data->args->opmc.mcyhist[1], data->args->opmc.mcnstat, data->controls->opmccontrols->window);
}

static void op_mc_hist_Hit(Data *data)
{
    show_histogram("H_IT", "MPa", data->args->opmc.mcxhist[2], data->args->opmc.mcyhist[2], data->args->opmc.mcnstat, data->controls->opmccontrols->window);
}

static void op_mc_hist_Er(Data *data)
{
    show_histogram("E_r", "GPa", data->args->opmc.mcxhist[3], data->args->opmc.mcyhist[3], data->args->opmc.mcnstat, data->controls->opmccontrols->window);
}

static void op_mc_hist_Eit(Data *data)
{
    show_histogram("E_IT", "GPa", data->args->opmc.mcxhist[4], data->args->opmc.mcyhist[4], data->args->opmc.mcnstat, data->controls->opmccontrols->window);
}

static void op_mc_hist_S(Data *data)
{
    show_histogram("S", "mN/nm", data->args->opmc.mcxhist[5], data->args->opmc.mcyhist[5], data->args->opmc.mcnstat, data->controls->opmccontrols->window);
}

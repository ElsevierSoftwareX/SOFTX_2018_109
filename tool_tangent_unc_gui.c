#include "tool_tangent_unc_gui.h"
#include "tool_tangent_unc.h"

#include "controls.h"
#include "gui-utils.h"
#include "niget-common.h"

#include "niget-gtk.h" /* controls definitions */

#include <gtk/gtk.h>

#define MC_WHEN_WARN 1e6

static void mc_iterations_changed(Data *data);
static void unc_depth_changed(Data *data);
static void unc_load_changed(Data *data);
static void unc_nu_changed(Data *data);
static void unc_nui_changed(Data *data);
static void unc_Ei_changed(Data *data);

static void tangent_uncertainty_save_data(Args *args);

static void tangent_mc_hist_hc(Data *data);
static void tangent_mc_hist_Aphc(Data *data);
static void tangent_mc_hist_Hit(Data *data);
static void tangent_mc_hist_Er(Data *data);
static void tangent_mc_hist_Eit(Data *data);
static void tangent_mc_hist_S(Data *data);


void tangent_propagate_uncertainties(TangentUncControls *ctrls, Tangentdata *tgdata, TangentUncdata *unc, FDdata *fddata, Instdata *instdata, Area *area)
{
    if (!tgdata->has_fit) {
        return;
    }

    if (unc->Sc == NULL || unc->Ec == NULL || unc->Hc == NULL) {
        return;
    }

    // calculate the uncertainties
    tangent_propagate_uncertainties_calc(tgdata, unc, fddata, instdata, area);

    /* show */
    // total uncertainties
    label_set_gdouble_format(ctrls->uS, unc->uS, SLOPE);
    label_set_gdouble_format(ctrls->uhc, unc->uhc, DEPTH);
    label_set_gdouble_format(ctrls->uA, unc->uA, AREA);
    label_set_gdouble_format(ctrls->uHit, unc->uHit, HARD);
    label_set_gdouble_format(ctrls->uEr, unc->uEr, MODUL);
    label_set_gdouble_format(ctrls->uEit, unc->uEit, MODUL);

    // other contributions
    label_set_gdouble_format(ctrls->uEitunu, unc->uEitunu, MODUL);
    label_set_gdouble_format(ctrls->uEitunui, unc->uEitunui, MODUL);
    label_set_gdouble_format(ctrls->uEituEi, unc->uEituEi, MODUL);

    /*table_contact doesn't change with any of the adjustables nu, nui or Ei, uh, uF, unu, uniu, uEi
    no need to update */
}

static void mc_iterations_changed(Data *data)
{
    data->args->tgunc.Nmc = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->tgunccontrols->Nmc));

    if (data->args->tgunc.Nmc > MC_WHEN_WARN)
        label_set_gint_format(data->controls->tgunccontrols->mc_warn, data->args->tgunc.Nmc,
                              "Warning: Monte Carlo calculation for %d iterations may take very long!");
    else {
        label_clear(data->controls->tgunccontrols->mc_warn);
    }
}

static void unc_depth_changed(Data *data)
{
    data->args->tgunc.uh = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->tgunccontrols->uh));
    tangent_propagate_uncertainties(data->controls->tgunccontrols, &(data->args->tgdata), &(data->args->tgunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

static void unc_load_changed(Data *data)
{
    data->args->tgunc.uF = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->tgunccontrols->uF));
    tangent_propagate_uncertainties(data->controls->tgunccontrols, &(data->args->tgdata), &(data->args->tgunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

static void unc_nu_changed(Data *data)
{
    data->args->instdata.unu = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->tgunccontrols->unu));
    tangent_propagate_uncertainties(data->controls->tgunccontrols, &(data->args->tgdata), &(data->args->tgunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

static void unc_nui_changed(Data *data)
{
    data->args->instdata.unui = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->tgunccontrols->unui));
    tangent_propagate_uncertainties(data->controls->tgunccontrols, &(data->args->tgdata), &(data->args->tgunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

static void unc_Ei_changed(Data *data)
{
    data->args->instdata.uEi = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->tgunccontrols->uEi)) * 1e9;
    tangent_propagate_uncertainties(data->controls->tgunccontrols, &(data->args->tgdata), &(data->args->tgunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));
}

static void tangent_uncertainty_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (tangent_unc_has_results(&(args->tgunc))) {
	buffer = tangent_uncertainty_export_data(&(args->tgdata), &(args->tgunc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
    }
}

void tangent_unc_close_window(Data *data)
{
    g_free(data->args->tgunc.Sc);
    g_free(data->args->tgunc.Ec);
    g_free(data->args->tgunc.Hc);
    data->args->tgunc.Sc = NULL;
    data->args->tgunc.Ec = NULL;
    data->args->tgunc.Hc = NULL;
    tangent_mc_close_window(data);
    gtk_widget_destroy(data->controls->tgunccontrols->window);
}

void tangent_uncertainty(Data *data)
{
    GtkWidget *label, *spin;
    gint i, j, k;
    TangentUncControls *ctrls;
    TangentUncdata *unc;
    Instdata *inst;
    gchar str[50];

    if ((data->args->tgunc.Sc != NULL) && (data->args->tgunc.Ec != NULL) && (data->args->tgunc.Hc != NULL)) {
        gtk_window_present(GTK_WINDOW(data->controls->tgunccontrols->window));
        return;
    }

    unc = &(data->args->tgunc);
    ctrls =	data->controls->tgunccontrols;
    inst = &(data->args->instdata);


    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Uncertainty in tangent analysis");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(tangent_unc_close_window), data);
    gtk_window_set_transient_for(GTK_WINDOW(ctrls->window), GTK_WINDOW(data->controls->maincontrols->window_main));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(ctrls->window), TRUE);


    /* Uncertainties due to uncertainties in input values */

    ctrls->frame_unc_input = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_unc_input),
                               GTK_WIDGET(gwy_label_new_header("Uncertainties due to uncertainties in input values")));

    ctrls->table_unc_input = gtk_table_new(6, 10, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_unc_input), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_unc_input), ctrls->table_unc_input);

    i = 0;
    j = 3;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(S)/mN.nm<sup>-1</sup>"), j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    j++;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(h<sub>c</sub>)/nm"), j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    j++;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(A<sub>p</sub>)/nm<sup>2</sup>"), j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    j++;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(H<sub>it</sub>)/MPa"), j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    j++;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(E<sub>r</sub>)/GPa"), j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    j++;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(E<sub>it</sub>)/GPa"), j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

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

    j = 3;
    ctrls->uS = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uS), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uS, j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    j++;

    ctrls->uhc = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uhc), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uhc, j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    j++;

    ctrls->uA = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uA), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uA, j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    j++;

    ctrls->uHit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uHit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uHit, j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    j++;

    ctrls->uEr = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEr), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEr, j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    j++;

    ctrls->uEit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEit, j, j + 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);


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
    j = 0;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("Contact point shift"), j, j + 1, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    j++;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("S/mN.nm<sup>-1</sup>"), j, j + 1, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    j++;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("E<sub>r</sub>/GPa"), j, j + 1, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    j++;
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("H<sub>it</sub>/GPa"), j, j + 1, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    j++;
    // gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
    // gwy_label_new_header("weight"), 2, 3, i,i+1,  GTK_EXPAND | GTK_FILL, GTK_FILL, 0, 0);

    i++;

    // Ec is Er, Hc is Hit, Sc is S
    unc->Ec = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
    unc->Hc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
    unc->Sc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));

    /* values */

    for (j = -CONTACTMAX; j <= CONTACTMAX; j++)	{

        /* vertical alignments: j: 1.0 (bottom), others: 0.0 (top); looks quite nice, still probably incorrect */

        k = 0;
        label = gtk_label_new(NULL);
        label_set_gint_format(label, j, "%d\n");
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 1.0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, k, k + 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
        k++;

        tangent_fit_shift_contact(&(data->args->tgdata), &(data->args->fddata), &(data->args->instdata), &(data->args->area), unc->Sc, unc->Ec, unc->Hc, j);

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
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, k, k + 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
        k++;


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
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, k, k + 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
        k++;

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
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, k, k + 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        i++;

    }

    tangent_propagate_uncertainties(data->controls->tgunccontrols, &(data->args->tgdata), &(data->args->tgunc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));

    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(tangent_uncertainty_save_data), data->args);

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
    data->controls->tgmccontrols->has_window = FALSE;
    ctrls->button_mc = gtk_button_new_with_mnemonic("_Monte Carlo");
    gtk_button_set_image(GTK_BUTTON(ctrls->button_mc), gtk_image_new_from_stock(GWY_STOCK_POLYNOM, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_mc, "clicked", G_CALLBACK(tangent_uncertainty_montecarlo), data);

    ctrls->vbox = gtk_vbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_unc_input, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_unc_contact, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_save, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_mc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_mc, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->vbox);
    gtk_widget_show_all(ctrls->window);
}

static void tangent_mc_hist_hc(Data *data)
{
    show_histogram("h_c", "nm", data->args->tgmc.mcxhist[0], data->args->tgmc.mcyhist[0], data->args->tgmc.mcnstat, data->controls->tgmccontrols->window);
}

static void tangent_mc_hist_Aphc(Data *data)
{
    show_histogram("A_p(h_c)", "nm^2", data->args->tgmc.mcxhist[1], data->args->tgmc.mcyhist[1], data->args->tgmc.mcnstat, data->controls->tgmccontrols->window);
}

static void tangent_mc_hist_Hit(Data *data)
{
    show_histogram("H_IT", "MPa", data->args->tgmc.mcxhist[2], data->args->tgmc.mcyhist[2], data->args->tgmc.mcnstat, data->controls->tgmccontrols->window);
}

static void tangent_mc_hist_Er(Data *data)
{
    show_histogram("E_r", "GPa", data->args->tgmc.mcxhist[3], data->args->tgmc.mcyhist[3], data->args->tgmc.mcnstat, data->controls->tgmccontrols->window);
}

static void tangent_mc_hist_Eit(Data *data)
{
    show_histogram("E_IT", "GPa", data->args->tgmc.mcxhist[4], data->args->tgmc.mcyhist[4], data->args->tgmc.mcnstat, data->controls->tgmccontrols->window);
}

static void tangent_mc_hist_S(Data *data)
{
    show_histogram("S", "(mN/nm)", data->args->tgmc.mcxhist[5], data->args->tgmc.mcyhist[5], data->args->tgmc.mcnstat, data->controls->tgmccontrols->window);
}

static void tangent_mc_save_data(Args *args)
{
    gchar *buffer = NULL;

    /* add check for mc results */
    // if (tangent_mc_has_results(&(args->tgmc))) {
    buffer = tangent_mc_export_data(&(args->tgdata), &(args->tgunc), &(args->tgmc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN);
	buffer_to_file_dialog(buffer);
	g_free(buffer);
	// }
}

void tangent_mc_close_window(Data *data)
{
    gint j;

    if (data->controls->tgmccontrols->has_window) {
        gtk_widget_destroy(data->controls->tgmccontrols->window);

        destroy_TangentMCdata(&(data->args->tgmc));
        data->controls->tgmccontrols->has_window = FALSE;

        for (j = 0; j < NPARAM_TG; j++) {
            if (GTK_IS_WIDGET(data->controls->tgmccontrols->hist_window[j])) {
                gtk_widget_destroy(data->controls->tgmccontrols->hist_window[j]);
            }
        }

        g_free(data->controls->tgmccontrols->hist_window);
        g_free(data->controls->tgmccontrols->button_hist);
    }
}

void tangent_uncertainty_montecarlo(Data *data)
{
    gint i, j;
    GtkWidget *label;
    GdkColor color;

    TangentMCControls *ctrls;
    TangentMCdata *tgmc;

    gdk_color_parse("red", &color);

    //if there was a previous window, close it and any histograms, free data
    tangent_mc_close_window(data);

    tgmc = &(data->args->tgmc);

    // initialize for number of iterations which was chosen in the Uncertainties dialog
    init_TangentMCdata(tgmc, data->args->tgunc.Nmc);

    // run calculation
    tangent_uncertainty_montecarlo_run_calc(&(data->args->tgdata), &(data->args->tgunc), &(data->args->tgmc), &(data->args->fddata), &(data->args->instdata), &(data->args->area));

    //create window and display results
    ctrls = data->controls->tgmccontrols;

    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Monte Carlo calculation of uncertainties in tangent analysis");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(tangent_mc_close_window), data);
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
    label_set_gdouble_format(label, data->args->tgunc.uh, "%g");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /*          u(F)             */

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgunc.uF, "%g");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /*          Number of iteration */

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("Number of iterations "), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gint_format(label, tgmc->Nmc, "%d");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    // Results
    ctrls->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_results),
                               GTK_WIDGET(gwy_label_new_header("Results of Monte Carlo simulation")));

    ctrls->table_results = gtk_table_new(NPARAM_TG + 1, 4, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_results), ctrls->table_results);

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("mean"), 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("stdev"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    ctrls->button_hist = (GtkWidget **)g_malloc(NPARAM_TG * sizeof(GtkWidget *));
    ctrls->hist_window = (GtkWidget **)g_malloc(NPARAM_TG * sizeof(GtkWidget *));

    for (j = 0; j < NPARAM_TG; j++) {
        ctrls->hist_window[j] = NULL;
    }

    /*                       h_c           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("h<sub>c</sub>/nm"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcavg[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcstd[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(tangent_mc_hist_hc), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       A_p(h_c)           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("A<sub>p</sub>/nm<sup>2</sup>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcavg[i - 1], AREA);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcstd[i - 1], AREA);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(tangent_mc_hist_Aphc), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       H_IT           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("H<sub>it</sub>/MPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcavg[i - 1], HARD);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcstd[i - 1], HARD);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(tangent_mc_hist_Hit), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       E_r           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("E<sub>r</sub>/GPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcavg[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcstd[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(tangent_mc_hist_Er), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       E_IT          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("E<sub>it</sub>/GPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcavg[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcstd[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(tangent_mc_hist_Eit), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       S          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("S/(mN/nm)"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcavg[i - 1], SLOPE);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->tgmc.mcstd[i - 1], SLOPE);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(tangent_mc_hist_S), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    if (100.*data->args->tgmc.skipmc / tgmc->Nmc > MAXBAD) {
        label = gtk_label_new(NULL);
        label_set_gdouble_format(label, 100. * data->args->tgmc.skipmc / tgmc->Nmc, "Warning: Approx %.1f %% of the results had to be removed!");
        gtk_widget_modify_fg(label, GTK_STATE_NORMAL, &color);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 0, 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    }

    // save data
    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(tangent_mc_save_data), data->args);

    ctrls->vbox = gtk_vbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_input, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_results, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_save, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->vbox);
    gtk_widget_show_all(ctrls->window);
    ctrls->has_window = TRUE;
}

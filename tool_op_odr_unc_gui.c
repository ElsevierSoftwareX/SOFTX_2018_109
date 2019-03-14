#ifndef NOFORTRAN

#include "tool_op_odr_unc_gui.h"
#include "tool_op_odr_unc.h"

#include "controls.h"
#include "datatypes.h"
#include "gui-utils.h"
#include "niget-common.h"

#include "niget-gtk.h" /* controls definitions */

#include <gtk/gtk.h>

#define MC_WHEN_WARN 5e2

static GtkWidget *op_odr_uncertainty_standard(Data *data, OPODRUnc2Controls *ctrls);
static GtkWidget *op_odr_uncertainty_mc(Data *data, OPODRUnc2Controls *ctrls);
static GtkWidget *op_odr_uncertainty_cp(Data *data, OPODRUnc2Controls *ctrls);

static void unc_depth_changed(Data *data);
static void unc_load_changed(Data *data);
static void unc_nu_changed(Data *data);
static void unc_nui_changed(Data *data);
static void unc_Ei_changed(Data *data);

static void op_odr_uncertainty_montecarlo(Data *data);
static void op_odr_uncertainty_montecarlo2(Data *data);

static gboolean mc_iterations_changed(Data *data);

static void op_odr_mc_hist_hc(Data *data);
static void op_odr_mc_hist_Aphc(Data *data);
static void op_odr_mc_hist_Hit(Data *data);
static void op_odr_mc_hist_Er(Data *data);
static void op_odr_mc_hist_Eit(Data *data);
static void op_odr_mc_hist_S(Data *data);
static void op_odr_mc_hist_m(Data *data);
static void op_odr_mc_hist_hp(Data *data);
static void op_odr_mc_hist_eps(Data *data);

static void op_odr_uncertainty_save_data(Args *args);
static void op_odr_mc_save_data(Args *args);

void op_odr_uncertainty(Data *data)
{
    GtkWidget *label, *spin;
    gint i, j;
    OPODRUncControls *ctrls;
    OPODRUncdata *unc;
    gchar str[500];

    if ((data->args->opodrunc.Ec != NULL) && (data->args->opodrunc.Hc != NULL) && (data->args->opodrunc.Ac != NULL) && (data->args->opodrunc.hc != NULL)) {
        gtk_window_present(GTK_WINDOW(data->controls->opodrunccontrols->window));
        return;
    }

    unc = &(data->args->opodrunc);
    ctrls = data->controls->opodrunccontrols;

    /*copy niget-wide Instdata to local Instdata (only in this function) */
    unc->instdata = data->args->instdata;

    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Uncertainty in Oliver Pharr analysis (ODR fit)");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(op_odr_unc_close_window), data);
    gtk_window_set_transient_for(GTK_WINDOW(ctrls->window), GTK_WINDOW(data->controls->maincontrols->window_main));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(ctrls->window), TRUE);

    /* Uncertainties due to uncertainties in input values */

    ctrls->frame_unc_input = gtk_frame_new(NULL);
    //	gtk_frame_set_label_widget (GTK_FRAME (ctrls->frame_unc_input),
    //				    GTK_WIDGET (gwy_label_new_header("Uncertainties due to uncertainties in input values")));
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_unc_input),
                               GTK_WIDGET(gwy_label_new_header("Uncertainties in input values")));

    ctrls->table_unc_input = gtk_table_new(6, 11, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_unc_input), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_unc_input), ctrls->table_unc_input);

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(m)"), 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(S)/(mN/nm)"), 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(h<sub>c</sub>)/nm"), 5, 6, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(A<sub>p</sub>)/nm<sup>2</sup>"), 6, 7, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(H<sub>it</sub>)/MPa"), 7, 8, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(E<sub>r</sub>)/GPa"), 8, 9, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(E<sub>it</sub>)/GPa"), 9, 10, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    i++;

    /* u(h) */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(h)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->uh = gtk_adjustment_new(unc->instdata.sigma_h, 0.0, 1e3, 0.01, 0.1, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uh), 0.001, 3);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->uh, "value-changed", G_CALLBACK(unc_depth_changed), data);

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /* u(F) */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->uF = gtk_adjustment_new(unc->instdata.sigma_F, 0.0, 1e3, 0.01, 0.1, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uF), 0.000001, 6);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->uF, "value-changed", G_CALLBACK(unc_load_changed), data);

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /* u(nu) */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(ν)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->unu = gtk_adjustment_new(unc->instdata.unu, 0, 0.5, 0.01, 0.1, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->unu), 0.01, 2);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->unu, "value-changed", G_CALLBACK(unc_nu_changed), data);

    ctrls->uEitunu = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEitunu), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEitunu, 9, 10, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    i++;

    /* u(nu_i) */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(ν<sub>i</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->unui = gtk_adjustment_new(unc->instdata.unui, 0, 0.5, 0.01, 0.1, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->unui), 0.01, 2);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->unui, "value-changed", G_CALLBACK(unc_nui_changed), data);

    ctrls->uEitunui = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEitunui), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEitunui, 9, 10, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    i++;

    /* u(E_i) */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(E<sub>i</sub>)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->uEi = gtk_adjustment_new(unc->instdata.uEi * 1e-9, 0, 1e5, 1, 10, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uEi), 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->uEi, "value-changed", G_CALLBACK(unc_Ei_changed), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->uEituEi = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEituEi), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEituEi, 9, 10, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    i++;

    /* total */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("Total"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->um = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->um), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->um, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uS = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uS), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uS, 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uhc = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uhc), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uhc, 5, 6, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uA = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uA), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uA, 6, 7, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uHit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uHit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uHit, 7, 8, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uEr = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEr), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEr, 8, 9, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uEit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEit, 9, 10, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);


    /* Uncertainties due to contact point */

    ctrls->frame_unc_contact = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_unc_contact),
                               GTK_WIDGET(gwy_label_new_header("Uncertainties due to choice of contact point")));

    ctrls->table_unc_contact = gtk_table_new(2 * CONTACTMAX + 1 + 1, 5, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_unc_contact), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_unc_contact), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_unc_contact), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_unc_contact), ctrls->table_unc_contact);


    // Ec is Er, Hc is Hit, Ac is Aphc
    unc->Ec = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
    unc->Hc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
    unc->Ac = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
    unc->hc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));


    /* table is filled in propagate_uncertainties, because it is updated at any change of uh, uF */
    op_odr_propagate_uncertainties(data->controls->opodrunccontrols, &(data->args->opodrdata), &(data->args->opodrunc),
                                   &(data->args->fddata), &(data->args->area));

    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(op_odr_uncertainty_save_data), data->args);

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

    ctrls->spin_Nmc = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->Nmc), 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(ctrls->spin_Nmc), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_mc), ctrls->spin_Nmc, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->Nmc, "value-changed", G_CALLBACK(mc_iterations_changed), data);
    i++;

    ctrls->mc_warn = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(ctrls->table_mc), ctrls->mc_warn, 0, 3, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    if (unc->Nmc > MC_WHEN_WARN) {
        label_set_gint_format(ctrls->mc_warn, unc->Nmc,
                              "Warning: Monte Carlo calculation for %d iterations may take very long!");
    }

    // no Monte Carlo window without the uncertainty dialog
    data->controls->opodrmccontrols->has_window = FALSE;
    ctrls->button_mc = gtk_button_new_with_mnemonic("_Monte Carlo");
    gtk_button_set_image(GTK_BUTTON(ctrls->button_mc), gtk_image_new_from_stock(GWY_STOCK_POLYNOM, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_mc, "clicked", G_CALLBACK(op_odr_uncertainty_montecarlo), data);

    ctrls->vbox = gtk_vbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_unc_input, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_unc_contact, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_save, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_mc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_mc, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->vbox);
    gtk_widget_show_all(ctrls->window);
}

#define MAXBAD 5 //maximum percentage of nonsense results before warning is triggered
static GtkWidget *op_odr_uncertainty_standard(Data *data, OPODRUnc2Controls *ctrls)
{
    gint i;
    GtkWidget *spin;
    OPODRUncdata *unc;
    Instdata *inst;

    ctrls->vbox_prop = gtk_vbox_new(FALSE, CTRLS_SPACING);

    inst = &(data->args->instdata);
    unc = &(data->args->opodrunc);

    ctrls->table_unc_input = gtk_table_new(7, 11, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_unc_input), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_ROW_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox_prop), ctrls->table_unc_input, FALSE, FALSE, CTRLS_SPACING);

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(m)"), 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(S)/(mN/nm)"), 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(h<sub>c</sub>)/nm"), 5, 6, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(A<sub>p</sub>)/nm<sup>2</sup>"), 6, 7, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(H<sub>it</sub>)/MPa"), 7, 8, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(E<sub>r</sub>)/GPa"), 8, 9, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(E<sub>it</sub>)/GPa"), 9, 10, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    i++;

    /* u(h) */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(h)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->uh = gtk_adjustment_new(unc->instdata.sigma_h, 0.0, 1e3, 0.01, 0.1, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uh), 0.001, 3);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), spin, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->uh, "value-changed", G_CALLBACK(unc_depth_changed), data);

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /* u(F) */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->uF = gtk_adjustment_new(unc->instdata.sigma_F, 0.0, 1e3, 0.01, 0.1, 0.);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->uF), 0.000001, 6);
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
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEitunu, 9, 10, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

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
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEitunui, 9, 10, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

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
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEituEi, 9, 10, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    i++;

    /* total */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("Total"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->um = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->um), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->um, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uS = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uS), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uS, 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uhc = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uhc), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uhc, 5, 6, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uA = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uA), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uA, 6, 7, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uHit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uHit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uHit, 7, 8, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uEr = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEr), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEr, 8, 9, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uEit = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uEit), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uEit, 9, 10, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    op_odr_propagate_uncertainties(data->controls->opodrunccontrols, &(data->args->opodrdata), &(data->args->opodrunc),
                                   &(data->args->fddata), &(data->args->area));

    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(op_odr_uncertainty_save_data), data->args);
    return NULL; /* shut-up warning */
}

/* podobne jako hlavni gui - separatni funkce poskytujici jednotlive taby */
void op_odr_uncertainty2(Data *data)
{
    OPODRUnc2Controls *ctrls;

    if ((data->args->opodrunc.Ec != NULL) && (data->args->opodrunc.Hc != NULL) && (data->args->opodrunc.Ac != NULL) && (data->args->opodrunc.hc != NULL)) {
        gtk_window_present(GTK_WINDOW(data->controls->opodrunccontrols->window));
        return;
    }

    //ctrls =	data->controls->opodrunccontrols;
    ctrls = g_malloc(sizeof(OPODRUnc2Controls));

    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Uncertainty in Oliver Pharr analysis (ODR fit)");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(op_odr_unc_close_window), data);
    gtk_window_set_transient_for(GTK_WINDOW(ctrls->window), GTK_WINDOW(data->controls->maincontrols->window_main));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(ctrls->window), TRUE);

    ctrls->notebook = gtk_notebook_new();
    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->notebook);

    ctrls->vbox_prop = op_odr_uncertainty_standard(data, ctrls);
    gtk_notebook_append_page(GTK_NOTEBOOK(ctrls->notebook), ctrls->vbox_prop, gtk_label_new("Instrument (unc. propagation)"));

    ctrls->vbox_mc = op_odr_uncertainty_mc(data, ctrls);
    gtk_notebook_append_page(GTK_NOTEBOOK(ctrls->notebook), ctrls->vbox_mc, gtk_label_new("Instrument (Monte Carlo)"));

    ctrls->vbox_cp = op_odr_uncertainty_cp(data, ctrls);
    gtk_notebook_append_page(GTK_NOTEBOOK(ctrls->notebook), ctrls->vbox_cp, gtk_label_new("Contact point choice"));

    gtk_widget_show_all(ctrls->window);
}

#define NPARAM 9  //hc, Aphc, Hit, Er, Eit, S, m, hp, eps
static GtkWidget *op_odr_uncertainty_mc(Data *data, OPODRUnc2Controls *ctrls)
{
    gint i, j;
    GtkWidget *label;
    OPODRUncdata *unc;
    OPODRMCdata *opodrmc;
    GdkColor color;

    gdk_color_parse("red", &color);

    unc = &(data->args->opodrunc);
    opodrmc = &(data->args->opodrmc);

    ctrls->vbox_mc = gtk_vbox_new(FALSE, CTRLS_SPACING);

    ctrls->frame_input = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_input),
                               GTK_WIDGET(gwy_label_new_header("Monte Carlo calculation of uncertainties")));
    gtk_box_pack_start(GTK_BOX(ctrls->vbox_mc), ctrls->frame_input, FALSE, FALSE, CTRLS_SPACING);

    ctrls->table_input = gtk_table_new(4, 4, FALSE); //rows, cols
    gtk_container_add(GTK_CONTAINER(ctrls->frame_input), ctrls->table_input);
    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_input), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_input), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_input), CTRLS_TABLE_ROW_SPACING);

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("Number of iterations"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->N_mc = gtk_adjustment_new(unc->Nmc, 1, 1e9, 10, 100, 0);
    ctrls->spin_N_mc = gtk_spin_button_new(GTK_ADJUSTMENT(ctrls->N_mc), 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(ctrls->spin_N_mc), TRUE);
    gtk_table_attach(GTK_TABLE(ctrls->table_input), ctrls->spin_N_mc, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(ctrls->N_mc, "value-changed", G_CALLBACK(mc_iterations_changed), data);

    i++;

    ctrls->mc_warn = gtk_label_new(NULL);
    gtk_table_attach(GTK_TABLE(ctrls->table_input), ctrls->mc_warn, 0, 3, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    if (unc->Nmc > MC_WHEN_WARN) {
        label_set_gint_format(ctrls->mc_warn, unc->Nmc,
                              "Warning: Monte Carlo calculation for %d iterations may take very long!");
    }

    // no Monte Carlo window without the uncertainty dialog
    data->controls->opodrmccontrols->has_window = FALSE;
    ctrls->button_run = gtk_button_new_with_mnemonic("_Run");
    gtk_button_set_image(GTK_BUTTON(ctrls->button_run), gtk_image_new_from_stock(GWY_STOCK_POLYNOM, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_run, "clicked", G_CALLBACK(op_odr_uncertainty_montecarlo2), data);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox_mc), ctrls->button_run, FALSE, FALSE, CTRLS_SPACING);


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
    label_set_gdouble_format(label, data->args->opodrunc.instdata.sigma_h, "%g");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /*          u(F)             */
    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrunc.instdata.sigma_F, "%g");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /*          Number of iteration */
    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("Number of iterations "), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gint_format(label, opodrmc->Nmc, "%d");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    // Results
    ctrls->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_results),
                               GTK_WIDGET(gwy_label_new_header("Results of Monte Carlo simulation")));

    ctrls->table_results = gtk_table_new(NPARAM + 1, 4, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_results), ctrls->table_results);

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("mean"), 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("stdev"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    ctrls->button_hist = (GtkWidget **)g_malloc(NPARAM * sizeof(GtkWidget *));
    ctrls->hist_window = (GtkWidget **)g_malloc(NPARAM * sizeof(GtkWidget *));

    for (j = 0; j < NPARAM; j++) {
        ctrls->hist_window[j] = NULL;
    }

    /*                       h_c           */
    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("h<sub>c</sub>/nm"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_hc), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       A_p(h_c)           */
    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("A<sub>p</sub>/nm<sup>2</sup>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], AREA);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], AREA);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_Aphc), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       H_IT           */
    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("H<sub>it</sub>/MPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], HARD);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], HARD);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_Hit), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       E_r           */
    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("E<sub>r</sub>/GPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_Er), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       E_IT          */
    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("E<sub>it</sub>/GPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_Eit), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       S          */
    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("S/(mN/nm)"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], SLOPE);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], SLOPE);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_S), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       m          */
    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("m"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], NUMBER);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], NUMBER);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_m), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       hp          */
    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("hp"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_hp), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       eps          */
    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("eps"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], NUMBER);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], NUMBER);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_eps), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    if (100.*data->args->opodrmc.skipmc / opodrmc->Nmc > MAXBAD) {
        label = gtk_label_new(NULL);
        label_set_gdouble_format(label, 100. * data->args->opodrmc.skipmc / opodrmc->Nmc,
                                 "Warning: Approx %.1f %% of the results had to be removed!");
        gtk_widget_modify_fg(label, GTK_STATE_NORMAL, &color);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 0, 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    }

    // save data
    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(op_odr_mc_save_data), data->args);

    ctrls->vbox_mc = gtk_vbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox_mc), ctrls->frame_input, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox_mc), ctrls->frame_results, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox_mc), ctrls->button_save, FALSE, FALSE, 0);

    return NULL; /* shut-up warning */
}

static GtkWidget *op_odr_uncertainty_cp(Data *data, OPODRUnc2Controls *ctrls)
{
    gint i, j;
    gchar str[256];
    GtkWidget *label;
    OPODRUncdata *unc;

    gint shift_range;

    unc = &(data->args->opodrunc);

    ctrls->vbox_cp = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* ctrls->frame_unc_contact = gtk_frame_new(NULL); */
    /* gtk_frame_set_label_widget (GTK_FRAME (ctrls->frame_unc_contact), */
    /* 			    GTK_WIDGET (gwy_label_new_header("Uncertainties due to choice of contact point"))); */

    shift_range = CONTACTMAX; /* protoze se pouiziva i ve volanych funkcich; eliminovat ve prospech sdilene promenne */

    ctrls->table_unc_contact = gtk_table_new(2 * shift_range + 1 + 1, 4, FALSE);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox_cp), ctrls->table_unc_contact, FALSE, FALSE, CTRLS_SPACING);
    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_unc_contact), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_unc_contact), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_unc_contact), CTRLS_TABLE_ROW_SPACING);

    /* table headers */

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("Contact point shift"), 0, 1, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("E<sub>r</sub>/GPa"), 1, 2, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("H<sub>it</sub>/GPa"), 2, 3, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("A<sub>p</sub>(h<sub>c</sub>)/nm<sup>2</sup>"), 3, 4, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    // gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
    // gwy_label_new_header("weight"), 2, 3, i,i+1,  GTK_EXPAND | GTK_FILL, GTK_FILL, 0, 0);

    i++;

    // Ec is Er, Hc is Hit, Ac is Aphc
    unc->Ec = (gdouble *)g_malloc((2 * shift_range + 1) * sizeof(gdouble));
    unc->Hc = (gdouble *)g_malloc((2 * shift_range + 1) * sizeof(gdouble));
    unc->Ac = (gdouble *)g_malloc((2 * shift_range + 1) * sizeof(gdouble));
    unc->hc = (gdouble *)g_malloc((2 * shift_range + 1) * sizeof(gdouble));

    /* values */

    for (j = -shift_range; j <= shift_range; j++) {

        /* vertical alignments: j: 1.0 (bottom), others: 0.0 (top); looks quite nice, still probably incorrect */

        label = gtk_label_new(NULL);
        label_set_gint_format(label, j, "%d\n");
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 0, 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        op_odr_fit_shift_contact(&(data->args->opodrdata), &(data->args->fddata), &(data->args->instdata),
                                 &(data->args->area), unc->Ec, unc->Hc, unc->Ac, unc->hc, j);

        //fill in Er
        label = gtk_label_new(NULL);

        if (unc->Ec[j + shift_range] > 0) {
            g_snprintf(str, sizeof(str), MODUL, unc->Ec[j + shift_range]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        /*
          label = gtk_label_new(NULL);
          g_snprintf(str, sizeof(str), "1/%d", 2*shift_range + 1);
          gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
          gtk_table_attach (GTK_TABLE(ctrls->table_unc_contact), label, 2, 3, i, i+1, GTK_FILL, GTK_FILL, 0, 0);
        */

        //fill in Hit
        label = gtk_label_new(NULL);

        if (unc->Hc[j + shift_range] > 0) {
            g_snprintf(str, sizeof(str), HARD, unc->Hc[j + shift_range]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);


        //fill in Aphc
        label = gtk_label_new(NULL);

        if (unc->Ac[j + shift_range] > 0) {
            g_snprintf(str, sizeof(str), AREA, unc->Ac[j + shift_range]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        //fill in hc
        label = gtk_label_new(NULL);

        if (unc->hc[j + shift_range] > 0) {
            g_snprintf(str, sizeof(str), DEPTH, unc->hc[j + shift_range]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 4, 5, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);


        i++;
    }

    return NULL; /* shut-up warning */
}

static void unc_depth_changed(Data *data)
{
    data->args->opodrunc.instdata.sigma_h = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opodrunccontrols->uh));
    op_odr_propagate_uncertainties(data->controls->opodrunccontrols, &(data->args->opodrdata), &(data->args->opodrunc),
                                   &(data->args->fddata), &(data->args->area));
}

static void unc_load_changed(Data *data)
{
    data->args->opodrunc.instdata.sigma_F = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opodrunccontrols->uF));
    op_odr_propagate_uncertainties(data->controls->opodrunccontrols, &(data->args->opodrdata), &(data->args->opodrunc),
                                   &(data->args->fddata), &(data->args->area));
}

static void unc_nu_changed(Data *data)
{
    data->args->opodrunc.instdata.unu = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opodrunccontrols->unu));
    op_odr_propagate_uncertainties(data->controls->opodrunccontrols, &(data->args->opodrdata), &(data->args->opodrunc),
                                   &(data->args->fddata), &(data->args->area));
}

static void unc_nui_changed(Data *data)
{
    data->args->opodrunc.instdata.unui = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opodrunccontrols->unui));
    op_odr_propagate_uncertainties(data->controls->opodrunccontrols, &(data->args->opodrdata), &(data->args->opodrunc),
                                   &(data->args->fddata), &(data->args->area));
}

static void unc_Ei_changed(Data *data)
{
    data->args->opodrunc.instdata.uEi = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opodrunccontrols->uEi)) * 1e9;
    op_odr_propagate_uncertainties(data->controls->opodrunccontrols, &(data->args->opodrdata), &(data->args->opodrunc),
                                   &(data->args->fddata), &(data->args->area));
}

static void op_odr_uncertainty_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (op_odr_unc_has_results(&(args->opodrunc))) {
        buffer = op_odr_uncertainty_export_data(&(args->opodrdata), &(args->opodrunc), &(args->fddata),
                                                &(args->area), EXPORT_FORMAT_PLAIN);
        buffer_to_file_dialog(buffer);
        g_free(buffer);
    }
}

void op_odr_unc_close_window(Data *data)
{
    g_free(data->args->opodrunc.Ec);
    g_free(data->args->opodrunc.Hc);
    g_free(data->args->opodrunc.Ac);
    g_free(data->args->opodrunc.hc);
    data->args->opodrunc.Ec = NULL;
    data->args->opodrunc.Hc = NULL;
    data->args->opodrunc.Ac = NULL;
    data->args->opodrunc.hc = NULL;
    op_odr_mc_close_window(data);
    gtk_widget_destroy(data->controls->opodrunccontrols->window);
}

static void op_odr_uncertainty_montecarlo(Data *data)
{
    gint i, j;
    GtkWidget *label;
    GdkColor color;

    OPODRMCControls *ctrls;
    OPODRMCdata *opodrmc;

    gdk_color_parse("red", &color);

    //if there was a previous window, close it and any histograms, free data
    op_odr_mc_close_window(data);

    opodrmc = &(data->args->opodrmc);

    // initialize for number of iterations which was chosen in the Uncertainties dialog
    init_OPODRMCdata(opodrmc, data->args->opodrunc.Nmc);

    // run calculation
    op_odr_uncertainty_montecarlo_run_calc(&(data->args->opodrdata), &(data->args->opodrunc), &(data->args->opodrmc), &(data->args->fddata), &(data->args->area));

    //create window and display results
    ctrls = data->controls->opodrmccontrols;

    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Monte Carlo calculation of uncertainties in Oliver Pharr analysis (ODR fit)");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(op_odr_mc_close_window), data);
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
    label_set_gdouble_format(label, data->args->opodrunc.instdata.sigma_h, "%g");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /*          u(F)             */

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrunc.instdata.sigma_F, "%g");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /*          Number of iteration */

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("Number of iterations "), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gint_format(label, opodrmc->Nmc, "%d");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    // Results
    ctrls->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_results),
                               GTK_WIDGET(gwy_label_new_header("Results of Monte Carlo simulation")));

    ctrls->table_results = gtk_table_new(NPARAM + 1, 4, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_results), ctrls->table_results);

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("mean"), 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("stdev"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    ctrls->button_hist = (GtkWidget **)g_malloc(NPARAM * sizeof(GtkWidget *));
    ctrls->hist_window = (GtkWidget **)g_malloc(NPARAM * sizeof(GtkWidget *));

    for (j = 0; j < NPARAM; j++) {
        ctrls->hist_window[j] = NULL;
    }

    /*                       h_c           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("h<sub>c</sub>/nm"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_hc), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       A_p(h_c)           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("A<sub>p</sub>/nm<sup>2</sup>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], AREA);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], AREA);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_Aphc), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       H_IT           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("H<sub>it</sub>/MPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], HARD);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], HARD);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_Hit), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       E_r           */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("E<sub>r</sub>/GPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_Er), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       E_IT          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("E<sub>it</sub>/GPa"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], MODUL);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_Eit), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       S          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("S/(mN/nm)"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], SLOPE);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], SLOPE);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_S), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       m          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("m"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], NUMBER);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], NUMBER);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_m), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       hp          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("hp"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_m), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       eps          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("eps"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcavg[i - 1], NUMBER);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->opodrmc.mcstd[i - 1], NUMBER);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(op_odr_mc_hist_m), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    if (100.*data->args->opodrmc.skipmc / opodrmc->Nmc > MAXBAD) {
        label = gtk_label_new(NULL);
        label_set_gdouble_format(label, 100. * data->args->opodrmc.skipmc / opodrmc->Nmc,
                                 "Warning: Approx %.1f %% of the results had to be removed!");
        gtk_widget_modify_fg(label, GTK_STATE_NORMAL, &color);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 0, 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    }

    // save data
    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(op_odr_mc_save_data), data->args);

    ctrls->vbox = gtk_vbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_input, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_results, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_save, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->vbox);
    gtk_widget_show_all(ctrls->window);
    ctrls->has_window = TRUE;
}

static void op_odr_mc_hist_hc(Data *data)
{
    show_histogram("h_c", "nm", data->args->opodrmc.mcxhist[0], data->args->opodrmc.mcyhist[0], data->args->opodrmc.mcnstat, data->controls->opodrmccontrols->window);
}

static void op_odr_mc_hist_Aphc(Data *data)
{
    show_histogram("A_p(h_c)", "nm^2", data->args->opodrmc.mcxhist[1], data->args->opodrmc.mcyhist[1], data->args->opodrmc.mcnstat, data->controls->opodrmccontrols->window);
}

static void op_odr_mc_hist_Hit(Data *data)
{
    show_histogram("H_IT", "MPa", data->args->opodrmc.mcxhist[2], data->args->opodrmc.mcyhist[2], data->args->opodrmc.mcnstat, data->controls->opodrmccontrols->window);
}

static void op_odr_mc_hist_Er(Data *data)
{
    show_histogram("E_r", "GPa", data->args->opodrmc.mcxhist[3], data->args->opodrmc.mcyhist[3], data->args->opodrmc.mcnstat, data->controls->opodrmccontrols->window);
}

static void op_odr_mc_hist_Eit(Data *data)
{
    show_histogram("E_IT", "GPa", data->args->opodrmc.mcxhist[4], data->args->opodrmc.mcyhist[4], data->args->opodrmc.mcnstat, data->controls->opodrmccontrols->window);
}

static void op_odr_mc_hist_S(Data *data)
{
    show_histogram("S", "mN/nm", data->args->opodrmc.mcxhist[5], data->args->opodrmc.mcyhist[5], data->args->opodrmc.mcnstat, data->controls->opodrmccontrols->window);
}

static void op_odr_mc_hist_m(Data *data)
{
    show_histogram("m", "", data->args->opodrmc.mcxhist[6], data->args->opodrmc.mcyhist[6], data->args->opodrmc.mcnstat, data->controls->opodrmccontrols->window);
}

static void op_odr_mc_hist_hp(Data *data)
{
    show_histogram("hp", "nm", data->args->opodrmc.mcxhist[7], data->args->opodrmc.mcyhist[7], data->args->opodrmc.mcnstat, data->controls->opodrmccontrols->window);
}

static void op_odr_mc_hist_eps(Data *data)
{
    show_histogram("eps", "", data->args->opodrmc.mcxhist[8], data->args->opodrmc.mcyhist[8], data->args->opodrmc.mcnstat, data->controls->opodrmccontrols->window);
}

static void op_odr_mc_save_data(Args *args)
{
    gchar *buffer = NULL;

    /* add check for mc results */
    // if (op_odr_mc_has_results(&(args->opodrmc))) {
    buffer = op_odr_mc_export_data(&(args->opodrdata), &(args->opodrunc), &(args->opodrmc), &(args->fddata),
                                   EXPORT_FORMAT_PLAIN);
    buffer_to_file_dialog(buffer);
    g_free(buffer);
    // }
}

/* this is not a callback, therefore no need to pass whole Data */
void op_odr_propagate_uncertainties(OPODRUncControls *ctrls, const OPODRdata *opodrdata, OPODRUncdata *unc,
                                    const FDdata *fddata, const Area *area)
{
    gchar str[30];
    gint i, j;
    GtkWidget *label;

    if (!opodrdata->has_fit) {
        return;
    }

    if (unc->Ec == NULL || unc->Hc == NULL || unc->Ac == NULL || unc->hc == NULL) {
        return;
    }

    op_odr_propagate_uncertainties_calc(opodrdata, unc, fddata, area);

    /* show */

    // uncertainties from uh and uF
    label_set_gdouble_format(ctrls->um, unc->um, NUMBER);
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

    /* contact point uncertainties */

    //destroy the labels of the old table_unc_contact
    remove_table_content(ctrls->table_unc_contact);

    //create new table
    /* table headers */

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("Contact point shift"), 0, 1, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("E<sub>r</sub>/GPa"), 1, 2, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("H<sub>it</sub>/GPa"), 2, 3, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("A<sub>p</sub>(h<sub>c</sub>)/nm<sup>2</sup>"), 3, 4, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("h<sub>c</sub>/nm"), 4, 5, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    // gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
    // gwy_label_new_header("weight"), 2, 3, i,i+1,  GTK_EXPAND | GTK_FILL, GTK_FILL, 0, 0);

    i++;

    /* values */

    for (j = -CONTACTMAX; j <= CONTACTMAX; j++)	{

        /* vertical alignments: j: 1.0 (bottom), others: 0.0 (top); looks quite nice, still probably incorrect */

        label = gtk_label_new(NULL);
        label_set_gint_format(label, j, "%d\n");
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 1.0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 0, 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        op_odr_fit_shift_contact(opodrdata, fddata, &(unc->instdata),
                                 area, unc->Ec, unc->Hc, unc->Ac, unc->hc, j);

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
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

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
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        //fill in Aphc
        label = gtk_label_new(NULL);

        if (unc->Ac[j + CONTACTMAX] > 0) {
            g_snprintf(str, sizeof(str), AREA, unc->Ac[j + CONTACTMAX]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        //fill in hc
        label = gtk_label_new(NULL);

        if (unc->hc[j + CONTACTMAX] > 0) {
            g_snprintf(str, sizeof(str), DEPTH, unc->hc[j + CONTACTMAX]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 4, 5, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        i++;
    }

    gtk_widget_show_all(ctrls->table_unc_contact);
}

static gboolean mc_iterations_changed(Data *data)
{
    data->args->opodrunc.Nmc = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->opodrunccontrols->Nmc));

    if (data->args->opodrunc.Nmc > MC_WHEN_WARN) {
        label_set_gint_format(data->controls->opodrunccontrols->mc_warn, data->args->opodrunc.Nmc,
                              "Warning: Monte Carlo calculation for %d iterations may take very long!");
    }
    else {
        label_clear(data->controls->opodrunccontrols->mc_warn);
    }

    return FALSE;
}

void op_odr_uncertainty_montecarlo2(Data *data)
{
    OPODRMCControls *ctrls;
    OPODRMCdata *opodrmc;

    //if there was a previous window, close it and any histograms, free data
    op_odr_mc_close_window(data);

    opodrmc = &(data->args->opodrmc);

    // initialize for number of iterations which was chosen in the Uncertainties dialog
    init_OPODRMCdata(opodrmc, data->args->opodrunc.Nmc);

    // run calculation
    op_odr_uncertainty_montecarlo_run_calc(&(data->args->opodrdata), &(data->args->opodrunc), &(data->args->opodrmc), &(data->args->fddata), &(data->args->area));

    //create window and display results
    ctrls = data->controls->opodrmccontrols;

    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->vbox);
    gtk_widget_show_all(ctrls->window);
    ctrls->has_window = TRUE;
}

void op_odr_mc_close_window(Data *data)
{
    gint j;

    if (data->controls->opodrmccontrols->has_window) {
        gtk_widget_destroy(data->controls->opodrmccontrols->window);

        destroy_OPODRMCdata(&(data->args->opodrmc));
        data->controls->opodrmccontrols->has_window = FALSE;

        for (j = 0; j < NPARAM; j++) {
            if (GTK_IS_WIDGET(data->controls->opodrmccontrols->hist_window[j])) {
                gtk_widget_destroy(data->controls->opodrmccontrols->hist_window[j]);
            }
        }

        g_free(data->controls->opodrmccontrols->hist_window);
        g_free(data->controls->opodrmccontrols->button_hist);
    }
}

#endif

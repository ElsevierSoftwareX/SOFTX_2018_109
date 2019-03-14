#ifndef NOFORTRAN

#include "tool_stiffness_unc_gui.h"
#include "tool_stiffness_unc.h"
#include "tool_stiffness.h"

#include "controls.h"
#include "datatypes.h"
#include "fddata.h"
#include "fit-utils.h"
#include "file-utils.h"
#include "gui-utils.h"
#include "niget-common.h"

#include "niget-gtk.h" /* controls definitions */

#include <math.h>
#include <string.h>

#include <gtk/gtk.h>

#include <app/file.h>

#define MC_WHEN_WARN 5e2

static gboolean mc_iterations_changed(Data *data);
static void unc_depth_changed(Data *data);
static void unc_load_changed(Data *data);

static void stiffness_mc_hist_kl(Data *data);
static void stiffness_mc_hist_ql(Data *data);
static void stiffness_mc_hist_ku(Data *data);
static void stiffness_mc_hist_qu(Data *data);
static void stiffness_uncertainty_save_data(Args *args);
static void stiffness_mc_save_data(Args *args);


void stiffness_unc_close_window(Data *data)
{
    g_free(data->args->stiffnessunc.klc);
    g_free(data->args->stiffnessunc.qlc);
    g_free(data->args->stiffnessunc.kuc);
    g_free(data->args->stiffnessunc.quc);
    data->args->stiffnessunc.klc = NULL;
    data->args->stiffnessunc.qlc = NULL;
    data->args->stiffnessunc.kuc = NULL;
    data->args->stiffnessunc.quc = NULL;
    stiffness_mc_close_window(data);
    gtk_widget_destroy(data->controls->stiffnessunccontrols->window);
}

void stiffness_uncertainty(Data *data)
{
    GtkWidget *label, *spin;
    gint i, j;
    StiffnessUncControls *ctrls;
    StiffnessUncdata *unc;
    gchar str[500];

    if ((data->args->stiffnessunc.klc != NULL) && (data->args->stiffnessunc.qlc != NULL)
            && (data->args->stiffnessunc.kuc != NULL) && (data->args->stiffnessunc.quc != NULL)) {
        gtk_window_present(GTK_WINDOW(data->controls->stiffnessunccontrols->window));
        return;
    }

    unc = &(data->args->stiffnessunc);
    ctrls = data->controls->stiffnessunccontrols;

    /*copy niget-wide Instdata to local Instdata (only in this function) */
    unc->instdata = data->args->instdata;


    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Uncertainty in stiffness analysis (ODR fit)");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(stiffness_unc_close_window), data);
    gtk_window_set_transient_for(GTK_WINDOW(ctrls->window), GTK_WINDOW(data->controls->maincontrols->window_main));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(ctrls->window), TRUE);

    /* Uncertainties due to uncertainties in input values */

    ctrls->frame_unc_input = gtk_frame_new(NULL);
    //	gtk_frame_set_label_widget (GTK_FRAME (ctrls->frame_unc_input),
    //				    GTK_WIDGET (gwy_label_new_header("Uncertainties due to uncertainties in input values")));
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_unc_input),
                               GTK_WIDGET(gwy_label_new_header("Uncertainties in input values")));

    ctrls->table_unc_input = gtk_table_new(6, 13, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_unc_input), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_unc_input), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_unc_input), ctrls->table_unc_input);

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(k<sub>load</sub>)/(mN/nm)"), 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(q<sub>load</sub>)/nm"), 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(k<sub>unload</sub>)/(mN/nm)"), 5, 6, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input),
                     gwy_label_new_header("u(q<sub>unload</sub>)/nm"), 6, 7, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

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

    /* total */

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), gwy_label_new_header("Total"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    ctrls->ukl = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->ukl), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->ukl, 3, 4, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uql = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uql), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uql, 4, 5, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uku = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uku), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uku, 5, 6, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    ctrls->uqu = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(ctrls->uqu), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_input), ctrls->uqu, 6, 7, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);



    /* Uncertainties due to contact point */

    ctrls->frame_unc_contact = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_unc_contact),
                               GTK_WIDGET(gwy_label_new_header("Uncertainties due to choice of contact point")));

    ctrls->table_unc_contact = gtk_table_new(2 * CONTACTMAX + 1 + 1, 5, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_unc_contact), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_unc_contact), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_unc_contact), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_unc_contact), ctrls->table_unc_contact);

    unc->klc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
    unc->qlc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
    unc->kuc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
    unc->quc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));


    /* table is filled in propagate_uncertainties, because it is updated at any change of uh, uF */
    stiffness_propagate_uncertainties(data->controls->stiffnessunccontrols, &(data->args->stiffnessdata), &(data->args->stiffnessunc), &(data->args->fddata));

    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(stiffness_uncertainty_save_data), data->args);

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
        label_set_gint_format(ctrls->mc_warn, unc->Nmc, "Warning: Monte Carlo calculation for %d iterations may take very long!");
    }

    // no Monte Carlo window without the uncertainty dialog
    data->controls->stiffnessmccontrols->has_window = FALSE;
    ctrls->button_mc = gtk_button_new_with_mnemonic("_Monte Carlo");
    gtk_button_set_image(GTK_BUTTON(ctrls->button_mc), gtk_image_new_from_stock(GWY_STOCK_POLYNOM, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_mc, "clicked", G_CALLBACK(stiffness_uncertainty_montecarlo), data);

    ctrls->vbox = gtk_vbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_unc_input, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_unc_contact, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_save, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_mc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_mc, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->vbox);
    gtk_widget_show_all(ctrls->window);
}

static void stiffness_uncertainty_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (stiffness_unc_has_results(&(args->stiffnessunc))) {
        buffer = stiffness_uncertainty_export_data(&(args->stiffnessdata), &(args->stiffnessunc), &(args->fddata), EXPORT_FORMAT_PLAIN);
        buffer_to_file_dialog(buffer);
        g_free(buffer);
    }
}

static gboolean mc_iterations_changed(Data *data)
{
    data->args->stiffnessunc.Nmc = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->stiffnessunccontrols->Nmc));

    if (data->args->stiffnessunc.Nmc > MC_WHEN_WARN) {
        label_set_gint_format(data->controls->stiffnessunccontrols->mc_warn, data->args->stiffnessunc.Nmc,
                              "Warning: Monte Carlo calculation for %d iterations may take very long!");
    }
    else {
        label_clear(data->controls->stiffnessunccontrols->mc_warn);
    }

    return FALSE;
}

static void unc_depth_changed(Data *data)
{
    data->args->stiffnessunc.instdata.sigma_h = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->stiffnessunccontrols->uh));
    stiffness_propagate_uncertainties(data->controls->stiffnessunccontrols, &(data->args->stiffnessdata), &(data->args->stiffnessunc),
                                      &(data->args->fddata));
}

static void unc_load_changed(Data *data)
{
    data->args->stiffnessunc.instdata.sigma_F = gtk_adjustment_get_value(GTK_ADJUSTMENT(data->controls->stiffnessunccontrols->uF));
    stiffness_propagate_uncertainties(data->controls->stiffnessunccontrols, &(data->args->stiffnessdata), &(data->args->stiffnessunc),
                                      &(data->args->fddata));
}


void stiffness_uncertainty_montecarlo(Data *data)
{
    gint i, j;
    GtkWidget *label;
    GdkColor color;

    StiffnessMCControls *ctrls;
    StiffnessMCdata *stiffnessmc;
    Instdata *instdata;

    instdata = &(data->args->stiffnessunc.instdata);
    gdk_color_parse("red", &color);

    //if there was a previous window, close it and any histograms, free data
    stiffness_mc_close_window(data);

    stiffnessmc = &(data->args->stiffnessmc);

    // initialize for number of iterations which was chosen in the Uncertainties dialog
    init_StiffnessMCdata(stiffnessmc, data->args->stiffnessunc.Nmc);

    // run calculation
    stiffness_uncertainty_montecarlo_run_calc(&(data->args->stiffnessdata), &(data->args->stiffnessunc), &(data->args->stiffnessmc), &(data->args->fddata));

    //create window and display results
    ctrls = data->controls->stiffnessmccontrols;

    ctrls->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctrls->window), "Monte Carlo calculation of uncertainties in two stiffness analysis (ODR fit)");
    gtk_window_set_default_size(GTK_WINDOW(ctrls->window), 400, 300);
    g_signal_connect_swapped(ctrls->window, "delete_event", G_CALLBACK(stiffness_mc_close_window), data);
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
    label_set_gdouble_format(label, instdata->sigma_h, "%g");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /*          u(F)             */

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("u(F)"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, instdata->sigma_F, "%g");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    gtk_table_attach(GTK_TABLE(ctrls->table_input), label_new_with_markup_left("normal distribution"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    /*          Number of iteration */

    gtk_table_attach(GTK_TABLE(ctrls->table_input), gwy_label_new_header("Number of iterations "), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gint_format(label, stiffnessmc->Nmc, "%d");
    gtk_table_attach(GTK_TABLE(ctrls->table_input), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    // Results
    ctrls->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(ctrls->frame_results),
                               GTK_WIDGET(gwy_label_new_header("Results of Monte Carlo simulation")));

    ctrls->table_results = gtk_table_new(NPARAM_STIFFNESS + 1, 4, FALSE);

    gtk_container_set_border_width(GTK_CONTAINER(ctrls->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(ctrls->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(ctrls->frame_results), ctrls->table_results);

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("mean"), 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), gwy_label_new_header("stdev"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    i++;

    ctrls->button_hist = (GtkWidget **)g_malloc(NPARAM_STIFFNESS * sizeof(GtkWidget *));
    ctrls->hist_window = (GtkWidget **)g_malloc(NPARAM_STIFFNESS * sizeof(GtkWidget *));

    for (j = 0; j < NPARAM_STIFFNESS; j++) {
        ctrls->hist_window[j] = NULL;
    }

    /*                       kload          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("kl/(mN/nm)"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->stiffnessmc.mcavg[i - 1], SLOPE);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->stiffnessmc.mcstd[i - 1], SLOPE);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(stiffness_mc_hist_kl), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       qload          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("ql/nm"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->stiffnessmc.mcavg[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->stiffnessmc.mcstd[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(stiffness_mc_hist_ql), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;

    /*                       kunload          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("ku/(mN/nm)"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->stiffnessmc.mcavg[i - 1], SLOPE);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->stiffnessmc.mcstd[i - 1], SLOPE);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(stiffness_mc_hist_ku), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;


    /*                       qunload          */

    gtk_table_attach(GTK_TABLE(ctrls->table_results),
                     gwy_label_new_header("qu/nm"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0); //TODO units

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->stiffnessmc.mcavg[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 1, 2, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    label = gtk_label_new(NULL);
    label_set_gdouble_format(label, data->args->stiffnessmc.mcstd[i - 1], DEPTH);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

    ctrls->button_hist[i - 1] = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(ctrls->button_hist[i - 1]), gtk_image_new_from_stock(GWY_STOCK_GRAPH, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(ctrls->button_hist[i - 1], "clicked", G_CALLBACK(stiffness_mc_hist_qu), data);
    gtk_table_attach(GTK_TABLE(ctrls->table_results), ctrls->button_hist[i - 1], 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    i++;


    if (100.*data->args->stiffnessmc.skipmc / stiffnessmc->Nmc > MAXBAD) {
        label = gtk_label_new(NULL);
        label_set_gdouble_format(label, 100. * data->args->stiffnessmc.skipmc / stiffnessmc->Nmc,
                                 "Warning: Approx %.1f %% of the results had to be removed!");
        gtk_widget_modify_fg(label, GTK_STATE_NORMAL, &color);
        gtk_table_attach(GTK_TABLE(ctrls->table_results), label, 0, 1, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);
    }

    // save data
    ctrls->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(ctrls->button_save, "clicked", G_CALLBACK(stiffness_mc_save_data), data->args);

    ctrls->vbox = gtk_vbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_input, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->frame_results, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(ctrls->vbox), ctrls->button_save, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(ctrls->window), ctrls->vbox);
    gtk_widget_show_all(ctrls->window);
    ctrls->has_window = TRUE;
}

static void stiffness_mc_save_data(Args *args)
{
    gchar *buffer = NULL;

    /* add mc has results */
    //if (stiffness_mc_has_results(&(args->stiffnessmc))) {
    buffer = stiffness_mc_export_data(&(args->stiffnessdata), &(args->stiffnessunc), &(args->stiffnessmc), &(args->fddata), EXPORT_FORMAT_PLAIN);
    buffer_to_file_dialog(buffer);
    g_free(buffer);
    //}
}

void stiffness_propagate_uncertainties(StiffnessUncControls *ctrls, const Stiffnessdata *stiffnessdata, StiffnessUncdata *unc, const FDdata *fddata)
{
    gchar str[30];
    gint i, j;
    GtkWidget *label;

    if (!stiffnessdata->has_fit_load || !stiffnessdata->has_fit_unload) {
        return;
    }

    if (unc->klc == NULL || unc->qlc == NULL || unc->kuc == NULL || unc->quc == NULL) {
        return;
    }

    // calculate the uncertainties
    stiffness_propagate_uncertainties_calc(stiffnessdata, unc, fddata);

    // show them
    // uncertainties from uh and uF
    label_set_gdouble_format(ctrls->ukl, unc->ukl, SLOPE);
    label_set_gdouble_format(ctrls->uku, unc->uku, SLOPE);
    label_set_gdouble_format(ctrls->uql, unc->uql, "%.6f");
    label_set_gdouble_format(ctrls->uqu, unc->uqu, "%.6f");

    /* contact point uncertainties */

    //destroy the labels of the old table_unc_contact
    remove_table_content(ctrls->table_unc_contact);

    //create new table
    /* table headers */

    i = 0;

    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("Contact point shift"), 0, 1, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("k<sub>load</sub>/GPa"), 1, 2, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("q<sub>load</sub>/GPa"), 2, 3, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("k<sub>load</sub>/GPa"), 3, 4, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact),
                     gwy_label_new_header("q<sub>unload</sub>/GPa"), 4, 5, i, i + 1,  GTK_EXPAND | GTK_FILL, 0, 0, 0);
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

        stiffness_fit_shift_contact(stiffnessdata, fddata, &(unc->instdata), unc->klc, unc->qlc, unc->kuc, unc->quc, j);

        //fill in kl
        label = gtk_label_new(NULL);

        if (unc->klc[j + CONTACTMAX] > -G_MAXINT) {
            g_snprintf(str, sizeof(str), SLOPE, unc->klc[j + CONTACTMAX]);
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

        //fill in ql
        label = gtk_label_new(NULL);

        if (unc->qlc[j + CONTACTMAX] > -G_MAXINT) {
            g_snprintf(str, sizeof(str), DEPTH, unc->qlc[j + CONTACTMAX]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 2, 3, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        //fill in ku
        label = gtk_label_new(NULL);

        if (unc->kuc[j + CONTACTMAX] > -G_MAXINT) {
            g_snprintf(str, sizeof(str), SLOPE, unc->kuc[j + CONTACTMAX]);
        }
        else {
            g_snprintf(str, sizeof(str), "not defined");
        }

        gtk_label_set_text(GTK_LABEL(label), str);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
        gtk_table_attach(GTK_TABLE(ctrls->table_unc_contact), label, 3, 4, i, i + 1, GTK_FILL, GTK_FILL, 0, 0);

        //fill in qu
        label = gtk_label_new(NULL);

        if (unc->quc[j + CONTACTMAX] > -G_MAXINT) {
            g_snprintf(str, sizeof(str), DEPTH, unc->quc[j + CONTACTMAX]);
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

void stiffness_mc_close_window(Data *data)
{
    gint j;

    if (data->controls->stiffnessmccontrols->has_window) {
        gtk_widget_destroy(data->controls->stiffnessmccontrols->window);

        destroy_StiffnessMCdata(&(data->args->stiffnessmc));
        data->controls->stiffnessmccontrols->has_window = FALSE;

        for (j = 0; j < NPARAM_STIFFNESS; j++) {
            if (GTK_IS_WIDGET(data->controls->stiffnessmccontrols->hist_window[j])) {
                gtk_widget_destroy(data->controls->stiffnessmccontrols->hist_window[j]);
            }
        }

        g_free(data->controls->stiffnessmccontrols->hist_window);
        g_free(data->controls->stiffnessmccontrols->button_hist);
    }
}

static void stiffness_mc_hist_kl(Data *data)
{
    show_histogram("kl", "mN/nm", data->args->stiffnessmc.mcxhist[0], data->args->stiffnessmc.mcyhist[0], data->args->stiffnessmc.mcnstat, data->controls->stiffnessmccontrols->window);
}

static void stiffness_mc_hist_ql(Data *data)
{
    show_histogram("ql", "nm", data->args->stiffnessmc.mcxhist[1], data->args->stiffnessmc.mcyhist[1], data->args->stiffnessmc.mcnstat, data->controls->stiffnessmccontrols->window);
}

static void stiffness_mc_hist_ku(Data *data)
{
    show_histogram("ku", "mN/nm", data->args->stiffnessmc.mcxhist[2], data->args->stiffnessmc.mcyhist[2], data->args->stiffnessmc.mcnstat, data->controls->stiffnessmccontrols->window);
}


static void stiffness_mc_hist_qu(Data *data)
{
    show_histogram("qu", "nm", data->args->stiffnessmc.mcxhist[3], data->args->stiffnessmc.mcyhist[3], data->args->stiffnessmc.mcnstat, data->controls->stiffnessmccontrols->window);
}

#endif

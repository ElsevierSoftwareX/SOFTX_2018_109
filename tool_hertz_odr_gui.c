#ifndef NOFORTRAN

#include "tool_hertz_odr_gui.h"
#include "tool_hertz_odr.h"
#include "tool_hertz_odr_unc.h"
#include "tool_hertz_odr_unc_gui.h"

#include "controls.h"
#include "datatypes.h"
#include "fit-utils.h"
#include "gui-utils.h"
#include "niget-common.h"

#include <math.h>

#include <stdlib.h>
#include <string.h>

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>
#include <app/file.h>

static void hertz_odr_remove_all_fit_labels(HertzODRControls *hzodrcontrols);
static void hertz_odr_remove_all_fit_curves(HertzODRControls *hzodrcontrols);

static void hz_odr_zoom_in(Data *data);
static void hz_odr_zoom_out(Data *data);


void hertz_odr_remove_fit_results_labels(Data *data)
{
    gwy_selection_clear(data->controls->hertzodrcontrols->selection);
    hertz_odr_remove_all_fit_results(&(data->args->hertzodrdata));
    hertz_odr_remove_all_fit_labels(data->controls->hertzodrcontrols);
    hertz_odr_remove_all_fit_curves(data->controls->hertzodrcontrols);
}

static void hertz_odr_remove_all_fit_labels(HertzODRControls *hzodrcontrols)
{
    /* remove all labels that result from a fit */
    /*  entry_clear(hzodrcontrols->from);
        entry_clear(hzodrcontrols->to); */

    label_clear(hzodrcontrols->gamma);
    label_clear(hzodrcontrols->h0);
    label_clear(hzodrcontrols->n);
    label_clear(hzodrcontrols->SSres);
    label_clear(hzodrcontrols->R2);
    label_clear(hzodrcontrols->R2adj);
    label_clear(hzodrcontrols->chi2);
    label_clear(hzodrcontrols->Erres);
    label_clear(hzodrcontrols->Eitres);
    label_clear(hzodrcontrols->range);
}

static void hertz_odr_remove_all_fit_curves(HertzODRControls *hzodrcontrols)
{
    gwy_graph_model_remove_curve_by_description(hzodrcontrols->graph_model, "Fit");
}

static void recalculate_Eit(HertzODRControls *hertzodrcontrols, HertzODRdata *hertzodrdata, const Instdata *instdata)
{
    if (hertzodrdata->has_fit) {
        hertzodrdata->Eit = calc_Eit(hertzodrdata->Er * 1e9, instdata->nu, instdata->Ei, instdata->nui) * 1e-9;
        label_set_gdouble_format(hertzodrcontrols->Eitres, hertzodrdata->Eit, MODUL);
    }
}

static void recalculate_Er(HertzODRControls *hertzodrcontrols, HertzODRdata *hertzodrdata)
{
    hertzodrdata->Er = calc_Er_Hertz(hertzodrdata->gamma, hertzodrdata->radius * 1e9) * 1e6;
    label_set_gdouble_format(hertzodrcontrols->Erres, hertzodrdata->Er, MODUL);
}

static void recalculate_radius(HertzODRControls *hertzodrcontrols, HertzODRdata *hertzodrdata)
{
    hertzodrdata->radius = calc_R_Hertz(hertzodrdata->gamma, hertzodrdata->Er) * 1e3;
    label_set_gdouble_format(hertzodrcontrols->radiusres, hertzodrdata->radius * 1e9, MODUL);
}

static void radius_changed(Data *data)
{
    data->args->hertzodrdata.radius = gtk_adjustment_get_value(GTK_ADJUSTMENT(gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(data->controls->hertzodrcontrols->radiusinp)))) * 1e-9;

    if (data->args->hertzodrdata.has_fit) {
        recalculate_Er(data->controls->hertzodrcontrols, &(data->args->hertzodrdata));
        recalculate_Eit(data->controls->hertzodrcontrols, &(data->args->hertzodrdata), &(data->args->instdata));
    }

    hertz_odr_propagate_uncertainties(data->controls->hertzodrunccontrols, &(data->args->hertzodrdata), &(data->args->hertzodrunc),
                                      &(data->args->fddata));
}

static void Er_changed(Data *data)
{
    data->args->hertzodrdata.Er = gtk_adjustment_get_value(GTK_ADJUSTMENT(gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(data->controls->hertzodrcontrols->Erinp))));

    if (data->args->hertzodrdata.has_fit) {
        recalculate_radius(data->controls->hertzodrcontrols, &(data->args->hertzodrdata));
    }

    hertz_odr_propagate_uncertainties(data->controls->hertzodrunccontrols, &(data->args->hertzodrdata), &(data->args->hertzodrunc),
                                      &(data->args->fddata));
}

static void Eit_changed(Data *data)
{
    gdouble nu, nui, Ei;

    data->args->hertzodrdata.Eit = gtk_adjustment_get_value(GTK_ADJUSTMENT(gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(data->controls->hertzodrcontrols->Eitinp))));

    nu = data->args->instdata.nu;
    nui = data->args->instdata.nui;
    Ei = data->args->instdata.Ei;

    data->args->hertzodrdata.Er = calc_Er(data->args->hertzodrdata.Eit, nu, Ei * 1e-9, nui);

    if (data->args->hertzodrdata.has_fit) {
        recalculate_radius(data->controls->hertzodrcontrols, &(data->args->hertzodrdata));
    }

    hertz_odr_propagate_uncertainties(data->controls->hertzodrunccontrols, &(data->args->hertzodrdata), &(data->args->hertzodrunc),
                                      &(data->args->fddata));
}

static void hertz_odr_radio_callback(GtkWidget *widget, Data *data)
{
    HertzODRControls *hzodrcontrols;

    hzodrcontrols = data->controls->hertzodrcontrols;

    if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget))) {
        gtk_label_set_text(GTK_LABEL(data->controls->hertzodrcontrols->radiusres), NULL);
        gtk_label_set_text(GTK_LABEL(data->controls->hertzodrcontrols->Erres), NULL);
        gtk_label_set_text(GTK_LABEL(data->controls->hertzodrcontrols->Eitres), NULL);
        return;
    }

    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(hzodrcontrols->radio_R))) {
        data->args->hertzodrdata.mode = R_MODE;
        gtk_widget_set_sensitive(hzodrcontrols->radiusinp, TRUE);
        gtk_widget_set_sensitive(hzodrcontrols->Erinp, FALSE);
        gtk_widget_set_sensitive(hzodrcontrols->Eitinp, FALSE);
        radius_changed(data);
    }
    else if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(hzodrcontrols->radio_Er))) {
        data->args->hertzodrdata.mode = ER_MODE;
        gtk_widget_set_sensitive(hzodrcontrols->radiusinp, FALSE);
        gtk_widget_set_sensitive(hzodrcontrols->Erinp, TRUE);
        gtk_widget_set_sensitive(hzodrcontrols->Eitinp, FALSE);
        Er_changed(data);
    }
    else if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(hzodrcontrols->radio_Eit))) {
        data->args->hertzodrdata.mode = EIT_MODE;
        gtk_widget_set_sensitive(hzodrcontrols->radiusinp, FALSE);
        gtk_widget_set_sensitive(hzodrcontrols->Erinp, FALSE);
        gtk_widget_set_sensitive(hzodrcontrols->Eitinp, TRUE);
        Eit_changed(data);
    }
    else {
        g_printerr("You should never get here! Something weird is going on \n");
    }
}

static void toggle_radial_changed(GtkWidget *button, Data *data)
{
    HertzODRdata *hzodr;
    HertzODRControls *hzctrls;

    hzodr = &(data->args->hertzodrdata);
    hzctrls = data->controls->hertzodrcontrols;

    hzodr->radial_corr = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));

    hertz_odr_remove_all_fit_labels(hzctrls);
    //hertz_odr_remove_all_fit_results(hzodr);
    hertz_odr_remove_all_fit_curves(hzctrls);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(hzctrls->fixgamma), TRUE);
}

static void fixgamma_changed(GtkWidget *button, Data *data)
{
    /*fix gamma 1-fitted, 0-fixed */
    HertzODRdata *hzodr;
    HertzODRControls *hzctrls;

    hzodr = &(data->args->hertzodrdata);
    hzctrls = data->controls->hertzodrcontrols;

    hzodr->ifixb[2] = 1 - gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));

    hertz_odr_remove_all_fit_labels(hzctrls);
    //hertz_odr_remove_all_fit_results(hzodr);
    hertz_odr_remove_all_fit_curves(hzctrls);

    if (verbose) {
        if (hzodr->ifixb[2]) {
            g_print("Gamma fitted \n");
        }
        else {
            g_print("Gamma fixed \n");
        }
    }
}

void hertz_odr_redraw(Data *data)
{
    // update loading curve and Fmax and hmax labels
    if (data->args->fddata.hload) {
        gwy_graph_curve_model_set_data(data->controls->hertzodrcontrols->cmodelload,
                                       data->args->fddata.hload->data, data->args->fddata.Fload->data,
                                       data->args->fddata.hload->res);

        label_set_gdouble_format(data->controls->hertzodrcontrols->hmax, data->args->fddata.hmax, DEPTH);
        label_set_gdouble_format(data->controls->hertzodrcontrols->Fmax, data->args->fddata.Fmax, FORCE);

        if (gwy_graph_model_get_n_curves(data->controls->hertzodrcontrols->graph_model) == 2) {
            gwy_graph_curve_model_set_data(data->controls->hertzodrcontrols->cmodelfit,
                                           data->args->hertzodrdata.xfit->data, data->args->hertzodrdata.yfit->data,  data->args->hertzodrdata.xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(data->controls->hertzodrcontrols->cmodelfit);
#endif
        }
    }
}

static void hz_odr_range_from_changed(Data *data)
{
    HertzODRControls *fcontrols;
    HertzODRdata *fdata;

    fcontrols = data->controls->hertzodrcontrols;
    fdata = &(data->args->hertzodrdata);

    range_from_changed(fcontrols->from, fcontrols->to, fcontrols->selection, &(fdata->from), &(fdata->from));
}

static void hz_odr_range_to_changed(Data *data)
{
    HertzODRControls *tcontrols;
    HertzODRdata *tdata;

    tcontrols = data->controls->hertzodrcontrols;
    tdata = &(data->args->hertzodrdata);

    range_to_changed(tcontrols->from, tcontrols->to, tcontrols->selection, &(tdata->to), &(tdata->to));
}

static void hertz_odr_graph_selected(GwySelection *selection, gint i, Data *data)
{
    gdouble range[2];

    g_return_if_fail(i <= 0);

    //take data from selection and put them in the "from" and "to" values
    gwy_selection_get_object(selection, 0, range);

    data->args->hertzodrdata.from = MIN(range[0], range[1]);
    data->args->hertzodrdata.to = MAX(range[0], range[1]);

    //update labels
    entry_set_gdouble_format(data->controls->hertzodrcontrols->from, data->args->hertzodrdata.from, DEPTH);
    entry_set_gdouble_format(data->controls->hertzodrcontrols->to, data->args->hertzodrdata.to, DEPTH);
}

static void hertz_odr_save_data(Args *args)
{
    gchar *buffer = NULL;

    if (hertz_odr_has_results(&(args->hertzodrdata))) {
        buffer = hertz_odr_export_data(&(args->hertzodrdata), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN);
        buffer_to_file_dialog(buffer);
        g_free(buffer);
    }
}

static void hertz_odr_show_log(Data *data)
{
    gint infolog;

    infolog = data->args->hertzodrdata.infolog;

    /* check if there are logs */
    switch (infolog) {
    case 0:
        show_log(data->controls->hertzodrcontrols->log_error_window, data->args->hertzodrdata.logfnm, ".err", "Fit error log");
        show_log(data->controls->hertzodrcontrols->log_report_window, data->args->hertzodrdata.logfnm, ".rpt", "Fit report log");
        break;

    case 1:
        show_log(data->controls->hertzodrcontrols->log_error_window,  data->args->hertzodrdata.logfnm, ".err", "Fit error log");
        show_warning("Fit report lpg could not be created.");
        break;

    case 2:
        show_log(data->controls->hertzodrcontrols->log_report_window, data->args->hertzodrdata.logfnm, ".rpt", "Fit report log");
        show_warning("Fit error log could not be created.");
        break;

    case 3:
        show_warning("Fit error log could not be created.");
        show_warning("Fit report log could not be created.");
        break;

    default:
        g_printerr("Should not get here.\n");
        break;
    }
}

void hertz_odr_fit_and_update(Data *data)
{
    gchar str[30];
    gint nc, i;
    gint istart, iend, ndata;
    gint nremove;
    gint fitinfo;
    gboolean fit = TRUE;
    GwyDataLine *x, *y;
    gdouble b, w;

    HertzODRdata *hzodrdata;
    FDdata *fddata;
    Instdata *instdata;
    HertzODRControls *hzodrcontrols;
    GdkColor color;

    gdk_color_parse("red", &color);

    // TEMPORARY update GtkEntries for ranges
    /*
    range_from_changed(data);
    range_to_changed(data);
    */

    //create pointers for easier reading
    hzodrdata = &(data->args->hertzodrdata);
    hzodrcontrols = data->controls->hertzodrcontrols;
    fddata = &(data->args->fddata);
    instdata = &(data->args->instdata);

    //if there was a previous fit, invalidate previous fitted curve
    if (hzodrdata->has_fit) {
        g_object_unref(hzodrdata->xfit);
        g_object_unref(hzodrdata->yfit);
    }

    // check range
    if (hzodrdata->from == hzodrdata->to) {
        fit = FALSE;
    }

    //find appropriate data, which should be fitted
    range_to_indices(hzodrdata->from, hzodrdata->to, fddata->hload, FALSE, &istart, &iend, &ndata);
    // if there are not enough data in the selected range, return fail
    ndata = MIN(ndata, fddata->hload->res - istart);

    if (ndata < 2) {
        fit = FALSE;
    }

    hzodrdata->nfitdata = ndata;

    if (fit) {
        //get data to be fitted
        x = gwy_data_line_part_extract(fddata->hload, istart, ndata);
        y = gwy_data_line_part_extract(fddata->Fload, istart, ndata);

        //remove negative depths from data
        nremove = filter_negative_data(x, y);
        hzodrdata->nfitdata -= nremove;

        //fit the data and evaluate
        fitinfo = hertz_odr_fit(x, y, hzodrdata, fddata, instdata, NULL);

        if (hzodrdata->radial_corr) {
            radial_correction_hertz(hzodrdata->mode, &(hzodrdata->Er), &(hzodrdata->Eit), &(hzodrdata->radius), hzodrdata->to, &(data->args->instdata));
        }

        if (verbose) {
            g_print("fitinfo pre %d infolog %d \n", fitinfo, hzodrdata->infolog);
        }

        if (fitinfo >= 1e6) {
            hzodrdata->infolog = fitinfo / 1000000;
            fitinfo = fitinfo % 1000000;
        }

        if (verbose) {
            g_print("fitinfo post %d infolog %d \n", fitinfo, hzodrdata->infolog);
        }

        if (fitinfo >= 5) {
            gtk_label_set_text(GTK_LABEL(hzodrcontrols->fitinfo), "Warning: Questionable results \n or fatal errors in fitting!");
            gtk_widget_modify_fg(hzodrcontrols->fitinfo, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(hzodrcontrols->fitinfo), TRUE);
        }
        else if (fitinfo == 4) {
            gtk_label_set_text(GTK_LABEL(hzodrcontrols->fitinfo), "Warning: Iteration limit reached!");
            gtk_widget_modify_fg(hzodrcontrols->fitinfo, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(hzodrcontrols->fitinfo), TRUE);
        }
        else if (hzodrdata->n < 1 && !hzodrdata->radial_corr) {
            gtk_label_set_text(GTK_LABEL(hzodrcontrols->fitinfo), "Warning: Unphysical fit results!");
            gtk_widget_modify_fg(hzodrcontrols->fitinfo, GTK_STATE_NORMAL, &color);
            gtk_label_set_line_wrap(GTK_LABEL(hzodrcontrols->fitinfo), TRUE);
        }
        else {
            label_clear(hzodrcontrols->fitinfo);
        }

        //update range labels
        entry_set_gdouble_format(hzodrcontrols->from, hzodrdata->from, DEPTH);
        entry_set_gdouble_format(hzodrcontrols->to, hzodrdata->to, DEPTH);
        g_snprintf(str, sizeof(str), "[%.3f, %.3f]", hzodrdata->from, hzodrdata->to);
        gtk_label_set_text(GTK_LABEL(hzodrcontrols->range), str);

        label_set_gdouble_format(hzodrcontrols->gamma, hzodrdata->gamma, "%.4g");
        label_set_gdouble_format(hzodrcontrols->h0, hzodrdata->h0, DEPTH);
        label_set_gdouble_format(hzodrcontrols->n, hzodrdata->n, NUMBER);

        label_set_gdouble_format(hzodrcontrols->SSres, hzodrdata->SSres, "%.4g");
        label_set_gdouble_format(hzodrcontrols->R2, hzodrdata->R2, "%.4g");
        label_set_gdouble_format(hzodrcontrols->R2adj, hzodrdata->R2adj, "%.4g");
        label_set_gdouble_format(hzodrcontrols->chi2, hzodrdata->chi2, "%.4g");

        if (hzodrdata->mode == R_MODE) {
            label_set_gdouble_format(hzodrcontrols->Eitres, hzodrdata->Eit, MODUL);
            label_set_gdouble_format(hzodrcontrols->Erres, hzodrdata->Er, MODUL);
        }
        else {
            label_set_gdouble_format(hzodrcontrols->radiusres, hzodrdata->radius * 1e9, DEPTH);

        }

        // create fit curve data
        hzodrdata->xfit = gwy_data_line_duplicate(x);
        hzodrdata->yfit = gwy_data_line_part_extract(fddata->hload, 0, hzodrdata->nfitdata);

        if (hzodrdata->radial_corr) {
            b = 2. / 3. / M_PI * (1 - 2 * instdata->nu) / (1 - instdata->nu) / sqrt(hzodrdata->radius * 1e9); /* must be in nanometers*/

            for (i = 0; i < hzodrdata->xfit->res; i++) {
                w = sqrt(hzodrdata->xfit->data[i] - hzodrdata->h0);
                hzodrdata->yfit->data[i] =  hzodrdata->gamma * (1 + b * w) * w * w * w;
            }
        }
        else {
            for (i = 0; i < hzodrdata->xfit->res; i++) {
                hzodrdata->yfit->data[i] =  hzodrdata->gamma * pow(hzodrdata->xfit->data[i] - hzodrdata->h0, hzodrdata->n);
            }
        }

        //clean up
        g_object_unref(x);
        g_object_unref(y);

        // update fitted curve, if necessary create it.
        nc = gwy_graph_model_get_n_curves(hzodrcontrols->graph_model);

        if (nc == 2) {
            gwy_graph_curve_model_set_data(hzodrcontrols->cmodelfit, hzodrdata->xfit->data, hzodrdata->yfit->data, hzodrdata->xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(hzodrcontrols->cmodelfit);
#endif
        }
        else {
            gwy_graph_curve_model_set_data(hzodrcontrols->cmodelfit, hzodrdata->xfit->data, hzodrdata->yfit->data, hzodrdata->xfit->res);
#if (GWY_VERSION_MAJOR == 2) && (GWY_VERSION_MINOR < 45)
            gwy_graph_curve_model_enforce_order(hzodrcontrols->cmodelfit);
#endif
            gwy_graph_model_add_curve(hzodrcontrols->graph_model, hzodrcontrols->cmodelfit);
            g_object_set(hzodrcontrols->cmodelfit,
                         "mode", GWY_GRAPH_CURVE_LINE,
                         "description", "Fit",
                         "color", CURVE_COLOR_FIT,
                         "line-style", GDK_LINE_SOLID,
                         "line-width", 3,
                         NULL);
        }
    }
    else {
        nc = gwy_graph_model_get_n_curves(hzodrcontrols->graph_model);
        hertz_odr_remove_all_fit_labels(hzodrcontrols);
        hertz_odr_remove_all_fit_results(hzodrdata);
        hzodrdata->has_fit = FALSE;

        if (nc == 2) {
            gwy_graph_model_remove_curve_by_description(hzodrcontrols->graph_model, "Fit");
        }
    }

    //
    // log of fitting procedure is available
    if (hzodrdata->has_fit) {
        gtk_widget_set_sensitive(hzodrcontrols->button_log, TRUE);
    }
    else {
        gtk_widget_set_sensitive(hzodrcontrols->button_log, FALSE);
    }

    // update previous uncertainties if existing, set button
    if (hzodrdata->has_fit) {
        hertz_odr_propagate_uncertainties(data->controls->hertzodrunccontrols, &(data->args->hertzodrdata), &(data->args->hertzodrunc),
                                          &(data->args->fddata));
        gtk_widget_set_sensitive(hzodrcontrols->button_unc, TRUE);
    }
    else {
        gtk_widget_set_sensitive(hzodrcontrols->button_unc, FALSE);
    }
}

GtkWidget *tool_hertz_odr_create(Data *data)
{
    gchar str[30];
    GtkWidget *label;
    GtkObject *adj;
    GwyGraphArea *area;
    Args *args;
    HertzODRControls *hzodrcontrols;
    gint i;

    args = data->args;
    hzodrcontrols = data->controls->hertzodrcontrols;

    /* Log, uncertainties, save buttons */

    hzodrcontrols->button_log = gtk_button_new_with_mnemonic("Show _log");
    gtk_button_set_image(GTK_BUTTON(hzodrcontrols->button_log), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(hzodrcontrols->button_log, "clicked", G_CALLBACK(hertz_odr_show_log), data);
    gtk_widget_set_sensitive(hzodrcontrols->button_log, FALSE);
    hzodrcontrols->log_error_window = NULL;
    hzodrcontrols->log_report_window = NULL;

    hzodrcontrols->button_save = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    g_signal_connect_swapped(hzodrcontrols->button_save, "clicked", G_CALLBACK(hertz_odr_save_data), args);

    hzodrcontrols->button_unc = gtk_button_new_with_mnemonic("_Uncertainties");
    gtk_button_set_image(GTK_BUTTON(hzodrcontrols->button_unc), gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(hzodrcontrols->button_unc, "clicked", G_CALLBACK(hertz_odr_uncertainty), data);
    gtk_widget_set_sensitive(hzodrcontrols->button_unc, FALSE);

    // GraphArea, two curves will be displayed: loading curve (F-d) and the resulting fit
    hzodrcontrols->graph_model = gwy_graph_model_new();
    g_object_set(hzodrcontrols->graph_model, "title", "Load", NULL);
    gwy_graph_model_set_axis_label(hzodrcontrols->graph_model, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(hzodrcontrols->graph_model, GTK_POS_LEFT, "F / mN");

    hzodrcontrols->zoomstack = NULL;
    hzodrcontrols->zoomstack = g_slist_prepend(hzodrcontrols->zoomstack, hzodrcontrols->graph_model);

    hzodrcontrols->hbox_graph = gtk_hbox_new(FALSE, CTRLS_SPACING);

    hzodrcontrols->graph = gwy_graph_new(hzodrcontrols->graph_model);
    /* gtk_widget_set_size_request(GTK_WIDGET(hzodrcontrols->graph), 600, 400); */

    gwy_graph_enable_user_input(GWY_GRAPH(hzodrcontrols->graph), FALSE);
    gwy_graph_set_status(GWY_GRAPH(hzodrcontrols->graph), GWY_GRAPH_STATUS_XSEL);

    area = GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(hzodrcontrols->graph)));
    hzodrcontrols->selection = gwy_graph_area_get_selection(area, GWY_GRAPH_STATUS_XSEL);
    gwy_selection_set_max_objects(hzodrcontrols->selection, 1);
    g_signal_connect(hzodrcontrols->selection, "changed", G_CALLBACK(hertz_odr_graph_selected), data);

    hzodrcontrols->cmodelload = gwy_graph_curve_model_new();
    hzodrcontrols->cmodelfit = gwy_graph_curve_model_new();

    g_object_set(hzodrcontrols->cmodelload,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", "Loading",
                 "point-size", 0,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_LOAD,
                 "line-style", GDK_LINE_SOLID,
                 "line-width", 2,
                 NULL);
    gwy_graph_model_add_curve(hzodrcontrols->graph_model, hzodrcontrols->cmodelload);


    /* Ctrls */

    hzodrcontrols->ctrls = gtk_vbox_new(FALSE, CTRLS_SPACING);

    /* Info section */

    hzodrcontrols->frame_info = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(hzodrcontrols->frame_info), GTK_WIDGET(gwy_label_new_header("Info")));
    hzodrcontrols->table_info = gtk_table_new(2, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(hzodrcontrols->table_info), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(hzodrcontrols->table_info), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(hzodrcontrols->table_info), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(hzodrcontrols->frame_info), hzodrcontrols->table_info);

    i = 0;

    /* h_max */
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_info), label_new_with_markup_left("h<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    hzodrcontrols->hmax = gtk_label_new(NULL);
    label_set_gdouble_format(hzodrcontrols->hmax, args->fddata.hmax, DEPTH);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->hmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_info), hzodrcontrols->hmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_info), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    /* F_max */
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_info), label_new_with_markup_left("F<sub>max</sub>"), 0, 1, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    hzodrcontrols->Fmax = gtk_label_new(NULL);
    label_set_gdouble_format(hzodrcontrols->Fmax, args->fddata.Fmax, FORCE);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->Fmax), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_info), hzodrcontrols->Fmax, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_info), label_new_with_markup_left("mN"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);

    /* Parameters section */

    hzodrcontrols->frame_parameters = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(hzodrcontrols->frame_parameters), GTK_WIDGET(gwy_label_new_header("Parameters")));
    hzodrcontrols->table_parameters = gtk_table_new(4, 4, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(hzodrcontrols->table_parameters), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(hzodrcontrols->table_parameters), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(hzodrcontrols->table_parameters), CTRLS_TABLE_ROW_SPACING);
    gtk_container_add(GTK_CONTAINER(hzodrcontrols->frame_parameters), hzodrcontrols->table_parameters);

    i = 0;

    hzodrcontrols->radio_R   = gtk_radio_button_new_with_label(NULL, "Radius");
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), hzodrcontrols->radio_R, 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(G_OBJECT(hzodrcontrols->radio_R), "toggled", G_CALLBACK(hertz_odr_radio_callback), data);
    /*hzodrcontrols->radiusinp = gtk_adjustment_new(args->hertzodrdata.radius*1e9, 1e-9, 1e6, 1, 10,0);
    gtk_widget_set_sensitive(hzodrcontrols->radiusinp, TRUE);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(hzodrcontrols->radiusinp), 0.1, 1);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), spin, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    */
    adj = gtk_adjustment_new(args->hertzodrdata.radius * 1e9, 1e-9, 1e6, 1, 10, 0);
    hzodrcontrols->radiusinp = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 0.1, 1);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(hzodrcontrols->radiusinp), TRUE);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), hzodrcontrols->radiusinp, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_widget_set_sensitive(hzodrcontrols->radiusinp, TRUE);
    data->args->hertzodrdata.mode = R_MODE;
    g_signal_connect_swapped(hzodrcontrols->radiusinp, "value-changed", G_CALLBACK(radius_changed), data);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    hzodrcontrols->radio_Er  = gtk_radio_button_new_from_widget(GTK_RADIO_BUTTON(hzodrcontrols->radio_R));
    label = label_new_with_markup("E<sub>r</sub>");
    gtk_container_add(GTK_CONTAINER(hzodrcontrols->radio_Er), label);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), hzodrcontrols->radio_Er, 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(G_OBJECT(hzodrcontrols->radio_Er), "toggled", G_CALLBACK(hertz_odr_radio_callback), data);
    /*
    hzodrcontrols->Erinp = gtk_adjustment_new(args->hertzodrdata.Er*1e-9, 1e-9, 1e6, 1, 10,0);
    gtk_widget_set_sensitive(hzodrcontrols->Erinp, FALSE);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(hzodrcontrols->Erinp), 0.1, 1);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), spin, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(hzodrcontrols->Erinp, "value-changed", G_CALLBACK(Er_changed), data);
    */
    adj = gtk_adjustment_new(10, 1e-3, 1e6, 1, 10, 0);
    hzodrcontrols->Erinp = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 0.001, 3);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(hzodrcontrols->Erinp), TRUE);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), hzodrcontrols->Erinp, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_widget_set_sensitive(hzodrcontrols->Erinp, FALSE);
    g_signal_connect_swapped(hzodrcontrols->Erinp, "value-changed", G_CALLBACK(Er_changed), data);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    hzodrcontrols->radio_Eit  = gtk_radio_button_new_from_widget(GTK_RADIO_BUTTON(hzodrcontrols->radio_R));
    label = label_new_with_markup("E<sub>IT</sub>");
    gtk_container_add(GTK_CONTAINER(hzodrcontrols->radio_Eit), label);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), hzodrcontrols->radio_Eit, 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(G_OBJECT(hzodrcontrols->radio_Eit), "toggled", G_CALLBACK(hertz_odr_radio_callback), data);
    /*
    hzodrcontrols->Eitinp = gtk_adjustment_new(args->hertzodrdata.Eit*1e-9, 1e-9, 1e6, 1, 10,0);
    gtk_widget_set_sensitive(GTK_WIDGET(hzodrcontrols->Eitinp), FALSE);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(hzodrcontrols->Eitinp), 0.1, 1);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), spin, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(hzodrcontrols->Eitinp, "value-changed", G_CALLBACK(Eit_changed), data);
    */
    adj = gtk_adjustment_new(10, 1e-3, 1e6, 1, 10, 0);
    hzodrcontrols->Eitinp = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 0.001, 3);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(hzodrcontrols->Eitinp), TRUE);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), hzodrcontrols->Eitinp, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_widget_set_sensitive(hzodrcontrols->Eitinp, FALSE);
    g_signal_connect_swapped(hzodrcontrols->Eitinp, "value-changed", G_CALLBACK(Eit_changed), data);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), label_new_with_markup_left("Fit range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->from = gtk_entry_new();
    g_object_set_data(G_OBJECT(hzodrcontrols->from), "id", (gpointer)"from");
    gtk_entry_set_width_chars(GTK_ENTRY(hzodrcontrols->from), 8);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), hzodrcontrols->from, 1, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(hzodrcontrols->from, TRUE);
    g_signal_connect_swapped(hzodrcontrols->from, "activate", G_CALLBACK(hz_odr_range_from_changed), data);
    /*	g_signal_connect_swapped(hzodrcontrols->from, "focus-out-event", G_CALLBACK(range_from_changed),data); */

    hzodrcontrols->to = gtk_entry_new();
    g_object_set_data(G_OBJECT(hzodrcontrols->to), "id", (gpointer)"to");
    gtk_entry_set_width_chars(GTK_ENTRY(hzodrcontrols->to), 8);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), hzodrcontrols->to, 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(hzodrcontrols->to, TRUE);
    g_signal_connect_swapped(hzodrcontrols->to, "activate", G_CALLBACK(hz_odr_range_to_changed), data);
    /*	g_signal_connect_swapped(hzodrcontrols->to, "focus-out-event", G_CALLBACK(range_to_changed),data);*/

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), label_new_with_markup_left("nm"), 3, 4, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    g_snprintf(str, sizeof(str), "Fix exponent at %g", N_FIX_HZODR);
    hzodrcontrols->fixgamma = gtk_check_button_new_with_label(str);
    //	hzodrcontrols->fixgamma = gtk_check_button_new_with_label("Fix exponent at %g");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(hzodrcontrols->fixgamma), FALSE);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), hzodrcontrols->fixgamma, 0, 2, i, i + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(hzodrcontrols->fixgamma, "toggled", G_CALLBACK(fixgamma_changed), data);
    i++;

    /* hzodrcontrols->toggle_radial = gtk_check_button_new_with_label("Radial correction (fixed exponent)"); */
    /* gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(hzodrcontrols->toggle_radial), FALSE); */
    /* gtk_table_attach(GTK_TABLE(hzodrcontrols->table_parameters), hzodrcontrols->toggle_radial, 0, 4, i, i + 1, GTK_FILL, 0, 0, 0); */
    /* g_signal_connect(hzodrcontrols->toggle_radial, "toggled", G_CALLBACK(toggle_radial_changed), data); */

    /* Fit button */
    hzodrcontrols->button_fit = gtk_button_new_with_mnemonic(gwy_sgettext("verb|_Fit"));
    gtk_button_set_image(GTK_BUTTON(hzodrcontrols->button_fit), gtk_image_new_from_stock(GTK_STOCK_EXECUTE, GTK_ICON_SIZE_BUTTON));
    g_signal_connect_swapped(hzodrcontrols->button_fit, "clicked", G_CALLBACK(hertz_odr_fit_and_update), data);

    /* Results section */

    hzodrcontrols->frame_results = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(hzodrcontrols->frame_results), GTK_WIDGET(gwy_label_new_header("Results")));

    hzodrcontrols->scrolledwindow_results = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(hzodrcontrols->scrolledwindow_results), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(hzodrcontrols->frame_results), hzodrcontrols->scrolledwindow_results);

    hzodrcontrols->table_results = gtk_table_new(8, 3, FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(hzodrcontrols->table_results), CTRLS_BORDER);
    gtk_table_set_col_spacings(GTK_TABLE(hzodrcontrols->table_results), CTRLS_TABLE_COL_SPACING);
    gtk_table_set_row_spacings(GTK_TABLE(hzodrcontrols->table_results), CTRLS_TABLE_ROW_SPACING);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(hzodrcontrols->scrolledwindow_results), hzodrcontrols->table_results);

    i = 0;

    hzodrcontrols->fitinfo = gtk_label_new(NULL);
    gtk_label_set_line_wrap(GTK_LABEL(hzodrcontrols->fitinfo), TRUE);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->fitinfo, 0, 3, i, i + 2, GTK_FILL, GTK_FILL, 0, 0);
    i += 2;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("γ"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->gamma = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->gamma), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->gamma, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("h<sub>0</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->h0 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->h0), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->h0, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("n"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->n = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->n), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->n, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    /*
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("Tip radius"), 0, 1, i, i+1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->radius = gtk_adjustment_new(args->hertzodrdata.radius*1e9, 1e-9, 1e6, 1, 10,0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(hzodrcontrols->radius), 0.1, 1);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), spin, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(hzodrcontrols->radius, "value-changed", G_CALLBACK(radius_changed), data);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i+1, GTK_FILL, 0, 0, 0);
    i++;
    */

    //	gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("χ<sup>2</sup> fit"), 0, 1, i, i+1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->chi2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->chi2), 0.0, 0.5);
    //	gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->chi2, 1, 2, i, i+1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("SS<sub>res</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->SSres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->SSres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->SSres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("R<sup>2</sup>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->R2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->R2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->R2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("R<sup>2</sup><sub>adj</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->R2adj = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->R2adj), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->R2adj, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("χ<sup>2</sup> fit"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->chi2 = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->chi2), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->chi2, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("Tip radius"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->radiusres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->radiusres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->radiusres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("E<sub>r</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->Erres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->Erres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->Erres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("E<sub>IT</sub>"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->Eitres = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->Eitres), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->Eitres, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("GPa"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);
    i++;

    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("Fitted range"), 0, 1, i, i + 1, GTK_FILL, 0, 0, 0);
    hzodrcontrols->range = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(hzodrcontrols->range), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), hzodrcontrols->range, 1, 2, i, i + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_table_attach(GTK_TABLE(hzodrcontrols->table_results), label_new_with_markup_left("nm"), 2, 3, i, i + 1, GTK_FILL, 0, 0, 0);


    //final assembly

    gtk_box_pack_start(GTK_BOX(hzodrcontrols->ctrls), hzodrcontrols->frame_info, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hzodrcontrols->ctrls), hzodrcontrols->frame_parameters, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hzodrcontrols->ctrls), hzodrcontrols->button_fit, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hzodrcontrols->ctrls), hzodrcontrols->frame_results, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(hzodrcontrols->ctrls), hzodrcontrols->button_log, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hzodrcontrols->ctrls), hzodrcontrols->button_unc, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hzodrcontrols->ctrls), hzodrcontrols->button_save, FALSE, FALSE, 0);


    /* Zooming */

    hzodrcontrols->vbox_zoom_buttons = gtk_vbox_new(TRUE, 0);
    hzodrcontrols->vbox_zoom_outer = gtk_vbox_new(FALSE, 0);

    hzodrcontrols->button_zoom_in = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(hzodrcontrols->button_zoom_in), gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(hzodrcontrols->vbox_zoom_buttons), hzodrcontrols->button_zoom_in, TRUE, TRUE, 0);
    g_signal_connect_swapped(hzodrcontrols->button_zoom_in, "clicked", G_CALLBACK(hz_odr_zoom_in), data);

    hzodrcontrols->button_zoom_out = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(hzodrcontrols->button_zoom_out), gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(hzodrcontrols->vbox_zoom_buttons), hzodrcontrols->button_zoom_out, TRUE, TRUE, 0);
    g_signal_connect_swapped(hzodrcontrols->button_zoom_out, "clicked", G_CALLBACK(hz_odr_zoom_out), data);

    hzodrcontrols->button_zoom_restore = gtk_button_new();
    gtk_button_set_image(GTK_BUTTON(hzodrcontrols->button_zoom_restore), gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT, GTK_ICON_SIZE_BUTTON));
    gtk_box_pack_start(GTK_BOX(hzodrcontrols->vbox_zoom_buttons), hzodrcontrols->button_zoom_restore, TRUE, TRUE, 0);
    g_signal_connect_swapped(hzodrcontrols->button_zoom_restore, "clicked", G_CALLBACK(hz_odr_zoom_restore), data);

    gtk_box_pack_start(GTK_BOX(hzodrcontrols->vbox_zoom_outer), hzodrcontrols->vbox_zoom_buttons, TRUE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(hzodrcontrols->hbox_graph), hzodrcontrols->graph, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(hzodrcontrols->hbox_graph), hzodrcontrols->vbox_zoom_outer, FALSE, FALSE, 0);

    //zoom and unzoom buttons
    /* hzodrcontrols->button_zoom_in = gtk_button_new(); */
    /* gtk_button_set_image (GTK_BUTTON (hzodrcontrols->button_zoom_in), gtk_image_new_from_stock (GWY_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON)); */
    /* g_signal_connect_swapped(hzodrcontrols->button_zoom_in, "clicked", G_CALLBACK (zoom_in), data); */

    /* hzodrcontrols->button_zoom_out = gtk_button_new(); */
    /* gtk_button_set_image (GTK_BUTTON (hzodrcontrols->button_zoom_out), gtk_image_new_from_stock (GWY_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON)); */
    /* g_signal_connect_swapped(hzodrcontrols->button_zoom_out, "clicked", G_CALLBACK (zoom_out), data); */

    /* vbox = gtk_vbox_new(FALSE,FALSE); */
    /* gtk_box_pack_start (GTK_BOX(vbox), hzodrcontrols->graph, TRUE, TRUE, 0); */
    /* gtk_box_pack_start (GTK_BOX(vbox), hzodrcontrols->button_zoom_in, FALSE, FALSE, 0); */
    /* gtk_box_pack_start (GTK_BOX(vbox), hzodrcontrols->button_zoom_out, FALSE, FALSE, 0); */

    hzodrcontrols->hertz_odr_gui = gtk_hbox_new(FALSE, CTRLS_SPACING);

    gtk_box_pack_start(GTK_BOX(hzodrcontrols->hertz_odr_gui), hzodrcontrols->ctrls, FALSE, FALSE, 0);   /* expand, fill, padding */
    gtk_box_pack_start(GTK_BOX(hzodrcontrols->hertz_odr_gui), hzodrcontrols->hbox_graph, TRUE, TRUE, 0);

    gtk_widget_show(hzodrcontrols->hertz_odr_gui);

    return hzodrcontrols->hertz_odr_gui;
}

void hertz_odr_recalculate(Data *data)
{
    if (data->args->hertzodrdata.mode == R_MODE) {
        recalculate_Eit(data->controls->hertzodrcontrols, &(data->args->hertzodrdata), &(data->args->instdata));
    }
    else if (data->args->hertzodrdata.mode == EIT_MODE) {
        data->args->hertzodrdata.Er = 1. / ((1 - data->args->instdata.nui * data->args->instdata.nui) / (data->args->instdata.Ei * 1e-9) +
                                            (1 - data->args->instdata.nu * data->args->instdata.nu) / data->args->hertzodrdata.Eit);
        recalculate_radius(data->controls->hertzodrcontrols, &(data->args->hertzodrdata));
    }

    hertz_odr_propagate_uncertainties(data->controls->hertzodrunccontrols, &(data->args->hertzodrdata), &(data->args->hertzodrunc),
                                      &(data->args->fddata));
}

static void update_models_hz(HertzODRControls *hzodrcontrols)
{
    GwyGraphCurveModel *cmodel;

    hzodrcontrols->graph_model = gwy_graph_get_model(GWY_GRAPH(hzodrcontrols->graph));
    cmodel = gwy_graph_model_get_curve_by_description(hzodrcontrols->graph_model, "Loading");

    if (cmodel) {
        hzodrcontrols->cmodelload = cmodel;
    }

    cmodel = gwy_graph_model_get_curve_by_description(hzodrcontrols->graph_model, "Fit");

    if (cmodel) {
        hzodrcontrols->cmodelfit = cmodel;
    }
}

static void hz_odr_zoom_in(Data *data)
{
    HertzODRControls *hzodrcontrols;

    hzodrcontrols = data->controls->hertzodrcontrols;

    rs_zoom_in(GWY_GRAPH(hzodrcontrols->graph), &(hzodrcontrols->zoomstack), FALSE);   /* do not rescale y */
    update_models_hz(hzodrcontrols);
}

static void hz_odr_zoom_out(Data *data)
{
    HertzODRControls *hzodrcontrols;

    hzodrcontrols = data->controls->hertzodrcontrols;

    rs_zoom_out(GWY_GRAPH(hzodrcontrols->graph), &(hzodrcontrols->zoomstack));
    update_models_hz(hzodrcontrols);
}

void hz_odr_zoom_restore(Data *data)
{
    HertzODRControls *hzodrcontrols;

    hzodrcontrols = data->controls->hertzodrcontrols;

    rs_zoom_restore(GWY_GRAPH(hzodrcontrols->graph), &(hzodrcontrols->zoomstack));
    update_models_hz(hzodrcontrols);
}

#endif

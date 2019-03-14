#include "gui-utils.h"

#include "file-utils.h"
#include "niget-common.h"

#include <stdlib.h>
#include <string.h>

#include <gtk/gtk.h>
#include <app/gwyapp.h>
#include <libgwydgets/gwydgets.h>

/* adapted from Gwydgets */
GtkWidget *label_new_left(const gchar *text)
{
    GtkWidget *label;

    label = gtk_label_new(text);
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);

    return label;
}

GtkWidget *label_new_with_markup(const gchar *text)
{
    GtkWidget *label;

    label = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(label), text);

    return label;
}

GtkWidget *label_new_with_markup_with_mnemonic(const gchar *text)
{
    GtkWidget *label;

    label = gtk_label_new(NULL);
    gtk_label_set_markup_with_mnemonic(GTK_LABEL(label), text);

    return label;
}

GtkWidget *label_new_with_markup_left(const gchar *text)
{
    GtkWidget *label;

    label = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(label), text);
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);

    return label;
}

void remove_table_content(GtkWidget *table)
{
    GList *children, *iter;

    children = gtk_container_get_children(GTK_CONTAINER(table));

    for (iter = children; iter != NULL; iter = g_list_next(iter)) {
        gtk_widget_destroy(GTK_WIDGET(iter->data));
        iter->data = NULL;
    }

    g_list_free(children);
}

void entry_clear(GtkWidget *entry)
{
    if (gtk_entry_get_text(GTK_ENTRY(entry))) {
        gtk_entry_set_text(GTK_ENTRY(entry), "");
    }
}

void label_clear(GtkWidget *label)
{
    if (gtk_label_get_text(GTK_LABEL(label))) {
        gtk_label_set_text(GTK_LABEL(label), NULL);
    }
}

void entry_set_gdouble_format(GtkWidget *entry, const gdouble num, const gchar *format)
{
    gchar str[256];

    g_snprintf(str, 256, format, num);
    gtk_entry_set_text(GTK_ENTRY(entry), str);
}

void entry_set_gint_format(GtkWidget *entry, const gint num, const gchar *format)
{
    gchar str[256];

    g_snprintf(str, 256, format, num);
    gtk_entry_set_text(GTK_ENTRY(entry), str);
}

void entry_set_text_format(GtkWidget *entry, const gchar *text, const gchar *format)
{
    gchar str[256];

    g_snprintf(str, 256, format, text);
    gtk_entry_set_text(GTK_ENTRY(entry), str);
}

void label_set_gdouble_format(GtkWidget *label, const gdouble num, const gchar *format)
{
    gchar str[256];

    g_snprintf(str, 256, format, num);
    gtk_label_set_text(GTK_LABEL(label), str);
}

void label_set_gint_format(GtkWidget *label, const gint num, const gchar *format)
{
    gchar str[256];

    g_snprintf(str, 256, format, num);
    gtk_label_set_text(GTK_LABEL(label), str);
}

void label_set_text_format(GtkWidget *label, const gchar *text, const gchar *format)
{
    gchar str[256];

    g_snprintf(str, 256, format, text);
    gtk_label_set_text(GTK_LABEL(label), str);
}

/* TODO: message windows should have parents */
void show_error(const gchar *msg)
{
    GtkWidget *error;

    error = gtk_message_dialog_new(NULL, GTK_DIALOG_DESTROY_WITH_PARENT, GTK_MESSAGE_ERROR, GTK_BUTTONS_CLOSE, "%s", msg);
    gtk_dialog_run(GTK_DIALOG(error));
    gtk_widget_destroy(error);
}

void show_save_error(FileSaveStatus filesavestatus, const gchar *filename)
{
    gchar buffer[500];

    switch (filesavestatus) {
    case SAVE_FILE_NO_DATA:
        show_error("No data to save.");
        break;

    case SAVE_FILE_COULD_NOT_OPEN:
        g_snprintf(buffer, sizeof(buffer), "Could not open file %s for saving.", filename);
        show_error(buffer);
        break;

    default:
	break;
    }
}

/* TODO: message windows should have parents */
void show_warning(const gchar *msg)
{
    GtkWidget *warning;

    warning = gtk_message_dialog_new(NULL, GTK_DIALOG_DESTROY_WITH_PARENT, GTK_MESSAGE_WARNING, GTK_BUTTONS_OK, "%s", msg);
    gtk_dialog_run(GTK_DIALOG(warning));
    gtk_widget_destroy(warning);
}

#define CURVE_COLOR_HISTO gwy_graph_get_preset_color(3)
void show_histogram(const gchar *quantity, const gchar *unit, const gdouble *xdata, const gdouble *ydata, gint ndata, GtkWidget *parent)
{
    GtkWidget *window, *graph;
    GwyGraphModel *graph_model;
    GwyGraphCurveModel *cmodel;

    gchar str[256];

    g_snprintf(str, sizeof(str), "show %s histogram", quantity);

    if (verbose) {
        g_print("Show histogram: %s\n", str);
    }

    window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(parent));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(window), TRUE);

    g_snprintf(str, sizeof(str), "Probability distribution of %s", quantity);
    gtk_window_set_title(GTK_WINDOW(window), str);
    gtk_window_set_default_size(GTK_WINDOW(window), 300, 200);
    g_signal_connect(window, "delete_event", G_CALLBACK(gtk_widget_destroy), NULL);   // CHECK

    graph_model = gwy_graph_model_new();
    g_snprintf(str, sizeof(str), "PDF %s", quantity);
    g_object_set(graph_model, "title", str, NULL);

    g_snprintf(str, sizeof(str), "%s / %s", quantity, unit);
    gwy_graph_model_set_axis_label(graph_model, GTK_POS_BOTTOM, str);
    gwy_graph_model_set_axis_label(graph_model, GTK_POS_LEFT, "count / a.u.");

    graph = gwy_graph_new(graph_model);
    gtk_widget_set_size_request(GTK_WIDGET(graph), 200, 150);
    gwy_graph_enable_user_input(GWY_GRAPH(graph), FALSE);
    gwy_graph_set_status(GWY_GRAPH(graph), GWY_GRAPH_STATUS_PLAIN);

    cmodel = gwy_graph_curve_model_new();
    g_snprintf(str, sizeof(str), "PDF(%s)", quantity);
    g_object_set(cmodel,
                 "mode", GWY_GRAPH_CURVE_LINE_POINTS,
                 "description", str,
                 "point-size", 3,
                 "point-type", GWY_GRAPH_POINT_CIRCLE,
                 "color", CURVE_COLOR_HISTO,
                 "line-style", GDK_LINE_SOLID,
                 NULL);
    gwy_graph_model_add_curve(graph_model, cmodel);
    gwy_graph_curve_model_set_data(cmodel, xdata, ydata, ndata);

    gtk_container_add(GTK_CONTAINER(window), graph);
    gtk_widget_show_all(window);
}

void show_log(GtkWidget *log_window, const gchar *logfnm, const gchar *logext, const gchar *window_title)
{
    GtkWidget *scrolled_window;
    GtkWidget *textview;
    GtkTextBuffer *buffer;
    gchar *log;
    gchar errmsg[500];
    gint len;
    FILE *logfile;
    gchar *logfnm_full;

    /*open log */
    logfnm_full = g_malloc(sizeof(gchar) * (strlen(logfnm) + 5));
    strcpy(logfnm_full, logfnm);
    strcat(logfnm_full, logext);
    logfile = fopen(logfnm_full, "r");

    if (!logfile) { /* should not get here */
	g_snprintf(errmsg, sizeof(errmsg), "Could not open log file '%s'.", logfnm_full);
	show_error(errmsg);
	g_free(logfnm_full);
	return;
    }
    else {
	/*find length and read file */
	fseek(logfile, 0, SEEK_END);
	len = ftell(logfile);
	log = (gchar *) g_malloc((len + 1) * sizeof(gchar));
	rewind(logfile);
	fread(log, len, 1, logfile);
	fclose(logfile);
	log[len] = '\0';
	g_free(logfnm_full);

	/*create buffer, fill it and create textview */
	buffer = gtk_text_buffer_new(NULL);
	gtk_text_buffer_set_text(buffer, log, -1);
	textview = gtk_text_view_new_with_buffer(buffer);
	gtk_widget_modify_font(textview, pango_font_description_from_string("monospace"));
	gtk_text_view_set_editable(GTK_TEXT_VIEW(textview), FALSE);
	g_free(log);

	scrolled_window = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_window), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
	gtk_container_add(GTK_CONTAINER(scrolled_window), textview);

	if (log_window != NULL && GTK_IS_WIDGET(log_window)) {
	    gtk_widget_destroy(log_window);
	}

	log_window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(log_window), window_title);
	gtk_widget_set_size_request(GTK_WIDGET(log_window), 450, 400);
	g_signal_connect(log_window, "delete_event", G_CALLBACK(gtk_widget_destroy), NULL);
	gtk_container_add(GTK_CONTAINER(log_window), scrolled_window);
	gtk_widget_show_all(log_window);
    }
}

void buffer_to_file_dialog(const gchar *buffer)
{
    GtkWidget *dialog;
    gchar *filename;
    gboolean file_ok;
    gint response;
    FileSaveStatus filestatus;

    if (!buffer) {
	show_error("No data to save.");
	return;
    }

    dialog = gtk_file_chooser_dialog_new("Save File",
                                         NULL,
                                         GTK_FILE_CHOOSER_ACTION_SAVE,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_SAVE, GTK_RESPONSE_OK,
                                         NULL);

    gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), gwy_app_get_current_directory());

    file_ok = FALSE;
    response = gtk_dialog_run(GTK_DIALOG(dialog));

    while (response == GTK_RESPONSE_OK && !file_ok) {
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        if (gwy_app_file_confirm_overwrite(dialog)) {
	    buffer_to_file(filename, buffer, &filestatus);
            file_ok = (filestatus == SAVE_FILE_OK);
        }

        g_free(filename);

        if (!file_ok) {
            show_save_error(filestatus, filename);
            response = gtk_dialog_run(GTK_DIALOG(dialog));
        }
    }

    gtk_widget_destroy(dialog);
}

#define DSWAP(x, y) GWY_SWAP(gdouble, x, y)

void rs_zoom_in(GwyGraph *graph, GSList **stack, gboolean free_y_rescale)
{
    GwyGraphModel *newgmodel, *oldgmodel;
    GwyGraphCurveModel *newcmodel, *oldcmodel;
    GwySelection *selection;
    GtkAllocation widget_alloc;

    gint ci, i, istart, iend, ndata, maxndata, oldndata, pointsize;
    gdouble ymin, ymax, avg, rng;
    const gdouble *oldxdata, *oldydata;
    gdouble *newxdata, *newydata;
    gdouble range[2];
    gboolean descending;

    oldgmodel = gwy_graph_get_model(graph);
    newgmodel = gwy_graph_model_new_alike(oldgmodel);

    selection = gwy_graph_area_get_selection(GWY_GRAPH_AREA(gwy_graph_get_area(GWY_GRAPH(graph))), GWY_GRAPH_STATUS_PLAIN);
    /* here PLAIN stands for the current selection mode */

    if (gwy_selection_get_object(selection, 0, range)) {

        if (range[0] > range[1]) {
            DSWAP(range[0], range[1]);
        }

        maxndata = 0;

        for (ci = 0; ci < gwy_graph_model_get_n_curves(oldgmodel); ci++) {

            oldcmodel = gwy_graph_model_get_curve(oldgmodel, ci);
            newcmodel = gwy_graph_curve_model_new_alike(oldcmodel);  /* LEAK */

            oldxdata = gwy_graph_curve_model_get_xdata(oldcmodel);
            oldydata = gwy_graph_curve_model_get_ydata(oldcmodel);
            oldndata = gwy_graph_curve_model_get_ndata(oldcmodel);

            descending = FALSE;

            if (oldndata > 0)
                if (oldxdata[0] > oldxdata[oldndata - 1]) {
                    descending = TRUE;
                }

            range_to_indices_array(range[0], range[1], oldxdata, oldndata, descending, &istart, &iend, &ndata);

            newxdata = (gdouble *)g_malloc(ndata * sizeof(gdouble));
            newydata = (gdouble *)g_malloc(ndata * sizeof(gdouble));

            for (i = 0; i < ndata; i++) {
                newxdata[i] = oldxdata[istart + i];
                newydata[i] = oldydata[istart + i];
            }

            if (free_y_rescale) {

                ymin = DBL_MAX;
                ymax = -DBL_MAX;

                for (i = 0; i < ndata; i++) {
                    if (newydata[i] < ymin) {
                        ymin = newydata[i];
                    }

                    if (newydata[i] > ymax) {
                        ymax = newydata[i];
                    }
                }

                avg = (ymin + ymax) / 2;
                rng = ymax - ymin;

                if (ymin != ymax) { /* rescale only nonconstant curves */
                    /* if constant, both ymin and ymax would be assigned from the same value, so == is correct */
                    for (i = 0; i < ndata; i++) {
                        newydata[i] = (newydata[i] - avg) / rng;
                    }
                }
                else {
                    ndata = 0; /* do not transfer data to the new model */
                }
            }

            if (ndata > 0) {
                gwy_graph_curve_model_set_data(newcmodel, newxdata, newydata, ndata);
                gwy_graph_model_add_curve(newgmodel, newcmodel);

                if (ndata > maxndata) {
                    maxndata = ndata;
                }
            }

            g_free(newxdata);
            g_free(newydata);
        }

        /* set point size, roughly based on the number of points in new curves */
        if (maxndata > 0) {
            gtk_widget_get_allocation(gwy_graph_get_area(GWY_GRAPH(graph)), &widget_alloc);
            pointsize = widget_alloc.width / (3 * maxndata);
            pointsize = CLAMP(pointsize, 0, 8);

            if (pointsize >= 3) {
                for (ci = 0; ci < gwy_graph_model_get_n_curves(newgmodel); ci++) {
                    newcmodel = gwy_graph_model_get_curve(newgmodel, ci);
                    g_object_set(newcmodel, "point-size", pointsize, NULL);
                }
            }
        }

        gwy_graph_set_model(graph, newgmodel);
        *stack = g_slist_prepend(*stack, newgmodel);
        gwy_selection_clear(selection);
    }
}

void rs_zoom_out(GwyGraph *graph, GSList **stack)
{
    if (g_slist_length(*stack) >= 2) {
        gwy_graph_set_model(graph, (GwyGraphModel *)g_slist_nth_data(*stack, 1));
        g_object_unref(g_slist_nth_data(*stack, 0));
        *stack = g_slist_remove(*stack, g_slist_nth_data(*stack, 0));
    }
}

void rs_zoom_restore(GwyGraph *graph, GSList **stack)
{
    guint i, n;

    n = g_slist_length(*stack);

    if (n >= 2) {
        for (i = 0; i < n - 1; i++) {
            g_object_unref(g_slist_nth_data(*stack, 0));
            *stack = g_slist_remove(*stack, g_slist_nth_data(*stack, 0));
        }

        gwy_graph_set_model(graph, (GwyGraphModel *)g_slist_nth_data(*stack, 0));
    }
}

gboolean beta_changed(GtkWidget *wbeta, gdouble *beta)
{
    const gchar *sbeta;
    gchar buffer[24];
    gdouble beta_prev;

    beta_prev = *beta;
    sbeta = gtk_entry_get_text(GTK_ENTRY(wbeta));

    if (strlen(sbeta)) {
        *beta = g_strtod(sbeta, NULL);

        /* if beta reasonable */
        if (*beta > 0) {
            g_snprintf(buffer, sizeof(buffer), "%.03f", *beta);
            gtk_entry_set_text(GTK_ENTRY(wbeta), buffer);
            return TRUE;
        }

        /* otherwise keep old beta */
        *beta = beta_prev;
    }

    g_snprintf(buffer, sizeof(buffer), "%.03f", *beta);
    gtk_entry_set_text(GTK_ENTRY(wbeta), buffer);
    return FALSE;
}

/* tyhle dve from,to mozna spojit do jedne */
void range_from_changed(GtkWidget *wfrom, GtkWidget *wto, GwySelection *selection, gdouble *datafrom, gdouble *data2from)
{
    gdouble from, to;
    gdouble range[2];
    gchar buffer[24];
    const gchar *sfrom, *sto;

    sfrom = gtk_entry_get_text(GTK_ENTRY(wfrom));
    sto   = gtk_entry_get_text(GTK_ENTRY(wto));

    if (strlen(sfrom)) {
        from = g_strtod(sfrom, NULL);

        if (strlen(sto)) {
            to = g_strtod(sto, NULL);

            if (from > to) {
                to = from;
            }

            range[0] = from;
            range[1] = to;
            gwy_selection_set_object(selection, 0, range);
        }

        *datafrom = from;
        *data2from = from;
        g_snprintf(buffer, sizeof(buffer), DEPTH, *datafrom);
        gtk_entry_set_text(GTK_ENTRY(wfrom), buffer);
    }
}

void range_to_changed(GtkWidget *wfrom, GtkWidget *wto, GwySelection *selection, gdouble *datato, gdouble *data2to)
{
    gdouble from, to;
    gdouble range[2];
    gchar buffer[24];
    const gchar *sfrom, *sto;

    sfrom = gtk_entry_get_text(GTK_ENTRY(wfrom));
    sto   = gtk_entry_get_text(GTK_ENTRY(wto));

    if (strlen(sto)) {
        to  = g_strtod(sto, NULL);

        if (strlen(sfrom)) {
            from = g_strtod(sfrom, NULL);

            if (from > to) {
                from = to;
            }

            range[0] = from;
            range[1] = to;
            gwy_selection_set_object(selection, 0, range);
        }

        *datato = to;
        *data2to = to;
        g_snprintf(buffer, sizeof(buffer), DEPTH, *datato);
        gtk_entry_set_text(GTK_ENTRY(wto), buffer);
    }
}

/* tyhle dve from,to pct Fmax mozna spojit do jedne, pripadne s predchozimi */
gboolean range_from_pct_Fmax_changed(GtkWidget *from_pct_Fmax, GtkWidget *to_pct_Fmax, GwySelection *selection,
                                     const gdouble Fmax, GwyDataLine *hdataline, GwyDataLine *Fdataline,
                                     gdouble *datafrom, gdouble *data2from, gdouble *datato, gdouble *data2to,
                                     gboolean *Finput, gboolean *Finput2)
{
    gdouble from_pct, to_pct;
    gdouble range[2];
    gint istart, iend, ndata;
    gdouble from_F, to_F;
    gchar buffer[24];
    const gchar *sfrom, *sto;

    sfrom = gtk_entry_get_text(GTK_ENTRY(from_pct_Fmax));
    sto   = gtk_entry_get_text(GTK_ENTRY(to_pct_Fmax));

    if (strlen(sfrom)) {
        from_pct = g_strtod(sfrom, NULL);

        if (strlen(sto)) {
            to_pct  = g_strtod(sto, NULL);

            if (from_pct > to_pct) {
                to_pct = from_pct;
            }

            from_F = from_pct * 0.01 * Fmax;
            to_F = to_pct * 0.01 * Fmax;

            range_to_indices(from_F, to_F, Fdataline, TRUE, &istart, &iend, &ndata);

            range[0] = hdataline->data[iend];
            range[1] = hdataline->data[istart];

            gwy_selection_set_object(selection, 0, range);

            *datato = to_pct;
            *data2to = to_pct;
            g_snprintf(buffer, sizeof(buffer), DEPTH, *datato);
            gtk_entry_set_text(GTK_ENTRY(to_pct_Fmax), buffer);
            *Finput = TRUE;
            *Finput2 = TRUE;
        }

        *datafrom = from_pct;
        *data2from = from_pct;
        g_snprintf(buffer, sizeof(buffer), DEPTH, *datafrom);
        gtk_entry_set_text(GTK_ENTRY(from_pct_Fmax), buffer);
    }

    return FALSE;
}

gboolean range_to_pct_Fmax_changed(GtkWidget *from_pct_Fmax, GtkWidget *to_pct_Fmax, GwySelection *selection,
                                   const gdouble Fmax, GwyDataLine *hdataline, GwyDataLine *Fdataline,
                                   gdouble *datafrom, gdouble *data2from, gdouble *datato, gdouble *data2to,
                                   gboolean *Finput, gboolean *Finput2)
{
    gdouble from_pct, to_pct;
    gdouble range[2];
    gint istart, iend, ndata;
    gdouble from_F, to_F;
    gchar buffer[24];
    const gchar *sfrom, *sto;

    sfrom = gtk_entry_get_text(GTK_ENTRY(from_pct_Fmax));
    sto   = gtk_entry_get_text(GTK_ENTRY(to_pct_Fmax));

    if (strlen(sto)) {
        to_pct = g_strtod(sto, NULL);

        if (strlen(sfrom)) {
            from_pct = g_strtod(sfrom, NULL);

            if (from_pct > to_pct) {
                from_pct = to_pct;
            }

            from_F = from_pct * 0.01 * Fmax;
            to_F = to_pct * 0.01 * Fmax;

            range_to_indices(from_F, to_F, Fdataline, TRUE, &istart, &iend, &ndata);

            range[0] = hdataline->data[iend];
            range[1] = hdataline->data[istart];

            gwy_selection_set_object(selection, 0, range);

            *datafrom = from_pct;
            *data2from = from_pct;
            g_snprintf(buffer, sizeof(buffer), DEPTH, *datafrom);
            gtk_entry_set_text(GTK_ENTRY(from_pct_Fmax), buffer);
            *Finput = TRUE;
            *Finput2 = TRUE;
        }

        *datato = to_pct;
        *data2to = to_pct;
        g_snprintf(buffer, sizeof(buffer), DEPTH, *datato);
        gtk_entry_set_text(GTK_ENTRY(to_pct_Fmax), buffer);
    }

    return FALSE;
}

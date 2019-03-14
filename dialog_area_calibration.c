#include "dialog_area_calibration.h"

#include "controls.h" /* until AreaControls are defined locally */
#include "datatypes.h"
#include "gui-utils.h"
#include "file-utils.h"
#include "niget-common.h"

#include "niget-gtk.h" /* controls definitions */
#include "tool_contact.h" /* controls definitions */

#ifdef G_OS_WIN32
#include "getline.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h> /* memset */

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>
#include <libprocess/gwyprocess.h>


#define XAREAMAX 300
#define NPOLYDATA 100

static void area_open_file_dialog(Data *data, G_GNUC_UNUSED GtkWidget *parent);
static gboolean load_area_file(Data *data, const gchar *filename);
static gboolean read_coefficient_file(Data *data, const gchar *filename);
static void area_polynom_callback(GtkWidget *widget, Data *data);
static void area_data_callback(GtkWidget *widget, Data *data);
static void polynom_coeff_changed(Data *data);
static void area_open_coeff_file_dialog(Data *data, G_GNUC_UNUSED GtkWidget *parent);

static void area_open_file_dialog(Data *data, G_GNUC_UNUSED GtkWidget *parent)
{
    GtkWidget *dialog;
    gchar *filename, *basename;
    gboolean fileok;

    dialog = gtk_file_chooser_dialog_new("Open File",
                                         GTK_WINDOW(data->controls->areacontrols->dialog),
                                         GTK_FILE_CHOOSER_ACTION_OPEN,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                         NULL);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
        fileok = load_area_file(data, filename);

        if (fileok) {
            if (data->args->area.xdata->res > 0) {
                data->args->area.filename = g_strdup(filename);
                basename = g_path_get_basename(data->args->area.filename);
                gtk_entry_set_text(GTK_ENTRY(data->controls->areacontrols->area_entry), basename);
                g_free(basename);
            }
        }

        g_free(filename);
    }

    gtk_widget_destroy(dialog);
}

static gboolean load_area_file(Data *data, const gchar *filename)
{
    /*
    FILE *file;
    gint numdata, i, n;
    gchar *line = NULL;
    size_t size;
    gdouble x, y;
    */

    AreaFileReadStatus status;

    status = read_area_dat_file(filename, &data->args->area);

    switch (status) {
    case AREA_FILE_COULD_NOT_OPEN:
        show_error("Could not open area data file.");
        return FALSE;
        break;

    case AREA_FILE_NOT_ENOUGH_ENTRIES:
        show_error("File does not contain enough entries with specified data format.");
        return FALSE;
        break;

    case AREA_FILE_OK:
        gwy_graph_curve_model_set_data(data->controls->areacontrols->cmodel, data->args->area.xdata->data, data->args->area.ydata->data, data->args->area.xdata->res);
        return TRUE;
        break;

    default:
        show_error("Shouldn't get here.");
        return FALSE;
        break;
    }

    /*
    file = fopen(filename, "r");

    if (!file) {
        show_error("Could not open area calibration file.");
        return FALSE;
    }

    numdata = 0;
    while (getline(&line, &size, file) != EOF)
    	numdata++;
    rewind(file);

    gwy_data_line_resample(data->args->area.xdata, numdata+1, GWY_INTERPOLATION_NONE);
    gwy_data_line_resample(data->args->area.ydata, numdata+1, GWY_INTERPOLATION_NONE);

    data->args->area.xdata->data[0] = 0; // TODO je to nutne a rozumne?
    data->args->area.ydata->data[0] = 0;

    for (i = 0; i < numdata; i++) {
    	getline(&line, &size, file);
    	localize_decimal_points(line);
    	n = sscanf(line, "%lf %lf", &x, &y);
    	// printf("x %g y %g \n", x, y);
    	if (n < 2) {
    	    show_error("Need two columns of data.");
    	    return FALSE;
    	}
    	data->args->area.xdata->data[i+1] = x;
    	data->args->area.ydata->data[i+1] = y;
    }
    fclose(file);
    g_free(line);
    */
}

void area_calibration(GtkWidget *widget, Data *data)
{
    AreaControls *areacontrols;
    Area *area;

    gint response;
    gint i;
    gint row;
    gboolean fileok;

    GSList *group;
    gchar str[300];
    gchar *str2;
    gchar *basename;
    GtkWidget *label;
    gboolean nok = TRUE;

    areacontrols = data->controls->areacontrols;
    area = &(data->args->area);

    // poresit, aby se pri kazdem volani nealokovalo cele nove gui

    areacontrols->dialog = gtk_dialog_new_with_buttons("Area Function",
                           GTK_WINDOW(data->controls->maincontrols->window_main),
                           GTK_DIALOG_DESTROY_WITH_PARENT,
                           GTK_STOCK_APPLY, GTK_RESPONSE_OK,
                           GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                           NULL);

    areacontrols->hbox = gtk_hbox_new(FALSE, 2);
    gtk_container_add(GTK_CONTAINER(gtk_dialog_get_content_area(GTK_DIALOG(areacontrols->dialog))), areacontrols->hbox);

    areacontrols->area_table = gtk_table_new(10, 4, FALSE);

    gtk_table_set_row_spacings(GTK_TABLE(areacontrols->area_table), CTRLS_TABLE_ROW_SPACING);
    gtk_table_set_col_spacings(GTK_TABLE(areacontrols->area_table), CTRLS_TABLE_COL_SPACING);
    gtk_box_pack_start(GTK_BOX(areacontrols->hbox), areacontrols->area_table, FALSE, FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(areacontrols->area_table), CTRLS_BORDER);

    /* OP polynomial button */
    areacontrols->opp_button = gtk_radio_button_new_with_label(NULL, "OP polynomial");
    gtk_table_attach(GTK_TABLE(areacontrols->area_table), areacontrols->opp_button, 1, 2, 0, 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(G_OBJECT(areacontrols->opp_button), "toggled", G_CALLBACK(area_polynom_callback), data);

    if (area->mode == AREA_DATA) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(areacontrols->opp_button), FALSE);
    }
    else {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(areacontrols->opp_button), TRUE);
    }

    /* OP Polynomial */
    areacontrols->area_polynom = (GtkWidget **)g_malloc(area->npoly * sizeof(GtkWidget *));
    row = 1;

    for (i = 0; i < area->npoly; i ++) {
        g_snprintf(str, sizeof(str), "x<sup>%s</sup>", area->polylabel[i]);
        label = label_new_with_markup_left(str);
        gtk_table_attach(GTK_TABLE(areacontrols->area_table), label, 2, 3, row, row + 1, GTK_FILL, 0, 0, 0);
        areacontrols->area_polynom[i] = gtk_entry_new();
        gtk_entry_set_width_chars(GTK_ENTRY(areacontrols->area_polynom[i]), 12);
        g_snprintf(str, sizeof(str), "%g", area->polycoef[i]);
        gtk_entry_set_text(GTK_ENTRY(areacontrols->area_polynom[i]), str);
        gtk_table_attach(GTK_TABLE(areacontrols->area_table), areacontrols->area_polynom[i], 1, 2, row, row + 1, GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(areacontrols->area_polynom[i], "activate", G_CALLBACK(polynom_coeff_changed), data);
        g_signal_connect_swapped(areacontrols->area_polynom[i], "focus-out-event", G_CALLBACK(polynom_coeff_changed), data);

        if (area->mode == AREA_DATA) {
            gwy_table_hscale_set_sensitive(GTK_OBJECT(areacontrols->area_polynom[i]), FALSE);
        }
        else {
            gwy_table_hscale_set_sensitive(GTK_OBJECT(areacontrols->area_polynom[i]), TRUE);
        }

        row++;
    }

    g_snprintf(str, sizeof(str), "x<sup>%g</sup>", area->polypower[0]);

    /* Polynomial file name entry */
    areacontrols->coeff_entry = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(areacontrols->coeff_entry), 18);

    gtk_table_attach(GTK_TABLE(areacontrols->area_table), areacontrols->coeff_entry, 1, 2, row, row + 1, GTK_FILL, 0, 0, 0);

    if (area->coeff_filename) {
        if (area->mode == AREA_POLYNOM) {
            fileok = read_coefficient_file(data, area->coeff_filename);

            if (fileok) {
                basename = g_path_get_basename(area->coeff_filename);
                gtk_entry_set_text(GTK_ENTRY(areacontrols->coeff_entry), basename);
                g_free(basename);

                for (i = 0; i < area->npoly; i ++) {
                    g_snprintf(str, sizeof(str), "%g", area->polycoef[i]);
                    gtk_entry_set_text(GTK_ENTRY(areacontrols->area_polynom[i]), str);
                }
            }
        }
    }

    gtk_widget_set_sensitive(areacontrols->coeff_entry, FALSE);
    row++;

    /* Polynomial open file button */
    areacontrols->coeff_file_button = gtk_button_new_with_mnemonic("Open _coeff. file");
    gtk_table_attach(GTK_TABLE(areacontrols->area_table), areacontrols->coeff_file_button, 1, 2, row, row + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(areacontrols->coeff_file_button, "clicked", G_CALLBACK(area_open_coeff_file_dialog), data);
    row++;

    /* Area entry */
    areacontrols->area_entry = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(areacontrols->area_entry), 18);

    gtk_table_attach(GTK_TABLE(areacontrols->area_table), areacontrols->area_entry, 3, 4, 1, 2, GTK_FILL, 0, 0, 0);

    if (area->filename) {
        if (area->mode == AREA_DATA) {
            fileok = load_area_file(data, area->filename);

            if (fileok) {
                basename = g_path_get_basename(area->filename);
                gtk_entry_set_text(GTK_ENTRY(areacontrols->area_entry), basename);
                g_free(basename);
            }
        }
    }

    gtk_widget_set_sensitive(areacontrols->area_entry, FALSE);

    /* Area file button */
    areacontrols->area_file_button = gtk_button_new_from_stock(GTK_STOCK_OPEN);
    gtk_table_attach(GTK_TABLE(areacontrols->area_table), areacontrols->area_file_button, 3, 4, 2, 3, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(areacontrols->area_file_button, "clicked", G_CALLBACK(area_open_file_dialog), data);

    if (area->mode == AREA_DATA) {
        fileok = load_area_file(data, area->filename);

        if (fileok) {
            basename = g_path_get_basename(area->filename);
            gtk_entry_set_text(GTK_ENTRY(areacontrols->area_entry), basename);
            g_free(basename);
            gtk_widget_set_sensitive(areacontrols->area_file_button, TRUE);
        }
    }
    else {
        gtk_widget_set_sensitive(areacontrols->area_file_button, FALSE);
    }

    /* Raw data button */
    //must be created after area_file_button and area_entry because its callback needs them?
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(areacontrols->opp_button));
    areacontrols->rd_button = gtk_radio_button_new_with_label(group, "Raw data");
    gtk_table_attach(GTK_TABLE(areacontrols->area_table), areacontrols->rd_button, 3, 4, 0, 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(G_OBJECT(areacontrols->rd_button), "toggled", G_CALLBACK(area_data_callback), data);

    if (area->mode == AREA_DATA) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(areacontrols->rd_button), TRUE);
    }
    else {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(areacontrols->rd_button), FALSE);
    }

    /*Graph */
    areacontrols->gmodel = gwy_graph_model_new();
    g_object_set(areacontrols->gmodel, "title", "Area Function", NULL);

    areacontrols->area_graph = gwy_graph_new(areacontrols->gmodel);
    gtk_widget_set_size_request(areacontrols->area_graph, 400, 300);

    gwy_graph_enable_user_input(GWY_GRAPH(areacontrols->area_graph), FALSE);
    gtk_box_pack_start(GTK_BOX(areacontrols->hbox), areacontrols->area_graph, TRUE, TRUE, 0);
    gwy_graph_set_status(GWY_GRAPH(areacontrols->area_graph), GWY_GRAPH_STATUS_PLAIN);

    areacontrols->cmodel = gwy_graph_curve_model_new();

    if (area->mode != AREA_DATA) {
        area->xdata = gwy_data_line_new(NPOLYDATA, NPOLYDATA, FALSE);
        area->ydata = gwy_data_line_new(NPOLYDATA, NPOLYDATA, FALSE);

        for (i = 0; i < NPOLYDATA; i++) {
            area->xdata->data[i] = XAREAMAX / NPOLYDATA * i;
            area->ydata->data[i] = polynom_eval(area, area->xdata->data[i]);
        }

        gwy_graph_curve_model_set_data(areacontrols->cmodel, area->xdata->data, area->ydata->data, area->xdata->res);
    }
    else {
        gwy_graph_curve_model_set_data(areacontrols->cmodel, area->xdata->data, area->ydata->data, area->xdata->res);
    }

    g_object_set(areacontrols->cmodel,
                 "mode", GWY_GRAPH_CURVE_LINE,
                 "color", gwy_graph_get_preset_color(1),
                 "line-style", GDK_LINE_SOLID,
                 NULL);
    gwy_graph_model_add_curve(areacontrols->gmodel, areacontrols->cmodel);
    gwy_graph_model_set_axis_label(areacontrols->gmodel, GTK_POS_BOTTOM, "h / nm");
    gwy_graph_model_set_axis_label(areacontrols->gmodel, GTK_POS_LEFT, "A / nm<sup>2</sup>");

    gtk_widget_show_all(areacontrols->dialog);

    response = gtk_dialog_run(GTK_DIALOG(areacontrols->dialog));

    while (nok) {
        nok = FALSE;

        switch (response) {
        case GTK_RESPONSE_CANCEL:
        case GTK_RESPONSE_DELETE_EVENT:
            break;

        case GTK_RESPONSE_NONE:
            return;
            break;

        case GTK_RESPONSE_OK:

            // pass
            // check if in case AREA_DATA was chosen a filename was provided

            if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(areacontrols->rd_button)) && (area->filename == NULL)) {
                area->mode = AREA_TRIVIAL;
                show_error("For the option 'Raw Data' an area calibration data file must be provided.");
                response = gtk_dialog_run(GTK_DIALOG(areacontrols->dialog));
                break;
            }

            if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(areacontrols->rd_button))) {
                area->mode = AREA_DATA;
            }
            else {
                area->mode = AREA_POLYNOM;

                if (area->filename != NULL) {
                    g_free(area->filename);
                    area->filename = NULL;
                }
            }

            str2 = write_area_label(area);
            gtk_label_set_markup(GTK_LABEL(data->controls->contactcontrols->area_function_label), str2);
            g_free(str2);

            data->functions->remove_all_tool_fit_results_labels(data, TRUE);
            break;

        default:
            g_warning("Unhandled dialog response");
            return;
            break;
        }
    }

    gtk_widget_destroy(areacontrols->dialog);
}

static void area_polynom_callback(GtkWidget *widget, Data *data)
{
    gint i;

    if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget))) {
        return;
    }

    gtk_widget_set_sensitive(data->controls->areacontrols->area_file_button, FALSE);
    gtk_widget_set_sensitive(data->controls->areacontrols->coeff_file_button, TRUE);

    for (i = 0; i < data->args->area.npoly; i++) {
        gwy_table_hscale_set_sensitive(GTK_OBJECT(data->controls->areacontrols->area_polynom[i]), TRUE);
    }

    gwy_data_line_resample(data->args->area.xdata, NPOLYDATA, GWY_INTERPOLATION_NONE);
    gwy_data_line_resample(data->args->area.ydata, NPOLYDATA, GWY_INTERPOLATION_NONE);

    for (i = 0; i < NPOLYDATA; i++) {
        data->args->area.xdata->data[i] = XAREAMAX / NPOLYDATA * i;
        data->args->area.ydata->data[i] = polynom_eval(&(data->args->area), data->args->area.xdata->data[i]);
    }

    gwy_graph_curve_model_set_data(data->controls->areacontrols->cmodel, data->args->area.xdata->data, data->args->area.ydata->data, data->args->area.xdata->res);
}

static void area_data_callback(GtkWidget *widget, Data *data)
{
    gint i;
    gboolean fileok = FALSE;

    if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget))) {
        return;
    }

    gtk_widget_set_sensitive(data->controls->areacontrols->area_file_button, TRUE);
    gtk_widget_set_sensitive(data->controls->areacontrols->coeff_file_button, FALSE);

    for (i = 0; i < data->args->area.npoly; i++) {
        gwy_table_hscale_set_sensitive(GTK_OBJECT(data->controls->areacontrols->area_polynom[i]), FALSE);
    }

    if (data->args->area.filename != NULL) {
        fileok = load_area_file(data, data->args->area.filename);
    }

    if (fileok) {
        gwy_graph_curve_model_set_data(data->controls->areacontrols->cmodel, data->args->area.xdata->data, data->args->area.ydata->data, data->args->area.xdata->res);
    }
    else {
        gwy_graph_curve_model_set_data(data->controls->areacontrols->cmodel, NULL, NULL, 0);
    }
}

static void polynom_coeff_changed(Data *data)
{
    gint i;
    const gchar *str;

    for (i = 0; i < data->args->area.npoly; i++) {
        str = gtk_entry_get_text(GTK_ENTRY(data->controls->areacontrols->area_polynom[i]));
        data->args->area.polycoef[i] = g_strtod(str, NULL);
    }

    gwy_data_line_resample(data->args->area.xdata, NPOLYDATA, GWY_INTERPOLATION_NONE);
    gwy_data_line_resample(data->args->area.ydata, NPOLYDATA, GWY_INTERPOLATION_NONE);

    for (i = 0; i < NPOLYDATA; i++) {
        data->args->area.xdata->data[i] = XAREAMAX / NPOLYDATA * i;
        data->args->area.ydata->data[i] = polynom_eval(&(data->args->area), data->args->area.xdata->data[i]);
    }

    gwy_graph_curve_model_set_data(data->controls->areacontrols->cmodel, data->args->area.xdata->data, data->args->area.ydata->data, NPOLYDATA);
}

static void area_open_coeff_file_dialog(Data *data, G_GNUC_UNUSED GtkWidget *parent)
{
    GtkWidget       *dialog;
    GtkFileFilter   *filter;
    gchar           *filename;
    gchar           *basename;
    gboolean        fileok;
    gchar str[300];
    gint i;
    AreaControls *areacontrols;
    Area *area;

    areacontrols = data->controls->areacontrols;
    area = &(data->args->area);

    dialog = gtk_file_chooser_dialog_new("Open File",
                                         GTK_WINDOW(areacontrols->dialog),
                                         GTK_FILE_CHOOSER_ACTION_OPEN,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                         NULL);

    filter = gtk_file_filter_new();
    gtk_file_filter_set_name(filter, "Coefficient file (*.ara)");
    gtk_file_filter_add_pattern(filter, "*.ara");
    gtk_file_filter_set_name(filter, "Coefficient file (*.ind)");
    gtk_file_filter_add_pattern(filter, "*.ind");
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filter);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        basename = g_path_get_basename(filename);
        fileok = read_coefficient_file(data, filename);

        if (fileok) {
            area->coeff_filename = g_strdup(basename);
            gtk_entry_set_text(GTK_ENTRY(areacontrols->coeff_entry), basename);

            for (i = 0; i < area->npoly; i ++) {
                g_snprintf(str, sizeof(str), "%g", area->polycoef[i]);
                gtk_entry_set_text(GTK_ENTRY(areacontrols->area_polynom[i]), str);
            }

            polynom_coeff_changed(data);
        }

        g_free(filename);
        g_free(basename);
    }

    gtk_widget_destroy(dialog);
}


static gboolean read_coefficient_file(Data *data, const gchar *filename)
{
    AreaFileReadStatus filestatus;
    gchar **splitlist;
    gint i;

    filestatus = AREA_FILE_COULD_NOT_OPEN;

    splitlist = g_strsplit(filename, ".", 0);
    i = 0;
    
    while (splitlist[i]) {
	i++;
    }
    
    i--;
    
    if (i > 0) { /* filename is not NULL and contains a suffix */
        if (0 == g_ascii_strncasecmp(splitlist[i], "ara", 3)) {
            filestatus = read_area_hysitron_file(filename, &data->args->area);
        }
        else if (0 == g_ascii_strncasecmp(splitlist[i], "ind", 3)) {
            filestatus = read_area_csm_file(filename, &data->args->area);
        }
        else if (0 == g_ascii_strncasecmp(splitlist[i], "ani", 3)) {
            filestatus = read_area_niget_file(filename, &data->args->area);
        }
    }

    g_strfreev(splitlist);

    switch (filestatus) {
    case AREA_FILE_COULD_NOT_OPEN:
        show_error("Could not open file.");
        return FALSE;
        break;

    case AREA_FILE_NOT_ENOUGH_ENTRIES:
        show_error("File does not contain enough entries with specified data format.");
        return FALSE;
        break;

    case AREA_FILE_FORMAT_MISMATCH:
        show_error("File does not have the correct format.");
        return FALSE;
        break;

    case AREA_FILE_POLYNOM:
        show_error("Currently only Oliver Pharr polynoms up to order 5 are supported.");
        return FALSE;
        break;

    case AREA_FILE_TIP:
        show_error("Currently only Berkovich tips with coefficient 24.5 are implemented.");
        return FALSE;
        break;

    case AREA_FILE_OK:
	
	if (verbose) {
	    for (i = 0; i < data->args->area.npoly; i++) {
		printf(" %d %g \n", i, data->args->area.polycoef[i]);
	    }
	}

        return TRUE;
        break;

    default:
        show_error("Shouldn't get here.");
        return FALSE;
        break;
    }
}

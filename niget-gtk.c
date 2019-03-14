#include "niget-gtk.h"

#include "controls.h"
#include "datatypes.h"
#include "fddata.h"
#include "gui-utils.h"
#include "niget-common.h"
#include "settings.h"

/* #ifndef G_OS_WIN32 */
/* #include "resources.h" */
/* #endif */

#include "dialog_about.h"
#include "dialog_area_calibration.h" /* controls definitions */

#include "tool_contact.h"

#include "tool_op.h"
#include "tool_op_gui.h"
#include "tool_op_unc.h"
#include "tool_op_unc_gui.h" /* controls definitions */

#ifndef NOFORTRAN

#include "tool_op_odr.h"
#include "tool_op_odr_gui.h"
#include "tool_op_odr_unc.h"
#include "tool_op_odr_unc_gui.h" /* controls definitions */

#include "tool_twoslopes.h"
#include "tool_twoslopes_gui.h"
#include "tool_twoslopes_unc.h"
#include "tool_twoslopes_unc_gui.h" /* controls definitions */

#include "tool_hertz_odr.h"
#include "tool_hertz_odr_gui.h"
#include "tool_hertz_odr_unc.h"
#include "tool_hertz_odr_unc_gui.h" /* controls definitions */

#include "tool_stiffness.h"
#include "tool_stiffness_gui.h"
#include "tool_stiffness_unc.h"
#include "tool_stiffness_unc_gui.h" /* controls definitions */

#endif

#include "tool_tangent.h"
#include "tool_tangent_gui.h"
#include "tool_tangent_unc.h"
#include "tool_tangent_unc_gui.h" /* controls definitions */

#include "tool_hertz.h"
#include "tool_hertz_gui.h"
#include "tool_hertz_unc.h"
#include "tool_hertz_unc_gui.h" /* controls definitions */

#include "tool_Ph2.h"
#include "tool_Ph2_gui.h"

#include "tool_apopins.h"
#include "tool_apopins_gui.h"

#include "tool_epwork.h"
#include "tool_epwork_gui.h"
#include "tool_epwork_unc.h"
#include "tool_epwork_unc_gui.h" /* controls definitions */


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gtk/gtk.h>

#include <libgwydgets/gwydgets.h>

#if defined(G_OS_WIN32) && (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 54)
/* glib 2.54 proved working, older versions not available to check */
#include <windows.h>
#endif

/* place somewhere else if possible */
/* keep items consitent with the enum DataFileFormat */
static const char *data_file_format_names[] = {"Niget: Preprocessed and split data",
                                               "UNHT: Time / Depth (nm) / Load (mN)",
                                               "UNHT: Time / Depth (nm) / Load (uN)",
                                               "Hysitron: Depth (nm) / Load (uN) / ...",
                                               "Other: Depth (nm) / Load (mN) / ...",
                                               "Other: Load (mN) / Depth (nm) / ...",
                                               "Other: Load (uN) / Depth (nm) / ..."
                                              };

/* zmenit na "niget.conf" nebo tak neco */
static const gchar *settings_filename = "niget_settings.cfg";

extern const gchar *point_sel_status[4];

enum DataFileFormat lastfileformat;

static void init_args(Args *args);
static void load_data(GtkFileChooserButton *widget, Data *data);
static void load_data_gui(const gchar *filename, enum DataFileFormat fileformat, Data *data);
static gboolean open_file_dialog(gchar **filename, enum DataFileFormat *fileformat, GtkWidget *parent);
static gboolean file_confirm_overwrite(gchar *fnm);

static void save_tool_export_data(gchar *basedir, gchar *basefile, const gchar *suffix, gchar *buffer);

static void save_data(Data *data);
static void show_help(void);
static void notebook_set_sensitive_all_pages(MainControls *maincontrols, gboolean val);
static void restore_zoom_all_tool_graphs(Data *data);
static void remove_all_tool_fit_results_labels(Data *data, gboolean areadependent);
static void redraw_all_tools(Data *data);

static void init_args(Args *args)
{
    init_FDdata(&(args->fddata));
    init_Area(&(args->area));
    init_Instdata(&(args->instdata));
    init_OPdata(&(args->opdata));
    init_OPUncdata(&(args->opunc));
#ifndef NOFORTRAN
    init_OPODRdata(&(args->opodrdata));
    init_OPODRUncdata(&(args->opodrunc), &(args->instdata));
    init_Slopesdata(&(args->slopesdata));
    init_SlopesUncdata(&(args->slopesunc), &(args->instdata));
    init_HertzODRdata(&(args->hertzodrdata));
    init_HertzODRUncdata(&(args->hertzodrunc), &(args->instdata));
    init_Stiffnessdata(&(args->stiffnessdata));
    init_StiffnessUncdata(&(args->stiffnessunc), &(args->instdata));
#endif
    init_Tangentdata(&(args->tgdata));
    init_TangentUncdata(&(args->tgunc));
    init_Hertzdata(&(args->hertzdata));
    init_HertzUncdata(&(args->hertzunc));
    init_Ph2data(&(args->ph2data));
    init_Apopindata(&(args->apopdata));
    init_Workdata(&(args->workdata));
    init_WorkUncdata(&(args->workunc));
}

static gchar *build_tool_filename(gchar *basefile, gchar *basedir, const gchar *prefix)
{
    gchar toolprefix[100];
    gchar *toolfilename;
    gchar *fnm;

    g_snprintf(toolprefix, sizeof(toolprefix), "%s", prefix);
    toolfilename = g_malloc(1 + strlen(basefile) + strlen(toolprefix));;
    strcpy(toolfilename, basefile);
    strcat(toolfilename, toolprefix);
    fnm = g_build_filename(basedir, toolfilename, NULL);

    g_free(toolfilename);
    return fnm;
}

static gboolean file_confirm_overwrite(gchar *fnm)
{
    gboolean write  = TRUE;
    GtkWidget *dialog2;

    if (g_file_test(fnm, G_FILE_TEST_EXISTS)) {
        dialog2 = gtk_message_dialog_new(NULL, GTK_DIALOG_MODAL,
                                         GTK_MESSAGE_ERROR, GTK_BUTTONS_YES_NO,
                                         _("File `%s' already exists. Replace?"),
                                         fnm);
        write = (gtk_dialog_run(GTK_DIALOG(dialog2)) == GTK_RESPONSE_YES);
        gtk_widget_destroy(dialog2);
    }

    return write;
}

static void save_tool_export_data(gchar *basedir, gchar *basefile, const gchar *suffix, gchar *buffer)
{
    gchar *fnm;
    FileSaveStatus filestatus;

    fnm = build_tool_filename(basefile, basedir, suffix);

    if (file_confirm_overwrite(fnm)) {
        buffer_to_file(fnm, buffer, &filestatus);

        if (filestatus != SAVE_FILE_OK) {
            show_save_error(filestatus, fnm);
        }
    }

    g_free(fnm);
}

static void save_data(Data *data)
{
    GtkWidget *dialog;
    gchar *suggestdir, *suggestname;
    gchar *basefilename;
    gchar *basedir, *basefile;
    gchar *buffer;
    Args *args;
    Controls *ctrls;

    args = data->args;
    ctrls = data->controls;

    if (args->fddata.ndata == 0) {
        show_error("No data loaded.");
        return;
    }

    dialog = gtk_file_chooser_dialog_new("Save File",
                                         NULL,
                                         GTK_FILE_CHOOSER_ACTION_SAVE,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_SAVE, GTK_RESPONSE_OK,
                                         NULL);

    suggestdir = g_path_get_dirname(data->args->fddata.filename);
    suggestname = g_path_get_basename(data->args->fddata.filename);

    gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), suggestdir); // TODO davat ke zdrojovym datum nebo do aktualniho adresare?
    gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(dialog), suggestname);

    // TODO zjistit nejdriv jestli vubec existuji nejaka data k ulozeni, nebo nejdriv jmeno souboru?
    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK) {
        basefilename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
        basefile = g_path_get_basename(basefilename);
        basedir = g_path_get_dirname(basefilename);

        // data
        buffer = fddata_export_data(&(args->fddata), EXPORT_FORMAT_PLAIN); /* */
        save_tool_export_data(basedir, basefile, "_data.txt", buffer); /* */
        g_free(buffer);

        // OP
        if (op_has_results(&args->opdata)) {
            buffer = op_export_data(&(args->opdata), &(args->fddata), &(args->instdata), &(args->area), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_op.txt", buffer); /* */
            g_free(buffer);
        }

        // OP Uncertainties & MC
        if (op_unc_has_results(&args->opunc)) {
            buffer = op_uncertainty_export_data(&(args->opdata), &(args->opunc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_op_unc.txt", buffer); /* */
            g_free(buffer);

            if (ctrls->opmccontrols->has_window) {
                buffer = op_mc_export_data(&(args->opdata), &(args->opunc), &(args->opmc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN); /* */
                save_tool_export_data(basedir, basefile, "_op_mc.txt", buffer); /* */
                g_free(buffer);
            }
        }

        // Tangent
        if (tangent_has_results(&args->tgdata)) {
            buffer = tangent_export_data(&(args->tgdata), &(args->fddata), &(args->instdata), &(args->area), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_tangent.txt", buffer); /* */
            g_free(buffer);
        }

        // Tangent Uncertainties & MC
        if (tangent_unc_has_results(&args->tgunc)) {
            buffer = tangent_uncertainty_export_data(&(args->tgdata), &(args->tgunc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_tangent_unc.txt", buffer); /* */
            g_free(buffer);

            if (ctrls->tgmccontrols->has_window) {
                buffer = tangent_mc_export_data(&(args->tgdata), &(args->tgunc), &(args->tgmc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN); /* */
                save_tool_export_data(basedir, basefile, "_tangent_mc.txt", buffer); /* */
                g_free(buffer);
            }
        }

        // Hertz
        if (hertz_has_results(&args->hertzdata)) {
            buffer = hertz_export_data(&(args->hertzdata), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_hertz.txt", buffer); /* */
            g_free(buffer);
        }

        // Hertz Uncertainties & MC
        if (hertz_unc_has_results(&args->hertzunc)) {
            buffer = hertz_uncertainty_export_data(&(args->hertzdata), &(args->hertzunc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_hertz_unc.txt", buffer); /* */
            g_free(buffer);

            if (ctrls->hertzmccontrols->has_window) {
                buffer = hertz_mc_export_data(&(args->hertzdata), &(args->hertzunc), &(args->hertzmc), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN); /* */
                save_tool_export_data(basedir, basefile, "_hertz_mc.txt", buffer); /* */
                g_free(buffer);
            }
        }

        // Diff Hard.
        if (Ph2_has_results(&args->ph2data)) {
            buffer = Ph2_export_data(&(args->ph2data), &(args->fddata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_Ph2.txt", buffer); /* */
            g_free(buffer);
        }

        // Pop-in
        if (apopin_has_results(&args->apopdata)) {
            buffer = apopin_export_data(&(args->apopdata), &(args->fddata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_popin.txt", buffer); /* */
            g_free(buffer);
        }

        // Elastic-plastic Work
        if (work_has_results(&args->workdata)) {
            buffer = work_export_data(&(args->workdata), &(args->fddata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_work.txt", buffer); /* */
            g_free(buffer);
        }

        // Elastic- plastic Work Uncertainties
        if (work_unc_has_results(&args->workunc)) {
            buffer = work_uncertainty_export_data(&(args->workunc), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_work_unc.txt", buffer); /* */
            g_free(buffer);
        }

#ifndef NOFORTRAN

        // OP ODR
        if (op_odr_has_results(&args->opodrdata)) {
            buffer = op_odr_export_data(&(args->opodrdata), &(args->fddata), &(args->instdata), &(args->area), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_op_odr.txt", buffer); /* */
            g_free(buffer);
        }

        // OP ODR Uncertainties & MC
        if (op_odr_unc_has_results(&args->opodrunc)) {
            buffer = op_odr_uncertainty_export_data(&(args->opodrdata), &(args->opodrunc), &(args->fddata), &(args->area), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_op_odr_unc.txt", buffer); /* */
            g_free(buffer);

            if (ctrls->opodrmccontrols->has_window) {
                buffer = op_odr_mc_export_data(&(args->opodrdata), &(args->opodrunc), &(args->opodrmc), &(args->fddata), EXPORT_FORMAT_PLAIN); /* */
                save_tool_export_data(basedir, basefile, "_op_odr_mc.txt", buffer); /* */
                g_free(buffer);
            }
        }

        // Two Slopes
        if (slopes_has_results(&args->slopesdata)) {
            buffer = slopes_export_data(&(args->slopesdata), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_slopes.txt", buffer); /* */
            g_free(buffer);
        }

        // Two Slopes Uncertainties & MC
        if (slopes_unc_has_results(&args->slopesunc)) {
            buffer = slopes_uncertainty_export_data(&(args->slopesdata), &(args->slopesunc), &(args->fddata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_slopes_unc.txt", buffer); /* */
            g_free(buffer);

            if (ctrls->slopesmccontrols->has_window) {
                buffer = slopes_mc_export_data(&(args->slopesdata), &(args->slopesunc), &(args->slopesmc), &(args->fddata), EXPORT_FORMAT_PLAIN); /* */
                save_tool_export_data(basedir, basefile, "_slopes_mc.txt", buffer); /* */
                g_free(buffer);
            }
        }

        // Hertz ODR
        if (hertz_odr_has_results(&args->hertzodrdata)) {
            buffer = hertz_odr_export_data(&(args->hertzodrdata), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_hertz_odr.txt", buffer); /* */
            g_free(buffer);
        }

        // Hertz ODR Uncertainties & MC
        if (hertz_odr_unc_has_results(&args->hertzodrunc)) {
            buffer = hertz_odr_uncertainty_export_data(&(args->hertzodrdata), &(args->hertzodrunc), &(args->fddata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_hertz_odr_unc.txt", buffer); /* */
            g_free(buffer);

            if (ctrls->hertzodrmccontrols->has_window) {
                buffer = hertz_odr_mc_export_data(&(args->hertzodrdata), &(args->hertzodrunc), &(args->hertzodrmc), &(args->fddata), EXPORT_FORMAT_PLAIN); /* */
                save_tool_export_data(basedir, basefile, "_hertz_odr_mc.txt", buffer); /* */
                g_free(buffer);
            }
        }

        // Stiffness ODR
        if (stiffness_has_results(&args->stiffnessdata)) {
            buffer = stiffness_export_data(&(args->stiffnessdata), &(args->fddata), &(args->instdata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_stiffness.txt", buffer); /* */
            g_free(buffer);
        }

        // Stiffness ODR Uncertainties & MC
        if (stiffness_unc_has_results(&args->stiffnessunc)) {
            buffer = stiffness_uncertainty_export_data(&(args->stiffnessdata), &(args->stiffnessunc), &(args->fddata), EXPORT_FORMAT_PLAIN); /* */
            save_tool_export_data(basedir, basefile, "_stiffness_unc.txt", buffer); /* */
            g_free(buffer);

            if (ctrls->stiffnessmccontrols->has_window) {
                buffer = stiffness_mc_export_data(&(args->stiffnessdata), &(args->stiffnessunc), &(args->stiffnessmc), &(args->fddata), EXPORT_FORMAT_PLAIN); /* */
                save_tool_export_data(basedir, basefile, "_stiffness_mc.txt", buffer); /* */
                g_free(buffer);
            }
        }

#endif

        g_free(basefilename);
        g_free(basefile);
        g_free(basedir);
    }

    g_free(suggestdir);
    g_free(suggestname);
    gtk_widget_destroy(dialog);
}

static void load_data(GtkFileChooserButton *widget, Data *data)
{
    gchar *filename;
    enum DataFileFormat fileformat;

    filename = NULL;
    fileformat = lastfileformat;

    if (open_file_dialog(&filename, &fileformat, NULL)) {
        load_data_gui(filename, fileformat, data);
    }
    else if (verbose) {
        g_print("No data loaded \n");
    }
}

static void load_data_gui(const gchar *filename, enum DataFileFormat fileformat, Data *data)
{
    FDdata *fddata;
    ContactControls *contactcontrols;
    MainControls *maincontrols;

    char str[1024];
    FileReadStatus filereadstatus;
    gboolean splitok;

    fddata = &(data->args->fddata);
    contactcontrols = data->controls->contactcontrols;
    maincontrols = data->controls->maincontrols;

    filereadstatus = load_data_nogui(filename, fileformat, fddata, &splitok);

    switch (filereadstatus) {
    case READ_FILE_OK:

        fddata->fileformat = fileformat;
        lastfileformat = fileformat;

        /* adjust window title */
        g_snprintf(str, 1024, "%s - %s", g_get_application_name(), filename);
        gtk_window_set_title(GTK_WINDOW(maincontrols->window_main), str);

        if (!splitok) {
            show_warning("Could not split data into loading and unloading curves. Check your data.");
        }

        /* spravit: veci s dilcimi krivkami jenom je-li split OK */

        /* f-d */
        /* contact_zoom_restore_fd(data); */
        rs_zoom_restore(GWY_GRAPH(contactcontrols->graph_fd), &(contactcontrols->zoomstack_fd));
        update_models_fd(contactcontrols);

        update_fd_plot(contactcontrols, fddata);
        gwy_selection_clear(contactcontrols->selection_fd);

        /* f(t), d(t) */
        fddata->htime = gwy_data_line_duplicate(fddata->horig); /* LEAK */
        fddata->Ftime = gwy_data_line_duplicate(fddata->Forig); /* LEAK */
        gwy_data_line_multiply(fddata->htime, 1 / fddata->hmax);
        gwy_data_line_multiply(fddata->Ftime, 1 / fddata->Fmax);

        rs_zoom_restore(GWY_GRAPH(contactcontrols->graph_fd_time), &(contactcontrols->zoomstack_fd_time));
        update_models_fd_time(contactcontrols);

        update_fd_time_plot(contactcontrols, fddata);
        update_point_marklines(contactcontrols, fddata);
        gwy_selection_clear(contactcontrols->selection_fd_time);

        /* prepare for work */
        gtk_label_set_text(GTK_LABEL(contactcontrols->status_contact_load), point_sel_status[fddata->status_contact_load]);
        gtk_label_set_text(GTK_LABEL(contactcontrols->status_load_hold), point_sel_status[fddata->status_load_hold]);
        gtk_label_set_text(GTK_LABEL(contactcontrols->status_hold_unload), point_sel_status[fddata->status_hold_unload]);
        gtk_label_set_text(GTK_LABEL(contactcontrols->status_end_unload), point_sel_status[fddata->status_end_unload]);

        update_max_labels(contactcontrols, fddata);
        data->functions->remove_all_tool_fit_results_labels(data, FALSE);
        restore_zoom_all_tool_graphs(data);
        data->functions->redraw_all_tools(data);

        contact_switch_to_fd(data);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(contactcontrols->radio_fd), TRUE);
        notebook_set_sensitive_all_pages(maincontrols, TRUE);
        gtk_notebook_set_current_page(GTK_NOTEBOOK(maincontrols->notebook_tools), CONTACT);
        gtk_widget_set_sensitive(contactcontrols->button_save_data, TRUE);
        break;

    case READ_FILE_COULD_NOT_OPEN:
        show_error("Could not open file.");
        break;

    case READ_FILE_NOT_ENOUGH_ENTRIES:
        show_error("File does not contain enough entries with specified data format.");
        break;

    case READ_FILE_FORMAT_MISMATCH:
        show_error("Data do not match the chosen file format.");
        break;

    case READ_FILE_CONTAINS_NAN:
        show_error("File contains NaN values.");
        break;

    case READ_FILE_CONTAINS_INF:
        show_error("File contains infinite values.");
        break;
    } /* switch filestatus */
}

static gboolean open_file_dialog(gchar **filename, enum DataFileFormat *fileformat, GtkWidget *parent)
{
    GtkWidget *dialog, *extra_hbox, *cbox_fileformat;
    gchar *fn;
    enum DataFileFormat ff;
    FileReadStatus filestatus;
    gint i, response;
    gboolean file_ok;

    dialog = gtk_file_chooser_dialog_new("Open File",
                                         GTK_WINDOW(parent),
                                         GTK_FILE_CHOOSER_ACTION_OPEN,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                         NULL);

    cbox_fileformat = gtk_combo_box_text_new();

    for (i = 0; i < ndatafileformats; i++) {
        gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(cbox_fileformat), data_file_format_names[i]);
    }

    gtk_combo_box_set_active(GTK_COMBO_BOX(cbox_fileformat), *fileformat);

    extra_hbox = gtk_hbox_new(FALSE, CTRLS_SPACING);
    gtk_box_pack_end(GTK_BOX(extra_hbox), cbox_fileformat, FALSE, FALSE, 0);
    gtk_box_pack_end(GTK_BOX(extra_hbox), gtk_label_new("Input data format"), FALSE, FALSE, 0);

    gtk_file_chooser_set_extra_widget(GTK_FILE_CHOOSER(dialog), extra_hbox);
    gtk_widget_show_all(extra_hbox);

    file_ok = FALSE;
    response = gtk_dialog_run(GTK_DIALOG(dialog));

    while ((response == GTK_RESPONSE_ACCEPT) && !file_ok) {
        fn = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
        ff = gtk_combo_box_get_active(GTK_COMBO_BOX(cbox_fileformat));

        if (verbose) {
            g_print("Picked filename %s, fileformat %d (%s)\n", fn, ff, data_file_format_names[ff]);
        }

        file_ok = file_contains_data(fn, ff, &filestatus);

        if (file_ok) {
            if (*filename != NULL) {
                g_free(*filename);
            }

            *filename = g_strdup(fn);
            *fileformat = ff;
        }
        else {
            switch (filestatus) {
            case READ_FILE_COULD_NOT_OPEN:
                show_error("Could not open file.");
                break;

            case READ_FILE_FORMAT_MISMATCH:
                show_error("File header does not conform to Niget data file format.");
                break;

            case READ_FILE_NOT_ENOUGH_ENTRIES:
                show_error("File does not contain enough entries with specified data format.");
                break;

            default:
                break;
            }

            response = gtk_dialog_run(GTK_DIALOG(dialog));
        }

        g_free(fn);
    }

    gtk_widget_destroy(dialog);

    return file_ok;
}

static void remove_all_tool_fit_results_labels(Data *data, gboolean areadependent)
{
    op_remove_fit_results_labels(data);

#ifndef NOFORTRAN
    op_odr_remove_fit_results_labels(data);

    if (!areadependent) {
        slopes_remove_fit_results_labels(data);
        hertz_odr_remove_fit_results_labels(data);
        stiffness_remove_fit_results_labels(data);
    }

#endif

    tangent_remove_fit_results_labels(data);

    if (!areadependent) {
        hertz_remove_fit_results_labels(data);
        Ph2_remove_fit_results_labels(data);
        apopin_remove_fit_results_labels(data);
        work_remove_fit_results_labels(data);
    }
}

static void redraw_all_tools(Data *data)
{
    op_redraw(data);
#ifndef NOFORTRAN
    op_odr_redraw(data);
#endif
    tangent_redraw(data);
    hertz_redraw(data);
#ifndef NOFORTRAN
    hertz_odr_redraw(data);
    stiffness_redraw(data);
    slopes_redraw(data);
#endif
    Ph2_redraw(data);
    apopin_redraw(data);
    work_redraw(data);
}

static void quit(Data *data)
{
    save_settings(g_build_filename(g_get_user_config_dir(), settings_filename, NULL), data->args);
    gtk_main_quit();
}

static void show_help(void)
{
#if defined(G_OS_WIN32) && (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 54)
    /* glib 2.54 proved working, older versions not available to check */

    /* Should be fixed in glib */
    if (!((INT_PTR)ShellExecute(NULL, "open", "http://nanometrologie.cz/niget/doc/html4/", NULL, NULL, SW_SHOWNORMAL) > 32)) {
        show_error("Failed to open online help.");
    }

#else
    GError *error = NULL;

    g_app_info_launch_default_for_uri("http://nanometrologie.cz/niget/doc/html4/", NULL, &error);

    if (error) {
        g_printerr("Could not open help file, error code: %d, message %s\n", error->code, error->message);
    }

#endif
}

static void notebook_set_sensitive_all_pages(MainControls *maincontrols, gboolean val)
{
    gtk_widget_set_sensitive(GTK_WIDGET(maincontrols->tab_contact), val);
    gtk_widget_set_sensitive(GTK_WIDGET(maincontrols->tab_op), val);
#ifndef NOFORTRAN
    gtk_widget_set_sensitive(GTK_WIDGET(maincontrols->tab_op_odr), val);
    gtk_widget_set_sensitive(GTK_WIDGET(maincontrols->tab_slopes), val);
    gtk_widget_set_sensitive(GTK_WIDGET(maincontrols->tab_hertz_odr), val);
    gtk_widget_set_sensitive(GTK_WIDGET(maincontrols->tab_stiffness), val);
#endif
    gtk_widget_set_sensitive(GTK_WIDGET(maincontrols->tab_tangent), val);
    gtk_widget_set_sensitive(GTK_WIDGET(maincontrols->tab_hertz), val);
    gtk_widget_set_sensitive(GTK_WIDGET(maincontrols->tab_ph2), val);
    gtk_widget_set_sensitive(GTK_WIDGET(maincontrols->tab_apopins), val);
    gtk_widget_set_sensitive(GTK_WIDGET(maincontrols->tab_epwork), val);
}

static void restore_zoom_all_tool_graphs(Data *data)
{
    op_zoom_restore(data);
    tg_zoom_restore(data);
    hz_zoom_restore(data);
    ph2_zoom_restore_linked(data);
    apop_zoom_restore_linked(data);
    work_zoom_restore(data);
#ifndef NOFORTRAN
    op_odr_zoom_restore(data);
    slopes_zoom_restore(data);
    hz_odr_zoom_restore(data);
#endif
}

int main(int argc, char **argv)
{
    Args args;
    Controls controls;
    struct Functions functions;

    Data data;

    MainControls maincontrols;
    ContactControls contactcontrols;
    AreaControls areacontrols;
    OPControls opcontrols;
    OPUncControls opunccontrols;
    OPMCControls opmccontrols;
#ifndef NOFORTRAN
    OPODRControls opodrcontrols;
    OPODRUncControls opodrunccontrols;
    OPODRMCControls opodrmccontrols;
    SlopesControls slopescontrols;
    SlopesUncControls slopesunccontrols;
    SlopesMCControls slopesmccontrols;
    HertzODRControls hertzodrcontrols;
    HertzODRUncControls hertzodrunccontrols;
    HertzODRMCControls hertzodrmccontrols;
    StiffnessControls stiffnesscontrols;
    StiffnessUncControls stiffnessunccontrols;
    StiffnessMCControls stiffnessmccontrols;
#endif
    TangentControls tgcontrols;
    TangentUncControls tgunccontrols;
    TangentMCControls tgmccontrols;
    HertzControls hzcontrols;
    HertzUncControls hzunccontrols;
    HertzMCControls hzmccontrols;
    Ph2Controls ph2controls;
    ApopinControls apopcontrols;
    WorkControls workcontrols;
    WorkUncControls workunccontrols;

    GdkScreen *screen;

    static gchar *fformat = NULL;
    static gchar **filenames = NULL;
    GError *error = NULL;
    GOptionContext *optioncontext;

    enum DataFileFormat fileformat;
    gboolean cmdlineok, fileformatok;

    /* be non-verbose by default */
    verbose = FALSE;

    static GOptionEntry entries[] = {
        { "fileformat", 'f', 0, G_OPTION_ARG_STRING, &fformat, "Data file format: niget (default), unht, hysitron", "FORMAT" },
        { "verbose", 'v', 0, G_OPTION_ARG_NONE, &verbose, "Be verbose", NULL },
        { G_OPTION_REMAINING, 0, 0, G_OPTION_ARG_FILENAME_ARRAY, &filenames, "Data file names", "filename1 filename2 …" },
        { NULL }
    };

    optioncontext = g_option_context_new("(only 1 file supported at present)");
    g_option_context_add_main_entries(optioncontext, entries, NULL);
    g_option_context_add_group(optioncontext, gtk_get_option_group(TRUE));
    g_option_context_set_summary(optioncontext, NIGET_LONGNAME ", version " NIGET_VERSION "\n"
                                 NIGET_URL);
    g_option_context_set_description(optioncontext, NIGET_LICENSE "\n"
                                     "Authors: " NIGET_AUTHORS "\n"
                                     "Please report any bugs to rslesinger@cmi.cz\n");

    cmdlineok = g_option_context_parse(optioncontext, &argc, &argv, &error);
    g_option_context_free(optioncontext);

    if (!cmdlineok) {
        g_print("Failed to parse commandline options: %s\n", error->message);
        exit(1);
    }

    fileformat = DATA_FILE_FORMAT_NIGET; /* default */

    if (fformat) {
        if (g_ascii_strcasecmp(fformat, "niget") == 0) {
            fileformatok = TRUE;
        }
        else if (g_ascii_strcasecmp(fformat, "unht") == 0) {
            fileformat = DATA_FILE_FORMAT_UNHT_D_NM_F_MN;
            fileformatok = TRUE;
        }
        else if (g_ascii_strcasecmp(fformat, "hysitron") == 0) {
            fileformat = DATA_FILE_FORMAT_HYSITRON;
            fileformatok = TRUE;
        }
        else {
            fileformatok = FALSE;
        }

        if (!fileformatok) {
            g_printerr("File type not recognized. Currently available options are: niget, unht, hysitron\n");
            exit(1);
        }
    }

    /* Initialize the widget set */
    gtk_init(&argc, &argv);
    g_set_application_name(NIGET_SHORTNAME);
    gwy_stock_register_stock_items();  // Gwyddion stock icons

    /* #ifndef G_OS_WIN32 */
    /* handle resources properly */
    /* g_resources_register(resources_get_resource()); */
    /* #endif */

    controls.maincontrols = &maincontrols;
    controls.contactcontrols = &contactcontrols;
    controls.areacontrols = &areacontrols;
    controls.opcontrols = &opcontrols;
    controls.opunccontrols = &opunccontrols;
    controls.opmccontrols = &opmccontrols;
#ifndef NOFORTRAN
    controls.opodrcontrols = &opodrcontrols;
    controls.opodrunccontrols = &opodrunccontrols;
    controls.opodrmccontrols = &opodrmccontrols;
    controls.slopescontrols = &slopescontrols;
    controls.slopesunccontrols = &slopesunccontrols;
    controls.slopesmccontrols = &slopesmccontrols;
    controls.hertzodrcontrols = &hertzodrcontrols;
    controls.hertzodrunccontrols = &hertzodrunccontrols;
    controls.hertzodrmccontrols = &hertzodrmccontrols;
    controls.stiffnesscontrols = &stiffnesscontrols;
    controls.stiffnessunccontrols = &stiffnessunccontrols;
    controls.stiffnessmccontrols = &stiffnessmccontrols;
#endif
    controls.tgcontrols = &tgcontrols;
    controls.tgunccontrols = &tgunccontrols;
    controls.tgmccontrols = &tgmccontrols;
    controls.hertzcontrols = &hzcontrols;
    controls.hertzunccontrols = &hzunccontrols;
    controls.hertzmccontrols = &hzmccontrols;
    controls.ph2controls = &ph2controls;
    controls.apopcontrols = &apopcontrols;
    controls.workcontrols = &workcontrols;
    controls.workunccontrols = &workunccontrols;

    init_args(&args);

    functions.redraw_all_tools = redraw_all_tools;
    functions.remove_all_tool_fit_results_labels = remove_all_tool_fit_results_labels;

    data.args = &args;
    data.controls = &controls;
    data.functions = &functions;

    load_settings(g_build_filename(g_get_user_config_dir(), settings_filename, NULL), data.args);


    /* Main window */
    maincontrols.window_main = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(maincontrols.window_main), NIGET_SHORTNAME " - " NIGET_LONGNAME);

    screen = gtk_window_get_screen(GTK_WINDOW(maincontrols.window_main));
    gtk_window_set_default_size(GTK_WINDOW(maincontrols.window_main),
                                (gint)(0.9 * gdk_screen_get_width(screen)), (gint)(0.9 * gdk_screen_get_height(screen)));
    gtk_window_set_position(GTK_WINDOW(maincontrols.window_main), GTK_WIN_POS_CENTER);

    g_signal_connect_swapped(G_OBJECT(maincontrols.window_main), "destroy", G_CALLBACK(quit), &data);

    /* Main notebook and its contents */
    maincontrols.notebook_tools = gtk_notebook_new();
    g_object_set(maincontrols.notebook_tools, "scrollable", TRUE, NULL);

    maincontrols.tab_contact = tool_contact_create(&data);
    maincontrols.tab_contact_label = gtk_label_new("Data & instrument");
    gtk_notebook_append_page(GTK_NOTEBOOK(maincontrols.notebook_tools), maincontrols.tab_contact, maincontrols.tab_contact_label);

    maincontrols.tab_op = tool_op_create(&data);
    maincontrols.tab_op_label = gtk_label_new("Oliver-Pharr analysis");
    gtk_notebook_append_page(GTK_NOTEBOOK(maincontrols.notebook_tools), maincontrols.tab_op, maincontrols.tab_op_label);

#ifndef NOFORTRAN
    maincontrols.tab_op_odr = tool_op_odr_create(&data); //
    maincontrols.tab_op_odr_label = gtk_label_new("Oliver-Pharr analysis with ODR");
    gtk_notebook_append_page(GTK_NOTEBOOK(maincontrols.notebook_tools), maincontrols.tab_op_odr, maincontrols.tab_op_odr_label);
#endif

    maincontrols.tab_tangent = tool_tangent_create(&data); //
    maincontrols.tab_tangent_label = gtk_label_new("Tangent analysis");
    gtk_notebook_append_page(GTK_NOTEBOOK(maincontrols.notebook_tools), maincontrols.tab_tangent, maincontrols.tab_tangent_label);

    maincontrols.tab_hertz = tool_hertz_create(&data); //
    maincontrols.tab_hertz_label = gtk_label_new("Hertz model");
    gtk_notebook_append_page(GTK_NOTEBOOK(maincontrols.notebook_tools), maincontrols.tab_hertz, maincontrols.tab_hertz_label);

#ifndef NOFORTRAN
    maincontrols.tab_hertz_odr = tool_hertz_odr_create(&data); //
    maincontrols.tab_hertz_odr_label = gtk_label_new("Hertz model ODR");
    gtk_notebook_append_page(GTK_NOTEBOOK(maincontrols.notebook_tools), maincontrols.tab_hertz_odr, maincontrols.tab_hertz_odr_label);

    maincontrols.tab_stiffness = tool_stiffness_create(&data); //
    maincontrols.tab_stiffness_label = gtk_label_new("Stiffness ODR");
    gtk_notebook_append_page(GTK_NOTEBOOK(maincontrols.notebook_tools), maincontrols.tab_stiffness, maincontrols.tab_stiffness_label);

    maincontrols.tab_slopes = tool_slopes_create(&data); //
    maincontrols.tab_slopes_label = gtk_label_new("Two slopes model");
    gtk_notebook_append_page(GTK_NOTEBOOK(maincontrols.notebook_tools), maincontrols.tab_slopes, maincontrols.tab_slopes_label);
#endif

    maincontrols.tab_ph2 = tool_Ph2_create(&data); //
    maincontrols.tab_ph2_label = label_new_with_markup("P-h<sup>2</sup> analysis");
    gtk_notebook_append_page(GTK_NOTEBOOK(maincontrols.notebook_tools), maincontrols.tab_ph2, maincontrols.tab_ph2_label);

    maincontrols.tab_apopins = tool_apopins_create(&data); //
    maincontrols.tab_apopins_label = gtk_label_new("Pop-in detection");
    gtk_notebook_append_page(GTK_NOTEBOOK(maincontrols.notebook_tools), maincontrols.tab_apopins, maincontrols.tab_apopins_label);

    maincontrols.tab_epwork = tool_epwork_create(&data); //
    maincontrols.tab_epwork_label = gtk_label_new("Elastic/plastic work");
    gtk_notebook_append_page(GTK_NOTEBOOK(maincontrols.notebook_tools), maincontrols.tab_epwork, maincontrols.tab_epwork_label);


    /* About */
    about_dialog_create(&data);

    /* Buttons */

    maincontrols.button_open = gtk_button_new_from_stock(GTK_STOCK_OPEN);
    gtk_widget_set_tooltip_text(maincontrols.button_open, "Open file containing force-depth indentation data");
    g_signal_connect(maincontrols.button_open, "clicked", G_CALLBACK(load_data), &data);
    lastfileformat = DATA_FILE_FORMAT_NIGET;

    //maincontrols.button_save = gtk_button_new_from_stock (GTK_STOCK_SAVE);
    maincontrols.button_save = gtk_button_new_with_mnemonic("Export _all data & results");
    gtk_button_set_image(GTK_BUTTON(maincontrols.button_save), gtk_image_new_from_stock(GTK_STOCK_SAVE, GTK_ICON_SIZE_BUTTON));
    gtk_widget_set_tooltip_text(maincontrols.button_save, "Export all data and results in a human-readable text format");
    g_signal_connect_swapped(maincontrols.button_save, "clicked", G_CALLBACK(save_data), &data);

    //maincontrols.button_help = gtk_button_new_from_stock (GTK_STOCK_HELP);
    maincontrols.button_help = gtk_button_new_with_mnemonic("H_elp");
    gtk_button_set_image(GTK_BUTTON(maincontrols.button_help), gtk_image_new_from_stock(GTK_STOCK_HELP, GTK_ICON_SIZE_BUTTON));
    gtk_widget_set_tooltip_text(maincontrols.button_help, "Show HTML documentation in a web browser");
    g_signal_connect(maincontrols.button_help, "clicked", G_CALLBACK(show_help), NULL);

    //maincontrols.button_about = gtk_button_new_from_stock (GTK_STOCK_ABOUT);
    maincontrols.button_about = gtk_button_new_with_mnemonic("Abou_t");
    gtk_button_set_image(GTK_BUTTON(maincontrols.button_about), gtk_image_new_from_stock(GTK_STOCK_ABOUT, GTK_ICON_SIZE_BUTTON));
    gtk_widget_set_tooltip_text(maincontrols.button_about, "Show basic information on the software");
    g_signal_connect(maincontrols.button_about, "clicked", G_CALLBACK(show_about), &data);

    maincontrols.button_quit = gtk_button_new_from_stock(GTK_STOCK_QUIT);
    g_signal_connect_swapped(maincontrols.button_quit, "clicked", G_CALLBACK(quit), &data);

    maincontrols.vbox_main = gtk_vbox_new(FALSE, CTRLS_SPACING);
    maincontrols.hbox_buttons = gtk_hbox_new(FALSE, CTRLS_SPACING);

    gtk_box_pack_start(GTK_BOX(maincontrols.hbox_buttons), maincontrols.button_open, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(maincontrols.hbox_buttons), maincontrols.button_save, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(maincontrols.hbox_buttons), maincontrols.button_help, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(maincontrols.hbox_buttons), maincontrols.button_about, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(maincontrols.hbox_buttons), maincontrols.button_quit, FALSE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(maincontrols.vbox_main), maincontrols.notebook_tools, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(maincontrols.vbox_main), maincontrols.hbox_buttons, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(maincontrols.window_main), maincontrols.vbox_main);

    notebook_set_sensitive_all_pages(&maincontrols, FALSE);
    gtk_widget_show_all(maincontrols.window_main);
    gtk_notebook_set_current_page(GTK_NOTEBOOK(maincontrols.notebook_tools), CONTACT);  /* until we find what causes OP tool to be active by default */

    if (filenames) {
        if (filenames[0]) { /* only 1 file can be loaded at present */
            load_data_gui(filenames[0], fileformat, &data);
        }

        g_strfreev(filenames);
    }

    gtk_main();

    return 0;
}

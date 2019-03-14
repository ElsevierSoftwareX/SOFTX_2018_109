#ifndef GUI_UTILS_H
#define GUI_UTILS_H

#include "controls.h"
#include "file-utils.h"

#include <gtk/gtk.h>
#include <libgwydgets/gwydgets.h>

#define CURVE_COLOR_END_UNLOAD gwy_graph_get_preset_color(0)
#define CURVE_COLOR_LOAD gwy_graph_get_preset_color(1)
#define CURVE_COLOR_HOLD gwy_graph_get_preset_color(2)
#define CURVE_COLOR_UNLOAD gwy_graph_get_preset_color(3)

#define CURVE_COLOR_H gwy_graph_get_preset_color(4)
#define CURVE_COLOR_F gwy_graph_get_preset_color(5)
#define CURVE_COLOR_FC gwy_graph_get_preset_color(6)

#define CURVE_COLOR_FIT gwy_graph_get_preset_color(4)
#define CURVE_COLOR_FIT2 gwy_graph_get_preset_color(5)

/* tool_apopins */
#define CURVE_COLOR_DER gwy_graph_get_preset_color(1)
#define CURVE_COLOR_SMOOTH_LOAD gwy_graph_get_preset_color(2)
#define CURVE_COLOR_TH gwy_graph_get_preset_color(2)
#define CURVE_COLOR_TH2 gwy_graph_get_preset_color(3)
/* tool_Ph2 */
#define CURVE_COLOR_PH2 gwy_graph_get_preset_color(1)
#define CURVE_COLOR_DIFF_PH2 gwy_graph_get_preset_color(2)


#define CTRLS_BORDER 3
#define CTRLS_SPACING 5
#define CTRLS_TABLE_COL_SPACING 10
#define CTRLS_TABLE_ROW_SPACING 3

struct Functions {
    void (*redraw_all_tools)(Data *data);
    void (*remove_all_tool_fit_results_labels)(Data *data, gboolean areadependent);
};

GtkWidget *label_new_left(const gchar *text);
GtkWidget *label_new_with_markup(const gchar *text);
GtkWidget *label_new_with_markup_with_mnemonic(const gchar *text);
GtkWidget *label_new_with_markup_left(const gchar *text);
void entry_clear(GtkWidget *entry);
void label_clear(GtkWidget *label);
void entry_set_gdouble_format(GtkWidget *entry, const gdouble num, const gchar *format);
void entry_set_gint_format(GtkWidget *entry, const gint num, const gchar *format);
void entry_set_text_format(GtkWidget *entry, const gchar *text, const gchar *format);
void label_set_gdouble_format(GtkWidget *label, const gdouble num, const gchar *format);
void label_set_gint_format(GtkWidget *label, const gint num, const gchar *format);
void label_set_text_format(GtkWidget *label, const gchar *text, const gchar *format);

void remove_table_content(GtkWidget *table);

void show_error(const gchar *msg);
void show_warning(const gchar *msg);
void show_save_error(FileSaveStatus filesavestatus, const gchar *filename);

void show_histogram(const gchar *quantity, const gchar *unit, const gdouble *xdata, const gdouble *ydata, gint ndata, GtkWidget *parent);
void show_log(GtkWidget *log_window, const gchar *logfnm, const gchar *logext, const gchar *window_title);

void buffer_to_file_dialog(const gchar *buffer);

void rs_zoom_in(GwyGraph *graph, GSList **stack, gboolean free_y_rescale);
void rs_zoom_out(GwyGraph *graph, GSList **stack);
void rs_zoom_restore(GwyGraph *graph, GSList **stack);

gboolean beta_changed(GtkWidget *wbeta, gdouble *beta);
void range_from_changed(GtkWidget *wfrom, GtkWidget *wto, GwySelection *selection, gdouble *datafrom, gdouble *data2from);
void range_to_changed(GtkWidget *wfrom, GtkWidget *wto, GwySelection *selection, gdouble *datato, gdouble *data2to);
gboolean range_from_pct_Fmax_changed(GtkWidget *from_pct_Fmax, GtkWidget *to_pct_Fmax, GwySelection *selection,
                                     const gdouble Fmax, GwyDataLine *hdataline, GwyDataLine *Fdataline,
                                     gdouble *datafrom, gdouble *data2from, gdouble *datato, gdouble *data2to,
                                     gboolean *Finput, gboolean *Finput2);
gboolean range_to_pct_Fmax_changed(GtkWidget *from_pct_Fmax, GtkWidget *to_pct_Fmax, GwySelection *selection,
                                   const gdouble Fmax, GwyDataLine *hdataline, GwyDataLine *Fdataline,
                                   gdouble *datafrom, gdouble *data2from, gdouble *datato, gdouble *data2to,
                                   gboolean *Finput, gboolean *Finput2);

#endif

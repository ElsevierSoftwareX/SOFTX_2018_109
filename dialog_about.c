#include "dialog_about.h"

#include "controls.h"
#include "datatypes.h"
#include "niget-common.h"

#include <gtk/gtk.h>

GtkWidget *dialog_about;

void show_about(GtkWidget *widget, Data *data)
{
    gtk_dialog_run(GTK_DIALOG(dialog_about));
    gtk_widget_hide(GTK_WIDGET(dialog_about));
}

void about_dialog_create(Data *data)
{
    /* #ifndef G_OS_WIN32 */
    /*     GdkPixbuf *logo_pixbuf; */
    /*     GError *error = NULL; */
    /* #endif */
    const gchar *authors[] = {"Anna Charvátová Campbell",
                              "Radek Šlesinger",
                              "Petr Grolich",
                              "",
                              "ODR function fitting:",
                              "J. W. Zwolak, P. T. Boggs, L. T. Watson: Algorithm 869: ODRPACK95: A weighted orthogonal distance regression code with bound constraints. ACM Trans. Math. Softw. (2007)",
                              NULL
                             };

    dialog_about = gtk_about_dialog_new();
    gtk_about_dialog_set_program_name(GTK_ABOUT_DIALOG(dialog_about), NIGET_SHORTNAME);
    gtk_about_dialog_set_comments(GTK_ABOUT_DIALOG(dialog_about), NIGET_LONGNAME);
    gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(dialog_about), NIGET_VERSION);
    gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG(dialog_about), NIGET_AUTHORS "\n" NIGET_CMI);
    gtk_about_dialog_set_license(GTK_ABOUT_DIALOG(dialog_about), NIGET_LICENSE);
    /* gtk_about_dialog_set_website requires gspawn-win32-helper-console.exe on Windows to launch the browser */
    gtk_about_dialog_set_website(GTK_ABOUT_DIALOG(dialog_about), NIGET_URL);
    gtk_about_dialog_set_website_label(GTK_ABOUT_DIALOG(dialog_about), NIGET_URL);
    gtk_about_dialog_set_authors(GTK_ABOUT_DIALOG(dialog_about), authors);
    gtk_window_set_icon_name(GTK_WINDOW(dialog_about), GTK_STOCK_ABOUT);

    /* #ifndef G_OS_WIN32 */
    /*     logo_pixbuf = gdk_pixbuf_new_from_resource ("/resources/logo.jpg", &error); */
    /*     gtk_about_dialog_set_logo (GTK_ABOUT_DIALOG(dialog_about), logo_pixbuf); */
    /* #endif */
}

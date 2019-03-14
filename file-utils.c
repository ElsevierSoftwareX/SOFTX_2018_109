#include "file-utils.h"

#include "niget-common.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <libprocess/gwyprocess.h>

#ifdef G_OS_WIN32
#include "getline.h"
#include <glib/gstdio.h>
#endif

static gchar *get_file_format_string(enum DataFileFormat fileformat, gint *nvars);

// unifikovat, dotahnout
gboolean file_contains_data(const gchar *filename, enum DataFileFormat fileformat, FileReadStatus *filestatus)
{
    FILE *infile;
    gchar *line = NULL;
    gchar *formatstring;
    size_t size;
    gint numdata, nvars, n;
    gdouble hlp1, hlp2, hlp3;

    infile = fopen(filename, "r");

    if (!infile) {
        *filestatus = READ_FILE_COULD_NOT_OPEN;
        return FALSE;
    }

    if (fileformat == DATA_FILE_FORMAT_NIGET) {
        getline(&line, &size, infile);

        if (sscanf(line, DATAFORMAT_HEADER_FORMAT, &n, &n, &n, &n) != 4) { /* do not care about the values */
            *filestatus = READ_FILE_FORMAT_MISMATCH;
            return FALSE;
        }
    }

    formatstring = get_file_format_string(fileformat, &nvars);

    numdata = 0;

    /* 3 points define an elementary loading-unloading curve */
    while ((getline(&line, &size, infile) != EOF) && (numdata < 3)) {
        localize_decimal_points(line);
        n = sscanf(line, formatstring, &hlp1, &hlp2, &hlp3);

        if (n == nvars) {
            numdata++;
        }
    }

    fclose(infile);
    g_free(formatstring);

    if (numdata == 3) {
        if (verbose) {
            g_print("File %s contains at least 3 valid data items, OK\n", filename);
        }

        *filestatus = READ_FILE_OK;
        return TRUE;
    }
    else {
        g_printerr("File %s contains only %d valid data items, not OK\n", filename, numdata);
        *filestatus = READ_FILE_NOT_ENOUGH_ENTRIES;
        return FALSE;
    }
}

/*
read_file_xxx (dle typu souboru - niget, plaintext, ...) jenom primitivni, zaplni tri zadane datalines - t,F,d
dalsi ukony jako uklid, inicializaci, split, zajisti load_data (non-gui verze)
*/

gint num_data_in_file(FILE *infile, enum DataFileFormat fileformat)
{
    gint numdata, n, nvars;
    size_t size;
    gdouble x1, x2, x3;
    gchar *formatstring, *line;

    formatstring = get_file_format_string(fileformat, &nvars);

    line = NULL;
    numdata = 0;

    while (getline(&line, &size, infile) != EOF) {
        localize_decimal_points(line);
        n = sscanf(line, formatstring, &x1, &x2, &x3);

        if (n == nvars) {
            numdata++;
        }
    }

    rewind(infile);

    if (verbose) {
        g_print("Number of data points in file: %d\n", numdata);
    }

    g_free(formatstring);

    return numdata;
}

/* have separate function for Niget files to allow future format changes */
FileReadStatus read_file_niget2(FILE *infile, gint numdata, GwyDataLine *time, GwyDataLine *hdata, GwyDataLine *Fdata,
                                gint *i_contact_load, gint *i_load_hold, gint *i_hold_unload, gint *i_end_unload)
{
    gchar *line = NULL;
    size_t size;

    getline(&line, &size, infile);
    localize_decimal_points(line);
    sscanf(line, DATAFORMAT_HEADER_FORMAT,
           i_contact_load, i_load_hold, i_hold_unload, i_end_unload);

    return read_file_plaintext2(infile, DATA_FILE_FORMAT_NIGET, numdata, time, hdata, Fdata);
}

FileReadStatus read_file_plaintext2(FILE *infile, enum DataFileFormat fileformat,
                                    gint numdata, GwyDataLine *time, GwyDataLine *hdata, GwyDataLine *Fdata)
{
    gint i, n, nvars;
    gchar *formatstring, *line;
    size_t size;
    gdouble x1, x2, x3;

    formatstring = get_file_format_string(fileformat, &nvars);
    line = NULL;
    i = 0;

    while (getline(&line, &size, infile) != EOF) {
        localize_decimal_points(line);
        n = sscanf(line, formatstring, &x1, &x2, &x3);

        if (n == nvars) {
            switch (fileformat) {
            case (DATA_FILE_FORMAT_UNHT_D_NM_F_MN):
                time->data[i] = i;
                hdata->data[i] = x2;
                Fdata->data[i] = x3;
                break;

            case (DATA_FILE_FORMAT_UNHT_D_NM_F_UN):
                time->data[i] = i;
                hdata->data[i] = x2;
                Fdata->data[i] = x3 / 1000;
                break;

            case (DATA_FILE_FORMAT_HYSITRON):
                time->data[i] = i;
                hdata->data[i] = x1;
                Fdata->data[i] = x2 / 1000;
                break;

            case (DATA_FILE_FORMAT_NIGET):
            case (DATA_FILE_FORMAT_D_NM_F_MN):
                time->data[i] = i;
                hdata->data[i] = x1;
                Fdata->data[i] = x2;
                break;

            case (DATA_FILE_FORMAT_F_MN_D_NM):
                time->data[i] = i;
                hdata->data[i] = x2;
                Fdata->data[i] = x1;
                break;

            case (DATA_FILE_FORMAT_F_UN_D_NM):
                time->data[i] = i;
                hdata->data[i] = x2;
                Fdata->data[i] = x1 / 1000;
                break;

            default:
                break;
            }

            i++;
        }
    }

    g_free(formatstring);

    if (i == numdata) {
        return READ_FILE_OK;
    }
    else {
        return READ_FILE_NOT_ENOUGH_ENTRIES;
    }
}

gboolean buffer_to_file(const gchar *filename, const gchar *buffer, FileSaveStatus *filestatus)
{
    FILE *file;

    if (!buffer) {
	g_printerr("No data in the buffer to save to file %s\n", filename);
	*filestatus = SAVE_FILE_NO_DATA;
	return FALSE;
    }
    
    file = fopen(filename, "w");

    if (!file) {
	g_printerr("Could not open file %s for writing\n", filename);
        *filestatus = SAVE_FILE_COULD_NOT_OPEN;
        return FALSE;
    }

    fputs(buffer, file);
    fclose(file);
    *filestatus = SAVE_FILE_OK;
    return TRUE;
}

/* static functions */

static gchar *get_file_format_string(enum DataFileFormat fileformat, gint *nvars)
{
    switch (fileformat) {
    case (DATA_FILE_FORMAT_UNHT_D_NM_F_MN):
    case (DATA_FILE_FORMAT_UNHT_D_NM_F_UN):
        *nvars = 3;
        return g_strdup("%lf %lf %lf");
        break;

    case (DATA_FILE_FORMAT_NIGET):
    case (DATA_FILE_FORMAT_HYSITRON):
    case (DATA_FILE_FORMAT_D_NM_F_MN):
    case (DATA_FILE_FORMAT_F_MN_D_NM):
    case (DATA_FILE_FORMAT_F_UN_D_NM):
        *nvars = 2;
        return g_strdup("%lf %lf");
        break;

    default:
        return NULL;
        break;
    }
}

/* Creates full log filename path from filename.
   Returned value must be deallocated.
*/
gchar *create_log_filename_path(const gchar *fname)
{
    gchar *res = NULL;

#ifdef G_OS_WIN32
    GStatBuf st;

    res = g_malloc(512 * sizeof(gchar));
    strcpy(res, getenv("APPDATA"));
    strncat(res, "\\Niget\\", strlen("\\Niget\\"));

    if (g_stat(res, &st) == -1) {
        g_mkdir(res, 0700);
    }

    strncat(res, fname, strlen(fname));
#else
    gint n;

    n = strlen(fname) ;
    res = g_malloc((n + 1) * sizeof(gchar));
    strncpy(res, fname, n + 1);
#endif

    return res;
}

gint get_number_lines(FILE *infile)
{
    gint n = 0 ;

    gchar *line = NULL;
    size_t size;

    while (getline(&line, &size, infile) > -1) {
        n++;
    }

    rewind(infile);
    return n;
}

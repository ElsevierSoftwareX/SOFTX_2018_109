#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include "datatypes.h"

#include <stdio.h>

typedef enum {
    READ_FILE_OK = 0,
    READ_FILE_COULD_NOT_OPEN,
    READ_FILE_NOT_ENOUGH_ENTRIES,
    READ_FILE_FORMAT_MISMATCH,
    READ_FILE_CONTAINS_NAN,
    READ_FILE_CONTAINS_INF
} FileReadStatus;

typedef enum {
    SAVE_FILE_OK = 0,
    SAVE_FILE_COULD_NOT_OPEN,
    SAVE_FILE_NO_DATA
} FileSaveStatus;

#define DATAFORMAT_HEADER_FORMAT "#NIGET %d %d %d %d"
#define DATAFORMAT_HEADER_FORMAT_NEWLINE "#NIGET %d %d %d %d\n"

gboolean file_contains_data(const gchar *filename, enum DataFileFormat fileformat, FileReadStatus *filestatus);
gint num_data_in_file(FILE *infile, enum DataFileFormat fileformat);

gboolean buffer_to_file(const gchar *filename, const gchar *buffer, FileSaveStatus *filestatus);

gchar *create_log_filename_path(const gchar *fname);
gint get_number_lines(FILE *infile);

FileReadStatus read_file_niget2(FILE *infile, gint numdata, GwyDataLine *time, GwyDataLine *hdata, GwyDataLine *Fdata,
				gint *i_contact_load, gint *i_load_hold, gint *i_hold_unload, gint *i_end_unload);
FileReadStatus read_file_plaintext2(FILE *infile, enum DataFileFormat fileformat,
				    gint numdata, GwyDataLine *time, GwyDataLine *hdata, GwyDataLine *Fdata);

#endif

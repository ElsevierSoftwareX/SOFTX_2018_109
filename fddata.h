#ifndef FDDATA_H
#define FDDATA_H

#include "datatypes.h"
#include "file-utils.h"
#include <libprocess/gwyprocess.h>

void init_FDdata(FDdata *fddata);
void destroy_FDdata(FDdata *fddata);

gboolean split_hF_data(FDdata *fddata);
gboolean split_hF_data_auto(FDdata *fddata);
gboolean split_hF_data_spec(FDdata *fddata, gint i_contact_load, gint i_load_hold, gint i_hold_unload, gint i_end_unload);

void dataline_clear(GwyDataLine *dataline);
void clear_hF_data(FDdata *fddata);

void fddata_get_hmax_Fmax(GwyDataLine *h, GwyDataLine *F, gdouble *hmax, gdouble *Fmax);
void fddata_find_contact(GwyDataLine *hdata, GwyDataLine *Fdata, gint *i_contact_load, gdouble *hc, gdouble *Fc);

gdouble create_shifted_FDdata(const FDdata *src, FDdata *dest, gint j);
gdouble create_shifted_FDdata_old(const FDdata *src, FDdata *dest, gint j);

FileReadStatus load_data_nogui(const gchar *filename, enum DataFileFormat fileformat, FDdata *fddata, gboolean *splitok);
gchar* fddata_export_data(const FDdata *fddata, enum ExportFormat exportformat);

#endif

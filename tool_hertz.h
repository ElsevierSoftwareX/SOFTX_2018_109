#ifndef TOOL_HERTZ_H
#define TOOL_HERTZ_H

#include "datatypes.h"
#include "file-utils.h"

void init_Hertzdata(Hertzdata *hzdata);
void hertz_fit(GwyDataLine *x, GwyDataLine *y, Hertzdata *hzdata, const Instdata *instdata);
gchar* hertz_export_data(const Hertzdata *hzdata, const FDdata *fddata, const Instdata *instdata, enum ExportFormat exportformat);
gboolean hertz_has_results(const Hertzdata *hzdata);
void hertz_remove_all_fit_results(Hertzdata *hzdata);

#endif

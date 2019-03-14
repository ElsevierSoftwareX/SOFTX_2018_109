#ifndef TOOL_PH2_H
#define TOOL_PH2_H

#include "datatypes.h"
#include "file-utils.h"

void init_Ph2data(Ph2data *ph2data);
GwyDataLine *copy_data_remove_ind(GwyDataLine *source, GwyDataLine *indices);
gchar* Ph2_export_data(const Ph2data *ph2data, const FDdata *fddata, enum ExportFormat exportformat);
gboolean Ph2_has_results(const Ph2data *ph2data);
void Ph2_calc(const FDdata *fddata, Ph2data *ph2data);

#endif

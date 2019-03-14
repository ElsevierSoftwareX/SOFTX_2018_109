#ifndef TOOL_EPWORK_H
#define TOOL_EPWORK_H

#include "datatypes.h"
#include "file-utils.h"

void init_Workdata(Workdata *workdata);
void work_calc(Workdata *workdata, const FDdata *fddata);
gchar* work_export_data(const Workdata *workdata, const FDdata *fddata, enum ExportFormat exportformat);
gboolean work_has_results(const Workdata *workdata);
void work_remove_results(Workdata *workdata);

#endif

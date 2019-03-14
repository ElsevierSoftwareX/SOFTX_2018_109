#ifndef TOOL_EPWORK_UNC_H
#define TOOL_EPWORK_UNC_H

#include "datatypes.h"
#include "file-utils.h"

#define CONTACTMAX 3

void init_WorkUncdata(WorkUncdata *workunc);
gboolean work_unc_has_results(const WorkUncdata *workuncdata);
gchar* work_uncertainty_export_data(const WorkUncdata *wunc, enum ExportFormat exportformat);
void work_run_shift_contact(const Workdata *wdata, const FDdata *fddata, gdouble *We, gdouble *Wp, gint j);

#endif

#ifndef NOFORTRAN

#ifndef TOOL_STIFFNESS_H
#define TOOL_STIFFNESS_H

#include "datatypes.h"
#include "file-utils.h"

void init_Stiffnessdata(Stiffnessdata *stdata);
gint stiffness_fit_unload(GwyDataLine *x, GwyDataLine *y, Stiffnessdata *stdata, const FDdata *fddata, const Instdata *instdata, const gdouble estimate[3]);
gint stiffness_fit_load(GwyDataLine *x, GwyDataLine *y, Stiffnessdata *stdata, const FDdata *fddata, const Instdata *instdata, const gdouble estimate[3]);
gchar* stiffness_export_data(const Stiffnessdata *stdata, const FDdata *fddata, const Instdata *instdata, enum ExportFormat exportformat);
void stiffness_combine_fit_results(Stiffnessdata *stdata, const FDdata *fddata, const Instdata *instdata);
gboolean stiffness_has_results(const Stiffnessdata *stiffness);
void stiffness_remove_all_fit_results_unload(Stiffnessdata *stdata);
void stiffness_remove_all_fit_results_load(Stiffnessdata *stdata);

#endif

#endif

#ifndef NOFORTRAN

#ifndef TOOL_SLOPES_H
#define TOOL_SLOPES_H

#include "datatypes.h"
#include "file-utils.h"

#define BETA_DEFAULT_SLOPES 1.0
#define N_FIX_SLOPES 2.0

void init_Slopesdata(Slopesdata *sldata);
gint slopes_fit_unload(GwyDataLine *x, GwyDataLine *y, Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata, const gdouble estimate[3]);
gint slopes_fit_load(GwyDataLine *x, GwyDataLine *y, Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata, const gdouble estimate[3]);
gchar* slopes_export_data(const Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata, enum ExportFormat exportformat) ;
void slopes_combine_fit_results(Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata);
gboolean slopes_has_results(const Slopesdata *slopes) ;
void slopes_remove_all_fit_results_load(Slopesdata *sldata);
void slopes_remove_all_fit_results_unload(Slopesdata *sldata);

#endif

#endif

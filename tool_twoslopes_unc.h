#ifndef NOFORTRAN

#ifndef TOOL_SLOPES_UNC_H
#define TOOL_SLOPES_UNC_H

#include "datatypes.h"
#include "file-utils.h"
#include "instrument.h"

#define NPARAM_SLOPES 8  //  Sload, Sunload,Hit, Er, Eit, Aphc, m, n
#define MAXBAD 5 //maximum percentage of nonsense results before warning is triggered

#define CONTACTMAX 3

void init_SlopesUncdata(SlopesUncdata *slopesunc, Instdata *instdata);
void slopes_propagate_uncertainties_calc(const Slopesdata *sldata, SlopesUncdata *unc, const FDdata *fd);
gboolean slopes_unc_has_results(const SlopesUncdata *slopesunc);
gchar *slopes_uncertainty_export_data(const Slopesdata *sldata, const SlopesUncdata *slunc, const FDdata *fddata, enum ExportFormat exportformat);
void slopes_fit_shift_contact(const Slopesdata *sldata, const FDdata *fddata, const Instdata *instdata, gdouble *Ec, gdouble *Hc, gdouble *Ac, gint j);

gchar *slopes_mc_export_data(const Slopesdata *slopesdata, const SlopesUncdata *slopesunc, const SlopesMCdata *slopesmc, const FDdata *fddata, enum ExportFormat exportformat);
void slopes_uncertainty_montecarlo_run_calc(const Slopesdata *slopesdata, SlopesUncdata *slopesunc, SlopesMCdata *slopesmc, FDdata *fddata);
void init_SlopesMCdata(SlopesMCdata *slopesmc, gint Nmc);
void destroy_SlopesMCdata(SlopesMCdata *slopesmc);

#endif

#endif

#ifndef TOOL_TANGENT_UNC_H
#define TOOL_TANGENT_UNC_H

#include "datatypes.h"
#include "file-utils.h"
#include "instrument.h"

#define CONTACTMAX 3
#define NPARAM_TG 6  // hc, Aphc, Hit, Er, Eit,S
#define MAXBAD 5 //maximum percentage of nonsense results before warning is triggered

void init_TangentUncdata(TangentUncdata *tgunc);

void tangent_propagate_uncertainties_calc(Tangentdata *tg, TangentUncdata *unc, FDdata *fd, Instdata *inst, Area *area);
gboolean tangent_unc_has_results(TangentUncdata *tgunc);
gchar* tangent_uncertainty_export_data(Tangentdata *tgdata, TangentUncdata *tgunc, FDdata *fddata, Instdata *instdata, enum ExportFormat exportformat);
void tangent_fit_shift_contact(Tangentdata *tgdata, FDdata *fddata, Instdata *inst, Area *area, gdouble *Sc, gdouble *Ec, gdouble *Hc, gint j);

gchar* tangent_mc_export_data(Tangentdata *tgdata, TangentUncdata *tgunc, TangentMCdata *tgmc, FDdata *fddata, Instdata *instdata, enum ExportFormat exportformat);
void tangent_uncertainty_montecarlo_run_calc(Tangentdata *tgdata, TangentUncdata *tgunc, TangentMCdata *tgmc, FDdata *fddata, Instdata *instdata, Area *area);
void init_TangentMCdata(TangentMCdata *tgmc, gint Nmc);
void destroy_TangentMCdata(TangentMCdata *tgmc);

#endif

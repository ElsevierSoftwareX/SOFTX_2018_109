#ifndef TOOL_HERTZ_UNC_H
#define TOOL_HERTZ_UNC_H

#include "datatypes.h"
#include "file-utils.h"
#include "instrument.h"

#define CONTACTMAX 3
#define MAXBAD 5 //maximum percentage of nonsense data
#define NPARAM_HZ 3  // Er, Eit or Radius

void init_HertzUncdata(HertzUncdata *hzunc);

void hertz_propagate_uncertainties_calc(const Hertzdata *hz, HertzUncdata *unc, const Instdata *inst);

gboolean hertz_unc_has_results(const HertzUncdata *hzunc);
gchar* hertz_uncertainty_export_data(const Hertzdata *hzdata, const HertzUncdata *hzunc, const FDdata *fddata,
				     const Instdata *instdata, enum ExportFormat exportformat);
void hertz_fit_shift_contact(const Hertzdata *hzdata, const FDdata *fddata, const Instdata *instdata, gdouble *ERc, gint j);

gchar* hertz_mc_export_data(const Hertzdata *hertzdata, const HertzUncdata *hertzunc, const HertzMCdata *hertzmc,
			    const FDdata *fddata, const Instdata *instdata, enum ExportFormat exportformat);
void hertz_uncertainty_montecarlo_run_calc(const Hertzdata *hertzdata, const HertzUncdata *hertzunc, HertzMCdata *hertzmc,
					   const FDdata *fddata, const Instdata *instdata);
void init_HertzMCdata(HertzMCdata *hzmc, gint Nmc);
void destroy_HertzMCdata(HertzMCdata *hzmc);

#endif

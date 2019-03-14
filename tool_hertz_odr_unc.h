#ifndef TOOL_HERTZ_ODR_UNC_H
#define TOOL_HERTZ_ODR_UNC_H

#include "datatypes.h"
#include "file-utils.h"
#include "instrument.h"

#define CONTACTMAX 3
#define MC_WHEN_WARN 1e2
#define NPARAM_HZODR 3  // Er, Eit or Radius

void init_HertzODRUncdata(HertzODRUncdata *hzodrunc, Instdata *inst);
void hertz_odr_propagate_uncertainties_calc(HertzODRdata *hzodr, HertzODRUncdata *unc, const FDdata *fd);
gboolean hertz_odr_unc_has_results(const HertzODRUncdata *hertzodrunc);
gchar *hertz_odr_uncertainty_export_data(const HertzODRdata *hzodrdata, const HertzODRUncdata *hzodrunc,
        const FDdata *fddata, enum ExportFormat exportformat);
void hertz_odr_fit_shift_contact(const HertzODRdata *hzodrdata, const FDdata *fddata, const Instdata *instdata, gdouble *ERc, gint j);

void init_HertzODRMCdata(HertzODRMCdata *hzodrmc, gint Nmc);
void destroy_HertzODRMCdata(HertzODRMCdata *hzodrmc);
void hertz_odr_uncertainty_montecarlo_run_calc(const HertzODRdata *hertzodrdata, HertzODRUncdata *hertzodrunc, HertzODRMCdata *hertzodrmc,
        const FDdata *fddata);
gchar *hertz_odr_mc_export_data(const HertzODRdata *hertzodrdata, const HertzODRUncdata *hertzodrunc, const HertzODRMCdata *hertzodrmc,
                                const FDdata *fddata, enum ExportFormat exportformat);

#endif

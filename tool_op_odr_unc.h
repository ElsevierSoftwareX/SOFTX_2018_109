#ifndef NOFORTRAN

#ifndef TOOL_OPODR_UNC_H
#define TOOL_OPODR_UNC_H

#include "datatypes.h"
#include "file-utils.h"
#include "instrument.h"

void init_OPODRUncdata(OPODRUncdata *opodrunc, Instdata *instdata);

void op_odr_propagate_uncertainties_calc(const OPODRdata *opodr, OPODRUncdata *unc, const FDdata *fd, const Area *area);

gboolean op_odr_unc_has_results(const OPODRUncdata *opodrunc);
gchar *op_odr_uncertainty_export_data(const OPODRdata *opodrdata, const OPODRUncdata *opodrunc, const FDdata *fddata, const Area *area,
                                      enum ExportFormat exportformat);
void op_odr_fit_shift_contact(const OPODRdata *opodrdata, const FDdata *fddata, const Instdata *instdata, const Area *area,
                              gdouble *Ec, gdouble *Hc, gdouble *Ac, gdouble *hc, gint j);

void init_OPODRMCdata(OPODRMCdata *opodrmc, gint Nmc);
void destroy_OPODRMCdata(OPODRMCdata *opodrmc);

gchar *op_odr_mc_export_data(const OPODRdata *opodrdata, const OPODRUncdata *opodrunc, const OPODRMCdata *opodrmc,
                             const FDdata *fddata, enum ExportFormat exportformat);
void op_odr_uncertainty_montecarlo_run_calc(const OPODRdata *opodrdata, OPODRUncdata *opodrunc, OPODRMCdata *opodrmc,
        FDdata *fddata, const Area *area);

#define CONTACTMAX 3

#endif

#endif

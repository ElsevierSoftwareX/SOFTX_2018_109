#ifndef TOOL_OP_UNC_H
#define TOOL_OP_UNC_H

#include "datatypes.h"
#include "file-utils.h"
#include "instrument.h"

void init_OPUncdata(OPUncdata *opunc);

void op_propagate_uncertainties_calc(const OPdata *op, OPUncdata *unc, const FDdata *fd, const Instdata *inst, const Area *area);
void op_fit_shift_contact(const OPdata *opdata, const FDdata *fddata, const Instdata *instdata, const Area *area,
			  gdouble *Sc, gdouble *Ec, gdouble *Hc, gint j);

gchar* op_uncertainty_export_data(const OPdata *opdata, const OPUncdata *opunc, const FDdata *fddata, const Instdata *instdata,
				  enum ExportFormat exportformat);
gchar* op_mc_export_data(const OPdata *opdata, const OPUncdata *opunc, const OPMCdata *opmc, const FDdata *fddata, const Instdata *instdata,
			 enum ExportFormat exportformat);

void op_uncertainty_montecarlo_run_calc(const OPdata *opdata, const OPUncdata *opunc, OPMCdata *opmc,
					FDdata *fddata, const Instdata *instdata, const Area *area);
void init_OPMCdata(OPMCdata *opmc, gint Nmc);
void destroy_OPMCdata(OPMCdata *opmc);

gboolean op_unc_has_results(const OPUncdata *opuncdata);

#define NPARAM_OP 6  //hc, Aphc, Hit, Er, Eit, S
#define MAXBAD 5 //maximum percentage of nonsense results before warning is triggered

#define CONTACTMAX 3

#endif

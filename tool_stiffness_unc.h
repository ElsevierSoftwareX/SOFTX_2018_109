#ifndef NOFORTRAN

#ifndef TOOL_STIFFNESS_UNC_H
#define TOOL_STIFFNESS_UNC_H

#include "datatypes.h"
#include "file-utils.h"
#include "instrument.h"

#define NPARAM_STIFFNESS 4  //  kload, qload, kunload, qunload
#define MAXBAD 5 //maximum percentage of nonsense results before warning is triggered

#define CONTACTMAX 3

void init_StiffnessUncdata(StiffnessUncdata *stiffnessunc, Instdata *instdata);
void stiffness_propagate_uncertainties_calc(const Stiffnessdata *stdata, StiffnessUncdata *unc, const FDdata *fd);
gboolean stiffness_unc_has_results(const StiffnessUncdata *stiffnessunc);
gchar *stiffness_uncertainty_export_data(const Stiffnessdata *stdata, const StiffnessUncdata *stunc, const FDdata *fddata, enum ExportFormat exportformat);
void stiffness_fit_shift_contact(const Stiffnessdata *stdata, const FDdata *fddata, const Instdata *instdata,
                                 gdouble *klc, gdouble *qlc, gdouble *kuc, gdouble *quc, gint j);

gchar *stiffness_mc_export_data(const Stiffnessdata *stiffnessdata, const StiffnessUncdata *stiffnessunc, const StiffnessMCdata *stiffnessmc,
                                const FDdata *fddata, enum ExportFormat exportformat);
void stiffness_uncertainty_montecarlo_run_calc(const Stiffnessdata *stiffnessdata, StiffnessUncdata *stiffnessunc, StiffnessMCdata *stiffnessmc,
        FDdata *fddata);
void init_StiffnessMCdata(StiffnessMCdata *stiffnessmc, gint Nmc);
void destroy_StiffnessMCdata(StiffnessMCdata *stiffnessmc);

#endif

#endif

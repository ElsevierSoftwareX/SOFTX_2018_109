#ifndef TOOL_OP_H
#define TOOL_OP_H

#include "datatypes.h"
#include "file-utils.h"

#define FIT_S_LR  TRUE // TRUE - Linear REgreeession, FALSE - Nonlinear Regression

#define BETA_DEFAULT_OP 1.034

void init_OPdata(OPdata *opdata);

void op_remove_all_fit_results_hp(OPdata *opdata);
void op_remove_all_fit_results_S(OPdata *opdata);

void op_fit_hp(GwyDataLine *x, GwyDataLine *y, OPdata *opdata, const FDdata *fddata);
void op_fit_S_nl(GwyDataLine *x, GwyDataLine *y, OPdata *opdata, const FDdata *fddata, const Instdata *instdata, const Area *area);
void op_fit_S_lr(GwyDataLine *x, GwyDataLine *y, OPdata *opdata, const FDdata *fddata, const Instdata *instdata, const Area *area);

gchar* op_export_data(const OPdata *opdata, const FDdata *fddata, const Instdata *instdata, const Area *area, enum ExportFormat exportformat);

gboolean op_has_results(const OPdata *opdata);

#endif

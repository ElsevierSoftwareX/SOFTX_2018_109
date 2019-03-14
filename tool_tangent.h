#ifndef TOOL_TANGENT_H
#define TOOL_TANGENT_H

#include "datatypes.h"
#include "file-utils.h"

#define BETA_DEFAULT_TG 1.034
#define EPS 0.75

void init_Tangentdata(Tangentdata *tgdata);
void tangent_remove_all_fit_results(Tangentdata *tgdata);
void tangent_fit(GwyDataLine *x, GwyDataLine *y, Tangentdata *tgdata, FDdata *fddata, Instdata *instdata, Area *area, gdouble ux, gdouble uy);
gchar* tangent_export_data(Tangentdata *tgdata, FDdata *fddata, Instdata *instdata, Area *area, enum ExportFormat exportformat);
gboolean tangent_has_results(Tangentdata *tgdata);

#define FIT_TLS TRUE
#define CORRECT FALSE
#define RESCALE FALSE
#define DELT  TRUE

#endif

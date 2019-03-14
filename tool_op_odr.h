#ifndef NOFORTRAN

#ifndef TOOL_OP_ODR_H
#define TOOL_OP_ODR_H

#include "datatypes.h"
#include "file-utils.h"

#define BETA_DEFAULT_OPODR 1.0
#define RANGE_FROM_DEFAULT_OPODR 20.0
#define RANGE_TO_DEFAULT_OPODR 98.0
#define RADIAL_ANGLE_DEFAULT_OPODR 70.3 /* cone equivalent for Berkovich */

void init_OPODRdata(OPODRdata *opodrdata);
void op_odr_remove_all_fit_results(OPODRdata *opodrdata);
gint op_odr_fit(GwyDataLine *x, GwyDataLine *y, OPODRdata *opodrdata, const FDdata *fddata, const Instdata *instdata, const Area *area, const gdouble estimate[3]);
void radial_correction_op(gdouble alpha, gdouble *H, gdouble *Er, gdouble *Eit, gdouble *Aphc, const Instdata *instdata);
gchar* op_odr_export_data(const OPODRdata *opodrdata, const FDdata *fddata, const Instdata *instdata, const Area *area, enum ExportFormat exportformat);
gboolean op_odr_has_results(const OPODRdata *opodr);

#endif

#endif

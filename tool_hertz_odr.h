#ifndef NOFORTRAN

#ifndef TOOL_HERTZ_ODR_H
#define TOOL_HERTZ_ODR_H

#include "datatypes.h"
#include "file-utils.h"

#define N_FIX_HZODR 1.5

void init_HertzODRdata(HertzODRdata *hzodrdata);
gint hertz_odr_fit(GwyDataLine *x, GwyDataLine *y, HertzODRdata *hzodrdata, const FDdata *fddata, const Instdata *instdata, const gdouble estimate[3]);
void radial_correction_hertz(gint mode, gdouble *Er, gdouble *Eit, gdouble *radius, gdouble to, const Instdata *instdata); /* z mode udelat enum */
gchar* hertz_odr_export_data(const HertzODRdata *hzodrdata, const FDdata *fddata, const Instdata *instdata, enum ExportFormat exportformat);
gboolean hertz_odr_has_results(const HertzODRdata *hzodr);
void hertz_odr_remove_all_fit_results(HertzODRdata *hzodrdata);

#endif

#endif

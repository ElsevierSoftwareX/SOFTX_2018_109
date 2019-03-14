#ifndef FIT_UTILS_H
#define FIT_UTILS_H

#include "datatypes.h"
#include <libprocess/gwyprocess.h>

/* LinReg should be moved here after datatypes.h cleanup */

gdouble least_squares_line_get_correction(GwyDataLine *x, gdouble ux);
void total_least_squares_line_fit__(GwyDataLine *x, GwyDataLine *y, LinReg *reg, gdouble ux, gdouble uy, gboolean correct, gboolean rescale, gboolean delta);
void total_least_squares_line_fit_delta(GwyDataLine *x, GwyDataLine *y, LinReg *reg, gdouble delta);
void total_least_squares_line_fit_nonconst(GwyDataLine *x, GwyDataLine *y, LinReg *reg, GwyDataLine *ux, GwyDataLine *uy);
void total_least_squares_line_fit_const(GwyDataLine *x, GwyDataLine *y, LinReg *reg, gdouble ux, gdouble uy);
void ordinary_least_squares_line_fit__(GwyDataLine *x, GwyDataLine *y, LinReg *reg, gdouble ux, gdouble uy, gboolean correct, gboolean rescale);
void ordinary_least_squares_line_fit(GwyDataLine *x, GwyDataLine *y, LinReg *reg);
void ordinary_least_squares_line_fit_unc_const(GwyDataLine *x, GwyDataLine *y, LinReg *reg, gdouble ux, gdouble uy);
void ordinary_least_squares_line_fit_unc_general(GwyDataLine *x, GwyDataLine *y, LinReg *reg, GwyDataLine *ux, GwyDataLine *uy);

void least_squares_32power_fit(GwyDataLine *x, GwyDataLine *y, gdouble ux, gdouble uy, LinReg *reg) ;
void least_squares_32power_fit_two(GwyDataLine *x, GwyDataLine *y, LinReg *reg);
void nonlinear_least_squares_power_fit(GwyDataLine *x, GwyDataLine *y, gdouble *param);

gint filter_negative_data(GwyDataLine *x, GwyDataLine *y);

gint odr_run_fit(GwyDataLine *data_x, GwyDataLine *data_y, const gint np, gdouble *beta, const gint *ifix, gchar *fnm, gdouble ux, gdouble uy, gdouble *std_beta, gdouble *covariance, gdouble *SSres, gdouble *R2, gdouble *R2adj, gdouble *chi2, gint what_to_fit);

#endif

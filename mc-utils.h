#ifndef MC_UTILS_H
#define MC_UTILS_H

#include <libprocess/gwyprocess.h>

gdouble gaussian_random_number(GRand *rng);
void generate_uncorrelated_noise(GwyDataLine *noise, gdouble unc, GRand *rng);
gboolean sum_data_lines(GwyDataLine *result, GwyDataLine *x, GwyDataLine *y);
gboolean create_histogram(const gdouble *data, gint ndata, gdouble *xhist, gdouble *yhist, gint nstat);

#endif

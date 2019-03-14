#include "mc-utils.h"

#include <math.h>
#include <stdio.h>

#include <libprocess/gwyprocess.h>


gdouble gaussian_random_number(GRand *rng)
{
    gdouble x, y, w;

    /* to avoid four random number generations for a signle gaussian number,
     * we use g_rand_int() instead of g_rand_double() and use the second
     * gaussian number to dd noise to the lower bits */
    do {
        x = -1.0 + 2.0 / 4294967295.0 * g_rand_int(rng);
        y = -1.0 + 2.0 / 4294967295.0 * g_rand_int(rng);
        w = x * x + y * y;
    }
    while (w >= 1.0 || w == 0);

    w = sqrt(-2.0 * log(w) / w);

    return (x + y / 4294967295.0) * w;
}

void generate_uncorrelated_noise(GwyDataLine *noise, gdouble unc, GRand *rng)
{
    // generates dataset of mutually independent Gaussian random data
    gint i;
    gdouble *gt;
    gdouble z;

    if (unc < 1e-6) {
        gwy_data_line_clear(noise);
    }
    else {
        gt = gwy_data_line_get_data(noise);

        for (i = noise->res; i; i--, gt++) {
            z = gaussian_random_number(rng);
            *gt = unc * z; // TODO check if pointer arithmetic is faster
        }
    }

    return ;
}

gboolean sum_data_lines(GwyDataLine *result, GwyDataLine *x, GwyDataLine *y)
{
    gint i;
    gdouble *r, *p, *q;

    p = x->data;
    q = y->data;
    r = result->data;

    if ((result->res != x->res) || (result->res != y->res)) {
        g_printerr("Sum data lines: cannot sum data lines with different sizes!!!\n");
        return FALSE;
    }

    for (i = result->res; i; i--, p++, q++, r++) {
        *r = *p + *q;
    }

    return TRUE;
}

gboolean create_histogram(const gdouble *data, gint ndata, gdouble *xhist, gdouble *yhist, gint nstat)
{
    gint i, k;
    gdouble min = 1e20, max = -1e20;

    for (i = 0; i < ndata; i++) {
        min = (data[i] < min ? data[i] : min);
        max = (data[i] > max ? data[i] : max);

    }

    for (i = 0; i < nstat; i++) {
        xhist[i] = min + i * (max - min) / nstat;
        yhist[i] = 0;
    }

    for (i = 0; i < ndata; i++) {
        k = (gint)((data[i] - min) / (max - min) * nstat);

        if (G_UNLIKELY(k >= nstat)) {
            k = nstat - 1;
        }

        if (G_UNLIKELY(k < 0)) {
            k = 0;
        }

        yhist[k]++;
    }

    /*
       printf("hist \n");
       for (i=0;i< nstat;i++)
       printf("%g %g \n",xhist[i], yhist[i]);
       */
    return TRUE;
}

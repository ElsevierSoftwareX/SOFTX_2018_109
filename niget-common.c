#include "niget-common.h"

#include <locale.h>
#include <math.h>

#include <libprocess/gwyprocess.h>

/* shared variable for verbosity control */
gboolean verbose;

static gdouble psi(gdouble z);


/* opsane z Octavi verze */
static gdouble psi(gdouble z)
{
    gdouble h;

    h = 1e-9;

    if (z > 0) {
        return (lgamma(z + h) - lgamma(z - h)) / (2 * h);
    }
    else if (z < 0) {
        return (lgamma(1 - z + h) - lgamma(1 - z - h)) / (2 * h) + M_PI / tan(M_PI * (1 - z));
    }
    else {
        return  -1e30;
    }
}

gdouble epsilon(gdouble m)
{
    return m * (1 - 2 * (m - 1) * tgamma(m / 2 / (m - 1)) / tgamma(1.0 / 2.0 / (m - 1)) / sqrt(M_PI));
}

gdouble der_epsilon(gdouble m)
{
    gdouble z, eps, deps;

    z = m / (m - 1) / 2;
    eps = epsilon(m);
    deps = eps / m - m / (m - 1) * (1 - eps / m) - m / 2 / (m - 1) / (m - 1) * (eps / m - 1) * (psi(z) - psi(z - 1. / 2));
    return  deps;
}

gdouble sq(gdouble x)
{
    return x*x;
}

gdouble calc_Hit(gdouble Fmax, gdouble Aphc)
{
    return Fmax / Aphc;
}

gdouble calc_S(gdouble hmax, gdouble Fmax, gdouble hp, gdouble m)
{
    return m * Fmax / (hmax - hp);
}

gdouble calc_hc(gdouble hmax, gdouble Fmax, gdouble S, gdouble eps)
{
    return hmax - eps * Fmax / S;
}

gdouble calc_Er_OP(gdouble Aphc, gdouble S, gdouble beta)
{
    return sqrt(M_PI / Aphc) * S / (2 * beta);
}

gdouble calc_Eit(gdouble Er, gdouble nu, gdouble Ei, gdouble nui)
{
    return (1 - sq(nu)) / (1 / Er - (1 - sq(nui)) / Ei);
}

gdouble calc_Er_Hertz(gdouble a, gdouble R)
{
    return 0.75 * a / sqrt(R);
}

gdouble calc_Er(gdouble Eit, gdouble nu, gdouble Ei, gdouble nui)
{
    return 1.0 / ((1 - sq(nui)) / Ei + (1 - sq(nu)) / Eit);
}

gdouble calc_R_Hertz(gdouble a, gdouble Er)
{
    return 9.0 / 16.0 * sq(a) / sq(Er);
}


GwyDataLine *moving_average(GwyDataLine *x, gint nmove) //TODO prepsat hezceji, rychleji, obecnejsi vahy, doplnit okraje
{
    GwyDataLine *data_line;
    gint i, j, s, n;

    if (nmove == 1 || nmove > (x->res / 2)) { /* temporary fix to avoid out of mem access */
        data_line = gwy_data_line_duplicate(x);
    }
    else {
        s = (nmove - 1) / 2; // nmove must be odd
        n = x->res;
        data_line = gwy_data_line_new_alike(x, TRUE); /* LEAK */

        for (i = s; i < n - s; i++) {
            for (j = -s; j <= s ; j++) {
                data_line->data[i] += x->data[i + j] / nmove;
            }
        }

        for (i = 0; i < s; i++) {
            /*			for (j = MAX(-s,-i); j <= s ; j++){
            			data_line->data[i] += x->data[i+j]/(nmove-s+i);
            			}
            			*/
            data_line->data[i] = x->data[i];
        }

        for (i = n - s; i < n; i++) {
            /*			for (j = -s; j <MIN( s,n-i) ; j++){
            			data_line->data[i] += x->data[i+j]/(nmove-s+n-1-i);
            			}
            			*/
            data_line->data[i] = x->data[i];
        }

    }

    return data_line;
}

/* not published in .h */
GwyDataLine *moving_average_rs(GwyDataLine *x, gint nmove) //TODO prepsat hezceji, rychleji, obecnejsi vahy, doplnit okraje
{
    GwyDataLine *data_line;
    gint i, j, l, u, w, s, n;

    if (nmove == 1) {
        data_line = gwy_data_line_duplicate(x);
    }
    else {
        s = (nmove - 1) / 2; // nmove must be odd
        n = x->res;
        data_line = gwy_data_line_new_alike(x, TRUE);

        for (i = 0; i < n; i++) {
            l = MAX(0,   i - s);
            u = MIN(n - 1, i + s);
            w = u - l + 1;

            for (j = l; j <= u; j++) {
                data_line->data[i] += x->data[j];
            }

            data_line->data[i] /= w;
        }
    }

    return data_line;
}

void range_to_indices(gdouble from, gdouble to, GwyDataLine *x, gboolean descending, gint *istart, gint *iend, gint *ndata)
{
    gint i = 0;

    *ndata = 0;

    if (descending) {
        while (i < x->res && x->data[i] > to) {
            i++;
        }

        *istart = i;

        while (i < x->res && x->data[i] > from) {
            i++;
        }

        *iend = i;
    }
    else {
        while (i < x->res && x->data[i] < from) {
            i++;
        }

        *istart = i;

        while (i < x->res && x->data[i] < to) {
            i++;
        }

        if (i == x->res) {
            i--;
        }

        *iend = i;
    }

    *ndata = *iend - *istart;
}

void range_to_indices_array(gdouble from, gdouble to, const gdouble *x, gint l, gboolean descending, gint *istart, gint *iend, gint *ndata)
{
    gint i;

    if (l > 0) {
        if (!descending) {
            i = 0;

            while (i < (l - 1) && x[i] < from) {
                i++;
            }

            *istart = i;

            i = l - 1;

            while (i > *istart && x[i] > to) {
                i--;
            }

            *iend = i;
        }
        else {
            i = l - 1;

            while (i > 0 && x[i] < from) {
                i--;
            }

            *iend = i;

            i = 0;

            while (i < *iend && x[i] > to) {
                i++;
            }

            *istart = i;
        }

        if ((*istart == *iend) && ((x[*istart] < from) || (x[*istart] > to))) {
            *istart = -1;
            *iend = -1;
            *ndata = 0;
        }
        else {
            *ndata = *iend - *istart + 1;
        }
    }
    else { /* empty array */
        *istart = -1;
        *iend = -1;
        *ndata = 0;
    }
}

void range_to_indices_dataline(gdouble from, gdouble to, GwyDataLine *x, gboolean descending, gint *istart, gint *iend, gint *ndata)
{
    range_to_indices_array(from, to, x->data, gwy_data_line_get_res(x), descending, istart, iend, ndata);
}

void localize_decimal_points(gchar *str)
{
    gchar *p;
    struct lconv *locale;

    if (str) {
        locale = localeconv();

        if (locale->decimal_point) {
            p = str;

            while (*p) {
                if ((*p == '.') || (*p == ',')) {
                    *p = *(locale->decimal_point);
                }

                p++;
            }
        }
    }
}

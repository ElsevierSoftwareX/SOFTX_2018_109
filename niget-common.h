#ifndef NIGET_COMMON_H
#define NIGET_COMMON_H

#include <libprocess/gwyprocess.h>

/* temporary fix for MSVC */
#ifndef VERBOSE
#define VERBOSE FALSE
#endif

#ifndef SAVEMC
#define SAVEMC FALSE
#endif

#ifndef SAVEMC2
#define SAVEMC2 FALSE
#endif


#define NIGET_SHORTNAME "Niget"
#define NIGET_LONGNAME "NanoIndentation General Evaluation Tool"
#define NIGET_AUTHORS "Anna Charvátová Campbell, Radek Šlesinger, Petr Grolich"
#define NIGET_CMI "Czech Metrology Institute"
#define NIGET_LICENSE "Licensed under the GNU General Public License version 2"
#define NIGET_URL "http://nanometrologie.cz/niget"

#define NIGET_VERSION "0.5" /* change before release! */


#define N_MC 10000 // number of Monte Carlo iterations

#ifndef NOFORTRAN
enum {
    POWER_3PARAM = 0,
    HERTZ_RADIAL,
    POWER_2PARAM,
    LINEAR,
};

enum {
    CONTACT = 0,
    OP,
    OPODR,
    TANGENT,
    HERTZ,
    HERTZODR,
    SLOPES,
    DIFF,
    APOP,
    WORK,
    CREEP,
    RELAX,
    PPOP,
};
#else
enum {
    CONTACT = 0,
    OP,
    TANGENT,
    HERTZ,
    DIFF,
    APOP,
    WORK,
    CREEP,
    RELAX,
    PPOP,
};
#endif

#define DEPTH "%.3f"
#define FORCE "%.3f"
#define HARD "%.2f"
#define MODUL "%.4f"
#define NUMBER "%.3f"
#define AREA  "%.0f"
#define SLOPE "%.5g"
#define EPWORK  "%.2f"
#define DELTA   "= σ<sub>F</sub><sup>2</sup>/σ<sub>h</sub><sup>2</sup> = %.03g"

extern gboolean verbose;

gdouble epsilon(gdouble m);
gdouble der_epsilon(gdouble m);
gdouble sq(gdouble x);
gdouble calc_Hit(gdouble Fmax, gdouble Aphc);
gdouble calc_S(gdouble hmax, gdouble Fmax, gdouble hp, gdouble m);
gdouble calc_hc(gdouble hmax, gdouble Fmax, gdouble S, gdouble eps);
gdouble calc_Er_OP(gdouble Aphc, gdouble S, gdouble beta);
gdouble calc_Er_Hertz(gdouble a, gdouble R);
gdouble calc_Eit(gdouble Er, gdouble nu, gdouble Ei, gdouble nui);
gdouble calc_Er(gdouble Eit, gdouble nu, gdouble Ei, gdouble nui);
gdouble calc_R_Hertz(gdouble a, gdouble Er);

GwyDataLine *moving_average(GwyDataLine *x, gint nmove);

void range_to_indices(gdouble from, gdouble to, GwyDataLine *x, gboolean descending, gint *istart, gint *iend, gint *ndata);
void range_to_indices_array(gdouble from, gdouble to, const gdouble *x, gint l, gboolean descending, gint *istart, gint *iend, gint *ndata);
void range_to_indices_dataline(gdouble from, gdouble to, GwyDataLine *x, gboolean descending, gint *istart, gint *iend, gint *ndata);

void localize_decimal_points(gchar *str);

#endif

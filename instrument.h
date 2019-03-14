#ifndef INSTRUMENT_H
#define INSTRUMENT_H

#include <libprocess/gwyprocess.h>

enum {
    AREA_TRIVIAL = 0,
    AREA_POLYNOM,
    AREA_DATA,
};

typedef enum {
    AREA_FILE_OK = 0,
    AREA_FILE_COULD_NOT_OPEN,
    AREA_FILE_NOT_ENOUGH_ENTRIES,
    AREA_FILE_FORMAT_MISMATCH,
    AREA_FILE_POLYNOM,
    AREA_FILE_TIP
} AreaFileReadStatus;

typedef struct {
    gchar *filename;
    GwyDataLine *xdata; /* in nm */
    GwyDataLine *ydata; /* in nm^2 */
    gint mode;
    gint npoly; // = 6 inicializovat a prozatim na to nesahat
    gdouble *polycoef;
    gdouble *polypower;
    gchar **polylabel;
    gdouble **covpolycoef;
    gdouble xmax;
    gchar *coeff_filename;
} Area;

typedef struct {
    gdouble nu;
    gdouble nui;
    gdouble Ei; /* Ei in GPa */
    gdouble unu;
    gdouble unui;
    gdouble uEi; /* u(Ei) in GPa */
    gdouble sigma_h, sigma_F, delta;
} Instdata;

void init_Area(Area *area);
void init_Instdata(Instdata *instdata);

gdouble polynom_eval(const Area *area, gdouble h);
gdouble eval_area(const Area *area, gdouble h);
gdouble eval_area_derivative(const Area *area, gdouble h);
gdouble eval_area_uncertainty(const Area *area, gdouble h);
gchar *write_area_label(const Area *area);
gchar *write_area_file(const Area *area);

AreaFileReadStatus read_area_hysitron_file(const gchar *fnm, Area *area);
AreaFileReadStatus read_area_csm_file(const gchar *fnm, Area *area);
AreaFileReadStatus read_area_niget_file(const gchar *fnm, Area *area);
AreaFileReadStatus read_area_dat_file(const gchar *fnm, Area *area);

#endif

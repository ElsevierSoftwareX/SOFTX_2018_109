#include "instrument.h"

#include "file-utils.h"
#include "niget-common.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>


#define NU_DEFAULT 0.3
#define NUI_DEFAULT 0.07
#define EI_DEFAULT 1141e9
#define UNC_NU_DEFAULT 0.1
#define UNC_NUI_DEFAULT 0.1
#define UNC_EI_DEFAULT 100e9
#define SIGMA_H_DEFAULT 0.1
#define SIGMA_F_DEFAULT 0.0001


static void fit_area_polynom(const Area *area, gdouble *a, gdouble *b);
static gdouble interpolate_area_function(const Area *area, gdouble h);
static gdouble interpolate_area_function_derivative(const Area *area, gdouble h);
static gdouble polynom_eval_derivative(const Area *area, gdouble h);
static gdouble polynom_eval_uncertainty(const Area *area, gdouble h);

static void fit_area_polynom(const Area *area, gdouble *a, gdouble *b)
{
    gdouble sx4, sx3, sx2, syx2, syx;
    gint n;
    gint i;
    gdouble *x, *y;

    x = area->xdata->data;
    y = area->ydata->data;
    n = area->xdata->res;

    sx4 = 0;
    sx3 = 0;
    sx2 = 0;
    syx2 = 0;
    syx = 0;

    for (i = 0; i < n; i++) {
        sx4 += x[i] * x[i] * x[i] * x[i];
        sx3 += x[i] * x[i] * x[i];
        sx2 += x[i] * x[i];
        syx2 += y[i] * x[i] * x[i];
        syx += y[i] * x[i];
    }

    *a = (syx2 * sx2 - syx * sx3) / (sx4 * sx2 - sx3 * sx3);
    *b = (syx * sx4 - syx2 * sx3) / (sx4 * sx2 - sx3 * sx3);
}

static gdouble interpolate_area_function(const Area *area, gdouble h)
{
    gint i, n;
    gboolean asc;
    gdouble *x;
    gdouble val;
    gdouble rest, a, b;

    x = area->xdata->data;
    n = area->xdata->res;

    asc = (x[0] < x[n - 1] ? TRUE : FALSE);

    if (h < 0) {
        return 0;
    }

    i = 0;

    if (asc) {
        while (i < n - 1 && x[i] < h) {
            i++;
        }
    }
    else {
        while (i < n - 1 && x[i] > h) {
            i++;
        }
    }

    if (h > x[n - 1]) {
        fit_area_polynom(area, &a, &b);
        val = a * h * h + b * h;
    }
    else {
        rest = (h - x[i - 1]) / (x[i] - x[i - 1]);
        val = gwy_data_line_get_dval(area->ydata, i - 1 + rest + 0.5, GWY_INTERPOLATION_BSPLINE);
    }

    return val;
}

static gdouble interpolate_area_function_derivative(const Area *area, gdouble h)
{
    gint i, n;
    gboolean asc;
    gdouble *x;
    gdouble dval;
    GwyDataLine *darea;
    gdouble rest, a, b;

    x = area->xdata->data;
    n = area->xdata->res;

    darea = gwy_data_line_new_alike(area->ydata, FALSE);

    for (i = 0; i < darea->res; i++) {
        darea->data[i] = gwy_data_line_get_der(area->ydata, i);
    }

    asc = (x[0] < x[n - 1] ? TRUE : FALSE);

    if (h < 0) {
        return 0;
    }

    i = 0;

    if (asc) {
        while (i < n - 1 && x[i] < h) {
            i++;
        }
    }
    else {
        while (i < n - 1 && x[i] > h) {
            i++;
        }
    }

    if (h > x[n - 1]) {
        fit_area_polynom(area, &a, &b);
        dval = a * 2 * h + b;
    }
    else {
        rest = (h - x[i - 1]) / (x[i] - x[i - 1]);
        dval = gwy_data_line_get_dval(darea, i - 1 + rest + 0.5, GWY_INTERPOLATION_BSPLINE);
    }

    g_object_unref(darea);

    return dval;
}

static gdouble polynom_eval_derivative(const Area *area, gdouble h)
{
    gint i;
    gdouble val = 0;

    for (i = 0; i < area->npoly; i++) {
        val += area->polycoef[i] * area->polypower[i] * pow(h, area->polypower[i] - 1);
    }

    return val;
}

void init_Area(Area *area)
{
    gint i;

    area->filename = NULL;

    area->xdata = NULL;
    area->ydata = NULL;
    area->xmax = 0;

    area->mode = AREA_TRIVIAL;

    area->npoly = 6;

    area->polycoef = (gdouble *)g_malloc0(area->npoly * sizeof(gdouble)); // neuvolni se nikdy, valgrind error, zatim nechat
    area->polycoef[0] = 24.5;

    area->polypower = (gdouble *)g_malloc(area->npoly * sizeof(gdouble)); // neuvolni se nikdy, valgrind error, zatim nechat

    for (i = 0;  i < area->npoly; i++) {
        area->polypower[i] = pow(2, 1 - i);
    }

    area->polylabel = (gchar **)g_malloc(area->npoly * sizeof(gchar *)); // neuvolni se nikdy, valgrind error, zatim nechat
    area->polylabel[0] = "2";
    area->polylabel[1] = "1";
    area->polylabel[2] = "1/2";
    area->polylabel[3] = "1/4";
    area->polylabel[4] = "1/8";
    area->polylabel[5] = "1/16";
    /* area->polylabel[6] = "1/32"; */

    area->coeff_filename = NULL;

    /* covariance matrix */
    area ->covpolycoef = (gdouble **) g_malloc(area->npoly * sizeof(gdouble *));

    for (i = 0;  i < area->npoly; i++) {
        area->covpolycoef[i] = (gdouble *)g_malloc0(area->npoly * sizeof(gdouble));
    }
}

void init_Instdata(Instdata *instdata)
{
    instdata->Ei = EI_DEFAULT;
    instdata->nui = NUI_DEFAULT;
    instdata->nu = NU_DEFAULT;
    instdata->uEi = UNC_EI_DEFAULT;
    instdata->unui = UNC_NUI_DEFAULT;
    instdata->unu = UNC_NU_DEFAULT;
    instdata->sigma_h = SIGMA_H_DEFAULT;
    instdata->sigma_F = SIGMA_F_DEFAULT;
    instdata->delta = SIGMA_F_DEFAULT * SIGMA_F_DEFAULT / SIGMA_H_DEFAULT / SIGMA_H_DEFAULT;
}

/* not static as it is presently used in the area dialog */
gdouble polynom_eval(const Area *area, gdouble h)
{
    gint i;
    gdouble val = 0;

    for (i = 0; i < area->npoly; i++) {
        val += area->polycoef[i] * pow(h, area->polypower[i]);
    }

    return val;
}

gdouble polynom_eval_uncertainty(const Area *area, gdouble h)
{
    gint i, j;
    gdouble val = 0;
    gdouble *der;

    der = (gdouble *)g_malloc0(area->npoly * sizeof(gdouble));

    for (i = 0; i < area->npoly; i++) {
        der[i] = pow(h, area->polypower[i]);
    }

    for (i = 0; i < area->npoly; i++) {
        for (j = 0; j < area->npoly; j++) {
            val += der[i] * der[j] * area->covpolycoef[i][j];
        }
    }

    g_free(der);

    return val;
}

gdouble eval_area(const Area *area, gdouble h)
{
    gdouble Aphc;

    if (h < 0) {
        return 0;
    }

    switch (area->mode) {
    case (AREA_TRIVIAL):
        Aphc = 24.5 * h * h;
        break;

    case (AREA_POLYNOM):
        Aphc = polynom_eval(area, h);
        break;

    case (AREA_DATA):
        Aphc = interpolate_area_function(area, h);
        break;
    }

    return Aphc;
}

gdouble eval_area_derivative(const Area *area, gdouble h)
{
    gdouble dAphc;

    if (h < 0) {
        return 0;
    }

    switch (area->mode) {
    case (AREA_TRIVIAL):
        dAphc = 24.5 * 2 * h;
        break;

    case (AREA_POLYNOM):
        dAphc = polynom_eval_derivative(area, h);
        break;

    case (AREA_DATA):
        dAphc = interpolate_area_function_derivative(area, h);
        break;
    }

    return dAphc;
}

gdouble eval_area_uncertainty(const Area *area, gdouble h)
{
    gdouble uncA;

    if (h < 0) {
        return 0;
    }

    switch (area->mode) {
    case (AREA_TRIVIAL):
        uncA = 0;
        break;

    case (AREA_POLYNOM):
        uncA = polynom_eval_uncertainty(area, h);
        break;

    case (AREA_DATA):
        uncA = 0 ; //TODO not yet implemented
        break;
    }

    return uncA;
}
/* zkonsolidovat s write_area_file */
gchar *write_area_label(const Area *area)
{
    gchar **strpoly;
    gdouble eps = 0.001;
    gint nsize = 64;
    gchar *str;
    gint i;

    if (area->mode == AREA_DATA) {
        str = g_strdup(area->filename);
    }
    else {
        strpoly = (gchar **)g_malloc((area->npoly + 1) * sizeof(gchar *));

        for (i = 0; i < area->npoly; i++) {
            strpoly[i] = (gchar *)g_malloc(nsize * sizeof(gchar));
        }

        strpoly[area->npoly] = NULL;

        for (i = 0; i < area->npoly; i++) {
            if (fabs(area->polycoef[i]) > eps) {
                if (i == 0) {
                    g_snprintf(strpoly[i], nsize, "%.3gx<sup>%s</sup>", area->polycoef[i], area->polylabel[i]);
                }
                else if (area->polycoef[i] > eps) {
                    g_snprintf(strpoly[i], nsize, " + %.3gx<sup>%s</sup>", area->polycoef[i], area->polylabel[i]);
                }
                else {
                    g_snprintf(strpoly[i], nsize, " - %.3gx<sup>%s</sup>", fabs(area->polycoef[i]), area->polylabel[i]);
                }
            }
            else {
                strpoly[i] = g_strdup("");
            }
        }

        str = g_strjoinv("", strpoly);
        g_strfreev(strpoly);
    }

    return str;
}

/* zkonsolidovat s write_area_label */
gchar *write_area_file(const Area *area)
{
    gchar **strpoly;
    gdouble eps = 0.001;
    gint nsize = 64;
    gchar *str;
    gint i;
    gchar buffer[1000];

    if (area->mode == AREA_DATA) {
        g_snprintf(buffer, sizeof(buffer), "#Raw data file: %s \n", area->filename);
        str = g_strdup(buffer);
    }
    else {
        strpoly = (gchar **)g_malloc((area->npoly + 1) * sizeof(gchar *));

        for (i = 0; i < area->npoly; i++) {
            strpoly[i] = (gchar *)g_malloc(nsize * sizeof(gchar));
        }

        strpoly[area->npoly] = NULL;

        for (i = 0; i < area->npoly; i++) {
            if (fabs(area->polycoef[i]) > eps) {
                if (i == 0) {
                    g_snprintf(strpoly[i], nsize, "%gx^%s", area->polycoef[i], area->polylabel[i]);
                }
                else if (area->polycoef[i] > eps) {
                    g_snprintf(strpoly[i], nsize, " + %gx^%s", area->polycoef[i], area->polylabel[i]);
                }
                else {
                    g_snprintf(strpoly[i], nsize, " - %gx^%s", fabs(area->polycoef[i]), area->polylabel[i]);
                }
            }
            else {
                strpoly[i] = g_strdup("");
            }
        }

        str = g_strjoinv("", strpoly);
        g_strfreev(strpoly);
        g_snprintf(buffer, sizeof(buffer), "#Polynomial area function:  %s \n", str);
        g_free(str);
        str = g_strdup(buffer);
    }

    return str;
}

AreaFileReadStatus read_area_niget_file(const gchar *fnm, Area *area)
{
    FILE *infile;
    gchar *line = NULL;
    size_t size;
    gint i, j, n ;
    gchar *pch;
    gdouble tmp;

    infile = fopen(fnm, "r");

    if (!infile) {
        return AREA_FILE_COULD_NOT_OPEN;
    }

    n = get_number_lines(infile);

    if (n < 5) {
        return AREA_FILE_NOT_ENOUGH_ENTRIES;
    }

    getline(&line, &size, infile); // header

    area->npoly = (n - 3);
    area->polycoef[0] = 24.5;

    if (verbose) {
        g_print("Read area NIGET:\n");
        g_print("  n %d \n", n);
        g_print("  npoly %d \n", area->npoly);
    }

    getline(&line, &size, infile); // coefficients
    localize_decimal_points(line);

    pch = strtok(line, " ");

    for (i = 0; i < area->npoly - 1; i++) {
        if (sscanf(pch, "%lf", &(area->polycoef[i + 1])) != 1) {
            return AREA_FILE_FORMAT_MISMATCH;
        }

        pch = strtok(NULL, " ");
    }

    // if the line was longer, there is a problem with the format
    if (sscanf(pch, "%lf", &tmp) == 1) {
        return AREA_FILE_FORMAT_MISMATCH;
    }


    getline(&line, &size, infile); //empty line
    getline(&line, &size, infile); // header

    //covariance matrix
    for (i = 0; i < area->npoly - 1; i++) {
        getline(&line, &size, infile);
        pch = strtok(line, " ");

        for (j = 0; j < area->npoly - 1; j++) {
            if (sscanf(pch, "%lf", &(area->covpolycoef[i + 1][j + 1])) != 1) {
                return AREA_FILE_FORMAT_MISMATCH;
            }

            pch = strtok(NULL, " ");
        }

        if (sscanf(pch, "%lf", &tmp) == 1) {
            return AREA_FILE_FORMAT_MISMATCH;
        }
    }

    /*
    if (verbose) {
    	g_print("polynom \n");
    	for (i = 0; i< area->npoly; i++){
    		g_print(" %g ", area->polycoef[i]);
    	}
    	g_print("\n");
    	g_print("covariance matrix \n");
    	for (i = 0; i< area->npoly; i++){
    		for (j = 0; j< area->npoly; j++){
    			g_print(" %g ", area->covpolycoef[i][j]);
    		}
    		g_print("\n");
    	}
    }
    */

    fclose(infile);
    g_free(line);

    area->mode = AREA_POLYNOM;
    return AREA_FILE_OK;
}

AreaFileReadStatus read_area_hysitron_file(const gchar *fnm, Area *area)
{
    FILE *infile;
    gchar *line = NULL;
    size_t size;
    gint i, n ;

    infile = fopen(fnm, "r");

    if (!infile) {
        return AREA_FILE_COULD_NOT_OPEN;
    }

    n = get_number_lines(infile);

    if (n < 14) {
        return AREA_FILE_NOT_ENOUGH_ENTRIES;
    }

    getline(&line, &size, infile); // header

    area->npoly = (n - 22 - 4) / 4;

    if (verbose) {
        g_print("Read area Hysitron:\n");
        g_print("  n %d \n", n);
        g_print("  npoly %d \n", area->npoly);
    }

    for (i = 0; i < area->npoly; i++) {
        getline(&line, &size, infile); // name of coefficient
        getline(&line, &size, infile); // coefficient
        localize_decimal_points(line);

        if (sscanf(line, "%lf", &(area->polycoef[i])) != 1) {
            return AREA_FILE_FORMAT_MISMATCH;
        }
    }

    fclose(infile);
    g_free(line);

    area->mode = AREA_POLYNOM;
    return AREA_FILE_OK;
}

AreaFileReadStatus read_area_csm_file(const gchar *fnm, Area *area)
{
    FILE *infile;
    gchar *line = NULL;
    size_t size;
    gint i, n ;
    gchar **c;
    gdouble p;

    infile = fopen(fnm, "r");

    if (!infile) {
        return AREA_FILE_COULD_NOT_OPEN;
    }

    n = get_number_lines(infile);

    if (n < 14) {
        return AREA_FILE_NOT_ENOUGH_ENTRIES;
    }

    // TODO add possibility to read h,A data pairs

    /* search for keyword FIT */
    getline(&line, &size, infile);

    while (strstr(line, "FIT") == NULL) {
        getline(&line, &size, infile);
    }

    getline(&line, &size, infile); // degree of polynome

    if (strstr(line, "POLYNOME DEGREE") == NULL) {
        return AREA_FILE_FORMAT_MISMATCH;
    }

    c = g_strsplit(line, "=", -1);

    if (verbose) {
        g_print("Read area CSM:\n");
#ifndef G_OS_WIN32
        g_print("  size c len c %zu %zu \n", sizeof(c), sizeof(c[0]));
#else
        g_print("  size c len c %Iu %Iu \n", sizeof(c), sizeof(c[0]));
#endif
    }

    area->npoly = atoi(c[1]) + 1;
    g_free(c);

    if (area->npoly > 6) {
        if (verbose) {
            g_print("  Currently only OP polynom up to order 5 is supported. \n");
        }

        return AREA_FILE_POLYNOM;
    }

    getline(&line, &size, infile); // OP polynome

    if (strstr(line, "OLIVER") == NULL) {
        if (verbose) {
            g_print("  Currently only OP polynome up to order 5 is supported. \n");
        }

        return AREA_FILE_POLYNOM;
    }

    getline(&line, &size, infile); // 24.5 coefficient

    if (strstr(line, "24.5") == NULL) {
        if (verbose) {
            g_print("  Currently only Berkovich tips with coefficient 24.5 are implemented. \n");
        }

        return AREA_FILE_TIP;
    }

    for (i = 1; i <  area->npoly; i++) {
        if (getline(&line, &size, infile) > -1) {
            c = g_strsplit(line, "=", -1);
            p = g_strtod(c[1], NULL);

            if (verbose) {
                g_print("  c %s p %g \n", c[1], p);
            }

            g_strfreev(c);
            area->polycoef[i] = p * pow(10, 9 * (pow(2, i) - 1) / pow(2, i - 1));

            if (verbose) {
                g_print("  %d  %g %g \n", i, p, area->polycoef[i]);
            }
        }
        else {
            g_printerr("  Number of coefficients does not correspond to order of polynom.\n");
            return AREA_FILE_FORMAT_MISMATCH;
        }
    }

    fclose(infile);
    g_free(line);

    area->mode = AREA_POLYNOM;
    return AREA_FILE_OK;
}

AreaFileReadStatus read_area_dat_file(const gchar *fnm, Area *area)
{
    FILE *infile;
    gchar *line = NULL;
    size_t size;
    gint i, n, nn ;
    gdouble  x, y;

    infile = fopen(fnm, "r");

    if (!infile) {
        return AREA_FILE_COULD_NOT_OPEN;
    }

    n = get_number_lines(infile);

    if (n < 10) {
        return AREA_FILE_NOT_ENOUGH_ENTRIES;
    }

    if (!area->xdata) {
        area->xdata = gwy_data_line_new(n + 1, n + 1, FALSE);
    }
    else {
        gwy_data_line_resample(area->xdata, n + 1, GWY_INTERPOLATION_NONE);
    }

    if (!area->ydata) {
        area->ydata = gwy_data_line_new(n + 1, n + 1, FALSE);
    }
    else {
        gwy_data_line_resample(area->ydata, n + 1, GWY_INTERPOLATION_NONE);
    }

    area->xdata->data[0] = 0; // TODO je to nutne a rozumne?
    area->ydata->data[0] = 0;

    for (i = 0; i < n; i++) {
        getline(&line, &size, infile);
        localize_decimal_points(line);
        nn = sscanf(line, "%lf %lf", &x, &y);

        if (nn < 2) {
            return AREA_FILE_FORMAT_MISMATCH;
        }

        area->xdata->data[i + 1] = x;
        area->ydata->data[i + 1] = y;
    }

    fclose(infile);
    g_free(line);

    area->mode = AREA_DATA;
    area->xmax = gwy_data_line_get_max(area->xdata);
    return AREA_FILE_OK;
}

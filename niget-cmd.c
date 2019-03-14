#include "fddata.h"
#include "file-utils.h"
#include "fit-utils.h"
#include "niget-common.h"
#include "instrument.h"

#ifdef G_OS_WIN32
#include "getline.h"
#endif

#include "tool_op.h"
#include "tool_op_unc.h"
#ifndef NOFORTRAN
#include "tool_op_odr.h"
#include "tool_op_odr_unc.h"
#include "tool_twoslopes.h"
#include "tool_twoslopes_unc.h"
#include "tool_hertz_odr.h"
#include "tool_hertz_odr_unc.h"
#include "tool_stiffness.h"
#include "tool_stiffness_unc.h"
#endif
#include "tool_tangent.h"
#include "tool_tangent_unc.h"
#include "tool_epwork.h"
#include "tool_epwork_unc.h"
#include "tool_hertz.h"
#include "tool_hertz_unc.h"
#include "tool_Ph2.h"
#include "tool_apopins.h"
#include "settings.h"

#include <locale.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define appname "Niget"

#define CONTACTMAX 3

void  print_help()
{
    printf(" -h    show this help\n");
    printf(" -f FILENAME    input file in NIGET format\n");
    printf(" -a AREA FILENAME    input area file in Hysitron (.ara) or CSM (.ind)  format\n");
    printf(" -U SETTINGS.FILENAME     classical calculation of uncertainties\n");
    printf(" -M NUMBER.OF.MONTE.CARLO.CYCLES     Monte Carlo calculation of uncertainties\n");
    printf(" -D      save loading curve in gnuplot friendly format\n");
    printf("\n");
    printf("Tools:\n");
#ifndef NOFORTRAN
    printf(" -d SETTINGS.FILENAME     Oliver Pharr (ODR) method\n");
    printf(" -s SETTINGS.FILENAME     Two Slopes method\n");
#endif
    printf(" -z SETTINGS.FILENAME     Hertz method\n");
#ifndef NOFORTRAN
    printf(" -t SETTINGS.FILENAME     Hertz (ODR) method\n");
    printf(" -e SETTINGS.FILENAME     Stiffness (ODR) \n");
#endif
    printf(" -p SETTINGS.FILENAME     popin detection\n");
    printf(" -P moving average width  P-h2 curve\n");
    printf(" -w                       elastic-plastic work\n");
}

gboolean read_area_file_cmd(gchar *areafnm, Area *area)
{
    AreaFileReadStatus status;

    status = AREA_FILE_COULD_NOT_OPEN;

    if (strstr(areafnm, ".ara") != NULL) {
        //		printf("hysitron format \n");
        status = read_area_hysitron_file(areafnm, area);
    }
    else if (strstr(areafnm, ".ind") != NULL) {
        //	printf("csm format \n");
        status = read_area_csm_file(areafnm, area);
    }
    else if (strstr(areafnm, ".ani") != NULL) {
        //	printf("cmi format \n");
        status = read_area_niget_file(areafnm, area);
    }
    else if (strstr(areafnm, ".dat") != NULL) {
        printf("raw data format \n");
        status = read_area_dat_file(areafnm, area);
        area->filename = g_strdup(areafnm);
    }
    else {
        printf(" Unknown file extension %s. \n", areafnm);
    }

    //printf("status %d \n", status);

    switch (status) {
    case AREA_FILE_COULD_NOT_OPEN:
        printf("Could not open file %s.\n", areafnm);
        return FALSE;
        break;

    case AREA_FILE_NOT_ENOUGH_ENTRIES:
        printf("Bad file format in file %s.\n", areafnm);
        return FALSE;
        break;

    case AREA_FILE_FORMAT_MISMATCH:
        printf("Bad file format in file %s.\n", areafnm);
        return FALSE;
        break;

    case AREA_FILE_POLYNOM:
        return FALSE;
        break;

    case AREA_FILE_TIP:
        return FALSE;
        break;

    case AREA_FILE_OK:
        return TRUE;
        break;

    default:
        return FALSE;
        break;
    }
}

gboolean read_settings_apopin(gchar *settings, Apopindata *apopdata)
{
    FILE *infile;
    gchar *line = NULL;
    size_t size;
    gboolean ok = TRUE;
    gint nlines;

    infile = fopen(settings, "r");

    if (!infile) {
        return FALSE;
    }

    nlines = get_number_lines(infile) ;

    if (nlines < 6) {
        ok = FALSE;
    }

    getline(&line, &size, infile); // header

    getline(&line, &size, infile); // nmove

    if (sscanf(line, "%d", &(apopdata->nmove)) != 1) {
        ok = FALSE;
    }

    getline(&line, &size, infile); //thresh
    localize_decimal_points(line);

    if (sscanf(line, "%lf", &(apopdata->thresh)) != 1) {
        ok = FALSE;
    }

    getline(&line, &size, infile); //thresh2
    localize_decimal_points(line);

    if (sscanf(line, "%lf", &(apopdata->thresh2)) != 1) {
        ok = FALSE;
    }

    getline(&line, &size, infile); //wpop
    localize_decimal_points(line);

    if (sscanf(line, "%d", &(apopdata->wpop)) != 1) {
        ok = FALSE;
    }

    getline(&line, &size, infile); //hpop
    localize_decimal_points(line);

    if (sscanf(line, "%lf", &(apopdata->hpop)) != 1) {
        ok = FALSE;
    }

    fclose(infile);
    return ok;
}

gboolean read_settings_op_odr(gchar *settings, OPODRdata *opodrdata, Instdata *instdata)
{
    FILE *infile;
    gchar *line = NULL;
    size_t size;
    gboolean ok = TRUE;
    gint nlines;

    infile = fopen(settings, "r");

    if (!infile) {
        return FALSE;
    }

    nlines = get_number_lines(infile) ;

    if (nlines < 5) {
        ok = FALSE;
    }

    getline(&line, &size, infile); // header

    getline(&line, &size, infile); // h or F mode

    if (sscanf(line, "%d", &(opodrdata->Finput)) != 1) {
        ok = FALSE;
    }

    if (opodrdata->Finput) {
        getline(&line, &size, infile); //range in F mode
        localize_decimal_points(line);

        if (sscanf(line, "%lf %lf ", &(opodrdata->from_pct_Fmax), &(opodrdata->to_pct_Fmax)) != 2) {
            ok = FALSE;
        }
    }
    else {
        getline(&line, &size, infile); //range in h mode
        localize_decimal_points(line);

        if (sscanf(line, "%lf %lf ",  &(opodrdata->from), &(opodrdata->to)) != 2) {
            ok = FALSE;
        }
    }

    getline(&line, &size, infile); //mechanical parameters: nu, nui, Ei;
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf %lf",  &(instdata->nu), &(instdata->nui), &(instdata->Ei)) != 3) {
        ok = FALSE;
    }

    getline(&line, &size, infile); //estimate of noise
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf ", &(instdata->sigma_h), &(instdata->sigma_F)) != 2) {
        ok = FALSE;
    }
    instdata->delta = instdata->sigma_F * instdata->sigma_F / instdata->sigma_h / instdata->sigma_h;

    // following lines are optional in the settings file, their presence must be checked
    // radial correction
    if ( getline(&line, &size, infile) != -1) { 
	    if (sscanf(line, "%d %lf", &(opodrdata->radial_corr), &(opodrdata->radial_angle)) != 2) {
		    ok = FALSE;
	    }
    }

    //beta
    if ( getline(&line, &size, infile) != -1){
	    if (sscanf(line, "%lf", &(opodrdata->beta)) != 1) {
		    ok = FALSE;
	    }
    }

    fclose(infile);
    return ok;

}

gboolean read_settings_twoslopes(gchar *settings, Slopesdata *sldata, Instdata *instdata)
{
    FILE *infile;
    gchar *line = NULL;
    size_t size;
    gint nlines;

    infile = fopen(settings, "r");

    if (!infile) {
        return FALSE;
    }

    nlines = get_number_lines(infile) ;

    if (nlines < 8) {
        return FALSE;
    }

    getline(&line, &size, infile); // header

    getline(&line, &size, infile); // h or F mode for loading

    if (sscanf(line, "%d", &(sldata->Finputload)) != 1) {
        return FALSE;
    }

    if (sldata->Finputload) {
        getline(&line, &size, infile); //range for loading in F mode
        localize_decimal_points(line);

        if (sscanf(line, "%lf %lf ", &(sldata->loadfrom_pct_Fmax), &(sldata->loadto_pct_Fmax)) != 2) {
            return FALSE;
        }
    }
    else {
        getline(&line, &size, infile); //range for loading in h mode
        localize_decimal_points(line);

        if (sscanf(line, "%lf %lf ",  &(sldata->loadfrom), &(sldata->loadto)) != 2) {
            return FALSE;
        }
    }

    getline(&line, &size, infile); //whether exponent should be fitted

    if (sscanf(line, "%d ",  &(sldata->loadifixb[2])) != 1) {
        return FALSE;
    }

    getline(&line, &size, infile); // h or F mode for unloading

    if (sscanf(line, "%d", &(sldata->Finputunload)) != 1) {
        return FALSE;
    }

    if (sldata->Finputunload) {
        getline(&line, &size, infile); //range for unloading in F mode
        localize_decimal_points(line);

        if (sscanf(line, "%lf %lf ",  &(sldata->unloadfrom_pct_Fmax), &(sldata->unloadto_pct_Fmax)) != 2) {
            return FALSE;
        }
    }
    else {
        getline(&line, &size, infile); //range for unloading in h mode
        localize_decimal_points(line);

        if (sscanf(line, "%lf %lf ",  &(sldata->unloadfrom), &(sldata->unloadto)) != 2) {
            return FALSE;
        }
    }

    getline(&line, &size, infile); //mechanical parameters: nu, nui, Ei;
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf %lf",  &(instdata->nu), &(instdata->nui), &(instdata->Ei)) != 3) {
        return FALSE;
    }

    getline(&line, &size, infile); //estimate of noise
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf ", &(instdata->sigma_h), &(instdata->sigma_F)) != 2) {
        return FALSE;
    }

    instdata->delta = instdata->sigma_F * instdata->sigma_F / instdata->sigma_h / instdata->sigma_h;


    fclose(infile);
    return TRUE;
}

gboolean read_settings_stiffness(gchar *settings, Stiffnessdata *stdata, Instdata *instdata)
{
    FILE *infile;
    gchar *line = NULL;
    size_t size;
    gint nlines;

    infile = fopen(settings, "r");

    if (!infile) {
        return FALSE;
    }

    nlines = get_number_lines(infile) ;

    if (nlines < 4) {
        return FALSE;
    }

    getline(&line, &size, infile); // header

    getline(&line, &size, infile); //range for loading in h mode
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf ",  &(stdata->loadfrom), &(stdata->loadto)) != 2) {
        return FALSE;
    }

    getline(&line, &size, infile); //range for unloading in h mode
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf ",  &(stdata->unloadfrom), &(stdata->unloadto)) != 2) {
        return FALSE;
    }


    getline(&line, &size, infile); //estimate of noise
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf ", &(instdata->sigma_h), &(instdata->sigma_F)) != 2) {
        return FALSE;
    }

    instdata->delta = instdata->sigma_F * instdata->sigma_F / instdata->sigma_h / instdata->sigma_h;


    fclose(infile);
    return TRUE;
}

gboolean read_settings_hertz(gchar *settings, Hertzdata *hzdata, Instdata *instdata)
{
    FILE *infile;
    gchar *line = NULL;
    size_t size;
    gdouble hlp;
    gint nlines;

    infile = fopen(settings, "r");

    if (!infile) {
        return FALSE;
    }

    nlines = get_number_lines(infile) ;

    if (nlines < 5) {
        return FALSE;
    }

    getline(&line, &size, infile); // header

    getline(&line, &size, infile); // R, Er or Eit mode for loading

    if (sscanf(line, "%d", &(hzdata->mode)) != 1) {
        return FALSE;
    }

    getline(&line, &size, infile); // radius for R-mode, contact modulus  for ER-mode, Young's modulus for Eit-mode
    localize_decimal_points(line);

    if (sscanf(line, "%lf",  &hlp) != 1) {
        return FALSE;
    }

    switch (hzdata->mode) {
    case R_MODE:
        hzdata->radius = hlp;
        break;

    case ER_MODE:
        hzdata->Er = hlp;
        break;

    case EIT_MODE:
        hzdata->Eit = hlp;
        break;
    }

    getline(&line, &size, infile); //range for loading
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf ",  &(hzdata->from), &(hzdata->to)) != 2) {
        return FALSE;
    }

    getline(&line, &size, infile); //mechanical parameters: nu, nui, Ei;
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf %lf",  &(instdata->nu), &(instdata->nui), &(instdata->Ei)) != 3) {
        return FALSE;
    }

    fclose(infile);
    return TRUE;
}

gboolean read_settings_hertz_odr(gchar *settings, HertzODRdata *hzodrdata, Instdata *instdata)
{
    FILE *infile;
    gchar *line = NULL;
    size_t size;
    gdouble hlp;
    gint nlines;

    infile = fopen(settings, "r");

    if (!infile) {
        printf("no file %s \n", settings);
        return FALSE;
    }

    nlines = get_number_lines(infile) ;

    if (nlines < 8) {
        printf("only %d lines\n", nlines);
        return FALSE;
    }

    getline(&line, &size, infile); // header

    getline(&line, &size, infile); // R, Er or Eit mode for loading

    if (sscanf(line, "%d", &(hzodrdata->mode)) != 1) {
        printf("bad line %s \n", line);
        return FALSE;
    }

    getline(&line, &size, infile); // radius for R-mode, contact modulus  for ER-mode, Young's modulus for Eit-mode
    localize_decimal_points(line);

    if (sscanf(line, "%lf",  &hlp) != 1) {
        printf("bad line %s \n", line);
        return FALSE;
    }

    switch (hzodrdata->mode) {
    case R_MODE:
        hzodrdata->radius = hlp;
        break;

    case ER_MODE:
        hzodrdata->Er = hlp;
        break;

    case EIT_MODE:
        hzodrdata->Eit = hlp;
        break;
    }

    getline(&line, &size, infile); //range for loading
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf ",  &(hzodrdata->from), &(hzodrdata->to)) != 2) {
        printf("bad line %s \n", line);
        return FALSE;
    }

    getline(&line, &size, infile); //whether exponent should be fitted

    if (sscanf(line, "%d ",  &(hzodrdata->ifixb[2])) != 1) {
        printf("bad line %s \n", line);
        return FALSE;
    }

    getline(&line, &size, infile); //mechanical parameters: nu, nui, Ei;
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf %lf",  &(instdata->nu), &(instdata->nui), &(instdata->Ei)) != 3) {
        printf("bad line %s \n", line);
        return FALSE;
    }

    getline(&line, &size, infile); //estimate of noise
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf ", &(instdata->sigma_h), &(instdata->sigma_F)) != 2) {
        printf("bad line %s \n", line);
        return FALSE;
    }

    instdata->delta = instdata->sigma_F * instdata->sigma_F / instdata->sigma_h / instdata->sigma_h;

    getline(&line, &size, infile); //whether radial correction should be used

    if (sscanf(line, "%d ",  &(hzodrdata->radial_corr)) != 1) {
        printf("bad line %s \n", line);
        return FALSE;
    }

    printf("aux %d \n", hzodrdata->radial_corr);

    if (hzodrdata->radial_corr) {
        hzodrdata->ifixb[2] = 0;
    }

    fclose(infile);
    return TRUE;
}

gboolean read_settings_uncertainties(gchar *settings, Instdata *instdata,
                                     OPODRUncdata *opodrunc, HertzUncdata *hertzunc, HertzODRUncdata *hertzodrunc, SlopesUncdata *slopesunc, StiffnessUncdata *stiffnessunc,
                                     gboolean calc_op_odr, gboolean calc_hertzodr, gboolean calc_slopes, gboolean calc_stiffness)
{
    FILE *infile;
    gchar *line = NULL;
    size_t size;
    gdouble uh, uF;
    infile = fopen(settings, "r");

    if (!infile) {
        return FALSE;
    }

    getline(&line, &size, infile); // header

    getline(&line, &size, infile); //  unc(h), unc(F)
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf ", &uh, &uF) != 2) {
        return FALSE;
    }

    hertzunc->uh = uh;
    hertzunc->uF = uF;


    if (calc_op_odr || calc_hertzodr || calc_slopes || calc_stiffness) {
        if (instdata->sigma_h != uh || instdata->sigma_F != uF) {
            printf("mismatch between settings for method and uncertainties\n");
            return FALSE;
        }
    }
    else {

        instdata->sigma_h = uh;
        instdata->sigma_F = uF;
        instdata->delta = uF * uF / uh / uh;
    }

    getline(&line, &size, infile); //mechanical parameters: unc(nu), unc(nui), unc(Ei);
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf %lf",  &(instdata->unu), &(instdata->unui), &(instdata->uEi)) != 3) {
        return FALSE;
    }

    getline(&line, &size, infile); //for Hertz only: unc(radius), unc(Er), unc(Eit);
    localize_decimal_points(line);

    if (sscanf(line, "%lf %lf %lf",  &(hertzunc->uradius), &(hertzunc->uEr), &(hertzunc->uEit)) != 3) {
        return FALSE;
    }

    if (calc_hertzodr) {
        hertzodrunc->uradius = hertzunc->uradius;
        hertzodrunc->uEr = hertzunc->uEr;
        hertzodrunc->uEit = hertzunc->uEit;
    }

    fclose(infile);
    return TRUE;
}

int main(int argc, char *argv[])
{
    gchar *fnm = NULL;
    gchar *areafnm = NULL;
    gchar fnm1[200];
    gchar *settings = NULL;
    gint istart, iend, ndata;
    GwyDataLine *x, *y;
    gint i, j;
    gint Nmc = 1e1;
    FileSaveStatus filestatus;
    gboolean splitok;
    gchar *buffer;

    FDdata fddata;
    Instdata instdata;
    Area area;

#ifndef NOFORTRAN
    gint fitinfo;
    gint nremove;

    OPODRdata opodrdata;
    OPODRMCdata opodrmc;

    Slopesdata sldata;
    SlopesMCdata slopesmc;

    HertzODRdata hzodrdata;
    HertzODRMCdata hertzodrmc;

    Stiffnessdata stiffdata;
    StiffnessMCdata stiffmc;
#endif

    /* FIXME - taken out from NOFORTRAN block only because of read_settings_uncertainties */
    OPODRUncdata opodrunc;
    SlopesUncdata slunc;
    HertzODRUncdata hzodrunc;
    StiffnessUncdata stiffunc;
    /* */
	    
    Hertzdata hzdata;
    HertzUncdata hzunc;
    HertzMCdata hertzmc;

    Ph2data ph2data;
    Workdata wdata;
    Apopindata apopdata;

    gint c;
    /* gboolean calc_op = FALSE; */
    gboolean calc_op_odr = FALSE;
    /* gboolean calc_tangent = FALSE; */
    gboolean calc_twoslopes = FALSE;
    gboolean calc_hertz = FALSE;
    gboolean calc_hertz_odr = FALSE;
    gboolean calc_stiffness = FALSE;
    gboolean calc_Ph2 = FALSE;
    gboolean calc_work = FALSE;
    gboolean calc_apopin = FALSE;
    gboolean calc_MC = FALSE;
    gboolean calc_unc = FALSE;
    gboolean write_data = FALSE;

#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 36)
    g_type_init();
#endif

    init_FDdata(&fddata);
    init_Instdata(&instdata);
    init_Area(&area);

#ifndef NOFORTRAN
    init_OPODRdata(&opodrdata);
    init_OPODRUncdata(&opodrunc, &instdata);

    init_Slopesdata(&sldata);
    init_SlopesUncdata(&slunc, &instdata);

    init_HertzODRdata(&hzodrdata);
    init_HertzODRUncdata(&hzodrunc, &instdata);

    init_Stiffnessdata(&stiffdata);
    init_StiffnessUncdata(&stiffunc, &instdata);

#endif

    init_Hertzdata(&hzdata);
    init_HertzUncdata(&hzunc);

    init_Ph2data(&ph2data);
    init_Apopindata(&apopdata);
    init_Workdata(&wdata);

    while ((c = getopt(argc, argv, "a:f:d:e:t:s:p:z:P:whDM:U:")) != -1) {
        switch (c) {
        case 'D':
            write_data = TRUE;
            break;

        case 'a':
            areafnm = optarg;

            if (!read_area_file_cmd(areafnm, &area)) {
                printf(" Incorrect area file. \n");
                exit(31);
            }

            break;

        case 'f':
            fnm = optarg;
            break;
            /* case 'o': */
            /* 	calc_op = TRUE; */
            /* 	break; */
#ifndef NOFORTRAN

        case 'd':
            calc_op_odr = TRUE;
            settings = optarg;

            if (! read_settings_op_odr(settings, &opodrdata, &instdata)) {
                printf(" Incorrect settings file for Oliver Pharr ODR method. \n");
                exit(32);
            }

            break;
#endif
            /* case 'g': */
            /* 	calc_tangent = TRUE; */
            /* 	break; */
#ifndef NOFORTRAN

        case 's':
            calc_twoslopes = TRUE;
            settings = optarg;

            if (! read_settings_twoslopes(settings, &sldata, &instdata)) {
                printf(" Incorrect settings file for slopes method. \n");
                exit(33);
            }

            break;
#endif

        case 'z':
            calc_hertz = TRUE;
            settings = optarg;

            if (!read_settings_hertz(settings, &hzdata, &instdata)) {
                printf(" Incorrect settings file for Hertz method. \n");
                exit(34);
            }

            break;
#ifndef NOFORTRAN

        case 't':
            calc_hertz_odr = TRUE;
            settings = optarg;

            if (!read_settings_hertz_odr(settings, &hzodrdata, &instdata)) {
                printf(" Incorrect settings file for Hertz (ODR) method. \n");
                exit(35);
            }

            break;
#endif

        case 'P':
            calc_Ph2 = TRUE;
            ph2data.nmove = atoi(optarg);
            break;

        case 'w':
            calc_work = TRUE;
            break;

        case 'p':
            calc_apopin = TRUE;
            settings = optarg;

            if (!read_settings_apopin(settings, &apopdata)) {
                printf(" Incorrect settings file for popin detection. \n");
                exit(37);
            }

            break;
#ifndef NOFORTRAN

        case 'e':
            calc_stiffness = TRUE;
            settings = optarg;

            if (!read_settings_stiffness(settings, &stiffdata, &instdata)) {
                printf(" Incorrect settings file for stiffness method. \n");
                exit(38);
            }

            break;
#endif

        case 'M':
            calc_MC = TRUE;
            Nmc = atoi(optarg);
            break;

        case 'U':
            calc_unc = TRUE;
            settings = optarg;

            if (!read_settings_uncertainties(settings, &instdata, &opodrunc, &hzunc, &hzodrunc, &slunc, &stiffunc,
                                             calc_op_odr, calc_hertz_odr, calc_twoslopes, calc_stiffness)) { // TODO FIX THIS SOON
                printf(" Incorrect settings file for uncertainties. \n");
                exit(36);
            }

            break;

        case 'h':
        default:
            print_help();
            exit(0);
            break;
        }
    }

    // check filename
    if (fnm == NULL) {
        printf("No input file. \n");
        exit(1);
    }

    // check NOFORTRAN
#ifdef NOFORTRAN

    if (calc_op_odr || calc_twoslopes || calc_hertz_odr) {
        printf("ODR not available without FORTRAN.\n");
        exit(1);
    }

#endif

    if (Nmc <= 0) {
        printf(" Number of Monte Carlo cycles must be positive \n");
        calc_MC = FALSE;
    }


    if (calc_MC && !calc_unc) {
        printf(" Monte Carlo evaluation needs uncertainties \n");
        exit(1);
    }

    load_data_nogui(fnm, DATA_FILE_FORMAT_NIGET, &fddata, &splitok);

    if (!splitok) {
        printf("Problem with splitting the data \n");
        exit(2);
    }

#ifndef NOFORTRAN

    if (calc_op_odr) {
        /* OP ODR */
        if (opodrdata.Finput)
            range_to_indices(opodrdata.from_pct_Fmax * 0.01 * fddata.Fmax, opodrdata.to_pct_Fmax * 0.01 * fddata.Fmax,
                             fddata.Funload, TRUE, &istart, &iend, &ndata);
        else {
            range_to_indices(opodrdata.from, opodrdata.to, fddata.hunload, TRUE, &istart, &iend, &ndata);
        }

        ndata = MIN(ndata, fddata.hunload->res - istart);

        if (ndata < 2) {
            if (opodrdata.Finput) {
                fprintf(stderr,  "Not enough data points for  range [%g, %g] %% from file %s\n", opodrdata.from_pct_Fmax, opodrdata.to_pct_Fmax, fddata.filename);
            }
            else {
                fprintf(stderr,  "Not enough data points for  range [%g, %g] nm from file %s\n", opodrdata.from, opodrdata.to, fddata.filename);
            }
        }
        else {

            x = gwy_data_line_part_extract(fddata.hunload, istart, ndata);
            y = gwy_data_line_part_extract(fddata.Funload, istart, ndata);
            opodrdata.nfitdata = ndata;

            fitinfo = op_odr_fit(x, y, &opodrdata, &fddata, &instdata, &area, NULL);

            if (opodrdata.radial_corr) {
                radial_correction_op(opodrdata.radial_angle, &(opodrdata.Hit), &(opodrdata.Er), &(opodrdata.Eit), &(opodrdata.Aphc), &instdata);
            }

            opodrdata.xfit = gwy_data_line_duplicate(x);
            opodrdata.yfit = gwy_data_line_new_alike(y, FALSE);

            for (i = 0; i < opodrdata.xfit->res; i++) {
                opodrdata.yfit->data[i] =  opodrdata.alpha * pow(opodrdata.xfit->data[i] - opodrdata.hp, opodrdata.m);
            }

            g_object_unref(x);
            g_object_unref(y);

            if (opodrdata.Finput) {
                sprintf(fnm1, "%s__op_odr_range_%.2g_%.2g_%%.dat", fnm, opodrdata.from_pct_Fmax, opodrdata.to_pct_Fmax);
            }
            else {
                sprintf(fnm1, "%s__op_odr_range_%.2g_%.2g_nm.dat", fnm, opodrdata.from, opodrdata.to);
            }

            buffer = op_odr_export_data(&opodrdata, &fddata, &instdata, &area, EXPORT_FORMAT_PLAIN); /* */
            buffer_to_file(fnm1, buffer, &filestatus);
            g_free(buffer);

            if (calc_unc) {
                /* uncertainties */

                /*copy instdata to a local copy */
                opodrunc.instdata = instdata;

                op_odr_propagate_uncertainties_calc(&opodrdata, &opodrunc, &fddata, &area);

                opodrunc.Ec = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
                opodrunc.Hc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
                opodrunc.Ac = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
                opodrunc.hc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));

                for (j = -CONTACTMAX; j <= CONTACTMAX ; j++) {
                    /* use instdata instead of opodrunc.instdata, because they are by definition the same*/
                    op_odr_fit_shift_contact(&opodrdata, &fddata, &instdata, &area, opodrunc.Ec, opodrunc.Hc, opodrunc.Ac, opodrunc.hc, j);
                }

                if (opodrdata.Finput) {
                    sprintf(fnm1, "%s__op_odr_unc_range_%.2g_%.2g_%%.dat", fnm, opodrdata.from_pct_Fmax, opodrdata.to_pct_Fmax);
                }
                else {
                    sprintf(fnm1, "%s__op_odr_unc_range_%.2g_%.2g_nm.dat", fnm, opodrdata.from, opodrdata.to);
                }

                buffer = op_odr_uncertainty_export_data(&opodrdata, &opodrunc, &fddata, &area, EXPORT_FORMAT_PLAIN); /* */
                buffer_to_file(fnm1, buffer, &filestatus);
                g_free(buffer);

                if (calc_MC) {
                    /* Monte Carlo */
                    printf(" Nmc %d \n", Nmc);
                    init_OPODRMCdata(&opodrmc, Nmc);
                    op_odr_uncertainty_montecarlo_run_calc(&opodrdata, &opodrunc, &opodrmc, &fddata,  &area);

                    if (opodrdata.Finput) {
                        sprintf(fnm1, "%s__op_odr_mc%d_range_%.2g_%.2g_%%.dat", fnm, Nmc, opodrdata.from_pct_Fmax, opodrdata.to_pct_Fmax);
                    }
                    else {
                        sprintf(fnm1, "%s__op_odr_mc%d_range_%.2g_%.2g_nm.dat", fnm, Nmc, opodrdata.from, opodrdata.to);
                    }

                    buffer = op_odr_mc_export_data(&opodrdata, &opodrunc, &opodrmc, &fddata, EXPORT_FORMAT_PLAIN); /* */
                    buffer_to_file(fnm1, buffer, &filestatus);
                    g_free(buffer);

                    destroy_OPODRMCdata(&opodrmc);
                }
            }
        }
    }

#endif

#ifndef NOFORTRAN

    if (calc_twoslopes) {
        /* TWOSLOPES */
        /*load */
        if (sldata.Finputload)
            range_to_indices(sldata.loadfrom_pct_Fmax * 0.01 * fddata.Fmax, sldata.loadto_pct_Fmax * 0.01 * fddata.Fmax,
                             fddata.Fload, FALSE, &istart, &iend, &ndata);
        else {
            range_to_indices(sldata.loadfrom, sldata.loadto, fddata.hload, FALSE, &istart, &iend, &ndata);
        }

        if (ndata < 2) {
            if (sldata.Finputload) {
                fprintf(stderr,  "Not enough data points for load range [%g, %g] %% from file %s\n",
                        sldata.loadfrom_pct_Fmax, sldata.loadto_pct_Fmax, fddata.filename);
            }
            else {
                fprintf(stderr,  "Not enough data points for load range [%g, %g] nm from file %s\n",
                        sldata.loadfrom, sldata.loadto, fddata.filename);
            }
        }
        else {
            ndata = MIN(ndata, fddata.hunload->res - istart);

            x = gwy_data_line_part_extract(fddata.hload, istart, ndata);
            y = gwy_data_line_part_extract(fddata.Fload, istart, ndata);

            if (ndata < 2) {
                if (sldata.Finputunload) {
                    fprintf(stderr,  "Not enough data points for unload range [%g, %g] %% from file %s\n",
                            sldata.unloadfrom_pct_Fmax, sldata.unloadto_pct_Fmax, fddata.filename);
                }
                else {
                    fprintf(stderr,  "Not enough data points for unload range [%g, %g] nm from file %s\n",
                            sldata.unloadfrom, sldata.unloadto, fddata.filename);
                }
            }
            else {
                sldata.nfitloaddata = ndata;

                fitinfo = slopes_fit_load(x, y, &sldata, &fddata, &instdata, NULL);
                sldata.xloadfit = gwy_data_line_duplicate(x);
                sldata.yloadfit = gwy_data_line_new_alike(y, FALSE);

                for (i = 0; i < sldata.xloadfit->res; i++) {
                    sldata.yloadfit->data[i] =  sldata.gamma * pow(sldata.xloadfit->data[i] - sldata.h0, sldata.n);
                }

                g_object_unref(x);
                g_object_unref(y);

                /*unload */
                if (sldata.Finputunload) {
                    range_to_indices(sldata.unloadfrom_pct_Fmax * 0.01 * fddata.Fmax, sldata.unloadto_pct_Fmax * 0.01 * fddata.Fmax,
                                     fddata.Funload, TRUE, &istart, &iend, &ndata);
                }
                else {
                    range_to_indices(sldata.unloadfrom, sldata.unloadto, fddata.hunload, TRUE, &istart, &iend, &ndata);
                }

                ndata = MIN(ndata, fddata.hunload->res - istart);

                x = gwy_data_line_part_extract(fddata.hunload, istart, ndata);
                y = gwy_data_line_part_extract(fddata.Funload, istart, ndata);
                sldata.nfitunloaddata = ndata;

                fitinfo = slopes_fit_unload(x, y, &sldata, &fddata, &instdata, NULL);
                sldata.xunloadfit = gwy_data_line_duplicate(x);
                sldata.yunloadfit = gwy_data_line_new_alike(y, FALSE);

                for (i = 0; i < sldata.xunloadfit->res; i++) {
                    sldata.yunloadfit->data[i] =  sldata.alpha * pow(sldata.xunloadfit->data[i] - sldata.hp, sldata.m);
                }

                g_object_unref(x);
                g_object_unref(y);

                if (sldata.Finputload) {
                    if (sldata.Finputunload) {
                        sprintf(fnm1, "%s__slopes_loadrange_%.2g_%.2g_%%_fit%d_unloadrange_%.2g_%.2g_%%.dat", fnm,
                                sldata.loadfrom_pct_Fmax, sldata.loadto_pct_Fmax, sldata.loadifixb[2], sldata.unloadfrom_pct_Fmax, sldata.unloadto_pct_Fmax);
                    }
                    else {
                        sprintf(fnm1, "%s__slopes_loadrange_%.2g_%.2g_%%_fit%d_unloadrange_%.2g_%.2g_nm.dat", fnm,
                                sldata.loadfrom_pct_Fmax, sldata.loadto_pct_Fmax, sldata.loadifixb[2], sldata.unloadfrom, sldata.unloadto);
                    }
                }
                else {
                    if (sldata.Finputunload) {
                        sprintf(fnm1, "%s__slopes_loadrange_%.2g_%.2g_nm_fit%d_unloadrange_%.2g_%.2g_%%.dat", fnm,
                                sldata.loadfrom, sldata.loadto, sldata.loadifixb[2], sldata.unloadfrom_pct_Fmax, sldata.unloadto_pct_Fmax);
                    }
                    else {
                        sprintf(fnm1, "%s__slopes_loadrange_%.2g_%.2g_nm_fit%d_unloadrange_%.2g_%.2g_nm.dat", fnm,
                                sldata.loadfrom, sldata.loadto, sldata.loadifixb[2], sldata.unloadfrom, sldata.unloadto);
                    }
                }

                buffer = slopes_export_data(&sldata, &fddata, &instdata, EXPORT_FORMAT_PLAIN); /* */
                buffer_to_file(fnm1, buffer, &filestatus);
                g_free(buffer);

                if (calc_unc) {

                    /*copy instdata to a local copy */
                    slunc.instdata = instdata;

                    /* uncertainties */
                    slopes_propagate_uncertainties_calc(&sldata, &slunc, &fddata);

                    slunc.Ec = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
                    slunc.Hc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
                    slunc.Ac = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));

                    for (j = -CONTACTMAX; j <= CONTACTMAX ; j++) {
                        /* use instdata instead of slunc.instdata, because they are by definition the same*/
                        slopes_fit_shift_contact(&sldata, &fddata, &instdata, slunc.Ec, slunc.Hc, slunc.Ac, j);
                    }

                    if (sldata.Finputload) {
                        if (sldata.Finputunload) {
                            sprintf(fnm1, "%s__slopes_unc_loadrange_%.2g_%.2g_%%_fit%d_unloadrange_%.2g_%.2g_%%.dat", fnm,
                                    sldata.loadfrom_pct_Fmax, sldata.loadto_pct_Fmax, sldata.loadifixb[2], sldata.unloadfrom_pct_Fmax, sldata.unloadto_pct_Fmax);
                        }
                        else {
                            sprintf(fnm1, "%s__slopes_unc_loadrange_%.2g_%.2g_%%_fit%d_unloadrange_%.2g_%.2g_nm.dat", fnm,
                                    sldata.loadfrom_pct_Fmax, sldata.loadto_pct_Fmax, sldata.loadifixb[2], sldata.unloadfrom, sldata.unloadto);
                        }
                    }
                    else {
                        if (sldata.Finputunload) {
                            sprintf(fnm1, "%s__slopes_unc_loadrange_%.2g_%.2g_nm_fit%d_unloadrange_%.2g_%.2g_%%.dat", fnm,
                                    sldata.loadfrom, sldata.loadto, sldata.loadifixb[2], sldata.unloadfrom_pct_Fmax, sldata.unloadto_pct_Fmax);
                        }
                        else {
                            sprintf(fnm1, "%s__slopes_unc_loadrange_%.2g_%.2g_nm_fit%d_unloadrange_%.2g_%.2g_nm.dat", fnm,
                                    sldata.loadfrom, sldata.loadto, sldata.loadifixb[2], sldata.unloadfrom, sldata.unloadto);
                        }
                    }

                    buffer = slopes_uncertainty_export_data(&sldata, &slunc, &fddata, EXPORT_FORMAT_PLAIN); /* */
                    buffer_to_file(fnm1, buffer, &filestatus);
                    g_free(buffer);

                    if (calc_MC) {
                        /* Monte Carlo */
                        init_SlopesMCdata(&slopesmc, Nmc);
                        slopes_uncertainty_montecarlo_run_calc(&sldata, &slunc, &slopesmc, &fddata);

                        if (sldata.Finputload) {
                            if (sldata.Finputunload) {
                                sprintf(fnm1, "%s__slopes_mc%d_loadrange_%.2g_%.2g_%%_unloadrange_%.2g_%.2g_%%.dat", fnm, Nmc,
                                        sldata.loadfrom_pct_Fmax, sldata.loadto_pct_Fmax, sldata.unloadfrom_pct_Fmax, sldata.unloadto_pct_Fmax);
                            }
                            else {
                                sprintf(fnm1, "%s__slopes_mc%d_range_%.2g_%.2g_%%_unloadrange_%.2g_%.2g_nm.dat", fnm, Nmc,
                                        sldata.loadfrom_pct_Fmax, sldata.loadto_pct_Fmax, sldata.unloadfrom, sldata.unloadto);
                            }
                        }
                        else {
                            if (sldata.Finputunload) {
                                sprintf(fnm1, "%s__slopes_mc%d_range_%.2g_%.2g_nm_unloadrange_%.2g_%.2g_%%.dat", fnm, Nmc,
                                        sldata.loadfrom, sldata.loadto, sldata.unloadfrom_pct_Fmax, sldata.unloadto_pct_Fmax);
                            }
                            else {
                                sprintf(fnm1, "%s__slopes_mc%d_range_%.2g_%.2g_nm_unloadrange_%.2g_%.2g_nm.dat", fnm, Nmc,
                                        sldata.loadfrom, sldata.loadto, sldata.unloadfrom, sldata.unloadto);
                            }
                        }

                        buffer = slopes_mc_export_data(&sldata, &slunc, &slopesmc, &fddata, EXPORT_FORMAT_PLAIN); /* */
                        buffer_to_file(fnm1, buffer, &filestatus);
                        g_free(buffer);

                        destroy_SlopesMCdata(&slopesmc);
                    }
                }
            }
        }
    }

#endif

    if (calc_hertz) {
        /* HERTZ */
        range_to_indices(hzdata.from, hzdata.to, fddata.hload, FALSE, &istart, &iend, &ndata);

        if (ndata < 2) {
            fprintf(stderr,  "Not enough data points for range [%g, %g] nm from file %s\n", hzdata.from, hzdata.to, fddata.filename);
        }
        else {
            x = gwy_data_line_part_extract(fddata.hload, istart, ndata);
            y = gwy_data_line_part_extract(fddata.Fload, istart, ndata);
            hertz_fit(x, y, &hzdata, &instdata);
            hzdata.xfit = gwy_data_line_duplicate(x);
            hzdata.yfit = gwy_data_line_part_extract(fddata.hload, 0, ndata);
            hzdata.nfitdata = ndata;

            for (i = 0; i < hzdata.xfit->res; i++) {
                hzdata.yfit->data[i] =  hzdata.reg.slope * pow(hzdata.xfit->data[i], 1.5) + hzdata.reg.intercept;
            }

            if (hzdata.mode == R_MODE) {
                sprintf(fnm1, "%s__hertz_range_%.2g_%.2g_nm_R_%.2g.dat", fnm, hzdata.from, hzdata.to, hzdata.radius);
            }
            else if (hzdata.mode == ER_MODE) {
                sprintf(fnm1, "%s__hertz_range_%.2g_%.2g_nm_Er_%.2g.dat", fnm, hzdata.from, hzdata.to, hzdata.Er);
            }
            else {
                sprintf(fnm1, "%s__hertz_range_%.2g_%.2g_nm_Eit_%.2g.dat", fnm, hzdata.from, hzdata.to, hzdata.Eit);
            }

            buffer = hertz_export_data(&hzdata, &fddata, &instdata, EXPORT_FORMAT_PLAIN); /* */
            buffer_to_file(fnm1, buffer, &filestatus);
            g_free(buffer);

            if (calc_unc) {
                /* uncertainties */
                hertz_propagate_uncertainties_calc(&hzdata, &hzunc, &instdata);

                hzunc.ERc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));

                for (j = -CONTACTMAX; j <= CONTACTMAX ; j++) {
                    hertz_fit_shift_contact(&hzdata, &fddata, &instdata, hzunc.ERc, j);
                }

                if (hzdata.mode == R_MODE) {
                    sprintf(fnm1, "%s__hertz_unc_range_%.2g_%.2g_nm_R_%.2g_uR_%.2g.dat", fnm, hzdata.from, hzdata.to, hzdata.radius, hzunc.uradius);
                }
                else if (hzdata.mode == ER_MODE) {
                    sprintf(fnm1, "%s__hertz_unc_range_%.2g_%.2g_nm_Er_%.2g_uEr_%.2g.dat", fnm, hzdata.from, hzdata.to, hzdata.Er, hzunc.uEr);
                }
                else {
                    sprintf(fnm1, "%s__hertz_unc_range_%.2g_%.2g_nm_Eit_%.2g_uEit_%.2g.dat", fnm, hzdata.from, hzdata.to, hzdata.Eit, hzunc.uEit);
                }

                buffer = hertz_uncertainty_export_data(&hzdata, &hzunc, &fddata, &instdata, EXPORT_FORMAT_PLAIN); /* */
                buffer_to_file(fnm1, buffer, &filestatus);
                g_free(buffer);

                if (calc_MC) {
                    /* Monte Carlo */
                    init_HertzMCdata(&hertzmc, Nmc);
                    hertz_uncertainty_montecarlo_run_calc(&hzdata, &hzunc, &hertzmc, &fddata, &instdata);

                    sprintf(fnm1, "%s__hertz_mc%d_range_%.2g_%.2g_nm_R_%.2g.dat", fnm, Nmc, hzdata.from, hzdata.to, hzdata.radius);

                    buffer = hertz_mc_export_data(&hzdata, &hzunc, &hertzmc, &fddata, &instdata, EXPORT_FORMAT_PLAIN); /* */
                    buffer_to_file(fnm1, buffer, &filestatus);
                    g_free(buffer);

                    destroy_HertzMCdata(&hertzmc);
                }
            }
        }
    }

#ifndef NOFORTRAN

    if (calc_hertz_odr) {
        double b, w;
        /* HERTZ ODR */
        range_to_indices(hzodrdata.from, hzodrdata.to, fddata.hload, FALSE, &istart, &iend, &ndata);

        if (ndata < 2) {
            fprintf(stderr,  "Not enough data points for range [%g, %g] nm from file %s\n", hzodrdata.from, hzodrdata.to, fddata.filename);
        }
        else {
            x = gwy_data_line_part_extract(fddata.hload, istart, ndata);
            y = gwy_data_line_part_extract(fddata.Fload, istart, ndata);
            hzodrdata.nfitdata = ndata;
            nremove = filter_negative_data(x, y);
            hzodrdata.nfitdata -= nremove;
            fitinfo = hertz_odr_fit(x, y, &hzodrdata, &fddata, &instdata, NULL);
            hzodrdata.xfit = gwy_data_line_duplicate(x);
            hzodrdata.yfit = gwy_data_line_part_extract(fddata.hload, 0, ndata);

            if (hzodrdata.radial_corr) {
                b = 2. / 3. / M_PI * (1 - 2 * instdata.nu) / (1 - instdata.nu) / sqrt(hzodrdata.radius * 1e9); /* must be in nanometers*/

                for (i = 0; i < hzodrdata.xfit->res; i++) {
                    w = sqrt(hzodrdata.xfit->data[i] - hzodrdata.h0);
                    hzodrdata.yfit->data[i] =  hzodrdata.gamma * (1 + b * w) * w * w * w;
                }
            }
            else {
                for (i = 0; i < hzodrdata.xfit->res; i++) {
                    hzodrdata.yfit->data[i] =  hzodrdata.gamma * pow(hzodrdata.xfit->data[i] - hzodrdata.h0, hzodrdata.n);
                }
            }

            g_object_unref(x);
            g_object_unref(y);

            if (hzodrdata.mode == R_MODE) {
                sprintf(fnm1, "%s__hertz_odr_range_%.2g_%.2g_nm_R_%.2g_fit%d_radial%d.dat", fnm,
                        hzodrdata.from, hzodrdata.to, hzodrdata.radius, hzodrdata.ifixb[2], hzodrdata.radial_corr);
            }
            else if (hzodrdata.mode == ER_MODE) {
                sprintf(fnm1, "%s__hertz_odr_range_%.2g_%.2g_nm_Er_%.2g_fit%d.dat", fnm,
                        hzodrdata.from, hzodrdata.to, hzodrdata.Er, hzodrdata.ifixb[2]);
            }
            else {
                sprintf(fnm1, "%s__hertz_odr_range_%.2g_%.2g_nm_Eit_%.2g_fit%d.dat", fnm,
                        hzodrdata.from, hzodrdata.to, hzodrdata.Eit, hzodrdata.ifixb[2]);
            }

            buffer = hertz_odr_export_data(&hzodrdata, &fddata, &instdata, EXPORT_FORMAT_PLAIN); /* */
            buffer_to_file(fnm1, buffer, &filestatus);
            g_free(buffer);

            if (calc_unc) {
                /* uncertainties */

                /*copy instdata to a local copy */
                hzodrunc.instdata = instdata;

                hertz_odr_propagate_uncertainties_calc(&hzodrdata, &hzodrunc, &fddata);

                hzodrunc.ERc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));

                for (j = -CONTACTMAX; j <= CONTACTMAX; j++) {
                    /* use instdata instead of hzodrunc.instdata, because they are by definition the same*/
                    hertz_odr_fit_shift_contact(&hzodrdata, &fddata, &instdata, hzodrunc.ERc, j);
                }

                if (hzodrdata.mode == R_MODE) {
                    sprintf(fnm1, "%s__hertz_odr_unc_range_%.2g_%.2g_nm_R_%.2g_fit%d_uR_%.2g.dat", fnm,
                            hzodrdata.from, hzodrdata.to, hzodrdata.radius, hzodrdata.ifixb[2], hzodrunc.uradius);
                }
                else if (hzodrdata.mode == ER_MODE) {
                    sprintf(fnm1, "%s__hertz_odr_unc_range_%.2g_%.2g_nm_Er_%.2g_fit%d_uEr_%.2g.dat", fnm,
                            hzodrdata.from, hzodrdata.to, hzodrdata.Er, hzodrdata.ifixb[2], hzodrunc.uEr);
                }
                else {
                    sprintf(fnm1, "%s__hertz_odr_unc_range_%.2g_%.2g_nm_Eit_%.2g_fit%d_uEit_%.2g.dat", fnm,
                            hzodrdata.from, hzodrdata.to, hzodrdata.Eit, hzodrdata.ifixb[2], hzodrunc.uEit);
                }

                buffer = hertz_odr_uncertainty_export_data(&hzodrdata, &hzodrunc, &fddata, EXPORT_FORMAT_PLAIN); /* */
                buffer_to_file(fnm1, buffer, &filestatus);
                g_free(buffer);

                if (calc_MC) {
                    /* Monte Carlo */
                    init_HertzODRMCdata(&hertzodrmc, Nmc);
                    hertz_odr_uncertainty_montecarlo_run_calc(&hzodrdata, &hzodrunc, &hertzodrmc, &fddata);

                    sprintf(fnm1, "%s__hertz_odr_mc%d_range_%.2g_%.2g_nm_R_%.2g.dat", fnm, Nmc, hzodrdata.from, hzodrdata.to, hzodrdata.radius);

                    buffer = hertz_odr_mc_export_data(&hzodrdata, &hzodrunc, &hertzodrmc, &fddata, EXPORT_FORMAT_PLAIN); /* */
                    buffer_to_file(fnm1, buffer, &filestatus);
                    g_free(buffer);

                    destroy_HertzODRMCdata(&hertzodrmc);
                }
            }
        }
    }

#endif

#ifndef NOFORTRAN

    if (calc_stiffness) {
        /* STIFFNESS */
        /*load */
        range_to_indices(stiffdata.loadfrom, stiffdata.loadto, fddata.hload, FALSE, &istart, &iend, &ndata);

        if (ndata < 2) {
            fprintf(stderr,  "Not enough data points for load range [%g, %g] nm from file %s\n",
                    stiffdata.loadfrom, stiffdata.loadto, fddata.filename);
        }
        else {
            ndata = MIN(ndata, fddata.hload->res - istart);
            stiffdata.nfitloaddata = ndata;

            x = gwy_data_line_part_extract(fddata.hload, istart, ndata);
            y = gwy_data_line_part_extract(fddata.Fload, istart, ndata);

            fitinfo = stiffness_fit_load(x, y, &stiffdata, &fddata, &instdata, NULL);
            stiffdata.xloadfit = gwy_data_line_duplicate(x);
            stiffdata.yloadfit = gwy_data_line_new_alike(y, FALSE);

            for (i = 0; i < stiffdata.xloadfit->res; i++) {
                stiffdata.yloadfit->data[i] =  stiffdata.kl * stiffdata.xloadfit->data[i] + stiffdata.ql;
            }

            g_object_unref(x);
            g_object_unref(y);
        }

        /*unload */
        range_to_indices(stiffdata.unloadfrom, stiffdata.unloadto, fddata.hunload, TRUE, &istart, &iend, &ndata);

        if (ndata < 2) {
            fprintf(stderr,  "Not enough data points for load range [%g, %g] nm from file %s\n",
                    stiffdata.loadfrom, stiffdata.loadto, fddata.filename);
        }

        else {
            ndata = MIN(ndata, fddata.hunload->res - istart);
            stiffdata.nfitunloaddata = ndata;

            x = gwy_data_line_part_extract(fddata.hunload, istart, ndata);
            y = gwy_data_line_part_extract(fddata.Funload, istart, ndata);

            fitinfo = stiffness_fit_unload(x, y, &stiffdata, &fddata, &instdata, NULL);
            stiffdata.xunloadfit = gwy_data_line_duplicate(x);
            stiffdata.yunloadfit = gwy_data_line_new_alike(y, FALSE);

            for (i = 0; i < stiffdata.xunloadfit->res; i++) {
                stiffdata.yunloadfit->data[i] =  stiffdata.ku * stiffdata.xunloadfit->data[i] + stiffdata.qu;
            }

            g_object_unref(x);
            g_object_unref(y);
        }

        if (stiffdata.has_fit_load && stiffdata.has_fit_unload) {
            stiffness_combine_fit_results(&stiffdata, &fddata, &instdata);
        }

        sprintf(fnm1, "%s__stiffness_loadrange_%.2g_%.2g_nm_unloadrange_%.2g_%.2g_nm.dat", fnm,
                stiffdata.loadfrom, stiffdata.loadto,  stiffdata.unloadfrom, stiffdata.unloadto);

        buffer = stiffness_export_data(&stiffdata, &fddata, &instdata, EXPORT_FORMAT_PLAIN); /* */
        buffer_to_file(fnm1, buffer, &filestatus);
        g_free(buffer);

        if (calc_unc) {
            /* uncertainties */

            /*copy instdata to a local copy */
            stiffunc.instdata = instdata;

            stiffness_propagate_uncertainties_calc(&stiffdata, &stiffunc, &fddata);

            stiffunc.klc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
            stiffunc.qlc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
            stiffunc.kuc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));
            stiffunc.quc = (gdouble *)g_malloc((2 * CONTACTMAX + 1) * sizeof(gdouble));

            for (j = -CONTACTMAX; j <= CONTACTMAX ; j++) {
                /* use instdata instead of stiffunc.instdata, because they are by definition the same*/
                stiffness_fit_shift_contact(&stiffdata, &fddata, &instdata, stiffunc.klc, stiffunc.qlc, stiffunc.kuc, stiffunc.quc, j);
            }

            sprintf(fnm1, "%s__stiffness_unc_loadrange_%.2g_%.2g_nm_unloadrange_%.2g_%.2g_nm.dat", fnm,
                    stiffdata.loadfrom, stiffdata.loadto, stiffdata.unloadfrom, stiffdata.unloadto);

            buffer = stiffness_uncertainty_export_data(&stiffdata, &stiffunc, &fddata, EXPORT_FORMAT_PLAIN); /* */
            buffer_to_file(fnm1, buffer, &filestatus);
            g_free(buffer);

            if (calc_MC) {
                /* Monte Carlo */
                init_StiffnessMCdata(&stiffmc, Nmc);
                stiffness_uncertainty_montecarlo_run_calc(&stiffdata, &stiffunc, &stiffmc, &fddata);

                sprintf(fnm1, "%s__stiffness_mc%d_range_%.2g_%.2g_nm_unloadrange_%.2g_%.2g_nm.dat", fnm, Nmc,
                        stiffdata.loadfrom, stiffdata.loadto, stiffdata.unloadfrom, stiffdata.unloadto);

                buffer = stiffness_mc_export_data(&stiffdata, &stiffunc, &stiffmc, &fddata, EXPORT_FORMAT_PLAIN); /* */
                buffer_to_file(fnm1, buffer, &filestatus);
                g_free(buffer);

                destroy_StiffnessMCdata(&stiffmc);
            }
        }
    }

#endif

    if (calc_Ph2) {
        /* P/h2 */
        if (ph2data.nmove % 2 == 0) {
            printf("Moving average cannot be taken over an even number of points. \n");
            exit(2);
        }

        Ph2_calc(&fddata, &ph2data);
        sprintf(fnm1, "%s__Ph2_%d.dat", fnm, ph2data.nmove);

        buffer = Ph2_export_data(&ph2data, &fddata, EXPORT_FORMAT_PLAIN); /* */
        buffer_to_file(fnm1, buffer, &filestatus);
        g_free(buffer);
    }

    if (calc_apopin) {
        apopin_create_aux_datalines(&apopdata, &fddata);
        apopin_calc(&apopdata, &fddata);
        sprintf(fnm1, "%s__popin.dat", fnm);
        printf("%s %d \n", fnm1, apopdata.npopin);

        buffer = apopin_export_data(&apopdata, &fddata, EXPORT_FORMAT_PLAIN); /* */
        buffer_to_file(fnm1, buffer, &filestatus);
        g_free(buffer);
    }

    if (calc_work) {
        work_calc(&wdata, &fddata);
        sprintf(fnm1, "%s__work.dat", fnm);

        buffer = work_export_data(&wdata, &fddata, EXPORT_FORMAT_PLAIN); /* */
        buffer_to_file(fnm1, buffer, &filestatus);
        g_free(buffer);
    }

    if (write_data) {
        sprintf(fnm1, "%s__data.dat", fnm);

        buffer = fddata_export_data(&fddata, EXPORT_FORMAT_PLAIN); /* */
        buffer_to_file(fnm1, buffer, &filestatus);
        g_free(buffer);
    }

    return 0;
}

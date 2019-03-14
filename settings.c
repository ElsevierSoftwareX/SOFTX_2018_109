#include "settings.h"

#include "niget-common.h"

#ifdef G_OS_WIN32
#include "getline.h"
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>


#ifndef NOFORTRAN
enum {
    FC = 0,
    OP_BETA,
    OP_UNC_H,
    OP_UNC_F,
    OP_UNC_NU, /* preserved for bw compatibility */
    OP_UNC_NUI, /* preserved for bw compatibility */
    OP_UNC_EI, /* preserved for bw compatibility */
    OP_NMC,
    OPODR_BETA,
    OPODR_UNC_H,
    OPODR_UNC_F,
    OPODR_UNC_NU, /* preserved for bw compatibility */
    OPODR_UNC_NUI, /* preserved for bw compatibility */
    OPODR_UNC_EI, /* preserved for bw compatibility */
    OPODR_NMC,
    TG_BETA,
    TG_UNC_H,
    TG_UNC_F,
    TG_UNC_NU, /* preserved for bw compatibility */
    TG_UNC_NUI, /* preserved for bw compatibility */
    TG_UNC_EI, /* preserved for bw compatibility */
    TG_NMC,
    HZ_UNC_H,
    HZ_UNC_F,
    HZ_UNC_NU, /* preserved for bw compatibility */
    HZ_UNC_NUI, /* preserved for bw compatibility */
    HZ_UNC_EI, /* preserved for bw compatibility */
    HZ_UNC_RADIUS,
    HZ_NMC,
    HZODR_UNC_H,
    HZODR_UNC_F,
    HZODR_UNC_NU, /* preserved for bw compatibility */
    HZODR_UNC_NUI, /* preserved for bw compatibility */
    HZODR_UNC_EI, /* preserved for bw compatibility */
    HZODR_UNC_RADIUS,
    HZODR_NMC,
    SLOPES_BETA,
    SLOPES_UNC_H,
    SLOPES_UNC_F,
    SLOPES_UNC_NU, /* preserved for bw compatibility */
    SLOPES_UNC_NUI, /* preserved for bw compatibility */
    SLOPES_UNC_EI, /* preserved for bw compatibility */
    SLOPES_NMC,
    APOP_THRESH,
    APOP_THRESH2,
    APOP_WPOP,
    APOP_HPOP,
    APOP_NMOVE,
    PH2_NMOVE,
    WORK_NMOVE,
    INST_NU,
    INST_NUI,
    INST_EI,
    INST_UNC_NU,
    INST_UNC_NUI,
    INST_UNC_EI,
    INST_SIGMA_H,
    INST_SIGMA_F,
    AREA_POLY, // MUST be the last item
};
#else
enum {
    FC = 0,
    OP_BETA,
    OP_UNC_H,
    OP_UNC_F,
    OP_UNC_NU, /* preserved for bw compatibility */
    OP_UNC_NUI, /* preserved for bw compatibility */
    OP_UNC_EI, /* preserved for bw compatibility */
    OP_NMC,
    TG_BETA,
    TG_UNC_H,
    TG_UNC_F,
    TG_UNC_NU, /* preserved for bw compatibility */
    TG_UNC_NUI, /* preserved for bw compatibility */
    TG_UNC_EI, /* preserved for bw compatibility */
    TG_NMC,
    HZ_UNC_H,
    HZ_UNC_F,
    HZ_UNC_NU, /* preserved for bw compatibility */
    HZ_UNC_NUI, /* preserved for bw compatibility */
    HZ_UNC_EI, /* preserved for bw compatibility */
    HZ_UNC_RADIUS,
    HZ_NMC,
    APOP_THRESH,
    APOP_THRESH2,
    APOP_WPOP,
    APOP_HPOP,
    APOP_NMOVE,
    PH2_NMOVE,
    WORK_NMOVE,
    INST_NU,
    INST_NUI,
    INST_EI,
    INST_UNC_NU,
    INST_UNC_NUI,
    INST_UNC_EI,
    INST_SIGMA_H,
    INST_SIGMA_F,
    AREA_POLY,// MUST be the last item
};
#endif

void save_settings(const gchar *filename, const Args *args)
{
    gint i;
    FILE *file;

    file = fopen(filename, "w");

    if (!file) {
        return;
    }

    fprintf(file, "Settings %s\n", g_get_application_name());

    fprintf(file, "%d  %g\n", FC, args->fddata.Fc);

    fprintf(file, "%d  %g\n", OP_BETA, args->opdata.beta);
    fprintf(file, "%d  %g\n", OP_UNC_H, args->opunc.uh);
    fprintf(file, "%d  %g\n", OP_UNC_F, args->opunc.uF);
    fprintf(file, "%d  %d\n", OP_NMC, args->opunc.Nmc);

#ifndef NOFORTRAN
    fprintf(file, "%d  %g\n", OPODR_BETA, args->opodrdata.beta);
    //   fprintf(file, "%d  %g\n", OPODR_UNC_H, args->opodrunc.uh);
    //  fprintf(file, "%d  %g\n", OPODR_UNC_F, args->opodrunc.uF);
    fprintf(file, "%d  %d\n", OPODR_NMC, args->opodrunc.Nmc);
#endif

    fprintf(file, "%d  %g\n", TG_BETA, args->tgdata.beta);
    fprintf(file, "%d  %g\n", TG_UNC_H, args->tgunc.uh);
    fprintf(file, "%d  %g\n", TG_UNC_F, args->tgunc.uF);
    fprintf(file, "%d  %d\n", TG_NMC, args->tgunc.Nmc);

    fprintf(file, "%d  %g\n", HZ_UNC_H, args->hertzunc.uh);
    fprintf(file, "%d  %g\n", HZ_UNC_F, args->hertzunc.uF);
    fprintf(file, "%d  %g\n", HZ_UNC_RADIUS, args->hertzunc.uradius);
    fprintf(file, "%d  %d\n", HZ_NMC, args->hertzunc.Nmc);

#ifndef NOFORTRAN
    // fprintf(file, "%d  %g\n", HZODR_UNC_H, args->hertzodrunc.uh);
    //fprintf(file, "%d  %g\n", HZODR_UNC_F, args->hertzodrunc.uF);
    fprintf(file, "%d  %g\n", HZODR_UNC_RADIUS, args->hertzodrunc.uradius);
    fprintf(file, "%d  %d\n", HZODR_NMC, args->hertzodrunc.Nmc);

    fprintf(file, "%d  %g\n", SLOPES_BETA, args->slopesdata.beta);
    // fprintf(file, "%d  %g\n", SLOPES_UNC_H, args->slopesunc.uh);
    //fprintf(file, "%d  %g\n", SLOPES_UNC_F, args->slopesunc.uF);
    fprintf(file, "%d  %d\n", SLOPES_NMC, args->slopesunc.Nmc);
#endif

    fprintf(file, "%d  %d\n", PH2_NMOVE, args->ph2data.nmove);

    fprintf(file, "%d  %g\n", APOP_THRESH, args->apopdata.thresh);
    fprintf(file, "%d  %g\n", APOP_THRESH2, args->apopdata.thresh2);
    fprintf(file, "%d  %d\n", APOP_WPOP, args->apopdata.wpop);
    fprintf(file, "%d  %g\n", APOP_HPOP, args->apopdata.hpop);
    fprintf(file, "%d  %d\n", APOP_NMOVE, args->apopdata.nmove);

    fprintf(file, "%d  %d\n", WORK_NMOVE, args->workdata.nmove);

    fprintf(file, "%d  %g\n", INST_NU, args->instdata.nu);
    fprintf(file, "%d  %g\n", INST_NUI, args->instdata.nui);
    fprintf(file, "%d  %g\n", INST_EI, args->instdata.Ei);
    fprintf(file, "%d  %g\n", INST_UNC_NU, args->instdata.unu);
    fprintf(file, "%d  %g\n", INST_UNC_NUI, args->instdata.unui);
    fprintf(file, "%d  %g\n", INST_UNC_EI, args->instdata.uEi);
    fprintf(file, "%d  %g\n", INST_SIGMA_H, args->instdata.sigma_h);
    fprintf(file, "%d  %g\n", INST_SIGMA_F, args->instdata.sigma_F);

    fprintf(file, "%d %d \n", AREA_POLY, args->area.npoly);

    for (i = 0 ; i < args->area.npoly; i++) {
        fprintf(file, "%d %g \n", AREA_POLY + 1 + i, args->area.polycoef[i]);
    }

    fclose(file);
}

void load_settings(const gchar *filename, Args *args)
{
    FILE *file;
    gchar *line = NULL;
    size_t size;
    gint key;
    gdouble value;
    gint n, npoly;
    gboolean polyok = TRUE;

    file = fopen(filename, "r");

    if (!file) {
        return;
    }

    getline(&line, &size, file);

    if (strstr(line, g_get_application_name()) == NULL) {
        return;
    }

    while (getline(&line, &size, file) != EOF) {
        localize_decimal_points(line);
        n = sscanf(line, " %d %lf", &key, &value);

        if (n == 2) {
            switch (key) {
            case (FC):
                args->fddata.Fc = MAX(value, 0);
                break;

            case (OP_BETA):
                args->opdata.beta = MAX(value, 0);
                break;

            case (OP_UNC_H):
                args->opunc.uh = MAX(value, 0);
                break;

            case (OP_UNC_F):
                args->opunc.uF = MAX(value, 0);
                break;

            case (OP_NMC):
                args->opunc.Nmc = MAX((gint)round(value), 1);
                break;
#ifndef NOFORTRAN

            case (OPODR_BETA):
                args->opodrdata.beta = MAX(value, 0);
                break;

            /*            case (OPODR_UNC_H):
                            args->opodrunc.uh = MAX(value, 0);
                            break;

                        case (OPODR_UNC_F):
                            args->opodrunc.uF = MAX(value, 0);
                            break;
                            */

            case (OPODR_NMC):
                args->opodrunc.Nmc = MAX((gint)round(value), 1);
                break;
#endif

            case (TG_BETA):
                args->tgdata.beta = MAX(value, 0);
                break;

            case (TG_NMC):
                args->tgunc.Nmc = MAX((gint)round(value), 1);
                break;

            case (TG_UNC_H):
                args->tgunc.uh = MAX(value, 0);
                break;

            case (TG_UNC_F):
                args->tgunc.uF = MAX(value, 0);
                break;

            case (HZ_UNC_H):
                args->hertzunc.uh = MAX(value, 0);
                break;

            case (HZ_UNC_F):
                args->hertzunc.uF = MAX(value, 0);
                break;

            case (HZ_UNC_RADIUS):
                args->hertzunc.uradius = MAX(value, 0);
                break;

            case (HZ_NMC):
                args->hertzunc.Nmc = MAX((gint)round(value), 1);
                break;
#ifndef NOFORTRAN

            /*
            case (HZODR_UNC_H):
            args->hertzodrunc.uh = MAX(value, 0);
            break;

            case (HZODR_UNC_F):
            args->hertzodrunc.uF = MAX(value, 0);
            break;
            */

            case (HZODR_UNC_RADIUS):
                args->hertzodrunc.uradius = MAX(value, 0);
                break;

            case (HZODR_NMC):
                args->hertzodrunc.Nmc = MAX((gint)round(value), 1);
                break;

            case (SLOPES_BETA):
                args->slopesdata.beta = MAX(value, 0);
                break;

            /*
            case (SLOPES_UNC_H):
            args->slopesunc.uh = MAX(value, 0);
            break;

            case (SLOPES_UNC_F):
            args->slopesunc.uF = MAX(value, 0);
            break;

            */

            case (SLOPES_NMC):
                args->slopesunc.Nmc = MAX((gint)round(value), 1);
                break;
#endif

            case (APOP_THRESH):
                args->apopdata.thresh = value;
                break;

            case (APOP_THRESH2):
                args->apopdata.thresh2 = value;
                break;

            case (APOP_WPOP):
                args->apopdata.wpop = MAX((gint)round(value), 0);
                break;

            case (APOP_HPOP):
                args->apopdata.hpop = value;
                break;

            case (APOP_NMOVE):
                args->apopdata.nmove = MAX((gint)round(value), 1);
                break;

            case (PH2_NMOVE):
                args->ph2data.nmove = MAX((gint)round(value), 1);
                break;

            case (WORK_NMOVE):
                args->workdata.nmove = MAX((gint)round(value), 1);
                break;

            case (INST_NU):
                args->instdata.nu = CLAMP(value, -1.0, 0.5);
                break;

            case (INST_NUI):
                args->instdata.nui = CLAMP(value, -1.0, 0.5);
                break;

            case (INST_EI):
                args->instdata.Ei = MAX(value, 0);
                break;

            case (INST_UNC_NU):
                args->instdata.unu = CLAMP(value, 0, 0.5);
                break;

            case (INST_UNC_NUI):
                args->instdata.unui = CLAMP(value, 0, 0.5);
                break;

            case (INST_UNC_EI):
                args->instdata.uEi = MAX(value, 0);
                break;

            case (INST_SIGMA_H):
                args->instdata.sigma_h = MAX(value, 0);
                args->instdata.delta = args->instdata.sigma_F * args->instdata.sigma_F / args->instdata.sigma_h / args->instdata.sigma_h;
                break;

            case (INST_SIGMA_F):
                args->instdata.sigma_F = MAX(value, 0);
                args->instdata.delta = args->instdata.sigma_F * args->instdata.sigma_F / args->instdata.sigma_h / args->instdata.sigma_h;
                break;

            case (AREA_POLY):
                npoly = (gint)round(value);
                polyok = (args->area.npoly >= npoly);
                break;

            default:
                if (polyok  && (key - AREA_POLY <= args->area.npoly)) {
                    args->area.polycoef[key - AREA_POLY - 1] = value;
                    break;
                }
                else {
                    g_printerr("Bad settings %d \n", key);
                }

                break;
            }
        }
    }
}

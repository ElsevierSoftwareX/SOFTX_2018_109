#include "tool_apopins.h"

#include "datatypes.h"
#include "niget-common.h"

#include <math.h>

#include <libprocess/gwyprocess.h>


void init_Apopindata(Apopindata *apopdata)
{
    apopdata->nmove = 1;
    apopdata->npopin = 0;
    apopdata->thresh = THRESH_DEFAULT;
    apopdata->thresh2 = THRESH2_DEFAULT;
    apopdata->wpop = WPOP_DEFAULT;
    apopdata->hpop = HPOP_DEFAULT;
    apopdata->Fpopin = NULL;
    apopdata->hpopin = NULL;
    apopdata->dhpopin = NULL;
    apopdata->ileft = NULL;
    apopdata->iright = NULL;
    apopdata->havg = NULL;
    apopdata->dh = NULL;
    apopdata->t = NULL;
    apopdata->dth = NULL;
    apopdata->dth2 = NULL;
}

gboolean apopin_has_results(const Apopindata *apopdata)
{
    return (apopdata->npopin > 0);
}

void apopin_remove_results(Apopindata *apopdata)
{
    // remove all results
    if (apopdata->npopin > 0) {
        g_free(apopdata->Fpopin);
        g_free(apopdata->hpopin);
        g_free(apopdata->dhpopin);
        g_free(apopdata->ileft);
        g_free(apopdata->iright);
        apopdata->Fpopin = NULL;
        apopdata->hpopin = NULL;
        apopdata->dhpopin = NULL;
        apopdata->ileft = NULL;
        apopdata->iright = NULL;
    }

    apopdata->npopin = 0;
}

void apopin_create_aux_datalines(Apopindata *apopdata, const FDdata *fddata)
{
    gint i;

    // check if pseudo-time data have been already created, if not create them
    //check and possibly create constant thresh, thresh2 datalines
    if (fddata->hload) {
        if (apopdata->t == NULL) {
            apopdata->t = gwy_data_line_new_alike(fddata->hload, TRUE); /* LEAK */

            for (i = 0; i < apopdata->t->res; i++) {
                apopdata->t->data[i] = i;
            }
        }
        else {
            if (apopdata->t->res != fddata->hload->res) {
                gwy_data_line_resample(apopdata->t, fddata->hload->res, GWY_INTERPOLATION_NONE);

                for (i = 0; i < apopdata->t->res; i++) {
                    apopdata->t->data[i] = i;
                }
            }
        }

        if (apopdata->dth == NULL) {
            apopdata->dth = gwy_data_line_new_alike(fddata->hload, TRUE); /* LEAK */
            gwy_data_line_add(apopdata->dth, apopdata->thresh);
        }
        else {
            gwy_data_line_resample(apopdata->dth, fddata->hload->res, GWY_INTERPOLATION_NONE);
            gwy_data_line_clear(apopdata->dth);
            gwy_data_line_add(apopdata->dth, apopdata->thresh);
        }

        if (apopdata->dth2 == NULL) {
            apopdata->dth2 = gwy_data_line_new_alike(fddata->hload, TRUE); /* LEAK */
            gwy_data_line_add(apopdata->dth2, apopdata->thresh2);
        }
        else {
            gwy_data_line_resample(apopdata->dth2, fddata->hload->res, GWY_INTERPOLATION_NONE);
            gwy_data_line_clear(apopdata->dth2);
            gwy_data_line_add(apopdata->dth2, apopdata->thresh2);
        }
    }
}

void apopin_calc(Apopindata *apopdata, const FDdata *fddata)
{
    GwyDataLine *x, *xavg;
    GwyDataLine *dx;
    gint n ;
    gint i, j;
    /* gdouble xm; */
    gdouble thresh, upvalue;

    gint prev;
    gint ngr;
    gint nind;
    gint *indices = NULL;
    gint **groups = NULL;
    gint *ngrind = NULL;

    gint wpop;
    gdouble hpop;
    gdouble thresh2;
    gint ileft;
    gint iright;

    gint npopin;

    x = fddata->hload;
    n = x->res;
    thresh = apopdata->thresh;
    thresh2 = apopdata->thresh2;
    wpop = apopdata->wpop;
    hpop = apopdata->hpop;

    // calculate sliding average of depth data
    xavg = moving_average(x, apopdata->nmove);

    //calculate derivative of smoothed depth
    dx = gwy_data_line_new_alike(x, FALSE);

    for (i = 0; i < n; i++) {
        dx->data[i] = gwy_data_line_get_der(xavg, i);
    }

    /* xm = gwy_data_line_get_avg(dx); */
    /* upvalue = xm + thresh; // works but would require plotting upvalue instead of thresh in apopdata->dth */
    upvalue = thresh;

    indices = (gint *)g_malloc(n * sizeof(gint));
    nind = 0;

    // find indices where the derivative is larger than a certain value (= average derivative + threshhold) or (=threshhold)
    for (i = 0; i < n; i++) {
        if (dx->data[i] > upvalue) {
            indices[nind++] = i;
        }
    }

    if (nind == 0) {
        ngr = 0;
        npopin = 0;
    }
    else {
        indices = (gint *)g_realloc(indices, nind * sizeof(gint));
        groups = (gint **) g_malloc(nind * sizeof(gint *));
        prev = indices[0];

        for (i = 0; i < nind; i++) {
            groups[i] = (gint *) g_malloc(nind * sizeof(gint));
        }

        if (verbose) {
            g_print(" nind %d \n", nind);
            g_print("indices \n");

            for (i = 1; i < nind; i++) {
                g_print("indices[%d] %d \n", i, indices[i]);
            }
        }

        ngrind = (gint *) g_malloc(nind * sizeof(gint));

        // split indices into groups with consecutive indices, ngr = number of these groups, ngrind = number of indices in each group, groups[igr] = indices belonging to group[igr]
        groups[0][0] = prev;
        j = 1;
        ngr = 0;

        if (verbose) {
            g_print("--------\n");
        }

        for (i = 1; i < nind; i++) {
            if (indices[i] == prev + 1) {
                groups[ngr][j] = indices[i];
                prev++;
                j++;

                if (verbose) {
                    g_print("add indices[%d]=%d to group %d \n", i, indices[i], ngr);
                }
            }
            else {
                if (verbose) {
                    g_print("finish previous group at ngr %d \n", ngr);
                }

                ngrind[ngr] = j;
                ngr++;

                if (verbose) {
                    g_print("start new group at ngr %d \n", ngr);
                }

                prev = indices[i];
                j = 0;
                groups[ngr][j] = indices[i];
                j++;
            }
        }

        ngrind[ngr] = j;

        if (verbose) {
            g_print("finish previous group at ngr %d \n", ngr);
        }

        ngr++;
    }

    if (verbose) {
        g_print(" ngr %d \n", ngr);

        for (i = 0; i < ngr; i++) {
            g_print("ngrind[%d] %d \n", i, ngrind[i]);

            for (j = 0; j < ngrind[i]; j++) {
                g_print(" groups[%d][%d] %d \n", i, j, groups[i][j]);
            }
        }
    }

    //copy local data to Data structure
    apopdata->havg =  gwy_data_line_duplicate(xavg);
    apopdata->dh =  gwy_data_line_duplicate(dx);

    // clean up old data if necessary
    if (apopdata->Fpopin != NULL) {
        g_free(apopdata->Fpopin);
    }

    if (apopdata->hpopin != NULL) {
        g_free(apopdata->hpopin);
    }

    if (apopdata->dhpopin != NULL) {
        g_free(apopdata->dhpopin);
    }

    if (apopdata->ileft != NULL) {
        g_free(apopdata->ileft);
    }

    if (apopdata->iright != NULL) {
        g_free(apopdata->iright);
    }

    //allocate new data
    apopdata->Fpopin = (gdouble *)g_malloc(ngr * sizeof(gdouble));
    apopdata->hpopin = (gdouble *)g_malloc(ngr * sizeof(gdouble));
    apopdata->dhpopin = (gdouble *)g_malloc(ngr * sizeof(gdouble));
    apopdata->ileft = (gint *)g_malloc(ngr * sizeof(gdouble));
    apopdata->iright = (gint *)g_malloc(ngr * sizeof(gdouble));

    // starting from a group of consecutive indices, add more indices at both ends so that the derivative is everywhere above threshhold2
    // if the width of this group is larger than wpop and the height jump associated with it larger than hpop, then this is a pop-in
    npopin = 0;

    for (i = 0; i < ngr; i++) {
        ileft = groups[i][0];

        while (ileft > 0 && dx->data[ileft] > thresh2 && dx->data[ileft] > dx->data[ileft - 1]) {
            ileft--;
        }

        iright = groups[i][ngrind[i] - 1];

        while (iright < n - 1 && dx->data[iright] > thresh2 && dx->data[iright] > dx->data[iright + 1]) {
            iright++;
        }

        if ((iright - ileft > wpop) & (x->data[iright] - x->data[ileft] > hpop)) {
            //it's a real pop-in, write the data to the Data structure
            apopdata->Fpopin[npopin] = gwy_data_line_get_avg(gwy_data_line_part_extract(fddata->Fload, ileft, iright - ileft + 1));
            apopdata->hpopin[npopin] = fddata->hload->data[ileft];
            apopdata->dhpopin[npopin] =  fddata->hload->data[iright] - fddata->hload->data[ileft];
            apopdata->ileft[npopin] = ileft;
            apopdata->iright[npopin] = iright;
            npopin++;
        }
    }

    //realloc data to match
    apopdata->Fpopin = (gdouble *)g_realloc(apopdata->Fpopin, npopin * sizeof(gdouble));
    apopdata->hpopin = (gdouble *)g_realloc(apopdata->hpopin, npopin * sizeof(gdouble));
    apopdata->dhpopin = (gdouble *)g_realloc(apopdata->dhpopin, npopin * sizeof(gdouble));
    apopdata->ileft = (gint *)g_realloc(apopdata->ileft, npopin * sizeof(gint));
    apopdata->iright = (gint *)g_realloc(apopdata->iright, npopin * sizeof(gint));
    apopdata->npopin = npopin;

    g_object_unref(xavg);
    g_object_unref(dx);

    if (indices != NULL) {
        g_free(indices);
    }

    if (ngrind != NULL) {
        g_free(ngrind);
    }

    if (groups != NULL) {
        for (i = 0 ; i < nind; i++) {
            g_free(groups[i]);
        }

        g_free(groups);
    }
}

gchar* apopin_export_data(const Apopindata *apopdata, const FDdata *fddata, enum ExportFormat exportformat)
{
    gint i, j;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#Popin detection  method A \n");
    g_string_append_printf(buf, "#Input parameters:\n");
    g_string_append_printf(buf, "#moving average window width  %d \n", apopdata->nmove);
    g_string_append_printf(buf, "#threshhold of depth derivative for identification of popin  %g \n", apopdata->thresh);
    g_string_append_printf(buf, "#threshhold of depth derivative for determination of popin width/height  %g \n", apopdata->thresh2);
    g_string_append_printf(buf, "#minimum popin width %d \n", apopdata->wpop);
    g_string_append_printf(buf, "#minimum popin height %g \n", apopdata->hpop);
    g_string_append_printf(buf, "\n\n");

    g_string_append_printf(buf, "\n#Popins found by method A : %d \n", apopdata->npopin);

    if (apopdata->npopin) {
        g_string_append_printf(buf, "#    load/mN     depth/mN     depth difference/mN    index left-end      index right-end \n");

        for (i = 0; i < apopdata->npopin; i++) {
            g_string_append_printf(buf, "#%d   %g    %g      %g    %d   %d \n", i + 1,
				   apopdata->Fpopin[i], apopdata->hpopin[i], apopdata->dhpopin[i], apopdata->ileft[i], apopdata->iright[i]);
        }

        g_string_append_printf(buf, "\n\n");

        g_string_append_printf(buf, "#Popin depth data \n");

        for (i = 0; i < apopdata->npopin; i++) {
            g_string_append_printf(buf, "#Popin #%d\n", i + 1);
            g_string_append_printf(buf, "# index i    hload/nm       havg/nm    Fload/nm   dh/di    threshold height   threshold width   \n");

            for (j = apopdata->ileft[i]; j <= apopdata->iright[i]; j++)
                g_string_append_printf(buf, "%g   %g %g   %g   %g %g %g\n",
				       apopdata->t->data[j], fddata->hload->data[j], apopdata->havg->data[j],
				       fddata->Fload->data[j], apopdata->dh->data[j], apopdata->thresh, apopdata->thresh2);

            g_string_append_printf(buf, "\n\n");
        }

        g_string_append_printf(buf, "\n\n");
        g_string_append_printf(buf, "\n\n");
        g_string_append_printf(buf, "# index i    hload/nm       havg/nm    Fload/nm   dh/di    threshold height   threshold width   \n");

        for (j = 0; j < apopdata->havg->res; j++)
            g_string_append_printf(buf, "%g   %g %g   %g   %g %g %g\n",
				   apopdata->t->data[j], fddata->hload->data[j], apopdata->havg->data[j],
				   fddata->Fload->data[j], apopdata->dh->data[j], apopdata->thresh, apopdata->thresh2);

    }

    return g_string_free(buf, FALSE);
}

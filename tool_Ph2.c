#include "tool_Ph2.h"

#include "datatypes.h"
#include "niget-common.h"

#include <math.h>

#include <libprocess/gwyprocess.h>


#define TINY  1e-10


void init_Ph2data(Ph2data *ph2data)
{
    ph2data->nmove = 1;
    ph2data->havg = NULL;
    ph2data->Favg = NULL;
    ph2data->Ph2 = NULL;
    ph2data->dPdh2 = NULL;
    ph2data->ind_remove = NULL;
}

gboolean Ph2_has_results(const Ph2data *ph2data)
{
    return (ph2data->Ph2 != NULL && ph2data->dPdh2 != NULL);
}

GwyDataLine *copy_data_remove_ind(GwyDataLine *source, GwyDataLine *indices)
{
    gint i, n, k;
    gint nremove;
    gdouble *gt;
    gint remain;
    gboolean in;
    GwyDataLine *dest;

    n = source->res;
    nremove = indices->res;

    if (nremove > n) {
        return NULL;
    }

    if (nremove == 1 && indices->data[0] == 0) {
        dest = gwy_data_line_part_extract(source, 1, n - 1);
    }
    else {
        dest = gwy_data_line_new(n - nremove, 1.0, TRUE);
        gt = gwy_data_line_get_data(dest);
        remain = nremove;

        for (i = 0; i < n ; i++) {
            in = FALSE;

            if (remain > 0) {
                in = FALSE;

                for (k = 0; k < nremove; k++) {
                    if (indices->data[k] == i) {
                        in = TRUE;
                    }
                }
            }

            if (in) {
                remain -= 1;
            }
            else {
                *gt = source->data[i];
                gt++;
            }
        }
    }

    return dest;
}

void Ph2_calc(const FDdata *fddata, Ph2data *ph2data)
{
    GwyDataLine *A;
    gdouble dF, dA;
    gint i, n;
    gdouble *ind_remove;
    gint num_remove = 0;
    gint k;
    gboolean in;
    gdouble *gt;

    n = fddata->hload->res;

    /* clear previous data if present */
    if (ph2data->havg != NULL) {
        g_object_unref(ph2data->havg);
    }

    if (ph2data->Favg != NULL) {
        g_object_unref(ph2data->Favg);
    }

    if (ph2data->Ph2 != NULL) {
        g_object_unref(ph2data->Ph2);
    }

    if (ph2data->dPdh2 != NULL) {
        g_object_unref(ph2data->dPdh2);
    }

    //calculate moving average (window with equal weights)
    ph2data->havg = moving_average(fddata->hload, ph2data->nmove);
    ph2data->Favg = moving_average(fddata->Fload, ph2data->nmove);

    //prepare data structures for P-h^2 and dP/d(h^2) curves
    ph2data->Ph2 = gwy_data_line_new_alike(fddata->hload, FALSE);
    ph2data->dPdh2 = gwy_data_line_new_alike(fddata->hload, FALSE);
    A = gwy_data_line_new_alike(fddata->hload, FALSE); /* A = h^2 */

    // it should not happen that A (= h^2) or dA (= dh^2) are zero, but we check and remove the corresponding data
    ind_remove = (gdouble *)g_malloc(n * sizeof(gdouble));
    ind_remove[0] = 0;
    num_remove = 1;

    for (i = 1; i < n; i++) {
        A->data[i] = ph2data->havg->data[i] * ph2data->havg->data[i];
        ph2data->Ph2->data[i] = ph2data->Favg->data[i] / A->data[i];

        if (ph2data->havg->data[i] < 0) {
            ind_remove[num_remove] = i;
            num_remove++;
        }
    }

    for (i = 0; i < n; i++) {
        dF = gwy_data_line_get_der(ph2data->Favg, i);
        dA = gwy_data_line_get_der(A, i);
        ph2data->dPdh2->data[i] = dF / dA;

        if (fabs(dA) < TINY) {
            in = FALSE;

            for (k = 0; k < n; k++)
                if (ind_remove[k] == i) {
                    in = TRUE;
                }

            if (!in) {
                ind_remove[num_remove] = i;
                num_remove++;
            }
        }
    }

    gwy_data_line_multiply(ph2data->dPdh2, 1e9); //TODO units MPa
    gwy_data_line_multiply(ph2data->Ph2, 1e9);//TODO units MPa

    if (ph2data->ind_remove != NULL) {
        g_object_unref(ph2data->ind_remove);
    }

    ph2data->ind_remove = gwy_data_line_new(num_remove, 1.0, TRUE);
    gt = gwy_data_line_get_data(ph2data->ind_remove);

    for (k = 0 ; k < num_remove; k++, gt++) {
        *gt = ind_remove[k];
    }

    g_free(ind_remove);
    g_object_unref(A);
}

gchar* Ph2_export_data(const Ph2data *ph2data, const FDdata *fddata, enum ExportFormat exportformat)
{
    gint i, k;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);
    g_string_append_printf(buf, "#raw data, smoothed depth and load data, P/h^2 and dP/d(h^2) \n");
    g_string_append_printf(buf, "#length of moving average window:  %d \n", ph2data->nmove);
    g_string_append_printf(buf, "\n\n");

    //Results

    if (ph2data->Ph2 != NULL) {
        g_string_append_printf(buf, "#h/nm  F/mN   havg/nm   Favg/mN    P/h^2/MPa   dP/dh^2/MPa \n");

        for (i = 0; i < ph2data->Ph2->res ; i ++) {
            for (k = 0; k < ph2data->ind_remove->res; k++)
                if (i == ph2data->ind_remove->data[k]) {
                    g_string_append_printf(buf, "#");
                }

            g_string_append_printf(buf, "%g %g %g %g %g %g \n",
				   fddata->hload->data[i], fddata->Fload->data[i], ph2data->havg->data[i], ph2data->Favg->data[i],
				   ph2data->Ph2->data[i], ph2data->dPdh2->data[i]);
        }
    }

    return g_string_free(buf, FALSE);
}

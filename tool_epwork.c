#include "tool_epwork.h"
#include "tool_epwork_unc.h"

#include "datatypes.h"
#include "file-utils.h"
#include "niget-common.h"

#include <math.h>

void init_Workdata(Workdata *workdata)
{
    workdata->nmove = 1;
    workdata->worktotal = 0;
    workdata->workelastic = 0;
    workdata->workplastic = 0;
    workdata->eta = 0;
    workdata->hloadavg = NULL;
    workdata->hunloadavg = NULL;
    workdata->hholdavg = NULL;
    workdata->Floadavg = NULL;
    workdata->Funloadavg = NULL;
    workdata->Fholdavg = NULL;
}

gboolean work_has_results(const Workdata *workdata)
{
    return (workdata->hloadavg != NULL && workdata->Floadavg != NULL &&
            workdata->hholdavg != NULL && workdata->Fholdavg != NULL &&
            workdata->hunloadavg != NULL && workdata->Funloadavg != NULL);
}

void work_calc(Workdata *workdata, const FDdata *fddata)
{
    gint i;
    gdouble wload, whold, wunload;

    //smooth loading curve and integrate it
    //calculate moving average (window with equal weights)
    workdata->hloadavg = moving_average(fddata->hload, workdata->nmove);
    workdata->Floadavg = moving_average(fddata->Fload, workdata->nmove);

    wload = 0;

    for (i = 1; i < fddata->hload->res; i++) {
        wload += 0.5 * (workdata->Floadavg->data[i] + workdata->Floadavg->data[i - 1]) * (workdata->hloadavg->data[i] - workdata->hloadavg->data[i - 1]);
    }

    //smooth holding curve and integrate it
    //calculate moving average (window with equal weights)
    workdata->hholdavg = moving_average(fddata->hhold, workdata->nmove);
    workdata->Fholdavg = moving_average(fddata->Fhold, workdata->nmove);

    whold = 0;

    for (i = 1; i < fddata->hhold->res; i++) {
        whold += 0.5 * (workdata->Fholdavg->data[i] + workdata->Fholdavg->data[i - 1]) * (workdata->hholdavg->data[i] - workdata->hholdavg->data[i - 1]);
    }

    //smooth unloading curve and integrate it
    //calculate moving average (window with equal weights)
    workdata->hunloadavg = moving_average(fddata->hunload, workdata->nmove);
    workdata->Funloadavg = moving_average(fddata->Funload, workdata->nmove);

    wunload = 0;

    for (i = 1; i < fddata->hunload->res; i++) {
        wunload += 0.5 * (workdata->Funloadavg->data[i] + workdata->Funloadavg->data[i - 1]) * (workdata->hunloadavg->data[i] - workdata->hunloadavg->data[i - 1]);
    }

    //units are pJ
    //mN * nm = 10^{-3} 10^{-9} Nm = pJ
    workdata->workelastic = fabs(wunload);
    workdata->workplastic = fabs(wload) + fabs(whold) - fabs(wunload);
    workdata->eta = workdata->workelastic / (workdata->workelastic + workdata->workplastic) * 100;
}

void work_remove_results(Workdata *workdata)
{
    // remove all results
    if (workdata->Floadavg != NULL) {
        g_object_unref(workdata->Floadavg);
        g_object_unref(workdata->Fholdavg);
        g_object_unref(workdata->Funloadavg);
        g_object_unref(workdata->hloadavg);
        g_object_unref(workdata->hholdavg);
        g_object_unref(workdata->hunloadavg);
        workdata->Floadavg = NULL;
        workdata->Fholdavg = NULL;
        workdata->Funloadavg = NULL;
        workdata->hloadavg = NULL;
        workdata->hholdavg = NULL;
        workdata->hunloadavg = NULL;
    }

    workdata->eta = 0;
    workdata->workelastic = 0;
    workdata->workplastic = 0;
}

gchar* work_export_data(const Workdata *workdata, const FDdata *fddata, enum ExportFormat exportformat)
{
    GString *buf;
    
    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);

    g_string_append_printf(buf, "#length of moving average window %d \n", workdata->nmove);
    g_string_append_printf(buf, "#elastic deformation work W_e: %g pJ\n", workdata->workelastic);
    g_string_append_printf(buf, "#plastic deformation work W_p: %g pJ\n", workdata->workplastic);
    g_string_append_printf(buf, "#elastic part of total work eta_IT: %g %%\n", workdata->eta);

    return g_string_free(buf, FALSE);
}

#include "tool_epwork_unc.h"
#include "tool_epwork.h"

#include "datatypes.h"
#include "file-utils.h"
#include "niget-common.h"

#include <math.h>

void init_WorkUncdata(WorkUncdata *workunc)
{
    workunc->Wp = NULL;
    workunc->We = NULL;
}

void work_run_shift_contact(const Workdata *wdata, const FDdata *fddata, gdouble *We, gdouble *Wp, gint j)
{
    Workdata w;
    FDdata fdc;
    gdouble dh, dF;
    gint iend;

    We[j + CONTACTMAX] = -1;
    Wp[j + CONTACTMAX] = -1;

    //	printf(" #-----------------------------j = %d ------------------------------------\n\n\n", j);
    // create shifted FDdata

    if (fddata->i_contact_load + j < 0) {
        return;
    }

    fdc.i_contact_load = fddata->i_contact_load + j;
    fdc.horig = gwy_data_line_duplicate(fddata->horig);
    fdc.Forig = gwy_data_line_duplicate(fddata->Forig);

    //	printf(" old contact %d new contact %d  j %d \n", fddata->i_contact_load, fdc.i_contact_load, j);

    fdc.hload = gwy_data_line_part_extract(fdc.horig, fdc.i_contact_load, fddata->hload->res - j);
    fdc.Fload = gwy_data_line_part_extract(fdc.Forig, fdc.i_contact_load, fddata->Fload->res - j);
    gwy_data_line_add(fdc.hload, -fdc.horig->data[fdc.i_contact_load]);
    gwy_data_line_add(fdc.Fload, -fdc.Forig->data[fdc.i_contact_load]);

    /*printf(" loading \n");
    printf("old contact \n");
    for (i = 0; i < fddata->hload->res; i++)
    	printf("%d %g \n", i, fddata->hload->data[i]);
    printf("new contact \n");
    for (i = 0; i < fdc.hload->res; i++)
    	printf("%d %g %g \n", i, fdc.hload->data[i], fddata->hload->data[i+j]-fdc.Fall_orig->data[fdc.i_contact_load]);
    	*/

    fdc.hunload = gwy_data_line_duplicate(fddata->hunload);
    fdc.Funload = gwy_data_line_duplicate(fddata->Funload);
    dh = fdc.horig->data[fdc.i_contact_load] - fddata->horig->data[fddata->i_contact_load];
    dF = fdc.Forig->data[fdc.i_contact_load] - fddata->Forig->data[fddata->i_contact_load];
    gwy_data_line_add(fdc.hunload, -dh);
    gwy_data_line_add(fdc.Funload, -dF);
    //	printf("dh %g dF %g \n", dh, dF);
    iend = 0;

    while (iend < fdc.hunload->res && fdc.hunload->data[iend] > 0 && fdc.Funload->data[iend] > 0) {
        iend++;
    }

    //	printf("iend %d res %d iend old %d \n", iend, fdc.hunload->res, fddata->i_end_unload);
    fdc.hunload = gwy_data_line_part_extract(fdc.hunload, 0, iend);
    fdc.Funload = gwy_data_line_part_extract(fdc.Funload, 0, iend);

    /*
    printf(" unloading \n");
    printf("old contact \n");
    for (i = 0; i < fddata->hunload->res; i++)
    	printf("%d %g %g \n", i, fddata->hunload->data[i], fddata->Funload->data[i]);
    printf("new contact \n");
    for (i = 0; i < fdc.hunload->res; i++)
    	printf("%d %g %g %g %g %g\n", i, fdc.hunload->data[i],fddata->hunload->data[i], fddata->hunload->data[i]-dh, fdc.Funload->data[i], fddata->Funload->data[i]);
    	*/

    fdc.hhold = gwy_data_line_duplicate(fddata->hhold);
    fdc.Fhold = gwy_data_line_duplicate(fddata->Fhold);
    gwy_data_line_add(fdc.hhold, -dh);
    gwy_data_line_add(fdc.Fhold, -dF);

    /*
    printf(" holding \n");
    printf("old contact \n");
    for (i = 0; i < fddata->hhold->res; i++)
    	printf("%d %g \n", i, fddata->hhold->data[i]);
    printf("new contact \n");
    for (i = 0; i < fdc.hhold->res; i++)
    	printf("%d %g %g \n", i, fdc.hhold->data[i], fddata->hhold->data[i]);
    	*/

    /*
    if (!create_shifted_FDdata(&fddata->&fdc, j))
    {
    	printf("could not create shifted FDdata\n");
    	return;
    }
    */

    //shifted Workdata, only nmove
    w.nmove = wdata->nmove;

    work_calc(&w, &fdc);

    We[j + CONTACTMAX] = w.workelastic;
    Wp[j + CONTACTMAX] = w.workplastic;
}

gchar* work_uncertainty_export_data(const WorkUncdata *wunc, enum ExportFormat exportformat)
{
    gint j;
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "#Changes of the elastic/plastic work corresponding to changes in the contact point \n");
    g_string_append_printf(buf, "# shift of contact point with respect to the globally chosen contact point (negative - to the left, positive - to the right)   ");
    g_string_append_printf(buf, "#W_elastic/pJ     W_plastic/pJ \n");

    for (j = -CONTACTMAX; j < CONTACTMAX + 1; j++) {
        g_string_append_printf(buf, "# %d   %g   %g \n", j, wunc->We[j + CONTACTMAX], wunc->Wp[j + CONTACTMAX]);
    }

    return g_string_free(buf, FALSE);
}

gboolean work_unc_has_results(const WorkUncdata *workuncdata)
{
    return (workuncdata->We != NULL && workuncdata->Wp != NULL);
}

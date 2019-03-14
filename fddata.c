#include "fddata.h"

#include "datatypes.h"
#include "file-utils.h"
#include "niget-common.h"

#include <string.h>

#include <libprocess/gwyprocess.h>

static void fddata_find_split_points(GwyDataLine *hdata, GwyDataLine *Fdata, gint i_contact_load, gint *i_load_hold, gint *i_hold_unload, gint *i_end_unload);
static void export_pair_datalines(GString *buf, GwyDataLine *x, GwyDataLine *y);

void init_FDdata(FDdata *fddata)
{
    if (verbose) {
        g_print("Initializing FD data\n");
    }

    fddata->horig = NULL;
    fddata->Forig = NULL;
    fddata->hload = NULL;
    fddata->Fload = NULL;
    fddata->hhold = NULL;
    fddata->Fhold = NULL;
    fddata->hunload = NULL;
    fddata->Funload = NULL;

    fddata->time = NULL;
    fddata->htime = NULL;
    fddata->Ftime = NULL;

    fddata->ndata = 0;
    fddata->hc = 0.0;
    fddata->Fc = 0.0;
    fddata->hmax = 0.0;
    fddata->Fmax = 0.0;
    fddata->i_contact_load = -1;
    fddata->i_load_hold = -1;
    fddata->i_hold_unload = -1;
    fddata->i_end_unload = -1;

    fddata->status_contact_load = POINT_STATUS_NOT_SET;
    fddata->status_load_hold = POINT_STATUS_NOT_SET;
    fddata->status_hold_unload = POINT_STATUS_NOT_SET;
    fddata->status_end_unload = POINT_STATUS_NOT_SET;

    fddata->filename = NULL;
}

void clear_hF_data(FDdata *fddata)
{
    dataline_clear(fddata->horig);
    dataline_clear(fddata->Forig);

    dataline_clear(fddata->time);
    dataline_clear(fddata->htime);
    dataline_clear(fddata->Ftime);

    dataline_clear(fddata->hload);
    dataline_clear(fddata->Fload);
    dataline_clear(fddata->hhold);
    dataline_clear(fddata->Fhold);
    dataline_clear(fddata->hunload);
    dataline_clear(fddata->Funload);
}

void destroy_FDdata(FDdata *fddata)
{
    if (verbose) {
        g_print("Freeing FD data\n");
    }

    clear_hF_data(fddata);

    if (fddata->filename) {
        g_free(fddata->filename);
        fddata->filename = NULL;
    }
}

void dataline_clear(GwyDataLine *dataline)
{
    if (G_IS_OBJECT(dataline)) { /* TOTO ASI BUDE PROBLEM, ted uz mozna ne */
        g_object_unref(dataline);

        /*dataline = NULL; */
        /* cppcheck: Assignment of function parameter has no effect outside the function. Did you forget dereferencing it? */
    }
}

void fddata_get_hmax_Fmax(GwyDataLine *h, GwyDataLine *F, gdouble *hmax, gdouble *Fmax)
{
    gint i;

    *Fmax = -G_MAXDOUBLE;
    *hmax = -G_MAXDOUBLE;

    for (i = 0; i < F->res; i++) {
        if (*Fmax < F->data[i]) {
            *Fmax = F->data[i];
            *hmax = h->data[i];
        }
    }
}


/* splitting */

void /* TODO */ fddata_find_contact(GwyDataLine *hdata, GwyDataLine *Fdata, gint *i_contact_load, gdouble *hc, gdouble *Fc)
{
    *i_contact_load = 0;

    *hc = hdata->data[*i_contact_load];
    *Fc = Fdata->data[*i_contact_load];
}

static void fddata_find_split_points(GwyDataLine *hdata, GwyDataLine *Fdata, gint i_contact_load,
                                     gint *i_load_hold, gint *i_hold_unload, gint *i_end_unload)
{
    gint i;
    GwyDataLine *F;
    gdouble hc, hmax, Fmax, uF;

    if (verbose) {
        g_print("Find split points:\n");
    }

    hc = hdata->data[i_contact_load];
    fddata_get_hmax_Fmax(hdata, Fdata, &hmax, &Fmax);

    uF = 0.0;

    if (i_contact_load > 8) { /* kde vzit sumovou cast u Hysitronu? */
        F = gwy_data_line_part_extract(Fdata, 0, i_contact_load);
        uF = gwy_data_line_get_rms(F);
        g_object_unref(F);
    }

    if (uF < DBL_EPSILON) {
        uF = 0.001 * Fmax;
    }

    if (verbose) {
        g_print("  uF: %g\n", uF);
    }

    i = i_contact_load;

    while ((i < hdata->res) && (Fdata->data[i] < (Fmax - 3 * uF))) {
        i++;
    }

    *i_load_hold = i; /* first point very close to Fmax */

    while ((i < hdata->res) && (hdata->data[i] >= hc)) {
        i++;
    }

    i--;
    *i_end_unload = i; /* last available point with nonnegative depth */

    while ((i > (*i_load_hold + 1)) && (Fdata->data[i] < (Fmax - 3 * uF))) {
        i--;
    }

    //i++; //TEST
    *i_hold_unload = i; /* last point very close to Fmax */

    if (verbose) {
        g_print("  split point indices: \n");
        g_print("   contact/load %d \n", i_contact_load);
        g_print("   load/hold %d \n", *i_load_hold);
        g_print("   hold/unload %d \n", *i_hold_unload);
        g_print("   unload/end %d \n", *i_end_unload);
    }
}

/* extract individual parts from orig data */
gboolean split_hF_data_spec(FDdata *fddata, gint i_contact_load, gint i_load_hold, gint i_hold_unload, gint i_end_unload)
{
    gdouble hc, Fc;
    GwyDataLine *hall, *Fall;

    /* check whether split points make sense */
    if ((0 <= i_contact_load) &&
            (i_contact_load <= i_load_hold) &&
            (i_load_hold <= i_hold_unload) &&
            (i_hold_unload <= i_end_unload) &&
            (i_end_unload < fddata->ndata)) {

        hc = fddata->horig->data[i_contact_load];
        Fc = fddata->Forig->data[i_contact_load];

        hall = gwy_data_line_duplicate(fddata->horig);
        Fall = gwy_data_line_duplicate(fddata->Forig);
        gwy_data_line_add(hall, -hc);
        gwy_data_line_add(Fall, -Fc);

        /* load */
        dataline_clear(fddata->hload);
        dataline_clear(fddata->Fload);
        fddata->hload = gwy_data_line_part_extract(hall, i_contact_load, i_load_hold - i_contact_load + 1);
        fddata->Fload = gwy_data_line_part_extract(Fall, i_contact_load, i_load_hold - i_contact_load + 1);
        gwy_data_line_set_real(fddata->hload, fddata->hload->res);
        gwy_data_line_set_real(fddata->Fload, fddata->Fload->res);

        /* hold */
        dataline_clear(fddata->hhold);
        dataline_clear(fddata->Fhold);
        fddata->hhold = gwy_data_line_part_extract(hall, i_load_hold, i_hold_unload - i_load_hold + 1);
        fddata->Fhold = gwy_data_line_part_extract(Fall, i_load_hold, i_hold_unload - i_load_hold + 1);
        gwy_data_line_set_real(fddata->hhold, fddata->hhold->res);
        gwy_data_line_set_real(fddata->Fhold, fddata->Fhold->res);

        /* unload */
        dataline_clear(fddata->hunload);
        dataline_clear(fddata->Funload);
        fddata->hunload = gwy_data_line_part_extract(hall, i_hold_unload, i_end_unload - i_hold_unload + 1);
        fddata->Funload = gwy_data_line_part_extract(Fall, i_hold_unload, i_end_unload - i_hold_unload + 1);
        gwy_data_line_set_real(fddata->hunload, fddata->hunload->res);
        gwy_data_line_set_real(fddata->Funload, fddata->Funload->res);

        dataline_clear(hall);
        dataline_clear(Fall);

        fddata->Fc = Fc;

        fddata->i_contact_load = i_contact_load;
        fddata->i_load_hold = i_load_hold;
        fddata->i_hold_unload = i_hold_unload;
        fddata->i_end_unload = i_end_unload;

        /* maximum values only from unloading curve to match CSM behavior */
        fddata_get_hmax_Fmax(fddata->hunload, fddata->Funload, &(fddata->hmax), &(fddata->Fmax));

        return TRUE;
    }
    else {
        return FALSE;
    }
}

gboolean split_hF_data_auto(FDdata *fddata)
{
    gint i_contact_load, i_load_hold, i_hold_unload, i_end_unload;
    gboolean isok, splitok;

    i_contact_load = fddata->i_contact_load;

    fddata_find_split_points(fddata->horig, fddata->Forig, i_contact_load, &i_load_hold, &i_hold_unload, &i_end_unload);

    fddata->status_load_hold = POINT_STATUS_SET_AUTO;
    fddata->status_hold_unload = POINT_STATUS_SET_AUTO;
    fddata->status_end_unload = POINT_STATUS_SET_AUTO;

    isok = TRUE;

    if (i_load_hold <= i_contact_load) {
        g_printerr("Warning: Could not split data between load and hold part.\n Proceed with care.\n");
        i_load_hold = i_contact_load;
        isok = FALSE;
    }

    if (i_hold_unload < i_load_hold) {
        g_printerr("Warning: Could not split data between hold and unload part.\n Proceed with care.\n");
        i_hold_unload = i_load_hold;
        isok = FALSE;
    }

    if (i_end_unload <= i_hold_unload) {
        g_printerr("Warning: Could not split data between unload and end part.\n Proceed with care.\n");
        i_end_unload = i_hold_unload;
        isok = FALSE;
    }

    splitok = split_hF_data_spec(fddata, i_contact_load, i_load_hold, i_hold_unload, i_end_unload);
    isok = isok && splitok;

    return isok;
}

/* shifting */

gdouble create_shifted_FDdata(const FDdata *src, FDdata *dest, gint j)
{
    gdouble dh;
    /* gdouble dF; */

    init_FDdata(dest);

    dest->i_contact_load = src->i_contact_load + j;
    dest->i_load_hold = src->i_load_hold;
    dest->i_hold_unload = src->i_hold_unload;
    dest->i_end_unload = src->i_end_unload;

    dest->horig = gwy_data_line_duplicate(src->horig);
    dest->Forig = gwy_data_line_duplicate(src->Forig);

    dest->ndata = src->ndata;

    split_hF_data_spec(dest, dest->i_contact_load, dest->i_load_hold, dest->i_hold_unload, dest->i_end_unload);

    //difference in shift between different contactpoints
    dh = dest->horig->data[dest->i_contact_load] - src->horig->data[src->i_contact_load];
    /* dF = dest->Forig->data[dest->i_contact_load] - src->Forig->data[src->i_contact_load]; */

    /*
    FILE *F;
    gchar str[100];
    gint i;

    sprintf(str, "create_shifted_FDdata_new_%d",j);
    F= fopen(str, "w");
    for (i=0 ; i< dest->hload->res; i++)
    	fprintf(F, "%g %g \n", dest->hload->data[i], dest->Fload->data[i]);
    fprintf(F,"\n\n");
    for (i=0 ; i< dest->hhold->res; i++)
    	fprintf(F, "%g %g \n", dest->hhold->data[i], dest->Fhold->data[i]);
    fprintf(F,"\n\n");
    for (i=0 ; i< dest->hunload->res; i++)
    	fprintf(F, "%g %g \n", dest->hunload->data[i], dest->Funload->data[i]);
    fprintf(F,"\n\n");
    fclose(F);
    */

    return dh;
}

FileReadStatus /* find a right place for this function */
load_data_nogui(const gchar *filename, enum DataFileFormat fileformat, FDdata *fddata, gboolean *splitok) {
    FileReadStatus readstatus;
    gint i, numdata, i_contact_load, i_load_hold, i_hold_unload, i_end_unload;
    FILE *infile;
    GwyDataLine *timedata, *hdata, *Fdata;

    /* we should check if fddata is not NULL ... */
    /* will be avoided when this function actually creates the fddata */
    /* so now we should clean up existing fddata to prevent leaks */
    /* ale data bych mel likvidovat az tehdy, kdyz vim, ze se neco fakt uspesne nacetlo */

    if (verbose) {
        g_print("Loading data from %s\n", filename);
    }

    infile = fopen(filename, "r");

    if (!infile) {
        return READ_FILE_COULD_NOT_OPEN;
    }

    numdata = num_data_in_file(infile, fileformat);

    timedata = gwy_data_line_new(numdata, numdata, FALSE);
    hdata = gwy_data_line_new(numdata, 1.0, FALSE);
    Fdata = gwy_data_line_new(numdata, 1.0, FALSE);

    switch (fileformat) { /* tady pridavat nove formaty */
    case DATA_FILE_FORMAT_NIGET:
        readstatus = read_file_niget2(infile, numdata, timedata, hdata, Fdata,
                                      &i_contact_load, &i_load_hold, &i_hold_unload, &i_end_unload);
        break;

    default:
        readstatus = read_file_plaintext2(infile, fileformat, numdata, timedata, hdata, Fdata);
        break;
    }

    fclose(infile);

    if (readstatus != READ_FILE_OK) { /* the right number of entries were not read */
        return readstatus;
    }

    /* check if data is not inf or nan */
    for (i = 0; i < numdata; i++) {
        if (isnan(hdata->data[i]) || isnan(Fdata->data[i])) {
            dataline_clear(timedata);
            dataline_clear(hdata);
            dataline_clear(Fdata);
            return READ_FILE_CONTAINS_NAN;
        }

        if (isinf(hdata->data[i]) || isinf(Fdata->data[i])) {
            dataline_clear(timedata);
            dataline_clear(hdata);
            dataline_clear(Fdata);
            return READ_FILE_CONTAINS_INF;
        }
    }

    /* everything OK, populate fddata (check whether all fields are complete!) */
    clear_hF_data(fddata);
    init_FDdata(fddata);

    fddata->ndata = numdata;
    fddata->time = gwy_data_line_duplicate(timedata);
    fddata->horig = gwy_data_line_duplicate(hdata);
    fddata->Forig = gwy_data_line_duplicate(Fdata);

    fddata->filename = strdup(filename);
    fddata->fileformat = fileformat;

    /* split */
    switch (fileformat) {
    case DATA_FILE_FORMAT_NIGET:

        fddata->i_contact_load = i_contact_load;
        fddata->i_load_hold = i_load_hold;
        fddata->i_hold_unload = i_hold_unload;
        fddata->i_end_unload = i_end_unload;

        *splitok = split_hF_data_spec(fddata,
                                      fddata->i_contact_load, fddata->i_load_hold,
                                      fddata->i_hold_unload, fddata->i_end_unload);

        fddata->status_contact_load = POINT_STATUS_FROM_FILE;
        fddata->status_load_hold = POINT_STATUS_FROM_FILE;
        fddata->status_hold_unload = POINT_STATUS_FROM_FILE;
        fddata->status_end_unload = POINT_STATUS_FROM_FILE;
        break;

    default: /* standard data formats */
        fddata->i_contact_load = 0; // musi se vynulovat i_load_hold i_hold_unload?
        fddata->i_load_hold = 0;
        fddata->i_hold_unload = 0;
        fddata->i_end_unload = numdata - 1;

        fddata_find_contact(fddata->horig, fddata->Forig, &(fddata->i_contact_load), &(fddata->hc), &(fddata->Fc));
        fddata->status_contact_load = POINT_STATUS_SET_AUTO;

        *splitok = split_hF_data_auto(fddata);
        break;
    }

    fddata->hc = hdata->data[fddata->i_contact_load];
    fddata->Fc = Fdata->data[fddata->i_contact_load];

    return READ_FILE_OK;
}

static void export_pair_datalines(GString *buf, GwyDataLine *x, GwyDataLine *y)
{
    gint i;
    
    if (x->res == y->res) {
	for (i = 0; i < x->res; i++) {
	    g_string_append_printf(buf, "%f\t%f\n", x->data[i], y->data[i]);
	}
    }

    g_string_append(buf, "\n\n");
}

gchar* fddata_export_data(const FDdata *fddata, enum ExportFormat exportformat)
{
    GString *buf;

    buf = g_string_new(NULL);

    g_string_append_printf(buf, "# %s \n", fddata->filename);

    if (fddata->horig && fddata->Forig) {
        g_string_append_printf(buf, "#fddata: horig, Forig\n");
	export_pair_datalines(buf, fddata->horig, fddata->Forig);
    }

    if (fddata->hload && fddata->Fload) {
        g_string_append_printf(buf, "#fddata: hload, Fload\n");
	export_pair_datalines(buf, fddata->hload, fddata->Fload);
    }

    if (fddata->hhold && fddata->Fhold) {
        g_string_append_printf(buf, "#fddata: hhold, Fhold\n");
	export_pair_datalines(buf, fddata->hhold, fddata->Fhold);
    }

    if (fddata->hunload && fddata->Funload) {
        g_string_append_printf(buf, "#fddata: hunload, Funload\n");
	export_pair_datalines(buf, fddata->hunload, fddata->Funload);
    }

    return g_string_free(buf, FALSE);
}

/* ================ */

# if 0
//TODO postupne smazat, zatim jeste nechat pro jistotu
gdouble create_shifted_FDdata_old(const FDdata *src, FDdata *dest, gint j)
{
    gdouble dh, dF;

    init_FDdata(dest);
    dest->i_contact_load = src->i_contact_load + j;

    dest->hall_orig = gwy_data_line_duplicate(src->hall_orig);
    dest->Fall_orig = gwy_data_line_duplicate(src->Fall_orig);
    dest->hall = gwy_data_line_duplicate(src->hall_orig);
    dest->Fall = gwy_data_line_duplicate(src->Fall_orig);

    dest->hload = gwy_data_line_duplicate(src->hload);
    dest->Fload = gwy_data_line_duplicate(src->Fload);
    dest->hhold = gwy_data_line_duplicate(src->hhold);
    dest->Fhold = gwy_data_line_duplicate(src->Fhold);
    dest->hunload = gwy_data_line_duplicate(src->hunload);
    dest->Funload = gwy_data_line_duplicate(src->Funload);

    //difference in shift between different contactpoints
    dh = dest->hall_orig->data[dest->i_contact_load] - src->hall_orig->data[src->i_contact_load];
    dF = dest->Fall_orig->data[dest->i_contact_load] - src->Fall_orig->data[src->i_contact_load];
    \

    if (verbose) {
        g_print("Create shifted FDdata OLD\n");
        g_print("dh %g dF %g \n", dh, dF);
    }

    gwy_data_line_add(dest->hload, -dh);
    gwy_data_line_add(dest->Fload, -dF);
    gwy_data_line_add(dest->hhold, -dh);
    gwy_data_line_add(dest->Fhold, -dF);
    gwy_data_line_add(dest->hunload, -dh);
    gwy_data_line_add(dest->Funload, -dF);

    dest->hmax = src->hmax - dh;
    dest->Fmax = src->Fmax - dF;

    FILE *F;
    gchar str[100];
    gint i;

    sprintf(str, "create_shifted_FDdata_old_%d", j);
    F = fopen(str, "w");

    for (i = 0 ; i < dest->hload->res; i++) {
        fprintf(F, "%g %g \n", dest->hload->data[i], dest->Fload->data[i]);
    }

    fprintf(F, "\n\n");

    for (i = 0 ; i < dest->hhold->res; i++) {
        fprintf(F, "%g %g \n", dest->hhold->data[i], dest->Fhold->data[i]);
    }

    fprintf(F, "\n\n");

    for (i = 0 ; i < dest->hunload->res; i++) {
        fprintf(F, "%g %g \n", dest->hunload->data[i], dest->Funload->data[i]);
    }

    fprintf(F, "\n\n");
    fclose(F);

    return dh;
}
#endif

/*
gboolean
create_shifted_FDdata(FDdata *FDold, FDdata *FDnew, gint j)
{

	gint i;
	gdouble dh,dF;
	// if there are no points far enough to the left, we cannot create the shifted data
	if (FDold->i_contact_load+j <0)
		return FALSE;


	FDnew->hall = gwy_data_line_duplicate(FDold->hall_orig);
	FDnew->Fall = gwy_data_line_duplicate(FDold->Fall_orig);

	FDnew->i_contact_load = FDold->i_contact_load+j;
	FDnew->Fc = FDold->Fc;

	// split_hF_data(FDnew);
	dh = FDnew->hall->data[FDnew->i_contact_load];
	dF = FDnew->Fall->data[FDnew->i_contact_load];

	gwy_data_line_add(FDnew->hall, -dh);
	gwy_data_line_add(FDnew->Fall, -dF);
	*/

/* refine_split_hF_data(FDnew); */
/*
	split_hF_data(FDnew);

	printf("\n\n\n ------------------Loading data-----------------------------\n");
	for (i=0; i< FDnew->hload->res; i++)
	{
		printf("%g %g \n", FDnew->hload->data[i], FDnew->Fload->data[i]);
	}
	printf("\n\n\n ------------------Hold data-----------------------------\n");
	for (i=0; i< FDnew->hhold->res; i++)
	{
		printf("%g %g \n", FDnew->hhold->data[i], FDnew->Fhold->data[i]);
	}
	printf("\n\n\n ------------------Unloading data-----------------------------\n");
	for (i=0; i< FDnew->hunload->res; i++)
	{
		printf("%g %g \n", FDnew->hunload->data[i], FDnew->Funload->data[i]);
	}
	printf(" \n \n \n -------------------------------------------\n");

	printf("prislusne indexy: \n");
	printf(" contact/load %d \n", FDnew->i_contact_load);
	printf(" load/hold %d \n", FDnew->i_load_hold);
	printf(" hold/unload %d \n", FDnew->i_hold_unload);
	printf(" unload/end %d \n", FDnew->i_end_unload);

	return TRUE;
}
*/

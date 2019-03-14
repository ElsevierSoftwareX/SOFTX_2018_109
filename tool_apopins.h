#ifndef TOOL_APOPINS_H
#define TOOL_APOPINS_H

#include "datatypes.h"
#include "file-utils.h"

#define THRESH_DEFAULT 2
#define THRESH2_DEFAULT 0.07
#define WPOP_DEFAULT 5
#define HPOP_DEFAULT 2.0

void init_Apopindata(Apopindata *apopdata);
gchar* apopin_export_data(const Apopindata *apopdata, const FDdata *fddata, enum ExportFormat exportformat);
gboolean apopin_has_results(const Apopindata *apopdata);
void apopin_remove_results(Apopindata *apopdata);
void apopin_calc(Apopindata *apopdata, const FDdata *fddata);
void apopin_create_aux_datalines(Apopindata *apopdata, const FDdata *fddata);

#endif

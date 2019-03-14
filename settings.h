#ifndef SETTINGS_H
#define SETTINGS_H

#include "datatypes.h"

void save_settings(const gchar *filename, const Args *args);
void load_settings(const gchar *filename, Args *args);

#endif

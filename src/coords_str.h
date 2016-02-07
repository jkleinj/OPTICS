/*==============================================================================
coords_str.h: object representation and related functions
Copyright (C) 2008-2015 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
==============================================================================*/

#ifndef COORDS_STR_H
#define COORDS_STR_H

#include <string.h>

#include "arg.h"
#include "config.h"
#include "safe.h"

/* input data */
typedef struct {
    char label[64];
    char string[512];
} String;

typedef struct {
	String *data; /* input data to order */
	int nData; /* number of data points */
} Dat;

/*____________________________________________________________________________*/
/* prototypes */
int get_data(char *inFileName, Dat *dat);
float calc_dist(Dat *dat, int i, int j, Arg *arg);
void print_header_object(FILE *outfile);
void print_object(FILE *outfile, Dat *dat, int index, int order, int cluster_id, float coreDist, float reachDist);

#endif

/*==============================================================================
coords_vec.h: object representation and related functions
Copyright (C) 2008-2015 Jens Kleinjung and Alessandro Pandini 
Read the COPYING file for license information.
==============================================================================*/

#ifndef COORDS_VEC_H
#define COORDS_VEC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "arg.h"
#include "config.h"
#include "safe.h"

/* input data */
typedef struct {
	float (*data)[512]; /* input data to order */
	int nData; /* number of data points */
	int lData; /* length of data vector */
} Dat;

/*____________________________________________________________________________*/
/* prototypes */
int get_data(char *inFileName, Dat *dat);
float calc_dist_asym(Dat *dat, int i, int j, Arg *arg);
float calc_dist_sym(Dat *dat, int i, int j, Arg *arg);
float calc_dist(Dat *dat, int i, int j, Arg *arg);
float calc_dist_ndim(Dat *dat, int i, int j, Arg *arg);
void print_header_object(FILE *outfile);
void print_object(FILE *outfile, Dat *dat, int index, int order, int cluster_id, float coreDist, float reachDist);

#endif

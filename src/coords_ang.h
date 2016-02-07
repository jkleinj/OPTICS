/*==============================================================================
coords_ang.h: object representation and related functions
Copyright (C) 2008-2015 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
==============================================================================*/

#ifndef COORDS_ANG_H
#define COORDS_ANG_H

#include <math.h>
#include <pcre.h>
#include <string.h>

#include "arg.h"
#include "config.h"
#include "safe.h"

/* input data */
typedef struct {
    float phi1; /* angle phi1 CA1 - CA2 - CA3 */
    float phi2; /* angle phi2 CA2 - CA3 - CA4 */
    float theta; /* torsion angle theta CA1 - CA2 - CA3 - CA4 */
} Ang;

typedef struct {
	Ang *data; /* input data to order */
	int nData; /* number of data points */
} Dat;

/*____________________________________________________________________________*/
/* prototypes */
int get_data(char *inFileName, Dat *dat);
float calc_dist(Dat *dat, int i, int j, Arg *arg);
void print_header_object(FILE *outfile);
void print_object(FILE *outfile, Dat *dat, int index, int order, int cluster_id, float coreDist, float reachDist);

#endif

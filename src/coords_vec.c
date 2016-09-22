/*==============================================================================
coords_vec.h : object representation and related functions
Copyright (C) 2008-2015 Jens Kleinjung and Alessandro Pandini 
Read the COPYING file for license information.
==============================================================================*/

#include "coords_vec.h"

/* define distance metric: KLD, SKLD, SHD, CTD, MHD */
#define MHD

/*____________________________________________________________________________*/
/* print header for point */
void print_header_object(FILE *outfile) {

    fprintf(outfile, "\tindex\torder\tclusterId\tCD\tRD\t");
}

/*____________________________________________________________________________*/
/* print data point */
void print_object(FILE *outfile, Dat *dat, int index, int order, int cluster_id, float coreDist, float reachDist) {
	unsigned int k;

    fprintf(outfile, "%10d %10d %10d %8.3f %8.3f ", 
            index, order, cluster_id,
            coreDist, reachDist);
	for (k = 0; k < dat->lData; ++ k)
		fprintf(outfile, "%8.3f ", dat->data[index][k]); 
    fprintf(outfile, "\n"); 
}

/*____________________________________________________________________________*/
/* calculate Kullback-Leibler distance (KLD) */
#ifdef KLD
float calc_dist(Dat *dat, int i, int j, Arg *arg) {
	unsigned int k;
    float dist = 0.;

	for (k = 0; k < dat->lData; ++ k)
		if ((dat->data[i][k] > 0.) && (dat->data[j][k] > 0.))
            dist += dat->data[i][k] * log(dat->data[i][k] / dat->data[j][k]);
    return dist;
}
#endif

/*____________________________________________________________________________*/
/* calculate symmetrised Kullback-Leibler distance (SKLD) */
#ifdef SKLD
float calc_dist(Dat *dat, int i, int j, Arg *arg) {
	unsigned int k;
    float dist = 0.;

	for (k = 0; k < dat->lData; ++ k)
		if ((dat->data[i][k] > 0.) && (dat->data[j][k] > 0.))
            dist += (.5 * dat->data[i][k] * log(dat->data[i][k] / dat->data[j][k])) +\
					(.5 * dat->data[j][k] * log(dat->data[j][k] / dat->data[i][k]));
    return dist;
}
#endif

/*____________________________________________________________________________*/
/* calculate SH distance (SHD) */
#ifdef SHD
float calc_dist(Dat *dat, int i, int j, Arg *arg) {
	unsigned int k;
    float dist = 0.;
	float s_x = 0.;
	float s_y = 0.;
	float s_xy = 0.;

	for (k = 0; k < dat->lData; ++ k) {
		if (dat->data[i][k] > 0)
			s_x += dat->data[i][k] * log2(dat->data[i][k]);
		if (dat->data[j][k] > 0)
			s_y += dat->data[j][k] * log2(dat->data[j][k]);
		if ((dat->data[i][k] > 0) && (dat->data[j][k] > 0))
			s_xy += (dat->data[i][k] + dat->data[j][k]) * log2(dat->data[i][k] + dat->data[j][k]);
	}

	dist = .5 * (s_xy - (s_x + s_y));

	return dist;
}
#endif

/*____________________________________________________________________________*/
/* calculate n-dimensional Cartesian distance (CTD)*/
#ifdef CTD
float calc_dist(Dat *dat, int i, int j, Arg *arg) {
	unsigned int k;
    float dist = 0.;
	float squaresum = 0.;

	for (k = 0; k < dat->lData; ++ k)
		squaresum += pow(dat->data[i][k] - dat->data[j][k], 2);

	dist = sqrt(squaresum / dat->lData);

    return dist;
}
#endif

/*____________________________________________________________________________*/
/* calculate n-dimensional Manhattan distance (MHD) */
#ifdef MHD
float calc_dist(Dat *dat, int i, int j, Arg *arg) {
	unsigned int k;
    float dist = 0.;

	for (k = 0; k < dat->lData; ++ k)
		dist += abs(dat->data[i][k] - dat->data[j][k]);

    return dist;
}
#endif

/*----------------------------------------------------------------------------*/
__inline__ int approximately_equal(float a, float b)
{
	return (int)(a == b || fabsf(a - b) <= .00001F * fabsf(a) + .00001F);
}

/*____________________________________________________________________________*/
/* read data points */
/* Data points are expected to be an array of vectors of arbitrary length;
	note that the maximal vector length can be adjusted (see below).
	The vectors must be of equal length. */
int get_data(char *inFileName, Dat *dat)
{
	FILE *inFile = 0;
	unsigned int allocated = 64;
	char line[512] = "";
	char cpline[512] = "";
	char *pch = 0;
	unsigned int k;
	/* allocate memory */
	/* The maximal vector length is defined in *data[512] 
		in the coord_vec.h header. If you modify the 512
		to a 2048 (for example), then set maxDim also to 2048. */
	const int maxDim = 512;
	dat->data = safe_malloc(allocated * sizeof(float [maxDim]));
	pch = safe_malloc(allocated * sizeof(char));

	/* read input data */
	inFile = safe_open(inFileName, "r");

	/* process input data */
	dat->nData = 0;
	while(fgets(line, maxDim, inFile) != NULL) {
		k = 0;
		/* read line */
#ifdef DEBUG
		fprintf(stderr, "%s", line);
#endif
		strcpy(&(cpline[0]), &(line[0]));
		/* split line and assign first vector element */
		pch = strtok(cpline, " ");
		dat->data[dat->nData][k ++] = atof(pch);
		/* assign remaining vector elements */
		while (pch != NULL) {
			pch = (strtok(NULL, " "));
			if (pch != NULL) {
				dat->data[dat->nData][k ++] = atof(pch); 
			}
			/* check for absolute vector length */
			assert((k < maxDim) && "Increase *data[] array in coord_vec.h!\n");
        }

		/* check for relative vector length */
		if (dat->nData == 0) {
			dat->lData = k;
		} else {
			assert((k == dat->lData) && "Input vectors must be of equal length!\n");
		}

		++ dat->nData;

		/* allocate additional memory */
		if (dat->nData == allocated) {
			allocated += 64;
			dat->data = safe_realloc(dat->data, allocated * sizeof(float [maxDim]));
		}
	}

	fclose(inFile);
	free(pch);
	return 0;
}


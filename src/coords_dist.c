/*==============================================================================
coords_dist.h : object representation and related functions
Copyright (C) 2008-2015 Jens Kleinjung and Alessandro Pandini 
Read the COPYING file for license information.
==============================================================================*/

#include "coords_dist.h"

/*____________________________________________________________________________*/
/* print header for point */
void print_header_object(FILE *outfile){
    fprintf(outfile, "     index      order  clusterId       CD       RD        dist\n");
}

/*____________________________________________________________________________*/
/* print data point */
void print_object(FILE *outfile, Dat *dat, int index, int order, int cluster_id, float coreDist, float reachDist){
    fprintf(outfile, "%10d %10d %10d %8.3f %8.3f\n", 
            index, order, cluster_id,
            coreDist, reachDist); 
}

/*____________________________________________________________________________*/
/* calculate distance */
float calc_dist(Dat *dat, int i, int j, Arg *arg){
    return dat->dist[i][j];
}

/*____________________________________________________________________________*/
/* read distances */
/* distances are expected to be in the form of a triangular matrix of floats */
int get_data(char *inFileName, Dat *dat)
{
	FILE *inFile = 0;
	unsigned int nd = 0; /* number of data points */
	unsigned int ndsym = 0;
	float dummy = 0.;
	float dist = 0.; /* point distance */
	unsigned int dim; /* dimension of input matrix */
	unsigned int i = 0;
	unsigned int j = 0;

	inFile = safe_open(inFileName, "r");

	/* count number of data points  */
	while(! feof(inFile)) {
		if (fscanf(inFile, "%f", &dummy) == 1) {
			++ nd;
		}
	}
	assert(nd > 1);
	fprintf(stdout, "%d distances in triangular input matrix\n", nd);

	/* assume triangular matrix of size n*(n-1)/2 */
	/* approximate square */
	ndsym = (2 * nd) + sqrt(2 * nd);
	dim = ceil(sqrt(ndsym));
	assert(nd == dim * (dim - 1) / 2 && "triangular matrix must be n * (n-1) / 2");

	/* allocate memory */
	dat->nData = dim;
	dat->dist = alloc_mat2D_float(dat->dist, dim, dim);
	init_mat2D_float(dat->dist, dim, dim, 0.);

	/* read distances */
	rewind(inFile);

	while(! feof(inFile)) {
		/* scan distance data */
		if (fscanf(inFile, "%f", &dist) == 1) {
			/* advance column */
			++ j;
			/* advance row and reset column for new row */
			if (j == dim) {
				++ i;
				j = i + 1;
			}
			/* assign distance values to symmetric matrix */
			dat->dist[i][j] = dist; 
			dat->dist[j][i] = dist;
		}
	}

	fclose(inFile);

#ifdef DEBUG
	for (i = 0; i < dim; ++ i) {
		for (j = 0; j < dim; ++ j) {
			fprintf(stderr, "%f ", dat->dist[i][j]);
			if (j == (dim - 1)) {
				fprintf(stderr, "\n");
			}
		}
	}
#endif

	return 0;
}


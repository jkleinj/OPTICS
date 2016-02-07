/*==============================================================================
coords_xyz.c : object representation and related functions
Copyright (C) 2008-2015 Jens Kleinjung and Alessandro Pandini 
Read the COPYING file for license information.
==============================================================================*/

#include "coords_xyz.h"

/*____________________________________________________________________________*/
/* print header for point */
void print_header_object(FILE *outfile){
    fprintf(outfile, "     index      order  clusterId       CD       RD        x        y        z\n");
}

/*____________________________________________________________________________*/
/* print data point */
void print_object(FILE *outfile, Dat *dat, int index, int order, int cluster_id, float coreDist, float reachDist){
    fprintf(outfile, "%10d %10d %10d %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
            index, order, cluster_id,
            coreDist, reachDist,
            dat->data[index].x, dat->data[index].y, dat->data[index].z); 
}

/*____________________________________________________________________________*/
/* calculate distance */
float calc_dist(Dat *dat, int i, int j, Arg *arg){

    float dist;

    dist = v_rmsd(&(dat->data[i]), &(dat->data[j]));

    return dist;
}

/*____________________________________________________________________________*/
/* read data points */
/* Data points are expected to be an array of vectors in 3D with
 * coordinates 'x y z'. */
int get_data(char *inFileName, Dat *dat)
{
	FILE *inFile = 0;
	unsigned int allocated = 64;

	/* allocate memory */
	dat->data = safe_malloc(allocated * sizeof(Vec));

	/* read data */
	inFile = safe_open(inFileName, "r");

	dat->nData = 0;
	while(! feof(inFile)) {
#ifdef DEBUG
		++ n;
		if (fscanf(inFile, "%f%f%f\n",
				&(dat->data[dat->nData].x), &(dat->data[dat->nData].y), &(dat->data[dat->nData].z)) != 3) {
			fprintf(stderr, "input data format has to be: [float] [float] [float]\n");
			fprintf(stderr, "format error in line %d\n", n);
			exit(1);
		}
#else
		assert(fscanf(inFile, "%f%f%f\n",
			&(dat->data[dat->nData].x), &(dat->data[dat->nData].y), &(dat->data[dat->nData].z)) == 3);
#endif
		++ dat->nData;

		if (dat->nData == allocated) {
			allocated += 64;
			dat->data = safe_realloc(dat->data, allocated * sizeof(Vec));
		}
	}

	assert(dat->nData > 1);

	fclose(inFile);
	return 0;
}


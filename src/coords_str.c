/*==============================================================================
coords_str.c: object representation and related functions
Copyright (C) 2008-2015 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
==============================================================================*/

#include "coords_str.h"

/*____________________________________________________________________________*/
/* print header for point */
void print_header_object(FILE *outfile){
    fprintf(outfile, "     index      order  clusterId       CD       RD string\n");
}

/*____________________________________________________________________________*/
/* print data point */
void print_object(FILE *outfile, Dat *dat, int index, int order, int cluster_id, float coreDist, float reachDist){
    fprintf(outfile, "%10d %10d %10d %8.3f %8.3f %s\n", 
            index, order, cluster_id,
            coreDist, reachDist,
            dat->data[index].string); 
}

/*____________________________________________________________________________*/
/* calculate distance */
float calc_dist(Dat *dat, int i, int j, Arg *arg){
    
    int dist = 0;
    char *iPtr = 0;
    char *jPtr = 0;

    iPtr = dat->data[i].string;
    jPtr = dat->data[j].string;

    if (strlen(iPtr) != strlen(jPtr)){
        fprintf(stderr, "ERROR: String %d (%d) and %d (%d) differ in length!", i, (int)strlen(iPtr), j, (int)strlen(jPtr));
        exit(1);
    }

    while(strlen(iPtr) >= arg->w){
        if (strncmp(iPtr, jPtr, arg->w) != 0)
            dist++;
        ++iPtr;
        ++jPtr;
    }

    return (float)dist;

}

/*____________________________________________________________________________*/
/* read data points */
/* Data points are expected to be an array of strings in a text file
 *     first column: label
 *     second column: data string. */
int get_data(char *inFileName, Dat *dat)
{
	FILE *inFile = 0;
	unsigned int allocated = 64;
	unsigned int n = 0;

	/* read data */
	inFile = safe_open(inFileName, "r");

	/* allocate memory */
	dat->data = safe_malloc(allocated * sizeof(String));

	dat->nData = 0;
	while(! feof(inFile)) {
		++ n;
#ifdef DEBUG
		if (fscanf(inFile, "%s %s\n", dat->data[dat->nData].label, dat->data[dat->nData].string) == 2) {
			fprintf(stderr, "input data format has to be: [string] [string]\n");
			fprintf(stderr, "format error in line %d\n", n);
			exit(1);
		}
#else
        assert(fscanf(inFile, "%s %s\n", dat->data[dat->nData].label, dat->data[dat->nData].string) == 2);
#endif
		++ dat->nData;

		if (dat->nData == allocated) {
			allocated += 64;
			dat->data = safe_realloc(dat->data, allocated * sizeof(String));
		}
	}

	assert(dat->nData > 1);

	fclose(inFile);
	return 0;
}


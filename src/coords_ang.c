/*==============================================================================
coords_ang.c: object representation and related functions
Copyright (C) 2008-2015 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
==============================================================================*/

#include "coords_ang.h"

#define OVECCOUNT 18

/*____________________________________________________________________________*/
/* print header for point */
void print_header_object(FILE *outfile){
    fprintf(outfile, "     index      order  clusterId       CD       RD     phi1     phi2    theta\n");
}

/*____________________________________________________________________________*/
/* print data point */
void print_object(FILE *outfile, Dat *dat, int index, int order, int cluster_id, float coreDist, float reachDist){
    fprintf(outfile, "%10d %10d %10d %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
            index, order, cluster_id,
            coreDist, reachDist, 
            dat->data[index].phi1, dat->data[index].phi2, dat->data[index].theta);  
}

/*____________________________________________________________________________*/
/* calculate distance */
float calc_dist(Dat *dat, int i, int j, Arg *arg){

    float diff_phi1, diff_phi2, diff_theta;
    float dist;

    diff_phi1 = dat->data[i].phi1 - dat->data[j].phi1;
    diff_phi2 = dat->data[i].phi2 - dat->data[j].phi2;
    diff_theta = dat->data[i].theta - dat->data[j].theta;

    if (abs(diff_theta) > 180)
        diff_theta = 360 - abs(diff_theta);

    dist = sqrt(diff_phi1 * diff_phi1 + diff_phi2 * diff_phi2 + diff_theta * diff_theta);

    return dist;
}

/*____________________________________________________________________________*/
/* read data points */
/* Data points are expected to be an array of angles in 3D with
 * coordinates phi, psi and theta. */
int get_data(char *inFileName, Dat *dat){

	FILE *inFile = 0; /* input file */
	char line[80]; /* line */
	unsigned int allocated = 64; /* allocation counter */
    pcre *re; /* regular expression */
    const char *error; /* error message string */
	int erroffset; /* error offset */
    int ovector[OVECCOUNT]; /* match vector */
    int rc; /* match return value */
    char **substring_list; /* substring list */

	/* initialise/allocate memory for set of (64) frag entries */
    dat->nData = 0;
	dat->data = safe_malloc(allocated * sizeof(Ang));

	/* read data */
	inFile = safe_open(inFileName, "r");

	/* compile regexp */
    re = pcre_compile("^.{7}\t.{4}\t.{6}\t(........)\t(........)\t(........)", 0, &error, &erroffset, NULL);

    /* count the number of models */
    while(fgets(line, 80, inFile) != NULL )
    {
        rc = pcre_exec(re, NULL, line, strlen(line), 0, 0, ovector, OVECCOUNT);
        if (rc == 4){
            pcre_get_substring_list(line, ovector, rc, (const char ***)&substring_list);
            dat->data[dat->nData].phi1   = atof(substring_list[1]);
            dat->data[dat->nData].phi2   = atof(substring_list[2]);
            dat->data[dat->nData].theta  = atof(substring_list[3]);
            ++dat->nData;

            /* free substring_list memory */
            pcre_free_substring_list((const char **)substring_list);
        }

		/* allocate more memory if needed */
		if (dat->nData == allocated) {
			allocated += 64;
			dat->data = safe_realloc(dat->data, allocated * sizeof(Ang));
		}
    }

	/* free regexp */
    pcre_free(re);

	assert(dat->nData > 1);

    /* close file handle */
	fclose(inFile);
	return 0;
}

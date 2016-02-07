/*==============================================================================
$Id: objects.c 174 2008-07-17 15:20:22Z apandini $
dist_angle.c: distance calculation based on angles
Copyright (c) 2008 Jens Kleinjung and Alessandro Pandini

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
==============================================================================*/

#include "dist_angle.h"

#define OVECCOUNT 18

/*____________________________________________________________________________*/
/* calculate distance */
__inline__ float calc_dist(Dat *dat, int i, int j) {

    float diff_phi1, diff_phi2, diff_theta;

    diff_phi1 = dat->data[i].phi1 - dat->data[j].phi1;
    diff_phi2 = dat->data[i].phi2 - dat->data[j].phi2;
    diff_theta = dat->data[i].theta - dat->data[j].theta;

    if (abs(diff_theta) > 180)
        diff_theta = 360 - abs(diff_theta);

    return sqrt(diff_phi1 * diff_phi1 + diff_phi2 * diff_phi2 + diff_theta * diff_theta);
}

/*____________________________________________________________________________*/
/* read data points */
/* The expected input file is the result of anathons analysis */
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

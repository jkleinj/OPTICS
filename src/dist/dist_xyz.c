/*==============================================================================
$Id: objects.c 169 2008-07-16 17:13:01Z apandini $
dist_xyz.c: distance calculation based on xyz coordinates
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

#include "dist_xyz.h"
/*____________________________________________________________________________*/
/* calculate distance */
__inline__ float calc_dist(Dat *dat, int i, int j) {

    return vec_rmsd(&(dat->data[i]), &(dat->data[j]));
}

/*____________________________________________________________________________*/
/* read data points */
/* The order points are expected to be an array of vectors in 3D with
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
		assert(fscanf(inFile, "%f%f%f\n",
			&(dat->data[dat->nData].x), &(dat->data[dat->nData].y), &(dat->data[dat->nData].z)) == 3);
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


/*==============================================================================
$Id: objects.h 169 2008-07-16 17:13:01Z apandini $
dist_xyz.h: distance calculation based on xyz coordinates
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

#if !defined DIST_XYZ_H
#define DIST_XYZ_H

#include "safe.h"
#include "vector.h"

/* input data */
typedef struct {
	Vec *data; /* input data to order */
	int nData; /* number of data points */
} Dat;

/*____________________________________________________________________________*/
/* prototypes */
int get_data(char *inFileName, Dat *dat);
__inline__ float calc_dist(Dat *dat, int i, int j);

#endif

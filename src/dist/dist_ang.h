/*==============================================================================
$Id: objects.h 173 2008-07-17 12:38:44Z apandini $
dist_angle.h: distance calculation based on angles
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

#if !defined DIST_ANGLE_H
#define DIST_ANGLE_H

#include <math.h>
#include <pcre.h>
#include <string.h>

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
__inline__ float calc_dist(Dat *dat, int i, int j);

#endif

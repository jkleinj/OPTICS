/*==============================================================================
optics.h : ordering/clustering via the OPTICS method
Copyright (C) 2008-2015 Jens Kleinjung and Alessandro Pandini
References :
	- Pandini et al., BMC Bioinformatics 11:97, 2010.
	- Ankerst et al., Proc.ACM SIGMOD'99 Int. Conf. on Management
		of Data, Philadelphia PA, 1999.
	- Daszykowski et al., J. Chem. Inf. Comput. Sci. 42:500-507, 2002.

Read the COPYING file for license information.
==============================================================================*/

#ifndef OPTICS_H
#define OPTICS_H

#include <assert.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "config.h"
#include "safe.h"
#include "sort.h"

/* optional compilation of program */
/* set DATAANG instead using 'configure --enable-data-ang */
#ifdef DATAANG
	#include "coords_ang.h"
#endif
/* set DATASTR instead using 'configure --enable-data-str */
#ifdef DATASTR 
	#include "coords_str.h"
#endif
/* set DATAVEC instead using 'configure --enable-data-vec */
#ifdef DATAVEC
	#include "coords_vec.h"
#endif
/* set DATAXYZ instead using 'configure --enable-data-xyz */
#ifdef DATAXYZ
	#include "coords_xyz.h"
#endif
/* set DATADIST instead using 'configure --enable-data-dist */
#ifdef DATADIST
	#include "coords_dist.h"
	#include "matrix.h"
#endif

/*____________________________________________________________________________*/
/* structures */

/* ordering parameters */
typedef struct {
	float eps; /* order: neighbourhood radius */
	int minPts; /* order: minimum number of objects in neighbourhood */
} Par;

/* epsilon neighbour */
typedef struct {
	int index; /* point reference = point index */
	float dist; /* distance between point i and neighbour j */
} Epsn;

/* point data */
typedef struct {
	int processed; /* flag indicating whether point has been processed */
    int index; /* index in the input data array */
    int order; /* index in the OPTICS ordered list */
	float coreDist; /* core-distance(i) = distance(i,j=MinPts), otherwsise undefined */
	float reachDist; /* reachability-distance(i,j) = max(core-distance(i), distance(i,j)) */
	Epsn *epsNeigh; /* epsilon-neighbourhood: list of neighbour points */
	int nEpsNeigh; /* epsilon-neighbourhood: number of neighbour points */
} Pt;

/* optics data */
typedef struct {
	Pt *pt; /* points */
	int nPt; /* number of points */
} OpticsDat;

/* cluster */
typedef struct {
	int start; /* index of starting point */
	int end; /* index of ending point */
	int parent; /* parent cluster */
    int size; /* cluster size */
	float minCD; /* min CD in the cluster */
	int minCDid; /* id of point with minCD */
} Cluster;

/* cluster list*/
typedef struct {
	Cluster *cluster; /* clusters */
	int nCluster; /* number of clusters */
} ClusterList;

#endif

/*=============================================================================
optics.c : OPTICS single-linkage point ordering
Rferences :
	- Pandini et al., BMC Bioinformatics 11:97, 2010.
	- Ankerst et al., Proc.ACM SIGMOD'99 Int. Conf. on Management
		of Data, Philadelphia PA, 1999.
	- Daszykowski et al., J. Chem. Inf. Comput. Sci. 42:500-507, 2002.

Copyright (C) 2008-2015 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
=============================================================================*/

#include "optics.h"

/*_____________________________________________________________________________*/
/** global parameters */
int silent = 0;

/*____________________________________________________________________________*/
/** min floats */
__inline__ static float min_float(float a, float b)
{
	return ((a < b ) ? a : b);
}

/*____________________________________________________________________________*/
/** max floats */
__inline__ static float max_float(float a, float b)
{
	return ((a > b ) ? a : b);
}

/*____________________________________________________________________________*/
/** compare distances */
__inline__ static int compare_distances(Epsn *a, Epsn *b)
{
	float arg_a = (float) a->dist;
	float arg_b = (float) b->dist;

	return ((arg_a < arg_b ) ? -1 : (arg_a > arg_b) ? 1 : 0);
}

/*____________________________________________________________________________*/
/** compare RD by dereferencing */
__inline__ static int compare_RD_deref(Pt **a, Pt **b)
{
	float arg_a = (float) (**a).reachDist;
	float arg_b = (float) (**b).reachDist;

    /* WARNING: decreasing */
	return ((arg_a < arg_b ) ? 1 : (arg_a > arg_b) ? -1 : 0);
}

/*____________________________________________________________________________*/
/* initialise parameters */
static int initialise(OpticsDat *opticsdat, Par *par)
{
	unsigned int i;

	for (i = 0; i < opticsdat->nPt; ++ i) {
		/* optics data point */
		opticsdat->pt[i].processed = 0;
		opticsdat->pt[i].index = -1;
        /* core distance is set to eps + 1
           to identify noise point easily */
		opticsdat->pt[i].coreDist = par->eps + 1.0;
        /* reachability set to eps */
        /* see Daszykowski et al. for details */
		opticsdat->pt[i].reachDist = par->eps;
		opticsdat->pt[i].nEpsNeigh = 0;
	}

	return 0;
}

/*____________________________________________________________________________*/
/* unset processed flag */
static int reset_processed_flag(Pt **pointer_to_data, int npoints)
{
	unsigned int i;

	for (i = 0; i < npoints; ++ i) {
		/* set process flag to 0 */
		(*pointer_to_data[i]).processed = 0;
	}

	return 0;
}

/*____________________________________________________________________________*/
/* epsilon-neighbourhood */
/* This defines the neighbourhood of point i: the number of neighbouring
 *	points within the boundary epsilon, a list of their 'index' and their 
 *	distance 'dist'. */
static void epsilon_neighbourhood(Dat *dat, Par *par, OpticsDat *opticsdat, int i, Arg *arg)
{
	unsigned int j;
	unsigned int allocated = 64;
	float dist = 0;

	/* allocate space to record neighbours */
	opticsdat->pt[i].epsNeigh = safe_malloc(allocated * sizeof(Epsn));

	for (j = 0, opticsdat->pt[i].nEpsNeigh = 0; j < dat->nData; ++ j) {
		if ((j == i) || (opticsdat->pt[j].processed > 0)) continue;

		/* record all points in epsilon neighbourhood */
		if ((dist = calc_dist(dat, i, j, arg)) <= par->eps) {
			/* assign neighbour point */
			opticsdat->pt[i].epsNeigh[opticsdat->pt[i].nEpsNeigh].index = j;
			/* assign neighbour distance */
			opticsdat->pt[i].epsNeigh[opticsdat->pt[i].nEpsNeigh].dist = dist;
			/* count neighbours */
			++ opticsdat->pt[i].nEpsNeigh;
		}

		/* allocate more space to neighbour array if needed */
		if (opticsdat->pt[i].nEpsNeigh == allocated) {
			allocated += 64;
			opticsdat->pt[i].epsNeigh = safe_realloc(opticsdat->pt[i].epsNeigh, allocated * sizeof(Epsn));
		}
	}

    if (! silent) 
        fprintf(stderr, "point %d has %d neighbours within %5.2e\n",
            i, opticsdat->pt[i].nEpsNeigh, par->eps);
}

/*____________________________________________________________________________*/
/* define core-distance of point i and reachability distance to all 
	epsilon neighbours j; update reachability-distance of point j 
	if the current one (to point i) is shorter than previously set/recorded */
static int cd_rd(Par *par, OpticsDat *opticsdat, int i)
{
	unsigned int j, k = 0;
	int (*fcmp)() = &compare_distances;
	float rd;
    /* minimum reachability distance(ij) set to eps */
    /* see Daszykowski et al. for details */
	float rd_min = par->eps; /* lowest RD */
	int i_rd_min = -1; /* index of neighbour j with lowest RD */
	int next = -1; /* next point: the one among j with minimum RD to point i */

	/* Definition 5 of Ref. 1.: 
		core-distance(i) = distance(i,j=MinPts) */
	/* if point i has enough neighbours */
	if (opticsdat->pt[i].nEpsNeigh >= par->minPts) {
		/* sort neighbours by neighbour distance */
		MergeSort((void *)opticsdat->pt[i].epsNeigh, opticsdat->pt[i].nEpsNeigh, sizeof(Epsn), fcmp);
		/* record neighbour distance of minPts point as core distance */
		opticsdat->pt[i].coreDist = opticsdat->pt[i].epsNeigh[par->minPts - 1].dist;
		/* Definition 6 of Ref. 2.:
			reachability-distance(i,j) = max(core-distance(i), distance(ij)) */
		for (j = 0; j < opticsdat->pt[i].nEpsNeigh; ++ j) {
            k = opticsdat->pt[i].epsNeigh[j].index; /* index of neighbourogh j in ord */
			rd = max_float(opticsdat->pt[i].coreDist, opticsdat->pt[i].epsNeigh[j].dist);
			opticsdat->pt[k].reachDist = min_float(rd, opticsdat->pt[k].reachDist); 
			/* define index of closest neighbour */
			if ((rd_min = min_float(rd_min, opticsdat->pt[k].reachDist)) == opticsdat->pt[k].reachDist)
				i_rd_min = k;
		}
		next = i_rd_min;
	} else {
		/* else not enough neighbours: set core distance and
			all reachability distances to 'undefined' */
        /* core distance is set to eps + 1
           to identify noise point easily */
		opticsdat->pt[i].coreDist = par->eps + 1.0;
		rd = FLT_MAX;
		for (j = 0; j < opticsdat->pt[i].nEpsNeigh; ++ j) {
            k = opticsdat->pt[i].epsNeigh[j].index; /* index of neighbour j in ord */
			opticsdat->pt[k].reachDist = min_float(rd, opticsdat->pt[k].reachDist); 
		}
		next  = -1;
	}
	
	return next;
}

/*____________________________________________________________________________*/
/* find next unprocessed point */
static int find_next(OpticsDat *opticsdat)
{
	unsigned int i = 0;

	while (i < opticsdat->nPt){
		if (opticsdat->pt[i].processed == 0)
			return i;
        ++ i;
    }

	return -1; /* when all points are already processed */
}

/*____________________________________________________________________________*/
/* order points */
/* Algorithm (from Ref. 2.):
	1. Select randomly a first object. It is called the current object.
		Its RD (= reachability distance) is undefined. Mark this object
		as processed, and plot it in the reachability plot in the first position.
	2. Calculate the RD of all objects with the current object.
	3. Select the object with the smallest RD to the current object and plot its
		RD in the reachability plot in the next position. Mark it processed.
		This object is now considered as the current object.
	4. Calculate the RD of all the remaining not processed objects with respect
		to the current object, and if the current RD for any object is smaller
		than the previous RD for that object, replace it with the current RD.
	5. Go to 3 and continue until all objects are processed. */

static int order(Dat *dat, Par *par, OpticsDat *opticsdat, int i, int *ptr_processed, Arg *arg)
{
	int next = -1;
    int perc;

    perc = (*ptr_processed  * 100 / dat->nData);

	/* the point is flagged up as processed */
	++ opticsdat->pt[i].processed;
    (*ptr_processed) ++;
    if ( ((*ptr_processed * 100 / dat->nData) != perc) & !silent){
        fprintf(stderr, "Processing... %3d%% completed\n", (*ptr_processed  * 100 / dat->nData));
    }

	/* and its index recorded */
	opticsdat->pt[i].index = i;

	/* Record unprocessed points in epsilon neighbourhood of point i. */
	epsilon_neighbourhood(dat, par, opticsdat, i, arg);
	/* Compute the CD (= core distance) of point i
		and update the RD of points in its epsilon neighbourhood.
		The returned point index 'next' is the closest neighbour of i. */
	next = cd_rd(par, opticsdat, i);

    /* empty the neighbour list to save memory */
	free(opticsdat->pt[i].epsNeigh);

	/* if CD and RD undefined, find next point to process */
	if (next < 0)
		next = find_next(opticsdat);

	return next;
}

/*____________________________________________________________________________*/
/* order data by RD starting from a OPTICS ordered vector of pointers 
        where first and last are indexes of the first and last element to 
        sort in the ordered_data array */
void order_by_RD(Pt **RD_ordered_data, Pt **ordered_data, int first, int last)
{
    int i, npoints;
	int (*fcmp)() = &compare_RD_deref;

    /* number of points to order */
    npoints = last - first + 1;

	/* initialise vector with pointers */
    for(i = 0; i < npoints; ++i){
        RD_ordered_data[i] = ordered_data[i + first];
    }

    /* sort ordered vector by RD */
    MergeSort((void *)RD_ordered_data, npoints, sizeof(Pt *), fcmp);
}

/*____________________________________________________________________________*/
/* extract root clusters */
void extract_root_clusters(ClusterList *cluster_list, Pt **RD_ordered_data, Pt **ordered_data, int npoints, Par *par, int *ptr_mem_allocated, int *ptr_pseudoClusterFlag)
{
    int i, j;
    Pt *clStart = NULL; /* start cluster point */
    float minCD; /*  min CD in the cluster */
    int minCDid; /* id of point with minCD */
    int totalProcessed = 0; /* processed counter */
    int cluster_size; /* size of root cluster */

    /* memory allocation for cluster list */
    *ptr_mem_allocated = 10; /* allocation counter */
    cluster_list->cluster = safe_malloc((*ptr_mem_allocated) * sizeof(Cluster));

    /* start with highest RD point */
    /* WARNING: sorting should be stable
     *          therefore the first highest point (RD = eps)
     *          would be the first ordered point in the Reachability plot */
    clStart = RD_ordered_data[0];

    for (i = 1; i < npoints; ++i) {
        if ((*RD_ordered_data[i]).processed == 1) {
            continue;
		}
        /* check MinPts condition */
        cluster_size = (*RD_ordered_data[i]).order - (*clStart).order;
        if (cluster_size > par->minPts ) {
            /* set minCD to highest value */
            minCD = par->eps + 1;
            minCDid = -1;

            for (j = (*clStart).order; j < (*RD_ordered_data[i]).order; ++j) {
                /* mark cluster point as processed */
                (*ordered_data[j]).processed = 1;
                totalProcessed ++;
                /* get minCD and minCDid */
                if ((*ordered_data[j]).coreDist < minCD) {
                    minCD = (*ordered_data[j]).coreDist;
                    minCDid = (*ordered_data[j]).order;
                }
            }

            /* add cluster to list */
            cluster_list->cluster[cluster_list->nCluster].start = (*clStart).order;
            /* last cluster point is set to i - 1 */
            cluster_list->cluster[cluster_list->nCluster].end = (*RD_ordered_data[i]).order - 1;
            cluster_list->cluster[cluster_list->nCluster].parent = -1;
            cluster_list->cluster[cluster_list->nCluster].size = cluster_size; 
            cluster_list->cluster[cluster_list->nCluster].minCD = minCD;
            cluster_list->cluster[cluster_list->nCluster].minCDid = minCDid;

            /* if cluster has no valid representative raise a flag */
            if (minCDid == -1) {
                (*ptr_pseudoClusterFlag) = 1;
			}

            /* increase cluster counts */
            cluster_list->nCluster ++;

            /* allocate more memory if needed */
            if (cluster_list->nCluster == (*ptr_mem_allocated)) {
                (*ptr_mem_allocated) += 10;
                cluster_list->cluster = safe_realloc(cluster_list->cluster, (*ptr_mem_allocated) * sizeof(Cluster));
            }

            /* new cluster can start */
            clStart = RD_ordered_data[i];
		/* cluster_size <= par->minPts */
        } else {
            /* mark cluster point as processed */
            (*RD_ordered_data[i]).processed = 1;
            totalProcessed ++;
            /* move on if points are levered */ 
            if ((*RD_ordered_data[i]).reachDist == (*RD_ordered_data[i - 1]).reachDist) {
                clStart = RD_ordered_data[i];
			}
        }
    }

    /* print number of processed point */
    if (! silent) 
        fprintf(stderr, "\n%d/%d points processed in root clustering\n", totalProcessed + 1, npoints);

    /* safe free if no cluster has been identified */
    if (cluster_list->nCluster == 0) {
        free(cluster_list->cluster);
    } else {
        /* add last point to right end side cluster if it is the only missing point */
        if ((npoints - cluster_list->cluster[cluster_list->nCluster - 1].end) == 2) {
            cluster_list->cluster[cluster_list->nCluster - 1].end ++;
		}
	}
}

/*____________________________________________________________________________*/
/* extract sub clusters */
int extract_sub_clusters(Cluster *cluster, ClusterList * cluster_list, Pt **ordered_data, int parent_id, Par *par, int *ptr_mem_allocated, int *ptr_pseudoClusterFlag)
{
    int i, j; /* counters */
    int start, end, npoints; /* start index, end index and number of points */
    Pt **cluster_RD_ordered_data; /* list of RD ordered data (cluster) */
    Pt *clStart = NULL; /* start cluster point */
    float minCD; /*  min CD in the cluster */
    int minCDid; /* id of point with minCD */
    float diff;
    int cluster_size; /* size of subcluster */
    int firstNewCluster = -1; /* index of first new cluster */

    start = cluster->start;
    end = cluster->end;
    npoints = end - start + 1;

    cluster_RD_ordered_data = safe_malloc(npoints * sizeof(Pt *));

    order_by_RD(cluster_RD_ordered_data, ordered_data, start, end);

    /* reset processed flag */
    reset_processed_flag(cluster_RD_ordered_data, npoints);

    /* start with highest RD point */
    /* WARNING: sorting should be stable
     *          therefore the first highest point (RD = eps)
     *          would be the first ordered point in the Reachability plot */
    clStart = cluster_RD_ordered_data[0];

	for (i = 1; i < npoints; ++ i) {
        if ((*cluster_RD_ordered_data[i]).processed == 1) {
            continue;
		}
        /* calculate RD diff between points */
        diff = (*cluster_RD_ordered_data[i]).reachDist - (*ordered_data[(*cluster_RD_ordered_data[i]).order - 1]).reachDist;
        /* cluster is present if increasing inflexion point */
        if (diff > 0) {
            /* check for MinPts conditions*/
            cluster_size = (*cluster_RD_ordered_data[i]).order - (*clStart).order;
            if (    /* cluster should include at least MinPts objects */
                    (cluster_size > par->minPts) &
                    /* cluster should be at least MinPts smaller than parent cluster */
                    ((cluster_list->cluster[parent_id].size - cluster_size ) >= par->minPts )) { 

                /*update first new cluster index */
                if (firstNewCluster == -1) {
                    firstNewCluster = cluster_list->nCluster;
				}

                /* set minCD to highest value */
                minCD = par->eps + 1; 
                minCDid = -1; 

                for (j = (*clStart).order; j < (*cluster_RD_ordered_data[i]).order; ++j) {
                    /* mark cluster point as processed */
                    (*ordered_data[j]).processed = 1;
                    /* get minCD and minCDid */
                    if ((*ordered_data[j]).coreDist < minCD) {
                        minCD = (*ordered_data[j]).coreDist;
                        minCDid = (*ordered_data[j]).order;
                    }
                }

                /* add cluster to list */
                cluster_list->cluster[cluster_list->nCluster].start = (*clStart).order;
                /* last cluster point is set to i - 1 */
                cluster_list->cluster[cluster_list->nCluster].end = (*cluster_RD_ordered_data[i]).order - 1;
                cluster_list->cluster[cluster_list->nCluster].size = cluster_size; 
                cluster_list->cluster[cluster_list->nCluster].parent = parent_id;
                cluster_list->cluster[cluster_list->nCluster].minCD = minCD;
                cluster_list->cluster[cluster_list->nCluster].minCDid = minCDid;

                /* if cluster has no valid representative raise a flag */
                if (minCDid == -1) {
                    (*ptr_pseudoClusterFlag) = 1;
				}

                /* increase cluster counts */
                cluster_list->nCluster ++;

                /* allocate more memory if needed */
                if (cluster_list->nCluster == (*ptr_mem_allocated)) {
                    (*ptr_mem_allocated) += 10;
                    cluster_list->cluster = safe_realloc(cluster_list->cluster, (*ptr_mem_allocated) * sizeof(Cluster));
                }
                /* next starting point is the actual left hand side */
                clStart = cluster_RD_ordered_data[i];
            }
        }
    }

    free(cluster_RD_ordered_data);

    return(firstNewCluster);
}

/*____________________________________________________________________________*/
int main(int argc, char *argv[])
{
	unsigned int i = 0; /* index */

    Arg arg; /** data structure for command line arguments */
	Dat dat; /* input data */
	Par par; /* ordering parameters */
	OpticsDat opticsdat; /* optics data points */
    int processed = 0; /*  of processed */
    Pt **ordered_data; /* list of ordered data */
    Pt **RD_ordered_data; /* list of RD ordered data */
	int next = 0; /* index of first point */
    ClusterList cluster_list; /* cluster list */
    int mem_allocated; /* memory counter for cluster list */
    int firstNewCluster = 0; /* index of first new cluster */
    int pseudoClusterFlag = 0; /* pseudo cluster flag */
    char pseudoClusterOutFileName[128] = "pseudo."; /* pseudo cluster file name */
    FILE *pseudoClusterOutFile = 0; /* pseudo cluster file handle */
	char outName[256] = "";
    
    /*____________________________________________________________________________*/
    /** parse command line arguments */
    parse_args(argc, &(argv[0]), &arg);

	/*____________________________________________________________________________*/
	/* parametrise */
	par.eps = arg.eps; /* order: neighbourhood radius */
	par.minPts = arg.minPts; /* order: minimum number of objects in neighbourhood */

	/*____________________________________________________________________________*/
	/** read input data */
	get_data(arg.dataInFileName, &dat);
    if (! silent) 
        fprintf(stderr, "%d points\n", dat.nData);
	opticsdat.nPt = dat.nData;

	/*____________________________________________________________________________*/
	/** perform ordering */
	/* allocate memory to points */
	opticsdat.pt = safe_malloc(opticsdat.nPt * sizeof(Pt));	
    ordered_data = safe_malloc(opticsdat.nPt * sizeof(Pt *));
    RD_ordered_data = safe_malloc(opticsdat.nPt * sizeof(Pt *));

	/* initialise optics data points */
	initialise(&opticsdat, &par);

	/* the RD of the first point is undefined and set to eps */
    /* see Daszykowski et al. for details */
	opticsdat.pt[next].reachDist = par.eps;

	/* order the points */
    while (next != -1) {
        /* the point is added to the order list */
        ordered_data[i] = &(opticsdat.pt[next]);
        /* and its order attribute is updated */
        opticsdat.pt[next].order = i;
        /* Compute CD and RD then return point index 'next' of closest neighbour */
        next = order(&dat, &par, &opticsdat, next, &processed, &arg);
        i ++; 
    }

	/* order the points by RD */
    order_by_RD(RD_ordered_data, ordered_data, 0, (opticsdat.nPt - 1));

    /* reset processed flag */
    reset_processed_flag(RD_ordered_data, opticsdat.nPt);

    /* initialise cluster list */
    cluster_list.nCluster = 0;

    /* extract root clusters */
    extract_root_clusters(&cluster_list, RD_ordered_data, ordered_data, opticsdat.nPt, &par, &mem_allocated, &pseudoClusterFlag);

    /* extract subclusters */
	if (cluster_list.nCluster > 0) {
		while(firstNewCluster != -1)
			for (i = firstNewCluster; i < cluster_list.nCluster; ++i)
				firstNewCluster = extract_sub_clusters(&(cluster_list.cluster[i]), &cluster_list, ordered_data, i, &par, &mem_allocated, &pseudoClusterFlag);
    } else {
        fprintf(stderr, "No clusters found! Try larger 'epsilon' or smaller 'minpts'\n");
        exit(1);
    }   

    /* open output files */
	/* concatenate path and file name */
	sprintf(outName, "%s/%s", arg.outPathName, arg.dataOutFileName);
    arg.dataOutFile = safe_open(outName, "w");
	sprintf(outName, "%s/%s", arg.outPathName, arg.clusterOutFileName);
    arg.clusterOutFile = safe_open(outName, "w");
	sprintf(outName, "%s/%s", arg.outPathName, arg.centerOutFileName);
    arg.centerOutFile = safe_open(outName, "w");
	sprintf(outName, "%s/%s", arg.outPathName, arg.uniqueOutFileName);
    arg.uniqueOutFile = safe_open(outName, "w");
    /* open pseudo cluster file if needed */
    if (pseudoClusterFlag != 0){
		sprintf(outName, "%s/%s%s",
			arg.outPathName, pseudoClusterOutFileName, arg.clusterOutFileName);
        pseudoClusterOutFile = safe_open(outName, "w");
    }

    /* print ordered points */
    if (! silent)
        fprintf(stderr, "Sorting Completed\n\n");
    fprintf(arg.dataOutFile, "    dataId       CD       RD\n");
	for (i = 0; i < opticsdat.nPt; ++ i){
        fprintf(arg.dataOutFile, "%10d %8.3f %8.3f\n", (*ordered_data[i]).index, (*ordered_data[i]).coreDist, (*ordered_data[i]).reachDist);
    }

    /* print clustering results */
    fprintf(arg.clusterOutFile, "  id   parent    start      end     size    minCD  minCDid\n");
    print_header_object(arg.centerOutFile);
    print_header_object(arg.uniqueOutFile);
    /* print header to pseudo cluster file if needed */
    if (pseudoClusterFlag != 0)
        fprintf(pseudoClusterOutFile, "  id   parent    start      end     size    minCD  minCDid\n");
    for (i = 0; i < cluster_list.nCluster;  ++i){
        /* check if cluster has a valid center.
         * i.e. at least one core point is in the cluster */
        if (cluster_list.cluster[i].minCDid != -1) {
            fprintf(arg.clusterOutFile, "%4d %8d %8d %8d %8d %8.3f %8d\n", i, cluster_list.cluster[i].parent, cluster_list.cluster[i].start, cluster_list.cluster[i].end, cluster_list.cluster[i].size, cluster_list.cluster[i].minCD, cluster_list.cluster[i].minCDid);
            print_object(arg.centerOutFile, &dat, (*ordered_data[cluster_list.cluster[i].minCDid]).index, (*ordered_data[cluster_list.cluster[i].minCDid]).order,  i, (*ordered_data[cluster_list.cluster[i].minCDid]).coreDist, (*ordered_data[cluster_list.cluster[i].minCDid]).reachDist);
            if (cluster_list.cluster[i].parent == -1)
                print_object(arg.uniqueOutFile, &dat, (*ordered_data[cluster_list.cluster[i].minCDid]).index, (*ordered_data[cluster_list.cluster[i].minCDid]).order,  i, (*ordered_data[cluster_list.cluster[i].minCDid]).coreDist, (*ordered_data[cluster_list.cluster[i].minCDid]).reachDist);
            else
                if (cluster_list.cluster[i].minCDid != cluster_list.cluster[cluster_list.cluster[i].parent].minCDid)
                    print_object(arg.uniqueOutFile, &dat, (*ordered_data[cluster_list.cluster[i].minCDid]).index, (*ordered_data[cluster_list.cluster[i].minCDid]).order,  i, (*ordered_data[cluster_list.cluster[i].minCDid]).coreDist, (*ordered_data[cluster_list.cluster[i].minCDid]).reachDist);
        }else{
            /* if the cluster has no valid center
             * print the info in a pseudo cluster file */
            fprintf(pseudoClusterOutFile, "%4d %8d %8d %8d %8d %8.3f %8d\n", i, cluster_list.cluster[i].parent, cluster_list.cluster[i].start, cluster_list.cluster[i].end, cluster_list.cluster[i].size, cluster_list.cluster[i].minCD, cluster_list.cluster[i].minCDid);
        }
    }

    /* close output files */
    fclose(arg.dataOutFile);
    fclose(arg.clusterOutFile);
    fclose(arg.centerOutFile);
    fclose(arg.uniqueOutFile);
    /* close pseudo cluster file if needed */
    if (pseudoClusterFlag != 0)
        fclose(pseudoClusterOutFile);

	/*____________________________________________________________________________*/
	/** free memory */
#ifdef DATADIST
	free_mat2D_float(dat.dist, dat.nData);
#else
	free(dat.data);
#endif
	free(opticsdat.pt);
	free(ordered_data);
	free(RD_ordered_data);
    if (cluster_list.nCluster > 0)
        free(cluster_list.cluster);

	/*____________________________________________________________________________*/
	/** terminate */
    if (! silent)
        fprintf(stderr, "Clean termination\n");
    return 0;
}

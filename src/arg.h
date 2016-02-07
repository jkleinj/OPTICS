/*==============================================================================
arg.h : parse command line arguments
Copyright (C) 2008-2015 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARG_H
#define ARG_H

#include <float.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

/*____________________________________________________________________________*/
/* structures */

/* non-parametric arguments */
typedef struct  
{
    FILE *dataInFile;
    char *dataInFileName;
    FILE *dataOutFile;
    char *dataOutFileName;
    FILE *clusterOutFile;
    char *clusterOutFileName;
    FILE *centerOutFile;
    char *centerOutFileName;
    FILE *uniqueOutFile;
    char *uniqueOutFileName;
    float eps;
    int minPts;
    int w;
	char *outPathName;
} Arg;

/*____________________________________________________________________________*/
/* prototypes */
int parse_args(int argc, char **argv, Arg *arg);

#endif

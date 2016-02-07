/*===============================================================================
sort.h : sort alogorithms implementation
(C) 2008-2013 by Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
================================================================================*/

#ifndef SORT_H 
#define SORT_H 

#include <stdio.h>
#include <string.h>

#include "safe.h"

/*____________________________________________________________________________*/
/* prototypes */
void MergeSort(void *array, size_t size, size_t esize, int (*compare) (const void *key1, const void *key2));

#endif


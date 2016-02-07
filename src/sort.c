/*=============================================================================
sort: sort algorithms implementation
(C) 2008-2013 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
==============================================================================*/

#include "sort.h"

/*_____________________________________________________________________________*/
/** merge function */
void merge(unsigned char *input, int left, int right, int length, int midpoint, size_t esize, unsigned char *scratch, int (*compare) (const void *key1, const void *key2)){

    int i = 0;
    int l, r; /* counters for left and right subarrays */

    l = left; /* left subarray starts at array start */ 
    r = left + midpoint; /* right subarray starts at array midpoint */

    /* operation are performed by copy esize chunks of memory
     * indexing is dealt consistently by bytes */

    /* merge the subarrays together */ 
    /* because (l+r == length) no risk of going out of array boundaries */
    for(i = 0; i < length; i++){
        /* check if left array has elements */
        if (l < left + midpoint){
            /* check if right array has elements */
            if (r < right){
                /* compare element in position l and r */
                if (compare(&input[l * esize], &input[r * esize]) <= 0){
                    memcpy(&scratch[i * esize],  &input[l * esize],  esize);
                    l++;
                }else{
                    memcpy(&scratch[i * esize],  &input[r * esize],  esize);
                    r++;
                }
            }else{
                /* otherwise add element from left array directly */
                memcpy(&scratch[i * esize],  &input[l * esize],  esize);
                l++;
            }
        }else{
            /* otherwise add element from right array directly */
            memcpy(&scratch[i * esize],  &input[r * esize],  esize);
            r++;
        }
    }

    /* update input array with sorted subarray from scratch */
    for(i = left; i < right; i++)
        memcpy(&input[i * esize], &scratch[(i - left) * esize], esize);
}

/*_____________________________________________________________________________*/
/** mergesort helper 
 * this is the function that actually perform the sorting where
 * left is the index of the leftmost element of the array
 * right is one past the index of the rightmost element */
void mergesort_hlp(unsigned char *input, int left, int right, size_t esize, unsigned char *scratch, int (*compare) (const void *key1, const void *key2))
{
    int length = right - left;
    int midpoint = length/2;

    /* if the array includes only one element return to parent call */
    if(length == 1)
        return;

    /* if the array includes more than one element sorting is performed */

    /* recursively sort each subarray */
    mergesort_hlp(input, left, left + midpoint, esize, scratch, compare);
    mergesort_hlp(input, left + midpoint, right, esize, scratch, compare);

    /* merge the subarrays */
    merge(input, left, right, length, midpoint, esize, scratch, compare);
}

/*_____________________________________________________________________________*/
/** mergesort */
void MergeSort(void *array, size_t size, size_t esize, int (*compare) (const void *key1, const void *key2))
{
    /* allocate memory for scratch array */
    unsigned char *scratch = safe_malloc(size * esize * sizeof(unsigned char));

    /* perform merge sorting */
    mergesort_hlp(array, 0, size, esize, scratch, compare);

    /* free memory from scratch array */
    free(scratch);
}


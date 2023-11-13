/* functions which operate on BAND_VECT structures */

#include "global_header.h"

/* returns the dot product of one banded vector with another thru the given range of indices */
double bv_dot(BAND_VECT * v1,   /* the first vector */
              BAND_VECT * v2,   /* the second vector */
              int idot_begin,   /* the beginning index for the dot product */
              int idot_end      /* the ending index for the dot product */
    )
{
    int i;                      /* loop counter */
    int ibegin;                 /* the beginning of the loop */
    int iend;                   /* the end of the loop */
    int isum;                   /* the number of entries to sum */
    double value;               /* the return value */
    double *a1, *a2;            /* pointers to the array values */

    /* calculate the beginning and end of the loop */
    ibegin = idot_begin;
    iend = idot_end;
    if (ibegin < v1->begin)
        ibegin = v1->begin;
    if (ibegin < v2->begin)
        ibegin = v2->begin;
    if (iend > v1->end)
        iend = v1->end;
    if (iend > v2->end)
        iend = v2->end;
    a1 = (v1->value) + (ibegin - v1->begin);
    a2 = (v2->value) + (ibegin - v2->begin);
    isum = iend - ibegin;

    /* calculates the dot product */
    for (i = 0, value = 0.0; i < isum; i++)
        value += a1[i] * a2[i];

    /* returns the calculated value */
    return (value);
}

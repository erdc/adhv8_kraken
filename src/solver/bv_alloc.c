/* functions which operate on BAND_VECT structures */

#include "global_header.h"

/* allocates a band vector */
void bv_alloc(BAND_VECT * bv_pntr   /* pointer to the block matrix to be initialized */
    )
{
    int i;                      /* loop counter */
    int new_size;               /* the new size of the band vector */

    /* calculate the new size of the band vector */
    new_size = bv_pntr->end - bv_pntr->begin;

    /* if the new size is larger than the current size, then allocate additional space */
    if (new_size > bv_pntr->size) {
        bv_pntr->value = (double *) tl_realloc(sizeof(double), new_size, bv_pntr->size, bv_pntr->value);
        for (i = bv_pntr->size; i < new_size; i++)
            bv_pntr->value[i] = 0.0;
        bv_pntr->size = new_size;
    }
}

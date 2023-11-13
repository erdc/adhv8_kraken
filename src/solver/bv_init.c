/* functions which operate on BAND_VECT structures */

#include "global_header.h"

/* initializes a band vector */
void bv_init(BAND_VECT * bv_pntr    /* pointer to the band vector to be initialized */
    )
{
    int i;                      /* loop counter */

    /* initialize the entries in the block matrix */
    bv_pntr->begin = 0;
    bv_pntr->end = 0;
    for (i = 0; i < bv_pntr->size; i++)
        bv_pntr->value[i] = 0.0;
}

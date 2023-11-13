/* functions which operate on BAND_VECT structures */

#include "global_header.h"

/* loads an entry in a band vector */
void bv_load(BAND_VECT * bv_pntr,   /* pointer to the band vector the entry goes to */
             double entry,      /* the value to be entered */
             int index          /* the index of the entry */
    )
{
    int i, j;                   /* loop counter */
    int iend;                   /* the end of the loop */
    int relative_index;         /* the index within the band vector */
    int new_size;               /* the size needed for the band vector */

    /* checks if the entry is before the begining of the vector and 
       enlarges the vector if needed */
    if (index < bv_pntr->begin) {
        /* checks to see if more space is needed */
        new_size = bv_pntr->end - index + 1;
        if (new_size > bv_pntr->size) {
            new_size += BV_BLOCK;
            bv_pntr->value = (double *) tl_realloc(sizeof(double), new_size, bv_pntr->size, bv_pntr->value);
            for (i = bv_pntr->size; i < new_size; i++)
                bv_pntr->value[i] = 0.0;
            bv_pntr->size = new_size;
        }

        /* shifts the band vector for the new begining */
        for (i = bv_pntr->end - bv_pntr->begin, j = bv_pntr->end - index; i >= 0; i--, j--)
            bv_pntr->value[j] = bv_pntr->value[i];
        for (i = 0, iend = bv_pntr->begin - index; i < iend; i++)
            bv_pntr->value[i] = 0.0;
        bv_pntr->begin = index;
    }

    /* checks if the entry is after the end of the vector and enlarges the 
       vector if needed */
    if (index >= bv_pntr->end) {
        /* checks to see if more space is needed */
        new_size = index - bv_pntr->begin + 1;
        if (new_size > bv_pntr->size) {
            new_size += BV_BLOCK;
            bv_pntr->value = (double *) tl_realloc(sizeof(double), new_size, bv_pntr->size, bv_pntr->value);
            for (i = bv_pntr->size; i < new_size; i++)
                bv_pntr->value[i] = 0.0;
            bv_pntr->size = new_size;
        }

        /* sets the new end to the band vector */
        bv_pntr->end = index + 1;
    }

    /* sets the entry in the vector */
    relative_index = index - bv_pntr->begin;
    bv_pntr->value[relative_index] = entry;
}

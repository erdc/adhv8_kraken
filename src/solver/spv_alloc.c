/* allocates space for the sparse vector */

#include "global_header.h"

void spv_alloc(SPARSE_VECT * sv,    /* the sparse vector */
               int max_nsys_sq)
{
    int new_max_size;           /* the new max size */
    int new_max_val_size;       /* the new max value size */
    int old_max_val_size;       /* the old max value size */
    int i;                      /* loop counter */

    /* checks to see if space is needed */
    if (sv->size > sv->max_size) {
        /* sets the new max size */
        new_max_size = sv->max_size;
        while (sv->size > new_max_size)
            new_max_size += SPV_BLOCK;
        old_max_val_size = sv->max_size * max_nsys_sq;
        new_max_val_size = new_max_size * max_nsys_sq;

        /* allocates space */
        sv->value = (double *) tl_realloc(sizeof(double), new_max_val_size, old_max_val_size, sv->value);
        sv->index = (int *) tl_realloc(sizeof(int), new_max_size, sv->max_size, sv->index);

        /* initializes the space */
        for (i = old_max_val_size; i < new_max_val_size; i++)
            sv->value[i] = 0.0;
        for (i = sv->max_size; i < new_max_size; i++)
            sv->index[i] = UNSET_INT;

        /* resets the max size */
        sv->max_size = new_max_size;
    }
    return;
}

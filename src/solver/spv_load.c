/* loads the given entry in the sparse vector */

#include "global_header.h"

void spv_load(SPARSE_VECT * sv, /* the sparse vector */
              int ind,          /* index to be added to the vector */
              double *value,    /* the value to be added to the entry */
              int p,            /* the number of equations being solved */
              int p2,           /* the number of equations being solved squared */
              int max_nsys_sq)
{
    int i = 0;                  /* loop counter */
    int j = 0;                  /* the value positions */


    /* loops over the indices to find the entry */
    i = 0;
    while (i < sv->size && ind != sv->index[i] && sv->index[i] != UNSET_INT)
        i++;

    /* adds the node if the index is UNSET_INT */
    if (i < sv->size && sv->index[i] == UNSET_INT)
        sv->index[i] = ind;

    if (p == 1) {
        /* adds the entry if it cannot be found */
        if (i == sv->size) {
            sv->size++;
            spv_alloc(sv, max_nsys_sq);
            sv->index[i] = ind;
            sv->value[i] = 0.0;
        }

        /* adds the value to the vector */
        sv->value[i] += value[0];
    }
    else if (p == 4) {
        /* sets the position for the values */
        j = i * p2;

        /* adds the entry if it cannot be found */
        if (i == sv->size) {
            sv->size++;
            spv_alloc(sv, max_nsys_sq);
            sv->index[i] = ind;
            sv->value[j] = 0.0;
            sv->value[j + 1] = 0.0;
            sv->value[j + 2] = 0.0;
            sv->value[j + 3] = 0.0;
            sv->value[j + 4] = 0.0;
            sv->value[j + 5] = 0.0;
            sv->value[j + 6] = 0.0;
            sv->value[j + 7] = 0.0;
            sv->value[j + 8] = 0.0;
            sv->value[j + 9] = 0.0;
            sv->value[j + 10] = 0.0;
            sv->value[j + 11] = 0.0;
            sv->value[j + 12] = 0.0;
            sv->value[j + 13] = 0.0;
            sv->value[j + 14] = 0.0;
            sv->value[j + 15] = 0.0;
        }

        /* adds the value to the vector */
        sv->value[j] += value[0];
        sv->value[j + 1] += value[1];
        sv->value[j + 2] += value[2];
        sv->value[j + 3] += value[3];
        sv->value[j + 4] += value[4];
        sv->value[j + 5] += value[5];
        sv->value[j + 6] += value[6];
        sv->value[j + 7] += value[7];
        sv->value[j + 8] += value[8];
        sv->value[j + 9] += value[9];
        sv->value[j + 10] += value[10];
        sv->value[j + 11] += value[11];
        sv->value[j + 12] += value[12];
        sv->value[j + 13] += value[13];
        sv->value[j + 14] += value[14];
        sv->value[j + 15] += value[15];
    }
    else if (p == 3) {
        /* sets the position for the values */
        j = i * p2;

        /* adds the entry if it cannot be found */
        if (i == sv->size) {
            sv->size++;
            spv_alloc(sv, max_nsys_sq);
            sv->index[i] = ind;
            sv->value[j] = 0.0;
            sv->value[j + 1] = 0.0;
            sv->value[j + 2] = 0.0;
            sv->value[j + 3] = 0.0;
            sv->value[j + 4] = 0.0;
            sv->value[j + 5] = 0.0;
            sv->value[j + 6] = 0.0;
            sv->value[j + 7] = 0.0;
            sv->value[j + 8] = 0.0;
        }

        /* adds the value to the vector */
        sv->value[j] += value[0];
        sv->value[j + 1] += value[1];
        sv->value[j + 2] += value[2];
        sv->value[j + 3] += value[3];
        sv->value[j + 4] += value[4];
        sv->value[j + 5] += value[5];
        sv->value[j + 6] += value[6];
        sv->value[j + 7] += value[7];
        sv->value[j + 8] += value[8];
    }

    else if (p == 2) {
        /* sets the position for the values */
        j = i * p2;

        /* adds the entry if it cannot be found */
        if (i == sv->size) {
            sv->size++;
            spv_alloc(sv, max_nsys_sq);
            sv->index[i] = ind;
            sv->value[j] = 0.0;
            sv->value[j + 1] = 0.0;
            sv->value[j + 2] = 0.0;
            sv->value[j + 3] = 0.0;
        }

        /* adds the value to the vector */
        sv->value[j] += value[0];
        sv->value[j + 1] += value[1];
        sv->value[j + 2] += value[2];
        sv->value[j + 3] += value[3];
    }

    else
        tl_error("Spv_load has not been written for the desired system of equations.");
    return;
}

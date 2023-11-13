/* computes the dot product between the given sparse vector and the given full vector */

#include "global_header.h"

void spv_dot(SPARSE_VECT sv,    /* the sparse vector */
             double *fv,        /* the full vector */
             double *product,   /* the product */
             int p              /* the number of equations being solved */
    )
{
    int i, j;                   /* loop counter */
    int id0, id11, id12, id13, id14, id21, id22, id23, id24, id31, id32, id33, id34, id41, id42, id43, id44;    /* indices for multiple equation multiplication */
    int i1, i2, i3, i4;         /* indices for multiple equation multiplication */

    /* splits over the number of equations */
    if (p == 1) {
        /* initializes the result */
        product[0] = 0.0;

        /* computes the dot product */
        for (i = 0; i < sv.size; i++)
            product[0] += sv.value[i] * fv[sv.index[i]];
    }
    else if (p == 4) {
        /* initializes the result */
        product[0] = 0.0;
        product[1] = 0.0;
        product[2] = 0.0;
        product[3] = 0.0;

        /* computes the dot product */
        for (j = 0, id0 = 0; j < sv.size; j++, id0 += 16) {
            i1 = (sv.index[j]) * p;
            i2 = i1 + 1;
            i3 = i1 + 2;
            i4 = i1 + 3;
            id11 = id0;
            id12 = id0 + 1;
            id13 = id0 + 2;
            id14 = id0 + 3;
            id21 = id0 + 4;
            id22 = id0 + 5;
            id23 = id0 + 6;
            id24 = id0 + 7;
            id31 = id0 + 8;
            id32 = id0 + 9;
            id33 = id0 + 10;
            id34 = id0 + 11;
            id41 = id0 + 12;
            id42 = id0 + 13;
            id43 = id0 + 14;
            id44 = id0 + 15;
            product[0] += sv.value[id11] * fv[i1] + sv.value[id12] * fv[i2] + sv.value[id13] * fv[i3] + sv.value[id14] * fv[i4];
            product[1] += sv.value[id21] * fv[i1] + sv.value[id22] * fv[i2] + sv.value[id23] * fv[i3] + sv.value[id24] * fv[i4];
            product[2] += sv.value[id31] * fv[i1] + sv.value[id32] * fv[i2] + sv.value[id33] * fv[i3] + sv.value[id34] * fv[i4];
            product[3] += sv.value[id41] * fv[i1] + sv.value[id42] * fv[i2] + sv.value[id43] * fv[i3] + sv.value[id44] * fv[i4];
        }
    }
    else if (p == 3) {
        /* initializes the result */
        product[0] = 0.0;
        product[1] = 0.0;
        product[2] = 0.0;

        /* computes the dot product */
        for (j = 0, id0 = 0; j < sv.size; j++, id0 += p * p) {
            i1 = (sv.index[j]) * p;
            i2 = i1 + 1;
            i3 = i1 + 2;
            id11 = id0;
            id12 = id0 + 1;
            id13 = id0 + 2;
            id21 = id0 + 3;
            id22 = id0 + 4;
            id23 = id0 + 5;
            id31 = id0 + 6;
            id32 = id0 + 7;
            id33 = id0 + 8;
            product[0] += sv.value[id11] * fv[i1] + sv.value[id12] * fv[i2] + sv.value[id13] * fv[i3];
            product[1] += sv.value[id21] * fv[i1] + sv.value[id22] * fv[i2] + sv.value[id23] * fv[i3];
            product[2] += sv.value[id31] * fv[i1] + sv.value[id32] * fv[i2] + sv.value[id33] * fv[i3];
        }
    }
    else if (p == 2) {
        /* initializes the result */
        product[0] = 0.0;
        product[1] = 0.0;

        /* computes the dot product */
        for (j = 0, id0 = 0; j < sv.size; j++, id0 += p * p) {
            i1 = (sv.index[j]) * p;
            i2 = i1 + 1;
            id11 = id0;
            id12 = id0 + 1;
            id21 = id0 + 2;
            id22 = id0 + 3;
            product[0] += sv.value[id11] * fv[i1] + sv.value[id12] * fv[i2];
            product[1] += sv.value[id21] * fv[i1] + sv.value[id22] * fv[i2];
        }
    }

}

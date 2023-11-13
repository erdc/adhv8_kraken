/*!
   \file spv_vscale.c
   \brief Scales the given sparse (SPV) vector

   Scales the given sparse (SPV) vector by the entries 
   in the given full vector and the given scalar 
 */

#include "global_header.h"

void spv_vscale(SPARSE_VECT * sv,   /* the sparse vector */
                double *fv,     /* the full vector */
                double *s,      /* the short vector */
                int p           /* the number of equations being solved */
    )
{
    int i, j;                   /* loop counters */
    int i1, i2, i3, i4;         /* indices for multiple equation multiplication */
    int id0, id11, id12, id13, id14, id21, id22, id23, id24, id31, id32, id33, id34, id41, id42, id43, id44;    /* indices for multiple equation multiplication */

    /* splits over the number of equations being solved */
    if (p == 1) {
        for (i = 0; i < sv->size; i++)
            sv->value[i] *= s[0] * fv[sv->index[i]];
    }
    else if (p == 4) {
        for (j = 0, id0 = 0; j < sv->size; j++, id0 += 16) {
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
            i1 = (sv->index[j]) * p;
            i2 = i1 + 1;
            i3 = i1 + 2;
            i4 = i1 + 3;
            sv->value[id11] *= s[0] * fv[i1];
            sv->value[id12] *= s[0] * fv[i2];
            sv->value[id13] *= s[0] * fv[i3];
            sv->value[id14] *= s[0] * fv[i4];
            sv->value[id21] *= s[1] * fv[i1];
            sv->value[id22] *= s[1] * fv[i2];
            sv->value[id23] *= s[1] * fv[i3];
            sv->value[id24] *= s[1] * fv[i4];
            sv->value[id31] *= s[2] * fv[i1];
            sv->value[id32] *= s[2] * fv[i2];
            sv->value[id33] *= s[2] * fv[i3];
            sv->value[id34] *= s[2] * fv[i4];
            sv->value[id41] *= s[3] * fv[i1];
            sv->value[id42] *= s[3] * fv[i2];
            sv->value[id43] *= s[3] * fv[i3];
            sv->value[id44] *= s[3] * fv[i4];
        }
    }
    else if (p == 3) {
        for (j = 0, id0 = 0; j < sv->size; j++, id0 += p * p) {
            id11 = id0;
            id12 = id0 + 1;
            id13 = id0 + 2;
            id21 = id0 + 3;
            id22 = id0 + 4;
            id23 = id0 + 5;
            id31 = id0 + 6;
            id32 = id0 + 7;
            id33 = id0 + 8;
            i1 = (sv->index[j]) * p;
            i2 = i1 + 1;
            i3 = i1 + 2;
            sv->value[id11] *= s[0] * fv[i1];
            sv->value[id12] *= s[0] * fv[i2];
            sv->value[id13] *= s[0] * fv[i3];
            sv->value[id21] *= s[1] * fv[i1];
            sv->value[id22] *= s[1] * fv[i2];
            sv->value[id23] *= s[1] * fv[i3];
            sv->value[id31] *= s[2] * fv[i1];
            sv->value[id32] *= s[2] * fv[i2];
            sv->value[id33] *= s[2] * fv[i3];
        }
    }
    else if (p == 2) {
        for (j = 0, id0 = 0; j < sv->size; j++, id0 += p * p) {
            id11 = id0;
            id12 = id0 + 1;
            id21 = id0 + 2;
            id22 = id0 + 3;
            i1 = (sv->index[j]) * p;
            i2 = i1 + 1;
            sv->value[id11] *= s[0] * fv[i1];
            sv->value[id12] *= s[0] * fv[i2];
            sv->value[id21] *= s[1] * fv[i1];
            sv->value[id22] *= s[1] * fv[i2];
        }
    }

    else
        tl_error("Spv_vscale has not been written for the desired system of equations.");
}

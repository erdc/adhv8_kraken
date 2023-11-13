/*!
   \file spv_jacobi.c
   \brief Implements point Jacobi preconditioning for SPV form

   implements point Jacobi preconditioning on 'matrix', 'diagonal', 
   B, and x this routine executes [S][A][S][y]=[S][B]
   where i
   [S] is a diagonal matrix of scale factors
   [A] is the coefficient matrix composed of 'diagonal' and 'matrix'
   [Y] = [Inverse S] [x]
   [B] = the residual
   we now are solving for [Y] instead of [x]
   this has to be undone to find [x] after the solve for [Y]  
 */

#include "global_header.h"
/* #include "solv_type.h" */

void spv_jacobi(SPARSE_VECT * matrix,   /* the matrix */
                double *diagonal,   /* the diagonal */
                double *rhs,    /* the right hand side */
                double *sol,    /* the initial guess */
                double *scale_vect, /* the scale vector */
                int n,          /* the number of degrees of freedom */
                int m,          /* the number of nodes */
                int p           /* the number of equations being solved */
#ifdef _MESSG
                ,SMPI *smpi
#endif
    )
{
    double dnorm[4];            /* the norm of the row */
    double row_factor[4];       /* the row factor */
    double factor;              /* coefficient to reduce memory fetches */
    double small = SOLV_TOL;    /* too small a number for a norm */
    double small_factor;        /* small^(-0.5) */
    int i, j;                   /* loop counter */
    int id0, id11, id12, id13, id14, id21, id22, id23, id24, id31, id32, id33, id34, id41, id42, id43, id44;    /* indices for multiple equation multiplication */
    int k1, k2, k3, k4;         /* indices for multiple equation multiplication */

    dnorm[0] = DZERO;
    dnorm[1] = DZERO;
    dnorm[2] = DZERO;
    dnorm[3] = DZERO;

    /* calculates the small factor */
    small_factor = 1.0 / sqrt(small);

    /* splits over the number of equations */
    if (p == 1) {

        /* constructs the scale vector */
        for (i = 0; i < m; i++)
            scale_vect[i] = fabs(diagonal[i]);
        for (i = 0; i < m; i++) {
            spv_norm(matrix[i], dnorm, p);
            if (dnorm[0] > scale_vect[i])
                scale_vect[i] = dnorm[0];
        }
#ifdef _MESSG
        comm_update_double(scale_vect, p, smpi);
#endif

        /* rolls the scale vector into the inverse of the square root of the scale vector */
        for (i = 0; i < n; i++)
            if (scale_vect[i] < small)
                scale_vect[i] = small_factor;
            else
                scale_vect[i] = 1.0 / sqrt(fabs(scale_vect[i]));

        /* scales the equations */
        for (i = 0; i < n; i++) {
            factor = scale_vect[i];
            rhs[i] *= factor;
            sol[i] /= factor;
        }
        for (i = 0; i < m; i++) {
            row_factor[0] = scale_vect[i];
            diagonal[i] *= row_factor[0] * row_factor[0];
            spv_vscale(matrix + i, scale_vect, row_factor, p);
        }
        return;
    }
    else if (p == 4) {
        /* constructs the scale vector */
        /* first find the largest term in the row of the diagonal block */
        for (i = 0, j = 0, id0 = 0; i < m; i++, j += 4, id0 += 16) {
            k1 = j;
            k2 = j + 1;
            k3 = j + 2;
            k4 = j + 3;
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
            scale_vect[k1] = DZERO;
            scale_vect[k2] = DZERO;
            scale_vect[k3] = DZERO;
            scale_vect[k4] = DZERO;
            if (fabs(diagonal[id11]) > scale_vect[k1])
                scale_vect[k1] = fabs(diagonal[id11]);
            if (fabs(diagonal[id12]) > scale_vect[k1])
                scale_vect[k1] = fabs(diagonal[id12]);
            if (fabs(diagonal[id13]) > scale_vect[k1])
                scale_vect[k1] = fabs(diagonal[id13]);
            if (fabs(diagonal[id14]) > scale_vect[k1])
                scale_vect[k1] = fabs(diagonal[id14]);
            if (fabs(diagonal[id21]) > scale_vect[k2])
                scale_vect[k2] = fabs(diagonal[id21]);
            if (fabs(diagonal[id22]) > scale_vect[k2])
                scale_vect[k2] = fabs(diagonal[id22]);
            if (fabs(diagonal[id23]) > scale_vect[k2])
                scale_vect[k2] = fabs(diagonal[id23]);
            if (fabs(diagonal[id24]) > scale_vect[k2])
                scale_vect[k2] = fabs(diagonal[id24]);
            if (fabs(diagonal[id31]) > scale_vect[k3])
                scale_vect[k3] = fabs(diagonal[id31]);
            if (fabs(diagonal[id32]) > scale_vect[k3])
                scale_vect[k3] = fabs(diagonal[id32]);
            if (fabs(diagonal[id33]) > scale_vect[k3])
                scale_vect[k3] = fabs(diagonal[id33]);
            if (fabs(diagonal[id34]) > scale_vect[k3])
                scale_vect[k3] = fabs(diagonal[id34]);
            if (fabs(diagonal[id41]) > scale_vect[k4])
                scale_vect[k4] = fabs(diagonal[id41]);
            if (fabs(diagonal[id42]) > scale_vect[k4])
                scale_vect[k4] = fabs(diagonal[id42]);
            if (fabs(diagonal[id43]) > scale_vect[k4])
                scale_vect[k4] = fabs(diagonal[id43]);
            if (fabs(diagonal[id44]) > scale_vect[k4])
                scale_vect[k4] = fabs(diagonal[id44]);
        }
        /* Now search the off-diagonal block for any terms in the row that are larger */
        for (i = 0, j = 0; i < m; i++, j += 4) {
            k1 = j;
            k2 = j + 1;
            k3 = j + 2;
            k4 = j + 3;
            spv_norm(matrix[i], dnorm, p);
            if (dnorm[0] > scale_vect[k1])
                scale_vect[k1] = dnorm[0];
            if (dnorm[1] > scale_vect[k2])
                scale_vect[k2] = dnorm[1];
            if (dnorm[2] > scale_vect[k3])
                scale_vect[k3] = dnorm[2];
            if (dnorm[3] > scale_vect[k4])
                scale_vect[k4] = dnorm[3];
        }
#ifdef _MESSG
        comm_update_double(scale_vect, p, smpi);
#endif

        /* rolls the scale vector into the inverse of the square root of the scale vector */
        for (i = 0; i < n; i++)
            if (scale_vect[i] < small)
                scale_vect[i] = small_factor;
            else
                scale_vect[i] = 1.0 / sqrt(fabs(scale_vect[i]));

        /* scales the equations */
        for (i = 0; i < n; i++) {
            factor = scale_vect[i];
            rhs[i] *= factor;
            sol[i] /= factor;
        }
        for (i = 0, j = 0, id0 = 0; i < m; i++, j += 4, id0 += 16) {
            k1 = j;
            k2 = j + 1;
            k3 = j + 2;
            k4 = j + 3;
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
            row_factor[0] = scale_vect[k1];
            row_factor[1] = scale_vect[k2];
            row_factor[2] = scale_vect[k3];
            row_factor[3] = scale_vect[k4];
            diagonal[id11] *= row_factor[0] * row_factor[0];
            diagonal[id12] *= row_factor[0] * row_factor[1];
            diagonal[id13] *= row_factor[0] * row_factor[2];
            diagonal[id14] *= row_factor[0] * row_factor[3];
            diagonal[id21] *= row_factor[1] * row_factor[0];
            diagonal[id22] *= row_factor[1] * row_factor[1];
            diagonal[id23] *= row_factor[1] * row_factor[2];
            diagonal[id24] *= row_factor[1] * row_factor[3];
            diagonal[id31] *= row_factor[2] * row_factor[0];
            diagonal[id32] *= row_factor[2] * row_factor[1];
            diagonal[id33] *= row_factor[2] * row_factor[2];
            diagonal[id34] *= row_factor[2] * row_factor[3];
            diagonal[id41] *= row_factor[3] * row_factor[0];
            diagonal[id42] *= row_factor[3] * row_factor[1];
            diagonal[id43] *= row_factor[3] * row_factor[2];
            diagonal[id44] *= row_factor[3] * row_factor[3];
            spv_vscale(matrix + i, scale_vect, row_factor, p);
        }
    }
    else if (p == 3) {
        /* constructs the scale vector */
        /* first find the largest term in the row of the diagonal block */
        for (i = 0, j = 0, id0 = 0; i < m; i++, j += p, id0 += p * p) {
            k1 = j;
            k2 = j + 1;
            k3 = j + 2;
            id11 = id0;
            id12 = id0 + 1;
            id13 = id0 + 2;
            id21 = id0 + 3;
            id22 = id0 + 4;
            id23 = id0 + 5;
            id31 = id0 + 6;
            id32 = id0 + 7;
            id33 = id0 + 8;
            scale_vect[k1] = DZERO;
            scale_vect[k2] = DZERO;
            scale_vect[k3] = DZERO;
            if (fabs(diagonal[id11]) > scale_vect[k1])
                scale_vect[k1] = fabs(diagonal[id11]);
            if (fabs(diagonal[id12]) > scale_vect[k1])
                scale_vect[k1] = fabs(diagonal[id12]);
            if (fabs(diagonal[id13]) > scale_vect[k1])
                scale_vect[k1] = fabs(diagonal[id13]);
            if (fabs(diagonal[id21]) > scale_vect[k2])
                scale_vect[k2] = fabs(diagonal[id21]);
            if (fabs(diagonal[id22]) > scale_vect[k2])
                scale_vect[k2] = fabs(diagonal[id22]);
            if (fabs(diagonal[id23]) > scale_vect[k2])
                scale_vect[k2] = fabs(diagonal[id23]);
            if (fabs(diagonal[id31]) > scale_vect[k3])
                scale_vect[k3] = fabs(diagonal[id31]);
            if (fabs(diagonal[id32]) > scale_vect[k3])
                scale_vect[k3] = fabs(diagonal[id32]);
            if (fabs(diagonal[id33]) > scale_vect[k3])
                scale_vect[k3] = fabs(diagonal[id33]);
        }
        /* Now search the off-diagonal block for any terms in the row that are larger */
        for (i = 0, j = 0; i < m; i++, j += p) {
            k1 = j;
            k2 = j + 1;
            k3 = j + 2;
            spv_norm(matrix[i], dnorm, p);
            if (dnorm[0] > scale_vect[k1])
                scale_vect[k1] = dnorm[0];
            if (dnorm[1] > scale_vect[k2])
                scale_vect[k2] = dnorm[1];
            if (dnorm[2] > scale_vect[k3])
                scale_vect[k3] = dnorm[2];
        }
#ifdef _MESSG
        comm_update_double(scale_vect, p, smpi);
#endif

        /* rolls the scale vector into the inverse of the square root of the scale vector */
        for (i = 0; i < n; i++)
            if (scale_vect[i] < small)
                scale_vect[i] = small_factor;
            else
                scale_vect[i] = 1.0 / sqrt(fabs(scale_vect[i]));

        /* scales the equations */
        for (i = 0; i < n; i++) {
            factor = scale_vect[i];
            rhs[i] *= factor;
            sol[i] /= factor;
        }
        for (i = 0, j = 0, id0 = 0; i < m; i++, j += p, id0 += p * p) {
            k1 = j;
            k2 = j + 1;
            k3 = j + 2;
            id11 = id0;
            id12 = id0 + 1;
            id13 = id0 + 2;
            id21 = id0 + 3;
            id22 = id0 + 4;
            id23 = id0 + 5;
            id31 = id0 + 6;
            id32 = id0 + 7;
            id33 = id0 + 8;
            row_factor[0] = scale_vect[k1];
            row_factor[1] = scale_vect[k2];
            row_factor[2] = scale_vect[k3];
            diagonal[id11] *= row_factor[0] * row_factor[0];
            diagonal[id12] *= row_factor[0] * row_factor[1];
            diagonal[id13] *= row_factor[0] * row_factor[2];
            diagonal[id21] *= row_factor[1] * row_factor[0];
            diagonal[id22] *= row_factor[1] * row_factor[1];
            diagonal[id23] *= row_factor[1] * row_factor[2];
            diagonal[id31] *= row_factor[2] * row_factor[0];
            diagonal[id32] *= row_factor[2] * row_factor[1];
            diagonal[id33] *= row_factor[2] * row_factor[2];
            spv_vscale(matrix + i, scale_vect, row_factor, p);
        }
    }
    else if (p == 2) {
        /* constructs the scale vector */
        /* first find the largest term in the row of the diagonal block */
        for (i = 0, j = 0, id0 = 0; i < m; i++, j += p, id0 += p * p) {
            k1 = j;
            k2 = j + 1;
            id11 = id0;
            id12 = id0 + 1;
            id21 = id0 + 2;
            id22 = id0 + 3;
            scale_vect[k1] = DZERO;
            scale_vect[k2] = DZERO;
            if (fabs(diagonal[id11]) > scale_vect[k1])
                scale_vect[k1] = fabs(diagonal[id11]);
            if (fabs(diagonal[id12]) > scale_vect[k1])
                scale_vect[k1] = fabs(diagonal[id12]);
            if (fabs(diagonal[id21]) > scale_vect[k2])
                scale_vect[k2] = fabs(diagonal[id21]);
            if (fabs(diagonal[id22]) > scale_vect[k2])
                scale_vect[k2] = fabs(diagonal[id22]);
        }
        /* Now search the off-diagonal block for any terms in the row that are larger */
        for (i = 0, j = 0; i < m; i++, j += p) {
            k1 = j;
            k2 = j + 1;
            spv_norm(matrix[i], dnorm, p);
            if (dnorm[0] > scale_vect[k1])
                scale_vect[k1] = dnorm[0];
            if (dnorm[1] > scale_vect[k2])
                scale_vect[k2] = dnorm[1];
        }
#ifdef _MESSG
        comm_update_double(scale_vect, p, smpi);
#endif

        /* rolls the scale vector into the inverse of the square root of the scale vector */
        for (i = 0; i < n; i++)
            if (scale_vect[i] < small)
                scale_vect[i] = small_factor;
            else
                scale_vect[i] = 1.0 / sqrt(fabs(scale_vect[i]));

        /* scales the equations */
        for (i = 0; i < n; i++) {
            factor = scale_vect[i];
            rhs[i] *= factor;
            sol[i] /= factor;
        }
        for (i = 0, j = 0, id0 = 0; i < m; i++, j += p, id0 += p * p) {
            k1 = j;
            k2 = j + 1;
            id11 = id0;
            id12 = id0 + 1;
            id21 = id0 + 2;
            id22 = id0 + 3;
            row_factor[0] = scale_vect[k1];
            row_factor[1] = scale_vect[k2];
            diagonal[id11] *= row_factor[0] * row_factor[0];
            diagonal[id12] *= row_factor[0] * row_factor[1];
            diagonal[id21] *= row_factor[1] * row_factor[0];
            diagonal[id22] *= row_factor[1] * row_factor[1];
            spv_vscale(matrix + i, scale_vect, row_factor, p);
        }
    }
    else
        tl_error("Solv_jacobi has not been written for the desired system of equations.");
}

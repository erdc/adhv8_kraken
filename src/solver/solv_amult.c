/* returns the matrix vector product */

#include "global_header.h"

void solv_amult(Solver_Info *solver,
                double *v,      /* the vector */
                double *av,     /* the matrix times vector */
                double *diagonal,   /* the diagonal */
                SPARSE_VECT * matrix,   /* the matrix */
                int m,          /* the number of nodes */
                int p           /* the number of equations being solved */
#ifdef _MESSG
                ,SMPI *smpi
#endif
    )
{
    int i, j;                   /* loop counter */
    int k1, k2, k3, k4;         /* indices for multiple equation multiplication */
    int id0, id11, id12, id13, id14, id21, id22, id23, id24, id31, id32, id33, id34, id41, id42, id43, id44;    /* indices for multiple equation multiplication */
    double intermediate_result[4];  /* the result of the row contribution */

    /* split over the number of equations being solved */
    if (p == 1) {
        /* update v */
#ifdef _MESSG
        comm_update_double(v, p, smpi);
#endif

        /* computes the diagonal contribution - the identity since diagonal scaling is used */
        for (i = 0; i < m; i++)
            av[i] = diagonal[i] * v[i];

        /* loops over the nodes computing the matrix vector product row by row */
        for (i = 0; i < m; i++) {
            spv_dot(matrix[i], v, intermediate_result, p);
            av[i] += intermediate_result[0];
        }

        /* update av - this may not be needed */
#ifdef _MESSG
        comm_update_double(av, p, smpi);
#endif
    }
    else if (p == 4) {
        /* update v */
#ifdef _MESSG
        comm_update_double(v, p, smpi);
#endif

        /* computes the diagonal contribution - the identity since diagonal scaling is used */
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
            av[k1] = diagonal[id11] * v[k1] + diagonal[id12] * v[k2] + diagonal[id13] * v[k3] + diagonal[id14] * v[k4];
            av[k2] = diagonal[id21] * v[k1] + diagonal[id22] * v[k2] + diagonal[id23] * v[k3] + diagonal[id24] * v[k4];
            av[k3] = diagonal[id31] * v[k1] + diagonal[id32] * v[k2] + diagonal[id33] * v[k3] + diagonal[id34] * v[k4];
            av[k4] = diagonal[id41] * v[k1] + diagonal[id42] * v[k2] + diagonal[id43] * v[k3] + diagonal[id44] * v[k4];
        }

        /* loops over the nodes computing the matrix vector product row by row */
        for (i = 0, j = 0; i < m; i++, j += 4) {
            spv_dot(matrix[i], v, intermediate_result, p);
            av[j] += intermediate_result[0];
            av[j + 1] += intermediate_result[1];
            av[j + 2] += intermediate_result[2];
            av[j + 3] += intermediate_result[3];
        }

        /* update av - this may not be needed */
#ifdef _MESSG
        comm_update_double(av, p, smpi);
#endif
    }
    else if (p == 3) {
        /* update v */
#ifdef _MESSG
        comm_update_double(v, p, smpi);
#endif

        /* computes the diagonal contribution - the identity since diagonal scaling is used */
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
            av[k1] = diagonal[id11] * v[k1] + diagonal[id12] * v[k2] + diagonal[id13] * v[k3];
            av[k2] = diagonal[id21] * v[k1] + diagonal[id22] * v[k2] + diagonal[id23] * v[k3];
            av[k3] = diagonal[id31] * v[k1] + diagonal[id32] * v[k2] + diagonal[id33] * v[k3];
        }

        /* loops over the nodes computing the matrix vector product row by row */
        for (i = 0, j = 0; i < m; i++, j += p) {
            spv_dot(matrix[i], v, intermediate_result, p);
            av[j] += intermediate_result[0];
            av[j + 1] += intermediate_result[1];
            av[j + 2] += intermediate_result[2];
        }

        /* update av - this may not be needed */
#ifdef _MESSG
        comm_update_double(av, p, smpi);
#endif
    }
    else if (p == 2) {
        /* update v */
#ifdef _MESSG
        comm_update_double(v, p, smpi);
#endif

        /* computes the diagonal contribution - the identity since diagonal scaling is used */
        for (i = 0, j = 0, id0 = 0; i < m; i++, j += p, id0 += p * p) {
            k1 = j;
            k2 = j + 1;
            id11 = id0;
            id12 = id0 + 1;
            id21 = id0 + 2;
            id22 = id0 + 3;
            av[k1] = diagonal[id11] * v[k1] + diagonal[id12] * v[k2];
            av[k2] = diagonal[id21] * v[k1] + diagonal[id22] * v[k2];
        }

        /* loops over the nodes computing the matrix vector product row by row */
        for (i = 0, j = 0; i < m; i++, j += p) {
            spv_dot(matrix[i], v, intermediate_result, p);
            av[j] += intermediate_result[0];
            av[j + 1] += intermediate_result[1];
        }

        /* update av - this may not be needed */
#ifdef _MESSG
        comm_update_double(av, p, smpi);
#endif
    }

    else
        tl_error("Solv_amult has not been written for the desired system of equations.");
}

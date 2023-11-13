/* sets up block Jacobi / additive Schwarz preconditioning */

#include "global_header.h"

void solv_blk_set_profile(Solver_Info *solver,
                          Profile_Info *prof,
                          SPARSE_VECT * matrix, /* the matrix */
                          double *diagonal, /* the diagonal */
                          int blk_my_nnode, /* the number of nodes I own */
                          int blk_nsys, /* the number of equations being solved */
                          int blk_nsys_sq   /* the number of equations being solved squared */
    )
{
    int i;                      /* loop counters */
    int new_size;               /* the new size of the array */

    /* allocate space if needed */
    new_size = blk_my_nnode * blk_nsys;
    if (new_size > prof->size) {
        prof->rows = (BAND_VECT *) tl_realloc(sizeof(BAND_VECT), new_size, prof->size, prof->rows);
        prof->cols = (BAND_VECT *) tl_realloc(sizeof(BAND_VECT), new_size, prof->size, prof->cols);
        for (i = prof->size; i < new_size; i++) {
            prof->rows[i].begin = UNSET_INT;
            prof->rows[i].end = UNSET_INT;
            prof->rows[i].size = 0;
            prof->rows[i].value = NULL;
            prof->cols[i].begin = UNSET_INT;
            prof->cols[i].end = UNSET_INT;
            prof->cols[i].size = 0;
            prof->cols[i].value = NULL;
        }
        prof->size = new_size;
    }

    /* initialize the arrays */
    for (i = 0; i < prof->size; i++) {
        bv_init(prof->rows + i);
        bv_init(prof->cols + i);
    }

    /* load the profile matrix */
    solv_blk_load_profile(solver->node_block, matrix, diagonal, prof->rows, prof->cols, blk_my_nnode, blk_nsys, blk_nsys_sq);

    /* calculate the LU factorization of the matrix */
    solv_blk_factor_profile(prof->rows, prof->cols, blk_nsys * blk_my_nnode);

}

/* replaces the profile matrix with the LU factorization of the matrix */

void solv_blk_factor_profile(BAND_VECT * rows,  /* the rows of the matrix */
                             BAND_VECT * cols,  /* the columns of the matrix */
                             int neq    /* the number of equations */
    )
{
    int i, j, k;                /* loop counters */
    int ip;                     /* the profile counter */
    double pivot;               /* the pivot the equation is scaled by */

    /* loop through the profiles */
    for (ip = 0; ip < neq; ip++) {
        /* calculate the lower triangular entries */
        for (j = rows[ip].begin, k = 0; j < rows[ip].end; j++, k++) {

            rows[ip].value[k] -= bv_dot(rows + ip, cols + j, 0, j);
            pivot = cols[j].value[cols[j].end - cols[j].begin - 1];
            if (fabs(pivot) > SMALL)
                rows[ip].value[k] /= pivot;
            else {
                rows[ip].value[k] /= SMALL;
            }
        }

        /* calculate the upper triangular entries */
        for (i = cols[ip].begin, k = 0; i < cols[ip].end; i++, k++) {
            cols[ip].value[k] -= bv_dot(rows + i, cols + ip, 0, i);
        }
    }
}

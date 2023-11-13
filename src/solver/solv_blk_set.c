/* sets up block Jacobi / additive Schwarz preconditioning */

#include "global_header.h"

void solv_blk_set(Solver_Info *solver,
                  SPARSE_VECT * matrix, /* the matrix */
                  double *diagonal, /* the diagonal */
                  int blk_my_nnode, /* the number of nodes I own */
                  int nsys,     /* the number of equations being solved */
                  int nsys_sq   /* the number of equations being solved squared */
#ifdef _MESSG
                  , MPI_Comm ADH_COMM
#endif
    )
{
    static int factored_flag = 0;   /* when this flag ==1, the sparse data structures will be freed before reloading and refactoring preconditioners */

    solver->UMFail = NO;

#ifdef _UMFPACK
    if (factored_flag == 1)
        solv_blk_free_sparse(); /* free fine grid */
#ifdef _MESSG
    solver->UMFail = solv_blk_set_sparse(solver->node_block, matrix, diagonal, blk_my_nnode, nsys, nsys_sq, ADH_COMM);
#else
    solver->UMFail = solv_blk_set_sparse(solver->node_block, matrix, diagonal, blk_my_nnode, nsys, nsys_sq);
#endif
#else
    solv_blk_set_profile(solver, &(solver->profile), matrix, diagonal, blk_my_nnode, nsys, nsys_sq);
#endif
    factored_flag = 1;

}

void solv_blk_set_clean(Profile_Info *prof)
{
    /* free memory allocated in prof.rows and prof.cols  */
    int i;

    if (prof->rows != NULL) {
        for (i = 0; i < prof->size; i++) {
            if (prof->rows[i].value != NULL) {
                prof->rows[i].value = (double *) tl_free(sizeof(double), prof->rows[i].size, prof->rows[i].value);
            }
        }
        prof->rows = (BAND_VECT *) tl_free(sizeof(BAND_VECT), prof->size, prof->rows);
    }

    if (prof->cols != NULL) {
        for (i = 0; i < prof->size; i++) {
            if (prof->cols[i].value != NULL) {
                prof->cols[i].value = (double *) tl_free(sizeof(double), prof->cols[i].size, prof->cols[i].value);
            }
        }
        prof->cols = (BAND_VECT *) tl_free(sizeof(BAND_VECT), prof->size, prof->cols);
    }

}

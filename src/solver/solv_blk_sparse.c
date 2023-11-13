#ifdef _UMFPACK

/* Version 1.0 */
/* sets up block Jacobi and additive Schwarz preconditioning */
/* uses UMFPACK */

#include "global_header.h"

/* UMFPACK data structures */

static long *ia = NULL;         /* ia[i] points to first coefficient in ith
                                 * row */
static long *ja = NULL;         /* ja[j] column index of coefficient */
static double *a = NULL;        /* a[j]  numeric value of coefficient */
static int nnz, nnz_alloc;
static int n;
static double *xtemp = NULL;

/* load and factor the matrix */

int solv_blk_set_sparse(int *node_block, SPARSE_VECT * matrix,  /* off-diagonal part of global sparse matrix */
                        double *diagonal,   /* diagonal of global sparse matrix */
                        int blk_my_nnode,   /* the number of nodes on pe */
                        int blk_nsys,   /* the number of equations per node  */
                        int blk_nsys_sq /* squared */
#ifdef _MESSG
                        , MPI_Comm ADH_COMM
#endif
    )
{
    int i, j, ii, jj;
    int UMFail = NO;
    SPARSE_VECT *my_sparse_row;

    n = blk_my_nnode * blk_nsys;

    /* find number of nonzeros */

    for (nnz = 0, i = 0; i < blk_my_nnode; i++) {
        my_sparse_row = matrix + i;
        for (ii = 0; ii < blk_nsys; ii++) {
            for (jj = 0; jj < blk_nsys; jj++)
                if (fabs(diagonal[i * blk_nsys_sq + ii * blk_nsys + jj]) > SMALL)
                    nnz++;
            for (j = 0; j < my_sparse_row->size; j++)
                if (node_block[my_sparse_row->index[j]] == node_block[i])
                    for (jj = 0; jj < blk_nsys; jj++)
                        if (fabs(my_sparse_row->value[j * blk_nsys_sq + ii * blk_nsys + jj]) > SMALL)
                            nnz++;
        }
    }
    /* allocate storage for the block matrix in triplet form */

    if (n < nnz)
        nnz_alloc = nnz;
    else
        nnz_alloc = n + 1;
    ia = (long *) tl_alloc(sizeof(long), nnz_alloc);
    ja = (long *) tl_alloc(sizeof(long), nnz_alloc);
    a = (double *) tl_alloc(sizeof(double), nnz_alloc);
    xtemp = (double *) tl_alloc(sizeof(double), nnz_alloc);

    /* load block matrix from global sparse matrix */
    /* on return, row indices in ia , column indices in ja */

    solv_blk_load_sparse(node_block, matrix, diagonal, blk_my_nnode, blk_nsys, blk_nsys_sq, n, nnz, ia, ja, a);

    /* LU factorization */
    /* on return, ia,ja,a in std compressed column format */
#ifdef _MESSG
    UMFail = solv_blk_UMFPACK_fact(n, nnz, ia, ja, a, "triplet", ADH_COMM);
#else
    UMFail = solv_blk_UMFPACK_fact(n, nnz, ia, ja, a, "triplet");
#endif
    return UMFail;
}

int solv_blk_solve_sparse       /* solves the LU factored profile
                                 * matrix */
    (double *sol,               /* the solution */
     int n                      /* the number of
                                 * equations */
#ifdef _MESSG
     , MPI_Comm ADH_COMM
#endif
    )
{
    int UMFail = NO;
#ifdef _MESSG
    UMFail = solv_blk_UMFPACK_solve(ia, ja, a, n, sol, xtemp, ADH_COMM);
#else
    UMFail = solv_blk_UMFPACK_solve(ia, ja, a, n, sol, xtemp);
#endif
    return UMFail;
}

void solv_blk_free_sparse()
{
    if (ia != NULL)
        ia = (long *) tl_free(sizeof(long), nnz_alloc, ia);
    if (ja != NULL)
        ja = (long *) tl_free(sizeof(long), nnz_alloc, ja);
    if (a != NULL)
        a = (double *) tl_free(sizeof(double), nnz_alloc, a);
    if (xtemp != NULL)
        xtemp = (double *) tl_free(sizeof(double), nnz_alloc, xtemp);
    solv_blk_UMFPACK_free();
}

#endif /* ifdef UMFPACK */

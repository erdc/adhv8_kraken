/* solves the LU factored profile matrix */

#include "global_header.h"

int solv_blk_solve(Profile_Info prof,
                    double *sol,    /* the solution */
                    int neq     /* the number of equations */
#ifdef _MESSG
                    , MPI_Comm ADH_COMM
#endif
    )
{
    int i, j, k;                /* loop counter */
    int iend;                   /* the end of the loop */
    double scale;               /* scale factor for back substitions */
    int UMFail = NO;

#ifdef _UMFPACK
#ifdef _MESSG
    UMFail = solv_blk_solve_sparse(sol, neq, ADH_COMM);
#else
    UMFail = solv_blk_solve_sparse(sol, neq);
#endif
    return UMFail;
#endif

    //printScreen_dble_array("solv_blk_solv", sol, neq, __FILE__, __LINE__);
    //printScreen_matrixProfile("matrix profile", &prof, __FILE__, __LINE__);

    /* forward substitution */
    for (k = 0; k < neq; k++) {
        /* find the end of the loop - caused by the matrix being profile rather than full */
        if (k > prof.rows[k].end)
            iend = prof.rows[k].end;
        else
            iend = k;

        /* loop thru the substitutions */
        for (i = prof.rows[k].begin, j = 0; i < iend; i++, j++)
            sol[k] -= prof.rows[k].value[j] * sol[i];
    }

    /* backward substitution - this may appear strange because the loops have 
       been reversed from the standard textbook (manual) method for back substitution - 
       this version accesses the elements of the columns sequentially in memory */
    for (k = neq - 1; k > -1; k--) {
        /* scale by the diagonal */
        sol[k] /= prof.cols[k].value[prof.cols[k].end - prof.cols[k].begin - 1];

        /* find the end of the loop - caused by the matrix being profile rather than full */
        if (k > prof.cols[k].end)
            iend = prof.cols[k].end;
        else
            iend = k;

        /* loops thru the substitutions */
        for (i = prof.cols[k].begin, scale = sol[k], j = 0; i < iend; i++, j++)
            sol[i] -= prof.cols[k].value[j] * scale;
    }

    return UMFail;
}

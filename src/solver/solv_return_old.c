/* strips off the right preconditioning diagonal */

#include "global_header.h"

void solv_return_old(Solver_Info *solver,    
                     double *solv_u,    /* the solution */
                     double *solv_u0,   /* the shift from the solver */
                     double *scale_vect,    /* the scale vector for the jacobi preconditioning */
                     int n,     /* the number of degrees of freedom */
                     int p      /* the number of equations being solved */
#ifdef _MESSG
                     ,SMPI *smpi
#endif
    )
{
    int i;                      /* loop counter */

    /* UNDO SCALING */
    if (solver->prec_value > 1)
        scale_vect = NULL;

    /* undoes the diagonal scaling and the shift */
    for (i = 0; i < n; i++) {
        solv_u[i] += solv_u0[i];
        /* solv_u[i] *= scale_vect[i]; UNDO SCALING */
    }
    if (solver->prec_value <= 1) {
        for (i = 0; i < n; i++)
            solv_u[i] *= scale_vect[i]; /* UNDO SCALING */
    }

    /* updates the solution */
#ifdef _MESSG
    comm_update_double(solv_u, p, smpi);
#else
    p = 0;
#endif
}

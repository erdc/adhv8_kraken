/*!
   \file solv_linear_system.c 
   \brief Actually solve the linear system.  

   This routine simply solves the linear system
   \f[A x = b\f]
   for \f$x\f$, where \f$x,b\f$ are vectors and \f$A\f$ 
   is a matrix.  The main procedure for doing this is the 
   iterative solver bicgstab, with preconditioning. 

   This needs to be called after solv_setup_system().
 */

#include "global_header.h"

int solv_linear_sys_solve(Solver_Info *solver,
                          int *bc_mask, /* the bc mask for each dof */
                          SPARSE_VECT * matrix, /* the matrix */
                          double *diagonal, /* the diagonal */
                          double *solv_b,   /* the right-hand-side vector */
                          double *solv_x,   /* the solution */
                          double *scale_vect,   /* temporary vector */
                          int mype_Ndof,    /* num. dof that my processor sees */
                          int mype_Ghost,   /* Ndof += ghosts, i.e. values on other pe's */
                          int Nsys_solving  /* the number of equations being solved */
#ifdef _MESSG
                          ,SMPI *smpi
#endif
    )
{
    int solv_flag = NO;
    int iterations = 0;

    if (solver->prec_value > 1) {
        tl_error("Preconditioner of 0 and 1 are only supported at this time.");
    }

  /************************************************************/
    /* Cases that call old BiCG-STAB */
    if ((solver->prec_value == 0) || (solver->prec_value == 1) || (solver->prec_value == 2) || (solver->prec_value == 3)) {
        solv_flag = solv_bcgstab(solver, bc_mask, matrix, diagonal, solv_b, solv_x, scale_vect, mype_Ghost, mype_Ndof, Nsys_solving
#ifdef _MESSG
                ,smpi
#endif
                );
    }
  /************************************************************/
    else {
    }
  /************************************************************/
  /************************************************************/
    /* Finish up Linear solver->Work */

    if (solver->max_lin_it == iterations) {
        if (solver->force_lin_it == YES)
            solv_flag = YES;
    }

    if (solv_flag == YES) {
        solver->it_count_lin += iterations;
        solver->lin_fail = NO;
    }
    else if (solv_flag == NO) {
        solver->lin_fail = YES;
        solver->it_count_lin_failed += iterations;
    }

    return (solv_flag);
}

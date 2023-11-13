/*!
 \file solv_setup_system.c
 \brief Setup the linear system before actual solve.
 
 This routine will actually set up the system before
 it is solved.  It should take care of:
 1. Conversion to/from different storage formats
 2. Getting the pre-conditioner set-up
 3. Possible "pre-solution" work, as might be done
 for LU factorizations.
 */

#include "global_header.h"

void solv_linear_sys_setup(Solver_Info *solver, int *bc_mask,  /* the bc mask for each dof */
                           SPARSE_VECT *matrix,    /* the matrix */
                           double *diagonal,    /* the diagonal */
                           double *solv_b,  /* the right-hand-side vector */
                           double *solv_x,  /* the solution */
                           double *scale_vect,  /* temporary vector */
                           int mype_Ndof,   /* num. dof that my processor sees */
                           int mype_Ghost,  /* Ndof += ghosts, i.e. values on other pe's */
                           int Nsys_solving /* the number of equations being solved */
#ifdef _MESSG
                           , SMPI *smpi
#endif
)
{
    if (solver->prec_value > 1) {
        tl_error("Preconditioner of 0 and 1 are only supported at this time.");
    }
    
    /************************************************************/
    /* Point Jacobi, Old Way */
    if (solver->prec_value == 0) {
        spv_jacobi(matrix, diagonal, solv_b, solv_x, scale_vect, Nsys_solving * mype_Ghost, mype_Ghost, Nsys_solving
#ifdef _MESSG
                   , smpi
#endif
                   );
    }
    /************************************************************/
    /* Block Jacobi, Old Way */
    else if (solver->prec_value == 1) {
        spv_jacobi(matrix, diagonal, solv_b, solv_x, scale_vect, Nsys_solving * mype_Ghost, mype_Ghost, Nsys_solving
#ifdef _MESSG
                   ,smpi
#endif
                   );
        
        /*
         printScreen_dble_array("solv_x", solv_x, mype_Ndof * Nsys_solving, __LINE__, __FILE__);
         printScreen_dble_array("solv_x", scale_vect, mype_Ndof * Nsys_solving, __LINE__, __FILE__);
         printScreen_resid("sw2 matrix",solv_b, mype_Ndof, Nsys_solving, __LINE__, __FILE__);
         printScreen_matrix("sw2 matrix", diagonal, matrix, mype_Ndof, Nsys_solving*Nsys_solving,__LINE__, __FILE__);
         */
        
        if (solver->refresh == YES) {
#ifdef _MESSG
            solv_blk_set(solver, matrix, diagonal, mype_Ndof, Nsys_solving, Nsys_solving * Nsys_solving, smpi->ADH_COMM);
#else
            solv_blk_set(solver, matrix, diagonal, mype_Ndof, Nsys_solving, Nsys_solving * Nsys_solving);
#endif
        }
        
        /*
         printScreen_matrixProfile("after blk set", &(solver->profile), __LINE__, __FILE__);
         exit(-1);
         */
    }
    /************************************************************/
    
    return;
}

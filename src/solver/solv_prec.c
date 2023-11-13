/* performs preconditioning operations on a vector */

/*  Ref:  Smith, B., P. Bjorstad, and W. Gropp, "Domain Decomposition:
   Parallel Multilevel Methods for Elliptic Partial Differential
   Equations", Cambridge University Press, Cambridge, 1996.

 */

#include "global_header.h"

void solv_prec(Solver_Info *solver,
               int *bc_mask,    /* the bc mask for each dof */
               SPARSE_VECT * matrix,    /* the matrix */
               double *diagonal,    /* the diagonal */
               double *v_input, /* the input vector */
               double *v_output,    /* the output vector */
               double *v_tmp,   /* temporary vector */
               int ndof_solv,   /* the number of degrees of freedom */
               int my_ndof_solv,    /* the number of my degrees of freedom */
               int my_nnode_solv,   /* the number of my nodes */
               int nsys_solv    /* the number of equations being solved */
#ifdef _MESSG
               , SMPI *smpi
#endif
    )
{
    int i;                      /* loop counter */

     //printf("prec_value: %d\n",solver->prec_value);

    /* BRANCH 1 */
    if (solver->prec_value == 0) {
        /* no preconditioning  */
        for (i = 0; i < my_ndof_solv; i++)
            v_output[i] = v_input[i];
        return;
        /* END OF BRANCH 1 */
    }

    if (solver->prec_value == 1) {
        /* BRANCH 2 */

        /* one-level additive Schwarz */
        /* The one level additive Schwarz preconditioner
           is just a block Jacobi preconditioner in our
           case.  We are using minimal overlap, so that
           subdomains do not share any of their "interior"
           nodes.  */

        for (i = 0; i < my_ndof_solv; i++) {
            v_output[i] = v_input[i];
        }

        /* calculate the block (subdomain) corrections */
        //printScreen_dble_array("v_output :: before solv_blk_solve", v_output, my_ndof_solv, __LINE__, __FILE__);
#ifdef _MESSG
        solver->UMFail = solv_blk_solve(solver->profile, v_output, my_ndof_solv, smpi->ADH_COMM);
#else
        solver->UMFail = solv_blk_solve(solver->profile, v_output, my_ndof_solv);
#endif
        //printScreen_dble_array("v_output :: after solv_blk_solve", v_output, my_ndof_solv, __LINE__, __FILE__);

        return;

        /* END OF BRANCH 2 */
    }

    if (solver->prec_value == 2) {
        /* BRANCH 3 */

        /* two-level additive Schwarz */
        /* This preconditioner adds a coarse grid correction 
           to the one-level Schwarz fine grid correction.

           Thus the form of the preconditioner M is

           M = B_c + B_f

           where B_c is the coarse grid preconditoner and
           B_f is the fine grid preconditioner (one-level
           additive Schwarz). */

        /* store input vectors for use in fine and coarse
           grid solves */
//        for (i = 0; i < my_ndof_solv; i++) {
//            v_tmp[i] = v_input[i];
//            v_output[i] = v_input[i];
//        }

        /* perform fine and coarse grid solves;
           add results together */
//        solv_coarse_solve(bc_mask, v_tmp, my_nnode_solv, nsys_solv);
//        solv_blk_solve(v_output, my_ndof_solv);
//        for (i = 0; i < my_ndof_solv; i++)
//            v_output[i] += v_tmp[i];
        return;

        /* END OF BRANCH 3 */
    }

    /* BRANCH 4 */

    /* two-level hybrid */
    /* The hybrid preconditioner is actually a mix of a one-level
       additive Schwarz preconditioner with one step of a 
       multiplicative Schwarz method.  The two "domains" are the
       fine domain and the coarse domain.  In this case, the 
       coarse grid preconditioner is applied first, and is followed
       by an application of the fine grid preconditioner in the
       multiplicative Schwarz step.  This is different from Mandel's 
       proposal, where the fine grid preconditioner is applied
       initially, followed by the coarse grid preconditioner in
       the multiplicative step. 

       The form of the preconditioner M is then

       M = B_f + (I - B_f*A_f)B_c

       where B_c is the coarse preconditioner, B_f is the fine
       preconditioner (one-level additive Schwarz), and A_f is 
       the fine grid matrix. */

    /* do the coarse grid solve first;
       store in a temporary value for later use */
//    for (i = 0; i < my_ndof_solv; i++) {
//        v_tmp[i] = v_input[i];
//    }
//    solv_coarse_solve(bc_mask, v_tmp, my_nnode_solv, nsys_solv);

    /* now do a matrix multiply and 
       subtract from identity;
       use this output as input for solv_blk_solve */
//    solv_amult(solver, v_tmp, v_output, diagonal, matrix, my_nnode_solv, nsys_solv);
//    for (i = 0; i < my_ndof_solv; i++)
//        v_output[i] = v_input[i] - v_output[i];

    /* fine grid solve */
//    solv_blk_solve(v_output, my_ndof_solv);

    /* add output from solv_blk_solve to 
       output from coarse solve */
//    for (i = 0; i < my_ndof_solv; i++)
//        v_output[i] += v_tmp[i];

    /* END OF BRANCH 4 */

    /* communicate the output - probably unnecessary */
#ifdef _MESSG
    comm_update_double(v_output, nsys_solv, smpi);
#endif
}

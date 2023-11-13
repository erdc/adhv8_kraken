#define SOLV_INIT
#include "global_header.h"

/* brief Initializes solver information, vectors, etc.
 * Should be called before any other solver routines are called.
 */
void solv_initialize(Solver_Info *solver)
{

    /* Non-Linear Tolerance */
    solver->tol_nonlin = SMALL;

    /* Linear Tolerance */
    solver->tol_lin = SMALL;

    /*
       if (mod->NS_FLOW) {
       solver.tol_lin = solver.tol_nonlin * 0.001;
       if (solver.tol_lin < SOLV_TOL)
       solver.tol_lin = SOLV_TOL;
       }
       if (mod->SW2_FLOW || mod->SW3_FLOW) {
       solver.tol_lin = solver.tol_nonlin * 0.0001;
       if (solver.tol_lin < SOLV_TOL)
       solver.tol_lin = SOLV_TOL;
       }
       if (mod->GW_FLOW) {
       solver.tol_lin = solver.tol_nonlin * 0.001;
       if (solver.tol_lin < SOLV_TOL)
       solver.tol_lin = SOLV_TOL;
       }
     */

    //printf("solver tol :: lin: %30.20e  nonlin: %30.20e\n",solver->tol_lin,solver->tol_nonlin);
    solver->tol_lin = solver->tol_nonlin * 0.0001;
    if(solver->tol_lin < SOLV_TOL) solver->tol_lin = SOLV_TOL;


    /* Preconditioner Choice */
    solver->prec_value = 21;
    /* Total Count of Linear Iterations */
    solver->it_count_lin = 0;
    /* Total Count of Non-Linear Iterations */
    solver->it_count_nonlin = 0;
    /* Total Count of Linear Iterations "Lost" to Failure */
    solver->it_count_lin_failed = 0;
    /* Total Count of Non-Linear Iterations "Lost" to Failure */
    solver->it_count_nonlin_failed = 0;
    /* Is Linear System currently in Failed State? */
    solver->lin_fail = NO;
    /* Is Non-Linear System currently in Failed State? */
    solver->nonlin_fail = NO;
    /* Maximum Number of Linear Iterations to Use */
    solver->max_lin_it = 500;
    /* Maximum Number of Non-Linear Iterations to Use */
    solver->max_nonlin_it = 20;
    /* Block Size to use for Preconditioner of Linear System */
    solver->blocksize = 8;
    /* Refresh the Preconditioner? */
    solver->refresh = YES;
    /* Force the acceptance of the Linear Iteration Step? */
    solver->force_lin_it = NO;
    /* Force the acceptance of the Non-Linear Iteration Step? */
    solver->force_nonlin_it = NO;
    /* Maximum Number of Line Search Cuts For a Non-Linear Step */
    solver->max_nonlin_linesearch_cuts = 5;
    /* Flag for if linear problem or not */
    solver->LINEAR_PROBLEM = NO;
    /* Which physics are we solving */
    solver->PRN_NEWTON = OFF;

    solver->profile.size = 0;      /* the size of the profile matrix */
    solver->profile.rows = NULL;
    solver->profile.cols = NULL;  
    solver->node_block = NULL;  

    return;
}

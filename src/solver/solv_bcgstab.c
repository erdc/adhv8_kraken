/* This is a biconjugate gradient squared stabilized matrix solver. 
 It is a Krylov-based iterative method used to solve nonsymmetric
 linear systems.  
 
 Ref:  H. A. vanderVorst, "Bi-CGStab:  A fast and smoothly converging 
 variant to Bi-CG for the solution of nonsymmetric systems",
 SIAM J. Sci. Statist. Comput., 13 (1992), pp. 631-644.
 
 or       C. T. Kelley, "Iterative Methods for Linear and Nonlinear
 Equations", SIAM, Philadelphia, PA, 1995.
 
 In general, the Bi-CGStab solver constructs bases for the Krylov 
 subspaces (A,b) and (A^T,c), where b is the right-hand-side vector
 and c is usually chosen so that c=b.  The method may break down, 
 although this does not seem to happen often in practice.
 
 Return: Linear Solver Failure?
 NO or YES
 */

#include "global_header.h"

static int isize = 0;           /* the size of the arrays */
static double *solv_r;          /* the linear solver residual */
static double *solv_p;          /* the search direction */
static double *solv_Ap;         /* AMp */
static double *solv_Mp;         /* Mp */
static double *solv_q;          /* the shadow residual */
static double *solv_s;          /* the second search direction for bcgstab */
static double *solv_As;         /* AMs */
static double *solv_Ms;         /* Ms */
static double *solv_Au;         /* Au */
static double *solv_u0;         /* solution = u+u0 - u0 is the shift */
static double *solv_ptmp;       /* temporary vector for the preconditioner */

/* Declaration */
void solv_print_my_matrix(SPARSE_VECT * matrix, double *, double *, int);

int solv_bcgstab(Solver_Info *solver,
                 int *bc_mask,  /* boundary condition mask */
                 SPARSE_VECT * matrix,  /* the matrix */
                 double *diagonal,  /* the diagonal */
                 double *solv_b,    /* the right hand side */
                 double *solv_u,    /* the solution */
                 double *scale_vect,    /* the scale vector for the jacobi preconditioning */
                 int nnode_solv,    /* the total number of nodes */
                 int my_nnode_solv, /* the number of my nodes */
                 int nsys_solv  /* the number of equations being solved */
#ifdef _MESSG
                 , SMPI *smpi
#endif
)
{
    int it;                     /* loop counter over the cg iterations */
    //int iapprox_update_flag;    /* flag for performing an update of the approximation */
    //int iresid_update_flag;     /* flag for performing an update of the residual */
    int isize_prev;             /* the previous size */
    int ndof_solv;              /* the total number of degrees of freedom */
    int my_ndof_solv;           /* the number of my degrees of freedom */
    double rnorm;               /* the norm of the residual */
    double alpha = 1.0;         /* scalar */
    double omega = 1.0;         /* scalar */
    double beta = 0.0;          /* scalar */
    double gamma = 0.0;         /* p dot ap */
    double rho = 1.0;           /* scalar */
    double rhop;                /* rho from the previous iteration */
    double bnorm;               /* the norm of the right hand side */
    //double snorm;               /* the norm of s */
    double conv_tol;            /* the convergence tolerance */
    double min_conv_tol;        /* the minimum convergence tolerance */
    
    /* sets the dof variables */
    ndof_solv = nnode_solv * nsys_solv;
    my_ndof_solv = my_nnode_solv * nsys_solv;
    
    /* allocates memory if needed */
    if (isize < ndof_solv) {
        isize_prev = isize;
        isize = ndof_solv;
        solv_r = (double *) tl_realloc(sizeof(double), isize, isize_prev, solv_r);
        solv_p = (double *) tl_realloc(sizeof(double), isize, isize_prev, solv_p);
        solv_Ap = (double *) tl_realloc(sizeof(double), isize, isize_prev, solv_Ap);
        solv_Mp = (double *) tl_realloc(sizeof(double), isize, isize_prev, solv_Mp);
        solv_q = (double *) tl_realloc(sizeof(double), isize, isize_prev, solv_q);
        solv_s = (double *) tl_realloc(sizeof(double), isize, isize_prev, solv_s);
        solv_As = (double *) tl_realloc(sizeof(double), isize, isize_prev, solv_As);
        solv_Ms = (double *) tl_realloc(sizeof(double), isize, isize_prev, solv_Ms);
        solv_Au = (double *) tl_realloc(sizeof(double), isize, isize_prev, solv_Au);
        solv_u0 = (double *) tl_realloc(sizeof(double), isize, isize_prev, solv_u0);
        solv_ptmp = (double *) tl_realloc(sizeof(double), isize, isize_prev, solv_ptmp);
    }
    
    /* zeroes the arrays */
    solv_copy(ndof_solv, solv_u, solv_u0);
    solv_init_dbl(ndof_solv, solv_u);
    solv_init_dbl(ndof_solv, solv_r);
    solv_init_dbl(ndof_solv, solv_q);
    solv_init_dbl(ndof_solv, solv_p);
    solv_init_dbl(ndof_solv, solv_Mp);
    solv_init_dbl(ndof_solv, solv_Ap);
    solv_init_dbl(ndof_solv, solv_s);
    solv_init_dbl(ndof_solv, solv_As);
    solv_init_dbl(ndof_solv, solv_Ms);
    
    /* Pre-Condition the RHS */
    /*
     printf("CUT TO HERE\n");
     printf("ndof_solv: %d \t my_ndof_solv: %d \t my_nnode_solv: %d\n",ndof_solv,my_ndof_solv,my_nnode_solv);
     printScreen_int_array("bc_mask", bc_mask, ndof_solv, __LINE__, __FILE__);
     printScreen_dble_array("scale_vect", scale_vect, ndof_solv, __LINE__, __FILE__);
     printScreen_dble_array("solv_p", solv_p, ndof_solv, __LINE__, __FILE__);
     printScreen_dble_array("solv_b", solv_b, ndof_solv, __LINE__, __FILE__);
     printScreen_dble_array("diagonal", diagonal, my_nnode_solv * nsys_solv * nsys_solv, __LINE__, __FILE__);
     printScreen_matrix("sw3 matrix", diagonal, matrix, my_nnode_solv, nsys_solv * nsys_solv, __LINE__, __FILE__);
     */
    
    solv_copy(ndof_solv, solv_b, solv_p);
    solv_init_dbl(ndof_solv, solv_b);

#ifdef _MESSG
    solv_prec(solver, bc_mask, matrix, diagonal, solv_p, solv_b, solv_ptmp, ndof_solv, my_ndof_solv, my_nnode_solv, nsys_solv, smpi);
#else
    solv_prec(solver, bc_mask, matrix, diagonal, solv_p, solv_b, solv_ptmp, ndof_solv, my_ndof_solv, my_nnode_solv, nsys_solv);
#endif
    solv_init_dbl(ndof_solv, solv_p);
    
    /*
     printScreen_dble_array("solv_p :: after bc_stab", solv_p, ndof_solv, __LINE__, __FILE__);
     printScreen_dble_array("solv_b :: after bc_stab", solv_b, ndof_solv, __LINE__, __FILE__);
     exit(-1);
     */
    
    /* Now Form Actual Residual */
    /* Assuming x=0 as the initial guess would have given r=b (the rhs). Easy, right? */
    /* However, we miss an opportunity if our preconditioner is meaningful. */
    /* This will basically take r = S b - S A x, where x = S b, as the initial guess. */
    solv_copy(ndof_solv, solv_b, solv_u);
    solv_amult(solver, solv_u, solv_Mp, diagonal, matrix, my_nnode_solv, nsys_solv
#ifdef _MESSG
               , smpi
#endif 
               );
    solv_prec(solver, bc_mask, matrix, diagonal, solv_Mp, solv_Ap, solv_ptmp, ndof_solv, my_ndof_solv, my_nnode_solv, nsys_solv
#ifdef _MESSG
              , smpi
#endif
              );
    solv_copy(ndof_solv, solv_b, solv_r);
    solv_daxpy(ndof_solv, NONE, solv_Ap, solv_r);
    solv_init_dbl(ndof_solv, solv_Mp);
    solv_init_dbl(ndof_solv, solv_Ap);
    solv_copy(ndof_solv, solv_r, solv_q);
#ifdef _MESSG
    double rnorm_max,bnorm_max;
    rnorm = solv_infty_norm(my_ndof_solv, solv_r, smpi->ADH_COMM);
    bnorm = solv_infty_norm(my_ndof_solv, solv_b, smpi->ADH_COMM);
    rnorm_max = messg_dmax(rnorm,smpi->ADH_COMM);
    bnorm_max = messg_dmax(bnorm,smpi->ADH_COMM);
    bnorm=bnorm_max;
    rnorm=rnorm_max;
#else 
    rnorm = solv_infty_norm(my_ndof_solv, solv_r);
    bnorm = solv_infty_norm(my_ndof_solv, solv_b);
#endif

    
    rhop = 1.0;
    alpha = 1.0;
    omega = 1.0;
    it = 0;
    if (solver->LINEAR_PROBLEM == NO)
        min_conv_tol = 1.E-10;
    if (solver->LINEAR_PROBLEM == YES)
        min_conv_tol = 1.E-10;
    if (solver->LINEAR_PROBLEM == NO)
        conv_tol = 1.E-5;
    if (solver->LINEAR_PROBLEM == YES)
        conv_tol = 1.E-5;
    rho = solv_dot(my_ndof_solv, solv_r, solv_q);
#ifdef _MESSG
    rho = messg_dsum(rho , smpi->ADH_COMM);
#endif
    
    //printf("rnorm: %20.10f \n",rnorm);
    //printf("bnorm: %20.10f \n",bnorm);
    //printf("rnorm: %30.20f \t\t tol: %30.20f \t\t it: %d \t\t max: %d\n",rnorm, min_conv_tol + bnorm * conv_tol, it,solver->max_lin_it);
    
    while (rnorm > (min_conv_tol + bnorm * conv_tol) && it < solver->max_lin_it) {
        /* (a) */
        it++;
        
        /* (b) */
        beta = (rho * alpha) / (rhop * omega);
        
        /* (c) */
        solv_daxpy(ndof_solv, -omega, solv_Ap, solv_p);
        solv_scal(ndof_solv, beta, solv_p);
        solv_daxpy(ndof_solv, DONE, solv_r, solv_p);
        
        /* (d) */
        solv_init_dbl(ndof_solv, solv_Mp);
        solv_init_dbl(ndof_solv, solv_Ap);
        solv_amult(solver, solv_p, solv_Mp, diagonal, matrix, my_nnode_solv, nsys_solv
#ifdef _MESSG
                   , smpi
#endif
                   );
        solv_prec(solver, bc_mask, matrix, diagonal, solv_Mp, solv_Ap, solv_ptmp, ndof_solv, my_ndof_solv, my_nnode_solv, nsys_solv
#ifdef _MESSG
                  , smpi
#endif
                  );
        
        /* (e) */
        gamma = solv_dot(my_ndof_solv, solv_q, solv_Ap);
#ifdef _MESSG
        gamma = messg_dsum(gamma, smpi->ADH_COMM);
#endif
        alpha = rho / gamma;
        
        /* (f) */
        solv_copy(ndof_solv, solv_r, solv_s);
        solv_daxpy(ndof_solv, -alpha, solv_Ap, solv_s);
        solv_init_dbl(ndof_solv, solv_Ms);
        solv_init_dbl(ndof_solv, solv_As);
        solv_amult(solver, solv_s, solv_Ms, diagonal, matrix, my_nnode_solv, nsys_solv
#ifdef _MESSG
                   , smpi
#endif
                   );
        solv_prec(solver, bc_mask, matrix, diagonal, solv_Ms, solv_As, solv_ptmp, ndof_solv, my_ndof_solv, my_nnode_solv, nsys_solv
#ifdef _MESSG
                  , smpi
#endif
                  );
        
        /* (g) */
        omega = solv_dot(my_ndof_solv, solv_As, solv_s);
        gamma = solv_dot(my_ndof_solv, solv_As, solv_As);
#ifdef _MESSG
        omega = messg_dsum(omega, smpi->ADH_COMM);
        gamma = messg_dsum(gamma, smpi->ADH_COMM);
#endif
        omega = omega / gamma;
        rhop = rho;
        rho = -1.0 * omega * solv_dot(my_ndof_solv, solv_As, solv_q);
#ifdef _MESSG
        rho = messg_dsum(rho, smpi->ADH_COMM);
#endif
        
        /* (h) */
        solv_daxpy(ndof_solv, alpha, solv_p, solv_u);
        solv_daxpy(ndof_solv, omega, solv_s, solv_u);
        
        /* (i) */
        solv_copy(ndof_solv, solv_s, solv_r);
        solv_daxpy(ndof_solv, -omega, solv_As, solv_r);
        /*rnorm = solv_l2_norm_scaled(my_ndof_solv, solv_r, (double)(global_nnode)); */
#ifdef _MESSG
        rnorm = solv_infty_norm(my_ndof_solv, solv_r, smpi->ADH_COMM);
        rnorm_max = messg_dmax(rnorm,smpi->ADH_COMM);
        rnorm = rnorm_max;
#else
       rnorm = solv_infty_norm(my_ndof_solv, solv_r);
#endif
        
#ifdef _DEBUG_KITCHEN_SINK
        //if (myid == 0) {
        //    printf("Resid[%d] = %16.15e\n", it, rnorm);
        //}
#endif
    }
#ifdef _MESSG 
    it = messg_imax(it, smpi->ADH_COMM);
    if (solver->PRN_NEWTON && smpi->myid <= 0)
#else
        if (solver->PRN_NEWTON)
#endif
            printf("MIT: %2d", it+1);
    if (it == solver->max_lin_it) {
        if (solver->force_lin_it == NO) {
            solv_return_old(solver, solv_u, solv_u0, scale_vect, my_ndof_solv, nsys_solv
#ifdef _MESSG
                            , smpi
#endif
                            );
            return (NO);
        }
        else {
            return (YES);
        }
    }
    /* the linear solver succeeded */
    
    solv_return_old(solver, solv_u, solv_u0, scale_vect, my_ndof_solv, nsys_solv
#ifdef _MESSG
                    , smpi
#endif
                    );
    return (YES);
}

void solv_bcgstab_clean(void)
{
    if (solv_r != NULL)
        solv_r = (double *) tl_free(sizeof(double), isize, solv_r);
    if (solv_p != NULL)
        solv_p = (double *) tl_free(sizeof(double), isize, solv_p);
    if (solv_Ap != NULL)
        solv_Ap = (double *) tl_free(sizeof(double), isize, solv_Ap);
    if (solv_Mp != NULL)
        solv_Mp = (double *) tl_free(sizeof(double), isize, solv_Mp);
    if (solv_q != NULL)
        solv_q = (double *) tl_free(sizeof(double), isize, solv_q);
    if (solv_s != NULL)
        solv_s = (double *) tl_free(sizeof(double), isize, solv_s);
    if (solv_As != NULL)
        solv_As = (double *) tl_free(sizeof(double), isize, solv_As);
    if (solv_Ms != NULL)
        solv_Ms = (double *) tl_free(sizeof(double), isize, solv_Ms);
    if (solv_Au != NULL)
        solv_Au = (double *) tl_free(sizeof(double), isize, solv_Au);
    if (solv_u0 != NULL)
        solv_u0 = (double *) tl_free(sizeof(double), isize, solv_u0);
    if (solv_ptmp != NULL)
        solv_ptmp = (double *) tl_free(sizeof(double), isize, solv_ptmp);
}

/*!
 \brief Print out the ENTIRE matrix for Debugging
 */
void solv_print_my_matrix(SPARSE_VECT * matrix, /* the matrix */
                          double *diagonal, /* the diagonal */
                          double *solv_b,   /* the rhs */
                          int mype_Ndof /* the number of my nodes */
)
{
    int ii, jj, kk;
    for (ii = 0; ii < mype_Ndof; ii++)
        printf("RHS[%5d]=%16.15e\n", ii, solv_b[ii]);
    for (ii = 0; ii < mype_Ndof; ii++) {
        printf("A[%5d][%5d]=%16.15e\n", ii, ii, diagonal[ii]);
        for (jj = 0; jj < matrix[ii].size; jj++) {
            kk = matrix[ii].index[jj];
            printf("A[%5d][%5d]=%16.15e\n", ii, kk, matrix[ii].value[jj]);
        }
    }
}

/*! \file  bcgstab_solver.c This file has functions responsible for solving linear system in split CSR format using BCG-stabilized method */
#include "adh.h"
//these could be macros or user-defined at run time
#define SOLV_TOL 1e-5
//keep umfpack calls in one script
static  void *Symbolic, *Numeric;
static int isize = 0;           /* the size of the arrays */
static double *r;          /* the linear solver residual */
static double *p;          /* the search direction */
static double *Ap;         /* AMp */
static double *Mp;         /* Mp */
static double *q;          /* the shadow residual */
static double *s;          /* the second search direction for bcgstab */
static double *As;         /* AMs */
static double *Ms;         /* Ms */
static double *x0;         /* solution = u+u0 - u0 is the shift */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Solves a system of already scaled equations in split CSR format using bcgstab
 *  this is also assuming there is some scaling prior to solve in scale_vect. If no scaling,
 *  then scale_vect needs to be all 1's. Preconditioning of block diagonal (local to process) via umfpack
 *  is employed
 * 
 *  Ref:  H. A. vanderVorst, "Bi-CGStab:  A fast and smoothly converging 
 *  variant to Bi-CG for the solution of nonsymmetric systems",
 *  SIAM J. Sci. Statist. Comput., 13 (1992), pp. 631-644.
 *
 *  or       C. T. Kelley, "Iterative Methods for Linear and Nonlinear
 *  Equations", SIAM, Philadelphia, PA, 1995.
 *
 *  In general, the Bi-CGStab solver constructs bases for the Krylov 
 *  subspaces (A,b) and (A^T,c), where b is the right-hand-side vector
 *  and c is usually chosen so that c=b.  The method may break down, 
 *  although this does not seem to happen often in practice.
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] x (double*) - solution to system Ax=b
 *  @param[in] indptr_diag (int*) - indeces of first/last index for each row in diagonal block (local to process) in split CSR matrix
 *  @param[in] cols_diag (int*) - column numbers (local to process) of diagonal block in split CSR matrix
 *  @param[in] vals_diag (double*) - array of doubles that are the nonzero entries of diagonal block (local to process) in split CSR matrix
 *  @param[in] indptr_off_diag (int*) - indeces of first/last index for each row in off-diagonal block (local to process) in split CSR matrix
 *  @param[in] cols_off_diag (int*) - column numbers (global) of off-diagonal block in split CSR matrix
 *  @param[in] vals_off_diag (double*) - array of doubles that are the nonzero entries of off-diagonal block (local to process) in split CSR matrix
 *  @param[in] b (double*) - scaled rhs of the linear system
 *  @param[in] scale_vect (double*) - the vector that scaled the system of equations prior to solve, used for rescaling after solve is finished
 *  @param[in] local_size (int) - number of equations (rows in the matrix) owned by the process
 *  @param[in] size_with_ghosts (int) - number of rows + number of ghost d.o.f.s present
 *  @param[in] rank (int) - MPI rank
 *  @param[in] ghosts (int*) - array of global dof numbers that are ghosts on the process
 *  @param[in] nghost (int) - number of ghost d.o.f.s on the process
 * \returns integer that is a code if the solve worked or not
 * \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int solve_linear_sys_bcgstab(double *x, int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag,double *vals_off_diag, double *b, 
  double *scale_vect, int local_size, int size_with_ghosts, int rank,
  int *ghosts, int nghost){
  int it;                     /* loop counter over the cg iterations */
  //int iapprox_update_flag;    /* flag for performing an update of the approximation */
  //int iresid_update_flag;     /* flag for performing an update of the residual */
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
  int isize_prev;
  /* allocates memory if needed */
  if (isize < size_with_ghosts) {
        isize_prev = isize;
        isize = size_with_ghosts;
        r = (double *) tl_realloc(sizeof(double), isize, isize_prev, r);
        p = (double *) tl_realloc(sizeof(double), isize, isize_prev, p);
        Ap = (double *) tl_realloc(sizeof(double), isize, isize_prev, Ap);
        Mp = (double *) tl_realloc(sizeof(double), isize, isize_prev, Mp);
        q = (double *) tl_realloc(sizeof(double), isize, isize_prev, q);
        s = (double *) tl_realloc(sizeof(double), isize, isize_prev, s);
        As = (double *) tl_realloc(sizeof(double), isize, isize_prev, As);
        Ms = (double *) tl_realloc(sizeof(double), isize, isize_prev, Ms);
        x0 = (double *) tl_realloc(sizeof(double), isize, isize_prev, x0);
  }
  printf("Ndof with ghosts = %d\n",size_with_ghosts);
  /* allocates memory if needed */
  // need to do different allocation later
  //breaks down for larger problems
  //double r[size_with_ghosts];
  //printf("allocated one doubles\n");
  //double p[size_with_ghosts];
  //printf("allocated two doubles\n");
  //double Ap[size_with_ghosts];
  //double Mp[size_with_ghosts];
  //double q[size_with_ghosts];
  //double s[size_with_ghosts];
  //double As[size_with_ghosts];
  //double Ms[size_with_ghosts];
  //double x0[size_with_ghosts];


  /* zeroes the arrays */
  //need to store this to shift solution after solve
  sarray_copy_dbl(x0, x, size_with_ghosts);

  //initialize arrays to 0
  sarray_init_dbl(x, size_with_ghosts);
  sarray_init_dbl(r, size_with_ghosts);
  sarray_init_dbl(q, size_with_ghosts);
  sarray_init_dbl(p, size_with_ghosts);
  sarray_init_dbl(Mp, size_with_ghosts);
  sarray_init_dbl(Ap, size_with_ghosts);
  sarray_init_dbl(s, size_with_ghosts);
  sarray_init_dbl(As, size_with_ghosts);
  sarray_init_dbl(Ms, size_with_ghosts);

  sarray_copy_dbl(p, b, size_with_ghosts);
  //in original routine but doesnt seem necessary
  sarray_init_dbl(b,size_with_ghosts);
  int i;
  //apply left preconditioning to RHS, this is stored in residual
  solve_umfpack(b, indptr_diag, cols_diag, vals_diag, p, local_size);
  //zero out rhs
  sarray_init_dbl(p, size_with_ghosts);
  //update ghosts

  /* Now Form Actual Residual */
  /* Assuming x=0 as the initial guess would have given r=b (the rhs). Easy, right? */
  /* However, we miss an opportunity if our preconditioner is meaningful. */
  /* This will basically take r = S b - S A x, where x = S b, as the initial guess. */
  //should it really be updating the ghosts?
  sarray_copy_dbl(x, b, size_with_ghosts);
  //for (i=0;i<local_size;i++){
  //  printf("Rank %d, x[%d] = %f\n",rank,i,vals_off_diag[i]);
  //}
  comm_update_double(x, local_size, 3, rank);
  split_CSR_mat_vec_mult(Mp, indptr_diag, cols_diag, vals_diag, indptr_off_diag, cols_off_diag, vals_off_diag,
   x, local_size,ghosts,nghost);
  //update matrix vector indices. is this necessary?
  comm_update_double(Mp, local_size, 3, rank);
  //for (i=0;i<local_size;i++){
  //  printf("Rank %d, Ax[%d] = %f\n",rank,i,Mp[i]);
  //}
  //next, apply preconditioner to Mp, save in Ap
  solve_umfpack(Ap, indptr_diag, cols_diag, vals_diag, Mp, local_size);
  //print vector
  //for (i=0;i<local_size;i++){
  //  printf("Rank %d, Ax[%d] = %f\n",rank,i,Ap[i]);
  //}
  //from b to r and set as r = b-Ax with prexondition
  sarray_copy_dbl(r, b, size_with_ghosts);
  sarray_y_plus_ax_dbl(r, -1, Ap, size_with_ghosts );
  sarray_init_dbl(Mp, size_with_ghosts);
  sarray_init_dbl(Ap, size_with_ghosts);
  sarray_copy_dbl(q, r, size_with_ghosts);
  //print vector

  //take linf norms of vectors
  rnorm = sarray_l_infty_norm(r, local_size);
  bnorm = sarray_l_infty_norm(b, local_size);
  //if it is in parallel. collect with allreduce
  rnorm = messg_dmax(rnorm);
  bnorm = messg_dmax(bnorm);

  printf("RNORM = %6.4e\n",rnorm);
  printf("BNORM = %6.4e\n",bnorm);
  rhop = 1.0;
  alpha = 1.0;
  omega = 1.0;
  it = 0;
  min_conv_tol = 1e-11;
  conv_tol = 1e-11;
  rho = sarray_dot_dbl(r, q, local_size);
  //if in parallel sum among processors
  rho = messg_dsum(rho);
  //printf("RHO = %f\n",rho);
  int max_it = 100;
  //start itertative loop
  while (rnorm > (min_conv_tol + bnorm * conv_tol) && it < max_it){
    // a
    it++;

    //(b)
    beta = (rho * alpha) / (rhop * omega);


    // c //
    sarray_y_plus_ax_dbl(p, -omega, Ap,size_with_ghosts );
    sarray_scale_replace_dbl(p,beta, size_with_ghosts);
    sarray_y_plus_ax_dbl(p, 1, r, size_with_ghosts);
    //for (i=0;i<local_size;i++){
    //printf("Rank %d, p[%d] = %f\n",rank,i,p[i]);
    //}
    //need comm_update_double befor matvec mult
    comm_update_double(p, local_size, 3, rank);

    // (d) //
    sarray_init_dbl(Mp, size_with_ghosts);
    sarray_init_dbl(Ap,size_with_ghosts);
    split_CSR_mat_vec_mult(Mp, indptr_diag, cols_diag, vals_diag, indptr_off_diag, cols_off_diag, vals_off_diag,
    p, local_size,ghosts,nghost);
    //if MPI
    //preconfitioner is local only
    //comm_update_double(Mp, local_size, 3, rank);
    //next, apply preconditioner to Mp, save in Ap
    solve_umfpack(Ap, indptr_diag, cols_diag, vals_diag, Mp, local_size);

    // (e) //
    gamma = sarray_dot_dbl(q, Ap, local_size);
    gamma = messg_dsum(gamma);
    alpha = rho / gamma;
    //printf("gamma, alpha %f, %f\n",gamma,alpha);

    /* (f) */
    sarray_copy_dbl(s, r, size_with_ghosts);
    sarray_y_plus_ax_dbl(s, -alpha, Ap, size_with_ghosts);
    sarray_init_dbl(Ms, size_with_ghosts);
    sarray_init_dbl(As, size_with_ghosts);
    //need comm_update_double befor matvec mult
    comm_update_double(s, local_size, 3, rank);
    split_CSR_mat_vec_mult(Ms, indptr_diag, cols_diag, vals_diag, indptr_off_diag, cols_off_diag, vals_off_diag,
    s, local_size,ghosts,nghost);
    solve_umfpack(As, indptr_diag, cols_diag, vals_diag, Ms, local_size);

    /* (g) */
    omega = sarray_dot_dbl(As, s, local_size);
    gamma = sarray_dot_dbl(As, As, local_size);
    //if MPI is on
    omega = messg_dsum(omega);
    gamma = messg_dsum(gamma);

    omega = omega / gamma;
    rhop = rho;
    rho = -1.0 * omega * sarray_dot_dbl(As, q, local_size);
    //only MPI
    rho = messg_dsum(rho);

    /* (h) */
    sarray_y_plus_ax_dbl(x, alpha, p, size_with_ghosts);
    sarray_y_plus_ax_dbl(x, omega, s, size_with_ghosts);
        
    /* (i) */
    sarray_copy_dbl(r, s, size_with_ghosts);
    sarray_y_plus_ax_dbl(r, -omega, As, size_with_ghosts);
    /*rnorm = solv_l2_norm_scaled(my_ndof_solv, r, (double)(global_nnode)); */
    rnorm = sarray_l_infty_norm(r, local_size);
    //only do this in parallel
    rnorm = messg_dmax(rnorm);
    //printf("Rnorm, %f\n",rnorm);

  }
  printf("BCGSTAB # IT = %d\n",it);

  //undo scaling (need to add in error checks)
  unscale_linear_system(x,x0,scale_vect,local_size);

//  for (i = 0; i < local_size; i++) {
//        x[i] += x0[i];
//        /* solv_u[i] *= scale_vect[i]; UNDO SCALING */
//    }
//  for (i = 0; i <local_size; i++){
//        x[i] *= scale_vect[i]; /* UNDO SCALING */
//    }
  comm_update_double(x, local_size, 3, rank);
 // for (i = 0; i <local_size; i++){
//    printf("Rank %d Final u[%d] = %f\n",rank,i,x[i]);
//  }

  //need to free after each solve? or not till end of simulation?
  umfpack_di_free_symbolic (&Symbolic) ;
  umfpack_di_free_numeric (&Numeric) ;
  //add later to see if converged or not
  int status = 0;
  return status;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Short helper function that clears umfpack structures
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  \note Clears structure that factors matrix, should this be done after each solve or only
 *  after refinement? 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void umfpack_clear(void){
  //clear memory
  umfpack_di_free_symbolic (&Symbolic) ;
  umfpack_di_free_numeric (&Numeric) ;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Function that pre-factors matrix, used for preconditioner, modifies the static Symbolic and Numeric structures
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in] indptr_diag (int*) - array of integers with first/last indeces of each row in CSR matrix
 *  @param[in] cols_diag (int*) - array of integers that are local (to process) column numbers
 *  @param[in] vals_diag (double*) - the values of the sparse matrix
 *  @param[in] nrow (int) - number of equations in system
 *  /returns (integer) status of umfpack routine
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int prep_umfpack(int *indptr_diag, int *cols_diag, double *vals_diag, int nrow){
  //load into umfpack, need to transpose because umfpack uses ccs format
  int status;
  double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];
  umfpack_di_defaults(Control);    /* default control parameters */
  status = umfpack_di_symbolic (nrow, nrow, indptr_diag, cols_diag, vals_diag, &Symbolic, Control, Info);
  status = umfpack_di_numeric (indptr_diag, cols_diag, vals_diag, Symbolic, &Numeric, Control, Info) ;
  return status;

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Function that solves system Ax=b using umfpack, provided A is in CSR format
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] x (double*) - solution array to system Ax=b
 *  @param[in] indptr_diag (int*) - array of integers with first/last indeces of each row in CSR matrix
 *  @param[in] cols_diag (int*) - array of integers that are local (to process) column numbers
 *  @param[in] vals_diag (double*) - the values of the sparse matrix
 *  @param[in] b (double*) - rhs of system Ax=b
 *  @param[in] nrow (int) - number of equations in system
 *  /returns (integer) status of umfpack routine
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int solve_umfpack(double *x, int *indptr_diag, int *cols_diag, double *vals_diag, double *b, int nrow){
  int status;
  double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];
  umfpack_di_defaults(Control);    /* default control parameters */
  //try this to account for CSR format, use solution as RHS
  status = umfpack_di_solve (UMFPACK_Aat, indptr_diag, cols_diag, vals_diag, x, b, Numeric, Control, Info);
  return status;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees all data structures allocate within this file
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void free_bcgstab(void){
  if (r != NULL)
        r = (double *) tl_free(sizeof(double), isize, r);
    if (p != NULL)
        p = (double *) tl_free(sizeof(double), isize, p);
    if (Ap != NULL)
        Ap = (double *) tl_free(sizeof(double), isize, Ap);
    if (Mp != NULL)
        Mp = (double *) tl_free(sizeof(double), isize, Mp);
    if (q != NULL)
        q = (double *) tl_free(sizeof(double), isize, q);
    if (s != NULL)
        s = (double *) tl_free(sizeof(double), isize, s);
    if (As != NULL)
        As = (double *) tl_free(sizeof(double), isize, As);
    if (Ms != NULL)
        Ms = (double *) tl_free(sizeof(double), isize, Ms);
    if (x0 != NULL)
        x0 = (double *) tl_free(sizeof(double), isize, x0);
}


//double get_global_max(double x){
//#ifdef _MESSG
//  double x_send = 0.0;        /* the variables for the pass */
//  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */
//  x_send = x;
//  ierr_code = MPI_Allreduce(&x_send, &x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//#endif
//  return x;
//}

//void scale_replace_dbl(double *v,  double alpha, int n){
//  int i;
//  //use blas maybe later
//  /* multiplies by a scalar */
//  for(i = 0; i < n; i++)
//        v[i] *= alpha;
//}

////\brief Performs daxpy operation for (double) vectors, y = \alpha x  + y$
//void y_plus_ax( double *y,     /* y vector (input and output) */
//                double alpha,     /* the scalar */
//                double *x,      /* x vector */
//                int n      /* the number of degrees of freedom */
//                
//)
//{
////can use blas or similar for this
//int i = 0;      /* loop counter */
//    
//    /* adds the arrays */
//    for(i = 0; i < n; i++)
//        y[i] += alpha * x[i];
//    return;
//}

/*!
 \brief Returns the \f$l_\infty\f (infinity) norm of a vector (maximum absolute value)
 
 \param n the length of the vectors (int)
 \param v1 pointer to the vector (double, length n)
 */
//double l_infty_norm(
//                       int n,     /* the number of degrees of freedom */
//                       double *v1     /* the vector */
//                       )
//{
//    double value = 0.0;   /* the partial sum */
//    int ii = 0;     /* loop counter */
//    
//    /* computes the value for this processor */
//    for(ii = 0; ii < n; ii++)
//    {
//        value = max_dbl(value, fabs(v1[ii]));
//    }
//    
//    /* returns the maximum */
//    return value;
//}

/*!
 \brief Returns the inner product, \f$(x \cdot y)\f$, of two (double) vectors
 
 Sums among all processors.
 
 \param n the length of the vectors (int)
 \param x pointer to the first vector (double, length \a n)
 \param y pointer to the second vector (double, length \a n)
 */
//double dot_dbl_array(
//                int n,      /* the number of degrees of freedom */
//                double *x,      /* the first vector */
//                double *y     /* the second vector */
//)
//{
//    double value = 0.0;   /* the partial sum */
//    int i = 0;      /* loop counter */
//    
//    for(i = 0; i < n; i++)
//        value += x[i] * y[i];//

//    /* returns the sum */
//    return (value);
//}






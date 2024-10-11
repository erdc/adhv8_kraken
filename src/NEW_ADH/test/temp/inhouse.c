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

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <umfpack.h>
 #define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define SOLV_TOL 1e-5
#define MAX_NODAL_DOF 4
void scale_linear_system(int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag,double *vals_off_diag, double *residual, 
  double *sol, double *scale_vect, int local_size, int size_with_ghosts, int rank,int *ghosts, int nghost);
int solve_linear_sys_bcgstab(double *x, int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag,double *vals_off_diag, double *b,
   double *scale_vect, int local_size, int size_with_ghosts, int rank,int *ghosts, int nghost);
double find_max_in_row(int *indptr, double *vals,int row_num);
void comm_update_double(double *vec,int size_v, int npe,int rank);
int global_to_local(int global_col_no,int local_size,int *ghosts, int nghost);
int prep_umfpack(int *indptr_diag, int *cols_diag, double *vals_diag, double *sol, double *resid, int nrow);
int solve_umfpack(double *x, int *indptr_diag, int *cols_diag, double *vals_diag, double *b, int nrow);
void umfpack_clear(void);
void split_CSR_mat_vec_mult(double *Ax, int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag, double *vals_off_diag,
  double *x, int nrows, int *ghosts, int nghost);
void copy_array_to_from(double *vec_to, double *vec_from, int n);
void zero_dbl_array(double *v, int n);
void y_plus_ax(double *y, double alpha, double *x,int n );
double get_global_max(double x);
double l_infty_norm(int n, double *v1);
double dot_dbl_array(int n, double *x,double *y);
double messg_dsum(double x);
void scale_dbl_array(double *v,  double alpha, int n);
double max_dbl(double a, double b);

//keep umfpack calls in one script
static  void *Symbolic, *Numeric;

int main(int argc, char **argv) {

  int status;
  int nghost = 0;
  int m=0;
  int n=0;
  int nnode=0;
  int nnz_diag=0;
  int nnz_off_diag=0;
  MPI_Init(NULL, NULL);
  int rank;
  int npe;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &npe);


  //  Memory allocates dynamically using malloc() 
  //indptr_diag = (int*)malloc(size * sizeof(int)); 


  //different processes have different rows
  if(rank==0){
    m = 3;
    n = 3;
    nnz_diag = 6;
    nnz_off_diag = 6;
    nghost = 4;
    nnode = 1;
  }else if(rank==1){
    m = 3;
    n = 3;
    nnz_diag = 9;
    nnz_off_diag = 4;
    nghost = 4;
    nnode = 2;
  }else if(rank == 2){
    m = 2;
    n = 2;
    nnz_diag = 2;
    nnz_off_diag = 8;
    nghost = 6;
    nnode = 1;
  }
  //allocate arrays
  int indptr_diag[m+1];
  int cols_diag[nnz_diag];
  double vals_diag[nnz_diag];
  int indptr_off_diag[m+1];
  int cols_off_diag[nnz_off_diag];
  double vals_off_diag[nnz_off_diag];
  int ghosts[nghost];
  double sol[n+nghost];
  double resid[n+nghost];
  double scale_vect[n+nghost];
  int ndof_node[nnode];
  //set values
  if(rank==0){
    indptr_diag[0] = 0; indptr_diag[1] =2; indptr_diag[2] = 4; indptr_diag[3] = 6;
    cols_diag[0] = 0; cols_diag[1] = 1; cols_diag[2] = 1; cols_diag[3] = 2;
    cols_diag[4] = 0; cols_diag[5] = 2;
    vals_diag[0] = 1; vals_diag[1] = 2; vals_diag[2] = 5; vals_diag[3] = 6;
    vals_diag[4] = 9; vals_diag[5] =10;
    indptr_off_diag[0] = 0; indptr_off_diag[1] = 2; indptr_off_diag[2] = 4; indptr_off_diag[3] = 6;
    cols_off_diag[0] = 4; cols_off_diag[1] = 7; cols_off_diag[2] =3;
    cols_off_diag[3] = 6; cols_off_diag[4] = 3; cols_off_diag[5] =6;
    vals_off_diag[0] = 3; vals_off_diag[1] = 4; vals_off_diag[2] = 7;
    vals_off_diag[3] = 8; vals_off_diag[4] = 11; vals_off_diag[5] = 12;
    //for now have them ordred, but need to try out of order too
    ghosts[0] = 3; ghosts[1] = 4; ghosts[2] = 6; ghosts[3] = 7;
    //set rhs which is resid
    resid[0] = 1; resid[1] = 1.5; resid[2] = .1; resid[3] = 1; resid[4] = 1;
    resid[5] = 1; resid[6] = 1;
    //set initial guess
    sol[0] = 0; sol[1] = 0; sol[2] = 0; sol[3] = 0; sol[4] = 0;
    sol[5] = 0; sol[6] = 0;
    //scale vector
    scale_vect[0] = 0; scale_vect[1] = 0; scale_vect[2] = 0; scale_vect[3] = 0; scale_vect[4] = 0;
    scale_vect[5] = 0; scale_vect[6] = 0;
    //for block size on preconditioner
    ndof_node[0] = 3;
  }else if (rank==1){
    indptr_diag[0] = 0; indptr_diag[1] = 3; indptr_diag[2] = 6; indptr_diag[3] = 9;
    cols_diag[0] = 0; cols_diag[1] = 1; cols_diag[2] =2; cols_diag[3] =0;
    cols_diag[4] = 1; cols_diag[5] = 2; cols_diag[6] = 0; cols_diag[7] =1; cols_diag[8] =2;
    vals_diag[0] = 15; vals_diag[1] = 16; vals_diag[2] = 17;
    vals_diag[3] = 19; vals_diag[4] = 20; vals_diag[5] = 21;
    vals_diag[6] = 22; vals_diag[7] = 23; vals_diag[8] = 1;
    indptr_off_diag[0] = 0; indptr_off_diag[1] = 2;
    indptr_off_diag[2] = 3; indptr_off_diag[3] = 4;
    cols_off_diag[0] = 0; cols_off_diag[1] = 2;
    cols_off_diag[2] = 1; cols_off_diag[3] = 6;
    vals_off_diag[0] = 13; vals_off_diag[1] = 14;
    vals_off_diag[2] = 18; vals_off_diag[3] = 24;
    ghosts[0] = 0; ghosts[1] = 1; ghosts[2] = 2; ghosts[3] = 6;
    resid[0] = 0.5; resid[1] = 0.2; resid[2] = 3.0; resid[3] = 1; resid[4] = 1;
    resid[5] = 1; resid[6] = 1;
    //set initial guess
    sol[0] = 0; sol[1] = 0; sol[2] = 0; sol[3] = 0; sol[4] = 0;
    sol[5] = 0; sol[6] = 0;
    //scale vector
    scale_vect[0] = 0; scale_vect[1] = 0; scale_vect[2] = 0; scale_vect[3] = 0; scale_vect[4] = 0;
    scale_vect[5] = 0; scale_vect[6] = 0;
    ndof_node[0] = 1; ndof_node[1] = 2;
  }else if (rank==2){
    indptr_diag[0] = 0; indptr_diag[1] = 1; indptr_diag[2] = 2;
    cols_diag[0] = 0; cols_diag[1] =1;
    vals_diag[0] = 29; vals_diag[1] =34;
    indptr_off_diag[0] = 0; indptr_off_diag[1] = 4; indptr_off_diag[2] = 8;
    cols_off_diag[0] = 0; cols_off_diag[1] = 1; cols_off_diag[2] = 2;
    cols_off_diag[3] = 5; cols_off_diag[4] = 0; cols_off_diag[5] = 3;
    cols_off_diag[6] = 4; cols_off_diag[7] = 5;
    vals_off_diag[0] = 25; vals_off_diag[1] = 26; vals_off_diag[2] = 27;
    vals_off_diag[3] = 28; vals_off_diag[4] = 30; vals_off_diag[5] = 31;
    vals_off_diag[6] = 32; vals_off_diag[7] =33;
    ghosts[0] = 0; ghosts[1] = 1; ghosts[2] = 2; ghosts[3] = 3;
    ghosts[4] = 4; ghosts[5] = 5;
    resid[0] = 1.5; resid[1] = 10; resid[2] = 1; resid[3] = 1; resid[4] = 1;
    resid[5] = 1; resid[6] = 1; resid[7] = 1;
    //set initial guess
    sol[0] = 0; sol[1] = 0; sol[2] = 0; sol[3] = 0; sol[4] = 0;
    sol[5] = 0; sol[6] = 0; sol[7] = 0;
    //scale vector
    scale_vect[0] = 0; scale_vect[1] = 0; scale_vect[2] = 0; scale_vect[3] = 0; scale_vect[4] = 0;
    scale_vect[5] = 0; scale_vect[6] = 0; scale_vect[7] = 0;
    ndof_node[0] = 2;
  }

  //update residual
  comm_update_double(resid,m,3,rank);
  //Global sizes
  //M = 8;
  //N = 8;

  //create this array directly with split arrays
  /*
            1  2  0  |  0  3  0  |  0  4
    Proc0   0  5  6  |  7  0  0  |  8  0
            9  0 10  | 11  0  0  | 12  0
    -------------------------------------
           13  0 14  | 15 16 17  |  0  0
    Proc1   0 18  0  | 19 20 21  |  0  0
            0  0  0  | 22 23  1  | 24  0
    -------------------------------------
    Proc2  25 26 27  |  0  0 28  | 29  0
           30  0  0  | 31 32 33  |  0 34

  */
  
//  for(k=0;k<nnz_off_diag;k++){
//    printf("Rank %d, unscaled vals_off_diag[%d] = %f\n",rank,k,vals_off_diag[k]);
//  }

  //now create routine that does the preconcitioning and then solve
  //scales matrix, rhs(residual), solution, and sets scale_vec
  scale_linear_system(indptr_diag, cols_diag, vals_diag, indptr_off_diag, cols_off_diag,
     vals_off_diag, resid, sol, scale_vect, m, n+nghost,rank, ghosts, nghost);

  //factor the matrix preconditioner
  status = prep_umfpack(indptr_diag,cols_diag,vals_diag, sol, resid, m);
  //and factors with umfpack
  //for(k=0;k<nnz_diag;k++){
  //  printf("Rank %d, vals_diag[%d] = %f\n",rank,k,vals_diag[k]);
  //}

  //for(k=0;k<nnz_off_diag;k++){
  //  printf("Rank %d, vals_off_diag[%d] = %f\n",rank,k,vals_off_diag[k]);
  //}

  //actually solve matrix
  status = solve_linear_sys_bcgstab(sol, indptr_diag, cols_diag, vals_diag, indptr_off_diag, cols_off_diag,
     vals_off_diag, resid, scale_vect, m, n+nghost,rank, ghosts, nghost);

  //print scale vec to see
  //for(k=0;k<m+nghost;k++){
  //  printf("Rank %d, scale_vec[%d] = %f\n",rank,k,scale_vect[k]);
  //}
  ////print out scaled matrix to see what happened
  //for(k=0;k<nnz_diag;k++){
  //  printf("Rank %d, scaled nnz_diag[%d] = %f\n",rank,k,vals_diag[k]);
  //}

  //now umfpack?
  MPI_Finalize();
  return 0;
}


void scale_linear_system(int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag,double *vals_off_diag, double *residual, 
  double *sol, double *scale_vect, int local_size, int size_with_ghosts, int rank,
  int *ghosts, int nghost){
  
  int i,j;
  double max_val_diag, max_val_off_diag;

  //first step is find max(abs(row)) and store in scale_vect, includes ghosts using MPI
  for (i=0;i<local_size;i++){
    max_val_diag = find_max_in_row(indptr_diag,vals_diag,i);
    max_val_off_diag = find_max_in_row(indptr_off_diag,vals_off_diag,i);
    
    if (max_val_diag > max_val_off_diag){
      scale_vect[i] = max_val_diag;
    }else{
      scale_vect[i] = max_val_off_diag;
    }

  }

  //ghosts of scale vect probably needed to be communicated
  //previously used comm_update_double
  //hardcoded for now
  comm_update_double(scale_vect, local_size, 3,rank);

  //now use scale factor as 1/sqrt(scale_vect), have check for small numbrs (may want warning statement)
  double small = SOLV_TOL;    /* too small a number for a norm */
  double small_factor;        /* small^(-0.5) */
  small_factor = 1.0 / sqrt(small);

  for (i = 0; i < size_with_ghosts; i++)
    if (scale_vect[i] < small)
      scale_vect[i] = small_factor;
    else
      scale_vect[i] = 1.0 / sqrt(fabs(scale_vect[i]));

  /* scales the equations */
  
  double factor;
  for (i = 0; i < size_with_ghosts; i++) {
    factor = scale_vect[i];
    residual[i] *= factor;
    sol[i] /= factor;
  }


  int col_no, global_col_no;


  //this is using nodal blocks but what I am seeing this is just multiplying by scale_vec[i]*scale_vec[j]
  int ind_start, ind_end, nentry_row;
  for (i=0;i<local_size;i++){
    ind_start = indptr_diag[i];ind_end = indptr_diag[i+1];
    nentry_row = ind_end-ind_start;
    for(j=0;j<nentry_row;j++){
      col_no = cols_diag[ind_start+j];
      vals_diag[ind_start+j]*=scale_vect[i]*scale_vect[col_no];
    }
  }
  //now off diag
  //remember this gives global column numbers
  for (i=0;i<local_size;i++){
    ind_start = indptr_off_diag[i];ind_end = indptr_off_diag[i+1];
    nentry_row = ind_end-ind_start;
    for(j=0;j<nentry_row;j++){
      //this is global column number, need global to local mapping
      //this is easy with ghost numbers
      global_col_no = cols_off_diag[ind_start+j];
      //this is easy with ghost numbers
      //search ghost no's and find the number then add on local_size
      col_no = global_to_local(global_col_no,local_size,ghosts,nghost);
      //printf("Rank %d, local row %d, row entry %d, global col_no %d, col_no %d\n", rank, i, j, global_col_no,col_no);
      vals_off_diag[ind_start+j]*=scale_vect[i]*scale_vect[col_no];
    }
  }

}



int prep_umfpack(int *indptr_diag, int *cols_diag, double *vals_diag, double *sol, double *resid, int nrow){
  //load into umfpack, need to transpose because umfpack uses ccs format
  int status;

  double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];

  status = umfpack_di_symbolic (nrow, nrow, indptr_diag, cols_diag, vals_diag, &Symbolic, Control, Info);
  status = umfpack_di_numeric (indptr_diag, cols_diag, vals_diag, Symbolic, &Numeric, Control, Info) ;
  //Adh Stops here, how do we want to call other stuff?
  return status;

}

int solve_umfpack(double *x, int *indptr_diag, int *cols_diag, double *vals_diag, double *b, int nrow){
  int status;
  double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];
  //try this to account for CSR format, use solution as RHS for some reason
  status = umfpack_di_solve (UMFPACK_Aat, indptr_diag, cols_diag, vals_diag, x, b, Numeric, Control, Info);
  return status;
}

//now call bigcstab, this is also assuming there is some scaling prior to solve
//if no scaling, the scale_vect should be all 1s
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

  /* allocates memory if needed */
  // need to do different allocation later
  double r[size_with_ghosts];
  double p[size_with_ghosts];
  double Ap[size_with_ghosts];
  double Mp[size_with_ghosts];
  double q[size_with_ghosts];
  double s[size_with_ghosts];
  double As[size_with_ghosts];
  double Ms[size_with_ghosts];
  double x0[size_with_ghosts];


  /* zeroes the arrays */
  //need to store this to shift solution after solve
  copy_array_to_from(x0, x, size_with_ghosts);

  //initialize arrays to 0
  zero_dbl_array(x, size_with_ghosts);
  zero_dbl_array(r, size_with_ghosts);
  zero_dbl_array(q, size_with_ghosts);
  zero_dbl_array(p, size_with_ghosts);
  zero_dbl_array(Mp, size_with_ghosts);
  zero_dbl_array(Ap, size_with_ghosts);
  zero_dbl_array(s, size_with_ghosts);
  zero_dbl_array(As, size_with_ghosts);
  zero_dbl_array(Ms, size_with_ghosts);

  copy_array_to_from(p, b, size_with_ghosts);
  //in original routine but doesnt seem necessary
  zero_dbl_array(b,size_with_ghosts);
  int i;
  //apply left preconditioning to RHS, this is stored in residual
  solve_umfpack(b, indptr_diag, cols_diag, vals_diag, p, local_size);
  //zero out rhs
  zero_dbl_array(p, size_with_ghosts);
  //update ghosts

  /* Now Form Actual Residual */
  /* Assuming x=0 as the initial guess would have given r=b (the rhs). Easy, right? */
  /* However, we miss an opportunity if our preconditioner is meaningful. */
  /* This will basically take r = S b - S A x, where x = S b, as the initial guess. */
  //should it really be updating the ghosts?
  copy_array_to_from(x, b, size_with_ghosts);
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
  copy_array_to_from(r, b, size_with_ghosts);
  y_plus_ax(r, -1, Ap, size_with_ghosts );
  zero_dbl_array(Mp, size_with_ghosts);
  zero_dbl_array(Ap, size_with_ghosts);
  copy_array_to_from(q, r, size_with_ghosts);
  //print vector

  //take linf norms of vectors
  rnorm = l_infty_norm(local_size, r);
  bnorm = l_infty_norm(local_size, b);
  //if it is in parallel. collect with allreduce
  rnorm = get_global_max(rnorm);
  bnorm = get_global_max(bnorm);


  rhop = 1.0;
  alpha = 1.0;
  omega = 1.0;
  it = 0;
  min_conv_tol = 1e-10;
  conv_tol = 1e-5;
  rho = dot_dbl_array(local_size, r, q);
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
    y_plus_ax(p, -omega, Ap,size_with_ghosts );
    scale_dbl_array(p,beta, size_with_ghosts);
    y_plus_ax(p, 1, r, size_with_ghosts);
    //for (i=0;i<local_size;i++){
    //printf("Rank %d, p[%d] = %f\n",rank,i,p[i]);
    //}
    //need comm_update_double befor matvec mult
    comm_update_double(p, local_size, 3, rank);

    // (d) //
    zero_dbl_array(Mp, size_with_ghosts);
    zero_dbl_array(Ap,size_with_ghosts);
    split_CSR_mat_vec_mult(Mp, indptr_diag, cols_diag, vals_diag, indptr_off_diag, cols_off_diag, vals_off_diag,
    p, local_size,ghosts,nghost);
    //if MPI
    //preconfitioner is local only
    //comm_update_double(Mp, local_size, 3, rank);
    //next, apply preconditioner to Mp, save in Ap
    solve_umfpack(Ap, indptr_diag, cols_diag, vals_diag, Mp, local_size);

    // (e) //
    gamma = dot_dbl_array(local_size, q, Ap);
    gamma = messg_dsum(gamma);
    alpha = rho / gamma;
    //printf("gamma, alpha %f, %f\n",gamma,alpha);

    /* (f) */
    copy_array_to_from(s, r, size_with_ghosts);
    y_plus_ax(s, -alpha, Ap, size_with_ghosts);
    zero_dbl_array(Ms, size_with_ghosts);
    zero_dbl_array(As, size_with_ghosts);
    //need comm_update_double befor matvec mult
    comm_update_double(s, local_size, 3, rank);
    split_CSR_mat_vec_mult(Ms, indptr_diag, cols_diag, vals_diag, indptr_off_diag, cols_off_diag, vals_off_diag,
    s, local_size,ghosts,nghost);
    solve_umfpack(As, indptr_diag, cols_diag, vals_diag, Ms, local_size);

    /* (g) */
    omega = dot_dbl_array(local_size, As, s);
    gamma = dot_dbl_array(local_size, As, As);
    //if MPI is on
    omega = messg_dsum(omega);
    gamma = messg_dsum(gamma);

    omega = omega / gamma;
    rhop = rho;
    rho = -1.0 * omega * dot_dbl_array(local_size, As, q);
    //only MPI
    rho = messg_dsum(rho);

    /* (h) */
    y_plus_ax(x, alpha, p, size_with_ghosts);
    y_plus_ax(x, omega, s, size_with_ghosts);
        
    /* (i) */
    copy_array_to_from(r, s, size_with_ghosts);
    y_plus_ax(r, -omega, As, size_with_ghosts);
    /*rnorm = solv_l2_norm_scaled(my_ndof_solv, r, (double)(global_nnode)); */
    rnorm = l_infty_norm(local_size, r);
    //only do this in parallel
    rnorm = get_global_max(rnorm);
    //printf("Rnorm, %f\n",rnorm);

  }
  printf("IT = %d\n\n",it);

  //undo scaling (need to add in error checks)
  for (i = 0; i < local_size; i++) {
        x[i] += x0[i];
        /* solv_u[i] *= scale_vect[i]; UNDO SCALING */
    }
  for (i = 0; i <local_size; i++){
        x[i] *= scale_vect[i]; /* UNDO SCALING */
    }
  comm_update_double(x, local_size, 3, rank);
  for (i = 0; i <local_size; i++){
    printf("Rank %d Final u[%d] = %f\n",rank,i,x[i]);
  }

  //add later to see if converged or not
  int status = 0;
  return status;
}


double get_global_max(double x){
  double x_send = 0.0;        /* the variables for the pass */
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */
  x_send = x;
  ierr_code = MPI_Allreduce(&x_send, &x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return x;
}

void scale_dbl_array(double *v,  double alpha, int n){
  int i;
  //use blas maybe later
  /* multiplies by a scalar */
  for(i = 0; i < n; i++)
        v[i] *= alpha;
}

void split_CSR_mat_vec_mult(double *Ax, int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag, double *vals_off_diag,
  double *x, int nrows, int *ghosts, int nghost){

  //consult with ppl to optimize later on
  //see discussion on https://scicomp.stackexchange.com/questions/27977/how-can-i-speed-up-this-code-for-sparse-matrix-vector-multiplication
  //does mat/vec multiplication between split CSR with vector x and returns Ax
  int i, j, k;
  double Ax_i;
  //loop over row
  for(i=0;i<nrows;i++){
    Ax_i = 0.0;
    for(j=indptr_diag[i];j<indptr_diag[i+1];j++){
      //diag blaock cols_diag are local numbers only so no issue here
      Ax_i += vals_diag[j]*x[cols_diag[j]];
    }

    for(k=indptr_off_diag[i];k<indptr_off_diag[i+1];k++){
      //off diag is a little more complex
      Ax_i += vals_off_diag[k]*x[global_to_local(cols_off_diag[k], nrows, ghosts, nghost)];
    }

    Ax[i] = Ax_i;

  }


}

void umfpack_clear(void){
  //clear memory
  umfpack_di_free_symbolic (&Symbolic) ;
  umfpack_di_free_numeric (&Numeric) ;
}

int global_to_local(int global_col_no,int local_size,int *ghosts, int nghost){
  //searches the array of ghosts and returns the local index
  //assumes on a process, dofs owned by process come first and then ghosts next
  //not assuming ghosts are ordered here
  int i;
  //printf("len of array, %d\n",len_array);
  for (i=0;i<nghost;i++){
    if(ghosts[i] == global_col_no){
      return i + local_size;
    }
  }
  return -1;

}

double find_max_in_row(int *indptr, double *vals,int row_num){
  int ind1 = indptr[row_num];
  int ind2 = indptr[row_num+1];
  int nind = ind2-ind1;
  int j;


  double max_val = 0;
  double temp;
  for(j=0;j<nind;j++){
      temp = fabs(vals[j+ind1]);
      if(temp > max_val){
        max_val = temp;
      }
  }

  return max_val;

}


void zero_dbl_array(
                   double *v,      /* the vector */
                   int n     /* the number of degrees of freedom */
)
{
    int ii = 0;     /* loop counter */
    
    /* zeroes the array */
    for(ii = 0; ii < n; ii++)
    {
        v[ii] = 0;
    }
    return;
}

void copy_array_to_from(double *vec_to, double *vec_from, int n){
  //better routine elsewhere
  int ii;
  for(ii=0;ii<n;ii++){
    vec_to[ii] = vec_from[ii];
  }
  return;
}

//\brief Performs daxpy operation for (double) vectors, y = \alpha x  + y$
void y_plus_ax( double *y,     /* y vector (input and output) */
                double alpha,     /* the scalar */
                double *x,      /* x vector */
                int n      /* the number of degrees of freedom */
                
)
{
//can use blas or similar for this
int i = 0;      /* loop counter */
    
    /* adds the arrays */
    for(i = 0; i < n; i++)
        y[i] += alpha * x[i];
    return;
}

/*!
 \brief Returns the \f$l_\infty\f (infinity) norm of a vector (maximum absolute value)
 
 \param n the length of the vectors (int)
 \param v1 pointer to the vector (double, length n)
 */
double l_infty_norm(
                       int n,     /* the number of degrees of freedom */
                       double *v1     /* the vector */
                       )
{
    double value = 0.0;   /* the partial sum */
    int ii = 0;     /* loop counter */
    
    /* computes the value for this processor */
    for(ii = 0; ii < n; ii++)
    {
        value = max_dbl(value, fabs(v1[ii]));
    }
    
    /* returns the maximum */
    return value;
}

double max_dbl(double a, double b)
{
    return a > b ? a : b;
}

/*!
 \brief Returns the inner product, \f$(x \cdot y)\f$, of two (double) vectors
 
 Sums among all processors.
 
 \param n the length of the vectors (int)
 \param x pointer to the first vector (double, length \a n)
 \param y pointer to the second vector (double, length \a n)
 */
double dot_dbl_array(
                int n,      /* the number of degrees of freedom */
                double *x,      /* the first vector */
                double *y     /* the second vector */
)
{
    double value = 0.0;   /* the partial sum */
    int i = 0;      /* loop counter */
    
    for(i = 0; i < n; i++)
        value += x[i] * y[i];

    /* returns the sum */
    return (value);
}


double messg_dsum(double x)
{
  double x_send = 0.0;        /* the variables for the pass */
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */

  x_send = x;
  ierr_code = MPI_Allreduce(&x_send, &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (x);
}

/*!
   \brief Update the Ghost Node Values (Get Info From Owning Processor)

   *hard coded for this example
 */
void comm_update_double(double *vec,  /* the vector to be updated */
                        int size_v,
                        int npe,
                        int rank //number of processors
  )
{
  //double *bpntr = NULL;         /* pointer to allow for de-referencing */
  int i_processor = 0, j_processor = 0, ii=0;  /* loop counter */

  //just send whole vector i guess?
  double temp_send[7],temp_recv[7];
  for (ii = 0;ii<7;ii++){
    if(ii<size_v){
      temp_send[ii] = vec[ii];
      temp_recv[ii] = 0;
    }else{temp_send[ii] = 0;
          temp_recv[ii] = 0;}
  }


  /* Post the Receives */
  for (i_processor = 0; i_processor < npe; i_processor++){
    if (rank ==i_processor){
      for (j_processor = 0; j_processor < npe; j_processor++){
        if (rank != j_processor){
          MPI_Send(temp_send, 7, MPI_DOUBLE, j_processor, j_processor, MPI_COMM_WORLD);
        }
      }
    }else{
       MPI_Recv(temp_recv, 7, MPI_DOUBLE, i_processor, rank, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
       //now move temp entries to appropriate vec!!!
       if (rank==0){
          if(i_processor == 1){
            vec[3] = temp_recv[0];
            vec[4] = temp_recv[1];
          }else if(i_processor == 2){
            vec[5] = temp_recv[0];
            vec[6] = temp_recv[1];

          }
       }else if(rank==1){
          if(i_processor == 0){
            vec[3] = temp_recv[0];
            vec[4] = temp_recv[1];
            vec[5] = temp_recv[2];
          }else if(i_processor == 2){
            vec[6] = temp_recv[0];

          }
       }else if(rank==2){
        if(i_processor == 0){
            vec[2] = temp_recv[0];
            vec[3] = temp_recv[1];
            vec[4] = temp_recv[2];
          }else if(i_processor == 1){
            vec[5] = temp_recv[0];
            vec[6] = temp_recv[1];
            vec[7] = temp_recv[2];

          }

       }
    }
  }



  return;
}

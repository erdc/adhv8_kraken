#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <umfpack.h>

#define SOLV_TOL 1e-5
#define MAX_NODAL_DOF 4
void solv_linear_sys_setup(int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag,double *vals_off_diag, double *residual, 
  double *sol, double *scale_vect, int local_size, int size_with_ghosts, int rank,int *ghosts);
void solv_linear_sys_bcgstab(int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag,double *vals_off_diag, double *residual, 
  double *sol, double *scale_vect, int local_size, int size_with_ghosts, int rank,int *ghosts);
double find_max_in_row(int *indptr, double *vals,int row_num);
void comm_update_double(double *vec,int size_v, int npe,int rank);
int global_to_local(int global_col_no,int local_size,int *ghosts);
void prep_umfpack(int *indptr_diag, int *cols_diag, double *vals_diag, double *sol, double *resid, int nrow);
void solve_umfpack(int *indptr_diag, int *cols_diag, double *vals_diag, double *sol, double *resid, int nrow);
void umfpack_clear();
void solv_init_dbl(int n,double *v);

//keep umfpack calls in one script
static  void *Symbolic, *Numeric;

int main(int argc, char **argv) {
  int k;
  int ierr;
  int nghost = 0;
  int m=0;
  int n=0;
  int nnode=0;
  int nnz_diag=0;
  int nnz_off_diag=0;
  int M,N;
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
    resid[0] = 1; resid[1] = 1; resid[2] = 1; resid[3] = 1; resid[4] = 1;
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
    resid[0] = 1; resid[1] = 1; resid[2] = 1; resid[3] = 1; resid[4] = 1;
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
    resid[0] = 1; resid[1] = 1; resid[2] = 1; resid[3] = 1; resid[4] = 1;
    resid[5] = 1; resid[6] = 1; resid[7] = 1;
    //set initial guess
    sol[0] = 0; sol[1] = 0; sol[2] = 0; sol[3] = 0; sol[4] = 0;
    sol[5] = 0; sol[6] = 0; sol[7] = 0;
    //scale vector
    scale_vect[0] = 0; scale_vect[1] = 0; scale_vect[2] = 0; scale_vect[3] = 0; scale_vect[4] = 0;
    scale_vect[5] = 0; scale_vect[6] = 0; scale_vect[7] = 0;
    ndof_node[0] = 2;
  }

  //Global sizes
  M = 8;
  N = 8;

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
  

  //now create routine that does the preconcitioning and then solve
  //scales matrix and factors with umfpack
  solv_linear_sys_setup(indptr_diag, cols_diag, vals_diag, indptr_off_diag, cols_off_diag,
     vals_off_diag, resid, sol, scale_vect, m, n+nghost,rank, ghosts);

  //actually solve matrix
  solv_linear_sys_bcgstab(indptr_diag, cols_diag, vals_diag, indptr_off_diag, cols_off_diag,
     vals_off_diag, resid, sol, scale_vect, m, n+nghost,rank, ghosts);

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


//preconditioning
void solv_linear_sys_setup(int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag,double *vals_off_diag, double *residual, 
  double *sol, double *scale_vect, int local_size, int size_with_ghosts, int rank,
  int *ghosts){
  
  int i,j,k,l;
  double max_val_diag, max_val_off_diag;
  double row_factor[MAX_NODAL_DOF];

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

  //ghosts of scale vect needs to be communicated
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
  int ctr=0,ctr2=0,ctr3=0;
  int nodal_block_size;
  int row_ind_start;
  bool finding_col;
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
      global_col_no = cols_diag[ind_start+j];
      //this is easy with ghost numbers
      //search ghost no's and find the number then add on local_size
      col_no = global_to_local(global_col_no,local_size,ghosts);
      vals_off_diag[ind_start+j]*=scale_vect[i]*scale_vect[col_no];
    }
  }
  

  //factor the matrix preconditioner
  prep_umfpack(indptr_diag,cols_diag,vals_diag, sol, residual, local_size);

  //show what happened to scalevec?
  //for(k=0;k<size_with_ghosts;k++){
  //  printf("Rank %d, scalevec[%d] = %f\n",rank,k,scale_vect[k]);
  //}

  //show what happened to residual?
  //for(k=0;k<local_size;k++){
  //  printf("Rank %d, scaled residual[%d] = %f\n",rank,k,residual[k]);
  //}
  //show what happened to nnz block?
  //for(k=0;k<4;k++){
  //  printf("Rank %d, scaled vals diag[%d] = %f\n",rank,k,vals_diag[k]);
  //}

  //does stuff get modified in sparse matrix too?
  //No
  //for(k=0;k< 7;k++){
  //  printf("Rank %d, scaled vals after umfpack[%d] = %f\n",rank,k,vals_diag[k]);
  //}

}


void prep_umfpack(int *indptr_diag, int *cols_diag, double *vals_diag, double *sol, double *resid, int nrow){
  //load into umfpack, need to transpose because umfpack uses ccs format
  int status, sys, nz;

  double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];

  status = umfpack_di_symbolic (nrow, nrow, indptr_diag, cols_diag, vals_diag, &Symbolic, Control, Info);
  status = umfpack_di_numeric (indptr_diag, cols_diag, vals_diag, Symbolic, &Numeric, Control, Info) ;
  //Adh Stops here, how do we want to call other stuff?

}

void solve_umfpack(int *indptr_diag, int *cols_diag, double *vals_diag, double *sol, double *resid, int nrow){
  int status, sys, nz;
  double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];
  //try this to account for CSR format, use solution as RHS for some reason
  status = umfpack_di_solve (UMFPACK_Aat, indptr_diag, cols_diag, vals_diag, sol, resid, Numeric, Control, Info);
}

//now call bigcstab
void solve_linear_sys_bcgstab(int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag,double *vals_off_diag, double *residual, 
  double *sol, double *scale_vect, int local_size, int size_with_ghosts, int rank,
  int *ghosts){

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
  isize_prev = isize;
  isize = ndof_solv;
  double solv_r[size_with_ghosts];
  double solv_p[size_with_ghosts];
  double solv_Ap[size_with_ghosts];
  double solv_Mp[size_with_ghosts];
  double solv_q[size_with_ghosts];
  double solv_s[size_with_ghosts];
  double solv_As[size_with_ghosts];
  double solv_Ms[size_with_ghosts];
  double solv_Au[size_with_ghosts];
  double solv_u0[size_with_ghosts];
  double solv_ptmp[size_with_ghosts];


  /* zeroes the arrays */
  solv_copy(size_with_ghosts, sol, solv_u0);
  solv_init_dbl(size_with_ghosts, sol);
  solv_init_dbl(size_with_ghosts, solv_r);
  solv_init_dbl(size_with_ghosts, solv_q);
  solv_init_dbl(size_with_ghosts, solv_p);
  solv_init_dbl(size_with_ghosts, solv_Mp);
  solv_init_dbl(size_with_ghosts, solv_Ap);
  solv_init_dbl(size_with_ghosts, solv_s);
  solv_init_dbl(size_with_ghosts, solv_As);
  solv_init_dbl(size_with_ghosts, solv_Ms);

  solv_copy(size_with_ghosts, residual, solv_p);
  //in original routine but doesnt seem necessary
  //solv_init_dbl(size_with_ghosts, residual);
  //why??
  //solv_copy(size_with_ghosts, solv_p, residual);

  //apply left preconditioning to RHS, this is stored in residual
  solve_umfpack(indptr_diag, cols_diag, vals_diag, residual, solv_p, local_size);
  //zero out rhs
  solv_init_dbl(ndof_solv, solv_p);

  /* Now Form Actual Residual */
  /* Assuming x=0 as the initial guess would have given r=b (the rhs). Easy, right? */
  /* However, we miss an opportunity if our preconditioner is meaningful. */
  /* This will basically take r = S b - S A x, where x = S b, as the initial guess. */
  //should it really be updating the ghosts?
  solv_copy(size_with_ghosts, residual, sol);


}


void umfpack_clear(){
  //clear memory
  umfpack_di_free_symbolic (&Symbolic) ;
  umfpack_di_free_numeric (&Numeric) ;
}

int global_to_local(int global_col_no,int local_size,int *ghosts){
  //searches the array of ghosts and returns the local index
  //assumes on a process, dofs owned by process come first and then ghosts next
  //not assuming ghosts are ordered here
  int i;
  int len_array = sizeof(&(ghosts))/sizeof(int);
  for (i=0;i<len_array;i++){
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


void solv_init_dbl(
                   int n,     /* the number of degrees of freedom */
                   double *v      /* the vector */
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

void solv_copy(int n, double *vec_from, double *vec_to){
  //better routine elsewhere
  int ii;
  for(ii=0;ii<n;ii++){
    vec_to[ii] = vec_from[ii];
  }
  return;
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
  double *bpntr = NULL;         /* pointer to allow for de-referencing */
  int i_processor = 0, j_processor = 0, ii = 0, jj = 0;  /* loop counter */
  int icnt = 0;                 /* location in the buffer */

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

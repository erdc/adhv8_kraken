/*! \file inhouse.c This file tests the BiCGStab solver for split CSR matrix */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function tests the BiCGStab solver
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int inhouse_test(int argc, char **argv) {
  int status;
  int nghost = 0;
  int m=0;
  int n=0;
  int nnode=0;
  int nnz_diag=0;
  int nnz_off_diag=0;
  int rank=0;
  
#ifdef _MPI
  int npe=1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &npe);
#endif

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
  //comm_update_double(resid,m,3,rank);
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
  status = prep_umfpack(indptr_diag,cols_diag,vals_diag, m);
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
  return status;
}




/*! \file  scale_linear_system.c This file has functions responsible for scaling the linear system in split CSR format based on max row entry */
#include "adh.h"
static double find_max_in_row(int *indptr, double *vals,int row_num);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Scales a linear of systems based on max row entry, modifies matrix values and saves a scale_vect for later reverse scaling
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in] indptr_diag (int*) - indeces of first/last index for each row in diagonal block (local to process) in split CSR matrix
 *  @param[in] cols_diag (int*) - column numbers (local to process) of diagonal block in split CSR matrix
 *  @param[in,out] vals_diag (double*) - array of doubles that are the nonzero entries of diagonal block (local to process) in split CSR matrix
 *  @param[in] indptr_off_diag (int*) - indeces of first/last index for each row in off-diagonal block (local to process) in split CSR matrix
 *  @param[in] cols_off_diag (int*) - column numbers (global) of off-diagonal block in split CSR matrix
 *  @param[in,out] vals_off_diag (double*) - array of doubles that are the nonzero entries of off-diagonal block (local to process) in split CSR matrix
 *  @param[in,out] b (double*) - scaled rhs of the linear system Ax=b
 *  @param[in,out] x (double*) - scaled solution of the linear system Ax=b
 *  @param[in] scale_vect (double*) - the vector that scaled the system of equations prior to solve, used for rescaling after solve is finished
 *  @param[in] local_size (int) - number of equations (rows in the matrix) owned by the process
 *  @param[in] size (int) - number of rows + number of ghost d.o.f.s present
 *  @param[in] rank (int) - MPI rank
 *  @param[in] ghosts (int*) - array of global dof numbers that are ghosts on the process
 *  @param[in] nghost (int) - number of ghost d.o.f.s on the process
 * \returns integer that is a code if the solve worked or not
 * \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void scale_linear_system(int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag,double *vals_off_diag, double *b, 
  double *x, double *scale_vect, int local_size, int size, int rank,
  int *ghosts, int nghost){
  
  int i,j;
  double max_val_diag, max_val_off_diag;

  //first step is find max(abs(row)) and store in scale_vect, includes ghosts using MPI
  if(indptr_off_diag!=NULL){
    for (i=0;i<local_size;i++){
      max_val_diag = find_max_in_row(indptr_diag,vals_diag,i);
      max_val_off_diag = find_max_in_row(indptr_off_diag,vals_off_diag,i);
    
      if (max_val_diag > max_val_off_diag){
        scale_vect[i] = max_val_diag;
      }else{
        scale_vect[i] = max_val_off_diag;
      }

    }
  }else{
    for (i=0;i<local_size;i++){
      max_val_diag = find_max_in_row(indptr_diag,vals_diag,i);
      scale_vect[i] = max_val_diag;
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

  for (i = 0; i < size; i++)
    if (scale_vect[i] < small)
      scale_vect[i] = small_factor;
    else
      scale_vect[i] = 1.0 / sqrt(fabs(scale_vect[i]));

  /* scales the equations */
  
  double factor;
  for (i = 0; i < size; i++) {
    factor = scale_vect[i];
    b[i] *= factor;
    x[i] /= factor;
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
  if (indptr_off_diag!=NULL){
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

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Rescales solution after a linear solve of system that was scaled using scale_linear_system
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] x (double*) - solution vector that comes from solving the prescaled Ax=b system
 *  @param[in] x0 (double*) - the initial solution passed to iterative solver (after scaling)
 *  @param[in] scale_vect (double*) - the vector that was used for scaling
 *  @param[in] local_size (int) - number of d.o.f owned by process
 * \note To be paired with a call to scale_linear_system BEFORE solve
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void unscale_linear_system(double *x, double *x0, double *scale_vect, int local_size){
  int i;
  for (i = 0; i < local_size; i++) {
      x[i] += x0[i];
      /* solv_u[i] *= scale_vect[i]; UNDO SCALING */
  }
  for (i = 0; i <local_size; i++){
      x[i] *= scale_vect[i]; /* UNDO SCALING */
  }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Short helper function that finds the max value in a CSR row, used for scaling matrix
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in] indptr (int*) - array of ints that are the row start/end indices of CSR matrix
 *  @param[in] vals (double*) - array of doubles that are the nonzero entries in CSR matrix
 *  @param[in] row_num (int) - the row number of interest
 * \returns (double) max value in the row_num row of the CSR matrix
 * \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
static double find_max_in_row(int *indptr, double *vals,int row_num){
  int ind1 = indptr[row_num];
  int ind2 = indptr[row_num+1];
  int j;
  double max_val = 0;
  double temp;
  for(j=ind1;j<ind2;j++){
      temp = fabs(vals[j]);
      if(temp > max_val){
        max_val = temp;
      }
  }
  return max_val;
}

/*! \file  matrix_utilities.c This file has various functions needed for solving sparse systems in split CSR format */
#include "adh.h" 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Small wrapper to initialize all structs pertaining to linear system
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] sm (SMODEL_SUPER*) - pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//void fe_allocate_initialize_linear_system(SMODEL_SUPER *sm) {//

//    if (sm->ndofs > sm->ndofs_old){    
//        allocate_adh_system(sm);
//#ifdef _PETSC
//        allocate_petsc_objects(sm);
//#endif
//            }
//}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Function to modify linear system to enforce Dirichlet conditions
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] sm (SMODEL_SUPER*) - pointer to an instant of the SuperModel struct
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void apply_Dirichlet_BC(SMODEL_SUPER *sm){
    //NOTE: THIS WILL BREAK IN PARALLEL, NEED TO ADD GHOSTS LATER

    //could illuciadate arguments, but will look at this later
    //modify linear system to have dirichlet conditions at
    //entries where bc_mask=YES
    //Method 1: doesnt require solution to have dirichlet conditions set
    //Method 2: assumes solution already has dirichlet condition values set

    //brute force for now, can find better way later
    int my_ndof = *(sm->my_ndofs);
    int row, row_start, row_end, row_entry,local_col_no;
    //Required for method 1 only
    //double temp;
    for (row=0;row<my_ndof;row++){
        row_start = sm->lin_sys->indptr_diag[row];
        row_end = sm->lin_sys->indptr_diag[row+1];
        //for every row, look for columns and subtract to RHS
        for(row_entry=row_start; row_entry<row_end; row_entry++){
            local_col_no = sm->lin_sys->cols_diag[row_entry];
            if(sm->bc_mask[local_col_no] == YES){
                //required for Method 1
                ////////////////////////////////////////////
                //temp = sm->lin_sys->vals_diag[row_entry];
                //adjust RHS
                //needed for method 1, for method 2 this should be unecceary
                //sm->residual[row] -= temp*(sm->dirichlet_data[local_col_no] - sm->sol[local_col_no]);
                ////////////////////////////////////////////

                //used in both method 1 and method 2
                //clear entry after moving
                sm->lin_sys->vals_diag[row_entry] = 0.0;
            }

        }

    }
    //go back through rows and clear rows that are no longer in equations
    //if it is dirichlet dof, wipe out row
    for (row=0;row<my_ndof;row++){
        if(sm->bc_mask[row] == YES){
            row_start = sm->lin_sys->indptr_diag[row];
            row_end = sm->lin_sys->indptr_diag[row+1];
            //assign dirichlet data to residual
            //Method 1
            //sm->residual[row] = sm->dirichlet_data[row] - sm->sol[row];
            //Method 2, requires sol to have bc already  there
            sm->lin_sys->residual[row] = 0.0;
            //wipe out row and set diagonal to 1 (rescaling diagonal shouldn't matter in method 2)
            for(row_entry=row_start; row_entry<row_end; row_entry++){
                local_col_no = sm->lin_sys->cols_diag[row_entry];
                //if non-diagonal set to 0, otherwise set diagonal to 1
                if(local_col_no == row){
                    sm->lin_sys->vals_diag[row_entry]= 1.0;
                }else{
                    sm->lin_sys->vals_diag[row_entry] = 0.0;
                }

            }

        }
    }

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Function to check for near zero entries and force dirichlet conditions there
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] sm (SMODEL_SUPER*) - pointer to an instant of the SuperModel struct
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void check_diag(SMODEL_SUPER *sm){
    //NOTE: THIS WILL BREAK IN PARALLEL, NEED TO ADD GHOSTS LATER

    //could illuciadate arguments, but will look at this later
    //modify linear system to have dirichlet conditions at
    //entries where bc_mask=YES
    //Method 1: doesnt require solution to have dirichlet conditions set
    //Method 2: assumes solution already has dirichlet condition values set

    //brute force for now, can find better way later
    int my_ndof = *(sm->my_ndofs);
    int row, row_start, row_end, index;
    //Required for method 1 only
    //double temp;
    //sarray_printScreen_int(sm->bc_mask, *(sm->ndofs), "bcmask old");
    for (row=0;row<my_ndof;row++){
        row_start = sm->lin_sys->indptr_diag[row];
        row_end = sm->lin_sys->indptr_diag[row+1];
        //find index of diagonal element in row
        index = binary_search_part(sm->lin_sys->cols_diag, row_start, row_end, row);
        //printf("Diagonal index of row[%d] %d (col %d)\n",row,index,sm->lin_sys->cols_diag[index]);
        //check if the entry is < some tolerance
        
        if (fabs(sm->lin_sys->vals_diag[index]) < SMALL){
            //(1) add to bc mask
            //(2) ensure dirichlet data is 0 (is this what we want?)
            sm->bc_mask[row] = YES;
            //sm->dirichlet_data[index] = 0.0;
        }

    }
    //sarray_printScreen_int(sm->bc_mask, *(sm->ndofs), "bcmask new");
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Performs mat/vec multiplication in split CSR format
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] Ax (double*) - the matrix-vector product Ax
 *  @param[in] indptr_diag (int*) - the first/last index of each row in CSR format (diagonal block wrt processor)
 *  @param[in] cols_diag (int*) - local (to process) column locations of each nonzero in CSR format (diagonal block wrt processor)
 *  @param[in] vals_diag (double*) - nonzero values in CSR matrix (diagonal block wrt processor)
 *  @param[in] indptr_off_diag (int*) - the first/last index of each row in CSR format (off-diagonal block wrt processor)
 *  @param[in] cols_off_diag (int*) - global column locations of each nonzero in CSR format (off-diagonal block wrt processor)
 *  @param[in] vals_off_diag (double*) - nonzero values in CSR matrix (off-diagonal block wrt processor)
 *  @param[in] x (double*) - the vector of the matrix-vector product
 *  @param[in] nrows (int) - the number of rows on this process (same as # owned d.o.fs)
 *  @param[in] ghosts (int*) - array of global dof numbers for ghost dofs
 *  @param[in] nghost (int) - length of ghosts array
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void split_CSR_mat_vec_mult(double *Ax, int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag, double *vals_off_diag,
  double *x, int nrows, int *ghosts, int nghost){

  //consult with ppl to optimize later on
  //see discussion on https://scicomp.stackexchange.com/questions/27977/how-can-i-speed-up-this-code-for-sparse-matrix-vector-multiplication
  //does mat/vec multiplication between split CSR with vector x and returns Ax
  int i, j, k;
  double Ax_i;

  if(indptr_off_diag!=NULL){
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
  }else{
    //loop over row
    for(i=0;i<nrows;i++){
      Ax_i = 0.0;
      for(j=indptr_diag[i];j<indptr_diag[i+1];j++){
        //diag blaock cols_diag are local numbers only so no issue here
        Ax_i += vals_diag[j]*x[cols_diag[j]];
      }
      
      Ax[i] = Ax_i;

    }
  }

}
//int* sarray_unique_int (int* first, int* last)
//{
//  if (first==last) return last;//
//  int* result = first;
//  while (++first != last)
//  {
//    if (!(*result == *first)) 
//      *(++result)=*first;
//  }
//  return ++result;
//}//
//int unique(int *arr, int size){
//    int sarray_unique_int_size = 1; // We'll start with the assumption that the first element is sarray_unique_int
//    for (int i = 1; i < size; i++) {
//        int is_sarray_unique_int = 1;
//        for (int j = 0; j < sarray_unique_int_size; j++) {
//            if (arr[i] == arr[j]) {
//                is_sarray_unique_int = 0; // We found a duplicate
//                break;
//            }
//        }
//        if (is_sarray_unique_int) {
//            arr[sarray_unique_int_size] = arr[i];
//            sarray_unique_int_size++;
//        }
//    }
//    return sarray_unique_int_size;
//}




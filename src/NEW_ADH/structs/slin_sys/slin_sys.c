/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  slin_sys.c This file collects methods of the SLIN_SYS structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes an AdH Linear System
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dmod           (SDMODEL **)  a pointer to an AdH design-level model
 * @param[in]  nSuperModels            (int) the total number of supermodels in the design model
 * @param[in]  nSubModels                (int*) the total number of submodels in each supermodel
 * @param[in]  nFluxInterfaces      (int*) the total number of flux interfaces between supermodels
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void slin_sys_alloc_array(SLIN_SYS **lin_sys, int nlin_sys){
	// allocate
    (*lin_sys) = (SLIN_SYS *) tl_alloc(sizeof(SLIN_SYS), nlin_sys);

    //some loop that allocates/initializes stuff in lin_sys
    slin_sys_init_array(*lin_sys, nlin_sys);

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes an AdH Linear System
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dmod           (SDMODEL **)  a pointer to an AdH design-level model
 * @param[in]  nSuperModels            (int) the total number of supermodels in the design model
 * @param[in]  nSubModels                (int*) the total number of submodels in each supermodel
 * @param[in]  nFluxInterfaces      (int*) the total number of flux interfaces between supermodels
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void slin_sys_init_array(SLIN_SYS *lin_sys, int nlin_sys){
    //set some loop that sets default values
    int isys;
    for(isys=0;isys<nlin_sys;isys++){
        lin_sys[isys].nnz_diag =0;
        lin_sys[isys].nnz_off_diag=0;
        lin_sys[isys].nnz_diag_old=0;
        lin_sys[isys].nnz_off_diag_old=0;
        lin_sys[isys].local_range[0]=0;
        lin_sys[isys].local_range[1]=0;
        lin_sys[isys].local_range_old[0]=0;
        lin_sys[isys].local_range_old[1]=0;
    }


}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Prints CSR format matrix to screen
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] Ax (double*) - the matrix-vector product Ax
 *  @param[in] indptr (int*) - the first/last index of each row in CSR format
 *  @param[in] cols (int*) - column locations of each nonzero in CSR format
 *  @param[in] vals (double*) - nonzero values in CSR matrix
 *  @param[in] nrows (int) - the number of rows on this process (same as # owned d.o.fs)
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void lin_sys_CSR_printScreen(SLIN_SYS *lin_sys){
	int *indptr = lin_sys->indptr_diag;
	int *cols = lin_sys->cols_diag;
	double *vals = lin_sys->vals_diag;
	int nrow = *(lin_sys->local_size);

    int i,j,nentry,row_start,row_end;
    for(i=0;i<nrow;i++){
        printf("CSR Matrix Row %d:",i);
        printf("(col,val): ");
        row_start = indptr[i];
        row_end = indptr[i+1];
        nentry = row_end-row_start;
        for(j=0;j<nentry;j++){
            printf(" (%d,%.17e)",cols[row_start+j],vals[row_start+j]);
        }
        printf("\n");
    }

}

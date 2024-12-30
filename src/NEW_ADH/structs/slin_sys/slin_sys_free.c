#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees and AdH grid
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] pgrid          (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void slin_sys_free_array(SLIN_SYS *lin_sys, int nlin_sys) {
    int i;
    for (i=0; i<nlin_sys; i++) {slin_sys_free(&(lin_sys[i]));}
    lin_sys = (SLIN_SYS *) tl_free(sizeof(SLIN_SYS), nlin_sys, lin_sys);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees and AdH grid
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] pgrid          (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void slin_sys_free(SLIN_SYS *lin_sys) {
#ifdef _PETSC
    MatDestroy(&(lin_sys->A));
    //PetscFree(&(lin_sys->A));
    KSPDestroy(&(lin_sys->ksp));
    VecDestroy(&(lin_sys->B)); //Petsc residual
    VecDestroy(&(lin_sys->X)); //Petsc residual //Persc solution
#endif
    //Free all pointers in Split CSR format
    if(lin_sys->indptr_diag!=NULL){
        lin_sys->indptr_diag= (int *) tl_free(sizeof(int), *(lin_sys->local_size)+1, lin_sys->indptr_diag);
    }
    if(lin_sys->cols_diag!=NULL){
        lin_sys->cols_diag= (int *) tl_free(sizeof(int), lin_sys->nnz_diag, lin_sys->cols_diag);
    }
    if(lin_sys->vals_diag!=NULL){
        lin_sys->vals_diag= (double *) tl_free(sizeof(double), lin_sys->nnz_diag, lin_sys->vals_diag);
    }
    if(lin_sys->indptr_off_diag!=NULL){
        lin_sys->indptr_off_diag= (int *) tl_free(sizeof(int), *(lin_sys->local_size)+1, lin_sys->indptr_off_diag);
    }
    if(lin_sys->cols_off_diag!=NULL){
        lin_sys->cols_off_diag= (int *) tl_free(sizeof(int), lin_sys->nnz_off_diag, lin_sys->cols_off_diag);
    }
    if(lin_sys->vals_off_diag!=NULL){
        lin_sys->vals_off_diag= (double *) tl_free(sizeof(double), lin_sys->nnz_off_diag, lin_sys->vals_off_diag);
    }
    if(lin_sys->residual!=NULL){
        lin_sys->residual= (double *) tl_free(sizeof(double), *(lin_sys->size), lin_sys->residual);
    }
    if(lin_sys->scale_vect!=NULL){
        lin_sys->scale_vect= (double *) tl_free(sizeof(double), *(lin_sys->size), lin_sys->scale_vect);
    }
    if(lin_sys->dsol!=NULL){
        lin_sys->dsol= (double *) tl_free(sizeof(double), *(lin_sys->size), lin_sys->dsol);
    }
    if(lin_sys->ghosts!=NULL){
        lin_sys->ghosts= (int *) tl_free(sizeof(int), lin_sys->nghost, lin_sys->ghosts);
    }
}


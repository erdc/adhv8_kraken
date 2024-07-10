 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file initalizes a SuperModel Newton Jacobian matrix
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_allocate_initialize_linear_system(SSUPER_MODEL *sm) {

    if (sm->ndofs > sm->ndofs_old){    
        #ifdef _PETSC
            allocate_petsc_system(sm);
        #else
            allocate_adh_system(sm);
        #endif
    }
}



 void allocate_adh_system(SSUPER_MODEL *sm){   
    // For now, we pre-allocate assuming that all nodes have the max number of equations attached.
    // This is not optimal, although the AdH storage is dumped to CCS and all extra zeros are removed before actually solving.
    // ndofs = nnodes * max_nsys (over-allocated for mixed dof systems)
    sm->residual = (double *) tl_realloc(sizeof(double), sm->ndofs, sm->ndofs_old, sm->residual);
    sm->sol =      (double *) tl_realloc(sizeof(double), sm->ndofs, sm->ndofs_old, sm->sol);
    sarray_init_dbl(sm->residual, sm->ndofs);
    sarray_init_dbl(sm->sol, sm->ndofs);
    
    // proprietary AdH matrix allocation
    int inode;
    sm->diagonal = (double *) tl_realloc(sizeof(double), sm->ndofs, sm->ndofs_old, sm->diagonal);
    sm->matrix = (SPARSE_VECT *) tl_realloc(sizeof(SPARSE_VECT), sm->ndofs, sm->ndofs_old, sm->matrix);
    for (inode = sm->ndofs_old; inode < sm->ndofs; inode++) {
            sm->matrix[inode].value           = (double *) tl_alloc(sizeof(double), SPV_BLOCK);
            sm->matrix[inode].index           = (int *)    tl_alloc(sizeof(int),    SPV_BLOCK); //this used to nodal, is this what we want?
            sm->matrix[inode].max_size        = SPV_BLOCK;
            sm->matrix[inode].size            = 0; //what is this?
        }
    
    init_adh_matrix(sm->ndofs, sm->matrix, sm->diagonal)
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Initializes the AdH Newton Jacobi matrix.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] matrix the sparse AdH matrix
 *  @param[in,out] diagonal the sparse AdH matrix diagonal blocks
 *  @param[in,out] nnodes the total number of grid nodes
 *  @param[in,out] max_nsys_sq the maximum number of equations in system squared
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void init_adh_matrix(int ndofs, SPARSE_VECT *matrix, double *diagonal) {
    
    int i, iend;

    for (i = 0, iend = ndofs; i < iend; i++) {
        diagonal[i] = 0.0;
    }
    
    for (i = 0; i < nnodes; i++) {
        spv_init(matrix[i], ndofs*ndofs);
    }
}

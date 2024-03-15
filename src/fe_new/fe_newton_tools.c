/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file fe_newton_tools.c This file  contains all routines for initialization, allocation, etc of the Jacobian and Residual */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
static int FILE_DEBUG = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file initalizes a SuperModel Newton residual
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct
 *  \returns ierr = error flag for initializing residual
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int fe_initialize_supermodel_residual(SSUPER_MODEL *sm) {
    ierr = 0;
#ifdef _PETSC
    ierr = VecZeroEntries(sm->residual);
#else
    sarray_init_dbl(sm->residual, sm->matrix_size);
#endif
    return ierr;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file assembles a local 3 degree-of-freedom element into the global SuperModel Newton residual
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in] SGRID *grid - the grid over which the monolithic residual resides
 *  @param[in] double *fmap - a map from the specific model dof to the supermodel dof
 *  @param[in] int *GnodeIDs - an array of length nnodes with the global node #'s of the local element nodes
 *  @param[in] int nnodes - the number of local nodes on the element
 *  @param[in] DOF_3 *elem_rhs - the local element right-hand-side
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct - adjusts the residual
 * \note elem_rhs[0] = x_eq, elem_rhs[1] = y_eq, elem_rhs[2] = c_eq,
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_elem_assemble(SSUPER_MODEL *sm, SGRID *grid, double *fmap, int *GnodeIDs, int nnodes, DOF_3 *elem_rhs) {
    int i,j;
    
#ifdef _PETSC
    PetscScalar values[3];
    for (i=0; i<nnodes; i++) {
        // No explicitly programmed ghost nodes so only add values from residential nodes
        if(grid->node[GNodeIDs[i]].resident_pe == grid->smpi->myid){
            j = fmap[GNodeIDs[i]];
            if (sm->max_ndof == 1) { // cjt :: why do we still use a 3dof matrix here?
                values[0] = -elem_rhs[i].c_eq;
                values[1] = 0;
                values[2] = 0;
            } else if (sm->max_ndof == 3) {
                values[0] = -elem_rhs[i].x_eq;
                values[1] = -elem_rhs[i].y_eq;
                values[2] = -elem_rhs[i].c_eq;
            } else {
                tl_error("max_ndof should be 1 or 3\n");
            }
            VecSetValuesBlockedLocal(sm->residual,1,&j,values,ADD_VALUES);
        }
    }
#else
    for (i=0; i<nnodes; i++) {
        if (sm->max_ndof == 1) {
            j = 1 * fmap[GNodeIDs[i]];
            sm->residual[j] = -elem_rhs[i].c_eq;
        } else if (sm->max_ndof == 3) {
            j = 3 * fmap[GNodeIDs[i]];
            sm->residual[j]     -= elem_rhs[i].x_eq;
            sm->residual[j + 1] -= elem_rhs[i].y_eq;
            sm->residual[j + 2] -= elem_rhs[i].c_eq;
        }
        else {
            tl_error("max_ndof should be 1 or 3\n");
        }
    }
#endif
    
}

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
void fe_allocate_initialize_supermodel_system(SSUPER_MODEL *sm) {
#ifdef _PETSC
    allocate_petsc_system(sm);
#else
    // For now, we pre-allocate assuming that all nodes have the max number of equations attached.
    // This is not optimal, although the AdH storage is dumped to CCS and all extra zeros are removed before actually solving.
    // ndofs = nnodes * max_nsys (over-allocated for mixed dof systems)
    sm->residual = (double *) tl_realloc(sizeof(double), sm->ndofs, sm->ndofs_old, sm->residual);
    sm->sol =      (double *) tl_realloc(sizeof(double), sm->ndofs_size, sm->ndofs_old, sm->sol);
    sarray_init_dbl(sm->residual, sm->ndofs);
    sarray_init_dbl(sm->sol, sm->ndofs);
    
    // proprietary AdH matrix allocation
    int inode;
    sm->diagonal = (double *) tl_realloc(sizeof(double), sm->ndofs, sm->ndofs_old, sm->diagonal);
    sm->matrix = (SPARSE_VECT *) tl_realloc(sizeof(SPARSE_VECT), sm->total_nnodes, sm->total_nnodes_matrix, sm->matrix);
    for (inode = sm->total_nnodes_matrix; inode < sm->total_nnodes; inode++) {
            sm->matrix[inode].value           = (double *) tl_alloc(sizeof(double), SPV_BLOCK * sm->max_nsys_sq);
            sm->matrix[inode].index           = (int *)    tl_alloc(sizeof(int),    SPV_BLOCK);
            sm->matrix[inode].max_size        = SPV_BLOCK;
            sm->matrix[inode].size            = 0;
        }
    fe_initialize_supermodel_matrix(sm);
#endif
}

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
 *  \returns ierr = error flag for initializing residual
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_initialize_supermodel_matrix(SSUPER_MODEL *sm) {
    ierr = 0;
#ifdef _PETSC
    ierr = MatZeroEntries(sm->A);
#else
    init_adh_matrix(sm->grid->nnodes, sm->max_nsys_sq, sm->matrix, sm->diagonal);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

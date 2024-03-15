/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This routine PETSC allocates a linear system
 *  \author    Sam Estes, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout]
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void allocate_petsc_system(SSUPER_MODEL *sm) {
    
    SGRID *sm->grid;
#ifdef _PETSC
    
    // Determine if the grid was updated (refined or unrefined) on any pe
#ifdef _MESSG
    int REF_FLAG_RECV=0, UNREF_FLAG_RECV=0;
    int REF_FLAG = grid->GRID_REFINED;
    int UNREF_FLAG = grid->GRID_UNREFINED;
    MPI_Allreduce(&REF_FLAG, &REF_FLAG_RECV, 1, MPI_INT, MPI_MAX, grid->smpi->ADH_COMM);
    MPI_Allreduce(&UNREF_FLAG, &UNREF_FLAG_RECV, 1, MPI_INT, MPI_MAX, grid->smpi->ADH_COMM);
#else
    REF_FLAG_RECV = grid->GRID_REFINED;
    UNREF_FLAG_RECV = grid->GRID_UNREFINED;
#endif
    
    // if the grid has been un/refined or the grid is unallocated, allocate it
    if (REF_FLAG_RECV==YES || UNREF_FLAG_RECV==YES || sm->A == PETSC_NULL){
        
        // total residential degrees of freedom over all physics in the monolithic solve
        int n = sm->my_ndofs;
        
        // Check if Jacobian, sol, and residual have already been created.
        // If so, destroy each of them before creating new PETSc objects.
        if(sm->A != PETSC_NULL){
            //printf("\n\nDestroying old matrix\n\n");
            ierr = MatDestroy(&(sm->A));
            ierr = KSPReset(sm->ksp);
            ierr = KSPSetFromOptions(sm->ksp);
        }
        if(sm->P != PETSC_NULL){
            ierr = MatDestroy(&(sm->P));
        }
        
        if(sm->sol != PETSC_NULL){
            ierr = VecDestroy(&(sm->sol));
        }
        if(sm->residual != PETSC_NULL){
            ierr = VecDestroy(&(sm->residual));
        }
        
        // Create Jacobian matrix
        ierr = MatCreate(PETSC_COMM_WORLD, &(sm->A));
        ierr = MatSetSizes(sm->A,n,n,PETSC_DETERMINE,PETSC_DETERMINE);
        ierr = MatSetFromOptions(sm->A);
        ierr = MatSetBlockSize(sm->A, sm->max_nsys); // CJT: COULD BE AN ISSUE
        ierr = MatSetUp(sm->A);
        
        // Create preallocator matrix
        if(sm->PREALLOC_FLAG == ON){
            ierr = MatCreate(PETSC_COMM_WORLD, &(sm->P));
            ierr = MatSetType(sm->P, MATPREALLOCATOR);
            ierr = MatSetSizes(sm->P,n,n,PETSC_DETERMINE,PETSC_DETERMINE);
            ierr = MatSetBlockSize(sm->P, sm->max_nsys); // CJT: COULD BE AN ISSUE
            ierr = MatSetUp(sm->P);
        }
        
        // Get ownership ranges of each proc
        ierr = MatGetOwnershipRange(sm->A, &(sm->Istart), &(sm->Iend));
        ierr = MatGetOwnershipRanges(sm->A, &(sm->ownership_range));
        
        // Create the Newtonian residual array
        ierr = VecCreate(PETSC_COMM_WORLD, &(sm->residual));
        ierr = VecSetSizes(sm->residual, n, PETSC_DETERMINE);
        ierr = VecSetFromOptions(sm->residual);
        ierr = VecSetBlockSize(sm->residual, sm->max_nsys); // CJT: COULD BE AN ISSUE
        ierr = VecSetUp(sm->residual);
        
        // Create the Newton solution (array of solution increments)
        ierr = VecCreate(PETSC_COMM_WORLD, &(sm->sol));
        ierr = VecSetSizes(sm->sol, n, PETSC_DETERMINE);
        ierr = VecSetFromOptions(sm->sol);
        ierr = VecSetUp(sm->sol);
        
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // CJT: NOT SURE HERE AND BELOW IS GOING TO WORK WITH PHYSICS COUPLING
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        // Local to global mapping object for convenient indexing when
        // setting values in residual/load/assemble functions
        ISLocalToGlobalMapping ltog_map;
        
        // TODO: Change the local to global mapping to work for
        // multiple models - SAM
        // Create local to global mapping for PETSc objects
        int petsc_global_indices[sm->nnodes];
        
        // Loop over the residential nodes first and add offset
        for(inode = 0; inode < sm->my_nnodes; inode++){
            petsc_global_indices[inode] = inode + (sm->Istart/sm->max_nsys);
        }
        // Loop over the ghost nodes and map them to the PETSc
        // index using information from the residential processors
        for(inode = sm->my_nnodes; inode < sm->nnodes; inode++){
            // Get the PETSc-global index of the non-residential node
            resident_id = grid->node[inode].resident_id;
            resident_pe = grid->node[inode].resident_pe;
            
            petsc_global_indices[inode] =
            resident_id + (sm->ownership_range[resident_pe]/sm->max_nsys);
        }
        
        // Create the local to global PETSc map
        ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, sm->max_nsys,
                                            sm->nnodes, petsc_global_indices, PETSC_COPY_VALUES, &(ltog_map));
        
        // Apply the mapping to A and residual
        ierr = MatSetLocalToGlobalMapping(sm->A, ltog_map, ltog_map);
        if(sm->PREALLOC_FLAG == ON){
            ierr = MatSetLocalToGlobalMapping(sm->P, ltog_map, ltog_map);
        }
        ierr = VecSetLocalToGlobalMapping(sm->residual, ltog_map);
        
        // Destroy mapping object now that it has been applied
        ierr = ISLocalToGlobalMappingDestroy(&(ltog_map));
        
        if(sm->PREALLOC_FLAG == ON){
            // Preallocate space for the matrix
            fe_preallocate(sm, 0);
            //ierr = MatView(sm->A,PETSC_VIEWER_STDOUT_WORLD);
            
            // Temporarily turn off allocation errors
            ierr = MatSetOption(sm->A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
        }
    }
    
#endif
}

#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function initalizes a the PETSc objects, should only be called if PETSC is
 *  active
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] sm (SMODEL_SUPER*) - pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void slin_sys_allocate_petsc_objects(SLIN_SYS *lin_sys){
    //set up anything PETSc after allocation but before load
 #ifdef _PETSC
    int ierr;   
    //PETSc solver will be stord in ksp object
    //ksp only will exist if PETSc is active
    int local_size = *(lin_sys->local_size);
    int global_size = *(lin_sys->global_size);
    //this is somehow not null 
    if(lin_sys->ksp == PETSC_NULLPTR){
        ierr = KSPCreate(PETSC_COMM_WORLD, &(lin_sys->ksp));
        ierr = KSPSetFromOptions(lin_sys->ksp);
    }

    // Check if Jacobian, sol, and residual have already been created.
    // If so, destroy each of them before creating new PETSc objects.
    //also somehow not null
    if(lin_sys->A != PETSC_NULLPTR){
        ierr = MatDestroy(&(lin_sys->A));
        ierr = KSPReset(lin_sys->ksp);
        ierr = KSPSetFromOptions(lin_sys->ksp);
    }
    //dont think we need preallocator matrix
//    if(sm->P != PETSC_NULLPTR){
//        ierr = MatDestroy(&(sm->P));
//    }
    if(lin_sys->X != PETSC_NULLPTR){
        ierr = VecDestroy(&(lin_sys->X));
    }
    if(lin_sys->B != PETSC_NULLPTR){
        ierr = VecDestroy(&(lin_sys->B));
    }
    // Create Jacobian matrix from CSR format
    //don't know if i need to do this or just call after load anyway
    //vals should update by reference
    if (lin_sys->indptr_off_diag == NULL){
        MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, local_size, local_size, lin_sys->indptr_diag, lin_sys->cols_diag, lin_sys->vals_diag, &(lin_sys->A));
        printf("PETSC Mat created with sequential array\n");
    }else{
        MatCreateMPIAIJWithSplitArrays(PETSC_COMM_WORLD, local_size, local_size, global_size, global_size, lin_sys->indptr_diag, lin_sys->cols_diag, lin_sys->vals_diag, lin_sys->indptr_off_diag, lin_sys->cols_off_diag, lin_sys->vals_off_diag, &(lin_sys->A));
        printf("PETSC Mat created with split arrays\n");
    }
    
    //ierr = MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, sm->my_ndofs, sm->my_ndofs, sm->macro_ndofs, sm->macro_ndofs, sm->indptr, sm->indices, sm->vals, &(sm->A))
    //ierr = MatCreate(PETSC_COMM_WORLD, &(sm->A));
    //ierr = MatSetSizes(sm->A,n,n,PETSC_DETERMINE,PETSC_DETERMINE);
    //ierr = MatSetFromOptions(sm->A);
    //ierr = MatSetBlockSize(sm->A, sm->max_nsys);
    //ierr = MatSetUp(sm->A);

    // Create preallocator matrix
    //dont need i dont think
    /*
    if(sm->PREALLOC_FLAG == ON){
        ierr = MatCreate(PETSC_COMM_WORLD, &(sm->P));
        ierr = MatSetType(sm->P, MATPREALLOCATOR);
        ierr = MatSetSizes(sm->P,n,n,PETSC_DETERMINE,PETSC_DETERMINE);
        ierr = MatSetBlockSize(sm->P, sm->max_nsys);
        ierr = MatSetUp(sm->P);
    }
    */
    //// Currently hard-coded for a single model - SAM
    //mod = &(sm->submodel[0]);

    // Get ownership ranges of each proc
    //dont think I need this either
    //ierr = MatGetOwnershipRange(sm->A, &(sm->Istart), &(sm->Iend));
    //ierr = MatGetOwnershipRanges(sm->A, &(sm->ownership_range));
        
    // Create the Newtonian residual array
    //lets see what we actually need to use of this
    //ierr = VecCreate(PETSC_COMM_WORLD, &(sm->B));
    //ierr = VecSetSizes(sm->B, sm->my_ndofs, PETSC_DETERMINE);
    //ierr = VecSetFromOptions(sm->B);
    //ierr = VecSetBlockSize(sm->residual, sm->max_nsys);
    //ierr = VecSetUp(sm->residual);
        
    // Create the Newton solution (array of solution increments)
    //ierr = VecCreate(PETSC_COMM_WORLD, &(sm->X));
    //ierr = VecSetSizes(sm->X, sm->my_ndofs, PETSC_DETERMINE);
    //ierr = VecSetFromOptions(sm->sol);
    //ierr = VecSetUp(sm->sol);
    ierr = VecCreateGhostWithArray(PETSC_COMM_WORLD, local_size, global_size, lin_sys->nghost, lin_sys->ghosts, lin_sys->residual, &(lin_sys->B));
    ierr = VecCreateGhostWithArray(PETSC_COMM_WORLD, local_size, global_size, lin_sys->nghost, lin_sys->ghosts, lin_sys->dsol, &(lin_sys->X));


    /*Maybe need to revisit for ghosts but skip for now

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
        resident_id = mod->grid->node[inode].resident_id;
        resident_pe = mod->grid->node[inode].resident_pe;
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
        
    // Below needs to be converted into PETSc vectors
    // TODO: Ask Corey about some of these
    //sm->solver_info.node_block = (int *) tl_realloc(sizeof(int), sm->total_nnodes, sm->total_nnodes_matrix, sm->solver_info.node_block);
    //sm->bc_mask = (int *) tl_realloc(sizeof(int), sm->total_nnodes * sm->max_nsys, sm->old_bc_mask_size, sm->bc_mask); // CJT :: prob not needed!
    sm->bc_mask = (int *) tl_realloc(sizeof(int), sm->total_nnodes * sm->max_nsys, sm->old_bc_mask_size, sm->bc_mask); // CJT :: prob not needed!
    sm->old_bc_mask_size = sm->total_nnodes * sm->max_nsys;
//#ifdef _MESSG
//        for (inode = sm->total_nnodes_matrix; inode < sm->total_nnodes; inode++) {
//            sm->solver_info.node_block[inode] = sm->supersmpi->myid;
//        }
//#endif
    // No longer use this for adaption condition above.
    // Could get rid of this line.
    sm->total_nnodes_matrix = sm->total_my_nnodes;
      */
#endif
}

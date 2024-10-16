 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function initalizes a SuperModel Newton Jacobian matrix
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_allocate_initialize_linear_system(SSUPER_MODEL *sm) {

    if (sm->ndofs > sm->ndofs_old){    
        allocate_adh_system(sm);
#ifdef _PETSC
        allocate_petsc_objects(sm);
#endif
            }
}




/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function initalizes a SuperModel Newton Jacobian matrix for suitesparse
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
 void allocate_adh_system(SSUPER_MODEL *sm){   
    // For now, we pre-allocate assuming that all nodes have the max number of equations attached.
    // This is not optimal, although the AdH storage is dumped to CCS and all extra zeros are removed before actually solving.
    // ndofs = nnodes * max_nsys (over-allocated for mixed dof systems)
    sm->residual = (double *) tl_realloc(sizeof(double), sm->ndofs, sm->ndofs_old, sm->residual);
    sm->sol =      (double *) tl_realloc(sizeof(double), sm->ndofs, sm->ndofs_old, sm->sol);
    sarray_init_dbl(sm->residual, sm->ndofs);
    sarray_init_dbl(sm->sol, sm->ndofs);
    
    // proprietary AdH matrix allocation
    // standard CSR format
    // indptr points to where column indeces for each row starts (length of nrows+1)
    // indices array of integers that contains column position of each nonzero
    // vals array of doubles containing matrix entries (same length as indices)
    
    //call routine to determine nnz? and where should we determine cols?
    //should we also keep nnz_row array for petsc? maybe not necessary

    //where does local_range get determined?
    //we now number of rows is sum(my_nnode*nodal_nvars)
    sm->diag_indptr= (int *) tl_realloc(sizeof(int), sm->local_range[1]-sm->local_range[0]+1, sm->local_range_old[1]-sm->local_range_old[0]+1, sm->diag_indptr);
    //if we split we have two structures
    sm->off_diag_indptr = (int *) tl_realloc(sizeof(int), sm->local_range[1]-sm->local_range[0]+1, sm->local_range_old[1]-sm->local_range_old[0]+1, sm->off_diag_indptr);
    
    //find nnz and store in indices
    //this does a single CSR format
    //create_sparsity_CSR(sm, grid);
    //this does a split CSR format, plan on doing this
    create_sparsity_split_CSR(sm, grid);
    //do we want to stor nnz? it is stored in sm->indptr[nrows]
    //do we want to store local_size = local_range[1]-local_range[0]
    sarray_init_dbl(sm->diag_vals, sm->diag_nnz);
    sarray_init_dbl(sm->off_diag_vals, sm->off_diag_nnz);
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function initalizes a SuperModel Newton Jacobian matrix for suitesparse
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void create_sparsity_split_CSR(sm, grid){
    //
    // We know that each row contains at least one non-zero entry, and the number of rows is known upfront. 
    // Thus, (i,j)-pairs can be represented as a vector of vectors
    //first allocate this by guessing a max size of nnz per row

    //guesses at max number of adjacent nodes in each dimension
    //multiply by max dof per node
    int nCon1d = 2*MAX_NVAR;
    int nCon2d = 7*3*MAX_NVAR;
    int nCon3d = 16*4*MAX_NVAR;

    int nrows = sm->local_range[1]-sm->local_range[0]

    // preallocate upfront with these dimensions. work out syntax later
    //my best idea for now, maybe way to avoid these but cant think of another way for now

    //temporarily stores column positions of diagonal and off diagonal blocks
    //diagonal blocks are local column numbers, off-diagonal are global
    int temp_cols_diag[nrows,nCon3d];
    int temp_cols_off_diag[nrows,nCon3d];
    int nnz_rows_diag[nrows];
    int nnz_rows_off_diag[nrows];
    sarray_init_dbl(nnz_rows_diag, nrows);
    sarray_init_dbl(nnz_rows_off_diag, nrows);

    //loop through all nelem3d to mark sparsity
    int row;
    int max_elem_dofs = MAX_NVAR*MAX_NNODE; 
    int dofs[max_elem_dofs],global_dofs[max_elem_dofs];
    int nnodes;
    int nvars_elem;
    int ndofs_ele;

    //loop thorugh each element and find sparsity
    for (j=0;j<grid->nelem3d;j++){
        nnodes = grid->elem3d[j].nnodes;
        //maybe this just comes from mat instead
        nvars_elem = sm->elem3d_nvars[j];
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem3d[j].nodes,nvars_elem,sm->elem3d_vars[j],sm->nodal_nvars, sm->nodal_vars);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    //differentiate  between diag and off diag blocks
                    if (row<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        temp_cols_diag[row,nnz_rows_diag[k]]=dofs[l];
                        nnz_rows_diag[k]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        temp_cols_off_diag[row,nnz_rows_off_diag[k]]=global_dofs[l];
                        nnz_rows_off_diag[k]+=1;

                    }
                }
            }

        }
    }

    //do 2d elements and then 1d
    for (j=0;j<grid->nelem2d;j++){
        nnodes = grid->elem2d[j].nnodes;
        //maybe this just comes from mat instead
        nvars_elem = sm->elem2d_nvars[j];
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem2d[j].nodes,nvars_elem,sm->elem2d_vars[j],sm->nodal_nvars, sm->nodal_vars);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    //differentiate  between diag and off diag blocks
                    if (row<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        temp_cols_diag[row,nnz_rows_diag[k]]=dofs[l];
                        nnz_rows_diag[k]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        temp_cols_off_diag[row,nnz_rows_off_diag[k]]=global_dofs[l];
                        nnz_rows_off_diag[k]+=1;

                    }
                }
            }

        }
    }
    //now 1d
    for (j=0;j<grid->nelem1d;j++){
        nnodes = grid->elem1d[j].nnodes;
        //maybe this just comes from mat instead
        nvars_elem = sm->elem1d_nvars[j];
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem1d[j].nodes,nvars_elem,sm->elem1d_vars[j],sm->nodal_nvars, sm->nodal_vars);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    //differentiate  between diag and off diag blocks
                    if (row<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        temp_cols_diag[row,nnz_rows_diag[k]]=dofs[l];
                        nnz_rows_diag[k]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        temp_cols_off_diag[row,nnz_rows_off_diag[k]]=global_dofs[l];
                        nnz_rows_off_diag[k]+=1;

                    }
                }
            }

        }
    }

    //now that all nnz_cols* and nnz* are filled, we need to sort and remove duplicates
    // counter for the total number of non-zero entries
    int NNZ_diag = 0, NNZ_off_diag=0;
    int *nnz_row_diag,*nnz_row_off_diag;
    //note nnz-rows is not actually nnz per row, includes duplicates
    //could be done with hybrid omp too?

    //do for diagonal and off diagonal blocks
    for (i=0;i<nrows;i++){
        // sort the column indices (j-entries)
        //use stdlib.h qsort
        qsort(temp_cols_diag[i], nnz_rows_diag[i], sizeof(int), compare_ints);
        qsort(temp_cols_off_diag[i], nnz_rows_off_diag[i], sizeof(int), compare_ints);
        //this should hopefully remove duplicates?
        nnz_row_diag = unique(temp_cols_diag[i], temp_cols_diag[i] + (sizeof(temp_cols_diag[i])/ sizeof(int)));
        nnz_row_off_diag = unique(temp_cols_off_diag[i], temp_cols_off_diag[i] + (sizeof(temp_cols_off_diag[i])/ sizeof(int)));
        //add nnz in a row to the NNZ
        //maybe overwrire nnz_rows[i] instead if we want this stored, and then sum it
        NNZ_diag+=*nnz_row_diag;
        NNZ_off_diag+=*nnz_row_off_diag;
    }
    sm->diag_nnz = NNZ; sm->off_diag_nnz= NNZ_off_diag;
    
    // resize as needed
    sm->diag_indices = (int *) tl_realloc(sizeof(int), sm->nnz_diag, sm->nnz_diag_old, sm->diag_indices);
    sm->diag_vals = (double *) tl_realloc(sizeof(double), sm->nnz_diag, sm->nnz_diag_old, sm->diag_vals);
    sm->global_off_diag_indices = (int *) tl_realloc(sizeof(int), sm->nnz_off_diag, sm->nnz_off_diag_old, sm->global_off_diag_indices);
    sm->off_diag_vals = (double *) tl_realloc(sizeof(double), sm->nnz_off_diag, sm->nnz_off_diag_old, sm->off_diag_vals);
    
    // enumerate the sorted indices to populate COL_INDEX and ROW_INDEX for diag and off diag
    int count_diag=0;
    int count_off_diag = 0;
    for(i=0;i<nrows;i++){
        sm->diag_indptr[i] = count_diag;
        sm->off_diag_indptr[i] = count_off_diag;
        //ask corey about syntax
        &rn_diag = *temp_cols_diag[i];
        &rn_off_diag = *temp_cols_diag[i];
        //think about how to do this since each temp_rows may be different size
        for(j=0;j<sizeof(rn_diag)/sizeof(int),j++){
            sm->diag_indices[count] = rn_diag[j];
            count_diag++;
        }
        for(j=0;j<sizeof(rn_off_diag)/sizeof(int),j++){
            sm->global_off_diag_indices[count] = rn_off_diag[j];
            count_off_diag++;
        }
    }
    //also last sm->indptr  needs to be the last entry
    sm->diag_indptr[nrows] = count_diag;
    sm->off_diag_indptr[nrows] = count_off_diag;
    //destroy temp_rows?
    free(temp_cols_diag);
    free(temp_cols_off_diag);
    free(nnz_rows);
    free(nnz_rows_diag);
}


int* unique (int* first, int* last)
{
  if (first==last) return last;

  int* result = first;
  while (++first != last)
  {
    if (!(*result == *first)) 
      *(++result)=*first;
  }
  return ++result;
}

int compare_ints(const void *a, const void *b) {
    return (*(int*)a - *(int*)b); // Ascending order
}




/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function initalizes a SuperModel Newton Jacobian matrix for suitesparse
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void allocate_petsc_objects(SSUPER_MODEL *sm){
    int ierr;
    //set up anything PETSc after allocation but before load

    
    //PETSc solver will be stord in ksp object
    //ksp only will exist if PETSc is active
    if(sm->ksp == PETSC_NULL){
        ierr = KSPCreate(PETSC_COMM_WORLD, &(sm->ksp));
        ierr = KSPSetFromOptions(sm->ksp);
    }
    // Check if Jacobian, sol, and residual have already been created.
    // If so, destroy each of them before creating new PETSc objects.
    if(sm->A != PETSC_NULL){
        printf("\n\nDestroying old matrix\n\n");
        ierr = MatDestroy(&(sm->A));
        ierr = KSPReset(sm->ksp);
        ierr = KSPSetFromOptions(sm->ksp);
    }
    //dont think we need preallocator matrix
    //if(sm->P != PETSC_NULL){
    //    ierr = MatDestroy(&(sm->P));
    //}
    if(sm->X != PETSC_NULL){
        ierr = VecDestroy(&(sm->X));
    }
    if(sm->B != PETSC_NULL){
            ierr = VecDestroy(&(sm->B));
    }

    // Create Jacobian matrix from CSR format
    //don't know if i need to do this or just call after load anyway
    //vals should update by reference
    MatCreateMPIAIJWithSplitArrays(PETSC_COMM_WORLD, sm->my_ndofs, sm->my_ndofs, sm->macro_ndofs, sm->macro_ndofs, sm->diag_indptr, sm->diag_indices, sm->diag_vals, sm->off_diag_indptr, sm->global_off_diag_indices, sm->off_diag_vals, &(sm->A));
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
    ierr = VecCreate(PETSC_COMM_WORLD, &(sm->B));
    ierr = VecSetSizes(sm->B, sm->my_ndofs, PETSC_DETERMINE);
    ierr = VecSetFromOptions(sm->B);
    //ierr = VecSetBlockSize(sm->residual, sm->max_nsys);
    //ierr = VecSetUp(sm->residual);
        
    // Create the Newton solution (array of solution increments)
    ierr = VecCreate(PETSC_COMM_WORLD, &(sm->X));
    ierr = VecSetSizes(sm->X, sm->my_ndofs, PETSC_DETERMINE);
    //ierr = VecSetFromOptions(sm->sol);
    //ierr = VecSetUp(sm->sol);


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
}
  

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function initalizes a SuperModel Newton Jacobian matrix for suitesparse
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void create_sparsity_CSR(sm, grid){
    //
    // We know that each row contains at least one non-zero entry, and the number of rows is known upfront. 
    // Thus, (i,j)-pairs can be represented as a vector of vectors
    //first allocate this by guessing a max size of nnz per row

    //guesses at max number of adjacent nodes in each dimension
    //multiply by max dof per node
    int nCon1d = 2*MAX_NVAR;
    int nCon2d = 7*3*MAX_NVAR;
    int nCon3d = 16*4*MAX_NVAR;

    int nrows = sm->local_range[1]-sm->local_range[0]

    // preallocate upfront with these dimensions. work out syntax later
    int temp_rows[nrows,nCon3d];
    int nnz_rows[nrows];
    sarray_init_dbl(nnz_rows, nrows);


    //loop through all nelem3d to mark sparsity

    int global_row;
    //loop thorugh each element and find sparsity
    for (j=0;j<grid->nelem3d;j++){
        nnodes = grid->elem3d[j].nnodes;
        //maybe this just comes from mat instead
        nvars_elem = sm->elem3d_nvars[j];
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem3d[j].nodes,nvars_elem,sm->elem3d_vars[j],sm->nodal_nvars, sm->nodal_vars);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        //keep sparsity pattern in temp_rows
        for (k=0;k<ndofs_ele;k++){
            global_row = global_dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
        if(global_row>=local_range[0] && global_row < local_range[1])
            //get row local to process
            row = global_row-local_range[0];
            for (l=0;l<ndofs_ele;l++){
                temp_rows[row,nnz_rows[k]]=global_dof[l];
                nnz_rows[k]+=1;
            }
        }

    }

    //do 2d elements and then 1d
    for (j=0;j<grid->nelem2d;j++){
        nnodes = grid->elem2d[j].nnodes;
        //maybe this just comes from mat instead
        nvars_elem = sm->elem2d_nvars[j];
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem2d[j].nodes,nvars_elem,sm->elem2d_vars[j],sm->nodal_nvars, sm->nodal_vars);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        //keep sparsity pattern in temp_rows
        for (k=0;k<ndofs_ele;k++){
            global_row = global_dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
        if(global_row>=local_range[0] && global_row < local_range[1])
            //get row local to process
            row = global_row-local_range[0];
            for (l=0;l<ndofs_ele;l++){
                temp_rows[row,nnz_rows[k]]=global_dof[l];
                nnz_rows[k]+=1;
            }
        }

    }

    for (j=0;j<grid->nelem1d;j++){
        nnodes = grid->elem1d[j].nnodes;
        //maybe this just comes from mat instead
        nvars_elem = sm->elem1d_nvars[j];
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem1d[j].nodes,nvars_elem,sm->elem1d_vars[j],sm->nodal_nvars, sm->nodal_vars);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        //keep sparsity pattern in temp_rows
        for (k=0;k<ndofs_ele;k++){
            global_row = global_dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
        if(global_row>=local_range[0] && global_row < local_range[1])
            //get row local to process
            row = global_row-local_range[0];
            for (l=0;l<ndofs_ele;l++){
                temp_rows[row,nnz_rows[k]]=global_dof[l];
                nnz_rows[k]+=1;
            }
        }

    }

    //now that all nnz_rows and temp_rows are filled, we need to sort and remove duplicates
    // counter for the total number of non-zero entries
    int NNZ = 0;
    int *nnz_row;
    //note nnz-rows is not actually nnz per row, includes duplicates
    //could be done with hybrid omp too?
    for (i=0;i<nrows;i++){
        // sort the column indices (j-entries)
        //use stdlib.h qsort
        qsort(temp_rows[i], nnz_rows[i], sizeof(int), compare_ints);
        //this should hopefully remove duplicates?
        nnz_row = unique(temp_rows[i], temp_rows[i] + (sizeof(temp_rows[i])/ sizeof(int)));
        //add nnz in a row to the NNZ
        //maybe overwrire nnz_rows[i] instead if we want this stored, and then sum it
        NNZ+=*nnz_row;
    }
    sm->nnz = NNZ;
    // resize as needed
    sm->indices = (int *) tl_realloc(sizeof(int), sm->nnz, sm->nnz_old, sm->indices);
    sm->vals = (double *) tl_realloc(sizeof(double), sm->nnz, sm->nnz_old, sm->vals);
    // enumerate the sorted indices to populate COL_INDEX and ROW_INDEX
    for(i=0;i<nrows;i++){
        sm->indptr[i] = count;
        //ask corey about syntax
        &rn = *temp_rows[i];
        //think about how to do this since each temp_rows may be different size
        for(j=0;j<sizeof(rn)/sizeof(int),j++){
            sm->indices[count] = rn[j];
            count++;
        }
    }
    //also last sm->indptr  needs to be the last entry
    sm->indptr[nrows] = count;
    //destroy temp_rows?
}




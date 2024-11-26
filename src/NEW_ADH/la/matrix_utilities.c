#include "adh.h" 
static int **temp_cols_diag = NULL;
static int **temp_cols_off_diag = NULL;
static int isize = 0;  /* the size of the nnz rows arrays */
static int max_elem_dofs = MAX_NVAR*MAX_NNODE; 
static int *nnz_rows_diag = NULL;
static int *nnz_rows_off_diag = NULL;
static int *nnz_rows_diag_no_duplicate = NULL;
static int *nnz_rows_off_diag_no_duplicate = NULL;
/*+++++++++++++++++++++++++++++++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function initalizes a SuperModel Newton Jacobian matrix for suitesparse
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SMODEL_SUPER *sm pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void create_sparsity_split_CSR(SMODEL_SUPER *sm, SGRID *grid){
    //temporarily stores column positions of diagonal and off diagonal blocks
    //diagonal blocks are local column numbers, off-diagonal are global


    bool has_off_diag = false;
    int row,col;
    int nnodes;
    int nvars_elem;
    int ndofs_ele;
    int i,j,k,l, mat_id;
    int elem_vars[MAX_NVAR];
    int dofs[max_elem_dofs],global_dofs[max_elem_dofs];
    int dum1,dum2;
    int isize_prev;
    int NNZ_diag = 0;
    int nnz_row_diag;
    int NNZ_off_diag=0;
    int nnz_row_off_diag;
    int count_diag=0;
    int *rn_diag;
    int count_off_diag = 0;
    int *rn_off_diag;
    // We know that each row contains at least one non-zero entry, and the number of rows is known upfront. 
    // Thus, (i,j)-pairs can be represented as a vector of vectors
    //first allocate this by guessing a max size of nnz per row
    int* fmap = sm->dof_map_local;
    int *local_range = sm->local_range;
    int nrows = local_range[1]-local_range[0];
    //int nnz_rows_diag[nrows];

    //see if we have off diagonal blocks to handle
    if(sm->indptr_off_diag!=NULL){
        has_off_diag=true;
    }

    //reallocate if necessary
    if (isize < nrows){
        isize_prev = isize;
        isize = nrows;
        nnz_rows_diag = (int *) tl_realloc(sizeof(int), isize, isize_prev, nnz_rows_diag);
        nnz_rows_diag_no_duplicate = (int *) tl_realloc(sizeof(int), isize, isize_prev, nnz_rows_diag_no_duplicate);
        if (has_off_diag){
            nnz_rows_off_diag = (int *) tl_realloc(sizeof(int), isize, isize_prev, nnz_rows_off_diag);
            nnz_rows_off_diag_no_duplicate = (int *) tl_realloc(sizeof(int), isize, isize_prev, nnz_rows_off_diag_no_duplicate);
        }
    }

    //printf("NROWS=%d\n",nrows);
    //initialize arrays with zeros
    sarray_init_int(nnz_rows_diag, nrows);
    if (has_off_diag){
        sarray_init_int(nnz_rows_off_diag, nrows);
    }

    //First set of loops is solely to establish how many nonzeros are there 
    //and how to dynamically allocate temporary sparsity arrays

    //loop through and find how much to allocate
    for (j=0;j<grid->nelems3d;j++){
        nnodes = grid->elem3d[j].nnodes;
        //pull all global information to local memory
        mat_id = sm->elem3d_physics_mat_id[j];
        //Get stuff from physics mat
        nvars_elem = sm->elem3d_physics_mat[mat_id].nvar;
        //get array of variables specified on element
        sarray_copy_int(elem_vars, sm->elem3d_physics_mat[mat_id].vars, nvars_elem);
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem3d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //get_cell_dofs_2(dofs, nnodes, grid->elem3d[j].nodes ,nvars_elem, elem_vars, sm->node_physics_mat, sm->node_physics_mat_id);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,sm->ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    //differentiate  between diag and off diag blocks by column number (local to process)
                    col = dofs[l];
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        nnz_rows_off_diag[row]+=1;

                    }
                }
            }

        }
    }

    //do 2d elements and then 1d
    //loop thorugh each element and find sparsity
    for (j=0;j<grid->nelems2d;j++){
        nnodes = grid->elem2d[j].nnodes;
        //pull all global information to local memory
        mat_id = sm->elem2d_physics_mat_id[j];
        //Get stuff from physics mat
        nvars_elem = sm->elem2d_physics_mat[mat_id].nvar;
        //get array of variables specified on element
        sarray_copy_int(elem_vars, sm->elem2d_physics_mat[mat_id].vars, nvars_elem);
        //printf("elemental variables copied\n");
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem2d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //get_cell_dofs_2(dofs, nnodes, grid->elem2d[j].nodes ,nvars_elem, elem_vars, sm->node_physics_mat, sm->node_physics_mat_id);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,sm->ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    col = dofs[l];
                    //differentiate  between diag and off diag blocks
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        nnz_rows_off_diag[row]+=1;

                    }
                }
            }

        }
    }

    //now 1d
    //loop thorugh each element and find sparsity
    for (j=0;j<grid->nelems1d;j++){
        nnodes = grid->elem1d[j].nnodes;
        //pull all global information to local memory
        mat_id = sm->elem1d_physics_mat_id[j];
        //Get stuff from physics mat
        nvars_elem = sm->elem1d_physics_mat[mat_id].nvar;
        //get array of variables specified on element
        sarray_copy_int(elem_vars, sm->elem1d_physics_mat[mat_id].vars, nvars_elem);
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem1d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //get_cell_dofs_2(dofs, nnodes, grid->elem1d[j].nodes ,nvars_elem, elem_vars, sm->node_physics_mat, sm->node_physics_mat_id);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,sm->ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    col = dofs[l];
                    //differentiate  between diag and off diag blocks
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        nnz_rows_off_diag[row]+=1;
                    }
                }
            }

        }
    }

//    printf("Found nnz\n");
//    for (j=0;j<nrows;j++){
//        printf("NNZ rows diag[%d]: %d,\n",j,nnz_rows_diag[j]);
//    }
    //use nnz_rows to dynamically allocate
    //int temp_cols_diag[nrows][nCon3d];
    temp_cols_diag = (int**) tl_alloc(sizeof(int*), nrows);
    for(j=0;j<nrows;j++){
        temp_cols_diag[j] = (int*) tl_alloc(sizeof(int), nnz_rows_diag[j]);
        for(k=0;k<nnz_rows_diag[j];k++){
            temp_cols_diag[j][k]=INT_MAX;
        }
    }

    if (has_off_diag){
            temp_cols_off_diag = (int**) tl_alloc(sizeof(int*), nrows);
        for(j=0;j<nrows;j++){
            temp_cols_off_diag[j] = (int*) tl_alloc(sizeof(int), nnz_rows_off_diag[j]);
            for(k=0;k<nnz_rows_off_diag[j];k++){
                temp_cols_off_diag[j][k]=INT_MAX;
            }
        }   
    }

    //Seems redundant but must reuse as indexing
    sarray_init_int(nnz_rows_diag, nrows);
    if (has_off_diag){
        sarray_init_int(nnz_rows_off_diag, nrows);
    }
    //loop thorugh each element and find sparsity
    for (j=0;j<grid->nelems3d;j++){
        nnodes = grid->elem3d[j].nnodes;
        //pull all global information to local memory
        mat_id = sm->elem3d_physics_mat_id[j];
        //Get stuff from physics mat
        nvars_elem = sm->elem3d_physics_mat[mat_id].nvar;
        //get array of variables specified on element
        sarray_copy_int(elem_vars, sm->elem3d_physics_mat[mat_id].vars, nvars_elem);
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem3d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //get_cell_dofs_2(dofs, nnodes, grid->elem3d[j].nodes ,nvars_elem, elem_vars, sm->node_physics_mat, sm->node_physics_mat_id);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,sm->ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    //differentiate  between diag and off diag blocks by column number (local to process)
                    col = dofs[l];
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        temp_cols_diag[row][nnz_rows_diag[row]]=col;
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        temp_cols_off_diag[row][nnz_rows_off_diag[row]]=global_dofs[l];
                        nnz_rows_off_diag[row]+=1;

                    }
                }
            }

        }
    }

    //do 2d elements and then 1d
    for (j=0;j<grid->nelems2d;j++){
        nnodes = grid->elem2d[j].nnodes;
        //pull all global information to local memory
        mat_id = sm->elem2d_physics_mat_id[j];
        //Get stuff from physics mat
        nvars_elem = sm->elem2d_physics_mat[mat_id].nvar;
        //get array of variables specified on element
        sarray_copy_int(elem_vars, sm->elem2d_physics_mat[mat_id].vars, nvars_elem);
        //printf("elemental variables copied\n");
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem2d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //get_cell_dofs_2(dofs, nnodes, grid->elem2d[j].nodes ,nvars_elem, elem_vars, sm->node_physics_mat, sm->node_physics_mat_id);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,sm->ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    col = dofs[l];
                    //differentiate  between diag and off diag blocks
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        temp_cols_diag[row][nnz_rows_diag[row]]=col;
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        temp_cols_off_diag[row][nnz_rows_off_diag[row]]=global_dofs[l];
                        nnz_rows_off_diag[row]+=1;
                    }
                }
            }

        }
    }

    //now 1d
    //loop thorugh each element and find sparsity
    for (j=0;j<grid->nelems1d;j++){
        nnodes = grid->elem1d[j].nnodes;
        //pull all global information to local memory
        mat_id = sm->elem1d_physics_mat_id[j];
        //Get stuff from physics mat
        nvars_elem = sm->elem1d_physics_mat[mat_id].nvar;
        //get array of variables specified on element
        sarray_copy_int(elem_vars, sm->elem1d_physics_mat[mat_id].vars, nvars_elem);
        ndofs_ele = nnodes*nvars_elem;
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem1d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //get_cell_dofs_2(dofs, nnodes, grid->elem1d[j].nodes ,nvars_elem, elem_vars, sm->node_physics_mat, sm->node_physics_mat_id);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,sm->ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    col = dofs[l];
                    //differentiate  between diag and off diag blocks
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        temp_cols_diag[row][nnz_rows_diag[row]]=col;
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        temp_cols_off_diag[row][nnz_rows_off_diag[row]]=global_dofs[l];
                        nnz_rows_off_diag[row]+=1;
                    }
                }
            }

        }
    }

//    printf("Found temp cols\n");
//    for (j=0;j<nrows;j++){
//        printf("cols diag[%d][0]: %d,\n",j,temp_cols_diag[j][0]);
//    }

    //now that all nnz_cols* and nnz* are filled, we need to sort and remove duplicates
    // counter for the total number of non-zero entries
    //note nnz-rows is not actually nnz per row, includes duplicates
    //could be done with hybrid omp too?
    //do for diagonal and off diagonal blocks
    for (i=0;i<nrows;i++){
        // sort the column indices (j-entries)
        //use stdlib.h qsort
        qsort(temp_cols_diag[i], nnz_rows_diag[i], sizeof(int), compare_ints);
        //this should hopefully remove duplicates?
        nnz_row_diag = unique(temp_cols_diag[i], nnz_rows_diag[i]);
        //overwrite nnz row with unique?
        nnz_rows_diag_no_duplicate[i] = nnz_row_diag;
        //add nnz in a row to the NNZ
        //maybe overwrire nnz_rows[i] instead if we want this stored, and then sum it
        NNZ_diag+=nnz_row_diag;
    }

    sm->nnz_diag = NNZ_diag; 

    //off diagonalblocks
    if (has_off_diag){

        for (i=0;i<nrows;i++){
            // sort the column indices (j-entries)
            //use stdlib.h qsort
            qsort(temp_cols_off_diag[i], nnz_rows_off_diag[i], sizeof(int), compare_ints);
            //this should hopefully remove duplicates?
            nnz_row_off_diag = unique(temp_cols_off_diag[i], nnz_rows_off_diag[i]);
            //overwrite nnz row with unique?
            nnz_rows_off_diag_no_duplicate[i] = nnz_row_off_diag;
            //add nnz in a row to the NNZ
            //maybe overwrire nnz_rows[i] instead if we want this stored, and then sum it
            NNZ_off_diag+=nnz_row_off_diag;
        }

        sm->nnz_off_diag = NNZ_off_diag;
    }

    // resize as needed
    sm->cols_diag = (int *) tl_realloc(sizeof(int), sm->nnz_diag, sm->nnz_diag_old, sm->cols_diag);
    sm->vals_diag = (double *) tl_realloc(sizeof(double), sm->nnz_diag, sm->nnz_diag_old, sm->vals_diag);
    if (has_off_diag){
        sm->cols_off_diag = (int *) tl_realloc(sizeof(int), sm->nnz_off_diag, sm->nnz_off_diag_old, sm->cols_off_diag);
        sm->vals_off_diag = (double *) tl_realloc(sizeof(double), sm->nnz_off_diag, sm->nnz_off_diag_old, sm->vals_off_diag);
    }

    //now use info to fill in indptr and cols
    for(i=0;i<nrows;i++){
        //printf("filling in index ptr and column entries, row %d\n",i);
        sm->indptr_diag[i] = count_diag;
        //printf("filling in index ptr and column entries, row %d\n",i);
        rn_diag = temp_cols_diag[i];
        //think about how to do this since each temp_rows may be different size
        for(j=0;j<nnz_rows_diag_no_duplicate[i];j++){
            sm->cols_diag[count_diag] = rn_diag[j];
            count_diag++;
        }
    }
    //also last sm->indptr  needs to be the last entry
    sm->indptr_diag[nrows] = count_diag;
    //assert(sm->indptr_diag[nrows] = sm->nnz_diag);

    
    if (has_off_diag){
        for(i=0;i<nrows;i++){
            //printf("filling in index ptr and column entries, row %d\n",i);
            sm->indptr_off_diag[i] = count_off_diag;
            //printf("filling in index ptr and column entries, row %d\n",i);
            rn_off_diag = temp_cols_off_diag[i];
            //think about how to do this since each temp_rows may be different size
            for(j=0;j<nnz_rows_off_diag_no_duplicate[i];j++){
                sm->cols_off_diag[count_off_diag] = rn_off_diag[j];
                count_off_diag++;
            }
        }
        sm->indptr_off_diag[nrows] = count_off_diag;
    }

    //destroy temp_rows?
    //not malloced so not necessary??
    //printf("Loop complete\n");
    //free(nnz_rows_diag);
    //free(nnz_rows_off_diag);
    //printf("Freed arrays");
    //for(i=0;i<nrows;i++){
    //    free(temp_cols_diag[i]);
    //    free(temp_cols_off_diag[i]);
    //}


    printf("About to free memory in create_sparsity_split_CSR\n");
//    for(j=0;j<nrows;j++){
//        tl_free(sizeof(int), nnz_temp_cols_diag[j]); 
//    }

    //for now free calls are in script, can separate to new routine later
    for (i=0; i<nrows; i++) {
        temp_cols_diag[i] = (int*) tl_free(sizeof(int), nnz_rows_diag[i], temp_cols_diag[i]);
    }
    temp_cols_diag = (int **) tl_free(sizeof(int *), nrows, temp_cols_diag);
    nnz_rows_diag= (int *) tl_free(sizeof(int), nrows, nnz_rows_diag);
    nnz_rows_diag_no_duplicate= (int *) tl_free(sizeof(int), nrows, nnz_rows_diag_no_duplicate);
    

    if (has_off_diag){
        for (i=0; i<nrows; i++) {
            temp_cols_off_diag[i] = (int*) tl_free(sizeof(int), nnz_rows_off_diag[i], temp_cols_off_diag[i]);
        }
        temp_cols_off_diag = (int **) tl_free(sizeof(int *), nrows, temp_cols_off_diag);
        nnz_rows_off_diag= (int *) tl_free(sizeof(int), nrows, nnz_rows_off_diag);
        nnz_rows_off_diag_no_duplicate= (int *) tl_free(sizeof(int), nrows, nnz_rows_off_diag_no_duplicate);
    }
    
    printf("CSR sparsity completed\n");

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
 *  @param[in,out]  SMODEL_SUPER *sm pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
 void allocate_adh_system(SMODEL_SUPER *sm){   
    // For now, we pre-allocate assuming that all nodes have the max number of equations attached.
    // This is not optimal, although the AdH storage is dumped to CCS and all extra zeros are removed before actually solving.
    // ndofs = nnodes * max_nsys (over-allocated for mixed dof systems)
    printf("allocating solution\n");
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
    printf("Made residual and solution\n");
    //where does local_range get determined?
    //we now number of rows is sum(my_nnode*nodal_nvars)
    sm->indptr_diag= (int *) tl_realloc(sizeof(int), sm->local_range[1]-sm->local_range[0]+1, sm->local_range_old[1]-sm->local_range_old[0]+1, sm->indptr_diag);
    //if we split we have two structures
    //if it is only one process this should be NULL!!!!
    //create conditional here
    sm->indptr_off_diag = (int *) tl_realloc(sizeof(int), sm->local_range[1]-sm->local_range[0]+1, sm->local_range_old[1]-sm->local_range_old[0]+1, sm->indptr_off_diag);
    printf("Made indptr\n");
    //find nnz and store in indices
    //this does a single CSR format
    //create_sparsity_CSR(sm, grid);
    //this does a split CSR format, plan on doing this
    create_sparsity_split_CSR(sm, sm->grid);
    //do we want to store nnz? it is stored in sm->indptr[nrows]
    //do we want to store local_size = local_range[1]-local_range[0]
    sarray_init_dbl(sm->vals_diag, sm->indptr_diag[sm->my_ndofs-1]);
    sarray_init_dbl(sm->vals_off_diag, sm->indptr_off_diag[sm->my_ndofs-1]);
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
 *  @param[in,out]  SMODEL_SUPER *sm pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void allocate_petsc_objects(SMODEL_SUPER *sm){
    int ierr;
    //set up anything PETSc after allocation but before load

 #ifdef _PETSC   
    //PETSc solver will be stord in ksp object
    //ksp only will exist if PETSc is active
    if(sm->ksp == PETSC_NULLPTR){
        ierr = KSPCreate(PETSC_COMM_WORLD, &(sm->ksp));
        ierr = KSPSetFromOptions(sm->ksp);
    }
    // Check if Jacobian, sol, and residual have already been created.
    // If so, destroy each of them before creating new PETSc objects.
    if(sm->A != PETSC_NULLPTR){
        printf("\n\nDestroying old matrix\n\n");
        ierr = MatDestroy(&(sm->A));
        ierr = KSPReset(sm->ksp);
        ierr = KSPSetFromOptions(sm->ksp);
    }
    //dont think we need preallocator matrix
    //if(sm->P != PETSC_NULLPTR){
    //    ierr = MatDestroy(&(sm->P));
    //}
    if(sm->X != PETSC_NULLPTR){
        ierr = VecDestroy(&(sm->X));
    }
    if(sm->B != PETSC_NULLPTR){
            ierr = VecDestroy(&(sm->B));
    }

    // Create Jacobian matrix from CSR format
    //don't know if i need to do this or just call after load anyway
    //vals should update by reference
    MatCreateMPIAIJWithSplitArrays(PETSC_COMM_WORLD, sm->my_ndofs, sm->my_ndofs, sm->macro_ndofs, sm->macro_ndofs, sm->indptr_diag, sm->cols_diag, sm->vals_diag, sm->indptr_off_diag, sm->cols_off_diag, sm->vals_off_diag, &(sm->A));
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
    ierr = VecCreateGhostWithArray(PETSC_COMM_WORLD, sm->my_ndofs, sm->macro_ndofs, sm->nghost, sm->ghosts, sm->residual, &(sm->B));
    ierr = VecCreateGhostWithArray(PETSC_COMM_WORLD, sm->my_ndofs, sm->macro_ndofs, sm->nghost, sm->ghosts, sm->sol, &(sm->X));


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
 *  @param[in,out]  SMODEL_SUPER *sm pointer to an instant of the SuperModel struct
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_allocate_initialize_linear_system(SMODEL_SUPER *sm) {

    if (sm->ndofs > sm->ndofs_old){    
        allocate_adh_system(sm);
#ifdef _PETSC
        allocate_petsc_objects(sm);
#endif
            }
}






void apply_Dirichlet_BC(SMODEL_SUPER *sm){
    //NOTE: THIS WILL BREAK IN PARALLEL, NEED TO ADD GHOSTS LATER

    //could illuciadate arguments, but will look at this later
    //modify linear system to have dirichlet conditions at
    //entries where bc_mask=YES

    //brute force for now, can find better way later
    int my_ndof = sm->my_ndofs;
    int row, row_start, row_end, row_entry,local_col_no;
    double temp;
    for (row=0;row<my_ndof;row++){
        row_start = sm->indptr_diag[row];
        row_end = sm->indptr_diag[row+1];
        //for every row, look for columns and subtract to RHS
        for(row_entry=row_start; row_entry<row_end; row_entry++){
            local_col_no = sm->cols_diag[row_entry];
            if(sm->bc_mask[local_col_no] == YES){
                temp = sm->vals_diag[row_entry];
                //adjust RHS
                //needed for method 1, for method 2 this should be unecceary
                //sm->residual[row] -= temp*(sm->dirichlet_data[local_col_no] - sm->sol[local_col_no]);
                //clear entry after moving
                sm->vals_diag[row_entry] = 0.0;
            }

        }

    }
    //go back through rows and clear rows that are no longer in equations
    //if it is dirichlet dof, wipe out row
    for (row=0;row<my_ndof;row++){
        if(sm->bc_mask[row] == YES){
            row_start = sm->indptr_diag[row];
            row_end = sm->indptr_diag[row+1];
            //assign dirichlet data to residual
            //method 1 appears to work
            //sm->residual[row] = sm->dirichlet_data[row] - sm->sol[row];
            //method 2, require sol to have bc already  there
            sm->residual[row] = 0.0;
            //wipe out row and set diagonal to 1 (rescaling diagonal shouldn't matter in method 2)
            for(row_entry=row_start; row_entry<row_end; row_entry++){
                local_col_no = sm->cols_diag[row_entry];
                //if non-diagonal set to 0, otherwise set diagonal to 1
                if(local_col_no == row){
                    sm->vals_diag[row_entry]= 1.0;
                }else{
                    sm->vals_diag[row_entry] = 0.0;
                }

            }

        }
    }

}



//int* unique (int* first, int* last)
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

int compare_ints(const void *a, const void *b) {
    return (*(int*)a - *(int*)b); // Ascending order
}



int unique(int *arr, int size){
    int unique_size = 1; // We'll start with the assumption that the first element is unique
    for (int i = 1; i < size; i++) {
        int is_unique = 1;
        for (int j = 0; j < unique_size; j++) {
            if (arr[i] == arr[j]) {
                is_unique = 0; // We found a duplicate
                break;
            }
        }
        if (is_unique) {
            arr[unique_size] = arr[i];
            unique_size++;
        }
    }
    return unique_size;
}


void Screen_print_CSR(int *indptr, int *cols, double *vals, int nrow){
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



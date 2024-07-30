 
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
    sm->indptr = (int *) tl_realloc(sizeof(int), sm->local_range[1]-sm->local_range[0]+1, sm->local_range_old[1]-sm->local_range_old[0]+1, sm->indptr);
    
    //find nnz and store in indices
    create_sparsity(sm, grid);
    //do we want to stor nnz? it is stored in sm->indptr[nrows]
    //do we want to store local_size = local_range[1]-local_range[0]
    sarray_init_dbl(sm->vals, sm->nnz);
}




void create_sparsity(sm, grid){
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
/*!
 *  \brief     Initializes the AdH Newton Jacobi matrix.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] SPARSE_VECT matrix, the sparse AdH matrix
 *  @param[in,out] double *diagonal the sparse AdH matrix diagonal blocks
 *  @param[in] int ndofs, the total number of equations to be solved in this linear system
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void init_adh_matrix(SPARSE_VECT *matrix, double *diagonal, int ndofs) {
    
    int i, iend;

    for (i = 0, iend = ndofs; i < iend; i++) {
        diagonal[i] = 0.0;
    }
    
    for (i = 0; i < nnodes; i++) {
        spv_init(matrix[i], ndofs*ndofs);
    }
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Initializes the sparse vec
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] SPARSE_VECT sv, the sparse matrix to be initializes
 *  @param[in,out] diagonal the sparse AdH matrix diagonal blocks
 *  @param[in,out] max_nsys_sq the maximum number of equations in system squared
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void spv_init(SPARSE_VECT sv,   /* the sparse vector */
              int max_nsys_sq)
{
    int i = 0;                  /* loop counter */
    int iend = 0;               /* the end of the loop */

    /* loops over the indices and initializes the values */
    for (i = 0, iend = sv.size * max_nsys_sq; i < iend; i++) {
        sv.value[i] = 0.0;
    }
    return;
}

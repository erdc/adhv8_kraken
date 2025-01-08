/*! \file  assemble_jacobian.c This file collections functions responsible for assembling the Jacobian matrix based on central finite difference  */
#include "adh.h"
static int DEBUG = OFF;
//static void sarray_init_double_2d_special(double to[MAX_ELEM_DOF][MAX_ELEM_DOF]);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function assembles the global FE Jacobian matrix elementwise, using F-D approximation to Jacobian
 *  designed to work within Newton method: \f$ \frac{\partial R (U^{i})}{\partial U^{i}} \delta U^{i} = -R(U^{i}) \f$
 *  With update formula: \f$ U^{i+1} = U^{i} + \delta U^{i}\f$.
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  sm (SMODEL_SUPER*) - pointer to an instant of the SMODEL_SUPER struct - contains pointers to sparse matrix
 * \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void assemble_jacobian(SMODEL_SUPER *sm) {
    SGRID *grid = sm->grid;
    int j,k;
    int* fmap = sm->dof_map_local;
    int* ghosts = sm->lin_sys->ghosts;
    int ndofs_ele;
    int nnodes, var_code;
    //maybe change this from 3d to 2d to 1d
    //in loop we dont want to have to alloc and free everytime, think about this
    double **elem_mat;
    elem_mat = (double**) tl_alloc(sizeof(double*), MAX_ELEM_DOF);
    for(j=0;j<MAX_ELEM_DOF;j++){
        elem_mat[j] = (double*) tl_alloc(sizeof(double), MAX_ELEM_DOF);
    }

    //better to do stack or heap?? 30x30 matrix, will keep on heap for now
    //double elem_mat[MAX_ELEM_DOF][MAX_ELEM_DOF];
    int dofs[MAX_ELEM_DOF],global_dofs[MAX_ELEM_DOF];

    //some alisases for convenience
    double *vals_diag = sm->lin_sys->vals_diag;
    double *vals_off_diag = sm->lin_sys->vals_off_diag;
    int *cols_diag = sm->lin_sys->cols_diag;
    int *cols_off_diag = sm->lin_sys->cols_off_diag;
    int *indptr_diag = sm->lin_sys->indptr_diag;
    int *indptr_off_diag = sm->lin_sys->indptr_off_diag;
    int nnz_diag = sm->lin_sys->nnz_diag;
    int nnz_off_diag = sm->lin_sys->nnz_off_diag;
    int local_range[2];
    local_range[0]= sm->lin_sys->local_range[0];
    local_range[1]= sm->lin_sys->local_range[1];
    //allocate 2d array, more memory than necessary
    //int vars_node[MAX_NNODE][MAX_NVAR];
//    int **vars_node;
//    vars_node = (int**) tl_alloc(sizeof(int*), MAX_NNODE);
//    for(j=0;j<MAX_NNODE;j++){
//        vars_node[j] = (int*) tl_alloc(sizeof(int), MAX_NVAR);
//        for(k=0;k<MAX_NVAR;k++){
//            vars_node[j][k]=0;
//        }
//    }


    //zero out stuff
    sarray_init_dbl(vals_diag, nnz_diag);
    //#ifdef? or avoid maybe
    if (indptr_off_diag !=NULL){
     sarray_init_dbl(vals_off_diag, nnz_off_diag);
    }



    int mat_id, nvars_elem, nphysics_models;
    int elem_vars[MAX_NVAR];
    //need sm->local_range whih is array of integers giving range (if process owned rows 0,1,2,3) then 
    // local range would be 0,4
    //loop through all nelem3d
    for (j=0;j<grid->nelems3d;j++){

        //pull all global information to local memory
        mat_id = sm->elem3d_physics_mat_id[j];
        //Get stuff from physics mat
        nvars_elem = sm->elem3d_physics_mat[mat_id].nvar;
        //get array of variables specified on element
        sarray_copy_int(elem_vars, sm->elem3d_physics_mat[mat_id].vars, nvars_elem);
        //number of phyics routines on element
        nphysics_models = sm->elem3d_physics_mat[mat_id].nSubmodels;
        //allows for mixed geometry quad/tri mesh
        nnodes = grid->elem3d[j].nnodes;
        ndofs_ele = nvars_elem*nnodes;
        //needs to be for 2d matrix
        sarray_init_double_2d(elem_mat,MAX_ELEM_DOF, MAX_ELEM_DOF);
        //sarray_init_double_2d_special(elem_mat);
        //maybe pull local nodal nvars and vars here too:
        //only necessary for CG i believe (or not idk)
        //if generalizing to higher dimension, this would be dof on basis functions
//        for (l=0;l<nnodes;l++){
//            node_id = grid->elem3d[j].nodes[l];
//            node_mat_id = sm->node_physics_mat_id[node_id];
//            nvar_node[l] = sm->node_physics_mat[node_mat_id].nvar;
//            for(m=0;m<nvar_node[l];m++){
//                vars_node[l][m] = sm->node_physics_mat[node_mat_id].vars[m];
//            }
//        }
        //need to loop over each variable and perturb it
        for (k=0;k<nvars_elem;k++){
            //use var code like PERTURB_H ,. ...
            //maybe get this from mat instead too!
            var_code = elem_vars[k];
            //want to loop over variables not necessarily submodels right?
            //need to think about this more
            //like sw2 is vector equation so that we would need 3 perturbations
            perturb_var(elem_mat, sm, sm->elem3d_physics_mat[mat_id].model, j, nnodes, nvars_elem, elem_vars,var_code, nphysics_models, k, grid->elem3d[j].nodes, DEBUG);
        }
        //store in global using 2 mappings
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //this gets dofs local to process
        //get_cell_dofs(dofs,fmap,nnodes,grid->elem3d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        get_cell_dofs(dofs,fmap,nnodes,grid->elem3d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat);
        //get_cell_dofs_2(dofs,nnodes,grid->elem3d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        //temp has local elemental matrix 
        //single CSR, we will keep for now but want to use split insted      
        //load_global_mat_CSR(sm->vals, sm->indptr, sm->indices, elem_mat, ndofs_ele, global_dofs, local_range);
        //split CSR
        load_global_mat_split_CSR(vals_diag, indptr_diag, cols_diag, vals_off_diag, indptr_off_diag, cols_off_diag, elem_mat, ndofs_ele, dofs, global_dofs, local_range);
    }

    //maybe change MAX_ELEM_DOF from element to element?
    //printf("Assembling jacobian\n");

    //2d loop
    for (j=0;j<grid->nelems2d;j++){
        //printf("Assembling element %d\n",j);
        //pull all global information to local memory
        mat_id = sm->elem2d_physics_mat_id[j];
        //Get stuff from physics mat
        nvars_elem = sm->elem2d_physics_mat[mat_id].nvar;
        //get array of variables specified on element
        sarray_copy_int(elem_vars, sm->elem2d_physics_mat[mat_id].vars, nvars_elem);
        //number of phyics routines on element
        nphysics_models = sm->elem2d_physics_mat[mat_id].nSubmodels;
        //allows for mixed geometry quad/tri mesh
        nnodes = grid->elem2d[j].nnodes;
        ndofs_ele = nvars_elem*nnodes;
        //needs to be for 2d matrix
        sarray_init_double_2d(elem_mat,MAX_ELEM_DOF, MAX_ELEM_DOF);
        //sarray_init_double_2d_special(elem_mat);
        //maybe pull local nodal nvars and vars here too:
        //only necessary for CG i believe (or not idk)
        //if generalizing to higher dimension, this would be dof on basis functions
//        for (l=0;l<nnodes;l++){
//            //printf("calling vars_node for node %d\n",l);
//            node_id = grid->elem2d[j].nodes[l];
//            node_mat_id = sm->node_physics_mat_id[node_id];
//            nvar_node[l] = sm->node_physics_mat[node_mat_id].nvar;
//            //printf("node mat id: %d, Nvar on this node%d\n",node_mat_id, nvar_node[l]);
//            for(m=0;m<nvar_node[l];m++){//

//                vars_node[l][m] = sm->node_physics_mat[node_mat_id].vars[m];
//                //printf("Node physics var assigned for node %d node var %d = %d\n",l,m,vars_node[l][m]);
//            }
//        //printf("Looping through eac variable on the elemnt\n");
//        }
        //need to loop over each variable and perturb it
        //printf("completed loop\n");
        for (k=0;k<nvars_elem;k++){
            //use var code like PERTURB_H ,. ...
            //maybe get this from mat instead too!
            var_code = elem_vars[k];
            //printf("elem vars [%d] = %d\n",k,var_code);
            //want to loop over variables not necessarily submodels right?
            //need to think about this more
            //like sw2 is vector equation so that we would need 3 perturbations
            perturb_var(elem_mat, sm, sm->elem2d_physics_mat[mat_id].model, j, nnodes, nvars_elem, elem_vars,var_code, nphysics_models, k, grid->elem2d[j].nodes, DEBUG);
            //printf("perturb var called\n");
        }
        //store in global using 2 mappings
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //this gets dofs local to process
        //get_cell_dofs(dofs,fmap,nnodes,grid->elem2d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        get_cell_dofs(dofs,fmap,nnodes,grid->elem2d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat);
        //get_cell_dofs_2(dofs,nnodes,grid->elem2d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        //temp has local elemental matrix 
        //single CSR, we will keep for now but want to use split insted      
        //load_global_mat_CSR(sm->vals, sm->indptr, sm->indices, elem_mat, ndofs_ele, global_dofs, local_range);
        //split CSR
        load_global_mat_split_CSR(vals_diag, indptr_diag, cols_diag, vals_off_diag, indptr_off_diag, cols_off_diag, elem_mat, ndofs_ele, dofs, global_dofs, local_range);
    }
    //loop through all nelem1d
    for (j=0;j<grid->nelems1d;j++){

        //pull all global information to local memory
        mat_id = sm->elem1d_physics_mat_id[j];
        //Get stuff from physics mat
        nvars_elem = sm->elem1d_physics_mat[mat_id].nvar;
        //get array of variables specified on element
        sarray_copy_int(elem_vars, sm->elem1d_physics_mat[mat_id].vars, nvars_elem);
        //number of phyics routines on element
        nphysics_models = sm->elem1d_physics_mat[mat_id].nSubmodels;
        //allows for mixed geometry quad/tri mesh
        nnodes = grid->elem1d[j].nnodes;
        ndofs_ele = nvars_elem*nnodes;
        //needs to be for 2d matrix
        sarray_init_double_2d(elem_mat,MAX_ELEM_DOF, MAX_ELEM_DOF);
        //sarray_init_double_2d_special(elem_mat);
        //maybe pull local nodal nvars and vars here too:
        //only necessary for CG i believe (or not idk)
        //if generalizing to higher dimension, this would be dof on basis functions
//        for (l=0;l<nnodes;l++){
//            node_id = grid->elem1d[j].nodes[l];
//            node_mat_id = sm->node_physics_mat_id[node_id];
//            nvar_node[l] = sm->node_physics_mat[node_mat_id].nvar;
//            for(m=0;m<nvar_node[l];m++){
//                vars_node[l][m] = sm->node_physics_mat[node_mat_id].vars[m];
//            }
//        }
        //need to loop over each variable and perturb it
        for (k=0;k<nvars_elem;k++){
            //use var code like PERTURB_H ,. ...
            //maybe get this from mat instead too!
            var_code = elem_vars[k];
            //want to loop over variables not necessarily submodels right?
            //need to think about this more
            //like sw2 is vector equation so that we would need 3 perturbations
            perturb_var(elem_mat, sm, sm->elem1d_physics_mat[mat_id].model, j, nnodes, nvars_elem, elem_vars,var_code, nphysics_models, k, grid->elem1d[j].nodes, DEBUG);
        }
        //store in global using 2 mappings
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //this gets dofs local to process
        //get_cell_dofs(dofs,fmap,nnodes,grid->elem1d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        get_cell_dofs(dofs,fmap,nnodes,grid->elem1d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat);
        //get_cell_dofs_2(dofs,nnodes,grid->elem1d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        //temp has local elemental matrix 
        //single CSR, we will keep for now but want to use split insted      
        //load_global_mat_CSR(sm->vals, sm->indptr, sm->indices, elem_mat, ndofs_ele, global_dofs, local_range);
        //split CSR
        load_global_mat_split_CSR(vals_diag, indptr_diag, cols_diag, vals_off_diag, indptr_off_diag, cols_off_diag, elem_mat, ndofs_ele, dofs, global_dofs, local_range);
    }
    //free memory, in time loop need to be smarter about this
    for(j=0;j<MAX_ELEM_DOF;j++){
        elem_mat[j] = (double *) tl_free(sizeof(double), MAX_ELEM_DOF, elem_mat[j]);
    }
    elem_mat = (double **) tl_free(sizeof(double *), MAX_ELEM_DOF, elem_mat);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function takes an elemental stiffness matrix and loads it to the full sparse matrix
 *  in split CSR format
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  vals (double*) - pointer to the double array containing all nonzero values in the diagonal block (local to process) part of the sparse matrix
 *  @param[in]  indptr (int*) - pointer to the array containing indeces to start/end of each row in the diagonal block (local to process) part of the sparse matrix
 *  @param[in]  indices (int*) - pointer to the array containing local column numbers of each nonzero entry in the diagonal block (local to process) part of the sparse matrix
 *  @param[in,out]  off_diag_vals (double*) - pointer to the values of the off-diagonal block (local to process) part of the sparse matrix
 *  @param[in]  off_diag_indptr (int*) - pointer to the array containing indeces to start/end of each row in the off-diagonal block (local to process) part of the sparse matrix
 *  @param[in]  off_diag_indices (int*) - pointer to the array containing global column numbers of each nonzero entry in the off-diagonal block (local to process) part of the sparse matrix
 *  @param[in] elem_mat (double**) - 2D array containing computed values from an elemental stiffness matrix
 *  @param[in] ndofs_ele (int) - number of degrees of freedom on the element
 *  @param[in] dofs (int*) - the degrees of freedom present in the element (local to the process)
 *  @param[in] global_dofs (int*) - the global degrees of freedom present in the element
 *  @param[in] local_range (int*) - an array of 2 integers that gives global start and end equation number
 * \note This function does depend on a binary search algorithm
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void load_global_mat_split_CSR(double *vals, int *indptr, int *indices, double *off_diag_vals, int *off_diag_indptr, int *off_diag_indices, double **elem_mat, int ndofs_ele, int *dofs, int *global_dofs, int *local_range){
    int i,j;
    //need local to global mapping stores in global_dofs
    //assume we have at least global row start that this PE owns
    int col_start_ind,col_end_ind,col_index;
    int row,col, global_col;
    int nrows = local_range[1]-local_range[0];


    /// assembles global residual
    // i and j are the row number and col number in elemental matrix
    //we only want rows local to process
    //this can easily be checked with global_dofs argument
    for (i=0; i<ndofs_ele; i++){
        //global row number
        row = dofs[i];
        //loop through each row of the local elemental matrix
        //we only need to fill in columns if the current row is local to process
        if(row<nrows){
            for (j=0; j<ndofs_ele;j++){
                //differentiate  between diag and off diag blocks
                if (row<nrows){
                    //then this is on diagonal block
                    //get column indices for this specific row
                    col_start_ind = indptr[row];
                    col_end_ind = indptr[row+1];
                    //get local column number
                    col = global_dofs[j];
                    //get offset of where to place value by binary search through part of the indices array
                    col_index = binary_search_part(indices, col_start_ind, col_end_ind, col);
                    //now that we have the index, add to values
                    vals[col_index] += elem_mat[i][j];
                }
                else{
                    //otherwise we are on off-diagonal part
                    //then this is on diagonal block
                    //get column indices for this specific row
                    col_start_ind = off_diag_indptr[row];
                    col_end_ind = off_diag_indptr[row+1];
                    //get global column number
                    global_col = global_dofs[j];
                    //get offset of where to place value by binary search through part of the indices array
                    col_index = binary_search_part(off_diag_indices, col_start_ind, col_end_ind, global_col);
                    //now that we have the index, add to values
                    off_diag_vals[col_index] +=elem_mat[i][j];

                }
            }
        }
    }

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function takes an elemental stiffness matrix and loads it to the full sparse matrix
 *  in standard CSR format
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  vals (double*) - pointer to the double array containing all nonzero values in the diagonal block (local to process) part of the sparse matrix
 *  @param[in]  indptr (int*) - pointer to the array containing indeces to start/end of each row in the diagonal block (local to process) part of the sparse matrix
 *  @param[in]  indices (int*) - pointer to the array containing local column numbers of each nonzero entry in the diagonal block (local to process) part of the sparse matrix
 *  @param[in] elem_mat (double**) - 2D array containing computed values from an elemental stiffness matrix
 *  @param[in] ndofs_ele (int) - number of degrees of freedom on the element
 *  @param[in] global_dofs (int*) - the global degrees of freedom present in the element
 *  @param[in] local_range (int*) - an array of 2 integers that gives global start and end equation number
 * \note This function does depend on a binary search algorithm
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void load_global_mat_CSR(double *vals, int *indptr, int *indices, double **elem_mat, int ndofs_ele, int *global_dofs, int *local_range){
    int i,j;
    //need local to global mapping stores in global_dofs
    //assume we have at least global row start that this PE owns
    int col_start_ind,col_end_ind,col_index;
    int global_row,row, global_col;

    /// assembles global residual
    // i and j are the row number and col number in elemental matrix
    //we only want rows local to process
    //this can easily be checked with global_dofs argument
    for (i=0; i<ndofs_ele; i++){
        //global row number
        global_row = global_dofs[i];
        //loop through each row of the local elemental matrix
        //we only need to fill in columns if the current row is local to process
        if(global_row>=local_range[0] && global_row < local_range[1]){
            //get row local to process
            row = global_row-local_range[0];
            //get column indices for this specific row
            col_start_ind = indptr[row];
            col_end_ind = indptr[row+1];
            for (j=0; j<ndofs_ele;j++){
                //get global column number
                global_col = global_dofs[j];
                //get offset of where to place value by binary search through part of the indices array
                col_index = binary_search_part(indices, col_start_ind, col_end_ind, global_col);
                //now that we have the index, add to values
                vals[col_index] += elem_mat[i][j];
            }
        }
    }

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes elemental Jacobian based on central F-D technique w.r.t one variable 
 * (thus filling out one column of the elemental Jacobian).
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] elem_mat (double**) - stores the Jacobian, elemental matrix
 *  @param[in] sm (SMODEL_SUPER*) - a pointer to the SMODEL_SUPER structure, needed for residual routines
 *  @param[in] model (SMODEL*) - an array of SMODEL structures, needed to get proper residual routines
 *  @param[in] ie (int) - the element number (local to process)
 *  @param[in] nodes_on_element (int) - the number of nodes on the element
 *  @param[in] nvar_ele (int) - the total number of variables active on element across all residual routines
 *  @param[in] elem_vars (int*) - the variable codes, an array of length nvar_ele
 *  @param[in] perturb_var_code (int) - the variable that is being differentiated
 *  @param[in] nsubModels (int) - number of residual routines to be called on this element
 *  @param[in] ele_var_no (int) - the index of the variable that is being differentiated (not the var code)
 *  @param[in] NodeIDs (int*) - the node numbers within the cell (local to the process)
 *  @param[in] DEBUG (int) - the debug code for more robust output
 *  @param[in] DEBUG a debug option
 *  \note The F-D spacing is based on magnitude of solution with lower limit hard coded as 1e-4. Catch is too small
 *  small of F-D space leads to large truncation error (division by very small number) but higher F-D space incurs
 *  large truncation error especially the more nonlinear the residual is
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void perturb_var(double **elem_mat, SMODEL_SUPER *sm, SMODEL *model, 
    int ie, int nodes_on_element, int nvar_ele, int *elem_vars ,int perturb_var_code, int nsubModels, int ele_var_no, int *NodeIDs, int DEBUG) {
    
    
    int i,j, NodeID = UNSET_INT;
    double epsilon = 0., epsilon2 = 0.;
    

    //pertrubation may depend on physics?
    //mark changed, the sqrt(SMALL) seems to be quite excessive
    double perturbation = 1e-4;//sqrt(SMALL);


    double temp_P[MAX_ELEM_DOF];
    double temp_M[MAX_ELEM_DOF];
    double elem_rhs_P[MAX_ELEM_DOF];    // +head perturbation initialized by elem_resid
    double elem_rhs_M[MAX_ELEM_DOF];    // -head perturbation initialized by elem_resid
    int eq_var_code,eq_var_code2,nvar_pde;
    int physics_vars[MAX_NVAR];
    //Mark add switch cases for all sol_var
    //switch case for which sm->sol_var
    double temp_sol;
//    if (perturb_var_code == PERTURB_H){
//        temp_sol = sm->head;
//    }else if (perturb_var_code == PERTURB_U){
//        temp_sol = sm->vel2d->x;
//    }else if (perturb_var_code == PERTURB_V){
//        temp_sol = sm->vel2d->y
//    }   
    for (i=0; i<nodes_on_element; i++) {
        NodeID = NodeIDs[i];
        //pull dof # to get info
        //for now use CG map in here

        //temp_sol = sm->sol[get_cg_dof_2(perturb_var_code, NodeID, sm->node_physics_mat, sm->node_physics_mat_id)];
        //temp_sol = sm->sol[get_cg_dof(perturb_var_code, NodeID, sm->dof_map_local, sm->node_physics_mat, sm->node_physics_mat_id)];
        temp_sol = sm->sol[get_cg_dof(perturb_var_code, NodeID, sm->dof_map_local, sm->node_physics_mat)];
        NUM_DIFF_EPSILON(epsilon, epsilon2, temp_sol, perturbation);    // calculates epsilon and 2*epsilon
        
        //epsilon = 1.0;
        //epsilon2 = 2.0;
        sarray_init_dbl(elem_rhs_P,MAX_NVAR*MAX_NNODE);
        sarray_init_dbl(elem_rhs_M,MAX_NVAR*MAX_NNODE);
        //safe assumption that all submodels on an element depend on all of the variables?
        for (j=0;j<nsubModels;j++){ 
            sarray_init_dbl(temp_P,MAX_NVAR*MAX_NNODE);
            sarray_init_dbl(temp_M,MAX_NVAR*MAX_NNODE);        
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) body perturbation of depth ++++++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: body pertubing +h for element %d :: perturbation: %20.10e\n",ie,epsilon);
#endif

            //this will give a local residual, in temp and will give code for model vars
            nvar_pde = model[j].nvar;
            sarray_init_int(physics_vars, nvar_pde);
            sarray_copy_int(physics_vars, model[j].physics_vars,nvar_pde);
            //eq_var_code = model[j].fe_resid(sm,temp_P,ie, epsilon,i, perturb_var_code, +1, DEBUG);
            eq_var_code = smodel_super_resid(sm,temp_P,ie, epsilon, i, perturb_var_code, +1, DEBUG, fe_resid[model[j].fe_resid]);
            add_replace_elem_rhs(elem_rhs_P,temp_P,nvar_ele,elem_vars,nvar_pde,physics_vars,nodes_on_element, 1);
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of depth ++++++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: body pertubing -h for element %d :: perturbation: %20.10e\n",ie,epsilon);
#endif
            //fe_sw2_body_resid(mod,elem_rhs_h_M,ie,epsilon,i,PERTURB_H,-1,DEBUG);
            //should always be same as var_code
            //eq_var_code2 = model[j].fe_resid(sm,temp_M,ie, epsilon,i, perturb_var_code, -1, DEBUG);
            eq_var_code2 = smodel_super_resid(sm,temp_M,ie, epsilon, i, perturb_var_code, -1, DEBUG, fe_resid[model[j].fe_resid]);
            add_replace_elem_rhs(elem_rhs_M,temp_M,nvar_ele,elem_vars,nvar_pde,physics_vars,nodes_on_element, 1);
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        }
        // calculate local contribution to a local elemental Jacobian matrix++++++++++++++++++++++++++++++
        // fills in one column
        elem_matrix_deriv(elem_mat,i, ele_var_no, nodes_on_element, nvar_ele, elem_rhs_P, elem_rhs_M, epsilon2); // gradient
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a scond order finite difference of a column within the Jacobian matrix elemental block
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] mat (double**) - the 2D array storing the values of the elemental stiffness matrix
 *  @param[in] node_no (int) - local node number (local to element)
 *  @param[in] var_no (int) - the solution variable number (local to element)
 *  @param[in] nnodes (int) - number of nodes on the element
 *  @param[in] elem_nvars (int) - number of solution variables active on the element
 *  @param[in] local1 (double*) - elemental residual vector with a (+) perturbation
 *  @param[in] local2 (double*) - elemental residual vector with a (-) perturbation
 *  @param[in] diff_ep (double) - two times the perturbations size
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void elem_matrix_deriv(double **mat, int node_no, int var_no, int nnodes, int elem_nvars, double *local1, double *local2, double diff_ep) {
    int inode,jvar, col_no, indx;
    col_no = node_no*elem_nvars + var_no;

    for (inode=0; inode<nnodes; inode++) {
        for (jvar=0; jvar<elem_nvars; jvar++) {
            indx = inode*elem_nvars + jvar;
            mat[indx][col_no]= (local1[indx] - local2[indx])/diff_ep;
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Short helper routine, seeing if using small array on stack helps performance, 
 *  doesnt seem to make difference
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] to (double[MAX_ELEM_DOF][MAX_ELEM_DOF]) - the 2D array storing the values of the elemental stiffness matrix
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//static void sarray_init_double_2d_special(double to[MAX_ELEM_DOF][MAX_ELEM_DOF]) {
//    int i,j;
//    for (i=0; i<MAX_ELEM_DOF; i++) {
//        for (j=0; j<MAX_ELEM_DOF; j++) {
//            to[i][j] = 0.0;
//        }
//    }
//}

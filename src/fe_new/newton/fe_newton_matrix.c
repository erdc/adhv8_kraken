/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function assembles the global FE Jacobian matrix elementwise
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct - adjusts the matrix
 *  @param[in] SGRID *grid - the grid
 *  @param[in] SMAT *mat - the set of materials
 * \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void assemble_matrix(SSUPER_MODEL *sm, SGRID *grid, SMAT *mat) {
    int j,k;
    double fmap = sm->fmap;
    int nsubMods;
    int ndofs_ele;
    int nvar_ele,nodes_on_ele, var_code;

    //maybe change this from 3 to 2 to 1
    int max_elem_dofs = MAX_NVAR*MAX_NNODE; 
    double elem_mat[max_elem_dofs,max_elem_dofs];
    int dofs[max_elem_dofs],global_dofs[max_elem_dofs];

    //zero out stuff
    sarray_init_dbl(sm->vals, sm->nnz_old);
    local_range = ??sm->local_range
    //need sm->local_range whih is array of integers giving range (if process owned rows 0,1,2,3) then 
    // local range would be 0,4
    //loop through all nelem3d
    for (j=0;j<grid->nelem3d;j++){

            //pull all global information to local memory

            //Get this from mat instead!!
            nvar_ele = sm->elem3d_nvars[j];
            //allows for mixed geometry quad/tri mesh
            nodes_on_ele = grid->elem3d[j].nnodes;
            ndofs_ele = nvar_ele*nodes_on_ele;
            nsubMods = sm->nSubMods3d[j];
            //needs to be for 2d matrix
            sarray_init_dbl(elem_mat,max_elem_dofs*max_elem_dofs);
            //maybe pull local nodal nvars and vars here too
            //need to loop over each variable and perturb it
            for (k=0;k<nvar_ele;k++){

                //use var code like PERTURB_H ,. ...
                //maybe get this from mat instead too!
                var_code = sm->elem3d_vars[j][k];
                //want to loop over variables not necessarily submodels right?
                //need to think about this more
                //like sw2 is vector equation so that we would need 3 perturbations
                perturb_var(elem_mat, sm, grid, mat, sm->elem3d_physics, j, nodes_on_ele, nvar_ele, sm->elem3d_vars[j],var_code, nsubMods, k, grid->elem3d[j].nodes, DEBUG);
            }
            //store in global using 2 mappings
            //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
            //this gets dofs local to process
            get_cell_dofs(dofs,fmap,nnodes,grid->elem3d[j].nodes,nvars_elem,sm->elem3d_vars[j],sm->nodal_nvars, sm->nodal_vars);
            //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
            local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
            //temp has local elemental matrix       
            load_global_mat(sm->vals, sm->indptr, sm->indices, elem_mat, ndofs_ele, global_dofs, local_range)
    }

    //maybe change max_elem_dofs from element to element?


    //loop through all nelem2d
    for (j=0;j<grid->nelem2d;j++){
        nvar_ele = sm->elem2d_nvars[j];
        nodes_on_ele = grid->elem2d[j].nnodes;
        ndofs_ele = nvar_ele*nodes_on_ele;
        nsubMods = sm->nSubMods2d[j];
        //needs to be for 2d matrix
        sarray_init_dbl(elem_mat,maxe_elem_dofs*max_elem_dofs);
        for (k=0;k<nvar_ele;k++){
            //use var code like PERTURB_H ,. ...
            var_code = sm->elem2d_vars[j][k];
            //want to loop over variables not necessarily submodels right?
            //need to think about this more
            //like sw2 is vector equation so that we would need 3 perturbations
            perturb_var(elem_mat, sm, grid, mat, sm->elem2d_physics, j, nodes_on_ele, nvar_ele, sm->elem2d_vars[j],var_code, nsubMods, k, grid->elem2d[j].nodes, DEBUG);
        }
        //store in global using 2 mappings
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //this gets dofs local to process
        get_cell_dofs(dofs,fmap,nnodes,grid->elem2d[j].nodes,nvars_elem,sm->elem2d_vars[j],sm->nodal_nvars, sm->nodal_vars);
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        //temp has local elemental matrix       
        load_global_mat(sm->vals, sm->indptr, sm->indices, elem_mat, ndofs_ele, global_dofs, local_range);
    
    }
    //loop through all nelem1d
    for (j=0;j<grid->nelem1d;j++){
        nvar_ele = sm->elem1d_nvars[j];
        nodes_on_ele = grid->elem1d[j].nnodes;
        ndofs_ele = nvar_ele*nodes_on_ele;
        nsubMods = sm->nSubMods1d[j];
        //needs to be for 2d matrix
        sarray_init_dbl(elem_mat,maxe_elem_dofs*max_elem_dofs);
        for (k=0;k<nvar_ele;k++){
            //use var code like PERTURB_H ,. ...
            var_code = sm->elem1d_vars[j][k];
            //want to loop over variables not necessarily submodels right?
            //need to think about this more
            //like sw2 is vector equation so that we would need 3 perturbations
            perturb_var(elem_mat, sm, grid, mat, sm->elem1d_physics, j, nodes_on_ele, nvar_ele, sm->elem1d_vars[j],var_code, nsubMods, k, grid->elem1d[j].nodes, DEBUG);
        }
        //store in global using 2 mappings
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //this gets dofs local to process
        get_cell_dofs(dofs,fmap,nnodes,grid->elem1d[j].nodes,nvars_elem,sm->elem3d_vars[j],sm->nodal_nvars, sm->nodal_vars);
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        //temp has local elemental matrix       
        load_global_mat(sm->vals, sm->indptr, sm->indices, elem_mat, ndofs_ele, global_dofs, local_range);
    }    
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function maps a local, dense matrix from an element and put into global sparse matrix
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  double *diagonal pointer to the diagonal  (block) part of the sparse matrix
 *  @param[in,out]  SPARSE_VECT *matrix pointer to the off-diagonal part of the sparse matrix
 *  @param[in] double **elem_mat - the local matrix that contains values we want to add to global sparse system
 *  @param[in] int nodes_on_ele - the number of nodes on element
 *  @param[in] int *fmap - a map from node number to row # in system of equations
 *  @param[in] in elem_nvars - number of distinct variables active on an element e.g. h,u,v would be 3
 *  @param[in] int *elem_vars - an array of length elem_nvars that contains the variable codes
 *  @param[in] in nodal_nvars - array of distinct variables active on each node e.g. h,u,v would be 3
 *  @param[in] int **nodal_nvars - an array of length nodal_nvars*nnodes that contains the variable codes
 * \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void local_dofs_to_global_dofs(int *global_dofs,int ndofs_on_ele,int *dofs,int *local_range,int *ghosts){
    
    int i;
    int local_size = local_range[1]-local_range[0];


    for(i=0;i<ndofs_on_ele;i++){
        if(dofs[i]<local_size){
            //residential dof
            global_dofs[i] = dofs[i] + local_range[0];
        }else{
            global_dofs[i] = ghosts[dofs[i]-local_size];
        }
    }


}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function maps a local, dense matrix from an element and put into global sparse matrix
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  double *diagonal pointer to the diagonal  (block) part of the sparse matrix
 *  @param[in,out]  SPARSE_VECT *matrix pointer to the off-diagonal part of the sparse matrix
 *  @param[in] double **elem_mat - the local matrix that contains values we want to add to global sparse system
 *  @param[in] int nodes_on_ele - the number of nodes on element
 *  @param[in] int *fmap - a map from node number to row # in system of equations
 *  @param[in] in elem_nvars - number of distinct variables active on an element e.g. h,u,v would be 3
 *  @param[in] int *elem_vars - an array of length elem_nvars that contains the variable codes
 *  @param[in] in nodal_nvars - array of distinct variables active on each node e.g. h,u,v would be 3
 *  @param[in] int **nodal_nvars - an array of length nodal_nvars*nnodes that contains the variable codes
 * \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void load_global_mat(double *vals, int *indptr, int *indices, double **elem_mat, int ndofs_ele, int *global_dofs, int *local_range){
    int i,j;
    //need local to global mapping stores in global_dofs
    //assume we have at least global row start that this PE owns
    int col_start_ind,col_end_ind,col_index;
    int global_row,row;

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
            row = global_row-local_range[0]
            //get column indices for this specific row
            col_start_ind = indptr[row];
            col_end_ind = indptr[row+1];
            for (j=0; j<ndofs_ele;j++){
                //get global column number
                global_col = global_dofs[j];
                //get offset of where to place value by binary search through part of the indices array
                col_index = binary_search_part(indices, col_start_ind, col_end_ind, global_col);
                //now that we have the index, add to values
                vals[col_index] +=elem_mat[i,j];
            }
        }
    }

}


//a small helper function, move this later
int binary_search_part(int *arr, int start, int end, int target) {
    int low = start;
    int high = end - 1;

    while (low <= high) {
        int mid = low + (high - low) / 2;

        if (arr[mid] == target) {
            return mid; // Target found
        } else if (arr[mid] < target) {
            low = mid + 1; // Search upper half
        } else {
            high = mid - 1; // Search lower half
        }
    }

    return -1; // Target not found
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds perturbation to a variable for each node in element for Newton Jacobian calculation.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] elem_mat stores the Jacobian, elemental matrix for the head perturbation
 *  @param[in]  sm a pointer to an AdH supermodel
 *  @param[in] nodes_on_element the number of nodes on the element
 *  @param[in] ie the element
 *  @param[in] dim the dimension of the resid
 *  @param[in] DEBUG a debug option
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void perturb_var(double **elem_mat, SSUPER_MODEL *sm, SGRID *grid, SMAT *mat, SELEM_PHYSICS **elem_physics, 
    int ie, int nodes_on_element, int nvar_ele, int *elem_vars ,int perturb_var_code, int nsubModels, int ele_var_no, int *GlobalNodeIDs, int DEBUG) {
    
    
    int i,j, GlobalNodeID = UNSET_INT;
    double epsilon = 0., epsilon2 = 0.;
    

    //pertrubation may depend on physics?
    double perturbation = sqrt(SMALL);


    double temp_P[MAX_NVAR*MAX_NNODE];
    double temp_M[MAX_NVAR*MAX_NNODE];
    double elem_rhs_P[MAX_NVAR*MAX_NNODE];    // +head perturbation initialized by elem_resid
    double elem_rhs_M[MAX_NVAR*MAX_NNODE];    // -head perturbation initialized by elem_resid
    int eq_var_code,eq_var_code2;

        

    //Mark add switch cases for all sol_var
    //switch case for which sm->sol_var
    double temp_sol[nnode];
    if (perturb_var_code == PERTURB_H){
        temp_sol = sm->head;
    }else if (perturb_var_code == PERTURB_U){
        temp_sol = sm->vel2d.x;
    }else if (perturb_var_code == PERTURB_V){
        temp_sol = sm->vel2d.y
    }

    
    for (i=0; i<nodes_on_element; i++) {
        GlobalNodeID = GlobalNodeIDs[i];
        NUM_DIFF_EPSILON(epsilon, epsilon2, temp_sol[GlobalNodeID], perturbation);    // calculates epsilon and 2*epsilon
        sarray_init_dbl(elem_rhs_P,MAX_NVAR*MAX_NNODE);
        sarray_init_dbl(elem_rhs_M,MAX_NVAR*MAX_NNODE);

        //safe assumption that all submodels on an element depend on all of the variables?
        for (j=0;j<nsubModels;j++){ 
            sarray_init_dbl(temp_P,MAX_NVAR*MAX_NNODE);
            sarray_init_dbl(temp_M,MAX_NVAR*MAX_NNODE);        

            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) body perturbation of depth ++++++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: body pertubing +h for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon);
#endif

            //this will give a local residual, in temp and will give code for model vars
            eq_var_code = elem_physics[ie][j]->fe_resid(sm,temp_P,grid,mat,ie, epsilon,i, perturb_var_code, +1, DEBUG);
            add_replace_elem_rhs(elem_rhs_P,temp_P,elem_vars,nvar_ele,eq_var_code,nodes_on_element);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of depth ++++++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: body pertubing -h for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon);
#endif
            //fe_sw2_body_resid(mod,elem_rhs_h_M,ie,epsilon,i,PERTURB_H,-1,DEBUG);
            //should always be same as var_code
            eq_var_code2 = elem_physics[ie][j]->fe_resid(sm,temp_M,grid,mat,ie, epsilon,i, perturb_var_code, -1, DEBUG);
            add_replace_elem_rhs(elem_rhs_M,temp_M,elem_vars,nvar_ele,eq_var_code2,nodes_on_element);
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        }
        // calculate local contribution to a local elemental Jacobian matrix++++++++++++++++++++++++++++++
        // fills in one column
        elem_matrix_deriv(i, ele_var_no, nodes_on_element, nvar_ele, elem_rhs_P, elem_rhs_M, elem_mat, epsilon2); // gradient
    }
    
    //DEBUG =  OFF;
}



/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a scond order finite difference of a Jacobi matrix elemental block
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] mat a 3 DOF matrix which stores the gradients
 *  @param[in] indx the block starting position
 *  @param[in] nnodes the total number of nodes on the element
 *  @param[in] local1 a 3 DOF residual with a (+) perturbation
 *  @param[in] local2 a 3 DOF residual with a (-) perturbation
 *  @param[in] diff_ep two times the perturbations size
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void elem_matrix_deriv(int node_no, int dof_no, int nnodes, int elem_nvars, double *local1, double *local2, double *mat, double diff_ep) {
    int inode,jvar, col_no, indx;
    col_no = node_no*elem_nvars + dof_no;

    for (inode=0; inode<nnodes; inode++) {
        for (jvar=0; jvar<elem_nvars; jvar++) {
            indx = inode*elem_nvars + jvar;
            mat[indx,col_no]= (local1[indx] - local2[indx])/diff_ep;
        }
    }
}
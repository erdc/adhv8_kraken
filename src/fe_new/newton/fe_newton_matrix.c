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

    //seems like the easiest way?
    //maybe think about this
    double fmap = sm->fmap;
    int nsubMods;
    int ndofs;
    int nvar_ele,nodes_on_ele, var_code;
    int max_elem_dofs = MAX_NVAR*MAX_NNODE; 
    double temp[max_elem_dofs,max_elem_dofs];

    //zero out stuff
    sarray_init_dbl(sm->vals, sm->nnz_old);
    

    //loop through all nelem3d
    for (j=0;j<grid->nelem3d;j++){

            //pull all global information to local memory
            nvar_ele = sm->elem3d_nvars[j];
            nodes_on_ele = grid->elem3d[j].nnodes;
            ndofs = nvar_ele*nodes_on_ele;
            nsubMods = sm->nSubMods3d[j];
            //needs to be for 2d matrix
            sarray_init_dbl(temp,maxe_elem_dofs*max_elem_dofs);
            //maybe pull local nodal nvars and vars here too
            //need to loop over each variable and perturb it
            for (k=0;k<nvar_ele;k++){

                //use var code like PERTURB_H ,. ...
                var_code = sm->elem3d_vars[j][k];
                //want to loop over variables not necessarily submodels right?
                //need to think about this more
                //like sw2 is vector equation so that we would need 3 perturbations
                perturb_var(temp, sm, grid, mat, sm->elem3d_physics, j, nodes_on_ele, nvar_ele, sm->elem3d_vars[j],var_code, nsubMods, k, grid->elem3d[j].nodes, DEBUG);
            }
            //store in global using 2 mappings
            //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
            //this gets dofs local to process
            get_cell_dofs(dofs,fmaplocal,nnodes,grid->elem3d[j].nodes,nvars_elem,sm->elem3d_vars[j],sm->nodal_nvars, sm->nodal_vars);
            //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
            local_dofs_to_global_dofs(global_dofs,ndofs,dofs,local_range,ghosts);
            //temp has local elemental matrix       
            load_global_mat(sm->vals, sm->indptr, sm->indices, temp, ndofd, global_dofs);
    }

    //loop through all nelem2d
    for (j=0;j<grid->nelem2d;j++){
        nvar_ele = sm->elem2d_nvars[j];
        nodes_on_ele = grid->elem2d[j].nnodes;
        nsubMods = sm->nSubMods2d[j];
        //needs to be for 2d matrix
        sarray_init_dbl(temp,maxe_elem_dofs*max_elem_dofs);
        for (k=0;k<nvar_ele;k++){
            //use var code like PERTURB_H ,. ...
            var_code = sm->elem2d_vars[j][k];
            //want to loop over variables not necessarily submodels right?
            //need to think about this more
            //like sw2 is vector equation so that we would need 3 perturbations
            perturb_var(temp, sm, grid, mat, sm->elem2d_physics, j, nodes_on_ele, nvar_ele, sm->elem2d_vars[j],var_code, nsubMods, k, grid->elem2d[j].nodes, DEBUG);
        }
        //store in global using the fmap
        //temp has local elemental matrix 
   
            load_global_mat(sm->diagonal, sm->matrix, temp, nodes_on_ele, nvar_ele, global_dofs);
    }
    //loop through all nelem1d
    for (j=0;j<grid->nelem1d;j++){
        nvar_ele = sm->elem1d_nvars[j];
        nodes_on_ele = grid->elem1d[j].nnodes;
        nsubMods = sm->nSubMods1d[j];
        //needs to be for 2d matrix
        sarray_init_dbl(temp,maxe_elem_dofs*max_elem_dofs);
        for (k=0;k<nvar_ele;k++){
            //use var code like PERTURB_H ,. ...
            var_code = sm->elem1d_vars[j][k];
            //want to loop over variables not necessarily submodels right?
            //need to think about this more
            //like sw2 is vector equation so that we would need 3 perturbations
            perturb_var(temp, sm, grid, mat, sm->elem1d_physics, j, nodes_on_ele, nvar_ele, sm->elem1d_vars[j],var_code, nsubMods, k, grid->elem1d[j].nodes, DEBUG);
        }
        //store in global using the fmap
        //temp has local elemental matrix 
       
        load_global_mat(sm->diagonal, sm->matrix, temp, nodes_on_ele, grid->elem1d[j].nodes, fmap, nvar_ele, var_ele, sm->nodal_nvars, sm->nodal_vars);

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
void local_dofs_to_global_dofs(int *global_dofs,int ndofs,int *dofs,int *local_range,int *ghosts){
    
    int i;
    int local_size = local_range[1]-local_range[0];


    for(i=0;i<ndofs;i++){
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
void load_global_mat(double *vals, int *indptr, int *indices, double **elem_mat, int nodes_on_ele, int *GnodeIDs, int *fmap, int elem_nvars, int *elem_vars, int *nodal_nvars, int **nodal_vars){
    int i,j,k,l,m,n,p;
    int elem_row,elem_col;
    int save_m,save_n;
    //need local to global mapping
    //assume we have at least global row start that this PE owns
    int row_start;

    /// assembles global residual
    // i and j are the row number
    for (i=0; i<nodes_on_ele; i++){
        //get current nodal row vars
        current_nodal_row_vars = nodal_vars[GnodeIDs[i]];
        for (j=0;j<elem_nvars;j++){
            elem_row = i*elem_nvars + j;
            m=0;
            notFound_row=TRUE;
            current_var_row = elem_vars[j];
            while (notFound_row){
                if(current_var_row == current_nodal_row_vars[m]){
                    notFound_row=FALSE;
                    save_m = m;
                    }
                m+=1;
            }
            //local row index to go in
            local_row_index = fmap[GnodeIDs[i]]+save_m;


            //k and l are column number
            for(k=0;k<nodes_on_ele;k++){
                current_nodal_col_vars = nodal_vars[GnodeIDs[k]];
                for(l=0;l<elem_nvars;l++){
                    elem_col = k*elem_nvars + l;
                    //map the current var from the residual to the correct var number in global residual
                    n=0;
                    notFound_col=TRUE;
                    current_var_col = elem_vars[l];
                    while (notFound_col){
                        if(current_var_col == current_nodal_col_vars[m]){
                            notFound=FALSE;
                            save_n = n;
                        }
                    n+=1;
                    }
                    //local column index to go in
                    local_col_index = fmap[GnodeIDs[k]] + save_n;

                    
                    //both conditions should be true or both should be false I think
                    //given a local row and column number, fit into CSR
                    index_start = indptr[row_index]
                    index_end = indptr[row_index+1]
                    //use these to find correct column index
                    for(p=0;p<index_end-index_start;p++){

                    }

                    


                }
            }

        }
    }

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
void perturb_var(double *elem_mat, SSUPER_MODEL *sm, SGRID *grid, SMAT *mat, SELEM_PHYSICS **elem_physics, 
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
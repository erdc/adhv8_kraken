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
    int nvar_ele,nnode, var_code;
    int max_elem_dofs = MAX_NVAR*MAX_NNODE; 
    double temp[max_elem_dofs,max_elem_dofs];

    //zero out stuff
#ifdef _PETSC
        // Zero the PETSc matrix
        ierr = MatZeroEntries(sm->A);   // TODO: DO WE NEED THIS - SAM
#else
        init_adh_matrix(grid->nnodes, mod->max_nsys_sq, sm->matrix, sm->diagonal);
#endif
    

    //loop through all nelem3d
    for (j=0;j<grid->nelem3d;j++){

            //pull all global information to local memory
            nvar_ele = grid->elem3d[j].nvars;
            var_ele = sm->elem3d_vars[j]
            nnode = grid->elem3d[j].nnodes;
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
                perturb_var(temp, sm, grid, mat, 3, j, nodes_on_element, nvar_ele, var_code, nsubMods, k, DEBUG);

            }
            //store in global using the fmap
            //temp has local elemental matrix 
#ifdef _PETSC
            //in progress
#else           
            load_global_mat(sm->diagonal, sm->matrix, temp, nnode, fmap, nvar_ele, var_ele, sm->nodal_nvars, sm->nodal_vars);
#endif
    }

    //loop through all nelem2d
    for (j=0;j<grid->nelem2d;j++){
        nvar_ele = grid->elem2d[j].nvars;
        nnode = grid->elem2d[j].nnodes;
        nsubMods = sm->nSubMods2d[j];
        //needs to be for 2d matrix
        sarray_init_dbl(temp,maxe_elem_dofs*max_elem_dofs);
        for (k=0;k<sm.nSubMods2d[j];k++){
            //use var code like PERTURB_H ,. ...
            var_code = sm->elem2d_vars[j][k];
            //want to loop over variables not necessarily submodels right?
            //need to think about this more
            //like sw2 is vector equation so that we would need 3 perturbations
            perturb_var(temp, sm, grid, mat, 2, j, nodes_on_element, nvar_ele,  var_code, nsubMods, k, DEBUG);
        }
        //store in global using the fmap
        //temp has local elemental matrix 
#ifdef _PETSC
            //in progress
#else           
            load_global_mat(sm->diagonal, sm->matrix, temp, nnode, fmap, nvar_ele, var_ele, sm->nodal_nvars, sm->nodal_vars);
#endif
    }
    //loop through all nelem1d
    for (j=0;j<grid->nelem1d;j++){
        nvar_ele = grid->elem1d[j].nvars;
        nnode = grid->elem1d[j].nnodes;
        nsubMods = sm->nSubMods1d[j];
        //needs to be for 2d matrix
        sarray_init_dbl(temp,maxe_elem_dofs*max_elem_dofs);
        for (k=0;k<sm.nSubMods1d[j];k++){
            //use var code like PERTURB_H ,. ...
            var_code = sm->elem2d_vars[j][k];
            //want to loop over variables not necessarily submodels right?
            //need to think about this more
            //like sw2 is vector equation so that we would need 3 perturbations
            perturb_var(temp, sm, grid, mat, 2, j, nodes_on_element, nvar_ele,  var_code, nsubMods, k, DEBUG);
        }
        //store in global using the fmap
        //temp has local elemental matrix 
#ifdef _PETSC
        //in progress
#else           
        load_global_mat(sm->diagonal, sm->matrix, temp, nnode, fmap, nvar_ele, var_ele, sm->nodal_nvars, sm->nodal_vars);
#endif
    }    
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file maps a local, dense matrix from an element and put into global sparse matrix
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
 *  @param[in] int *nodal_nvars - an array of length nodal_nvars*nnodes that contains the variable codes
 * \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void load_global_mat(double *diagonal, SPARSE_VECT *matrix, double **elem_mat, int nodes_on_ele, int *fmap, int elem_nvars, int *elem_vars, int *nodal_nvars, int **nodal_vars){
    int i,j,k,l,m,n;
    int elem_row,elem_col;
    int save_m,save_n;
    /// assembles global residual
    // i and j are the row number
    for (i=0; i<nodes_on_ele; i++){
        //get current nodal row vars
        current_nodal_row_vars = nodal_vars[i];
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
            row_index = fmap[node_no[i]]+save_m;


            //k and l are column number
            for(k=0;k<nodes_on_ele;k++){
                current_nodal_col_vars = nodal_vars[k];
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
                    
                    col_index = fmap[node_no[k]] + save_n;

                    //both conditions should be true or both should be false I think
                    //maybe check when we run a test
                    if  (elem_row==elem_col ){
                        diagonal[row_index]+=elem_mat[elem_row,elem_col];
                    }else{
                        matrix[??]+=elem_mat[elem_row,elem_col];
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
 *  \author    Mark Loveland
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
void perturb_var(double *elem_mat, SSUPER_MODEL *sm, SGRID *grid, SMAT *mat, int dim, int ie, int nodes_on_element, int nvar_ele,  int var_code, int nsubModels, int dof_no, int DEBUG) {
    
#ifdef _DEBUG
    assert(dim == 1 || dim == 2 || dim==3);
#endif
    
    int i,j, GlobalNodeID = UNSET_INT;
    double epsilon = 0., epsilon2 = 0.;
    

    //pertrubation may depend on physics?
    double perturbation = sqrt(SMALL);


    double temp_P[MAX_NVAR*MAX_NNODE];
    double temp_M[MAX_NVAR*MAX_NNODE];
    double elem_rhs_P[MAX_NVAR*MAX_NNODE];    // +head perturbation initialized by elem_resid
    double elem_rhs_M[MAX_NVAR*MAX_NNODE];    // -head perturbation initialized by elem_resid
    int eq_var_code,eq_var_code2;

    //maybe pass as argument instead?
    int **elem_vars;
    
    //DEBUG = ON;

    if (dim == 1) {
         SELEM_1D elem;
         elem = grid->elem1d[ie];
         elem_vars = sm->elem1d_vars[ie];
    }
    else if (dim == 2) {
        SELEM_2D elem;
        elem = grid->elem2d[ie];
        elem_vars = sm->elem2d_vars[ie];
    }
    else if (dim==3) {
        SELEM_3D elem;
        elem = grid->elem3d[ie];
        elem_vars = sm->elem3d_vars[ie];
    }


    

    //Mark add switch cases for all sol_var
    //switch case for which sm->sol_var
    double temp_sol[nnode];
    if (var_code == PERTURB_H){
        temp_sol = sm->head;
    }else if (var_code == PERTURB_U){
        temp_sol = sm->vel2d[].x;
    }else if (var_code == PERTURB_V){
        temp_sol = sm->vel2d.y
    }

    
    for (i=0; i<nodes_on_element; i++) {
        GlobalNodeID = elem.nodes[i];
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
            eq_var_code = sm.elem2d_physics[ie][j]->fe_resid(sm,temp_P,grid,mat,ie, epsilon,i, var_code, +1, DEBUG);
            add_replace_elem_rhs(elem_rhs_P,temp_P,elem_vars,nvar_ele,var_code,nodes_on_element);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of depth ++++++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: body pertubing -h for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon);
#endif
            //fe_sw2_body_resid(mod,elem_rhs_h_M,ie,epsilon,i,PERTURB_H,-1,DEBUG);
            //should always be same as var_code
            eq_var_code2 = sm.elem2d_physics[ie][j]->fe_resid(sm,temp_M,grid,mat,ie, epsilon,i, var_code, -1, DEBUG);
            add_replace_elem_rhs(elem_rhs_M,temp_M,elem_vars,nvar_ele,var_code,nodes_on_element);
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        }
        // calculate local contribution to a local elemental Jacobian matrix++++++++++++++++++++++++++++++
        // fills in one column
        elem_matrix_deriv(i, dof_no, nodes_on_element, nvar_ele, elem_rhs_P, elem_rhs_M, elem_mat, epsilon2); // gradient
    }
    
    //DEBUG =  OFF;
}



/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a 3 DOF first order finite difference of a Jacobi matrix elemental block
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
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
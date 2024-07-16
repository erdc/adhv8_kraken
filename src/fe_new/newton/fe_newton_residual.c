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
void assemble_residual(SSUPER_MODEL *sm, SGRID *grid, SMAT *mat) {
    int j,k;

    //seems like the easiest way?
    //maybe think about this
    double fmap = sm->fmap; 
    //zero out stuff
    #ifdef _PETSC
        ierr = VecZeroEntries(sm->residual);
    #else
        sarray_init_dbl(sm->residual, sm->ndofs);
    #endif

    //create array which is the max_nvar in the supermodel
    //and max_nnode of the grid
    //will need to define these properties later
    int nentry = MAX_NVAR*MAX_NNODE;     
    double elem_rhs[nentry];
    //also create a temporary variable which will recieve the residuals from individual fe_resid routines
    double temp[nentry];
    int nnodes;

    int var_code;

    //loop through all nelem3d
    for (j=0;j<grid->nelem3d;j++){
        sarray_init_dbl(elem_rhs,nentry);
        nnodes = grid->elem3d[j].nnodes;

        for (k=0;k<sm.nSubMods3d[j];k++){
            sarray_init_dbl(temp,nentry);
            //would this work for vector functions?
            //would it be better to put global residual outside inner loop if possible?
            //convention for filling temp will be:
            // for i in node (for j in nvar temp[nnode*i + j] = result)

            //either modify resid routines to return int or add an argument, can decide later
            //var_code will contain ordered digits of each equation, same codes used in elem*_vars[]
            var_code = sm.elem3d_physics[j][k]->fe_resid(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

            //add temp to elem_rhs
            //in order to do this we will need elemental vars and info about fe_resid routine
            add_replace_elem_rhs(elem_rhs,temp,sm->elem3d_vars[j],sm->elem3d_nvars[j],var_code,nnodes);
        }

        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->residual, elem_rhs, nnodes, j, fmap, sm->elem3d_nvars[j], sm->elem3d_vars[j], sm->nodal_nvars, sm->nodal_vars,3);

    }

    //loop through all nelem2d
    for (j=0;j<grid->nelem2d;j++){
        sarray_init_dbl(temp,nentry);
        nnodes = grid->elem2d[j].nnodes;
        for (k=0;k<sm.nSubMods2d[j];k++){

            var_code = sm.elem2d_physics[j][k]->fe_resid(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

            add_replace_elem_rhs(elem_rhs,temp,sm->elem2d_vars[j],sm->elem2d_nvars[j],var_code,nnodes);
        }

        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->residual, elem_rhs, nnodes, j, fmap, sm->elem2d_nvars[j], sm->elem2d_vars[j], sm->nodal_nvars, sm->nodal_vars,2);

    }


    //loop through all nelem1d
    for (j=0;j<grid->nelem1d;j++){
        sarray_init_dbl(temp,max_elem_dofs);
        nnodes = grid->elem1d[j].nnodes;
        for (k=0;k<sm.nSubMods1d[j];k++){

            var_code = sm.elem1d_physics[j][k]->fe_resid(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

            add_replace_elem_rhs(elem_rhs,temp,sm->elem1d_vars[j],sm->elem1d_nvars[j],var_code,nnodes);
        }

        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->residual, elem_rhs, nnodes, j, fmap, sm->elem1d_nvars[j], sm->elem1d_vars[j], sm->nodal_nvars, sm->nodal_vars,1);
    }


    //Dirichlet stuff?    
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
void add_replace_elem_rhs(double *elem_rhs, double *eq_rhs, int *elem_vars, int elem_nvars, int eq_vars, int nnodes){

    //for each node, place the rhs entries of a specific pde residual routine in the correct slots of an elemental rhs

    int inode,eq_var, neq_vars;
    int current_var;
    bool notFound = TRUE;
    int k,save_k;
    //number of digits in eq_vars will be the number of variables in this residual
    neq_vars = count_digits(eq_vars);

    for (inode=0;inode<nnodes,inode++){

        for (eq_var=0;eq_var<neq_vars,eq_var++){

            //start with first digit, then second and so on
            current_var = eq_vars/pow(10,neq_vars-eq_var-1);

            //map the current var from the residual to the correct var number in elem_vars
            k=0;
            notFound=TRUE;
            while (notFound){
                if(current_var == elem_vars[k]){
                    notFound=FALSE;
                    save_k = k;
                }
                k+=1;
            }

            //now we know what the current_var is in the whole element, put those entries into the elem_rhs
            elem_rhs[inode*elem_nvars + save_k] += eq_rhs[inode*neq_vars+eq_var];

        }
    
    } 


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
int count_digits(int num) {
  int count = 0;
  while (num != 0) {
    num /= 10;
    count++;
  }
  return count;
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
void load_global_resid(double *residual, double *elem_rhs, int nnodes, int ie, int nvars, int *fmap, int elem_nvars, int elem_vars, int nodal_nvars, int *nodal_vars, int dim) {
    int i,j,k,save_k,index;
    int current_var;
    bool notFound;


    //get global node numbers
    int node_no[nnodes];

    if(dim==1){
        for (i=0;i<nnodes;i++){
            node_no[i] = grid->elem1d[ie].nodes[i];
        }
    }else if (dim ==2){
        for (i=0;i<nnodes;i++){
            node_no[i] = grid->elem2d[ie].nodes[i];
        }
    }else if (dim==3){
        for (i=0;i<nnodes;i++){
            node_no[i] = grid->elem3d[ie].nodes[i];
        }
    }
    


    /// assembles global residual
    for (i=0; i<nnodes; i++) {
        //not sure if this is correct or not
        current_nodal_vars = nodal_vars[node_no[i]];
        for (j=0; j<elem_nvars; j++){
#ifdef _PETSC
            //Mark need to figure this out
#else
            //map the current var from the residual to the correct var number in global residual
            k=0;
            notFound=TRUE;
            current_var = elem_vars[j];
            while (notFound){
                if(current_var == current_nodal_vars[k]){
                    notFound=FALSE;
                    save_k = k;
                }
                k+=1;
            }
                    
            index = fmap[node_no[i]] + save_k;
            //does minus convention hold for all residuals?
            residual[index]     -= elem_rhs[i*elem_nvars+j];
#endif
                    
        }
    }

}

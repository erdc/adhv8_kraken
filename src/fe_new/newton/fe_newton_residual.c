/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function assembles the global residual
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] SSUPER_MODEL *sm - the super model where the residual resides
 *  @param[in] SGRID *grid - the grid over which the monolithic residual resides
 *  @param[in] SMAT *mat - the materials structure
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void assemble_residual(SSUPER_MODEL *sm, SGRID *grid, SMAT *mat) {
    int j,k;

    //seems like the easiest way?
    //maybe think about this
    //we want/need a local and global mapping (local as in local to PE)
    //maybe encapsulate as a function instead of hardcoding as vector?
    int fmap = sm->fmap;
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
    double eq_rhs[nentry];
    //dofs means to the process in this case, not cell
    //for a given elemental rhs, this will give index local to process where we should put entries
    int dofs[nentry];
    int nnodes;

    int var_code;

    int nvars_elem;
    //loop through all nelem3d
    for (j=0;j<grid->nelem3d;j++){
        sarray_init_dbl(elem_rhs,nentry);
        nnodes = grid->elem3d[j].nnodes;
        //maybe this just comes from mat instead
        nvars_elem = sm->elem3d_nvars[j];
        for (k=0;k<sm.nSubMods3d[j];k++){
            sarray_init_dbl(eq_rhs,nentry);
            //convention for filling temp will be:
            // for i in node (for j in nvar temp[nnode*i + j] = result)

            //either modify resid routines to return int or add an argument, can decide later
            //var_code will contain ordered digits of each equation, same codes used in elem*_vars[]
            //maybe elem_physics comes from mat as well? instead of elem number, get mat number
            var_code = sm.elem3d_physics[j][k]->fe_resid(sm,eq_rhs,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

            //add eq_rhs to elem_rhs
            //in order to do this we will need elemental vars and info about fe_resid routine
            add_replace_elem_rhs(elem_rhs,eq_rhs,nvars_elem,sm->elem3d_vars[j],var_code,nnodes);
        }
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem3d[j].nodes,nvars_elem,sm->elem3d_vars[j],sm->nodal_nvars, sm->nodal_vars);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->residual, elem_rhs, nnodes, nvars_elem, dofs);
    }

    //loop through all nelem2d
    for (j=0;j<grid->nelem2d;j++){
        sarray_init_dbl(elem_rhs,nentry);
        nnodes = grid->elem2d[j].nnodes;
        nvars_elem = sm->elem2d_nvars[j];
        for (k=0;k<sm.nSubMods2d[j];k++){
            sarray_init_dbl(temp,nentry);
            var_code = sm.elem2d_physics[j][k]->fe_resid(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);
            add_replace_elem_rhs(elem_rhs,temp,nvars_elem,sm->elem2d_vars[j],var_code,nnodes);
        }
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        get_cell_dofs(dofs,fmap,nnodes,grid->elem2d[j].nodes,nvars_elem,sm->elem2d_vars[j],sm->nodal_nvars, sm->nodal_vars);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->residual, elem_rhs, nnodes, nvars_elem, dofs);
    }


    //loop through all nelem1d
    for (j=0;j<grid->nelem1d;j++){
        sarray_init_dbl(elem_rhs,nentry);
        nnodes = grid->elem1d[j].nnodes;
        nvars_elem = sm->elem1d_nvars[j];
        for (k=0;k<sm.nSubMods1d[j];k++){
            sarray_init_dbl(temp,nentry);
            var_code = sm.elem1d_physics[j][k]->fe_resid(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);
            add_replace_elem_rhs(elem_rhs,temp,nvars_elem, sm->elem1d_vars[j],var_code,nnodes);
        }
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        get_cell_dofs(dofs,fmap,nnodes,grid->elem1d[j].nodes,nvars_elem,sm->elem1d_vars[j],sm->nodal_nvars, sm->nodal_vars);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->residual, elem_rhs, nnodes, nvars_elem, dofs);
    }


    //Dirichlet stuff?    
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function adds the residual of one PDE into the elemental residual
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  double *elem_rhs array containing the elemental residual
 *  @param[in]      double *eq_rhs array containg the residual from a single PDE routine
 *  @param[in]      int *elem_vars array containg the unique variables defined on an element
 *  @param[in]      int elem_nvars integer that is the length of *elem_nvars
 *  @param[in]      int eq_vars integer that contains the unique variables defined by an equation on each digit
 *  @param[in]      int nnodes  integer that is the number of nodes on one element
 * 
 *  \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void add_replace_elem_rhs(double *elem_rhs, double *eq_rhs, int elem_nvars, int *elem_vars,  int eq_vars, int nnodes){

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
 *  \brief     Short routine that gives the degrees of freedom local to the current process
 *  \\ sepearated as function so that maybe other functions can be constructed in future
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] local_dofs - an array of integers that will give the equation numbers for a
 *  given cell local to the process
 *  @param[in]  int ie - the local (to process) cell number
 *  @param[in]  int * - the map that takes in local node number and produces the start equation number
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void get_cell_dofs(int *local_dofs, int *fmaplocal, int nnodes, int *local_node_ids ,int elem_nvars, int *elem_vars, int *nodal_nvars, int **nodal_vars){

    int i,j,k,save_k,index,temp,ctr;
    int current_var;
    ctr =0;



    for (i=0; i<nnodes; i++){
        //on this node get nodal vars
        //need to iron this out, this is an array of ints
        current_nodal_vars = nodal_vars[local_node_ids[i]];
        temp = fmaplocal[local_node_ids[i]];
        for (j=0; j<elem_nvars; j++){
            //map the current var from the residual to the correct var number in global residual
            k=0;
            notFound=TRUE;
            current_var = elem_vars[j];
            while (notFound){
                //note current_nodal_vars depends on i
                if(current_var == current_nodal_vars[k]){
                    notFound=FALSE;
                    save_k = k;
                }
                k+=1;
            }
            local_dofs[ctr] =  temp + save_k;
            ctr+=1
        }
    }

}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Short routine that couns the number of digits in a given integer
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in] int num - the integer we want to find the number of digits
 *  \returns   int count - the number of digits
 *  \note assumes num > 0
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
 *  \brief     This function takes an elemental residual and loads into gloval residual
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] double *residual - the global residual
 *  @param[in] double *elem_rhs - the local element right-hand-side
 *  @param[in] int nnodes - the number of local nodes on the element
 *  @param[in] int *GnodeIDs - an array of length nnodes with the global node #'s of the local element nodes
 *  @param[in] int *fmap - a map from the specific global node ID to the first supermodel dof on that node
 *  @param[in] int elem_nvars - the number of active variables on the element
 *  @param[in] int *elem_vars - an array of length elem_nvars the codes of the active variable
 *  @param[in] int nodal_nvars - the number of active variables on the node (must always be >= elem_nvars)
 *  @param[in] int *nodal_vars - an array of length nodal_nvars the codes of the active variable on a node
 * \note elem_rhs[0] = x_eq, elem_rhs[1] = y_eq, elem_rhs[2] = c_eq,
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void load_global_resid(double *residual, double *elem_rhs, int nnodes, int elem_nvars, int *local_dofs) {
    int index;
    /// assembles global residual
    for (i=0; i<nnodes*elem_nvars; i++) {

            //map the current var from the residual to the correct var number in processor residual
            //using local_dofs map

            index = local_dofs[i];
            //does minus convention hold for all residuals?
            residual[index]     -= elem_rhs[i];
                    
        }

}

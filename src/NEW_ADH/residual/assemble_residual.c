#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function assembles the global residual
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] SMODEL_SUPER *sm - the super model where the residual resides
 *  @param[in] SGRID *grid - the grid over which the monolithic residual resides
 *  @param[in] SMAT *mat - the materials structure
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void assemble_residual(SMODEL_SUPER *sm, SGRID *grid) {
    int j,k,l,m;

    //seems like the easiest way?
    //maybe think about this
    //we want/need a local and global mapping (local as in local to PE)
    //maybe encapsulate as a function instead of hardcoding as vector?
    int* fmap = sm->dof_map_local;
    //zero out stuff
    sarray_init_dbl(sm->residual, sm->ndofs);
    printf("zeroed residual structure\n");

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
    int node_id, node_mat_id;

    int elem_vars[MAX_NVAR];
    int physics_vars[MAX_NVAR];
    int var_code;
    int nvar_node[MAX_NNODE];

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

    sarray_init_int(elem_vars,MAX_NVAR);

    //printf("Beggining 3d,2d,1d loops\n");

    int nvars_elem, nphysics_models, mat_id, nvar_pde;

    //printf("nelem3d = %d\n",grid->nelems3d);
    //loop through all nelem3d
    for (j=0;j<grid->nelems3d;j++){
        sarray_init_dbl(elem_rhs,nentry);
        nnodes = grid->elem3d[j].nnodes;
        //maybe this just comes from mat instead
        //nvars_elem = sm->elem3d_nvars[j];
        mat_id = sm->elem3d_physics_mat_id[j];
        nvars_elem = sm->elem3d_physics_mat[mat_id].nvar;
        nphysics_models = sm->elem3d_physics_mat[mat_id].nSubmodels;
        sarray_copy_int(elem_vars, sm->elem3d_physics_mat[mat_id].vars, nvars_elem);

        //pull out nodal variables
        //only necessary for CG i believe
        //for (l=0;l<nnodes;l++){
//            node_id = grid->elem3d[j].nodes[l];
//            node_mat_id = sm->node_physics_mat_id[node_id];
//            nvar_node[l] = sm->node_physics_mat[node_mat_id].nvar;
//            for(m=0;m<nvar_node[l];m++){
//                vars_node[l][m] = sm->node_physics_mat[node_mat_id].vars[m];
//            }
//        }


        for (k=0;k<nphysics_models;k++){
            sarray_init_dbl(eq_rhs,nentry);
            
            nvar_pde = sm->elem3d_physics[mat_id][k].nvar;
            sarray_init_int(physics_vars, nvar_pde);
            sarray_copy_int(physics_vars, sm->elem3d_physics[mat_id][k].physics_vars,nvar_pde);
            //convention for filling temp will be:
            // for i in node (for j in nvar temp[nnode*i + j] = result)

            //either modify resid routines to return int or add an argument, can decide later
            //var_code will contain ordered digits of each equation, same codes used in elem*_vars[]
            //maybe elem_physics comes from mat as well? instead of elem number, get mat number
            var_code = sm->elem3d_physics[mat_id][k].fe_resid(sm,eq_rhs,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

            //add eq_rhs to elem_rhs
            //in order to do this we will need elemental vars and info about fe_resid routine
            add_replace_elem_rhs(elem_rhs,eq_rhs,nvars_elem,elem_vars,nvar_pde,physics_vars,nnodes,-1.0);
        }
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem3d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //get_cell_dofs_2(dofs,nnodes,grid->elem3d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->residual, elem_rhs, nnodes, nvars_elem, dofs);
    }
    //printf("Passed 3d grid\n");
    
    //printf("number of 2d elements : %d\n",grid->nelems2d);
    //loop through all nelem2d, same thing but difficult to write one function because elems are different structures
    for (j=0;j<grid->nelems2d;j++){
        //printf("2d element %d\n",j);
        sarray_init_dbl(elem_rhs,nentry);
        nnodes = grid->elem2d[j].nnodes;
        //printf("2d element nodes = {%d,%d,%d}\n",grid->elem2d[j].nodes[0],grid->elem2d[j].nodes[1],grid->elem2d[j].nodes[2]);
        //maybe this just comes from mat instead
        //nvars_elem = sm->elem3d_nvars[j];


        mat_id = sm->elem2d_physics_mat_id[j];

        nvars_elem = sm->elem2d_physics_mat[mat_id].nvar;
        //printf("NVARS ELEM %d\n",nvars_elem);
        nphysics_models = sm->elem2d_physics_mat[mat_id].nSubmodels;
        sarray_copy_int(elem_vars, sm->elem2d_physics_mat[mat_id].vars, nvars_elem);
        //pull out nodal variables
        //only necessary for CG i believe

//        for (l=0;l<nnodes;l++){
//            node_id = grid->elem2d[j].nodes[l];
//            
//            node_mat_id = sm->node_physics_mat_id[node_id];
//            nvar_node[l] = sm->node_physics_mat[node_mat_id].nvar;
//            for(m=0;m<nvar_node[l];m++){
//                vars_node[l][m] = sm->node_physics_mat[node_mat_id].vars[m];
//            }
//        }

        for (k=0;k<nphysics_models;k++){
            sarray_init_dbl(eq_rhs,nentry);
            
            nvar_pde = sm->elem2d_physics[mat_id][k].nvar;
            sarray_init_int(physics_vars, nvar_pde);
            sarray_copy_int(physics_vars, sm->elem2d_physics[mat_id][k].physics_vars,nvar_pde);
            //convention for filling temp will be:
            // for i in node (for j in nvar temp[nnode*i + j] = result)

            //either modify resid routines to return int or add an argument, can decide later
            //var_code will contain ordered digits of each equation, same codes used in elem*_vars[]
            //maybe elem_physics comes from mat as well? instead of elem number, get mat number
            var_code = sm->elem2d_physics[mat_id][k].fe_resid(sm,eq_rhs,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);
            //printf("Element: %d, Physics module: %d, Physics Residual = {%f,%f,%f,%f,%f,%f,%f,%f,%f}\n",j,k,eq_rhs[0],eq_rhs[1],eq_rhs[2],eq_rhs[3],eq_rhs[4],eq_rhs[5],eq_rhs[6],eq_rhs[7],eq_rhs[8]);
            //add eq_rhs to elem_rhs
            //in order to do this we will need elemental vars and info about fe_resid routine
            add_replace_elem_rhs(elem_rhs,eq_rhs,nvars_elem,elem_vars,nvar_pde,physics_vars,nnodes,-1.0);
            //printf("Element: %d, Physics module: %d, Elemental Residual = {%f,%f,%f,%f,%f,%f,%f,%f,%f}\n",j,k,elem_rhs[0],elem_rhs[1],elem_rhs[2],elem_rhs[3],elem_rhs[4],elem_rhs[5],elem_rhs[6],elem_rhs[7],elem_rhs[8]);
        }
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem2d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //alternative, computes fmaplocal instead of storing
        //get_cell_dofs_2(dofs, nnodes, grid->elem2d[j].nodes ,nvars_elem, elem_vars, sm->node_physics_mat, sm->node_physics_mat_id);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->residual, elem_rhs, nnodes, nvars_elem, dofs);
    }


    for (j=0;j<grid->nelems1d;j++){
        sarray_init_dbl(elem_rhs,nentry);
        nnodes = grid->elem1d[j].nnodes;
        //maybe this just comes from mat instead
        //nvars_elem = sm->elem3d_nvars[j];
        mat_id = sm->elem1d_physics_mat_id[j];
        nvars_elem = sm->elem1d_physics_mat[mat_id].nvar;
        nphysics_models = sm->elem1d_physics_mat[mat_id].nSubmodels;
        sarray_copy_int(elem_vars, sm->elem1d_physics_mat[mat_id].vars, nvars_elem);
        //pull out nodal variables
        //only necessary for CG i believe
//        for (l=0;l<nnodes;l++){
//            node_id = grid->elem1d[j].nodes[l];
//            node_mat_id = sm->node_physics_mat_id[node_id];
//            nvar_node[l] = sm->node_physics_mat[node_mat_id].nvar;
//            for(m=0;m<nvar_node[l];m++){
//                vars_node[l][m] = sm->node_physics_mat[node_mat_id].vars[m];
//            }
//        }

        for (k=0;k<nphysics_models;k++){
            sarray_init_dbl(eq_rhs,nentry);
            
            nvar_pde = sm->elem1d_physics[mat_id][k].nvar;
            sarray_init_int(physics_vars, nvar_pde);
            sarray_copy_int(physics_vars, sm->elem1d_physics[mat_id][k].physics_vars,nvar_pde);
            //convention for filling temp will be:
            // for i in node (for j in nvar temp[nnode*i + j] = result)

            //either modify resid routines to return int or add an argument, can decide later
            //var_code will contain ordered digits of each equation, same codes used in elem*_vars[]
            //maybe elem_physics comes from mat as well? instead of elem number, get mat number
            var_code = sm->elem1d_physics[mat_id][k].fe_resid(sm,eq_rhs,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

            //add eq_rhs to elem_rhs
            //in order to do this we will need elemental vars and info about fe_resid routine
            add_replace_elem_rhs(elem_rhs,eq_rhs,nvars_elem,elem_vars,nvar_pde,physics_vars,nnodes,-1.0);
        }
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        get_cell_dofs(dofs,fmap,nnodes,grid->elem1d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //get_cell_dofs_2(dofs,nnodes,grid->elem1d[j].nodes,nvars_elem,elem_vars,sm->node_physics_mat, sm->node_physics_mat_id);
        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->residual, elem_rhs, nnodes, nvars_elem, dofs);
    }


//    for(j=0;j<MAX_NNODE;j++){
//        free(vars_node[j]);
//    } 
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
void add_replace_elem_rhs(double *elem_rhs, double *eq_rhs, int elem_nvars, int *elem_vars, int eq_nvars, int *eq_vars, int nnodes, double scale){

    //for each node, place the rhs entries of a specific pde residual routine in the correct slots of an elemental rhs

    int inode,eq_var; //eq_nvars;
    int current_var;
    bool notFound = TRUE;
    int k,save_k;
    //number of digits in eq_vars will be the number of variables in this residual
    //eq_nvars = count_digits(eq_vars);

    for (inode=0;inode<nnodes;inode++){

        for (eq_var=0;eq_var<eq_nvars;eq_var++){

            //start with first digit, then second and so on
            current_var = eq_vars[eq_var];
            //current_var = eq_vars/pow(10,eq_nvars-eq_var-1);

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
            elem_rhs[inode*elem_nvars + save_k] += scale*eq_rhs[inode*eq_nvars+eq_var];

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
 *  @param[in] int node_nvars - the number of active variables on the node (must always be >= elem_nvars)
 *  @param[in] int *node_vars - an array of length node_nvars the codes of the active variable on a node
 * \note elem_rhs[0] = x_eq, elem_rhs[1] = y_eq, elem_rhs[2] = c_eq,
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void load_global_resid(double *residual, double *elem_rhs, int nnodes, int elem_nvars, int *local_dofs) {
    int index,i;
    /// assembles global residual
    for (i=0; i<nnodes*elem_nvars; i++) {

            //map the current var from the residual to the correct var number in processor residual
            //using local_dofs map

            index = local_dofs[i];
            //does minus convention hold for all residuals?
            residual[index]     += elem_rhs[i];
                    
        }

}

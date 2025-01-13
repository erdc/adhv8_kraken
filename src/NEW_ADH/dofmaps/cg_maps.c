/*! \file  cg_maps.c This file collections functions responsible for finding order of finite element
 * degrees of freedom for CG elements  */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Routine that gives an array of degrees of freedom local to the current process for a CG element using fmaplocal array
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] local_dofs (int*) - an array of integers that will give the degree of freedom numbers (equation numbers) for a
 *  given element local to the process
 *  @param[in] fmaplocal (int*) - an array of integers that gives the lowest d.o.f at a given node
 *  @param[in] nnodes (int) - the number of nodes on the element
 *  @param[in] local_node_ids (int*) - array of length nnodes containing node numbers local to process
 *  @param[in] elem_nvars (int) - number of solution variables active on the element
 *  @param[in] elem_vars (int*) - array of length elem_nvars that has the integer code for each variable
 *  @param[in] node_physics_mat (SMAT_PHYSICS*) - an array of SMAT_PHYSICS structs that contains variable info for each node
 *  @param[in] nodal_physics_mat_id (int*) - array of integers that gives the nodal physics mat id
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void get_cell_dofs(int *local_dofs, int *fmaplocal, int nnodes, int *local_node_ids ,int elem_nvars, int *elem_vars, SMAT_PHYSICS **node_physics_mat){

    int i,j,k,save_k,ctr, nodeID, current_nodal_nvars,offset, temp_nodal_var;
    int current_var;
#ifdef _DEBUG
    bool isFound;
#endif
    ctr =0;
    save_k=-1;
    //fmap will only work for CG, need to rethink for DG or possibly mixed CG-DG materials

    //fmap will only work for CG, need to rethink for DG or possibly mixed CG-DG materials

    for (i=0; i<nnodes; i++){
        offset=0;
        nodeID=local_node_ids[i];
        //may need to adjust in parallel? shouldnt if node ids local to process
//        for(l=0;l<nodeID;l++){
//            nodal_mat_ID = nodal_physics_mat_id[l];
//            offset+=node_physics_mat[nodal_mat_ID].nvar;
//        }
        offset=fmaplocal[nodeID];
        //now for current node
        //nodal_mat_ID = nodal_physics_mat_id[nodeID];
        //on this node get nodal vars
        //need to iron this out, this is a 2d array of ints , nnodes x nvar
        //int current_nodal_vars = *nodal_vars[local_node_ids[i]];
        //current_nodal_nvars = node_physics_mat[nodal_mat_ID].nvar;
        current_nodal_nvars = node_physics_mat[nodeID]->nvar;
        for (j=0; j<elem_nvars; j++){
            //map the current var from the residual to the correct var number in global residual
#ifdef _DEBUG
            isFound=FALSE;
#endif
            current_var = elem_vars[j];

            //loop through the nodal vars to look for match
            for(k=0;k<current_nodal_nvars;k++){
                temp_nodal_var = node_physics_mat[nodeID]->vars[k];
                //note current_nodal_vars depends on i
                if(current_var == temp_nodal_var){
#ifdef _DEBUG
                    isFound=TRUE;
#endif
                    save_k = k;
                    break;
                }

            }
            //for debugging
#ifdef _DEBUG
            assert(isFound);
#endif
            local_dofs[ctr] =  offset + save_k;
            ctr+=1;
        }
    }

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Routine that gives an array of degrees of freedom local to the current process for a CG element fully implicitly (no fmap array) using redundant calculations
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] local_dofs (int*) - an array of integers that will give the degree of freedom numbers (equation numbers) for a
 *  given element local to the process
 *  @param[in] nnodes (int) - the number of nodes on the element
 *  @param[in] local_node_ids (int*) - array of length nnodes containing node numbers local to process
 *  @param[in] elem_nvars (int) - number of solution variables active on the element
 *  @param[in] elem_vars (int*) - array of length elem_nvars that has the integer code for each variable
 *  @param[in] node_physics_mat (SMAT_PHYSICS*) - an array of SMAT_PHYSICS structs that contains variable info for each node
 *  @param[in] nodal_physics_mat_id (int*) - array of integers that gives the nodal physics mat id
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void get_cell_dofs_2(int *local_dofs, int nnodes, int *local_node_ids ,int elem_nvars, int *elem_vars, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id){

    int i,j,k,l,save_k,ctr,nodal_mat_ID, nodeID, current_var,current_nodal_nvars,offset, temp_nodal_var;
#ifdef _DEBUG
    bool isFound;
#endif
    ctr =0;
    

    //fmap will only work for CG, need to rethink for DG or possibly mixed CG-DG materials

    for (i=0; i<nnodes; i++){
        offset=0;
        nodeID=local_node_ids[i];
        //may need to adjust in parallel? shouldnt if node ids local to process
        for(l=0;l<nodeID;l++){
            nodal_mat_ID = nodal_physics_mat_id[l];
            offset+=node_physics_mat[nodal_mat_ID].nvar;
        }
        //now for current node
        nodal_mat_ID = nodal_physics_mat_id[nodeID];
        //on this node get nodal vars
        //need to iron this out, this is a 2d array of ints , nnodes x nvar
        //int current_nodal_vars = *nodal_vars[local_node_ids[i]];
        current_nodal_nvars = node_physics_mat[nodal_mat_ID].nvar;
        for (j=0; j<elem_nvars; j++){
            //map the current var from the residual to the correct var number in global residual
#ifdef _DEBUG
            isFound=FALSE;
#endif
            current_var = elem_vars[j];

            //loop through the nodal vars to look for match
            for(k=0;k<current_nodal_nvars;k++){
                temp_nodal_var = node_physics_mat[nodal_mat_ID].vars[k];
                //note current_nodal_vars depends on i
                if(current_var == temp_nodal_var){
#ifdef _DEBUG
                    isFound=TRUE;
#endif
                    save_k = k;
                    break;
                }

            }
#ifdef _DEBUG
            assert(isFound);
#endif
            local_dofs[ctr] =  offset + save_k;
            ctr+=1;
        }
    }

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Routine that gives returns the degree of freedom local to the current process for a CG element using fmaplocal array
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in] var (int) - variable code
 *  @param[in] NodeID (int) - node number local to process
 *  @param[in] fmaplocal (int*) - an array of integers that gives the lowest d.o.f at a given node
 *  @param[in] node_physics_mat (SMAT_PHYSICS*) - an array of SMAT_PHYSICS structs that contains variable info for each node
 *  @param[in] nodal_physics_mat_id (int*) - array of integers that gives the nodal physics mat id
 *  \returns int - an integer that will give the degree of freedom number (equation number) for a
 *  given node and variable number
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int get_cg_dof(int var, int NodeID, int *fmaplocal, SMAT_PHYSICS **node_physics_mat){

    int k,save_k, current_nodal_nvars, temp_nodal_var;
#ifdef _DEBUG
    bool isFound;
#endif
    int offset=0;
    save_k=-1;
    
    //fmap will only work for CG, need to rethink for DG or possibly mixed CG-DG materials
    
    //may need to adjust in parallel? shouldnt if node ids local to process
    offset=fmaplocal[NodeID];
    //now current node
    current_nodal_nvars = node_physics_mat[NodeID]->nvar;
    //map the current var and node to the correct equation number in global residual
#ifdef _DEBUG
    isFound=FALSE;
#endif
    //loop through the nodal vars to look for match
    for(k=0;k<current_nodal_nvars;k++){
        temp_nodal_var = node_physics_mat[NodeID]->vars[k];
        if(var == temp_nodal_var){
#ifdef _DEBUG
            isFound=TRUE;
#endif
            save_k = k;
            break;
        }

    }
#ifdef _DEBUG
    assert(isFound);
#endif
    return offset + save_k;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Routine that gives returns the degree of freedom local to the current process for a CG element fully implicitly (not using fmaplocal array)
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in] var (int) - variable code
 *  @param[in] NodeID (int) - node number local to process
 *  @param[in] node_physics_mat (SMAT_PHYSICS*) - an array of SMAT_PHYSICS structs that contains variable info for each node
 *  @param[in] nodal_physics_mat_id (int*) - array of integers that gives the nodal physics mat id
 *  \returns int - an integer that will give the degree of freedom number (equation number) for a
 *  given node and variable number
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int get_cg_dof_2(int var, int NodeID, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id){

    int k,l,save_k, nodal_mat_ID, current_nodal_nvars, temp_nodal_var;
#ifdef _DEBUG
    bool isFound;
#endif
    int offset=0;
    save_k=-1;
    
    //fmap will only work for CG, need to rethink for DG or possibly mixed CG-DG materials
    
    //may need to adjust in parallel? shouldnt if node ids local to process
    for(l=0;l<NodeID;l++){
        nodal_mat_ID = nodal_physics_mat_id[l];
        offset+=node_physics_mat[nodal_mat_ID].nvar;
    }
    //now current node
    nodal_mat_ID = nodal_physics_mat_id[NodeID];
    current_nodal_nvars = node_physics_mat[nodal_mat_ID].nvar;
    //map the current var and node to the correct equation number in global residual
#ifdef _DEBUG
    isFound=FALSE;
#endif
    //loop through the nodal vars to look for match
    for(k=0;k<current_nodal_nvars;k++){
        temp_nodal_var = node_physics_mat[nodal_mat_ID].vars[k];
        if(var == temp_nodal_var){
#ifdef _DEBUG
            isFound=TRUE;
#endif
            save_k = k;
            break;
        }

    }
#ifdef _DEBUG
    assert(isFound);
#endif
    return offset + save_k;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     extracts sub-array of solution values for specific variable at a given set of nodes
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] local (double*) - the array containing local (to the element) values of a specific variable
 *  @param[in] global (double*) - the full array of the solution vector
 *  @param[in] nodes (int*) - array of node IDs local to process
 *  @param[in] nnodes (int) - the numer of nodes in the nodes array
 *  @param[in] var (int) - the variable code to be extracted
 *  @param[in] fmaplocal (int*) - array containing first dof number (local to process) at each node
 *  @param[in] node_physics_mat (SMAT_PHYSICS*) - array of SMAT_PHYSICS structs containing variable info at nodes
 *  @param[in] nodal_physics_mat_id (int*) - array of ints containing the nodal physics mat id at each node
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void global_to_local_dbl_cg(double *local, double *global, int *nodes, int nnodes, int var, int *fmaplocal, SMAT_PHYSICS **node_physics_mat) {
    int i=0;
    int temp;
    for (i=0; i<nnodes; i++) {
        temp = get_cg_dof(var, nodes[i], fmaplocal, node_physics_mat);
        local[i] = global[temp];
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     extracts sub-array of solution values for specific variable at a given set of nodes
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] local (double*) - the array containing local (to the element) values of a specific variable
 *  @param[in] global (double*) - the full array of the solution vector
 *  @param[in] nodes (int*) - array of node IDs local to process
 *  @param[in] nnodes (int) - the numer of nodes in the nodes array
 *  @param[in] var (int) - the variable code to be extracted
 *  @param[in] fmaplocal (int*) - array containing first dof number (local to process) at each node
 *  @param[in] node_physics_mat (SMAT_PHYSICS*) - array of SMAT_PHYSICS structs containing variable info at nodes
 *  @param[in] nodal_physics_mat_id (int*) - array of ints containing the nodal physics mat id at each node
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void global_to_local_SVECT2D_cg(SVECT2D *local, double *global, int *nodes, int nnodes, int varx, int vary, int *fmaplocal, SMAT_PHYSICS **node_physics_mat) {
    int i=0;
    int temp1,temp2;
    for (i=0; i<nnodes; i++) {
        temp1 = get_cg_dof(varx, nodes[i], fmaplocal, node_physics_mat);
        temp2 = get_cg_dof(vary, nodes[i], fmaplocal, node_physics_mat);
        local[i].x = global[temp1];
        local[i].y = global[temp2];
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     extracts sub-array of solution values for specific variable at a given set of nodes without fmaplocal
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] local (double*) - the array containing local (to the element) values of a specific variable
 *  @param[in] global (double*) - the full array of the solution vector
 *  @param[in] nodes (int*) - array of node IDs local to process
 *  @param[in] nnodes (int) - the numer of nodes in the nodes array
 *  @param[in] var (int) - the variable code to be extracted
 *  @param[in] node_physics_mat (SMAT_PHYSICS*) - array of SMAT_PHYSICS structs containing variable info at nodes
 *  @param[in] nodal_physics_mat_id (int*) - array of ints containing the nodal physics mat id at each node
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void global_to_local_dbl_cg_2(double *global, double *local, int *nodes, int nnodes, int var, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id) {
    int i=0;
    int temp;
    for (i=0; i<nnodes; i++) {
        temp = get_cg_dof_2(var, nodes[i], node_physics_mat, nodal_physics_mat_id);
        local[i] = global[temp];
    }
}



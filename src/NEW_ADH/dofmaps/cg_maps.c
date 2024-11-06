#include "adh.h"
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

    int i,j,k,save_k,index,temp,ctr,nvars_node;
    int current_var;
    bool isFound;
    ctr =0;

    //fmap will only work for CG, need to rethink for DG or possibly mixed CG-DG materials

    for (i=0; i<nnodes; i++){
        //on this node get nodal vars
        //need to iron this out, this is a 2d array of ints , nnodes x nvar
        //int current_nodal_vars = *nodal_vars[local_node_ids[i]];
        nvars_node=nodal_nvars[i];
        //an array that takes node id (local to process), and points to first dof for that node
        //(CG map only)
        temp = fmaplocal[local_node_ids[i]];
        for (j=0; j<elem_nvars; j++){
            //map the current var from the residual to the correct var number in global residual

            isFound=FALSE;
            current_var = elem_vars[j];

            //loop through the nodal vars to look for match
            for(k=0;k<nvars_node;k++){
                //note current_nodal_vars depends on i
                if(current_var == nodal_vars[i][k]){
                    isFound=TRUE;
                    save_k = k;
                    break;
                }

            }
            assert(isFound);
            local_dofs[ctr] =  temp + save_k;
            ctr+=1;
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
void get_cell_dofs_2(int *local_dofs, int nnodes, int *local_node_ids ,int elem_nvars, int *elem_vars, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id){

    int i,j,k,l,save_k,index,temp,ctr,nvars_node, nodal_mat_ID, nodeID, current_nodal_nvars,offset, temp_nodal_var;
    int current_var;
    bool isFound;
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
            isFound=FALSE;
            current_var = elem_vars[j];

            //loop through the nodal vars to look for match
            for(k=0;k<current_nodal_nvars;k++){
                temp_nodal_var = node_physics_mat[nodal_mat_ID].vars[k];
                //note current_nodal_vars depends on i
                if(current_var == temp_nodal_var){
                    isFound=TRUE;
                    save_k = k;
                    break;
                }

            }
            assert(isFound);
            local_dofs[ctr] =  offset + save_k;
            ctr+=1;
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
int get_cg_dof(int var, int NodeID, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id){

    int k,l,save_k, nodal_mat_ID, current_nodal_nvars, temp_nodal_var;
    int current_var;
    bool isFound;
    int offset=0;
    
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
    isFound=FALSE;
    //loop through the nodal vars to look for match
    for(k=0;k<current_nodal_nvars;k++){
        temp_nodal_var = node_physics_mat[nodal_mat_ID].vars[k];
        if(var == temp_nodal_var){
            isFound=TRUE;
            save_k = k;
            break;
        }

    }
    assert(isFound);
    return offset + save_k;
}



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



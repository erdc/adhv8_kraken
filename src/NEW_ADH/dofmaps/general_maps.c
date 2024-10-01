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


int global_to_local(int global_col_no,int local_size,int *ghosts, int nghost){
  //searches the array of ghosts and returns the local index
  //assumes on a process, dofs owned by process come first and then ghosts next
  //not assuming ghosts are ordered here
  int i;
  //printf("len of array, %d\n",len_array);
  for (i=0;i<nghost;i++){
    if(ghosts[i] == global_col_no){
      return i + local_size;
    }
  }
  return -1;

}

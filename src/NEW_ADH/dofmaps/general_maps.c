/*! \file  general_maps.c This file collections functions responsible for finding order of finite element
 * degrees of freedom  */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function takes an array of local (to process) degrees of freedom and gives the global d.o.f. numbers
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] global_dofs (int*) - array of global dofs corresponding to the local dofs
 *  @param[in] ndofs_on_ele (int) - number of degrees of freedom to be calculated
 *  @param[in] dofs (int*) - array of dofs local to the process
 *  @param[in] local_range (int*) -  array of 2 integers, first is global dof start and second is global dof end
 *  @param[in] ghosts (int*) -  array of global dof numbers for ghost dofs on the process
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
            //ghost dof (not owned by process)
            global_dofs[i] = ghosts[dofs[i]-local_size];
        }
    }


}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function takes a global degree of freedom and gives the local (to processor) d.o.f. number
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in] global_col_no (int) - global d.o.f. number corresponding to the local dofs
 *  @param[in] local_size (int) - number of degrees of freedom owned by processor
 *  @param[in] ghosts (int*) - array of d.o.f. numbers that aren't owned by processor
 *  @param[in] nghost (int) - number of ghost d.o.f
 *  \returns (int) - degree of freedom number local to process
 * \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
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

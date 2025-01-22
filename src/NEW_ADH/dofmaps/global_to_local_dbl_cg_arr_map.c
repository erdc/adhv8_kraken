/*! \file  global_to_local_dbl_cg_arr_map.c This file collections functions responsible for finding order of finite element
 * degrees of freedom  */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function takes a global array of values, and picks out the values for an array of nodes
 *  given a map array which takes nodeID -> dof index for the specific variable. (works for CG only)
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] local_vals (double*) - array of values corresponding to the node Ids
 *  @param[in] nodeIDs (int*) - node IDs (from the grid)
 *  @param[in] nnodes (int) - number of nodes requested to pull out from global array
 *  @param[in] map_array (int*) -  array that takes node ID and produces dof # for a specific variable
 *  @param[in] global_vals (double*) -  the global array of values that we want to pull from
 * \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void global_to_local_dbl_cg_arr_map(double *local_vals, int *nodeIDs, int nnodes, int *map_array, double *global_vals){
        for(int i =0;i<nnodes;i++){
                local_vals[i] = global_vals[map_array[nodeIDs[i]]];
        }
}

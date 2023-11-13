/* checks comm_update by comparing the ghost values to the actual values in the
   overlap region */

#include "global_header.h"

void comm_check_int(SGRID *grid)
{

  int ii = 0;                   /* loop counter */
  int *v=NULL;
 v = (int *) tl_alloc(sizeof(int), grid->nnodes);
  for (ii = 0; ii < grid->my_nnodes; ii++)
    {
      v[ii] = grid->smpi->myid;
    }
  for (ii = grid->my_nnodes; ii < grid->nnodes; ii++)
    {
      v[ii] = UNSET_INT;
    }
//tl_check_all_pickets(__FILE__,__LINE__);
  comm_update_int(v, 1, grid->smpi);
//  tl_check_all_pickets(__FILE__,__LINE__);
  for (ii = 0; ii < grid->nnodes; ii++)
    {
      if(v[ii] == UNSET_INT) printf("comm int error myid %d ii %d gid %d resident_pe %d \n", grid->smpi->myid, ii, grid->node[ii].gid, grid->node[ii].resident_pe);
}

printf("********myid %d int check done \n", grid->smpi->myid);
v = (int *) tl_free(sizeof(int), grid->nnodes, v);
return;
}

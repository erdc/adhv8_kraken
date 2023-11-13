/* this routine renumbers the tensor */
#include "global_header.h"

void node_renumber_tensor2d(int max_nnode,
                            STENSOR2D * the_array,   /* the array */
                            STENSOR2D * tmp_array,   /* the tmp array */
                            int *new_number /* the new node numbers */
  )
{
  int i;                        /* loop counter */

  /* set the current array to the tmp array */
  for (i = 0; i < max_nnode; i++) {
    tmp_array[i].xx = the_array[i].xx;
    tmp_array[i].xy = the_array[i].xy;
    tmp_array[i].yy = the_array[i].yy;
  }

  /* copies the tmp array back to the array and renumbers */
  for (i = 0; i < max_nnode; i++) {
    the_array[new_number[i]].xx = tmp_array[i].xx;
    the_array[new_number[i]].xy = tmp_array[i].xy;
    the_array[new_number[i]].yy = tmp_array[i].yy;
  }
}

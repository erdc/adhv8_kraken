/* this routine renumbers the VECT2D arrays */
#include "global_header.h"

void node_renumber_vect2d(int max_nnode,
                          SVECT2D * the_array,   /* the array */
                          SVECT2D * tmp_array,   /* the tmp array */
                          int *new_number,   /* the new node numbers */
                          int *order_tmp
  )
{
  int i;                        /* loop counter */

  /* copies the tmp array back to the array and renumbers */
  for (i = 0; i < max_nnode; i++) {
    if(new_number[i] != i) {
      tmp_array[new_number[i]].x = the_array[new_number[i]].x;
      tmp_array[new_number[i]].y = the_array[new_number[i]].y;
      if (order_tmp[i]==0) {
        the_array[new_number[i]].x = the_array[i].x;
        the_array[new_number[i]].y = the_array[i].y;
      }
      else {
        the_array[new_number[i]].x = tmp_array[i].x;
        the_array[new_number[i]].y = tmp_array[i].y;
      }
    }
  }
}

/* this routine renumbers the global node numbers */

#include "global_header.h"

#ifdef _MESSG_ADPT
void node_renumber_gn(int max_nnode,
                      GLOBAL_NODE * the_array,  /* the array */
                      GLOBAL_NODE * tmp_array,  /* the tmp array */
                      int *new_number,   /* the new node numbers */
                      int *order_tmp)
{
  int i;                        /* loop counter */

  /* copies the tmp array back to the array and renumbers */
  for (i = 0; i < max_nnode; i++) {
    if(new_number[i] != i) {
      tmp_array[new_number[i]].sd = the_array[new_number[i]].sd;
      tmp_array[new_number[i]].rnode = the_array[new_number[i]].rnode;
      if (order_tmp[i]==0) {
        the_array[new_number[i]].sd = the_array[i].sd;
        the_array[new_number[i]].rnode = the_array[i].rnode;
      }
      else {
        the_array[new_number[i]].sd = tmp_array[i].sd;
        the_array[new_number[i]].rnode = tmp_array[i].rnode;
      }
    }
  }
}
#else
void node_renumber_gn(void) {}
#endif

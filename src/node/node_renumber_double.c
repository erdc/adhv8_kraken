/* this routine renumbers the double arrays */
#include "global_header.h"

void node_renumber_double(int max_nnode,
                          double *the_array,    /* the array */
                          double *tmp_array,    /* the tmp array */
                          int *new_number,   /* the new node numbers */
                          int *order_tmp
  )
{
  int i;                        /* loop counter */
 
  for (i = 0; i < max_nnode; i++) {
    
    if(new_number[i] != i) { 
      tmp_array[new_number[i]] = the_array[new_number[i]];
      if (order_tmp[i]==0) {
        the_array[new_number[i]] = the_array[i];
      }
      else {
        the_array[new_number[i]] = tmp_array[i];
      }
    }
  }

}


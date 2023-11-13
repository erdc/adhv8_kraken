/* resets the sparse vector */

#include "global_header.h"

void spv_reset(
  SPARSE_VECT * sv		/* the sparse vector */
)
{
  int i = 0;			/* loop counter */

  /* loops over the indices and initializes the values */
  for(i = 0; i < sv->size; i++)
    sv->index[i] = UNSET_INT;
  return;
}

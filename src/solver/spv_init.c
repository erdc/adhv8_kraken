/* initializes the values of the sparse vector */

#include "global_header.h"

void spv_init(SPARSE_VECT sv,   /* the sparse vector */
              int max_nsys_sq)
{
    int i = 0;                  /* loop counter */
    int iend = 0;               /* the end of the loop */

    /* loops over the indices and initializes the values */
    for (i = 0, iend = sv.size * max_nsys_sq; i < iend; i++) {
        sv.value[i] = 0.0;
    }
    return;
}

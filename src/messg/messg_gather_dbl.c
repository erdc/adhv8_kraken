#include "global_header.h"

/* performs an all gather of an array */
void messg_gather_dbl(int root,
                      double *my_x, /* my part of the array */
                      double *x,    /* the resulting array */
                      int size  /* the size of x */
#ifdef _MESSG
                      ,MPI_Comm ADH_COMM
#endif
    )
{
#ifdef _MESSG
    int ierr_code = MPI_ERR_UNKNOWN;    /* the error code from an mpi call */

    ierr_code = MPI_Gather(my_x, size, MPI_DOUBLE, x, size, MPI_DOUBLE, root, ADH_COMM);
    if (ierr_code != MPI_SUCCESS) {
        messg_err(ierr_code);
    }
#else
    int ii;                     /* loop counter */

    /* copies the array */
    for (ii = 0; ii < size; ii++) {
        x[ii] = my_x[ii];
    }
#endif
    return;
}

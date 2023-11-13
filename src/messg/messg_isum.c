#include "global_header.h"

/*!
   \brief Sums group of int values with one contribution from each processor
 */
int messg_isum(int x      /* the value */
#ifdef _MESSG
                  , MPI_Comm ADH_COMM
#endif
  )
{
#ifdef _MESSG
  int x_send = UNSET_INT;        /* the variables for the pass */
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */

  x_send = x;
  ierr_code = MPI_Allreduce(&x_send, &x, 1, MPI_INT, MPI_SUM, ADH_COMM);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }
#endif
  return (x);
}

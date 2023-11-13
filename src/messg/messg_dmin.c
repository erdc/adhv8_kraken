#include "global_header.h"

/*!
   \brief Minimizes a single double value across processors

   \param x Value in
   \return Minimum Value Across Processors
 */
double messg_dmin(double x      /* the value */
#ifdef _MESSG
                  , MPI_Comm ADH_COMM
#endif
  )
{
#ifdef _MESSG
  double x_send = DZERO;        /* the variables for the pass */
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */

  x_send = x;
  ierr_code = MPI_Allreduce(&x_send, &x, 1, MPI_DOUBLE, MPI_MIN, ADH_COMM);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }
#endif
  return (x);
}

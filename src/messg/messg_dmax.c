#include "global_header.h"

/*!
   \brief Maximizes a single double value across processors

   \param x Value in
   \return Maximum Value Across Processors
 */
double messg_dmax(double x      /* the value */
#ifdef _MESSG
                  , MPI_Comm ADH_COMM
#endif
  )
{
#ifdef _MESSG
  double x_send = DZERO;        /* the variables for the pass */
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */

  x_send = x;
  ierr_code = MPI_Allreduce(&x_send, &x, 1, MPI_DOUBLE, MPI_MAX, ADH_COMM);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }
#endif
  return (x);
}

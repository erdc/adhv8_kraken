/*! \file  messg_imin.c This file has function that computes min integer value over all processors*/
#include "adh.h"
/*!
   \brief Minimizes a single integer value across processors
 */
int messg_imin(int x            /* the value */
#ifdef _MESSG
               , MPI_Comm ADH_COMM
#endif
  )
{
#ifdef _MESSG
  int x_send = 0;               /* the variables for the pass */
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */

  x_send = x;
  ierr_code = MPI_Allreduce(&x_send, &x, 1, MPI_INT, MPI_MIN, ADH_COMM);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }
#endif
  return (x);
}

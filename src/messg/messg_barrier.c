#include "global_header.h"

/*!
   \brief Places a barrier in the communications 
 */
void messg_barrier(MPI_Comm ADH_COMM)
{
#ifdef _MESSG
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */

  ierr_code = MPI_Barrier(ADH_COMM);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }
#endif
  return;
}

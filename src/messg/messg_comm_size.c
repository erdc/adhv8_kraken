#include "global_header.h"

/*!
   \brief Get the total number of processors 
 */
#ifdef _MESSG
int messg_comm_size(MPI_Comm ADH_COMM)
#else
int messg_comm_size(void)
#endif
{
  int size = 1;
#ifdef _MESSG
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */

  ierr_code = MPI_Comm_size(ADH_COMM, &size);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }
#endif
  return (size);
}

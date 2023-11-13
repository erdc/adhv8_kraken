#include "global_header.h"

/*!
   \brief Figure out which processor id I am 
 */

#ifdef _MESSG
    int messg_comm_rank(MPI_Comm comm)
#else
    int messg_comm_rank(void)
#endif
{
  int rank = 0;
#ifdef _MESSG
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */

  ierr_code = MPI_Comm_rank(comm, &rank);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }
#endif
  return (rank);
}

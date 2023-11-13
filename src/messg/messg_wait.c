#include "global_header.h"

/*!
   \brief Waits for all the outstanding asynchronous communications to complete 
 */
void messg_wait(SMPI *smpi)
{
#ifdef _MESSG
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */
  /*! Allocate enough memory for the status flags */
  if (smpi->nmsg_counter > smpi->max_nmsg_status)
    {
      smpi->msg_status = (MPI_Status *) tl_realloc(sizeof(MPI_Status), smpi->nmsg_counter, smpi->max_nmsg_status, smpi->msg_status);
      smpi->max_nmsg_status = smpi->nmsg_counter;
    }

  /*! Wait for the messages to come in */
  if (smpi->nmsg_counter > 0)
    {
      ierr_code = MPI_Waitall(smpi->nmsg_counter, smpi->msg_request, smpi->msg_status);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
    }
#endif
  smpi->nmsg_counter = 0;
  smpi->nmsg_request = 0;
  smpi->nmsg_status = 0;
  return;
}

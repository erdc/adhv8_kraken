#include "global_header.h"

/* blocking receive message */
int messg_precv(MESSG_BUFFER * buffer,  /* the message */
                int tag,         /* the message tag */
                SMPI *smpi
  )
{
  int ndat = 0;                 /* the number of data items */
#ifdef _MESSG
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */

  /* allocate more memory for the status flags */
  if (smpi->max_nmsg_status == 0)
    {
      smpi->msg_status = (MPI_Status *) tl_realloc(sizeof(MPI_Status), 1, smpi->max_nmsg_status, smpi->msg_status);
      smpi->max_nmsg_status = 1;
    }

  /* check to see if message from processor has arrived */
  ierr_code = MPI_Probe(buffer->sd, tag, smpi->ADH_COMM, smpi->msg_status);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }

  /* allocate memory for the message and receive the message */
  if (buffer->type == MESSG_INT)
    {
      ierr_code = MPI_Get_count(smpi->msg_status, MPI_INT, &ndat);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
      messg_buffer_alloc(ndat, sizeof(int), buffer);
      ierr_code = MPI_Recv(buffer->buffer, buffer->nitem, MPI_INT, buffer->sd, tag, smpi->ADH_COMM, smpi->msg_status);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
      return (ndat);
    }
  else if (buffer->type == MESSG_DOUBLE)
    {
      ierr_code = MPI_Get_count(smpi->msg_status, MPI_DOUBLE, &ndat);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
      messg_buffer_alloc(ndat, sizeof(double), buffer);
      ierr_code = MPI_Recv(buffer->buffer, buffer->nitem, MPI_DOUBLE, buffer->sd, tag, smpi->ADH_COMM, smpi->msg_status);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
      return (ndat);
    }
  else if (buffer->type == MESSG_PACKED)
    {
      ierr_code = MPI_Get_count(smpi->msg_status, MPI_PACKED, &ndat);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
      messg_buffer_alloc(ndat, ONE, buffer);
      ierr_code = MPI_Recv(buffer->buffer, buffer->nitem, MPI_PACKED, buffer->sd, tag, smpi->ADH_COMM, smpi->msg_status);
      /* reset the pointer for unpacking the messg later */
      buffer->pos = 0;
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
      return (ndat);
    }
  else
    {
      tl_error("Bad type specified in call to messg_precv.");
      return (UNSET_INT);
    }
#endif
  return (ndat);
}

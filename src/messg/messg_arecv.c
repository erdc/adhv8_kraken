#include "global_header.h"

/*!
   \brief Asynchronous receive message 
 */
void messg_arecv(MESSG_BUFFER *buffer, /* the message */
                 int tag_1,        /* the tag */
                 SMPI *smpi
  )
{

  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */
 
  /*! Do Some Checking */
  assert(buffer != NULL);

  /* get more request space if needed */
  if (smpi->nmsg_request == smpi->max_nmsg_request)
    {
      smpi->max_nmsg_request += MESSG_REQ_INC;
      smpi->msg_request = (MPI_Request *) tl_realloc(sizeof(MPI_Request), smpi->max_nmsg_request, smpi->nmsg_request, smpi->msg_request);
    }

  /* receive the message */
  if (buffer->type == MESSG_INT)
    {
      ierr_code = MPI_Irecv(buffer->buffer, buffer->nitem, MPI_INT, buffer->sd, tag_1, smpi->ADH_COMM, smpi->msg_request + smpi->nmsg_request);
    }
  else if (buffer->type == MESSG_DOUBLE)
    {
      ierr_code = MPI_Irecv(buffer->buffer, buffer->nitem, MPI_DOUBLE, buffer->sd, tag_1, smpi->ADH_COMM, smpi->msg_request + smpi->nmsg_request);
    }
  else if (buffer->type == MESSG_PACKED)
    {
      ierr_code = MPI_Irecv(buffer->buffer, buffer->nitem, MPI_PACKED, buffer->sd, tag_1, smpi->ADH_COMM, smpi->msg_request + smpi->nmsg_request);
      buffer->pos = 0;
    }
  else
    {
      tl_error("Bad type specified in call to messg_arecv.");
    }

  /* Check for Errors */
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }

  /* Increment the request and message counter */
  smpi->nmsg_counter++;
  smpi->nmsg_request++;

  return;
}

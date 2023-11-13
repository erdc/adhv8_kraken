#include "global_header.h"

/* unpack a given buffer into a target buffer - this allows messages to be lumped together */
void messg_unpack(MESSG_BUFFER * input, /* the buffer to be packed */
                  MESSG_BUFFER * target, /* the buffer to be packed into */
                  MPI_Comm ADH_COMM
  )
{
#ifdef _MESSG
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */
  int current_position;         /* the current position in the target buffer */
  int ibuff_info[2];            /* buffer that stores information about 
                                   the buffer being packed:  
                                   the type of data and the number of items 
                                   (this must be the same in messg_unpack) */
  char str[MAXLINE];

  /* checks that the input buffer is a packed buffer */
  if (input->type != MESSG_PACKED)
    tl_error("Tried to pack a buffer that was not of type MESSG_PACKED in messg_pack.");

  /* keeps the current position in the target buffer */
  current_position = input->pos;

  /* unpacks ibuff_info */
  ierr_code = MPI_Unpack(input->buffer, input->nitem, &current_position, ibuff_info, 2, MPI_INT, ADH_COMM);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }

  /* sets the type and number of the target buffer */
  target->type = ibuff_info[0];
  target->nitem = ibuff_info[1];

  /* allocates the space in the target buffer and 
     unpacks the packed buffer into the target buffer */
  if (target->type == MESSG_INT)
    {
      messg_buffer_alloc(target->nitem, sizeof(int), target);
      ierr_code = MPI_Unpack(input->buffer, input->nitem, &current_position, target->buffer, target->nitem, MPI_INT, ADH_COMM);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
    }
  else if (target->type == MESSG_DOUBLE)
    {
      messg_buffer_alloc(target->nitem, sizeof(double), target);
      ierr_code = MPI_Unpack(input->buffer, input->nitem, &current_position, target->buffer, target->nitem, MPI_DOUBLE, ADH_COMM);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
    }
  else if (target->type == MESSG_PACKED)
    tl_error("Tried to unpack a packed buffer from a packed buffer in messg_unpack.");
  else{
    sprintf(str, "Bad type specified in call to messg_unpack. Type %d", target->type);
    tl_error(str);
  }

  /* reset pos in the input buffer */
  input->pos = current_position;
#endif
  return;
}

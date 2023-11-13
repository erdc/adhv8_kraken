#include "global_header.h"

/* pack a given buffer into a target buffer - this allows messages to be lumped together */
void messg_pack(MESSG_BUFFER * input,   /* the buffer to be packed */
                MESSG_BUFFER * target,   /* the buffer to be packed into */
                MPI_Comm ADH_COMM
  )
{
#ifdef _MESSG
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */
  int current_position = 0;     /* the current position in the target buffer */
  int input_size = 0;           /* the size of the input data */
  int info_size = 0;            /* the size of the info data */
  int ibuff_info[2];            /* buffer that stores information about 
                                   the buffer being packed:  
                                   the type of data and the number of items 
                                   (this must be the same in messg_unpack) */

  /* checks that the target buffer is a packed buffer */
  assert(target->type == MESSG_PACKED);

  /* if this is the first buffer packed in the target, then reset nitem */
  if (target->pos == 0)
    {
      target->nitem = 0;
    }

  /* keeps the current position in the target buffer */
  current_position = target->nitem;

  /* sets ibuff_info */
  ibuff_info[0] = input->type;
  ibuff_info[1] = input->nitem;

  /* checks that there is sufficient space in the target buffer */
  if (input->type == MESSG_INT)
    {
      ierr_code = MPI_Pack_size(input->nitem, MPI_INT, ADH_COMM, &input_size);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
    }
  else if (input->type == MESSG_DOUBLE)
    {
      ierr_code = MPI_Pack_size(input->nitem, MPI_DOUBLE, ADH_COMM, &input_size);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
    }
  else if (input->type == MESSG_PACKED)
    {
      tl_error("Tried to pack a packed buffer into a second buffer in messg_pack.");
    }
  else
    {
      tl_error("Bad type specified in call to messg_pack.");
    }

  ierr_code = MPI_Pack_size(2, MPI_INT, ADH_COMM, &info_size);
  if (ierr_code != MPI_SUCCESS)
    {
      {
        messg_err(ierr_code);
      }
    }
  if (target->size < target->nitem + info_size + input_size)
    {
      tl_error("Buffer for packing was not sufficiently large.");
    }

  /* pack the info buffer */
  MPI_Pack(ibuff_info, 2, MPI_INT, target->buffer, target->size, &current_position, ADH_COMM);

  /* pack the buffer */
  if (input->type == MESSG_INT)
    {
      ierr_code = MPI_Pack(input->buffer, input->nitem, MPI_INT, target->buffer, target->size, &current_position, ADH_COMM);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
    }
  else if (input->type == MESSG_DOUBLE)
    {
      ierr_code = MPI_Pack(input->buffer, input->nitem, MPI_DOUBLE, target->buffer, target->size, &current_position, ADH_COMM);
      if (ierr_code != MPI_SUCCESS)
        {
          messg_err(ierr_code);
        }
    }
  else if (input->type == MESSG_PACKED)
    {
      tl_error("Tried to pack a packed buffer into a second buffer in messg_pack.");
    }
  else
    {
      tl_error("Bad type specified in call to messg_pack.");
    }

  /* reset nitem in the target buffer */
  target->nitem = current_position;
  target->pos = current_position;

#endif
  return;
}

#include "global_header.h"

/*! 
   \brief Allocate Enough Space for target to be packed in the future

   Only Allow input to be of type MESSG_INT, MESSG_DOUBLE

   On Entry:
   target->buffer is either allocated or not. If it is allocated,
   then target->size is up to date (i.e., it reflects the correct
   size of the allocated space in bytes). The input buffer has all of
   the information that needs to be packed (i.e., input->nitem 
   entries of type input->type, where input-type is INT/DOUBLE)

   On Exit: 
   target->buffer is re-allocated to reflect the space needed to
   pack in the data from the input buffer.
   target->size is up to date (includes integer buffer with information)
 */
void messg_pack_alloc(MESSG_BUFFER * input, /* the buffer to be packed */
                      MESSG_BUFFER * target, /* the buffer to be packed into */
                      MPI_Comm ADH_COMM
  )
{
#ifdef _MESSG
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */
  int data_size = 0;            /* the size of the input data */
  int information_size = 0;     /* the size of the info data */

  /*! Checks the Tags on both buffers for errors */
  assert(target != NULL);
  assert(input != NULL);
  assert(target->type == MESSG_PACKED);
  assert((input->type == MESSG_DOUBLE) || (input->type == MESSG_INT));

  /*! Get Pack Size of Information */
  ierr_code = MPI_Pack_size(2, MPI_INT, ADH_COMM, &information_size);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }

  /*! Get Pack Size of Data */
  if (input->type == MESSG_INT)
    {
      ierr_code = MPI_Pack_size(input->nitem, MPI_INT, ADH_COMM, &data_size);
    }
  else if (input->type == MESSG_DOUBLE)
    {
      ierr_code = MPI_Pack_size(input->nitem, MPI_DOUBLE, ADH_COMM, &data_size);
    }
  else if (input->type == MESSG_PACKED)
    {
      tl_error("Tried to pack a packed buffer into a second buffer in messg_pack_alloc.");
    }
  else
    {
      tl_error("Bad type specified in call to messg_pack_alloc.");
    }

  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }

  messg_buffer_alloc(target->size + information_size + data_size, ONE, target);
#endif
  return;
}

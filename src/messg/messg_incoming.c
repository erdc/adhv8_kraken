#include "global_header.h"

/*!
   \brief Count the number of messages coming into myid 

   Could be used to simply say "hey, I need to communicate with you"
   Or could be used for a normal All-To-All

   \param messg_flag - messages
   On Entry: Integer to send out to each PE from myid
   messg_flag[ii] goes to processor ii
   On Exit: Integer arriving from each PE
   messg_flag[ii] came from processor ii
 */
int messg_incoming(int *messg_flag,  /* messages to flag : return messages from flag */
                   SMPI *smpi
  )
{
  int num_total_message = 0;    /* Sum of all values */
#ifdef _MESSG
  int ii = 0;                   /* loop counter */
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */
  int *temp_flag = NULL;        /* temporary flag marking processors sending messages */

  /*! Do Some Checking */
  if (messg_flag == NULL)
    tl_error("comm/messg_incoming.c");

  /*! Allocate space for incoming messages */
  temp_flag = (int *) tl_alloc(sizeof(int), smpi->npes);

  /*! Send/Receive the All-To-All data */
  ierr_code = MPI_Alltoall(messg_flag, 1, MPI_INT, temp_flag, 1, MPI_INT, smpi->ADH_COMM);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }

  /*! Copy from Local Array */
  for (ii = 0; ii < smpi->npes; ii++)
    {
      messg_flag[ii] = temp_flag[ii];
      num_total_message += messg_flag[ii];
    }

  /*! Free the local array */
  temp_flag = (int *) tl_free(sizeof(int), smpi->npes, temp_flag);

#endif
  return (num_total_message);
}

#include "share_extern.h"

/*!
   \brief Update the Ghost Node Values (Get Info From Owning Processor)

   *) VECT2D Version with multiple entries
   *) Post Receives, then Post Sends.
 */
void comm_update_VECT2D_2(int in_dimension, /* Number of Values per Node */
                          VECT2D ** v   /* the vector to be updated */
  )
{
#ifdef _MPI
  double *bpntr = NULL;         /* pointer to the buffer location */
  int i_processor = 0, ii = 0, jj = 0, kk = 0;  /* loop counter */
  int icnt = 0;                 /* location in the buffer */
  int nvalues = 2;              /* the number of dimensions */

  /* Checking */
  assert(v != NULL);

  /* Post the Receives */
  for (i_processor = 0; i_processor < npes; i_processor++)
    {
      if ((recv_key[i_processor].end - recv_key[i_processor].begin) > 0)
        {
          messg_buffer_free(recv_msg + i_processor);
          messg_buffer_init(recv_msg + i_processor, i_processor);
          icnt = in_dimension * nvalues * (recv_key[i_processor].end - recv_key[i_processor].begin);
          messg_buffer_alloc(icnt, sizeof(double), recv_msg + i_processor);
          recv_msg[i_processor].type = MESSG_DOUBLE;
          messg_arecv(recv_msg + i_processor, TAG_UPDATE);
        }
    }

  /* Load and send the Send buffers */
  for (i_processor = 0; i_processor < npes; i_processor++)
    {
      if (send_key[i_processor].size > 0)
        {
          messg_buffer_free(send_msg + i_processor);
          messg_buffer_init(send_msg + i_processor, i_processor);
          icnt = in_dimension * nvalues * send_key[i_processor].size;
          messg_buffer_alloc(icnt, sizeof(double), send_msg + i_processor);
          send_msg[i_processor].type = MESSG_DOUBLE;
          bpntr = (double *) send_msg[i_processor].buffer;
          for (ii = 0, jj = 0; ii < send_key[i_processor].size; ii++)
            {
              for (kk = 0; kk < in_dimension; kk++)
                {
                  bpntr[jj++] = v[send_key[i_processor].key[ii]][kk].x;
                  bpntr[jj++] = v[send_key[i_processor].key[ii]][kk].y;
                }
            }
          messg_asend(send_msg + i_processor, TAG_UPDATE);
        }
    }

  /* Wait for the asynchronous communications to clear */
  messg_wait();

  /* Use the Receive buffers to update the ghost node values */
  for (i_processor = 0; i_processor < npes; i_processor++)
    {
      if ((recv_key[i_processor].end - recv_key[i_processor].begin) > 0)
        {
          bpntr = (double *) recv_msg[i_processor].buffer;
          for (ii = recv_key[i_processor].begin, jj = 0; ii < recv_key[i_processor].end; ii++)
            {
              for (kk = 0; kk < in_dimension; kk++)
                {
                  v[ii][kk].x = bpntr[jj++];
                  v[ii][kk].y = bpntr[jj++];
                }
            }
        }
    }
#endif
  return;
}

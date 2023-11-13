#include "global_header.h"

/*!
   \brief Update the Ghost Node Values (Get Info From Owning Processor)

   *) TENSOR2D Version
   *) Post Receives, then Post Sends.
 */
void comm_update_TENSOR2D(STENSOR2D * v,  /* the vector to be updated */
                          SMPI *smpi
  )
{
#ifdef _MESSG
  double *bpntr = NULL;         /* pointer to the buffer location */
  int i_processor = 0, ii = 0, jj = 0;  /* loop counter */
  int icnt = 0;                 /* location in the buffer */
  int nvalues = 3;              /* the number of dimensions */
  int start, end;

  /* Checking */
  assert(v != NULL);

  /* Post the Receives */
  for (i_processor = 0; i_processor < smpi->npes; i_processor++)
    {
      if (smpi->nrecv[i_processor] > 0)
        {
          messg_buffer_free(smpi->recv_msg + i_processor);
          messg_buffer_init(smpi->recv_msg + i_processor, i_processor);
          icnt = nvalues * (smpi->nrecv[i_processor]);
          messg_buffer_alloc(icnt, sizeof(double), smpi->recv_msg + i_processor);
          smpi->recv_msg[i_processor].type = MESSG_DOUBLE;
          messg_arecv(smpi->recv_msg + i_processor, TAG_UPDATE, smpi);
        }
    }

  /* Load and send the Send buffers */
  for (i_processor = 0; i_processor < smpi->npes; i_processor++)
    {
      if (smpi->send_key[i_processor].size > 0)
        {
          messg_buffer_free(smpi->send_msg + i_processor);
          messg_buffer_init(smpi->send_msg + i_processor, i_processor);
          messg_buffer_alloc(nvalues * smpi->send_key[i_processor].size, sizeof(double), smpi->send_msg + i_processor);
          smpi->send_msg[i_processor].type = MESSG_DOUBLE;
          bpntr = (double *) smpi->send_msg[i_processor].buffer;
          for (ii = 0, jj = 0; ii < smpi->send_key[i_processor].size; ii++)
            {
              bpntr[jj++] = v[smpi->send_key[i_processor].key[ii]].xx;
              bpntr[jj++] = v[smpi->send_key[i_processor].key[ii]].xy;
              bpntr[jj++] = v[smpi->send_key[i_processor].key[ii]].yy;
            }
          messg_asend(smpi->send_msg + i_processor, TAG_UPDATE, smpi);
        }
    }

  /* Wait for the asynchronous communications to clear */
  messg_wait(smpi);

  /* Use the Receive buffers to update the ghost node values */
  for (i_processor = 0; i_processor < smpi->npes; i_processor++)
    {
      if (smpi->nrecv[i_processor] > 0)
        {
          bpntr = (double *) smpi->recv_msg[i_processor].buffer;
          start = smpi->recv_init[i_processor];
          end = start + smpi->nrecv[i_processor];
          for (ii = start, jj = 0; ii < end; ii++)
            {
              v[ii].xx = bpntr[jj++];
              v[ii].xy = bpntr[jj++];
              v[ii].yy = bpntr[jj++];
            }
        }
    }
#endif
  return;
}

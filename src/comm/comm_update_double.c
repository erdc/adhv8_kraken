#include "global_header.h"

/*!
   \brief Update the Ghost Node Values (Get Info From Owning Processor)

   *) Double Version
   *) Post Receives, then Post Sends.
   *) This has been updated to allow for the receives to post data 
   directly to the actual location in memory where it will end
   up, instead of doing our own temporary buffer.
 */
void comm_update_double(double *v,  /* the vector to be updated */
                        int nvalues, /* number of values to be updated per ghost node */
                        SMPI *smpi
  )
{
#ifdef _MESSG
  double *bpntr = NULL;         /* pointer to allow for de-referencing */
  int i_processor = 0, ii = 0, jj = 0;  /* loop counter */
  int icnt = 0;                 /* location in the buffer */

  /* Checking */
  assert(v != NULL);
  assert(nvalues > 0);
  /* Post the Receives */
  for (i_processor = 0; i_processor < smpi->npes; i_processor++)
    {
      if (smpi->nrecv[i_processor] > 0)
        {
          
          messg_buffer_free(smpi->recv_msg + i_processor);
          /*recv_msg[i_processor].size = 0; Set Above */
          smpi->recv_msg[i_processor].nitem = nvalues * (smpi->nrecv[i_processor]);
          smpi->recv_msg[i_processor].type = MESSG_DOUBLE;
          smpi->recv_msg[i_processor].sd = i_processor;
          smpi->recv_msg[i_processor].pos = 0;
          smpi->recv_msg[i_processor].buffer = &(v[nvalues * smpi->recv_init[i_processor]]);
          messg_arecv(smpi->recv_msg + i_processor, TAG_UPDATE, smpi);
        }
    }

  /* Load and send the Send buffers */
  for (i_processor = 0; i_processor < smpi->npes; i_processor++)
    {
      if (smpi->send_key[i_processor].size > 0)
        {
          messg_buffer_alloc(nvalues * smpi->send_key[i_processor].size, sizeof(double), smpi->send_msg + i_processor);
          /* send_msg[i_processor].size =; Set Above */
          /* send_msg[i_processor].nitem =; Set Above */
          smpi->send_msg[i_processor].type = MESSG_DOUBLE;
          smpi->send_msg[i_processor].sd = i_processor;
          smpi->send_msg[i_processor].pos = 0;
          bpntr = (double *) smpi->send_msg[i_processor].buffer;
          for (ii = 0, icnt = 0; ii < smpi->send_key[i_processor].size; ii++)
            { 
              for (jj = 0; jj < nvalues; jj++)
                {
                  bpntr[icnt++] = v[smpi->send_key[i_processor].key[ii] * nvalues + jj];
                }
            }
          
          messg_asend(smpi->send_msg + i_processor, TAG_UPDATE, smpi);
        }
    }

  /* Wait for the asynchronous communications to clear */
  messg_wait(smpi);
#endif
  return;
}

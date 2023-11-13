#include "share_extern.h"

/*
   \brief Update "owner_pe" Info in Global_Node entries
 */
void comm_update_GN1(GLOBAL_NODE * v    /* the vector to be updated */
  )
{
#ifdef _MPI
  int i_processor = 0;          /* loop counter over the subdomains */
  int ii = 0, jj = 0;           /* loop counters */
  int local_id_check = 0;       /* Check on Global ID coming out of array send */
  int *bpntr = NULL;            /* pointer to the buffer location */

  /* Post the Receives */
  for (i_processor = 0; i_processor < npes; i_processor++)
    {
      if ((recv_key[i_processor].end - recv_key[i_processor].begin) > 0)
        {
          messg_buffer_free(recv_msg + i_processor);
          messg_buffer_init(recv_msg + i_processor, i_processor);
          messg_buffer_alloc(2 * (recv_key[i_processor].end - recv_key[i_processor].begin), sizeof(int), recv_msg + i_processor);
          recv_msg[i_processor].type = MESSG_INT;
          messg_arecv(recv_msg + i_processor, TAG_UPDATE);
        }
    }

  /* Load and send the send buffers */
  for (i_processor = 0; i_processor < npes; i_processor++)
    {
      if (send_key[i_processor].size != 0)
        {
          messg_buffer_free(send_msg + i_processor);
          messg_buffer_init(send_msg + i_processor, i_processor);
          messg_buffer_alloc(2 * send_key[i_processor].size, sizeof(int), send_msg + i_processor);
          send_msg[i_processor].type = MESSG_INT;
          bpntr = (int *) send_msg[i_processor].buffer;
          for (ii = 0, jj = 0; ii < send_key[i_processor].size; ii++)
            {
              bpntr[jj++] = node_map[send_key[i_processor].key[ii]].owner_pe;
              bpntr[jj++] = node_map[send_key[i_processor].key[ii]].global_id;
            }
          messg_asend(send_msg + i_processor, TAG_UPDATE);
        }
    }

  /* Wait for the asynchronous communications to clear */
  messg_wait();

  /* Update My Ghost Node Information */
  for (i_processor = 0; i_processor < npes; i_processor++)
    {
      if ((recv_key[i_processor].end - recv_key[i_processor].begin) > 0)
        {
          bpntr = (int *) recv_msg[i_processor].buffer;
          for (ii = recv_key[i_processor].begin, jj = 0; ii < recv_key[i_processor].end; ii++)
            {
              node_map[ii].owner_pe = bpntr[jj++];
              local_id_check = bpntr[jj++];
              assert(node_map[ii].global_id == local_id_check);
            }
        }
    }

#endif
  return;
}

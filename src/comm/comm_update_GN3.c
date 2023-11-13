#include "share_extern.h"

/*
   \brief Update "global_id" Info in node_map, node_ladj, and node_adj
 */
void comm_update_GN3(GLOBAL_NODE * v,   /* the vector to be updated */
                     GLOBAL_NODE * v2, GLOBAL_NODE * ladj, GLOBAL_NODE * radj)
{

#ifdef _MPI
  int i_processor = 0;          /* loop counter over the subdomains */
  int ii = 0, jj = 0, kk = 0;   /* loop counters */
  int *bpntr = NULL;            /* pointer to the buffer location */
  int recv_len;                 /* number of entries in recv list */
  int *rpntr = NULL;            /* pointer to map of old global_id to new global_id */

  /* Post the Receives */
  for (i_processor = 0; i_processor < npes; i_processor++)
    {
      messg_buffer_free(recv_msg + i_processor);
      if ((recv_key[i_processor].end - recv_key[i_processor].begin) > 0)
        {
          messg_buffer_init(recv_msg + i_processor, i_processor);
          messg_buffer_alloc(3 * (recv_key[i_processor].end - recv_key[i_processor].begin), sizeof(int), recv_msg + i_processor);
          recv_msg[i_processor].type = MESSG_INT;
          messg_arecv(recv_msg + i_processor, TAG_UPDATE);
        }
    }

  /* Load and send the send buffers */
  for (i_processor = 0; i_processor < npes; i_processor++)
    {
      messg_buffer_free(send_msg + i_processor);
      if (send_key[i_processor].size != 0)
        {
          messg_buffer_init(send_msg + i_processor, i_processor);
          messg_buffer_alloc(3 * send_key[i_processor].size, sizeof(int), send_msg + i_processor);
          send_msg[i_processor].type = MESSG_INT;
          bpntr = (int *) send_msg[i_processor].buffer;
          for (ii = 0, jj = 0; ii < send_key[i_processor].size; ii++)
            {
              bpntr[jj++] = v[send_key[i_processor].key[ii]].owner_pe;
              bpntr[jj++] = v2[send_key[i_processor].key[ii]].global_id;
              bpntr[jj++] = v[send_key[i_processor].key[ii]].global_id;
            }
          messg_asend(send_msg + i_processor, TAG_UPDATE);
        }
    }

  /* Wait for the asynchronous communications to clear */
  messg_wait();

  /* Update My Ghost Node Information */
  for (i_processor = 0; i_processor < npes; i_processor++)
    {
      recv_len = recv_key[i_processor].end - recv_key[i_processor].begin;

      if (recv_len > 0)
        {
          bpntr = (int *) recv_msg[i_processor].buffer;
          rpntr = (int *) tl_alloc(sizeof(int), 2 * recv_len);
          for (ii = recv_key[i_processor].begin, jj = 0, kk = 0; ii < recv_key[i_processor].end; ii++)
            {
              if (v[ii].owner_pe != bpntr[jj])
                {
                  printf("whoops\n");
                }
              /*              assert(v[ii].owner_pe == bpntr[jj]); */
              jj++;

              if ((v[ii].global_id != bpntr[jj]) && (v[ii].global_id != -999))
                {
                  printf("uh oh\n");
                }
              /*              assert((v[ii].global_id == bpntr[jj]) || (v[ii].global_id == -999)); */
              jj++;

              /* update global_id in v (which should be node_map!) */
              rpntr[kk] = v[ii].global_id;
              rpntr[kk + 1] = bpntr[jj];
              kk = kk + 2;

              v[ii].global_id = bpntr[jj];
              jj++;

            }
          /* Now, fix ladj and radj */
          for (ii = 0; ii < nnode; ii++)
            {
              for (kk = 0; kk < recv_len; kk++)
                {
                  if ((ladj[ii].global_id == rpntr[2 * kk]) && (ladj[ii].owner_pe == i_processor))
                    {
                      ladj[ii].global_id = rpntr[2 * kk + 1];
                      break;
                    }
                }
            }

          for (ii = 0; ii < nnode; ii++)
            {
              for (kk = 0; kk < recv_len; kk++)
                {
                  if ((radj[ii].global_id == rpntr[2 * kk]) && (radj[ii].owner_pe == i_processor))
                    {
                      radj[ii].global_id = rpntr[2 * kk + 1];
                      break;
                    }
                }
            }
          rpntr = (int *) tl_free(sizeof(int), 2 * recv_len, rpntr);
        }
    }

#endif
  return;
}

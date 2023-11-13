/* flag the edges for adaption */

#include "global_header.h"
#define TAG_FLAG 771
void comm_flag_edges(EDGE_LIST_ITEM ** edge_hashtab,    /* edge hash table */
                     SGRID* grid,
                     int *nrecv_edge, int *nsend_edge,
                     SMODEL *mod                   
      )
{
#ifdef _MESSG
  int isd;                      /* loop counter over the subdomains */
  int i;                        /* loop counter */
  int *bpntr;                   /* pointer to the buffer location */
  int *nsend_edge_flag = nrecv_edge;    /* send and receive are flipped for flags */
  int *nrecv_edge_flag = nsend_edge;
  MESSG_BUFFER *send_edge_flag_msg = grid->smpi->recv_edge_msg;
  MESSG_BUFFER *recv_edge_flag_msg = grid->smpi->send_edge_msg;
  EDGE_LIST_ITEM ***send_edge_flag_key = grid->smpi->recv_edge_key;
  EDGE_LIST_ITEM ***recv_edge_flag_key = grid->smpi->send_edge_key;
  int npes=grid->smpi->npes;


  /* checks that the buffers are large enough */
  for (isd = 0; isd < npes; isd++)
    {
      messg_buffer_alloc(nsend_edge_flag[isd], sizeof(int), send_edge_flag_msg + isd);
      messg_buffer_alloc(nrecv_edge_flag[isd], sizeof(int), recv_edge_flag_msg + isd);
    }

  /* post the receives */
  for (isd = 0; isd < npes; isd++)
    {
      if (nrecv_edge_flag[isd] != 0)
        {
          recv_edge_flag_msg[isd].type = MESSG_INT;
          messg_arecv(recv_edge_flag_msg + isd, TAG_FLAG, grid->smpi);
        }
    }

  /* load and send the send buffers */
  for (isd = 0; isd < npes; isd++)
    {
      if (nsend_edge_flag[isd] != 0)
        {
          bpntr = (int *) send_edge_flag_msg[isd].buffer;
          for (i = 0; i < nsend_edge_flag[isd]; i++)
            {
              if (send_edge_flag_key[isd][i] == NULL || send_edge_flag_key[isd][i]->new_node == UNSET_INT)
                bpntr[i] = NO;
              else
                bpntr[i] = YES;
            }
          send_edge_flag_msg[isd].type = MESSG_INT;
          messg_asend(send_edge_flag_msg + isd, TAG_FLAG, grid->smpi);
        }
    }

  /* wait for the asynchronous communications to clear */
  messg_wait(grid->smpi);

  /* use the recv buffers to update the ghost node values */
  for (isd = 0; isd < npes; isd++)
    if (nrecv_edge_flag[isd] != 0)
      {
        bpntr = (int *) recv_edge_flag_msg[isd].buffer;
        for (i = 0; i < nrecv_edge_flag[isd]; i++)
          if (recv_edge_flag_key[isd][i] == NULL)
            tl_error("Bad key found in comm_flag_edges.");
          else if (bpntr[i] == YES && recv_edge_flag_key[isd][i]->new_node == UNSET_INT)
            adpt_get_node(mod, recv_edge_flag_key[isd][i], edge_hashtab);

      }

 
#endif
}

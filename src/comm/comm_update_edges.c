/* updates the edges for adaption */

#include "global_header.h"

void comm_update_edges(EDGE_LIST_ITEM ** edge_hashtab,  /* edge hash table */
                     SGRID* grid,
                     int *nrecv_edge, int *nsend_edge,
                     SMODEL *mod
  )
{
#ifdef _MESSG
  int isd;                      /* loop counter over the subdomains */
  int i;                        /* loop counter */
  int *bpntr;                   /* pointer to the buffer location */
	int npes = grid->smpi->npes;


  /* checks that the buffers are large enough */
  for (isd = 0; isd < npes; isd++)
    {
      messg_buffer_alloc(nsend_edge[isd], sizeof(int), grid->smpi->send_edge_msg + isd);
      messg_buffer_alloc(nrecv_edge[isd], sizeof(int), grid->smpi->recv_edge_msg + isd);
    }

  /* post the receives */
  for (isd = 0; isd < npes; isd++)
    if (nrecv_edge[isd] != 0)
      {
        grid->smpi->recv_edge_msg[isd].type = MESSG_INT;
        messg_arecv(grid->smpi->recv_edge_msg + isd, TAG_UPDATE, grid->smpi);
      }

  /* load and send the send buffers */
  for (isd = 0; isd < npes; isd++)
    {
      if (nsend_edge[isd] != 0)
        {
          bpntr = (int *) grid->smpi->send_edge_msg[isd].buffer;
          for (i = 0; i < nsend_edge[isd]; i++)
            {
              if (grid->smpi->send_edge_key[isd][i] == NULL)
                tl_error("Bad key found in comm_update_edges.\n");
               bpntr[i] = grid->smpi->send_edge_key[isd][i]->new_node; 
            }
          grid->smpi->send_edge_msg[isd].type = MESSG_INT;
          messg_asend(grid->smpi->send_edge_msg + isd, TAG_UPDATE, grid->smpi);
        }
    }

  /* wait for the asynchronous communications to clear */
  messg_wait(grid->smpi);

  /* use the recv buffers to update the ghost node values */
  for (isd = 0; isd < npes; isd++)
    if (nrecv_edge[isd] != 0)
      {
        bpntr = (int *) grid->smpi->recv_edge_msg[isd].buffer;
        for (i = 0; i < nrecv_edge[isd]; i++)
          {
            if (grid->smpi->recv_edge_key[isd][i] != NULL && bpntr[i] != UNSET_INT)
              {
                if (grid->smpi->recv_edge_key[isd][i]->new_node == UNSET_INT)
                  {
                    grid->smpi->recv_edge_key[isd][i]->new_node = adpt_get_node(mod,grid->smpi->recv_edge_key[isd][i], edge_hashtab);
                  }
                grid->node[grid->smpi->recv_edge_key[isd][i]->new_node].resident_id = bpntr[i];
              }
          }
      }


#endif
}

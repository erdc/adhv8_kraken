/* communicates the element levels */

#include "global_header.h"
#define TAG_ELEM_LEVELS 775

void comm_elem_levels(ELEM_REF_LIST_ITEM ** send_3d_elem,   /* the linked list of 3D elements to send */
                      ELEM_REF_LIST_ITEM ** send_2d_elem,   /* the linked list of 2D elements to send */
                      ELEM_REF_LIST_ITEM ** send_1d_elem,   /* the linked list of 1D elements to send */
                      NODE_LIST_ITEM ** node_hashtab,   /* node hash table */
                      ELEM3D_LIST_ITEM ** elem3d_hashtab,   /* 3D element hash table */
                      ELEM2D_LIST_ITEM ** elem2d_hashtab,   /* 2D element hash table */
                      ELEM1D_LIST_ITEM ** elem1d_hashtab,   /* 1D element hash table */
                      SGRID *grid                           /* all the node and element information */
  )
{
#ifdef _MESSG
  int flag;                     /* flag indicating 1D, 2D, or 3D elements in message */
  int *ibuffer;                 /* integer form of buffer */
  int i;                        /* loop counter */
  int ielem;                    /* local element number */
  int ipos;                     /* position in a buffer */
  int isd;                      /* loop counter over the processors */
  int ndat;                     /* number of data items within a buffer */
  int nd1, nd2, nd3, nd4;       /* local nodes within an element */
  int lev1, lev2, lev3, lev4;   /* the node levels within an element */
  SNODE gn1, gn2, gn3, gn4;   /* global node numbers passed in a message */
  
  int nsend_3d_elem;            /* number of item sent for 3D element messages */
  int nsend_2d_elem;            /* number of item sent for 2D element messages */
  int nsend_1d_elem;            /* number of item sent for 1D element messages */
  ELEM3D_LIST_ITEM *elem3d_pntr;    /* temporary pointer to a 3D element */
  ELEM2D_LIST_ITEM *elem2d_pntr;    /* temporary pointer to a 2D element */
  ELEM1D_LIST_ITEM *elem1d_pntr;    /* temporary pointer to a 1D element */
  ELEM_REF_LIST_ITEM *elem_pntr;    /* element linked list counter */
  MESSG_BUFFER *send_elem_buff; /* send buffer for outgoing elements */
  int *nsend_elem;              /* number of elements to send */
  MESSG_BUFFER recv_elem_buff;  /* recv buffer for incoming elements */
  int *messg_flag;              /* yes/no flag to indicate if a node message is going 
                                   to the given processor */
  int npes = grid->smpi->npes;

  /* allocate and initialize buffers */
  send_elem_buff = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
  for (isd = 0; isd < npes; isd++)
    {
      messg_buffer_init(send_elem_buff + isd, isd);
      send_elem_buff[isd].type = MESSG_INT;
    }

  messg_buffer_init(&recv_elem_buff, UNSET_INT);
  recv_elem_buff.type = MESSG_INT;

  nsend_elem = (int *) tl_alloc(sizeof(int), npes);
  for (isd = 0; isd < npes; isd++)
    nsend_elem[isd] = 0;

  /* allocate the send buffers */
  for (isd = 0; isd < npes; isd++)
    {
      for (nsend_3d_elem = 0, elem_pntr = send_3d_elem[isd]; elem_pntr != NULL; elem_pntr = elem_pntr->next)
        nsend_3d_elem++;
      for (nsend_2d_elem = 0, elem_pntr = send_2d_elem[isd]; elem_pntr != NULL; elem_pntr = elem_pntr->next)
        nsend_2d_elem++;
      for (nsend_1d_elem = 0, elem_pntr = send_1d_elem[isd]; elem_pntr != NULL; elem_pntr = elem_pntr->next)
        nsend_1d_elem++;
      nsend_elem[isd] = nsend_3d_elem * 13 + nsend_2d_elem * 10 + nsend_1d_elem * 7;
      messg_buffer_alloc(nsend_elem[isd], sizeof(int), send_elem_buff + isd);
    }

  /* count the number of incoming messages */
  messg_flag = (int *) tl_alloc(sizeof(int), npes);
  for (isd = 0; isd < npes; isd++)
    if (nsend_elem[isd] > 0)
      messg_flag[isd] = 1;
    else
      messg_flag[isd] = 0;

  /* send the messages */
  for (isd = 0; isd < npes; isd++)
    if (nsend_elem[isd] != 0)
      {
        ipos = 0;
        ibuffer = (int *) send_elem_buff[isd].buffer;
        for (elem_pntr = send_3d_elem[isd]; elem_pntr != NULL; elem_pntr = elem_pntr->next)
          {
            ibuffer[ipos++] = COMM_3D_ELEM_LEVELS;
            ielem = elem_pntr->ielem;
            nd1 = grid->elem3d[ielem].nodes[0];
            nd2 = grid->elem3d[ielem].nodes[1];
            nd3 = grid->elem3d[ielem].nodes[2];
            nd4 = grid->elem3d[ielem].nodes[3];
            ibuffer[ipos++] = grid->node[nd1].resident_pe;
            ibuffer[ipos++] = grid->node[nd1].resident_id;
            ibuffer[ipos++] = grid->node[nd2].resident_pe;
            ibuffer[ipos++] = grid->node[nd2].resident_id;
            ibuffer[ipos++] = grid->node[nd3].resident_pe;
            ibuffer[ipos++] = grid->node[nd3].resident_id;
            ibuffer[ipos++] = grid->node[nd4].resident_pe;
            ibuffer[ipos++] = grid->node[nd4].resident_id;
            ibuffer[ipos++] = grid->elem3d[ielem].levels[0];
            ibuffer[ipos++] = grid->elem3d[ielem].levels[1];
            ibuffer[ipos++] = grid->elem3d[ielem].levels[2];
            ibuffer[ipos++] = grid->elem3d[ielem].levels[3];
          }
        for (elem_pntr = send_2d_elem[isd]; elem_pntr != NULL; elem_pntr = elem_pntr->next)
          {
            ibuffer[ipos++] = COMM_2D_ELEM_LEVELS;
            ielem = elem_pntr->ielem;
            nd1 = grid->elem2d[ielem].nodes[0];
            nd2 = grid->elem2d[ielem].nodes[1];
            nd3 = grid->elem2d[ielem].nodes[2];
            ibuffer[ipos++] = grid->node[nd1].resident_pe;
            ibuffer[ipos++] = grid->node[nd1].resident_id;
            ibuffer[ipos++] = grid->node[nd2].resident_pe;
            ibuffer[ipos++] = grid->node[nd2].resident_id;
            ibuffer[ipos++] = grid->node[nd3].resident_pe;
            ibuffer[ipos++] = grid->node[nd3].resident_id;
            ibuffer[ipos++] = grid->elem2d[ielem].levels[0];
            ibuffer[ipos++] = grid->elem2d[ielem].levels[1];
            ibuffer[ipos++] = grid->elem2d[ielem].levels[2];
            
          }
        for (elem_pntr = send_1d_elem[isd]; elem_pntr != NULL; elem_pntr = elem_pntr->next)
          {
            ibuffer[ipos++] = COMM_1D_ELEM_LEVELS;
            ielem = elem_pntr->ielem;
            nd1 = grid->elem1d[ielem].nodes[0];
            nd2 = grid->elem1d[ielem].nodes[1];
            ibuffer[ipos++] = grid->node[nd1].resident_pe;
            ibuffer[ipos++] = grid->node[nd1].resident_id;
            ibuffer[ipos++] = grid->node[nd2].resident_pe;
            ibuffer[ipos++] = grid->node[nd2].resident_id;
            ibuffer[ipos++] = grid->elem1d[ielem].levels[0];
            ibuffer[ipos++] = grid->elem1d[ielem].levels[1];
          }
        messg_asend(send_elem_buff + isd, TAG_ELEM_LEVELS, grid->smpi);
      }

  /* receive the messages */

  messg_incoming(messg_flag, grid->smpi);

  for (isd = 0; isd < npes; isd++)
    {
      if (messg_flag[isd] > 0)
        {
          /* receive the next message */
          recv_elem_buff.sd = isd;
          ndat = messg_precv(&recv_elem_buff, TAG_ELEM_LEVELS, grid->smpi);
          ibuffer = (int *) recv_elem_buff.buffer;
          /* unpack the buffer */
          for (i = 0; i < ndat;)
            {
              flag = ibuffer[i++];
              /* 3D elements */
              if (flag == COMM_3D_ELEM_LEVELS)
                {
                  gn1.resident_pe = ibuffer[i++];
                  gn1.resident_id = ibuffer[i++];
                  gn2.resident_pe = ibuffer[i++];
                  gn2.resident_id = ibuffer[i++];
                  gn3.resident_pe = ibuffer[i++];
                  gn3.resident_id = ibuffer[i++];
                  gn4.resident_pe = ibuffer[i++];
                  gn4.resident_id = ibuffer[i++];
                  lev1 = ibuffer[i++];
                  lev2 = ibuffer[i++];
                  lev3 = ibuffer[i++];
                  lev4 = ibuffer[i++];
                  nd1 = node_hash_lookup(gn1, node_hashtab, npes);
                  if (nd1 < 0){
                    tl_error("In comm_elem_levels node0 lookup failed for a 3D element.");}
                  nd2 = node_hash_lookup(gn2, node_hashtab, npes);
                  if (nd2 < 0){
                    tl_error("In comm_elem_levels node1 lookup failed for a 3D element.");}
                  nd3 = node_hash_lookup(gn3, node_hashtab, npes);
                  if (nd3 < 0){
                    tl_error("In comm_elem_levels node2 lookup failed for a 3D element.");}
                  nd4 = node_hash_lookup(gn4, node_hashtab, npes);
                  if (nd4 < 0){
                    tl_error("In comm_elem_levels node3 lookup failed for a 3D element.");}
                  elem3d_pntr = elem3d_hash_lookup(nd1, nd2, nd3, nd4, elem3d_hashtab);
                  if (elem3d_pntr == NULL){
                    tl_error("In comm_elem_levels 3D element lookup failed.");}
                  ielem = elem3d_pntr->ielem;
                  if (grid->elem3d[ielem].levels[0] == UNSET_INT){
                    if (lev1 != UNSET_INT)
                      {
                        if (grid->elem3d[ielem].nodes[0] == nd1)
                          grid->elem3d[ielem].levels[0] = lev1;
                        else if (grid->elem3d[ielem].nodes[0] == nd2)
                          grid->elem3d[ielem].levels[0] = lev2;
                        else if (grid->elem3d[ielem].nodes[0] == nd3)
                          grid->elem3d[ielem].levels[0] = lev3;
                        else if (grid->elem3d[ielem].nodes[0] == nd4)
                          grid->elem3d[ielem].levels[0] = lev4;
                        else
                          tl_error("Local nodes do not correspond to global nodes passed to me in comm_elem_levels for a 3D element.");
                        
                      }
                    else
                    {
                      tl_error("Bad levels found in comm_elem_levels for a 3D element.");
                    }
                  }
                  if (grid->elem3d[ielem].levels[1] == UNSET_INT){
                    if (lev2 != UNSET_INT)
                      {
                        if (grid->elem3d[ielem].nodes[1] == nd1)
                          grid->elem3d[ielem].levels[1] = lev1;
                        else if (grid->elem3d[ielem].nodes[1] == nd2)
                          grid->elem3d[ielem].levels[1] = lev2;
                        else if (grid->elem3d[ielem].nodes[1] == nd3)
                          grid->elem3d[ielem].levels[1] = lev3;
                        else if (grid->elem3d[ielem].nodes[1] == nd4)
                          grid->elem3d[ielem].levels[1] = lev4;
                        else{
                          tl_error("Local nodes do not correspond to global nodes passed to me in comm_elem_levels.");
                        }
                      }
                    else{
                      tl_error("Bad levels found in comm_elem_levels.");
                    }
                  }
                  if (grid->elem3d[ielem].levels[2] == UNSET_INT){
                    if (lev3 != UNSET_INT)
                      {
                        if (grid->elem3d[ielem].nodes[2] == nd1)
                          grid->elem3d[ielem].levels[2] = lev1;
                        else if (grid->elem3d[ielem].nodes[2] == nd2)
                          grid->elem3d[ielem].levels[2] = lev2;
                        else if (grid->elem3d[ielem].nodes[2] == nd3)
                          grid->elem3d[ielem].levels[2] = lev3;
                        else if (grid->elem3d[ielem].nodes[2] == nd4)
                          grid->elem3d[ielem].levels[2] = lev4;
                        else
                          tl_error("Local nodes do not correspond to global nodes passed to me in comm_elem_levels for a 3D element.");
                      }
                  }else
                      tl_error("Bad levels found in comm_elem_levels for a 3D element.");
                  if (grid->elem3d[ielem].levels[3] == UNSET_INT){
                    if (lev4 != UNSET_INT)
                      {
                        if (grid->elem3d[ielem].nodes[3] == nd1)
                          grid->elem3d[ielem].levels[3] = lev1;
                        else if (grid->elem3d[ielem].nodes[3] == nd2)
                          grid->elem3d[ielem].levels[3] = lev2;
                        else if (grid->elem3d[ielem].nodes[3] == nd3)
                          grid->elem3d[ielem].levels[3] = lev3;
                        else if (grid->elem3d[ielem].nodes[3] == nd4)
                          grid->elem3d[ielem].levels[3] = lev4;
                        else
                          tl_error("Local nodes do not correspond to global nodes passed to me in comm_elem_levels for a 3D element.");
                      }
                  }else
                      tl_error("Bad levels found in comm_elem_levels for a 3D element.");
                }

              /* 2D elements */
              if (flag == COMM_2D_ELEM_LEVELS)
                {
                  gn1.resident_pe = ibuffer[i++];
                  gn1.resident_id = ibuffer[i++];
                  gn2.resident_pe = ibuffer[i++];
                  gn2.resident_id = ibuffer[i++];
                  gn3.resident_pe = ibuffer[i++];
                  gn3.resident_id = ibuffer[i++];
                  lev1 = ibuffer[i++];
                  lev2 = ibuffer[i++];
                  lev3 = ibuffer[i++];
                  nd1 = node_hash_lookup(gn1, node_hashtab, npes);
                  if (nd1 < 0)
                    tl_error("In comm_elem_levels node0 lookup failed for a 2D element.");
                  nd2 = node_hash_lookup(gn2, node_hashtab, npes);
                  if (nd2 < 0)
                    tl_error("In comm_elem_levels node1 lookup failed for a 2D element.");
                  nd3 = node_hash_lookup(gn3, node_hashtab, npes);
                  if (nd3 < 0)
                    tl_error("In comm_elem_levels node2 lookup failed for a 2D element.");
                  elem2d_pntr = elem2d_hash_lookup(nd1, nd2, nd3, elem2d_hashtab);
                  if (elem2d_pntr == NULL)
                    tl_error("In comm_elem_levels 2D element lookup failed.");
                  ielem = elem2d_pntr->ielem;
                  
                  if (grid->elem2d[ielem].levels[0] == UNSET_INT){
                    if (lev1 != UNSET_INT)
                      {
                        if (grid->elem2d[ielem].nodes[0] == nd1)
                          grid->elem2d[ielem].levels[0] = lev1;
                        else if (grid->elem2d[ielem].nodes[0] == nd2)
                          grid->elem2d[ielem].levels[0] = lev2;
                        else if (grid->elem2d[ielem].nodes[0] == nd3)
                          grid->elem2d[ielem].levels[0] = lev3;
                        else
                          tl_error("Local nodes do not correspond to global nodes passed to me in comm_elem_levels for a 2D element.");
                      }else{
                    
                      tl_error("Bad levels found in comm_elem_levels for a 2D element.");
                    }
                  }
                  if (grid->elem2d[ielem].levels[1] == UNSET_INT){
                    if (lev2 != UNSET_INT)
                      {
                        if (grid->elem2d[ielem].nodes[1] == nd1)
                          grid->elem2d[ielem].levels[1] = lev1;
                        else if (grid->elem2d[ielem].nodes[1] == nd2)
                          grid->elem2d[ielem].levels[1] = lev2;
                        else if (grid->elem2d[ielem].nodes[1] == nd3)
                          grid->elem2d[ielem].levels[1] = lev3;
                        else
                          tl_error("Local nodes do not correspond to global nodes passed to me in comm_elem_levels for a 2D element.");
                      }
                    else
                      tl_error("Bad levels found in comm_elem_levels for a 2D element.");
                  }
                  if (grid->elem2d[ielem].levels[2] == UNSET_INT){
                    if (lev3 != UNSET_INT)
                      {
                        if (grid->elem2d[ielem].nodes[2] == nd1)
                          grid->elem2d[ielem].levels[2] = lev1;
                        else if (grid->elem2d[ielem].nodes[2] == nd2)
                          grid->elem2d[ielem].levels[2] = lev2;
                        else if (grid->elem2d[ielem].nodes[2] == nd3)
                          grid->elem2d[ielem].levels[2] = lev3;
                        else
                          tl_error("Local nodes do not correspond to global nodes passed to me in comm_elem_levels for a 2D element.");
                      }
                    else
                      tl_error("Bad levels found in comm_elem_levels for a 2D element.");
                  }
                }

              /* 1D elements */
              if (flag == COMM_1D_ELEM_LEVELS)
                {
                  gn1.resident_pe = ibuffer[i++];
                  gn1.resident_id = ibuffer[i++];
                  gn2.resident_pe = ibuffer[i++];
                  gn2.resident_id = ibuffer[i++];
                  lev1 = ibuffer[i++];
                  lev2 = ibuffer[i++];
                  nd1 = node_hash_lookup(gn1, node_hashtab, npes);
                  if (nd1 < 0)
                    tl_error("In comm_elem_levels node0 lookup failed for a 1D element.");
                  nd2 = node_hash_lookup(gn2, node_hashtab, npes);
                  if (nd2 < 0)
                    tl_error("In comm_elem_levels node1 lookup failed for a 1D element.");
                  elem1d_pntr = elem1d_hash_lookup(nd1, nd2, elem1d_hashtab);
                  if (elem1d_pntr == NULL)
                    tl_error("In comm_elem_levels 1D element lookup filed.");
                  ielem = elem1d_pntr->ielem;
                  if (grid->elem1d[ielem].levels[0] == UNSET_INT){
                    if (lev1 != UNSET_INT)
                      grid->elem1d[ielem].levels[0] = lev1;

                  }else{
                      tl_error("Bad levels found in comm_elem_levels for a 1D element.");
                  }
                  if (grid->elem1d[ielem].levels[1] == UNSET_INT){
                    if (lev2 != UNSET_INT)
                      grid->elem1d[ielem].levels[1] = lev2;
                  }else{
                      tl_error("Bad levels found in comm_elem_levels for a 1D element.");
                  }
                }
            }
        }
    }
  /* wait for sends to complete (resets msg_count to zero) */
  messg_wait(grid->smpi);
  

  /* free buffers */
  for (i = 0; i < npes; i++)
    {
      if (send_elem_buff[i].buffer != NULL)
        {
          send_elem_buff[i].buffer = tl_free(ONE, send_elem_buff[i].size, send_elem_buff[i].buffer);
        }
    }
  send_elem_buff = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), npes, send_elem_buff);
  if (recv_elem_buff.buffer != NULL)
    {
      recv_elem_buff.buffer = tl_free(ONE, recv_elem_buff.size, recv_elem_buff.buffer);
    }
  nsend_elem = (int *) tl_free(sizeof(int), npes, nsend_elem);
  messg_flag = (int *) tl_free(sizeof(int), npes, messg_flag);

#endif
}

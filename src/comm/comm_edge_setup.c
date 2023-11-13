/* sets up the keys for the edge update communication
   in the adaption */
#include "global_header.h"
#define INTERIOR 0  /* node owned by this processor which is not
                       a ghost on other processors */
#define BORDER 1    /* node owned by this processor which is a
                       ghost on another processor */
#define GHOST 2 /* node owned by another processor */
#define TAG_TMP_EDGE 770    /* message tag for temporary edge data message */

void comm_edge_setup(EDGE_LIST_ITEM** edge_hashtab,     /* hash table of the edges */
                     SGRID* grid,
                     int *nrecv_edge, int *nsend_edge  )
{
#ifdef _MESSG
  int i, j, k;                  /* loop counters */
  int i_processor;              /* loop counter over the processors */
  int *node_type;               /* label nodes as BORDER, INTERIOR or GHOSTS */
  int *edge_cnt;                /* number of communicating edges owned by given node */
  int nd1, nd2;                 /* nodes on an edge */
  EDGE_LIST_ITEM ***edge_list;     /* list of border edges attached to a given node */
  NODE_LIST_ITEM *node_hashtab[HASHSIZE];   /* hash table for nodes */
  int brdr_node;                /* temp holder for current border node */
  SNODE gn;               /* temp used to look up local node number of global node */
  EDGE_LIST_ITEM *next_edge;       /* pointer to the next edge */
  int *bpntr;                   /* points to the void buffer */
  int nnode = grid->nnodes;
  int npes = grid->smpi->npes;
  int myid = grid->smpi->myid;
  
  /* construct the hash table to convert node numbers
     from global node numbers to local node numbers */
  for (i = 0; i < HASHSIZE; i++)
    node_hashtab[i] = NULL;

  for (i = 0; i < nnode; i++)
    node_hash_add_entry(grid->node[i], node_hashtab, i, npes);

  /* IDENTIFY NODE TYPES */
  /* allocate node_type */
  node_type = (int *) tl_alloc(sizeof(int), nnode);

  /* initialize node types to GHOST or INTERIOR */
  for (i = 0; i < nnode; i++)
    if (grid->node[i].resident_pe == myid)
      node_type[i] = INTERIOR;
    else
      node_type[i] = GHOST;

  /* loop over edges and identify BORDER nodes */
  for (i = 0; i < HASHSIZE; i++)
    {
      next_edge = edge_hashtab[i];
      while (next_edge != NULL)
        {
          nd1 = next_edge->nd1;
          nd2 = next_edge->nd2;
          /* if(grid->node[nd1].resident_pe != myid && grid->node[nd2].resident_pe == myid) */
          if (grid->node[nd1].resident_pe != myid && grid->node[nd2].resident_pe == myid)
            {
              node_type[nd2] = BORDER;
            }
      //    else if (grid->node[nd2].resident_pe != myid && grid->node[nd1].resident_pe == myid)
       //     {
       //       node_type[nd1] = BORDER;
       //     }
          next_edge = next_edge->next;
        }
    }

  /* COUNT COMMUNICATING EDGES */
  /* allocate and initialize edge_cnt */
  edge_cnt = (int *) tl_alloc(sizeof(int), nnode);
  for (i = 0; i < nnode; i++)
    edge_cnt[i] = 0;


  /* loop over the edges and count the border edges
     DEFINITION:  a border edge is an edge for which one node is a border node and
     the other node is on the processor interface
     (either border or ghost - not interior)
   */

    for(i = 0; i < HASHSIZE; i++)
      {
        next_edge = edge_hashtab[i];
        while(next_edge != NULL)
        {
          nd1 = next_edge->nd1;
          nd2 = next_edge->nd2;
          /* is the edge a border edge owned by nd2 */
          if(node_type[nd1] != INTERIOR /* nd1 is a ghost or border node */
             && node_type[nd2] == BORDER /* nd2 is a border node */ )
            {
              /* does the edge belong to nd2 - yes if nd1 > nd2 judged by
                 global node number (grid->node)
                 NOTE:  this means first comparing pe numbers and if the nodes come
                 from the same pe, then comparing relative node numbers
               */
              if(grid->node[nd1].resident_pe > grid->node[nd2].resident_pe ||
                 (grid->node[nd1].resident_pe == grid->node[nd2].resident_pe &&
                  grid->node[nd1].resident_id > grid->node[nd2].resident_id))
                {
                  edge_cnt[nd2]++;
                }
            }
          next_edge = next_edge->next;
        }
   }
  
  
//  for (i = 0; i < HASHSIZE; i++)
//    {
//      next_edge = edge_hashtab[i];
//      while (next_edge != NULL)
//        {
//          nd1 = next_edge->nd1;
//          nd2 = next_edge->nd2;
//
//          if (node_type[nd1] != INTERIOR    /* nd1 is a ghost or border node */
//              && node_type[nd2] == BORDER)  /* nd2 is a border node */
//            {
//              if (grid->node[nd1].resident_pe > grid->node[nd2].resident_pe ||
//                 (grid->node[nd1].resident_id > grid->node[nd2].resident_id && grid->node[nd1].resident_pe == grid->node[nd2].resident_pe))
//                {
//                  edge_cnt[nd2]++;
//                }
//            }
//          else if (node_type[nd2] != INTERIOR   /* nd1 is a ghost or border node */
//                   && node_type[nd1] == BORDER) /* nd2 is a border node */
//            {
//              if (grid->node[nd2].resident_pe > grid->node[nd1].resident_pe ||
//                 (grid->node[nd2].resident_id > grid->node[nd1].resident_id && grid->node[nd2].resident_pe == grid->node[nd1].resident_pe))
//               {
//                  edge_cnt[nd1]++;
//                }
//            }
//          next_edge = next_edge->next;
//        }
//    }

  /* allocate list of edges for each node */
  edge_list = (EDGE_LIST_ITEM ***) tl_alloc(sizeof(EDGE_LIST_ITEM **), nnode);
  for (i = 0; i < nnode; i++)
    {
      if (node_type[i] == BORDER)
        {
          if (edge_cnt[i] != 0)
            {
              edge_list[i] = (EDGE_LIST_ITEM **) tl_alloc(sizeof(EDGE_LIST_ITEM *), edge_cnt[i]);
            }
        }
    }

  /* initialize the edge count */
  for (i = 0; i < nnode; i++)
    edge_cnt[i] = 0;

  /* set list of edges for each node */
  for (i = 0; i < HASHSIZE; i++)
    {
      next_edge = edge_hashtab[i];
      while (next_edge != NULL)
        {
          nd1 = next_edge->nd1;
          nd2 = next_edge->nd2;
          if (node_type[nd1] != INTERIOR    /* nd1 is a ghost or border node */
              && node_type[nd2] == BORDER)  /* nd2 is a border node */
            {
              if (grid->node[nd1].resident_pe > grid->node[nd2].resident_pe ||
                 (grid->node[nd1].resident_id > grid->node[nd2].resident_id && grid->node[nd1].resident_pe == grid->node[nd2].resident_pe))
                {
                  edge_list[nd2][edge_cnt[nd2]] = next_edge;    /* note: this is a pointer to the edge */
                  edge_cnt[nd2]++;
                }
            }
//          else if (node_type[nd2] != INTERIOR   /* nd2 is a ghost or border node */
//                   && node_type[nd1] == BORDER) /* nd1 is a border node */
//            {
//              if (grid->node[nd2].resident_pe > grid->node[nd1].resident_pe ||
//                 (grid->node[nd2].resident_id > grid->node[nd1].resident_id && grid->node[nd2].resident_pe == grid->node[nd1].resident_pe))
//                {
//                  edge_list[nd1][edge_cnt[nd1]] = next_edge;    /* note: this is a pointer to the edge */
//                  edge_cnt[nd1]++;
//                }
//            }
          next_edge = next_edge->next;
        }
    }

  /* update the edge count */
  comm_update_int(edge_cnt, 1, grid->smpi);

  /* SET EDGE KEYS */
  /* loop over the nodes and sum the ghost edges */
  for (i_processor = 0; i_processor < npes; i_processor++)
    nrecv_edge[i_processor] = 0;

  for (i = 0; i < nnode; i++)
    if (node_type[i] == GHOST)
      nrecv_edge[grid->node[i].resident_pe] += edge_cnt[i];

  /* allocate receive buffers */
  for (i_processor = 0; i_processor < npes; i_processor++)
    /* the if statement was added 11 Feb 2000 jph */
    if (nrecv_edge[i_processor] != 0)
      {
        messg_buffer_alloc(nrecv_edge[i_processor] * 4, sizeof(int), grid->smpi->recv_edge_msg + i_processor);
      }
  
  /* receive edge messages */
  for (i_processor = 0; i_processor < npes; i_processor++)
    if (nrecv_edge[i_processor] != 0)
      {
        grid->smpi->recv_edge_msg[i_processor].type = MESSG_INT;
        grid->smpi->recv_edge_msg[i_processor].sd = i_processor;
        messg_arecv(grid->smpi->recv_edge_msg + i_processor, TAG_TMP_EDGE, grid->smpi);

      }

  /* loop over the processors to load and send messages */
  for (i_processor = 0; i_processor < npes; i_processor++)
    {

      /* loop over node keys and sum the border edges */
      nsend_edge[i_processor] = 0;

      for (i = 0; i < grid->smpi->send_key[i_processor].size; i++)
        nsend_edge[i_processor] += edge_cnt[grid->smpi->send_key[i_processor].key[i]];

      /* allocate send_edge_msg and send_edge_key */
      if (nsend_edge[i_processor] != 0)
        {
          messg_buffer_alloc(nsend_edge[i_processor] * 4, sizeof(int), grid->smpi->send_edge_msg + i_processor);
          /* THE NEXT LINE IS MESSING UP EDGE_LIST VALUES */
          grid->smpi->send_edge_key[i_processor] = (EDGE_LIST_ITEM **) tl_alloc(sizeof(EDGE_LIST_ITEM *), nsend_edge[i_processor]);
        }

      /* loop over edges and set send_edge_msg and send_edge_key
         NOTE: we are sending global node numbers. */
      bpntr = (int *) grid->smpi->send_edge_msg[i_processor].buffer;

      for (i = 0, j = 0; i < grid->smpi->send_key[i_processor].size; i++)
        {
          brdr_node = grid->smpi->send_key[i_processor].key[i];
          for (k = 0; k < edge_cnt[brdr_node]; k++)
            {
              bpntr[4 * j] = grid->node[edge_list[brdr_node][k]->nd1].resident_pe;
              bpntr[4 * j + 1] = grid->node[edge_list[brdr_node][k]->nd1].resident_id;
              bpntr[4 * j + 2] = grid->node[edge_list[brdr_node][k]->nd2].resident_pe;
              bpntr[4 * j + 3] = grid->node[edge_list[brdr_node][k]->nd2].resident_id;
              grid->smpi->send_edge_key[i_processor][j] = edge_list[brdr_node][k];
              j++;
            }
        }

      /* allocate recv_edge_msg and recv_edge_key */
      if (nrecv_edge[i_processor] != 0)
        {
          grid->smpi->recv_edge_key[i_processor] = (EDGE_LIST_ITEM **) tl_alloc(sizeof(EDGE_LIST_ITEM *), nrecv_edge[i_processor]);
        }
      /* send send_edge_msg */
      if (nsend_edge[i_processor] != 0)
        {
          grid->smpi->send_edge_msg[i_processor].type = MESSG_INT;
          grid->smpi->send_edge_msg[i_processor].sd = i_processor;
          messg_asend(grid->smpi->send_edge_msg + i_processor, TAG_TMP_EDGE, grid->smpi);
        }
    }


  /* wait for messages to arrive */
  messg_wait(grid->smpi);

  /* look up and set recv_edge_key */
  for (i_processor = 0; i_processor < npes; i_processor++)
    {
      bpntr = (int *) grid->smpi->recv_edge_msg[i_processor].buffer;
      for (i = 0; i < nrecv_edge[i_processor]; i++)
        {
          /* looks the node up in the hash table */
          gn.resident_pe = bpntr[i * 4];
          gn.resident_id = bpntr[i * 4 + 1];
          nd1 = node_hash_lookup(gn, node_hashtab, npes);
          gn.resident_pe = bpntr[i * 4 + 2];
          gn.resident_id = bpntr[i * 4 + 3];
          nd2 = node_hash_lookup(gn, node_hashtab, npes);
          
          /* note: each key value is a pointer to the edge in the hash table */
          /* some edges are received that are not needed, and that have nodes
             that the processor does not see-- thus the check below */
          if (nd1 >= 0 && nd2 >= 0)
            grid->smpi->recv_edge_key[i_processor][i] = edge_hash_lookup(nd1, nd2, edge_hashtab);
          else
            grid->smpi->recv_edge_key[i_processor][i] = NULL;
        }
    }

  /* free memory */

  for (i = 0; i < nnode; i++)
    if (node_type[i] == BORDER)
      if (edge_cnt[i] != 0)
        {
          edge_list[i] = (EDGE_LIST_ITEM **)tl_free(sizeof(EDGE_LIST_ITEM *), edge_cnt[i], edge_list[i]);
        }
  edge_cnt = (int *) tl_free(sizeof(int), nnode, edge_cnt);
  node_type = (int *) tl_free(sizeof(int), nnode, node_type);

  edge_list = (EDGE_LIST_ITEM ***)tl_free(sizeof(EDGE_LIST_ITEM **), nnode, edge_list);
  tl_list_free_all(NODE_LIST);
#endif
  return;
}

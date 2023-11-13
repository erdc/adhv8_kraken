/* communicates the adjacent nodes */

#include "global_header.h"
#define TAG_UNREF 774

#ifdef _MESSG

void comm_adj_nodes(int *node_unref_flags,  /* flags nodes to be eliminated */
                    NODE_LIST_ITEM ** node_hashtab,  /* node hash table */
                    SMODEL *mod
  )
{
  int i, isd, i_processor;                        /* loop counter */
  int *nnode_out;              /* the number of nodes going to each processor */
  int **nodes_out;              /* the nodes going out to each processor */
  int *nnode_in;                /* number of nodes coming into processor */
  MESSG_BUFFER recv_nd_msg;     /* the incoming message buffer */
  MESSG_BUFFER *send_nd_msg;    /* the outgoing message buffers */
  MESSG_BUFFER node_ibuffer;    /* the buffer for the integer node part of the message */
  MESSG_BUFFER node_dbuffer;    /* the buffer for the double node part of the message */
  MESSG_BUFFER node_dbuffer_surface;    /* the buffer for the surface double node part of the message */
  MESSG_BUFFER node_dbuffer_bed;    /* the buffer for the bed double node part of the message */
  int npes = mod->grid->smpi->npes;
  int myid = mod->grid->smpi->myid;
  int nnode = mod->grid->nnodes;
  
  send_nd_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
  /* send adjacent node data */
  nnode_out = (int *)tl_alloc(sizeof(int), npes);
  nodes_out = (int **)tl_alloc(sizeof(int *), npes);
  nnode_in = (int *)tl_alloc(sizeof(int), npes);

  for(i = 0; i < npes; i++)
    nnode_out[i] = 0;
  messg_buffer_init(&recv_nd_msg, UNSET_INT);
  messg_buffer_init(&node_ibuffer, UNSET_INT);
  messg_buffer_init(&node_dbuffer, UNSET_INT);
  if(mod->grid->type == COLUMNAR){
    messg_buffer_init(&node_dbuffer_surface, UNSET_INT);
    messg_buffer_init(&node_dbuffer_bed, UNSET_INT);
  }
  /* set number of nodes to be sent according to receiving processors possibly only my_nnode?*/
  /* potentially expensive cache operation may want to set nodes_out to nnode size to prevent running this loop twice */
  
  
  for (i = 0; i < nnode; i++) {    
    if (node_unref_flags[i] == YES){ 
      if (mod->grid->node[i].parent_res_pe[0] == myid && mod->grid->node[i].parent_res_pe[1] != myid){
        nnode_out[mod->grid->node[i].parent_res_pe[1]]++;
      }
      else if (mod->grid->node[i].parent_res_pe[0] != myid && mod->grid->node[i].parent_res_pe[1] == myid){
        nnode_out[mod->grid->node[i].parent_res_pe[0]]++;
      }
    }
  }
  
 
  for(i = 0; i < npes; i++) {
    if(nnode_out[i] != 0) {
      nodes_out[i] = (int *)tl_alloc(sizeof(int), nnode_out[i]);
    }
  }
  
  for(i = 0; i < npes; i++)
    nnode_out[i] = 0;

  for (i = 0; i < nnode; i++) {
      if (node_unref_flags[i] == YES){
        if (mod->grid->node[i].parent_res_pe[0] == myid && mod->grid->node[i].parent_res_pe[1] != myid){
             i_processor = mod->grid->node[i].parent_res_pe[1];
             
            if((i_processor >= 0) && (i_processor < npes)) nodes_out[i_processor][nnode_out[i_processor]++] = mod->grid->node[i].parent_res_id[0];
            
          }
        else if (mod->grid->node[i].parent_res_pe[0] != myid && mod->grid->node[i].parent_res_pe[1] == myid)
          {
            i_processor = mod->grid->node[i].parent_res_pe[0];
            if((i_processor >= 0) && (i_processor < npes)) nodes_out[i_processor][nnode_out[i_processor]++] = mod->grid->node[i].parent_res_id[1];
            
          }
      }
    }
  

  for (isd = 0; isd < npes; isd++)
    {
      nnode_in[isd] = nnode_out[isd];
      if (nnode_out[isd] > 0)
        {
          /* initialize the auxilliary buffers */
          messg_buffer_free(&node_ibuffer);
          messg_buffer_free(&node_dbuffer);
          if(mod->grid->ndim == 3){
            messg_buffer_free(&node_dbuffer_surface);
            messg_buffer_free(&node_dbuffer_bed);
          }

          messg_buffer_init(&node_ibuffer, UNSET_INT);
          messg_buffer_init(&node_dbuffer, UNSET_INT);
       
          if(mod->grid->ndim ==  3){
            messg_buffer_init(&node_dbuffer_surface, UNSET_INT);
            messg_buffer_init(&node_dbuffer_bed, UNSET_INT);
          }

          /* load the node buffers */
          mesh_pack_node(nnode_out[isd], nodes_out[isd], &node_dbuffer, &node_dbuffer_surface, &node_dbuffer_bed, &node_ibuffer, mod);

          /* set the buffer type */
          messg_buffer_init(&send_nd_msg[isd], isd);
          send_nd_msg[isd].type = MESSG_PACKED;

          /* allocate space to pack the buffer */
          messg_pack_alloc(&node_ibuffer, send_nd_msg + isd, mod->grid->smpi->ADH_COMM);
          messg_pack_alloc(&node_dbuffer, send_nd_msg + isd, mod->grid->smpi->ADH_COMM);
          if(mod->grid->ndim ==  3){
            messg_pack_alloc(&node_dbuffer_surface, send_nd_msg + isd, mod->grid->smpi->ADH_COMM);
            messg_pack_alloc(&node_dbuffer_bed, send_nd_msg + isd, mod->grid->smpi->ADH_COMM);
          }
          messg_pack(&node_ibuffer, send_nd_msg + isd, mod->grid->smpi->ADH_COMM);
          messg_pack(&node_dbuffer, send_nd_msg + isd, mod->grid->smpi->ADH_COMM);
          if(mod->grid->ndim ==  3){
            messg_pack(&node_dbuffer_surface, send_nd_msg + isd, mod->grid->smpi->ADH_COMM);
            messg_pack(&node_dbuffer_bed, send_nd_msg + isd, mod->grid->smpi->ADH_COMM);
          }

          /* send the message */
          messg_asend(send_nd_msg + isd, TAG_NODE_DATA, mod->grid->smpi);

       }
    }
  
  /* receive messages */
  /* loop over the incoming messages */
  recv_nd_msg.type = MESSG_PACKED;
  messg_incoming(nnode_in, mod->grid->smpi);

  for (isd = 0; isd < npes; isd++)
    {

       if (nnode_in[isd] > 0)
        {

          /* initialize the auxilliary buffers */
          /* note:  since these buffers are reused, they must be freed before the buffer pointers
             are reset in messg_buffer_init */
          messg_buffer_free(&node_ibuffer);
          messg_buffer_free(&node_dbuffer);
          if(mod->grid->ndim ==  3){
            messg_buffer_free(&node_dbuffer_surface);
            messg_buffer_free(&node_dbuffer_bed);
          }
          messg_buffer_init(&node_ibuffer, isd);
          messg_buffer_init(&node_dbuffer, isd);
          if(mod->grid->ndim ==  3){
            messg_buffer_init(&node_dbuffer_surface, UNSET_INT);
            messg_buffer_init(&node_dbuffer_bed, UNSET_INT);
          }
          /* get a message */
          recv_nd_msg.sd = isd;
          messg_precv(&recv_nd_msg, TAG_NODE_DATA, mod->grid->smpi);
          /* unpacks the buffers from the recv buffer */
          messg_unpack(&recv_nd_msg, &node_ibuffer, mod->grid->smpi->ADH_COMM);
          messg_unpack(&recv_nd_msg, &node_dbuffer, mod->grid->smpi->ADH_COMM);
          if(mod->grid->ndim ==  3){
            messg_unpack(&recv_nd_msg, &node_dbuffer_surface, mod->grid->smpi->ADH_COMM);
            messg_unpack(&recv_nd_msg, &node_dbuffer_bed, mod->grid->smpi->ADH_COMM);
          }
          /* unpacks the node buffers */
          mesh_unpack_node(&node_dbuffer, &node_dbuffer_surface, &node_dbuffer_bed, &node_ibuffer, mod, node_hashtab, 1);
        }
    }

  /* wait for the asynchronous communication to complete */
  messg_wait(mod->grid->smpi);
  /* free the send buffers */
  for(i = 0; i < npes; i++) {
    if(nnode_out[i] != 0) {
      messg_buffer_free(send_nd_msg + i);
      nodes_out[i] = (int *)tl_free(sizeof(int), nnode_out[i], nodes_out[i]);
    }
  }
  
  
  nnode_in = (int *) tl_free(sizeof(int), npes, nnode_in);
  nnode_out = (int *)tl_free(sizeof(int), npes, nnode_out);
  nodes_out = (int **)tl_free(sizeof(int *), npes, nodes_out);
  send_nd_msg = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), npes, send_nd_msg);
  /* free the local buffer memory */
  messg_buffer_free(&recv_nd_msg);
  messg_buffer_free(&node_ibuffer);
  messg_buffer_free(&node_dbuffer);
  if(mod->grid->ndim ==  3){
    messg_buffer_free(&node_dbuffer_surface);
    messg_buffer_free(&node_dbuffer_bed);
  }

}
#else
void comm_adj_nodes(int *node_unref_flags,  /* flags nodes to be eliminated */
                    NODE_ENTRY ** node_hashtab  /* node hash table */
  )
{
}
#endif

/* move the departing node data */

#include "global_header.h"

#define NDIM_MSG 5  /* the number of items sent in the dimension message */

void comm_node_data(int *nnode_out, /* the number of nodes going to each processor */
                    int **nodes_out,    /* the nodes going out to each processor */
                    SNODE * old_node_pair,    /* node numbers before renumbering */
                    int old_nnode,   /* number of nodes before repartitioning */
                    SMODEL *mod
                    
                    )
{ 
    SGRID *g = mod->grid;
    int i;                        /* loop counters */
    int *nnode_in;                /* number of nodes coming into processor */
    int isd;                      /* loop counter over the processor */
    int ie;                       /* loop counter over the elements */
    int *pack_node;               /* flag for whether or not node is to be packed */
    int nnode = g->nnodes;
    int nelem3d = g->nelems3d;
    int nelem2d = g->nelems2d;
    int nelem1d = g->nelems1d;
    int npes = g->smpi->npes;
    
    MESSG_BUFFER recv_nd_msg;     /* the incoming message buffer */
    MESSG_BUFFER *send_nd_msg;    /* the outgoing message buffers */
    MESSG_BUFFER elem3d_ibuffer;   /* the buffer for the elem3d part of the message */
    MESSG_BUFFER elem2d_ibuffer;   /* the buffer for the elem2d part of the message */
    MESSG_BUFFER elem1d_ibuffer;   /* the buffer for the elem1d part of the message */
    MESSG_BUFFER elem3d_dbuffer;   /* the buffer for the elem3d part of the message */
    MESSG_BUFFER elem2d_dbuffer;   /* the buffer for the elem2d part of the message */
    MESSG_BUFFER elem1d_dbuffer;   /* the buffer for the elem1d part of the message */
    MESSG_BUFFER node_ibuffer;    /* the buffer for the integer node part of the message */
    MESSG_BUFFER node_dbuffer;    /* the buffer for the double node part of the message */
    MESSG_BUFFER node_dbuffer_surface;    /* the buffer for the surface double node part of the message */
    MESSG_BUFFER node_dbuffer_bed;    /* the buffer for the surface double node part of the message */
    NODE_LIST_ITEM *node_hashtab[HASHSIZE];   /* tops of the linked lists for the node hash table */
    SNODE gn;
    int nd_pntr;
    int loc;      /* the temporary pointer to an adjacent node */
    ELEM3D_LIST_ITEM *elem3d_hashtab[HASHSIZE];   /* tops of the linked lists for the 3d element hash table */
    ELEM2D_LIST_ITEM *elem2d_hashtab[HASHSIZE];   /* tops of the linked lists for the 2d element hash table */
    ELEM1D_LIST_ITEM *elem1d_hashtab[HASHSIZE];   /* tops of the linked lists for the 1d element hash table */
    
    /* initialize buffer sizes */
    elem3d_ibuffer.size = 0;
    elem2d_ibuffer.size = 0;
    elem1d_ibuffer.size = 0;
    elem3d_dbuffer.size = 0;
    elem2d_dbuffer.size = 0;
    elem1d_dbuffer.size = 0;
    node_ibuffer.size = 0;
    node_dbuffer.size = 0;
    node_dbuffer_surface.size = 0;
    node_dbuffer_bed.size = 0;
    
    /* construct the hash tables */
    for (i = 0; i < HASHSIZE; i++) {
        node_hashtab[i] = NULL;
        elem3d_hashtab[i] = NULL;
        elem2d_hashtab[i] = NULL;
        elem1d_hashtab[i] = NULL;
    }
    
    for (i = 0; i < nnode; i++) {
        node_hash_add_entry(g->node[i], node_hashtab, i, g->smpi->npes);
    }
    for (i = 0; i < old_nnode; i++) {
        nd_pntr = node_hash_lookup(old_node_pair[i], node_hashtab, g->smpi->npes);
        if (nd_pntr < 0)
            
            node_hash_add_entry(old_node_pair[i], node_hashtab, i, g->smpi->npes);
    }
    for (ie = 0; ie < nelem3d; ie++)
        elem3d_hash_add_entry(g->elem3d[ie].nodes[0], g->elem3d[ie].nodes[1], g->elem3d[ie].nodes[2], g->elem3d[ie].nodes[3], elem3d_hashtab, ie);
    for (ie = 0; ie < nelem2d; ie++)
        elem2d_hash_add_entry(g->elem2d[ie].nodes[0], g->elem2d[ie].nodes[1], g->elem2d[ie].nodes[2], elem2d_hashtab, ie);
    for (ie = 0; ie < nelem1d; ie++)
        elem1d_hash_add_entry(g->elem1d[ie].nodes[0], g->elem1d[ie].nodes[1], elem1d_hashtab, ie);
    
    /* correct adjacency information for all nodes using the new node numbers */
    for(i = 0; i < g->my_nnodes; i++){
        if(g->node[i].parent_res_pe[0] == UNSET_INT && g->node[i].parent_res_pe[1] != UNSET_INT)
            tl_error("One adjacent node found in comm_node_data.");
        else if(g->node[i].parent_res_pe[0] != UNSET_INT && g->node[i].parent_res_pe[1] == UNSET_INT)
            tl_error("One adjacent node found in comm_node_data.");
        else if(g->node[i].parent_res_pe[0] != UNSET_INT && g->node[i].parent_res_pe[1] != UNSET_INT){
            /* fix left adjacent node */
            gn.resident_id = g->node[i].parent_res_id[0];
            gn.resident_pe = g->node[i].parent_res_pe[0];
            nd_pntr = node_hash_lookup(gn, node_hashtab, g->smpi->npes);
            if(nd_pntr < 0)
                tl_error("Adjacent node lookup failed in comm_node_data.");
            loc = nd_pntr;
            g->node[i].parent_res_pe[0] = g->node[loc].resident_pe;
            g->node[i].parent_res_id[0] = g->node[loc].resident_id;
            
            /* fix right adjacent node */
            gn.resident_id = g->node[i].parent_res_id[1];
            gn.resident_pe = g->node[i].parent_res_pe[1];
            nd_pntr = node_hash_lookup(gn, node_hashtab, g->smpi->npes);
            if(nd_pntr < 0)
                tl_error("Adjacent node lookup failed in comm_node_data.");
            loc = nd_pntr;
            g->node[i].parent_res_pe[1] = g->node[loc].resident_pe;
            g->node[i].parent_res_id[1] = g->node[loc].resident_id;
        }
    }
    
    /* allocate space for the messages */
    send_nd_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
    for (isd = 0; isd < npes; isd++)
        messg_buffer_init(send_nd_msg + isd, isd);
    messg_buffer_init(&recv_nd_msg, UNSET_INT);
    
    /* allocate space for the pack node array */
    pack_node = (int *) tl_alloc(sizeof(int), nnode);
    
    /* loop over the subdomains, pack the appropriate data and send the message */
    nnode_in = (int *) tl_alloc(sizeof(int), npes);
    
    //tag(MPI_COMM_WORLD);
    
    for (isd = 0; isd < npes; isd++) {
        nnode_in[isd] = nnode_out[isd];
        if (nnode_out[isd] > 0) {
              //tag(MPI_COMM_WORLD);  MPI_Barrier(MPI_COMM_WORLD);
            
            /* initialize the auxilliary buffers */
            messg_buffer_free(&elem3d_ibuffer); //tag(MPI_COMM_WORLD);  //MPI_Barrier(MPI_COMM_WORLD);
            messg_buffer_free(&elem2d_ibuffer); //tag(MPI_COMM_WORLD);  //MPI_Barrier(MPI_COMM_WORLD);
            messg_buffer_free(&elem1d_ibuffer); //tag(MPI_COMM_WORLD);  //MPI_Barrier(MPI_COMM_WORLD);
            messg_buffer_free(&elem3d_dbuffer); //tag(MPI_COMM_WORLD);  //MPI_Barrier(MPI_COMM_WORLD);
            messg_buffer_free(&elem2d_dbuffer); //tag(MPI_COMM_WORLD);  //MPI_Barrier(MPI_COMM_WORLD);
            messg_buffer_free(&elem1d_dbuffer); //tag(MPI_COMM_WORLD);  //MPI_Barrier(MPI_COMM_WORLD);
            messg_buffer_free(&node_ibuffer); //tag(MPI_COMM_WORLD);  //MPI_Barrier(MPI_COMM_WORLD);
            messg_buffer_free(&node_dbuffer); //tag(MPI_COMM_WORLD);  //MPI_Barrier(MPI_COMM_WORLD);
            if(g->type == COLUMNAR){
                messg_buffer_free(&node_dbuffer_surface);
                messg_buffer_free(&node_dbuffer_bed);
            }
            //tag(MPI_COMM_WORLD);  MPI_Barrier(MPI_COMM_WORLD);
            
            messg_buffer_init(&elem3d_ibuffer, UNSET_INT);
            messg_buffer_init(&elem2d_ibuffer, UNSET_INT);
            messg_buffer_init(&elem1d_ibuffer, UNSET_INT);
            messg_buffer_init(&elem3d_dbuffer, UNSET_INT);
            messg_buffer_init(&elem2d_dbuffer, UNSET_INT);
            messg_buffer_init(&elem1d_dbuffer, UNSET_INT);
            messg_buffer_init(&node_ibuffer, UNSET_INT);
            messg_buffer_init(&node_dbuffer, UNSET_INT);
            if(g->type ==  COLUMNAR){
                messg_buffer_init(&node_dbuffer_surface, UNSET_INT);
                messg_buffer_init(&node_dbuffer_bed, UNSET_INT);
            }
            
            /* set the pack node array */
            for (i = 0; i < nnode; i++)
                pack_node[i] = NO;
            
            for (i = 0; i < nnode_out[isd]; i++)
                pack_node[nodes_out[isd][i]] = YES;
            /* load the element buffers */
            mesh_pack_elem3d(pack_node, &elem3d_ibuffer, &elem3d_dbuffer, g);
            mesh_pack_elem2d(pack_node, &elem2d_ibuffer, &elem2d_dbuffer, g);
            mesh_pack_elem1d(pack_node, &elem1d_ibuffer, &elem1d_dbuffer, g);
            /* load the node buffers */
            mesh_pack_node(nnode_out[isd], nodes_out[isd], &node_dbuffer, &node_dbuffer_surface, &node_dbuffer_bed, &node_ibuffer, mod);
            
            /* set the buffer type */
            send_nd_msg[isd].type = MESSG_PACKED;
            
            /* allocate space to pack the buffer */
            
            messg_pack_alloc(&node_ibuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack_alloc(&node_dbuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            if(g->type == COLUMNAR  && node_dbuffer_surface.nitem > 0){
                messg_pack_alloc(&node_dbuffer_surface, send_nd_msg + isd, g->smpi->ADH_COMM);
                messg_pack_alloc(&node_dbuffer_bed, send_nd_msg + isd, g->smpi->ADH_COMM);
            }
            messg_pack_alloc(&elem3d_ibuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack_alloc(&elem3d_dbuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack_alloc(&elem2d_ibuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack_alloc(&elem2d_dbuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack_alloc(&elem1d_ibuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack_alloc(&elem1d_dbuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            
            /* pack the buffer */
            messg_pack(&node_ibuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack(&node_dbuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            if(g->type == COLUMNAR && node_dbuffer_surface.nitem > 0){
                messg_pack(&node_dbuffer_surface, send_nd_msg + isd, g->smpi->ADH_COMM);
                messg_pack(&node_dbuffer_bed, send_nd_msg + isd, g->smpi->ADH_COMM);
            }
            messg_pack(&elem3d_ibuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack(&elem3d_dbuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack(&elem2d_ibuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack(&elem2d_dbuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack(&elem1d_ibuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            messg_pack(&elem1d_dbuffer, send_nd_msg + isd, g->smpi->ADH_COMM);
            
            /* send the message */
            messg_asend(send_nd_msg + isd, TAG_NODE_DATA, g->smpi);
            
        }
    }
    //tl_error("check\n");
    //tag(MPI_COMM_WORLD);
    
    /* free the pack node array */
    pack_node = (int *) tl_free(sizeof(int), nnode, pack_node);
    /* loop over the incoming messages */
    recv_nd_msg.type = MESSG_PACKED;
    messg_incoming(nnode_in, g->smpi);
    
    for (isd = 0; isd < npes; isd++)
    {
        
        if (nnode_in[isd] > 0)
        {
            
            /* initialize the auxilliary buffers */
            /* note:  since these buffers are reused, they must be freed before the buffer pointers
             are reset in messg_buffer_init */
            messg_buffer_free(&elem3d_ibuffer);
            messg_buffer_free(&elem2d_ibuffer);
            messg_buffer_free(&elem1d_ibuffer);
            messg_buffer_free(&elem3d_dbuffer);
            messg_buffer_free(&elem2d_dbuffer);
            messg_buffer_free(&elem1d_dbuffer);
            messg_buffer_free(&node_ibuffer);
            messg_buffer_free(&node_dbuffer);
            if(g->type == COLUMNAR){
                messg_buffer_free(&node_dbuffer_surface);
                messg_buffer_free(&node_dbuffer_bed);
            }
            messg_buffer_init(&elem3d_dbuffer, UNSET_INT);
            messg_buffer_init(&elem2d_dbuffer, UNSET_INT);
            messg_buffer_init(&elem1d_dbuffer, UNSET_INT);
            messg_buffer_init(&elem3d_ibuffer, UNSET_INT);
            messg_buffer_init(&elem2d_ibuffer, UNSET_INT);
            messg_buffer_init(&elem1d_ibuffer, UNSET_INT);
            messg_buffer_init(&node_ibuffer, UNSET_INT);
            messg_buffer_init(&node_dbuffer, UNSET_INT);
            if(g->type == COLUMNAR){
                messg_buffer_init(&node_dbuffer_surface, UNSET_INT);
                messg_buffer_init(&node_dbuffer_bed, UNSET_INT);
            }
            /* get a message */
            recv_nd_msg.sd = isd;
            messg_precv(&recv_nd_msg, TAG_NODE_DATA, g->smpi);
            /* unpacks the buffers from the recv buffer */
            messg_unpack(&recv_nd_msg, &node_ibuffer, g->smpi->ADH_COMM);
            messg_unpack(&recv_nd_msg, &node_dbuffer, g->smpi->ADH_COMM);
            if(g->type == COLUMNAR && node_dbuffer_surface.nitem > 0){
                messg_unpack(&recv_nd_msg, &node_dbuffer_surface, g->smpi->ADH_COMM);
                messg_unpack(&recv_nd_msg, &node_dbuffer_bed, g->smpi->ADH_COMM);
            }
            messg_unpack(&recv_nd_msg, &elem3d_ibuffer, g->smpi->ADH_COMM);
            messg_unpack(&recv_nd_msg, &elem3d_dbuffer, g->smpi->ADH_COMM);
            messg_unpack(&recv_nd_msg, &elem2d_ibuffer, g->smpi->ADH_COMM);
            messg_unpack(&recv_nd_msg, &elem2d_dbuffer, g->smpi->ADH_COMM);
            messg_unpack(&recv_nd_msg, &elem1d_ibuffer, g->smpi->ADH_COMM);
            messg_unpack(&recv_nd_msg, &elem1d_dbuffer, g->smpi->ADH_COMM);
            /* unpacks the node buffers */
            mesh_unpack_node(&node_dbuffer, &node_dbuffer_surface, &node_dbuffer_bed, &node_ibuffer, mod, node_hashtab, 0);
            /* unpacks the element buffers */
            mesh_unpack_elem3d(node_hashtab, elem3d_hashtab, &elem3d_ibuffer, &elem3d_dbuffer, mod);
            mesh_unpack_elem2d(node_hashtab, elem2d_hashtab, &elem2d_ibuffer, &elem2d_dbuffer, mod);
            mesh_unpack_elem1d(node_hashtab, elem1d_hashtab, &elem1d_ibuffer, &elem1d_dbuffer, mod);
            
        }
    }
    
    
    /* free the local buffer memory */
    //tl_error("check\n");
    messg_buffer_free(&recv_nd_msg); //tl_error("check\n");
    messg_buffer_free(&elem3d_ibuffer);
    messg_buffer_free(&elem2d_ibuffer);
    messg_buffer_free(&elem1d_ibuffer);
    messg_buffer_free(&elem3d_dbuffer);
    messg_buffer_free(&elem2d_dbuffer);
    messg_buffer_free(&elem1d_dbuffer);
    messg_buffer_free(&node_ibuffer);
    messg_buffer_free(&node_dbuffer);
    if(g->type == COLUMNAR){
        messg_buffer_free(&node_dbuffer_surface);
        messg_buffer_free(&node_dbuffer_bed);
    }
    
    /* wait for the asynchronous communication to complete */
    messg_wait(g->smpi);
    /* free the send buffers */
    for (isd = 0; isd < npes; isd++) messg_buffer_free(send_nd_msg + isd);
    nnode_in = (int *) tl_free(sizeof(int), npes, nnode_in);
    send_nd_msg = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), npes, send_nd_msg);
    
    /* free the hash tables */
    tl_list_free_all(NODE_LIST);
    tl_list_free_all(ELEM3D_LIST);
    tl_list_free_all(ELEM2D_LIST);
    tl_list_free_all(ELEM1D_LIST);
    
    
}

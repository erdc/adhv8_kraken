#include "global_header.h"

/* code to determine new node numbers for departing nodes */
void comm_new_node_nums(int *nnode_out, /* the number of nodes going to each processor */
                        int **nodes_out, /* the nodes going out to each processor */
                        int npes,
                        int myid,
                        SMODEL *mod
                        )
{
    SGRID *g = mod->grid;
    int i, j;                        /* loop counter */
    int isd;                      /* loop counter over the subdomains */
    int inode;                    /* the new node number */
    int *ibuffer;                 /* integer cast of the buffer */
    int *nnodes_in;               /* nodes coming into processor */
    MESSG_BUFFER *send_old_node_number;   /* send buffer for the old node numbers */
    MESSG_BUFFER recv_old_node_number;    /* recv buffer for the old node numbers */
    MESSG_BUFFER *send_new_node_number;   /* send buffer for the new node numbers */
    MESSG_BUFFER *recv_new_node_number;   /* recv buffer for the new node numbers */
    NODE_LIST_ITEM *node_hashtab[HASHSIZE];   /* tops of the linked lists for the node hash table */
    SNODE tmp_node;         /* temporary global node number used for lookup in hash table */
    
    /* allocations and initializations */
    nnodes_in = (int *) tl_alloc(sizeof(int), npes);
    send_old_node_number = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
    send_new_node_number = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
    recv_new_node_number = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
    
    for (isd = 0; isd < npes; isd++)
    {
        messg_buffer_init(send_old_node_number + isd, isd);
        messg_buffer_init(send_new_node_number + isd, isd);
        messg_buffer_init(recv_new_node_number + isd, isd);
        messg_buffer_alloc((nnode_out[isd]*3), sizeof(int), send_old_node_number + isd);
        messg_buffer_alloc(nnode_out[isd], sizeof(int), recv_new_node_number + isd);
    }
    
    messg_buffer_init(&recv_old_node_number, UNSET_INT);
    
    /* SEND THE OUT GOING NODE NUMBERS */
    for (isd = 0; isd < npes; isd++)
    {
        nnodes_in[isd] = nnode_out[isd];
        if (nnode_out[isd] > 0)
        {
            for (i = 0, ibuffer = (int *) send_old_node_number[isd].buffer; i < nnode_out[isd]; i++){
                j = (i*3);
                ibuffer[j] = nodes_out[isd][i];
                ibuffer[j+1] = g->node[nodes_out[isd][i]].global_surf_id;
                ibuffer[j+2] = g->node[nodes_out[isd][i]].global_bed_id;
            }
            send_old_node_number[isd].type = MESSG_INT;
            messg_asend(send_old_node_number + isd, TAG_NODE_OUT, g->smpi);             
        } 
    }
    
    //messg_barrier(MPI_COMM_WORLD);  // cjt :: not sure if this really needs to be here ... not in older code
    
    /* construct the hash table to convert ghost node numbers
     from global node numbers to local node numbers */
    for (i = 0; i < HASHSIZE; i++)
        node_hashtab[i] = NULL;
    
    for (i = g->my_nnodes; i < g->nnodes; i++)
        node_hash_add_entry(g->node[i], node_hashtab, i, npes);
    
    
    /* RECEIVE THE IN COMING NODE NUMBERS AND
     SEND THE NEW NUMBERS BACK */
    messg_incoming(nnodes_in, g->smpi);
    recv_old_node_number.type = MESSG_INT;
    
    for (isd = 0; isd < npes; isd++)
    {
        if (nnodes_in[isd] > 0)
        {
            /* receive the next message */
            recv_old_node_number.sd = isd;
            
            messg_precv(&recv_old_node_number, TAG_NODE_OUT, g->smpi);
            
            /* allocate memory for the reply */
            j=recv_old_node_number.nitem/3;
            messg_buffer_alloc(j, sizeof(int), send_new_node_number + recv_old_node_number.sd);
            
            /* parses the message */
            for (i = 0; i < j; i++)
            {
                /* load each incoming node in the temporary global node number */
                tmp_node.resident_id = *(((int *) recv_old_node_number.buffer) + (i*3));
                
                tmp_node.global_surf_id = *(((int *) recv_old_node_number.buffer) + ((i*3+1)));
                
                tmp_node.global_bed_id = *(((int *) recv_old_node_number.buffer) + ((i*3+2)));
                
                tmp_node.resident_pe = recv_old_node_number.sd;
                
                /* get a local node number for each incoming node */
                inode = node_get_local(tmp_node, node_hashtab, mod);
                g = mod->grid;
                
                /* fix the nodes global number on myid */
                g->node[inode].resident_id = inode;
                g->node[inode].resident_pe = myid;
                
                /* load the local node number in the nodes_in_num buffer */
                *(((int *) send_new_node_number[recv_old_node_number.sd].buffer) + i) = inode;
            }
            
            /* send node numbers back to sending PE */
            send_new_node_number[recv_old_node_number.sd].type = MESSG_INT;
            
            messg_asend(send_new_node_number + recv_old_node_number.sd, TAG_NODE_NUM, g->smpi);
            
        }
    }
    
    //messg_barrier(MPI_COMM_WORLD);   // cjt :: not sure if this really needs to be here ... not in older code
    
    /* RECEIVE THE NEW NUMBERS AND CHANGE THE GLOBAL NODE NUMBERS */
    for (isd = 0; isd < npes; isd++)
        if (nnode_out[isd] > 0)
        {
            recv_new_node_number[isd].type = MESSG_INT;
            messg_arecv(recv_new_node_number + isd, TAG_NODE_NUM, g->smpi);
        }
    
    /* wait for the asynchronous communication to complete */
    messg_wait(g->smpi);
    //messg_barrier(MPI_COMM_WORLD); // cjt :: not sure if this really needs to be here ... not in older code
    
    /* update nodes with new node numbers */
    for (isd = 0; isd < npes; isd++)
        for (i = 0, ibuffer = (int *) recv_new_node_number[isd].buffer; i < nnode_out[isd]; i++)
        {
            g->node[nodes_out[isd][i]].resident_id = ibuffer[i];
            g->node[nodes_out[isd][i]].resident_pe = isd;
        }
    
    //messg_barrier(MPI_COMM_WORLD);   // cjt :: not sure if this really needs to be here ... not in older code
    
    /* free the memory */
    for (isd = 0; isd < npes; isd++)
    {
        messg_buffer_free(send_old_node_number + isd);
        messg_buffer_free(send_new_node_number + isd);
        messg_buffer_free(recv_new_node_number + isd);
    }
    messg_buffer_free(&recv_old_node_number);
    nnodes_in = (int *) tl_free(sizeof(int), npes, nnodes_in);
    send_old_node_number = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), npes, send_old_node_number);
    send_new_node_number = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), npes, send_new_node_number);
    recv_new_node_number = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), npes, recv_new_node_number);
    tl_list_free_all(NODE_LIST);
    
}

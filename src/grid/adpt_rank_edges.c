/* ranks the edges */
#include "global_header.h"

#define TAG_RANK 773

/* function definition */
static int edge_comp(
                     const void *,
                     const void *
                     );

void adpt_rank_edges(
    EDGE_LIST_ITEM ** edge_hashtab,	/* the hash table of edges */
    SGRID *g                     /* grid nodes */
#ifdef _MESSG
    ,int *nrecv_edge, int *nsend_edge
#endif
) {

    int i, j;			/* loop counters */
    int nedge;			/* the number of edges */
    int *rank = NULL;		/* the ranks of the edges in hash table order */
    int *sorted_rank = NULL;	/* the ranks of the edges in sorted order */
    EDGE_RANK *edge_ranks = NULL;	/* the double lengths of the edges */
    EDGE_LIST_ITEM *next_edge;	/* pointer to the next edge */
    EDGE_LIST_ITEM *opp_order_edge;	/* same edges with opposite edge ordering */
    SNODE *node=g->node;                     /* grid nodes */
#ifdef _MESSG
    int isd;			/* loop counter over the subdomains */
    int *bpntr;			/* pointer to the buffer location */
    int myid = g->smpi->myid;
#endif
    
    /* initialize edge ranks in the edge hash table */
    for(i = 0; i < HASHSIZE; i++) {
        next_edge = edge_hashtab[i];
        while(next_edge != NULL) {
            next_edge->rank = UNSET_INT;
            next_edge = next_edge->next;
        }
    }
    
    /* count the # of unique edges and allocates space for the sort */
    for(i = 0, nedge = 0; i < HASHSIZE; i++) {
        next_edge = edge_hashtab[i];
        while(next_edge != NULL) {
#ifdef _MESSG
            if(node[next_edge->nd2].resident_pe== myid &&
               (node[next_edge->nd1].resident_pe > node[next_edge->nd2].resident_pe ||
                (node[next_edge->nd1].resident_pe == node[next_edge->nd2].resident_pe &&
                 node[next_edge->nd1].resident_id > node[next_edge->nd2].resident_id))) {
                    nedge++;
                }
            
#else
            if(next_edge->nd1 > next_edge->nd2) {
                nedge++;
            }
#endif
            next_edge = next_edge->next;
        }
    }
    
    /* allocate local arrays */
    if(nedge > 0) {
        edge_ranks = (EDGE_RANK *) tl_alloc(sizeof(EDGE_RANK), nedge);
        rank = (int *)tl_alloc(sizeof(int), nedge);
        sorted_rank = (int *)tl_alloc(sizeof(int), nedge);
    }
    
    /* form the edge lengths */
    for(i = 0, j = 0; i < HASHSIZE; i++) {
        next_edge = edge_hashtab[i];
        while(next_edge != NULL) {
#ifdef _MESSG
            if(node[next_edge->nd2].resident_pe == myid &&
               (node[next_edge->nd1].resident_pe > node[next_edge->nd2].resident_pe ||
                (node[next_edge->nd1].resident_pe == node[next_edge->nd2].resident_pe &&
                 node[next_edge->nd1].resident_id > node[next_edge->nd2].resident_id))) {
                    edge_ranks[j].length = DIST_3D(node[next_edge->nd1], node[next_edge->nd2]);
                    j++;
                }
#else
            if(next_edge->nd1 > next_edge->nd2) {
                edge_ranks[j].length = DIST_3D(node[next_edge->nd1], node[next_edge->nd2]);
                j++;
            }
#endif
            next_edge = next_edge->next;
        }
    }
    
    /* slide the indices into the edges */
    for(i = 0; i < nedge; i++) {
        edge_ranks[i].number = i;
    }
    
    /* sort the edges locally */
    qsort(edge_ranks, nedge, sizeof(EDGE_RANK), edge_comp);
    
    /* rank the edges globally */
    for(i = 0; i < nedge; i++) {
        sorted_rank[i] = -1;
    }
    
#ifdef _MESSG
    findrank(edge_ranks, nedge, sorted_rank, g->smpi);
#else
    for(i = 0; i < nedge; i++) {
        sorted_rank[i] = i;
    }
#endif
    
    /* compute the ranks */
    for(i = 0; i < nedge; i++)
        rank[edge_ranks[i].number] = sorted_rank[i];
    
    /* assign ranks to the unique edges */
    for(i = 0, j = 0; i < HASHSIZE; i++) {
        next_edge = edge_hashtab[i];
        while(next_edge != NULL) {
#ifdef _MESSG
            if(node[next_edge->nd2].resident_pe == myid &&
              (node[next_edge->nd1].resident_pe > node[next_edge->nd2].resident_pe ||
                (node[next_edge->nd1].resident_pe == node[next_edge->nd2].resident_pe &&
                 node[next_edge->nd1].resident_id > node[next_edge->nd2].resident_id))) {
                    next_edge->rank = rank[j];
                   j++;
                }
#else
            if(next_edge->nd1 > next_edge->nd2) {
                next_edge->rank = rank[j];
                j++;
            }
#endif
            next_edge = next_edge->next;
        }
    }
    
    /* update the edge ranks */
#ifdef _MESSG
    /* barrier - I think this can be removed */
//    messg_barrier();
    
    /* checks that the buffers are large enough */
    for(isd = 0; isd < g->smpi->npes; isd++) {
        messg_buffer_alloc(nsend_edge[isd], sizeof(int), g->smpi->send_edge_msg + isd);
        messg_buffer_alloc(nrecv_edge[isd], sizeof(int), g->smpi->recv_edge_msg + isd);
    }
    
    /* post the receives */
    for(isd = 0; isd < g->smpi->npes; isd++) {
        if(nrecv_edge[isd] != 0) {
            g->smpi->recv_edge_msg[isd].type = MESSG_INT;
            messg_arecv(g->smpi->recv_edge_msg + isd, TAG_RANK, g->smpi);
        }
    }
    
    /* load and send the send buffers */
    for(isd = 0; isd < g->smpi->npes; isd++) {
        if(nsend_edge[isd] != 0) {
            bpntr = (int *)g->smpi->send_edge_msg[isd].buffer;
            for(i = 0; i < nsend_edge[isd]; i++) {
                if(g->smpi->send_edge_key[isd][i] == NULL)
                    tl_error("Bad key found in adpt_rank_edges.\n");
                bpntr[i] = g->smpi->send_edge_key[isd][i]->rank;
            }
            g->smpi->send_edge_msg[isd].type = MESSG_INT;
            messg_asend(g->smpi->send_edge_msg + isd, TAG_RANK, g->smpi);
        }
    }
    
    /* wait for the asynchronous communications to clear */
    messg_wait(g->smpi);
    
    /* use the recv buffers to set the edge ranks */
    for(isd = 0; isd < g->smpi->npes; isd++)
        if(nrecv_edge[isd] != 0) {
            bpntr = (int *)g->smpi->recv_edge_msg[isd].buffer;
            for(i = 0; i < nrecv_edge[isd]; i++)
                if(g->smpi->recv_edge_key[isd][i] != NULL)
                    g->smpi->recv_edge_key[isd][i]->rank = bpntr[i];
        }
#endif
   int error_code=0; 
    /* assign ranks to the remaining edges */
    /* note :: two edges with same, but reverse nodes will have the same rank */
    for(i = 0; i < HASHSIZE; i++) {
        next_edge = edge_hashtab[i];
        while(next_edge != NULL) {
            if(next_edge->rank == UNSET_INT) {
                opp_order_edge = edge_hash_lookup(next_edge->nd2, next_edge->nd1, edge_hashtab);
                if(opp_order_edge == NULL) {
                    tl_error("Edge not found in adpt_rank_edges.\n");
                }
                if(opp_order_edge->rank == UNSET_INT) {
                
                  tl_error("Rank not set for either ordering of edge in adpt_rank_edges.");
                  
                }
                next_edge->rank = opp_order_edge->rank;
            }
            next_edge = next_edge->next;
        }
    }
    
    /* frees the edge rank memory */
    edge_ranks = (EDGE_RANK *) tl_free(sizeof(EDGE_RANK), nedge, edge_ranks);
    rank = (int *)tl_free(sizeof(int), nedge, rank);
    sorted_rank = (int *)tl_free(sizeof(int), nedge, sorted_rank);
}

/* an edge comparison function */
int edge_comp(
              const void *pntr1,		/* 1st pointer */
              const void *pntr2		/* 2nd pointer */
)
{
    EDGE_RANK *edge1 = (EDGE_RANK *) pntr1;	/* the first edge */
    EDGE_RANK *edge2 = (EDGE_RANK *) pntr2;	/* the second edge */
    if(edge1->length > edge2->length)
        return (1);
    else
        return (-1);
}

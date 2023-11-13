/* fixes the adjacent nodes */

#include "global_header.h"

void adpt_fix_adj(
                  SGRID *grid,              // the grid
                  int *node_unref_flags     // flags nodes to be eliminated */
#ifdef _MESSG
                  ,NODE_LIST_ITEM ** node_hashtab	// node hash table
#endif
) {
    int i;			/* loop counter */
    int nd1;			/* the 1st adjacent node */
    int nd2;			/* the 2nd adjacent node */
    int npes = grid->smpi->npes;
    SNODE nd_pntr;
    
    
    /* fix adjacent nodes */
    for(i = 0; i < grid->nnodes; i++) {
        if(node_unref_flags[i] == YES) {
            /* look up adjacent nodes */
#ifdef _MESSG
            nd_pntr.resident_pe = grid->node[i].parent_res_pe[0];
            nd_pntr.resident_id = grid->node[i].parent_res_id[0];
            nd1 = node_hash_lookup(nd_pntr, node_hashtab, npes);
            if(nd1 < 0) {
                tl_error("Left adjacent node lookup failed in adpt_unref.");
            }
            nd_pntr.resident_pe = grid->node[i].parent_res_pe[1];
            nd_pntr.resident_id = grid->node[i].parent_res_id[1];
            nd2 = node_hash_lookup(nd_pntr, node_hashtab, npes);
            if(nd2 < 0) {
                tl_error("Right adjacent node lookup failed in adpt_unref.");
            }
#else
            nd1 = grid->node[i].parent[0];
            nd2 = grid->node[i].parent[1];
#endif
            
            /* fix the adjacencies for nd1 */
#ifdef _MESSG
            if(grid->node[nd1].parent_res_pe[0] == grid->node[i].resident_pe &&
               grid->node[nd1].parent_res_id[0] == grid->node[i].resident_id) {
                grid->node[nd1].parent[0] = nd2;
                grid->node[nd1].parent_res_pe[0] = grid->node[nd2].resident_pe;
                grid->node[nd1].parent_res_id[0] = grid->node[nd2].resident_id;
            } else if(grid->node[nd1].parent_res_pe[1] == grid->node[i].resident_pe &&
                      grid->node[nd1].parent_res_id[1] == grid->node[i].resident_id) {
                grid->node[nd1].parent[1] = nd2;
                grid->node[nd1].parent_res_pe[1] = grid->node[nd2].resident_pe;
                grid->node[nd1].parent_res_id[1] = grid->node[nd2].resident_id;
            }
#else
            if(grid->node[nd1].parent[0] == i) {
                grid->node[nd1].parent[0] = nd2;
            } else if(grid->node[nd1].parent[1] == i) {
                grid->node[nd1].parent[1] = nd2;
            }
#endif
            
            /* fix the adjacencies for nd2 */
#ifdef _MESSG
            if(grid->node[nd2].parent_res_pe[0] == grid->node[i].resident_pe &&
               grid->node[nd2].parent_res_id[0] == grid->node[i].resident_id) {
                grid->node[nd2].parent[0] = nd1;
                grid->node[nd2].parent_res_pe[0] = grid->node[nd1].resident_pe;
                grid->node[nd2].parent_res_id[0] = grid->node[nd1].resident_id;
            } else if(grid->node[nd2].parent_res_pe[1] == grid->node[i].resident_pe &&
                      grid->node[nd2].parent_res_id[1] == grid->node[i].resident_id) {
                grid->node[nd2].parent[1] = nd1;
                grid->node[nd2].parent_res_pe[1] = grid->node[nd1].resident_pe;
                grid->node[nd2].parent_res_id[1] = grid->node[nd1].resident_id;
            }
#else
            if(grid->node[nd2].parent[0] == i) {
                grid->node[nd2].parent[0] = nd1;
            } else if(grid->node[nd2].parent[1] == i) {
                grid->node[nd2].parent[1] = nd1;
            }
#endif
            
        }
    }
    
    /* remove unrefined nodes */
    //printf("adpt_fix_adj1 :: node: %d \t string: %d \n",grid->node[291].id, grid->node[291].string);
    for(i = 0; i < grid->nnodes; i++) {
        if(node_unref_flags[i] == YES) {
            snode_init(&(grid->node[i]));
            grid->node[i].id = i;
        }
    }
    //printf("adpt_fix_adj2 :: node: %d \t string: %d \n",grid->node[291].id, grid->node[291].string);
    
    
}

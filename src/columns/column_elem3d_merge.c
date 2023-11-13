/* ADH Version 2.0.0 6-04 */
/* this routine merges a 3d elem */

#include "global_header.h"

void column_elem3d_merge(SGRID *grid,
                         int ielem, /* the right element to be merged */
                         int *new_local, /* the node of the element that is being removed */
                         int *node_unref_flags, NODE_LIST_ITEM ** node_hashtab, /* the node hashtable */
                         ELEM3D_LIST_ITEM ** elem3d_hashtab /* the 3D element hash table */
)
{
    int i;
    int adj11, adj12, adj21, adj22;/* local numbers of the adjacent nodes */
    int local_check_node1, local_check_node2;
    int new_node1, new_node2;     /* the newest node in the element */
    int old_node1, old_node2;     /* the old node replacing the newest node */
    ELEM3D_LIST_ITEM *elem_pntr;  /* pointer to the element in the hash table */
#ifdef _MESSG
    SNODE nd_pntr1, nd_pntr2;      /* pointer to the node in the hash table */
#endif
    SNODE nodes[NDONPRISM];

    if (grid->elem3d[ielem].nnodes == NDONTET){
    /**************************** TETS *******************************/
        new_node1 = grid->elem3d[ielem].nodes[new_local[0]];
    
    /* get adjacent nodes */
#ifdef _MESSG
        nd_pntr1.resident_pe=grid->node[new_node1].parent_res_pe[0];
        nd_pntr1.resident_id=grid->node[new_node1].parent_res_id[0];
        adj11 = -1;
        adj11 = node_hash_lookup(nd_pntr1, node_hashtab, grid->smpi->npes);
        if (adj11 < 0) tl_error("Adjacent node parent0 lookup failed in elem2d_merge!");
        nd_pntr1.resident_pe=grid->node[new_node1].parent_res_pe[1];
        nd_pntr1.resident_id=grid->node[new_node1].parent_res_id[1];
        adj12 = -1;
        adj12 = node_hash_lookup(nd_pntr1, node_hashtab, grid->smpi->npes);
        if (adj12 < 0) tl_error("Adjacent node parent1 lookup failed in elem2d_merge!");
#else
        adj11 = grid->node[new_node1].parent[0];
        adj12 = grid->node[new_node1].parent[1];
#endif
        
        old_node1 = UNSET_INT;
        for (i=0; i<grid->elem3d[ielem].nnodes; i++){
            if (adj11 == grid->elem3d[ielem].nodes[i]){
                old_node1=adj12;
                local_check_node1 = i;
                break;
            }
        }
        if (old_node1 == UNSET_INT){
            for (i=0; i<grid->elem3d[ielem].nnodes; i++){
                if (adj12 == grid->elem3d[ielem].nodes[i]){
                    old_node1=adj11;
                    local_check_node1 = i;
                    break;
                }
            }
        }
        if (old_node1 == UNSET_INT)
            tl_error("Neither adjacent node found in elem3d_merge.");

		
        /* reset the nodes */
        grid->elem3d[ielem].nodes[new_local[0]] = old_node1;
        grid->elem3d[ielem].levels[new_local[0]] = UNSET_INT; /* Gajanan gkc adding */
        
        node_unref_flags[grid->elem3d[ielem].nodes[new_local[0]]] = NO;
        
        /* look up element & return if duplicate -
         NOTE:  if we find the element then this is the second copy of the element to
         be unrefined, but if it is not there then this is the first copy of the element.
         We are using the hash table to remove duplicate copies */
        elem_pntr = elem3d_hash_lookup(grid->elem3d[ielem].nodes[0], grid->elem3d[ielem].nodes[1], grid->elem3d[ielem].nodes[2], grid->elem3d[ielem].nodes[3], elem3d_hashtab);
        /* second copy, initialize the element */
        if (elem_pntr != NULL) {
            grid->elem3d[elem_pntr->ielem].levels[local_check_node1] = grid->elem3d[ielem].levels[local_check_node1];
            selem3d_init(&grid->elem3d[ielem]);
            grid->elem_error[ielem]=UNSET_FLT;
        }
        /* first copy, fix it and enter it in the hash table */
        else {
            /* enter element in hash table */
            elem3d_hash_add_entry(grid->elem3d[ielem].nodes[0], grid->elem3d[ielem].nodes[1], grid->elem3d[ielem].nodes[2], grid->elem3d[ielem].nodes[3], elem3d_hashtab, ielem);
            
            /* fix element */
            grid->elem_error[ielem] = 0.0;
            for (i=0; i<NDONTET; i++) {
                snode_copy(&(nodes[i]), grid->node[grid->elem3d[ielem].nodes[i]]);
            }
            get_tet_linear_djac_gradPhi(&(grid->elem3d[ielem]), nodes, NULL);
        }
    }
    else if (grid->elem3d[ielem].nnodes == NDONPRISM){
    /*************************** PRISMS ******************************/
        new_node1 = grid->elem3d[ielem].nodes[new_local[0]];
        new_node2 = grid->elem3d[ielem].nodes[new_local[1]];
    
    /* get adjacent nodes */
#ifdef _MESSG
        nd_pntr1.resident_pe=grid->node[new_node1].parent_res_pe[0];
        nd_pntr1.resident_id=grid->node[new_node1].parent_res_id[0];
        adj11 = -1;
        adj11 = node_hash_lookup(nd_pntr1, node_hashtab, grid->smpi->npes);
        if (adj11 < 0) tl_error("Adjacent node parent0 lookup failed in elem2d_merge!");
        nd_pntr1.resident_pe=grid->node[new_node1].parent_res_pe[1];
        nd_pntr1.resident_id=grid->node[new_node1].parent_res_id[1];
        adj12 = -1;
        adj12 = node_hash_lookup(nd_pntr1, node_hashtab, grid->smpi->npes);
        if (adj12 < 0) tl_error("Adjacent node parent1 lookup failed in elem2d_merge!");

        nd_pntr2.resident_pe=grid->node[new_node2].parent_res_pe[0];
        nd_pntr2.resident_id=grid->node[new_node2].parent_res_id[0];
        adj21 = -1;
        adj21 = node_hash_lookup(nd_pntr2, node_hashtab, grid->smpi->npes);
        if (adj21 < 0) tl_error("Adjacent node parent0 lookup failed in elem2d_merge!");
        nd_pntr2.resident_pe=grid->node[new_node2].parent_res_pe[1];
        nd_pntr2.resident_id=grid->node[new_node2].parent_res_id[1];
        adj22 = -1;
        adj22 = node_hash_lookup(nd_pntr2, node_hashtab, grid->smpi->npes);
        if (adj22 < 0) tl_error("Adjacent node parent1 lookup failed in elem2d_merge!");


#else
        adj11 = grid->node[new_node1].parent[0];
        adj12 = grid->node[new_node1].parent[1];
        adj21 = grid->node[new_node2].parent[0];
        adj22 = grid->node[new_node2].parent[1];
#endif
        
        old_node1 = UNSET_INT;
        for (i=0; i<grid->elem3d[ielem].nnodes; i++){
            if (adj11 == grid->elem3d[ielem].nodes[i]){
                old_node1=adj12;
                local_check_node1 = i;
                break;
            }
        }
        if (old_node1 == UNSET_INT){
            for (i=0; i<grid->elem3d[ielem].nnodes; i++){
                if (adj12 == grid->elem3d[ielem].nodes[i]){
                    old_node1=adj11;
                    local_check_node1 = i;
                    break;
                }
            }
        }

        old_node2 = UNSET_INT;
        for (i=0; i<grid->elem3d[ielem].nnodes; i++){
            if (adj21 == grid->elem3d[ielem].nodes[i]){
                old_node2=adj22;
                local_check_node2 = i;
                break;
            }
        }
        if (old_node2 == UNSET_INT){
            for (i=0; i<grid->elem3d[ielem].nnodes; i++){
                if (adj22 == grid->elem3d[ielem].nodes[i]){
                    old_node2=adj21;
                    local_check_node2 = i;
                    break;
                }
            }
        }

        if (old_node1 == UNSET_INT || old_node2 == UNSET_INT)
            tl_error("Neither adjacent node found in elem3d_merge.");

        /* reset the nodes */
        grid->elem3d[ielem].nodes[new_local[0]] = old_node1;
        grid->elem3d[ielem].nodes[new_local[1]] = old_node2;
        grid->elem3d[ielem].levels[new_local[0]] = UNSET_INT; /* Gajanan gkc adding */
        grid->elem3d[ielem].levels[new_local[1]] = UNSET_INT; /* Gajanan gkc adding */
        
        node_unref_flags[grid->elem3d[ielem].nodes[new_local[0]]] = NO;
        node_unref_flags[grid->elem3d[ielem].nodes[new_local[1]]] = NO;
        
        /* look up element & return if duplicate -
         NOTE:  if we find the element then this is the second copy of the element to
         be unrefined, but if it is not there then this is the first copy of the element.
         We are using the hash table to remove duplicate copies */
        elem_pntr = elem3d_hash_lookup(grid->elem3d[ielem].nodes[0], grid->elem3d[ielem].nodes[1], grid->elem3d[ielem].nodes[2], grid->elem3d[ielem].nodes[3], elem3d_hashtab);
        /* second copy, initialize the element */
        if (elem_pntr != NULL) {
            grid->elem3d[elem_pntr->ielem].levels[local_check_node1] = grid->elem3d[ielem].levels[local_check_node1];
            grid->elem3d[elem_pntr->ielem].levels[local_check_node2] = grid->elem3d[ielem].levels[local_check_node2];
            selem3d_init(&grid->elem3d[ielem]);
            grid->elem_error[ielem]=UNSET_FLT;
        }
        /* first copy, fix it and enter it in the hash table */
        else {
            /* enter element in hash table */
            elem3d_hash_add_entry(grid->elem3d[ielem].nodes[0], grid->elem3d[ielem].nodes[1], grid->elem3d[ielem].nodes[2], grid->elem3d[ielem].nodes[3], elem3d_hashtab, ielem);
            
            /* fix element */
            grid->elem_error[ielem] = 0.0;
            for (i=0; i<NDONPRISM; i++) {
                snode_copy(&(nodes[i]), grid->node[grid->elem3d[ielem].nodes[i]]);
            }
            /* Gajanan gkc - I don't know what to do here. */
            //get_triPrism_linear_djac_gradPhi(.....);
        }
    }
}


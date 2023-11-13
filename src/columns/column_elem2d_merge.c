/* ADH Version 2.0.0 6-04 */
/* this routine merges a 2d elem */

#include "global_header.h"

int column_elem2d_merge(
#ifdef _MESSG
                        SGRID *grid,
                        int ielem, /* the right element to be merged */
                        int *new_local, /* the node of the element that is being removed */
                        int *new_elem_num_pntr,    /* the pointer to the merged element number */
                        NODE_LIST_ITEM ** node_hashtab,    /* the node hashtable */
                        ELEM2D_LIST_ITEM ** elem2d_hashtab /* the 2D element hash table */
#else
                        SGRID *grid,
                        int ielem, /* the right element to be merged */
                        int *new_local, /* the node of the element that is being removed */
                        int *new_elem_num_pntr,    /* the pointer to the merged element number */
                        ELEM2D_LIST_ITEM ** elem2d_hashtab /* the 2D element hash table */
#endif
)
{
    int i;
    int adj11, adj12, adj21, adj22;/* local numbers of the adjacent nodes */
    int check_node1, check_node2, check_node;               /* node to check if it needs this element */
    int hi_level = 0;             /* the highest level in the element */
    int current_level;            /* the level of the current node */
    int new_node1, new_node2;     /* the newest node in the element */
    int old_node1, old_node2;     /* the old node replacing the newest node */
    int local_check_node1, local_check_node2; /* the local number of the check node */
    ELEM2D_LIST_ITEM *elem_pntr;  /* pointer to the element in the hash table */
#ifdef _MESSG
    SNODE nd_pntr1, nd_pntr2;      /* pointer to the node in the hash table */
#endif
    
    if (grid->elem2d[ielem].nnodes == NDONTRI){
        new_node1 = grid->elem2d[ielem].nodes[new_local[0]];
    
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
        if (adj12 < 0) tl_error("Adjacent node parent0 lookup failed in elem2d_merge!");
     
#else
        adj11 = grid->node[new_node1].parent[0];
        adj12 = grid->node[new_node1].parent[1];
#endif
     
        if (adj11 == grid->elem2d[ielem].nodes[0]) {
            old_node1 = adj12;
            check_node1 = adj11;
            local_check_node1 = 0;
            check_node2 = grid->elem2d[ielem].nodes[3-(new_local[0]+local_check_node1)];
        }
        else if (adj11 == grid->elem2d[ielem].nodes[1]) {
            old_node1 = adj12;
            check_node1 = adj11;
            local_check_node1 = 1;
            check_node2 = grid->elem2d[ielem].nodes[3-(new_local[0]+local_check_node1)];
        }
        else if (adj11 == grid->elem2d[ielem].nodes[2]) {
            old_node1 = adj12;
            check_node1 = adj11;
            local_check_node1 = 2;
            check_node2 = grid->elem2d[ielem].nodes[3-(new_local[0]+local_check_node1)];
        }
        else if (adj12 == grid->elem2d[ielem].nodes[0]) {
            old_node1 = adj11;
            check_node1 = adj12;
            local_check_node1 = 0;
            check_node2 = grid->elem2d[ielem].nodes[3-(new_local[0]+local_check_node1)];
        }
        else if (adj12 == grid->elem2d[ielem].nodes[1]) {
            old_node1 = adj11;
            check_node1 = adj12;
            local_check_node1 = 1;
            check_node2 = grid->elem2d[ielem].nodes[3-(new_local[0]+local_check_node1)];
        }
        else if (adj12 == grid->elem2d[ielem].nodes[2]) {
            old_node1 = adj11;
            check_node1 = adj12;
            local_check_node1 = 2;
            check_node2 = grid->elem2d[ielem].nodes[3-(new_local[0]+local_check_node1)];
        }
        else {
            old_node1 = UNSET_INT;
            check_node = UNSET_INT;
            local_check_node1 = UNSET_INT;
            tl_error("Neither adjacent node found in elem2d_merge.");
        }
        
        /* reset the nodes and levels */
        grid->elem2d[ielem].nodes[new_local[0]] = old_node1;
        grid->elem2d[ielem].levels[new_local[0]] = UNSET_INT;
        
        /* look up element & return if duplicate -
         NOTE:  if we find the element then this is the second copy of the element to
         be unrefined, but if it is not there then this is the first copy of the element.
         We are using the hash table to remove duplicate copies */
        elem_pntr = elem2d_hash_lookup(grid->elem2d[ielem].nodes[0], grid->elem2d[ielem].nodes[1], grid->elem2d[ielem].nodes[2], elem2d_hashtab);
        /* second copy, initialize the element */
        if (elem_pntr != NULL) {
            grid->elem2d[elem_pntr->ielem].levels[local_check_node1] = grid->elem2d[ielem].levels[local_check_node1];
            *new_elem_num_pntr = elem_pntr->ielem;
            selem2d_init(&grid->elem2d[ielem]);
        }
        /* first copy, fix it and enter it in the hash table */
        else {
            /* enter element in hash table */
            elem2d_hash_add_entry(grid->elem2d[ielem].nodes[0], grid->elem2d[ielem].nodes[1], grid->elem2d[ielem].nodes[2], elem2d_hashtab, ielem);
            *new_elem_num_pntr = ielem;
            
            /* fix element */
            //grid->elem_error[ielem] = 0.0; // cjt :: took this out since we are only dealing with 3d element errors in 3d
            grid->elem2d[ielem].djac3d_fixed = 0.; // zero this so that lingrad recalculates
            
            get_triangle_linear_djac_nrml_gradPhi(&(grid->elem2d[ielem]), grid->node, NULL);
        }
    }
    else if (grid->elem2d[ielem].nnodes == NDONQUAD){
    
        new_node1 = grid->elem2d[ielem].nodes[new_local[0]];
        new_node2 = grid->elem2d[ielem].nodes[new_local[1]];
    
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
        if (adj12 < 0) tl_error("Adjacent node parent0 lookup failed in elem2d_merge!");

        nd_pntr2.resident_pe=grid->node[new_node2].parent_res_pe[0];
        nd_pntr2.resident_id=grid->node[new_node2].parent_res_id[0];
        adj21 = -1;
        adj21 = node_hash_lookup(nd_pntr2, node_hashtab, grid->smpi->npes);
        if (adj21 < 0) tl_error("Adjacent node parent0 lookup failed in elem2d_merge!");
        nd_pntr2.resident_pe=grid->node[new_node2].parent_res_pe[1];
        nd_pntr2.resident_id=grid->node[new_node2].parent_res_id[1];
        adj22 = -1;
        adj22 = node_hash_lookup(nd_pntr2, node_hashtab, grid->smpi->npes);
        if (adj22 < 0) tl_error("Adjacent node parent0 lookup failed in elem2d_merge!");

#else
        adj11 = grid->node[new_node1].parent[0];
        adj12 = grid->node[new_node1].parent[1];
        adj21 = grid->node[new_node2].parent[0];
        adj22 = grid->node[new_node2].parent[1];
#endif
    
        old_node1 = UNSET_INT;
        check_node1 = UNSET_INT;
        local_check_node1 = UNSET_INT;
        for (i=0; i<grid->elem2d[ielem].nnodes; i++){
            if (adj11 == grid->elem2d[ielem].nodes[i]) {
                old_node1 = adj12;
                check_node1 = adj11;
                local_check_node1 = i;
                break;
            }
        }
        if (old_node1 == UNSET_INT){
            for (i=0; i<grid->elem2d[ielem].nnodes; i++){
                if (adj12 == grid->elem2d[ielem].nodes[i]) {
                    old_node1 = adj11;
                    check_node1 = adj12;
                    local_check_node1 = i;
                    break;
                }
            }
        }

        old_node2 = UNSET_INT;
        check_node2 = UNSET_INT;
        local_check_node2 = UNSET_INT;
        for (i=0; i<grid->elem2d[ielem].nnodes; i++){
            if (adj21 == grid->elem2d[ielem].nodes[i]) {
                old_node2 = adj22;
                check_node2 = adj21;
                local_check_node2 = i;
                break;
            }
        }
        if (old_node2 == UNSET_INT){
            for (i=0; i<grid->elem2d[ielem].nnodes; i++){
                if (adj22 == grid->elem2d[ielem].nodes[i]) {
                    old_node2 = adj21;
                    check_node2 = adj22;
                    local_check_node2 = i;
                    break;
                }
            }
        }

        if (old_node1 == UNSET_INT || old_node2 == UNSET_INT)
            tl_error("Neither adjacent node found in elem2d_merge.");
        
        /* reset the nodes and levels */
        grid->elem2d[ielem].nodes[new_local[0]] = old_node1;
        grid->elem2d[ielem].levels[new_local[0]] = UNSET_INT;
        grid->elem2d[ielem].nodes[new_local[1]] = old_node2;
        grid->elem2d[ielem].levels[new_local[1]] = UNSET_INT;
        
        /* look up element & return if duplicate -
         NOTE:  if we find the element then this is the second copy of the element to
         be unrefined, but if it is not there then this is the first copy of the element.
         We are using the hash table to remove duplicate copies */
        elem_pntr = elem2d_hash_lookup(grid->elem2d[ielem].nodes[0], grid->elem2d[ielem].nodes[1], grid->elem2d[ielem].nodes[2], elem2d_hashtab);
        /* second copy, initialize the element */
        if (elem_pntr != NULL) {
            grid->elem2d[elem_pntr->ielem].levels[local_check_node1] = grid->elem2d[ielem].levels[local_check_node1];
            grid->elem2d[elem_pntr->ielem].levels[local_check_node2] = grid->elem2d[ielem].levels[local_check_node2];
            *new_elem_num_pntr = elem_pntr->ielem;
            selem2d_init(&grid->elem2d[ielem]);
        }
        /* first copy, fix it and enter it in the hash table */
        else {
            /* enter element in hash table */
            elem2d_hash_add_entry(grid->elem2d[ielem].nodes[0], grid->elem2d[ielem].nodes[1], grid->elem2d[ielem].nodes[2], elem2d_hashtab, ielem);
            *new_elem_num_pntr = ielem;
            
            /* fix element */
            //grid->elem_error[ielem] = 0.0; // cjt :: took this out since we are only dealing with 3d element errors in 3d
            grid->elem2d[ielem].djac3d_fixed = 0.; // zero this so that lingrad recalculates
            
            SVECT nds[grid->elem2d[ielem].nnodes];
            for (i=0; i<grid->elem2d[ielem].nnodes; i++) {
                nds[i].x = grid->node[ grid->elem2d[ielem].nodes[i] ].x;
                nds[i].y = grid->node[ grid->elem2d[ielem].nodes[i] ].y;
                nds[i].z = grid->node[ grid->elem2d[ielem].nodes[i] ].z;  // initial displacement should be added here
            }
            grid->elem2d[ielem].nrml = get_elem2d_normals(nds);
        }
    }
    
    /* indicate if elem is going out */
#ifdef _MESSG
    /* If this isn't a surface element we don't need to worry about updating the element levels */
    if(grid->node[new_node1].resident_pe == grid->smpi->myid && grid->node[check_node1].resident_pe != grid->smpi->myid ){
        return (grid->node[check_node1].resident_pe);
    }
    if (new_local[1] != UNSET_INT){
        if(grid->node[new_node2].resident_pe == grid->smpi->myid && grid->node[check_node2].resident_pe != grid->smpi->myid ){
            return (grid->node[check_node2].resident_pe);
        }
    }
    else {
        return (UNSET_INT);
    }
    tl_error("Code should never get here.");
#else
    return (UNSET_INT);
#endif
}

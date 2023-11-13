// CJT :: Merges two 3D element
// notes :: currently only for tetrahedrons although the algorithm is general
// notes :: will probably need some work for prisms, check column_elem3d_split
#include "global_header.h"

int elem3d_merge(SGRID *grid,                           /* the model grid */
                 int ielem,                            /* the element to be merged */
                 int *new_elem_num_pntr,               /* the pointer to the merged element number */
#ifdef _MESSG
                 NODE_LIST_ITEM ** node_hashtab,       /* the node hashtable */
#endif
                 ELEM3D_LIST_ITEM ** elem3d_hashtab,   /* the 3D element hash table */
                 int *node_unref_flags                 /* flags nodes to be eliminated */
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
    
    static int counter = 0;
    counter++;
    
    // find node to be removed
    int new_local = UNSET_INT;    /* the local # of the newest node in the element - if this is an original element it will be UNSET_INT */
    int hi_level = UNSET_INT;               // the highest level in the element
    int current_level = UNSET_INT;  // the level of the current node
    for (i=0; i<grid->elem3d[ielem].nnodes; i++){
        current_level = grid->elem3d[ielem].levels[i];
        if(current_level > hi_level) {
            hi_level = current_level;
            new_local = i;
        }
    }
    if(new_local == UNSET_INT) {tl_error("Tried to unref original element in elem3d_merge.");}
    
    
    if (grid->elem3d[ielem].nnodes == NDONTET){
        /**************************** TETS *******************************/
        new_node1 = grid->elem3d[ielem].nodes[new_local];
        
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
        if (old_node1 == UNSET_INT) {tl_error("Neither adjacent node found in elem3d_merge.");}
        
        /* reset the nodes */
        grid->elem3d[ielem].nodes[new_local] = old_node1;
        grid->elem3d[ielem].levels[new_local] = UNSET_INT; /* Gajanan gkc adding */
        node_unref_flags[grid->elem3d[ielem].nodes[new_local]] = NO;
        
        /* look up element & return if duplicate -
         NOTE:  if we find the element then this is the second copy of the element to
         be unrefined, but if it is not there then this is the first copy of the element.
         We are using the hash table to remove duplicate copies */
        elem_pntr = elem3d_hash_lookup(
                                       grid->elem3d[ielem].nodes[0],
                                       grid->elem3d[ielem].nodes[1],
                                       grid->elem3d[ielem].nodes[2],
                                       grid->elem3d[ielem].nodes[3], elem3d_hashtab);
        /* second copy, initialize the element */
        if (elem_pntr != NULL) {
            grid->elem3d[elem_pntr->ielem].levels[local_check_node1] = grid->elem3d[ielem].levels[local_check_node1];
            *new_elem_num_pntr = elem_pntr->ielem;
            selem3d_init(&grid->elem3d[ielem]);
            grid->elem_error[ielem]=0.0;
        }
        /* first copy, fix it and enter it in the hash table */
        else {
            /* enter element in hash table */
            elem3d_hash_add_entry(
                                  grid->elem3d[ielem].nodes[0],
                                  grid->elem3d[ielem].nodes[1],
                                  grid->elem3d[ielem].nodes[2],
                                  grid->elem3d[ielem].nodes[3], elem3d_hashtab, ielem);
            *new_elem_num_pntr = ielem;
            
            /* fix element */
            grid->elem_error[ielem] = 0.0;
            for (i=0; i<NDONTET; i++) {
                snode_copy(&(nodes[i]), grid->node[grid->elem3d[ielem].nodes[i]]);
            }
            get_tet_linear_djac_gradPhi(&(grid->elem3d[ielem]), nodes, NULL);
        }
    }
    else if (grid->elem3d[ielem].nnodes == NDONPRISM){
        tl_error("TRIANGULAR PRISMS NOT YET SUPPORTED IN ELEM3D_MERGE.C!");
    }
    
    /* indicate if elem is going out */
#ifdef _MESSG
    if(grid->node[new_node1].resident_pe == grid->smpi->myid && grid->node[local_check_node1].resident_pe != grid->smpi->myid){

    return (grid->node[local_check_node1].resident_pe);
  }
        return (UNSET_INT);
#else
    return (UNSET_INT);
#endif
}

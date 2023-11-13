#include "global_header.h"
/***********************************************************/

int elem2d_merge(
                 SGRID *grid,                            /* the model grid */
                 int ielem,                              /* the element to be merged */
                 int *new_elem_num_pntr,                 /* the pointer to the merged element number */
#ifdef _MESSG
                 NODE_LIST_ITEM **node_hashtab,          /* the node hashtable */
#endif
                 ELEM2D_LIST_ITEM **elem2d_hashtab,      /* the 2D element hash table */
                 int *node_unref_flags                   /* flags nodes to be eliminated */
)
{
    int i;                        /*loop counter*/
    int adj1, adj2;       /* local numbers of the adjacent nodes */
    int check_node1;      /* node to check if it needs this element */
    int check_node2;
    int current_level;        /* the level of the current node */
    int new_node;         /* the newest node in the element */
    int old_node;         /* the old node replacing the newest node */
    int local_check_node;     /* the local number of the check node */
    ELEM2D_LIST_ITEM *elem_pntr;  /* pointer to the element in the hash table */
#ifdef _MESSG
    SNODE nd_pntr;
#endif
    
    /* find node to be removed */
    int new_local = UNSET_INT;    /* the local # of the newest node in the element - if this is an original element it will be UNSET_INT */
    int hi_level = 0;     /* the highest level in the element */
    current_level = grid->elem2d[ielem].levels[0]; //printf("current_level: %d \t hi_level: %d\n",current_level,hi_level);
    if (current_level > hi_level) {
        hi_level = current_level;
        new_local = 0;
    }
    current_level = grid->elem2d[ielem].levels[1]; //printf("current_level: %d \t hi_level: %d\n",current_level,hi_level);
    if (current_level > hi_level) {
        hi_level = current_level;
        new_local = 1;
    }
    current_level = grid->elem2d[ielem].levels[2]; //printf("current_level: %d \t hi_level: %d\n",current_level,hi_level);
    if (current_level > hi_level) {
        hi_level = current_level;
        new_local = 2;
    }
    if(new_local == UNSET_INT) {
        printf("current_level: %d \t hi_level: %d \t new_local: %d \n",current_level, hi_level, new_local);
        printf("ielem: %di lvl1: %d \t lvl2: %d \t lvl3: %d \n",ielem,grid->elem2d[ielem].levels[0],grid->elem2d[ielem].levels[1],grid->elem2d[ielem].levels[2]);
        tl_error("Tried to unref original element in elem2d_merge.");
    }
    
    new_node = grid->elem2d[ielem].nodes[new_local];
    
    /* get adjacent nodes */
#ifdef _MESSG
    nd_pntr.resident_pe=grid->node[new_node].parent_res_pe[0];
    nd_pntr.resident_id=grid->node[new_node].parent_res_id[0];
    adj1 = -1;
    adj1 = node_hash_lookup(nd_pntr, node_hashtab, grid->smpi->npes);
    if (adj1 < 0) tl_error("Adjacent node parent0 lookup failed in elem2d_merge!");
    nd_pntr.resident_pe=grid->node[new_node].parent_res_pe[1];
    nd_pntr.resident_id=grid->node[new_node].parent_res_id[1];
    adj2 = -1;
    adj2 = node_hash_lookup(nd_pntr, node_hashtab, grid->smpi->npes);
    if (adj2 < 0) tl_error("Adjacent node parent0 lookup failed in elem2d_merge!");
#else
    adj1 = grid->node[new_node].parent[0];
    adj2 = grid->node[new_node].parent[1];
#endif
    
    if(adj1 == grid->elem2d[ielem].nodes[0]) {
        old_node = adj2;
        check_node1 = adj1;
        local_check_node = 0;
        check_node2 = grid->elem2d[ielem].nodes[3-(new_local+local_check_node)];
    } else if(adj1 == grid->elem2d[ielem].nodes[1]) {
        old_node = adj2;
        check_node1 = adj1;
        local_check_node = 1;
        check_node2 = grid->elem2d[ielem].nodes[3-(new_local+local_check_node)];
    } else if(adj1 == grid->elem2d[ielem].nodes[2]) {
        old_node = adj2;
        check_node1 = adj1;
        local_check_node = 2;
        check_node2 = grid->elem2d[ielem].nodes[3-(new_local+local_check_node)];
    } else if(adj2 == grid->elem2d[ielem].nodes[0]) {
        old_node = adj1;
        check_node1 = adj2;
        local_check_node = 0;
        check_node2 = grid->elem2d[ielem].nodes[3-(new_local+local_check_node)];
    } else if(adj2 == grid->elem2d[ielem].nodes[1]) {
        old_node = adj1;
        check_node1 = adj2;
        local_check_node = 1;
        check_node2 = grid->elem2d[ielem].nodes[3-(new_local+local_check_node)];
    } else if(adj2 == grid->elem2d[ielem].nodes[2]) {
        old_node = adj1;
        check_node1 = adj2;
        local_check_node = 2;
        check_node2 = grid->elem2d[ielem].nodes[3-(new_local+local_check_node)];
    } else {
        old_node = UNSET_INT;
        check_node1 = UNSET_INT;
        check_node2 = UNSET_INT;
        local_check_node = UNSET_INT;
        printf("\n ielem: %d adj1: %d adj2: %d element nodes: %d %d %d\n",ielem,adj1,adj2,grid->elem2d[ielem].nodes[0],grid->elem2d[ielem].nodes[1],grid->elem2d[ielem].nodes[2]);
        tl_error("Neither adjacent node found in elem2d_merge.");
    }
    
    /* reset the nodes and levels */
    grid->elem2d[ielem].nodes[new_local] = old_node;
    grid->elem2d[ielem].levels[new_local] = UNSET_INT;
    node_unref_flags[grid->elem2d[ielem].nodes[new_local]] = NO;
    
    /* look up element & return if duplicate -
     NOTE:  if we find the element then this is the second copy of the element to
     be unrefined, but if it is not there then this is the first copy of the element.
     We are using the hash table to remove duplicate copies */
    elem_pntr = elem2d_hash_lookup(grid->elem2d[ielem].nodes[0], grid->elem2d[ielem].nodes[1], grid->elem2d[ielem].nodes[2], elem2d_hashtab);
    /* second copy, initialize the element */
    if(elem_pntr != NULL) {
        grid->elem2d[elem_pntr->ielem].levels[local_check_node] = grid->elem2d[ielem].levels[local_check_node];
        *new_elem_num_pntr = elem_pntr->ielem;
        selem2d_init(&grid->elem2d[ielem]);
    }
    /* first copy, fix it and enter it in the hash table */
    else {
        /* enter element in hash table */
        elem2d_hash_add_entry(grid->elem2d[ielem].nodes[0], grid->elem2d[ielem].nodes[1], grid->elem2d[ielem].nodes[2], elem2d_hashtab, ielem);
        *new_elem_num_pntr = ielem;
        
        /* fix element */
        grid->elem_error[ielem] = 0.0;
        grid->elem2d[ielem].djac3d_fixed = 0.; // zero this so that lingrad recalculates
        
        // get 2D element djacs, normals and basis function gradients (elemental constants)
        if (grid->elem2d[ielem].nnodes == NDONTRI) {
            get_triangle_linear_djac_nrml_gradPhi(&(grid->elem2d[ielem]), grid->node, NULL);
        } else {
            SVECT nds[grid->elem2d[ielem].nnodes];
            for (i=0; i<grid->elem2d[ielem].nnodes; i++) {
                nds[i].x = grid->node[ grid->elem2d[ielem].nodes[i] ].x;
                nds[i].y = grid->node[ grid->elem2d[ielem].nodes[i] ].y;
                nds[i].z = grid->node[ grid->elem2d[ielem].nodes[i] ].z;  // initial displacement should be added here
            }
            grid->elem2d[ielem].nrml = get_elem2d_normals(nds);
        }
    }
    
#ifdef UNREF_CONSV
    *unref_node = new_node;
#endif
    
    /* indicate if elem is going out */
#ifdef _MESSG
    if(grid->node[new_node].resident_pe == grid->smpi->myid && grid->node[check_node1].resident_pe != grid->smpi->myid && grid->node[check_node1].resident_pe != grid->node[check_node2].resident_pe){
        
        return (grid->node[check_node1].resident_pe);
    }
    else
        return (UNSET_INT);
#else
    return (UNSET_INT);
#endif
}

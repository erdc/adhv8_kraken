#include "global_header.h"

/***********************************************************/
/* merge a 1d adapted element */
int elem1d_merge(
                 SGRID *grid,               /* the model grid */
                 int ielem,                 /* the element to be merged */
                 int *new_elem_num_pntr,    /* the pointer to the merged element number */
#ifdef _MESSG
                 NODE_LIST_ITEM ** node_hashtab,        /* the node hashtable */
#endif
                 ELEM1D_LIST_ITEM ** elem1d_hashtab,    /* the 1D element hash table */
                 int *node_unref_flags      /* flags nodes to be eliminated */

)
{
    int adj1, adj2;             /* local numbers of the adjacent nodes */
    int check_node;             /* node to check if it needs this element */
    int hi_level = 0;           /* the highest level in the element */
    int current_level;          /* the level of the current node */
    int new_local = UNSET_INT;  /* the local # of the newest node in the element - if this is an original element then this will maintain a value of UNSET_INT */
    int new_node;               /* the newest node in the element */
    int old_node;               /* the old node replacing the newest node */
    int local_check_node;       /* the local number of the check node */
    ELEM1D_LIST_ITEM *elem_pntr;/* pointer to the element in the hash table */
#ifdef _MESSG
    SNODE nd_pntr;    /* pointer to the node in the hash table */
#endif
    
    /* find node to be removed */
    current_level = grid->elem1d[ielem].levels[0];
    if(current_level > hi_level) {
        hi_level = current_level;
        new_local = 0;
    }
    current_level = grid->elem1d[ielem].levels[1];
    if(current_level > hi_level) {
        hi_level = current_level;
        new_local = 1;
    }
    if(new_local == UNSET_INT) {
        tl_error("Tried to unref original element in elem1d_merge.");
    }
    
    new_node = grid->elem1d[ielem].nodes[new_local];
    
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
  if (adj2 < 0) tl_error("Adjacent node parent1 lookup failed in elem2d_merge!");
#else
    adj1 = grid->node[new_node].parent[0];
    adj2 = grid->node[new_node].parent[1];
#endif
    
    if(adj1 == grid->elem1d[ielem].nodes[0])
    {
        old_node = adj2;
        check_node = adj1;
        local_check_node = 0;
    }
    else if(adj1 == grid->elem1d[ielem].nodes[1])
    {
        old_node = adj2;
        check_node = adj1;
        local_check_node = 1;
    }
    else if(adj2 == grid->elem1d[ielem].nodes[0])
    {
        old_node = adj1;
        check_node = adj2;
        local_check_node = 0;
    }
    else if(adj2 == grid->elem1d[ielem].nodes[1])
    {
        old_node = adj1;
        check_node = adj2;
        local_check_node = 1;
    }
    else
    {
        old_node = UNSET_INT;
        check_node = UNSET_INT;
        local_check_node = UNSET_INT;
        tl_error("Neither adjacent node found in elem1d_merge.");
    }
    
    /* reset the nodes and levels */
    grid->elem1d[ielem].nodes[new_local] = old_node;
    grid->elem1d[ielem].levels[new_local] = UNSET_INT;
    node_unref_flags[grid->elem1d[ielem].nodes[new_local]] = NO;
    
    /* look up element & return if duplicate -
     NOTE:  if we find the element then this is the second copy of the element to
     be unrefined, but if it is not there then this is the first copy of the element.
     We are using the hash table to remove duplicate copies */
    elem_pntr = elem1d_hash_lookup(grid->elem1d[ielem].nodes[0], grid->elem1d[ielem].nodes[1], elem1d_hashtab);
    /* second copy, initialize the element */
    if(elem_pntr != NULL) {
        grid->elem1d[elem_pntr->ielem].levels[local_check_node] = grid->elem1d[ielem].levels[local_check_node];
        *new_elem_num_pntr = elem_pntr->ielem;
        selem1d_init(&grid->elem1d[ielem]);
    }
    /* first copy, fix it and enter it in the hash table */
    else {
        /* enter element in hash table */
        elem1d_hash_add_entry(grid->elem1d[ielem].nodes[0], grid->elem1d[ielem].nodes[1], elem1d_hashtab, ielem);
        *new_elem_num_pntr = ielem;
        
        /* fix element */
        get_elem1d_linear_djac_gradPhi(grid, &(grid->elem1d[ielem]));
    }
    
    /* indicate if elem is going out */
#ifdef _MESSG
    if(grid->node[new_node].resident_pe == grid->smpi->myid && grid->node[check_node].resident_pe != grid->smpi->myid) 
      return (grid->node[check_node].resident_pe);
    else
        return (UNSET_INT);
#else
    return (UNSET_INT);
#endif
}

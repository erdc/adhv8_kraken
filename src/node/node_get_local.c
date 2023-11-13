/* this routine gets the local number of a node given the global node 
   number and the hashtable for the lookup */

#include "global_header.h"

#ifdef _MESSG
int node_get_local(SNODE gn,  /* the global node number */
                   NODE_LIST_ITEM ** node_hashtab,   /* hash table for the nodes */
                   SMODEL *mod
  )
{
  int nd_list_item; /* pointer to a node list item */
  int ind, ind_sur, ind_bed;                      /* the node */
  SGRID *g = mod->grid;
  int npes = g->smpi->npes;
   /* looks the node up in the hash table */
  nd_list_item = node_hash_lookup(gn, node_hashtab, npes);
  //printf("nd_list_item->local %d\n", nd_list_item->local);
  /* if the node was not found then create the node and 
     enter it in the hash table */
   
  if (nd_list_item < 0) {
   
    ind = node_new(mod);
    if((gn.global_surf_id != UNSET_INT) && (mod->grid->ndim==3)) {
      ind_sur = node_new_surface(mod);
      
      if (g->type == COLUMNAR) g->old_global_surf[ind_sur] = gn.global_surf_id;
      g->nodeID_2d_to_3d_sur[ind_sur] = ind;
      g->nodeID_3d_to_2d_sur[ind] = ind_sur;
    }
    if((gn.global_bed_id != UNSET_INT) && (mod->grid->ndim==3)) {
      ind_bed = node_new_bed(mod);
      
      if (g->type == COLUMNAR) g->old_global_bed[ind_bed] = gn.global_bed_id;
      g->nodeID_2d_to_3d_bed[ind_bed] = ind;
      g->nodeID_3d_to_2d_bed[ind] = ind_bed;
    }
    g->node[ind].global_surf_id = gn.global_surf_id;
    g->node[ind].global_bed_id = gn.global_bed_id;
    g->node[ind].id = ind;
    g->node[ind].string = gn.string;
    g->node[ind].edge_string = gn.edge_string;
    g->node[ind].original_id = gn.original_id;
    g->node[ind].parent[0] = gn.parent[0];
    g->node[ind].parent[1] = gn.parent[1];
    g->node[ind].level = gn.level;
    g->node[ind].block = gn.block;
    g->node[ind].els_flag = gn.els_flag;
    g->node[ind].x = gn.x;
    g->node[ind].y = gn.y;
    g->node[ind].z = gn.z;
    g->node[ind].gid = gn.gid;
    g->node[ind].resident_pe = gn.resident_pe;
    g->node[ind].resident_id = gn.resident_id;
    node_hash_add_entry(g->node[ind], node_hashtab, ind, g->smpi->npes);
  }
  /* if the node was found then the node is in the NODE_LIST_ITEM */
  else
    ind = nd_list_item;
  /* returns the local node number */
  return (ind);
}
#else
void node_get_local(void)
{
}
#endif

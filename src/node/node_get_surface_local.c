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
  int ind;                      /* the node */
  SGRID *g = mod->grid;
  int npes = g->smpi->npes;
   /* looks the node up in the hash table */
  nd_list_item = node_hash_lookup(gn, node_hashtab, npes);
 // printf("nd_list_item->local %d\n", nd_list_item->local);
  /* if the node was not found then create the node and 
     enter it in the hash table */
   
  if (nd_list_item < 0) {
   
    ind = node_new(mod);
    g = mod->grid;
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

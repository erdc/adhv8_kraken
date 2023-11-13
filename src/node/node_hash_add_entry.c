/* adds an entry to the node hash table */
#include "global_header.h"

#ifdef _MESSG
void node_hash_add_entry(SNODE global,    /* the global node pair for the entry */
                         NODE_LIST_ITEM ** node_hashtab,    /* the node hash table */
                         int local_node, /* the local node corresponding to the global node pair */
                         int npes
  )
{
  NODE_LIST_ITEM *np;           /* the pointer to the node hash entry */
  int node_hashval;             /* the index into the hash table */

  /* allocate space for the hash entry */
  np = (NODE_LIST_ITEM *) tl_list_alloc(NODE_LIST);

  /* store data in hash table */
  np->rnode = global.resident_id;
  np->sd = global.resident_pe;
  np->local = local_node;

  /* add structure to beginning of linked list */
  node_hashval = node_hash_index(global, npes);
  np->next = node_hashtab[node_hashval];
  node_hashtab[node_hashval] = np;
}
#else
void node_hash_add_entry(void)
{
}
#endif

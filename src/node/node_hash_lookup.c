/* returns pointer to structure containing global node info */

#include "global_header.h"

#ifdef _MESSG
int node_hash_lookup(SNODE global,    /* the global node number to be looked up */
                                 NODE_LIST_ITEM ** node_hashtab, /* the hash table */
                                 int npes
  )
{
  int node_hashval;             /* the index for the global node */
  int local = -1;
  NODE_LIST_ITEM *np;           /* the entry in the hash table */
  
  node_hashval = node_hash_index(global, npes);
  for (np = node_hashtab[node_hashval]; np != NULL; np = np->next){
    if ((global.resident_id == np->rnode) && (global.resident_pe == np->sd)) {
      local = np->local;
      break;    /* found */
    }
  }

  return local;  
}
#else
void node_hash_lookup(void)
{
}
#endif

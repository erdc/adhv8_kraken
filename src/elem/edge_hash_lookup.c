/* returns pointer to structure containing edge info */

#include "global_header.h"

EDGE_LIST_ITEM *edge_hash_lookup(
  int nd1,			/* the edge nodes to be looked up */
  int nd2,
  EDGE_LIST_ITEM ** edge_hashtab	/* the hash table */
)
{
  int edge_hashval;		/* the index for the edge */
  EDGE_LIST_ITEM *ep;		/* the entry in the hash table */

  edge_hashval = edge_hash_index(nd1, nd2);
  for(ep = edge_hashtab[edge_hashval]; ep != NULL; ep = ep->next)
    if((nd1 == ep->nd1) && (nd2 == ep->nd2))
      return ep;		/* found */
  return NULL;			/* not found */
}

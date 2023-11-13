/* returns pointer to structure containing the 1D element */

#include "global_header.h"

ELEM1D_LIST_ITEM *elem1d_hash_lookup(
  int nd1,			/* the first node */
  int nd2,			/* the second node */
  ELEM1D_LIST_ITEM ** elem1d_hashtab	/* the hash table */
)
{
  int elem1d_hashval;		/* the index for the element */
  ELEM1D_LIST_ITEM *ep;		/* the entry in the hash table */

  elem1d_hashval = elem_hash_index(nd1);
  for(ep = elem1d_hashtab[elem1d_hashval]; ep != NULL; ep = ep->next)
    if((nd1 == ep->nd1) && (nd2 == ep->nd2))
      return ep;		/* found */
  return NULL;			/* not found */
}

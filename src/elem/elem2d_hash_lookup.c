/* returns pointer to structure containing the 2D element */

#include "global_header.h"

ELEM2D_LIST_ITEM *elem2d_hash_lookup(
  int nd1,			/* the first node */
  int nd2,			/* the second node */
  int nd3,			/* the third node */
  ELEM2D_LIST_ITEM ** elem2d_hashtab	/* the hash table */
)
{
  int elem2d_hashval;		/* the index for the element */
  ELEM2D_LIST_ITEM *ep;		/* the entry in the hash table */

  elem2d_hashval = elem_hash_index(nd1);
  for(ep = elem2d_hashtab[elem2d_hashval]; ep != NULL; ep = ep->next)
    if((nd1 == ep->nd1) && (nd2 == ep->nd2) && (nd3 == ep->nd3))
      return ep;		/* found */
  return NULL;			/* not found */
}

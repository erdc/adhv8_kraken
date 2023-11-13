/* returns pointer to structure containing the 3D element */

#include "global_header.h"

ELEM3D_LIST_ITEM *elem3d_hash_lookup(
  int nd1,			/* the first node */
  int nd2,			/* the second node */
  int nd3,			/* the third node */
  int nd4,			/* the fourth node */
  ELEM3D_LIST_ITEM ** elem3d_hashtab	/* the hash table */
)
{
  int elem3d_hashval;		/* the index for the element */
  ELEM3D_LIST_ITEM *ep;		/* the entry in the hash table */

  elem3d_hashval = elem_hash_index(nd1);
  for(ep = elem3d_hashtab[elem3d_hashval]; ep != NULL; ep = ep->next)
    if((nd1 == ep->nd1) && (nd2 == ep->nd2) && (nd3 == ep->nd3) && (nd4 == ep->nd4))
      return ep;		/* found */
  return NULL;			/* not found */
}

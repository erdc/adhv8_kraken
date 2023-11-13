/* adds an entry to the 3d element hash table */
#include "global_header.h"

void elem3d_hash_add_entry(
  int nd1,			/* the first node */
  int nd2,			/* the second node */
  int nd3,			/* the third node */
  int nd4,			/* the fourth node */
  ELEM3D_LIST_ITEM ** elem3d_hashtab,	/* the hash table */
  int ielem			/* the element */
)
{
  ELEM3D_LIST_ITEM *ep;		/* the pointer to the element hash entry */
  int elem3d_hashval;		/* the index into the element hash table */

  /* allocate space for the hash entry */
  ep = (ELEM3D_LIST_ITEM *) tl_list_alloc(ELEM3D_LIST);

  /* store data in hash table */
  ep->nd1 = nd1;
  ep->nd2 = nd2;
  ep->nd3 = nd3;
  ep->nd4 = nd4;
  ep->ielem = ielem;

  /* add structure to beginning of linked list */
  elem3d_hashval = elem_hash_index(nd1);
  ep->next = elem3d_hashtab[elem3d_hashval];
  elem3d_hashtab[elem3d_hashval] = ep;
}

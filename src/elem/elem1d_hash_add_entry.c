/* adds an entry to the 1d element hash table */
#include "global_header.h"

void elem1d_hash_add_entry(
  int nd1,			/* the first node */
  int nd2,			/* the second node */
  ELEM1D_LIST_ITEM ** elem1d_hashtab,	/* the hash table */
  int ielem			/* the element */
)
{
  ELEM1D_LIST_ITEM *ep;		/* the pointer to the element hash entry */
  int elem1d_hashval;		/* the index into the element hash table */

  /* allocate space for the hash entry */
  ep = (ELEM1D_LIST_ITEM *) tl_list_alloc(ELEM1D_LIST);

  /* store data in hash table */
  ep->nd1 = nd1;
  ep->nd2 = nd2;
  ep->ielem = ielem;

  /* add structure to beginning of linked list */
  elem1d_hashval = elem_hash_index(nd1);
  ep->next = elem1d_hashtab[elem1d_hashval];
  elem1d_hashtab[elem1d_hashval] = ep;
}

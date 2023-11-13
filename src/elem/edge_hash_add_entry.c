/* adds an entry to the edge hash table */
#include "global_header.h"

void edge_hash_add_entry(
  int nd1,			/* the nodes on the edge to be added to the edge hash table */
  int nd2,
  EDGE_LIST_ITEM ** edge_hashtab	/* the edge hash table */
)
{
  EDGE_LIST_ITEM *ep;		/* the pointer to the edge hash entry */
  int edge_hashval;		/* the index into the edge hash table */

  /* check to see if edge is already entered */
  ep = edge_hash_lookup(nd1, nd2, edge_hashtab);
  if(ep != NULL)
    return;

  /* allocate space for the hash entry */
  ep = (EDGE_LIST_ITEM *) tl_list_alloc(EDGE_LIST);

  /* store data in hash table */
  ep->nd1 = nd1;
  ep->nd2 = nd2;
  ep->new_node = UNSET_INT;
  ep->rank = UNSET_INT;

  /* add structure to beginning of linked list */
  edge_hashval = edge_hash_index(nd1, nd2);
  ep->next = edge_hashtab[edge_hashval];
  edge_hashtab[edge_hashval] = ep;
}

/* compares two edges to determine the longest */
#include "global_header.h"

int tl_long_edge(
  SNODE *node,
  int e1nd1,			/* the first node on the first edge */
  int e1nd2,			/* the second node on the first edge */
  int e2nd1,			/* the first node on the second edge */
  int e2nd2,			/* the second node on the second edge */
  EDGE_LIST_ITEM ** edge_hashtab	/* hash table of the edges */
)
{
  int level1, level2;		/* the levels of the two edges */
  EDGE_LIST_ITEM *edge1, *edge2;	/* the two edges */

  /* compare the levels */
  if(node[e1nd1].level > node[e1nd2].level)
    level1 = node[e1nd1].level;
  else
    level1 = node[e1nd2].level;
  if(node[e2nd1].level > node[e2nd2].level)
    level2 = node[e2nd1].level;
  else
    level2 = node[e2nd2].level;
  if(level1 > level2)
    return (2);
  if(level2 > level1)
    return (1);

  /* look up the ranks of the edges - if the edge is not found, then 
     it is short */
  edge1 = edge_hash_lookup(e1nd1, e1nd2, edge_hashtab);
  edge2 = edge_hash_lookup(e2nd1, e2nd2, edge_hashtab);
  if(edge1 == NULL)
    return (2);
  if(edge2 == NULL)
    return (1);

  /* compare the ranks of the edges */
  if(edge1->rank > edge2->rank)
    return (1);
  else if(edge1->rank < edge2->rank)
    return (2);
  else if(edge1->rank == UNSET_INT)
    return (0);

  /* if you make it this far then there is a serious error */
  tl_error("Could not differentiate between the old edges in tl_long_edge.");
  return (0);
}

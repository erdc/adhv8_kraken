/* this routine splits a 1d elem */

#include "global_header.h"

void elem1d_split(
        SMODEL *mod,
        int ielem,			            /* the element to be split */
        EDGE_LIST_ITEM ** edge_hashtab	/* the hash table of edges */
)
{
  int new_node;			/* the new node */
  int new_node_level;		/* level for the new node */
  int new_elem;			/* the new element */
  int nd1, nd2;			/* the nodes on an edge */
  EDGE_LIST_ITEM *edge_pntr;	/* pointer to the edge */

  /* refines the edge if needed */
  nd1 = mod->grid->elem1d[ielem].nodes[0];
  nd2 = mod->grid->elem1d[ielem].nodes[1];
  edge_pntr = edge_hash_lookup(nd1, nd2, edge_hashtab);
  new_node = adpt_get_node(mod, edge_pntr, edge_hashtab);

  /* gets the new element */
  new_elem = elem1d_new(mod->grid, mod->nalloc_inc);

  /* copies the current element to the new element */
  selem1d_copy(&(mod->grid->elem1d[new_elem]), mod->grid->elem1d[ielem]);

  /* calculate level for the new node */
  new_node_level = elem1d_level(mod->grid->elem1d[ielem]) + 1;

  /* corrects the current elements data */
  mod->grid->elem1d[ielem].nodes[0] = new_node;
  mod->grid->elem1d[ielem].levels[0] = new_node_level;
  get_elem1d_linear_djac_gradPhi(mod->grid, &(mod->grid->elem1d[ielem]));

  /* corrects the new elements data */
  mod->grid->elem1d[new_elem].id = new_elem;
  mod->grid->elem1d[new_elem].nodes[1] = new_node;
  mod->grid->elem1d[new_elem].levels[1] = new_node_level;
  get_elem1d_linear_djac_gradPhi(mod->grid, &(mod->grid->elem1d[new_elem]));
}

/* Here we  find the connection to */
/* the 2d element.  This is a bit tricky.  The hash table is set up first.  It */
/* is setup twice, once with node[0] used to set the index (the first node) and */
/* once node[1] as the index.  The node numbers are sorted into ascending order */
/* and the hashtable is setup.  When the look up is done, the 1d element's nodes */
/* are put in ascending order as well.  This means at most we will have to make only */
/* two lookups. */

#include "global_header.h"

void elem1d_find_elem2d_init(SGRID *grid)  {

  ELEM2D_LIST_ITEM *ie_list_item;   /* pointer to an element list item */
  int sort_node[NDPRFC];        /* the sorted node numbers */
  int i;
  int ie;                       /* loop counter over the elements */
  int temp;                     /* a temporary variable */
  ELEM2D_LIST_ITEM *elem2d_hashtab[HASHSIZE];   /* tops of the linked lists for the 3d element hash table */

  /* this finds the 2d elements that share an edge with the 1d elements */
  /* this is presently only done for groundwater using overland flow */

  if (grid->nelems1d > 0) {
    for (i = 0; i < HASHSIZE; i++)
      elem2d_hashtab[i] = NULL;
    ie_list_item = NULL;
    /* setup the hash table for the 2d elements, do it with two different indices */
    for (ie = 0; ie < grid->nelems2d; ie++) {

#ifdef _ADH_GROUNDWATER
        if (grid->elem2d[ie].bflag == 0 || grid->elem2d[ie].bflag == UNSET_INT){
            //printf("\nHashing elem2d[%4i](%4i, %4i, %4i)",ie, grid->elem2d[ie].nodes[0], grid->elem2d[ie].nodes[1], grid->elem2d[ie].nodes[2]);
#endif
        /* sorts the nodes */
        ELEM2D_LOCAL_COPY(grid->elem2d[ie].nodes, sort_node);

        if (sort_node[0] > sort_node[1]) {
          temp = sort_node[0];
          sort_node[0] = sort_node[1];
          sort_node[1] = temp;
        }
        if (sort_node[1] > sort_node[2]) {
          temp = sort_node[1];
          sort_node[1] = sort_node[2];
          sort_node[2] = temp;
        }

        /* the first two nodes are now the two smallest nodes */
        elem2d_hash_add_entry(sort_node[0], sort_node[1], sort_node[2], elem2d_hashtab, ie);
        elem2d_hash_add_entry(sort_node[1], sort_node[2], sort_node[0], elem2d_hashtab, ie);

#ifdef _ADH_GROUNDWATER
        } //if (grid->elem2d[ie].bflag == 0){
#endif
    }

    /* now find the 2d element that shares a face with the 1d element */
    for (ie = 0; ie < grid->max_nelems1d; ie++) {

      if (grid->elem1d[ie].string != UNSET_INT) {
          ELEM1D_LOCAL_COPY(grid->elem1d[ie].nodes, sort_node);
          if (sort_node[0] < sort_node[1])
            ie_list_item = elem1d_find_elem2d(sort_node[0], sort_node[1], elem2d_hashtab);
          else
            ie_list_item = elem1d_find_elem2d(sort_node[1], sort_node[0], elem2d_hashtab);

          if (ie_list_item == NULL) {
            tl_error("elem1d_renumber: 2d element does not exist to match 1d element.");
          }
          /* if the element does exist */
          else {
            grid->elem1d[ie].elem2d = ie_list_item->ielem;
            grid->elem1d[ie].mat = grid->elem2d[ grid->elem1d[ie].elem2d ].mat;  /* set material type */
          }
      }
    }
    if (grid->nelems2d != 0)
      tl_list_free_all(ELEM2D_LIST);
  }
}

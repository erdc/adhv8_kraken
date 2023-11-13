/* this routine unpacks the integer information for the given element in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void elem1d_unpacki(
  void *int_buffer,			/* the buffer to be packed */
  void *dble_buffer,     /* the buffer to be packed */
  NODE_LIST_ITEM ** node_hashtab,	/* hash table for the nodes */
  ELEM1D_LIST_ITEM ** elem1d_hashtab,	/* hash table for the 1D elements */
  SMODEL *mod
)
{
  SNODE gn;		/* global node temp used to lookup the local 
				   node number */
  ELEM1D_LIST_ITEM *ie_list_item;	/* pointer to an element list item */
  int nd1, nd2;			/* the nodes in the element */
  int ie;			/* the element */
  int lev1, lev2;		/* the element node levels */
  int imat;			/* the element flag */
  int ibuff = 0, dbuff = 0;		/* position in the buffer */
  int *ibuffer;			/* the integer buffer */
  double *dbuffer;     /* the double buffer */
  int string;
  int elem2d;
  SVECT2D nrml;
  double djac;
  double grad_shp[NDPREDG];
  SGRID *grid =  mod->grid;

  snode_init(&gn);

  /* cast the buffer */
  ibuffer = (int *)int_buffer;
  dbuffer = (double *)dble_buffer;

  /* read the element flag */
  string = ibuffer[ibuff++];

  /* read and process the first node */
  gn.resident_pe = ibuffer[ibuff++];
  gn.resident_id = ibuffer[ibuff++];
  gn.gid = ibuffer[ibuff++];
  gn.global_surf_id = ibuffer[ibuff++];
  gn.global_bed_id = ibuffer[ibuff++];
  gn.string = ibuffer[ibuff++];
  gn.x = dbuffer[dbuff++];
  gn.y = dbuffer[dbuff++];
  gn.z = dbuffer[dbuff++];
  nd1 = node_get_local(gn, node_hashtab, mod);

  /* read and process the second node */
  gn.resident_pe = ibuffer[ibuff++];
  gn.resident_id = ibuffer[ibuff++];
  gn.gid = ibuffer[ibuff++];
  gn.global_surf_id = ibuffer[ibuff++];
  gn.global_bed_id = ibuffer[ibuff++];
  gn.string = ibuffer[ibuff++];
  gn.x = dbuffer[dbuff++];
  gn.y = dbuffer[dbuff++];
  gn.z = dbuffer[dbuff++];
  nd2 = node_get_local(gn, node_hashtab, mod);

  /* read the node levels */
  lev1 = ibuffer[ibuff++];
  lev2 = ibuffer[ibuff++];

  imat = ibuffer[ibuff++];
  elem2d = ibuffer[ibuff++];

  djac = dbuffer[dbuff++];
  grad_shp[0] = dbuffer[dbuff++];
  grad_shp[1] = dbuffer[dbuff++];
  nrml.x = dbuffer[dbuff++];
  nrml.y = dbuffer[dbuff++]; 

  /* look up the element */
  ie_list_item = elem1d_hash_lookup(nd1, nd2, elem1d_hashtab);

  /* if the element does not exist then create it 
     and enter it in the hash table */
  if(ie_list_item == NULL)
    {
      ie = elem1d_new(grid,mod->nalloc_inc);
      grid = mod->grid; 
      grid->elem1d[ie].nodes[0] = nd1;
      grid->elem1d[ie].nodes[1] = nd2;
      grid->elem1d[ie].levels[0] = lev1;
      grid->elem1d[ie].levels[1] = lev2;
      grid->elem1d[ie].mat = imat;
      grid->elem1d[ie].string = string;
      grid->elem1d[ie].elem2d = elem2d;
      grid->elem1d[ie].djac = djac;
      grid->elem1d[ie].grad_shp[0] = grad_shp[0];
      grid->elem1d[ie].grad_shp[1] = grad_shp[1];
      grid->elem1d[ie].nrml.x = nrml.x;
      grid->elem1d[ie].nrml.y = nrml.y;
      elem1d_hash_add_entry(grid->elem1d[ie].nodes[0], grid->elem1d[ie].nodes[1], elem1d_hashtab, ie);
    }
  /* if the element does exist then check the material flag */
  else
    {
      ie = ie_list_item->ielem;
      if(grid->elem1d[ie].mat != imat)
	tl_error("Inconsistent element types found in elem1d_unpacki.");
    }
}
#else
void elem1d_unpacki(
  void
)
{
}
#endif

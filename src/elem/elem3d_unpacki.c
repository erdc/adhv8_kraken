/* this routine unpacks the integer information for the given element in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void elem3d_unpacki(
  void *int_buffer,			/* the integer buffer to be unpacked */
  void *dble_buffer,      /* the double buffer to be packed */
  NODE_LIST_ITEM ** node_hashtab,	/* hash table for the nodes */
  ELEM3D_LIST_ITEM ** elem3d_hashtab,	/* hash table for the 3D elements */
  SMODEL *mod
)
{
  SNODE gn;		/* global node temp used to lookup the local 
				   node number */
  ELEM3D_LIST_ITEM *ie_list_item;	/* pointer to an element list item */
  int nd1, nd2, nd3, nd4;	/* the nodes in the element */
  int lev1, lev2, lev3, lev4;	/* the node levels in the element */
  int ie;			/* the element */
  int imat;			/* the element flag */
  int string; /*the element string number */
  int flx;      /* flux pointer*/
  int ibuff = 0;		/* position in the buffer */
  int dbuff = 0;
  int *ibuffer;			/* the integer buffer */
  double *dbuffer;     /* the double buffer */
  int gid;
  int id_orig;
  int surface2delement;
  double djac;              /* the jacobian of the element */
  double error;               /* max elemental error in continuity equations for the current ts */
  SVECT grad_shp[NDONTET];  /* the gradients of the shape functions */
  SGRID *grid = mod->grid;

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
  gn.string = ibuffer[ibuff++];
  gn.global_surf_id = ibuffer[ibuff++];
  gn.global_bed_id = ibuffer[ibuff++];
  gn.x = dbuffer[dbuff++];
  gn.y = dbuffer[dbuff++];
  gn.z = dbuffer[dbuff++];
  nd1 = node_get_local(gn, node_hashtab, mod);

  /* read and process the second node */
  gn.resident_pe = ibuffer[ibuff++];
  gn.resident_id = ibuffer[ibuff++];
  gn.gid = ibuffer[ibuff++];
  gn.string = ibuffer[ibuff++];
  gn.global_surf_id = ibuffer[ibuff++];
  gn.global_bed_id = ibuffer[ibuff++];
  gn.x = dbuffer[dbuff++];
  gn.y = dbuffer[dbuff++];
  gn.z = dbuffer[dbuff++];
  nd2 = node_get_local(gn, node_hashtab, mod);

  /* read and process the third node */
  gn.resident_pe = ibuffer[ibuff++];
  gn.resident_id = ibuffer[ibuff++];
  gn.gid = ibuffer[ibuff++];
  gn.string = ibuffer[ibuff++];
  gn.global_surf_id = ibuffer[ibuff++];
  gn.global_bed_id = ibuffer[ibuff++];
  gn.x = dbuffer[dbuff++];
  gn.y = dbuffer[dbuff++];
  gn.z = dbuffer[dbuff++];
  nd3 = node_get_local(gn, node_hashtab, mod);

  /* read and process the fourth node */
  gn.resident_pe = ibuffer[ibuff++];
  gn.resident_id = ibuffer[ibuff++];
  gn.gid = ibuffer[ibuff++];
  gn.string = ibuffer[ibuff++];
  gn.global_surf_id = ibuffer[ibuff++];
  gn.global_bed_id = ibuffer[ibuff++];
  gn.x = dbuffer[dbuff++];
  gn.y = dbuffer[dbuff++];
  gn.z = dbuffer[dbuff++];
  nd4 = node_get_local(gn, node_hashtab, mod);

  /* read the node levels */
  lev1 = ibuffer[ibuff++];
  lev2 = ibuffer[ibuff++];
  lev3 = ibuffer[ibuff++];
  lev4 = ibuffer[ibuff++];

  flx = ibuffer[ibuff++];
  id_orig = ibuffer[ibuff++];
 // surface2delement  = ibuffer[ibuff++];
  gid = ibuffer[ibuff++];
  imat = ibuffer[ibuff++];
  djac = dbuffer[dbuff++];
  error = dbuffer[dbuff++];
  grad_shp[0].x = dbuffer[dbuff++];
  grad_shp[0].y = dbuffer[dbuff++];
  grad_shp[0].z = dbuffer[dbuff++];
  grad_shp[1].x = dbuffer[dbuff++];
  grad_shp[1].y = dbuffer[dbuff++];
  grad_shp[1].z = dbuffer[dbuff++];
  grad_shp[2].x = dbuffer[dbuff++];
  grad_shp[2].y = dbuffer[dbuff++];
  grad_shp[2].z = dbuffer[dbuff++];
  grad_shp[3].x = dbuffer[dbuff++];
  grad_shp[3].y = dbuffer[dbuff++];
  grad_shp[3].z = dbuffer[dbuff++];


  /* look up the element */
  ie_list_item = elem3d_hash_lookup(nd1, nd2, nd3, nd4, elem3d_hashtab);

  /* if the element does not exist then create it 
     and enter it in the hash table */
  if(ie_list_item == NULL)
    {
      ie = elem3d_new(grid, mod->nalloc_inc, grid->elem3d[0].nnodes); //LP Hack here to get working
      grid = mod->grid;
      grid->elem3d[ie].id = ie;
      grid->elem3d[ie].nodes[0] = nd1;
      grid->elem3d[ie].nodes[1] = nd2;
      grid->elem3d[ie].nodes[2] = nd3;
      grid->elem3d[ie].nodes[3] = nd4;
      grid->elem3d[ie].levels[0] = lev1;
      grid->elem3d[ie].levels[1] = lev2;
      grid->elem3d[ie].levels[2] = lev3;
      grid->elem3d[ie].levels[3] = lev4;
      grid->elem3d[ie].string = string;
      grid->elem3d[ie].flux_ptr = flx;
      grid->elem3d[ie].mat = imat;
      grid->elem3d[ie].id_orig = id_orig;
      grid->elem3d[ie].djac = djac;
      grid->elem3d[ie].error = error;
      grid->elem3d[ie].grad_shp[0].x = grad_shp[0].x;
      grid->elem3d[ie].grad_shp[0].y = grad_shp[0].y;
      grid->elem3d[ie].grad_shp[0].z = grad_shp[0].z;
      grid->elem3d[ie].grad_shp[1].x = grad_shp[1].x;
      grid->elem3d[ie].grad_shp[1].y = grad_shp[1].y;
      grid->elem3d[ie].grad_shp[1].z = grad_shp[1].z;
      grid->elem3d[ie].grad_shp[1].x = grad_shp[2].x;
      grid->elem3d[ie].grad_shp[2].y = grad_shp[2].y;
      grid->elem3d[ie].grad_shp[2].z = grad_shp[2].z;
      grid->elem3d[ie].grad_shp[2].x = grad_shp[3].x;
      grid->elem3d[ie].grad_shp[2].y = grad_shp[3].y;
      grid->elem3d[ie].grad_shp[2].z = grad_shp[3].z;

     // grid->elem3d[ie].surface2delement = surface2delement;
      grid->elem3d[ie].gid = gid;
      elem3d_hash_add_entry(grid->elem3d[ie].nodes[0], grid->elem3d[ie].nodes[1], grid->elem3d[ie].nodes[2], grid->elem3d[ie].nodes[3], elem3d_hashtab, ie);
    }
  /* if the element does exist then check the material flag */
  else
    {
      ie = ie_list_item->ielem;
      if(grid->elem3d[ie].mat != imat){
      printf("ie %d imat %d elem_mat %d \n", ie, imat, grid->elem3d[ie].mat);
	tl_error("Inconsistent element types found in elem3d_unpacki.");
}
    }
}
#else
void elem3d_unpacki(
  void
)
{
}
#endif

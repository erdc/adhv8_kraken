/* this routine unpacks the integer information for the given element in the given buffer */
/* JRC: change the return type from void to int */

#include "global_header.h"

#ifdef _MESSG
void elem2d_unpacki(
  void *int_buffer,			/* the buffer to be packed */
  void *dble_buffer,     /* the buffer to be packed */
  NODE_LIST_ITEM ** node_hashtab,	/* hash table for the nodes */
  ELEM2D_LIST_ITEM ** elem2d_hashtab,	/* hash table for the 2D elements */
  SMODEL *mod
)
{
  SNODE gn;		/* global node temp used to lookup the local 
				   node number */
  ELEM2D_LIST_ITEM *ie_list_item;	/* pointer to an element list item */
  int nd1, nd2, nd3;		/* the nodes in the element */
  int lev1, lev2, lev3;		/* the node levels in the element */
  int ie;			/* the element */
  int imat;			/* the element material */
  int istrn;			/* the element flag */
  int ibuff = 0, dbuff = 0;		/* position in the buffer */
  int *ibuffer;			/* the integer buffer */
  double *dbuffer;     /* the integer buffer */
  int iflux; 
  int id_orig;
  int bflag;
  int id_3d;
  int gid;
  double djac;              /* jacobian given only 2D coordinates */
  double djac3d;            /* the 2D jacobian surface in 3D */
  double djac3d_fixed;      /* the area of the original element in 3D */
  SVECT2D grad_shp[NDONTRI]; /* the gradients of the shape functions */
  SVECT nrml;               /* the normal to the face */
  SGRID *grid = mod->grid;
  
  snode_init(&gn);
  /* cast the buffer */
  ibuffer = (int *)int_buffer;
  dbuffer = (double *)dble_buffer;
  
  /* read the element flag */
  istrn = ibuffer[ibuff++];
  imat = ibuffer[ibuff++];;

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
  

  /* read the node levels */
  lev1 = ibuffer[ibuff++];
  lev2 = ibuffer[ibuff++];
  lev3 = ibuffer[ibuff++];

  id_orig = ibuffer[ibuff++];
  bflag = ibuffer[ibuff++];
  id_3d = ibuffer[ibuff++];
  gid = ibuffer[ibuff++];

  djac = dbuffer[dbuff++];              /* jacobian given only 2D coordinates */
  djac3d = dbuffer[dbuff++];            /* the 2D jacobian surface in 3D */
  djac3d_fixed = dbuffer[dbuff++];      /* the area of the original element in 3D */
  grad_shp[0].x = dbuffer[dbuff++]; /* the gradients of the shape functions */
  grad_shp[0].y = dbuffer[dbuff++];
  grad_shp[1].x = dbuffer[dbuff++];
  grad_shp[1].y = dbuffer[dbuff++];
  grad_shp[2].x = dbuffer[dbuff++];
  grad_shp[2].y = dbuffer[dbuff++];
  nrml.x = dbuffer[dbuff++];
  nrml.y = dbuffer[dbuff++];
  nrml.z = dbuffer[dbuff++];
               /* the normal to the face */

  //iflux = ibuffer[ibuff++];

  /* look up the element */
  ie_list_item = elem2d_hash_lookup(nd1, nd2, nd3, elem2d_hashtab);

  /* if the element does not exist then create it 
     and enter it in the hash table */
  
  if(ie_list_item == NULL)
    { 
      ie = elem2d_new(grid, mod->nalloc_inc, grid->elem2d[0].nnodes); //LP Hacked this to work may need attention 
      grid = mod->grid;
      grid->elem2d[ie].id = ie;
      grid->elem2d[ie].nodes[0] = nd1;
      grid->elem2d[ie].nodes[1] = nd2;
      grid->elem2d[ie].nodes[2] = nd3;
      grid->elem2d[ie].levels[0] = lev1;
      grid->elem2d[ie].levels[1] = lev2;
      grid->elem2d[ie].levels[2] = lev3;
      grid->elem2d[ie].string = istrn;
      //elem2d_fptr[ie] = iflux;
      grid->elem2d[ie].mat = imat;
      grid->elem2d[ie].id_orig = id_orig;
      grid->elem2d[ie].bflag = bflag;
      grid->elem2d[ie].id_3d = id_3d;
      grid->elem2d[ie].gid = gid;
      grid->elem2d[ie].djac = djac;
      grid->elem2d[ie].djac3d = djac3d;
      grid->elem2d[ie].djac3d_fixed = djac3d_fixed;
      grid->elem2d[ie].grad_shp[0].x = grad_shp[0].x;
      grid->elem2d[ie].grad_shp[0].x = grad_shp[0].y;
      grid->elem2d[ie].grad_shp[0].x = grad_shp[1].x;
      grid->elem2d[ie].grad_shp[0].x = grad_shp[1].y;
      grid->elem2d[ie].grad_shp[0].x = grad_shp[2].x;
      grid->elem2d[ie].grad_shp[0].x = grad_shp[2].y;
      grid->elem2d[ie].nrml.x = nrml.x;
      grid->elem2d[ie].nrml.y = nrml.y;
      grid->elem2d[ie].nrml.z = nrml.z;

      elem2d_hash_add_entry(grid->elem2d[ie].nodes[0], grid->elem2d[ie].nodes[1], grid->elem2d[ie].nodes[2], elem2d_hashtab, ie);

      return;
    }
  /* if the element does exist then check the material flag */
  else
    {
      ie = ie_list_item->ielem;
      if(grid->elem2d[ie].mat != imat)
	      tl_error("Inconsistent element types found in elem2d_unpacki.");
      return;
    }
}
#else
int elem2d_unpacki(
  void
)
{
    return -1;
}
#endif

/* this routine packs the integer information for the given element in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void elem3d_packd(
  void *buffer,			/* the buffer to be packed */
  int ie,			/* the element to be packed */
  SGRID *grid
)
{
  double *ibuffer;			/* the integer buffer */
  int i=0; /*counter */
  /* cast the buffer */
  ibuffer = (double *)buffer;


  /* pack the buffer */
  
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[0]].x;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[0]].y;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[0]].z;
  
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[1]].x;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[1]].y;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[1]].z;
  
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[2]].x;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[2]].y;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[2]].z;
  
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[3]].x;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[3]].y;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[3]].z;

  ibuffer[i++] = grid->elem3d[ie].djac;
  ibuffer[i++] = grid->elem3d[ie].error;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[0].x;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[0].y;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[0].z;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[1].x;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[1].y;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[1].z;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[2].x;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[2].y;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[2].z;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[3].x;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[3].y;
  ibuffer[i++] = grid->elem3d[ie].grad_shp[3].z;

  
}
#else
void elem3d_packd(
  void
)
{
}
#endif

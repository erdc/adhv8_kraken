/* this routine packs the integer information for the given element in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void elem1d_packd(
  void *buffer,			/* the buffer to be packed */
  int ie,			/* the element to be packed */
  SGRID *grid
)
{
  double *ibuffer;			/* the integer buffer */

  /* cast the buffer */
  ibuffer = (double *)buffer;

  /* pack the buffer */
  
  ibuffer[0] = grid->node[grid->elem1d[ie].nodes[0]].x;
  ibuffer[1] = grid->node[grid->elem1d[ie].nodes[0]].y;
  ibuffer[2] = grid->node[grid->elem1d[ie].nodes[0]].z;
  
  ibuffer[3] = grid->node[grid->elem1d[ie].nodes[1]].x;
  ibuffer[4] = grid->node[grid->elem1d[ie].nodes[1]].y;
  ibuffer[5] = grid->node[grid->elem1d[ie].nodes[1]].z;

  ibuffer[6] = grid->elem1d[ie].djac;
  ibuffer[7] = grid->elem1d[ie].grad_shp[0];
  ibuffer[8] = grid->elem1d[ie].grad_shp[1];
  ibuffer[9] = grid->elem1d[ie].nrml.x;
  ibuffer[10] = grid->elem1d[ie].nrml.y;
  
}
#else
void elem1d_packd(
  void
)
{
}
#endif

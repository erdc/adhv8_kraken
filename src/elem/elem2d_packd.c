/* this routine packs the integer information for the given element in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void elem2d_packd(
  void *buffer,			/* the buffer to be packed */
  int ie,			/* the element to be packed */
  SGRID *grid
  
)
{
  double *ibuffer;			/* the integer buffer */

  /* cast the buffer */
  ibuffer = (double *)buffer;

  /* pack the buffer */
  
  
  ibuffer[0] = grid->node[grid->elem2d[ie].nodes[0]].x;
  ibuffer[1] = grid->node[grid->elem2d[ie].nodes[0]].y;
  ibuffer[2] = grid->node[grid->elem2d[ie].nodes[0]].z;
  
  ibuffer[3] = grid->node[grid->elem2d[ie].nodes[1]].x;
  ibuffer[4] = grid->node[grid->elem2d[ie].nodes[1]].y;
  ibuffer[5] = grid->node[grid->elem2d[ie].nodes[1]].z;
  
  ibuffer[6] = grid->node[grid->elem2d[ie].nodes[2]].x;
  ibuffer[7] = grid->node[grid->elem2d[ie].nodes[2]].y;
  ibuffer[8] = grid->node[grid->elem2d[ie].nodes[2]].z;

  ibuffer[9] = grid->elem2d[ie].djac;
  ibuffer[10] = grid->elem2d[ie].djac3d;
  ibuffer[11] = grid->elem2d[ie].djac3d_fixed;
  
  ibuffer[12] = grid->elem2d[ie].grad_shp[0].x;
  ibuffer[13] = grid->elem2d[ie].grad_shp[0].y;

  ibuffer[14] = grid->elem2d[ie].grad_shp[1].x;
  ibuffer[15] = grid->elem2d[ie].grad_shp[1].y;

  ibuffer[16] = grid->elem2d[ie].grad_shp[2].x;
  ibuffer[17] = grid->elem2d[ie].grad_shp[2].y;


  ibuffer[18] = grid->elem2d[ie].nrml.x;
  ibuffer[19] = grid->elem2d[ie].nrml.y;
  ibuffer[20] = grid->elem2d[ie].nrml.z;
  

}
#else
void elem2d_packd(
  void
)
{
}
#endif

/* this routine packs the integer information for the given element in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void elem1d_packi(
  void *buffer,			/* the buffer to be packed */
  int ie,			/* the element to be packed */
  SGRID *grid
)
{
  int *ibuffer;			/* the integer buffer */

  /* cast the buffer */
  ibuffer = (int *)buffer;

  /* pack the buffer */
  ibuffer[0] = grid->elem1d[ie].string;
  ibuffer[1] = grid->node[grid->elem1d[ie].nodes[0]].resident_pe;
  ibuffer[2] = grid->node[grid->elem1d[ie].nodes[0]].resident_id;
  ibuffer[3] = grid->node[grid->elem1d[ie].nodes[0]].gid;
  ibuffer[4] = grid->node[grid->elem1d[ie].nodes[0]].global_surf_id;
  ibuffer[5] = grid->node[grid->elem1d[ie].nodes[0]].global_bed_id;
  ibuffer[6] = grid->node[grid->elem1d[ie].nodes[0]].string;
  ibuffer[7] = grid->node[grid->elem1d[ie].nodes[1]].resident_pe;
  ibuffer[8] = grid->node[grid->elem1d[ie].nodes[1]].resident_id;
  ibuffer[9] = grid->node[grid->elem1d[ie].nodes[1]].gid;
  ibuffer[10] = grid->node[grid->elem1d[ie].nodes[1]].global_surf_id;
  ibuffer[11] = grid->node[grid->elem1d[ie].nodes[1]].global_bed_id;
  ibuffer[12] = grid->node[grid->elem1d[ie].nodes[1]].string;
  ibuffer[13] = grid->elem1d[ie].levels[0];
  ibuffer[14] = grid->elem1d[ie].levels[1];
  ibuffer[15] = grid->elem1d[ie].mat;
  ibuffer[16] = grid->elem1d[ie].elem2d;
}
#else
void elem1d_packi(
  void
)
{
}
#endif

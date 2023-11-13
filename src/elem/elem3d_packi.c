/* this routine packs the integer information for the given element in the given buffer */
// CJT :: HardWired for TETS!!
#include "global_header.h"

#ifdef _MESSG
void elem3d_packi(
  void *buffer,			/* the buffer to be packed */
  int ie,			/* the element to be packed */
  SGRID *grid
)
{
  int *ibuffer;			/* the integer buffer */
  int i=0; /*counter */
  /* cast the buffer */
  ibuffer = (int *)buffer;


  /* pack the buffer */
  ibuffer[i++] = grid->elem3d[ie].string;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[0]].resident_pe;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[0]].resident_id;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[0]].gid;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[0]].string;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[0]].global_surf_id;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[0]].global_bed_id;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[1]].resident_pe;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[1]].resident_id;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[1]].gid;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[1]].string;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[1]].global_surf_id;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[1]].global_bed_id;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[2]].resident_pe;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[2]].resident_id;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[2]].gid;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[2]].string;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[2]].global_surf_id;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[2]].global_bed_id;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[3]].resident_pe;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[3]].resident_id;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[3]].gid;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[3]].string;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[3]].global_surf_id;
  ibuffer[i++] = grid->node[grid->elem3d[ie].nodes[3]].global_bed_id;
  ibuffer[i++] = grid->elem3d[ie].levels[0];
  ibuffer[i++] = grid->elem3d[ie].levels[1];
  ibuffer[i++] = grid->elem3d[ie].levels[2];
  ibuffer[i++] = grid->elem3d[ie].levels[3];
  ibuffer[i++] = grid->elem3d[ie].flux_ptr;
  ibuffer[i++] = grid->elem3d[ie].id_orig;
//  ibuffer[15] = grid->elem3d[ie].surface2delement;
  ibuffer[i++] = grid->elem3d[ie].gid;
  ibuffer[i++] = grid->elem3d[ie].mat;
}
#else
void elem3d_packi(
  void
)
{
}
#endif

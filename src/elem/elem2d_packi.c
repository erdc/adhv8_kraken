/* this routine packs the integer information for the given element in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void elem2d_packi(
  void *buffer,			/* the buffer to be packed */
  int ie,			/* the element to be packed */
  SGRID *grid
  
)
{
  int *ibuffer;			/* the integer buffer */
  int ibuff = 0;

  /* cast the buffer */
  ibuffer = (int *)buffer;

  /* pack the buffer */
  ibuffer[ibuff++] = (int)grid->elem2d[ie].string;;
  ibuffer[ibuff++] = (int)grid->elem2d[ie].mat;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[0]].resident_pe;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[0]].resident_id;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[0]].gid;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[0]].string;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[0]].global_surf_id;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[0]].global_bed_id;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[1]].resident_pe;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[1]].resident_id;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[1]].gid;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[1]].string;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[1]].global_surf_id;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[1]].global_bed_id;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[2]].resident_pe;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[2]].resident_id;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[2]].gid;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[2]].string;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[2]].global_surf_id;
  ibuffer[ibuff++] = grid->node[grid->elem2d[ie].nodes[2]].global_bed_id;
  ibuffer[ibuff++] = grid->elem2d[ie].levels[0];
  ibuffer[ibuff++] = grid->elem2d[ie].levels[1];
  ibuffer[ibuff++] = grid->elem2d[ie].levels[2];
  ibuffer[ibuff++] = grid->elem2d[ie].id_orig;
  ibuffer[ibuff++] = grid->elem2d[ie].bflag;
  ibuffer[ibuff++] = grid->elem2d[ie].id_3d;
  ibuffer[ibuff++] = grid->elem2d[ie].gid;

}
#else
void elem2d_packi(
  void
)
{
}
#endif

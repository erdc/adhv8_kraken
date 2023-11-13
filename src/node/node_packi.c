/* this routine packs the integer information for the given node in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void node_packi(void *buffer,   /* the buffer to be packed */
                int ind,         /* the node to be packed */
                SGRID *grid
  )
{
  int *ibuffer;                 /* the integer buffer */
  int ibuff = 0;                /* current counter in the buffer */

  /* cast the buffer */
  ibuffer = (int *) buffer;

  /* packs the buffer */
  ibuffer[ibuff++] = grid->node[ind].string;
  ibuffer[ibuff++] = grid->node[ind].edge_string;
  ibuffer[ibuff++] = grid->node[ind].node_string;
  ibuffer[ibuff++] = grid->node[ind].original_id;
  ibuffer[ibuff++] = grid->node[ind].parent[0];
  ibuffer[ibuff++] = grid->node[ind].parent[1];
  ibuffer[ibuff++] = grid->node[ind].parent_res_pe[0];
  ibuffer[ibuff++] = grid->node[ind].parent_res_pe[1];
  ibuffer[ibuff++] = grid->node[ind].parent_res_id[0];
  ibuffer[ibuff++] = grid->node[ind].parent_res_id[1];
  ibuffer[ibuff++] = grid->node[ind].level;
  ibuffer[ibuff++] = grid->node[ind].block;
  ibuffer[ibuff++] = grid->node[ind].els_flag;
  ibuffer[ibuff++] = grid->node[ind].nelems_connected;  
  ibuffer[ibuff++] = grid->node[ind].gid;
  ibuffer[ibuff++] = grid->node[ind].global_surf_id;
  ibuffer[ibuff++] = grid->node[ind].global_bed_id;
  ibuffer[ibuff++] = grid->node[ind].resident_pe;
  ibuffer[ibuff++] = grid->node[ind].resident_id;
  
}
#else /*  */
void node_packi(void)
{
}
#endif

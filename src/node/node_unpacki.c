/* this routine unpacks the integer information for the given node in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void node_unpacki(void *buffer, /* the buffer to be unpacked */
                  int ind,      /* the local node number */
                  SGRID *grid   
    )
{
  int *ibuffer;                 /* the integer buffer */
  int ibuff = 0;
  /* cast the buffer */
  ibuffer = (int *) buffer;

  /* unpacks the buffer */
  grid->node[ind].string = ibuffer[ibuff++];
  grid->node[ind].edge_string = ibuffer[ibuff++];
  grid->node[ind].node_string = ibuffer[ibuff++];
  grid->node[ind].original_id = ibuffer[ibuff++];
  grid->node[ind].parent[0] = ibuffer[ibuff++];
  grid->node[ind].parent[1] = ibuffer[ibuff++];
  grid->node[ind].parent_res_pe[0] = ibuffer[ibuff++];
  grid->node[ind].parent_res_pe[1] = ibuffer[ibuff++];
  grid->node[ind].parent_res_id[0] = ibuffer[ibuff++];
  grid->node[ind].parent_res_id[1] = ibuffer[ibuff++];
  grid->node[ind].level = ibuffer[ibuff++];
  grid->node[ind].block = ibuffer[ibuff++];
  grid->node[ind].els_flag = ibuffer[ibuff++];
  grid->node[ind].nelems_connected = ibuffer[ibuff++];
  grid->node[ind].gid = ibuffer[ibuff++];
  grid->node[ind].global_surf_id = ibuffer[ibuff++];
  grid->node[ind].global_bed_id = ibuffer[ibuff++];
  grid->node[ind].resident_pe = ibuffer[ibuff++];
  grid->node[ind].resident_id = ibuffer[ibuff++];
  
}
#else /*  */
void node_unpacki(void)
{
}
#endif

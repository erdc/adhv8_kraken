#include "global_header.h"
/* This updates the ghost nodes after new elements have been added.
 * The nodal information may have been passed in the node_pack routines,
 * but new ghost nodes from the elemXX_unpack routines will not have this
 * information */

void comm_update_snode(SGRID *grid) {

  int i, j; //counter
  SVECT *coords;
  double *area;
  int *packed_ints;

  packed_ints = (int *)tl_alloc(sizeof(int), grid->nnodes * 16); //from snode.h but string and gid have been passed
  coords = (SVECT *)tl_alloc(sizeof(SVECT), grid->nnodes);
  area =  (double *)tl_alloc(sizeof(double), grid->nnodes); //might not need it but better safe than sorry
  
  /* pack the arrays */
  for(j=0,i=0;i<grid->nnodes;i++){
    packed_ints[j++] = grid->node[i].level;
    packed_ints[j++] = grid->node[i].block;
    packed_ints[j++] = grid->node[i].original_id;
    packed_ints[j++] = grid->node[i].resident_pe;
    packed_ints[j++] = grid->node[i].resident_id;
    packed_ints[j++] = grid->node[i].edge_string;
    packed_ints[j++] = grid->node[i].els_flag;
    packed_ints[j++] = grid->node[i].parent[0];
    packed_ints[j++] = grid->node[i].parent[1];
    packed_ints[j++] = grid->node[i].parent_res_pe[0];
    packed_ints[j++] = grid->node[i].parent_res_pe[1];
    packed_ints[j++] = grid->node[i].parent_res_id[0];
    packed_ints[j++] = grid->node[i].parent_res_id[1];
    packed_ints[j++] = grid->node[i].nelems_connected;
    packed_ints[j++] = grid->node[i].global_surf_id;
    packed_ints[j++] = grid->node[i].gid;
    coords[i].x = grid->node[i].x;
    coords[i].y = grid->node[i].y;
    coords[i].z = grid->node[i].z;
    area[i] = grid->node[i].area;
  }

  /*update the arrays */
  comm_update_int(packed_ints, 16, grid->smpi);
  comm_update_VECT(coords, grid->smpi);
  comm_update_double(area, 1, grid->smpi);

  /*update the grid struct (only ghost nodes to save time and be safe) */


  for(j=(16*grid->my_nnodes),i=grid->my_nnodes;i<grid->nnodes;i++){
    grid->node[i].level = packed_ints[j++];
    grid->node[i].block = packed_ints[j++];
    grid->node[i].original_id = packed_ints[j++];
    grid->node[i].resident_pe = packed_ints[j++];
    grid->node[i].resident_id = packed_ints[j++];
    grid->node[i].edge_string = packed_ints[j++];
    grid->node[i].els_flag = packed_ints[j++];
    grid->node[i].parent[0] = packed_ints[j++];
    grid->node[i].parent[1] = packed_ints[j++];
    grid->node[i].parent_res_pe[0] = packed_ints[j++];
    grid->node[i].parent_res_pe[1] = packed_ints[j++];
    grid->node[i].parent_res_id[0] = packed_ints[j++];
    grid->node[i].parent_res_id[1] = packed_ints[j++];
    grid->node[i].nelems_connected = packed_ints[j++];
    grid->node[i].global_surf_id = packed_ints[j++];
    grid->node[i].gid = packed_ints[j++];
    grid->node[i].x = coords[i].x;
    grid->node[i].y = coords[i].y;
    grid->node[i].z = coords[i].z;
    grid->node[i].area = area[i];
  }

  packed_ints = (int *)tl_free(sizeof(int), grid->nnodes * 16, packed_ints); //from snode.h but string and gid have been passed
  coords = (SVECT *)tl_free(sizeof(SVECT), grid->nnodes, coords);
  area =  (double *)tl_free(sizeof(double), grid->nnodes, area); //might not need it but better safe than sorry

  return;

}



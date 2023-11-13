#include "global_header.h"

/* This updates the SW3 struct. It needs to be called when ghost nodes have been added after element renumbering */

void comm_update_sw3_surface(SSW_3D *sw, SGRID *grid

    )
{

  int i, j, ind, jj, k; //counters
  int num_doubles_sw3;

  double *packed_sw3;

  num_doubles_sw3 = 11; // from ssw_3d.h no winds, waves, sediment, 2d or temp arrays;
  if (sw->waves != NULL) num_doubles_sw3 += 5;
  if (sw->winds != NULL) num_doubles_sw3 += 2;


  packed_sw3 = (double *)tl_alloc(sizeof(double), grid->nnodes_sur * num_doubles_sw3);

  /* pack the array */

  for(j=0,ind=0;ind<grid->nnodes_sur;ind++){
    packed_sw3[j++] = sw->depth[ind];
    packed_sw3[j++] = sw->old_depth[ind];
    packed_sw3[j++] = sw->bed_elevation[ind];
    packed_sw3[j++] = sw->depth_avg_vel[ind].x;
    packed_sw3[j++] = sw->depth_avg_vel[ind].y;
    packed_sw3[j++] = sw->old_depth_avg_vel[ind].x;
    packed_sw3[j++] = sw->old_depth_avg_vel[ind].y;
    packed_sw3[j++] = sw->surface_vel[ind].x;
    packed_sw3[j++] = sw->surface_vel[ind].y;
    packed_sw3[j++] = sw->bottom_vel[ind].x;
    packed_sw3[j++] = sw->bottom_vel[ind].y;

    if (sw->waves != NULL) {
      packed_sw3[j++] = sw->waves[ind].rads.xx;
      packed_sw3[j++] = sw->waves[ind].rads.yy;
      packed_sw3[j++] = sw->waves[ind].rads.xy;
      packed_sw3[j++] = sw->waves[ind].stress.x;
      packed_sw3[j++] = sw->waves[ind].stress.y;
    }

    if (sw->winds != NULL ) {
      packed_sw3[j++] = sw->winds[ind].stress.x;
      packed_sw3[j++] = sw->winds[ind].stress.y;
    }
  }


  /*update the array */
  comm_update_double_surf(packed_sw3, num_doubles_sw3, grid->smpi);
  

  /* update the sw3 struct (only ghost nodes to save time and be safe)*/
  for(j=0,ind=0;ind<grid->nnodes_sur;ind++){
    sw->depth[ind] = packed_sw3[j++];
    sw->old_depth[ind] = packed_sw3[j++];
    sw->bed_elevation[ind] = packed_sw3[j++];
    sw->depth_avg_vel[ind].x = packed_sw3[j++];
    sw->depth_avg_vel[ind].y = packed_sw3[j++];
    sw->old_depth_avg_vel[ind].x = packed_sw3[j++];
    sw->old_depth_avg_vel[ind].y = packed_sw3[j++];
    sw->surface_vel[ind].x = packed_sw3[j++];
    sw->surface_vel[ind].y = packed_sw3[j++];
    sw->bottom_vel[ind].x = packed_sw3[j++];
    sw->bottom_vel[ind].y = packed_sw3[j++];

    if (sw->waves != NULL) {
      sw->waves[ind].rads.xx = packed_sw3[j++];
      sw->waves[ind].rads.yy = packed_sw3[j++];
      sw->waves[ind].rads.xy = packed_sw3[j++];
      sw->waves[ind].stress.x = packed_sw3[j++];
      sw->waves[ind].stress.y = packed_sw3[j++];
    }

    if (sw->winds != NULL ) {
      sw->winds[ind].stress.x = packed_sw3[j++];
      sw->winds[ind].stress.y = packed_sw3[j++];
    }
  }

  packed_sw3 = (double *)tl_free(sizeof(double), grid->nnodes_sur * num_doubles_sw3, packed_sw3);

  return;

}

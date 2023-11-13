#include "global_header.h"

/* This updates the SW2 struct. It needs to be called when ghost nodes have been added after element renumbering */

void comm_update_sw2(SSW_2D *sw2, SGRID *grid
  ) {

  int ind,i,j, jj,k; //counters
  int num_doubles_sw2 = 13; // from ssw_2d.h no winds, waves or temp arrays;
  double *packed_sw2;

  if(sw2->waves != NULL) num_doubles_sw2 += 5;
  if(sw2->winds != NULL) num_doubles_sw2 += 2;

  packed_sw2 = (double *)tl_alloc(sizeof(double), grid->nnodes * num_doubles_sw2);


  /* pack the array */
  for(j=0,ind=0;ind<grid->nnodes;ind++){
    packed_sw2[j++] = sw2->head[ind];
    packed_sw2[j++] = sw2->old_head[ind];
    packed_sw2[j++] = sw2->older_head[ind];
    packed_sw2[j++] = sw2->vel[ind].x;
    packed_sw2[j++] = sw2->vel[ind].y;
    packed_sw2[j++] = sw2->old_vel[ind].x;
    packed_sw2[j++] = sw2->old_vel[ind].y;
    packed_sw2[j++] = sw2->older_vel[ind].x;
    packed_sw2[j++] = sw2->older_vel[ind].y;
    packed_sw2[j++] = sw2->density[ind];
    packed_sw2[j++] = sw2->bed_displacement[ind];
    packed_sw2[j++] = sw2->bed_elevation[ind];
    packed_sw2[j++] = sw2->error[ind];
    if (sw2->waves != NULL) {
      packed_sw2[j++] = sw2->waves[ind].rads.xx;
      packed_sw2[j++] = sw2->waves[ind].rads.xy;
      packed_sw2[j++] = sw2->waves[ind].rads.yy;
      packed_sw2[j++] = sw2->waves[ind].stress.x;
      packed_sw2[j++] = sw2->waves[ind].stress.y;
    }

    if (sw2->winds != NULL) {
      packed_sw2[j++] = sw2->winds[ind].stress.x;
      packed_sw2[j++] = sw2->winds[ind].stress.y;
    }
  }
  /*update the array */
  comm_update_double(packed_sw2, num_doubles_sw2, grid->smpi);
  /* update the sw2 struct (only ghost nodes to save time and be safe)*/
  for(j=(grid->my_nnodes*num_doubles_sw2),ind=grid->my_nnodes;ind<grid->nnodes;ind++){
    sw2->head[ind] = packed_sw2[j++];
    sw2->old_head[ind] = packed_sw2[j++];
    sw2->older_head[ind] = packed_sw2[j++];
    sw2->vel[ind].x = packed_sw2[j++];
    sw2->vel[ind].y = packed_sw2[j++];
    sw2->old_vel[ind].x = packed_sw2[j++];
    sw2->old_vel[ind].y = packed_sw2[j++];
    sw2->older_vel[ind].x = packed_sw2[j++];
    sw2->older_vel[ind].y = packed_sw2[j++];
    sw2->density[ind] = packed_sw2[j++];
    sw2->bed_displacement[ind] = packed_sw2[j++];
    sw2->bed_elevation[ind] = packed_sw2[j++];
    sw2->error[ind] = packed_sw2[j++];
    if (sw2->waves != NULL) {
      sw2->waves[ind].rads.xx = packed_sw2[j++];
      sw2->waves[ind].rads.xy = packed_sw2[j++];
      sw2->waves[ind].rads.yy = packed_sw2[j++];
      sw2->waves[ind].stress.x = packed_sw2[j++];
      sw2->waves[ind].stress.y = packed_sw2[j++];
    }

    if (sw2->winds != NULL) {
      sw2->winds[ind].stress.x = packed_sw2[j++];
      sw2->winds[ind].stress.y = packed_sw2[j++];
    }
  }

  packed_sw2 = (double *)tl_free(sizeof(double), grid->nnodes * num_doubles_sw2, packed_sw2);

  return;

}


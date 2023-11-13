#include "global_header.h"

/* This updates the SW3 struct. It needs to be called when ghost nodes have been added after element renumbering */

void comm_update_sw3(SSW_3D *sw, SGRID *grid) {

  int i,j; //counters
  int num_doubles_sw3;
  int num_doubles_comm;
  double *packed_sw3;

  num_doubles_sw3 = 25*grid->nnodes; // from ssw_3d.h no winds, waves, sediment, 2d or temp arrays;
  num_doubles_comm = 25;
#ifdef _SEDIMENT
  num_doubles_sw3 += (3*grid->nnodes);
  num_doubles_comm += 3;
#endif
  /* add 2d variables */
/*  num_doubles_sw3 += (11*grid->nnodes_sur);

  if(sw3->waves != NULL) num_doubles_sw3 += (5*grid->nnodes_sur);
  if(sw3->winds != NULL) num_doubles_sw3 += (2*grid->nnodes_sur);
*/
  packed_sw3 = (double *)tl_alloc(sizeof(double), num_doubles_sw3); 

  /* pack the array */
  
  for(j=0,i=0;i<grid->nnodes;i++){
      packed_sw3[j++] = sw->displacement[i];
      packed_sw3[j++] = sw->old_displacement[i];
      packed_sw3[j++] = sw->older_displacement[i];
      packed_sw3[j++] = sw->grid_speed[i];
      packed_sw3[j++] = sw->old_grid_speed[i];
      packed_sw3[j++] = sw->density[i];
      packed_sw3[j++] = sw->dpl_perturbation[i];
      packed_sw3[j++] = sw->prs[i];
      packed_sw3[j++] = sw->prs_plus[i];
      packed_sw3[j++] = sw->prs_minus[i];
      packed_sw3[j++] = sw->error[i];
      packed_sw3[j++] = sw->vertical_node_flux[i];
      packed_sw3[j++] = sw->hyd_viscosity[i]; 
      packed_sw3[j++] = sw->trn_diffusivity[i]; 

      packed_sw3[j++] = sw->vel[i].x;
      packed_sw3[j++] = sw->vel[i].y;
      packed_sw3[j++] = sw->vel[i].z;
      packed_sw3[j++] = sw->old_vel[i].x;
      packed_sw3[j++] = sw->old_vel[i].y;
      packed_sw3[j++] = sw->old_vel[i].z;
      packed_sw3[j++] = sw->older_vel[i].x;
      packed_sw3[j++] = sw->older_vel[i].y;
      packed_sw3[j++] = sw->older_vel[i].z;
      packed_sw3[j++] = sw->tanvec[i].x;
      packed_sw3[j++] = sw->tanvec[i].y;

#ifdef _SEDIMENT
      packed_sw3[j++] = sw->bed_displacement[i];
      packed_sw3[j++] = sw->old_bed_displacement[i];
      packed_sw3[j++] = sw->older_bed_displacement[i];
#endif
   /*   if(grid->nodeID_3d_to_2d_sur[i] != UNSET_INT){
        packed_sw3[j++] = sw->depth[grid->nodeID_3d_to_2d_sur[i]];
        packed_sw3[j++] = sw->old_depth[grid->nodeID_3d_to_2d_sur[i]];
        packed_sw3[j++] = sw->bed_elevation[grid->nodeID_3d_to_2d_sur[i]];
        packed_sw3[j++] = sw->depth_avg_vel[grid->nodeID_3d_to_2d_sur[i]].x;
        packed_sw3[j++] = sw->depth_avg_vel[grid->nodeID_3d_to_2d_sur[i]].y;
        packed_sw3[j++] = sw->old_depth_avg_vel[grid->nodeID_3d_to_2d_sur[i]].x;
        packed_sw3[j++] = sw->old_depth_avg_vel[grid->nodeID_3d_to_2d_sur[i]].y;
        packed_sw3[j++] = sw->surface_vel[grid->nodeID_3d_to_2d_sur[i]].x;
        packed_sw3[j++] = sw->surface_vel[grid->nodeID_3d_to_2d_sur[i]].y;
        packed_sw3[j++] = sw->bottom_vel[grid->nodeID_3d_to_2d_sur[i]].x;
        packed_sw3[j++] = sw->bottom_vel[grid->nodeID_3d_to_2d_sur[i]].y;
        if(sw->winds != NULL){
          packed_sw3[j++] = sw->winds[grid->nodeID_3d_to_2d_sur[i]].x;
          packed_sw3[j++] = sw->winds[grid->nodeID_3d_to_2d_sur[i]].y;
        }
        if(sw->waves != NULL){
          packed_sw3[j++] = sw->waves[grid->nodeID_3d_to_2d_sur[i]].rads.xx;
          packed_sw3[j++] = sw->waves[grid->nodeID_3d_to_2d_sur[i]].rads.xy;
          packed_sw3[j++] = sw->waves[grid->nodeID_3d_to_2d_sur[i]].rads.yy;
          packed_sw3[j++] = sw->waves[grid->nodeID_3d_to_2d_sur[i]].stress.x;
          packed_sw3[j++] = sw->waves[grid->nodeID_3d_to_2d_sur[i]].stress.y;
        }
      }
      */
  }

  /*update the array */
  comm_update_double(packed_sw3, num_doubles_sw3, grid->smpi);

  /* update the sw2 struct (only ghost nodes to save time and be safe)*/
  for(j=0,i=0;i<grid->nnodes;i++){
    if(i<grid->my_nnodes){
      j += 25;
#ifdef _SEDIMENT
      j += 3;
#endif
  /*    if(grid->nodeID_3d_to_2d_sur[i] != UNSET_INT){
        j += 11;
        if(sw->winds != NULL) j +=2;
        if(sw->winds != NULL) j +=5;
      }*/
    }
    else{
      sw->displacement[i] = packed_sw3[j++];
      sw->old_displacement[i] = packed_sw3[j++];
      sw->older_displacement[i] = packed_sw3[j++];
      sw->grid_speed[i] = packed_sw3[j++];
      sw->old_grid_speed[i] = packed_sw3[j++];
      sw->density[i] = packed_sw3[j++];
      sw->dpl_perturbation[i] = packed_sw3[j++];
      sw->prs[i] = packed_sw3[j++];
      sw->prs_plus[i] = packed_sw3[j++];
      sw->prs_minus[i] = packed_sw3[j++];
      sw->error[i] = packed_sw3[j++];
      sw->vertical_node_flux[i] = packed_sw3[j++];
      sw->hyd_viscosity[i] = packed_sw3[j++];
      sw->trn_diffusivity[i] = packed_sw3[j++];

      sw->vel[i].x = packed_sw3[j++];
      sw->vel[i].y = packed_sw3[j++];
      sw->vel[i].z = packed_sw3[j++];
      sw->old_vel[i].x = packed_sw3[j++];
      sw->old_vel[i].y = packed_sw3[j++];
      sw->old_vel[i].z = packed_sw3[j++];
      sw->older_vel[i].x = packed_sw3[j++];
      sw->older_vel[i].y = packed_sw3[j++];
      sw->older_vel[i].z = packed_sw3[j++];
      sw->tanvec[i].x = packed_sw3[j++];
      sw->tanvec[i].y = packed_sw3[j++];

#ifdef _SEDIMENT
      sw->bed_displacement[i] = packed_sw3[j++];
      sw->old_bed_displacement[i] = packed_sw3[j++];
      sw->older_bed_displacement[i] = packed_sw3[j++];
#endif
  /*    if(grid->nodeID_3d_to_2d_sur[i] != UNSET_INT){
        sw->depth[grid->nodeID_3d_to_2d_sur[i]] = packed_sw3[j++];
        sw->old_depth[grid->nodeID_3d_to_2d_sur[i]] = packed_sw3[j++];
        sw->bed_elevation[grid->nodeID_3d_to_2d_sur[i]] = packed_sw3[j++];
        sw->depth_avg_vel[grid->nodeID_3d_to_2d_sur[i]].x = packed_sw3[j++];
        sw->depth_avg_vel[grid->nodeID_3d_to_2d_sur[i]].y = packed_sw3[j++];
        sw->old_depth_avg_vel[grid->nodeID_3d_to_2d_sur[i]].x = packed_sw3[j++];
        sw->old_depth_avg_vel[grid->nodeID_3d_to_2d_sur[i]].y = packed_sw3[j++];
        sw->surface_vel[grid->nodeID_3d_to_2d_sur[i]].x = packed_sw3[j++];
        sw->surface_vel[grid->nodeID_3d_to_2d_sur[i]].y = packed_sw3[j++];
        sw->bottom_vel[grid->nodeID_3d_to_2d_sur[i]].x = packed_sw3[j++];
        sw->bottom_vel[grid->nodeID_3d_to_2d_sur[i]].y = packed_sw3[j++];
        if(sw->winds != NULL){
          sw->winds[grid->nodeID_3d_to_2d_sur[i]].x = packed_sw3[j++];
          sw->winds[grid->nodeID_3d_to_2d_sur[i]].y = packed_sw3[j++];
        }
        if(sw->waves != NULL){
          sw->waves[grid->nodeID_3d_to_2d_sur[i]].rads.xx = packed_sw3[j++];
          sw->waves[grid->nodeID_3d_to_2d_sur[i]].rads.xy = packed_sw3[j++];
          sw->waves[grid->nodeID_3d_to_2d_sur[i]].rads.yy = packed_sw3[j++];
          sw->waves[grid->nodeID_3d_to_2d_sur[i]].stress.x = packed_sw3[j++];
          sw->waves[grid->nodeID_3d_to_2d_sur[i]].stress.y = packed_sw3[j++];
        }
      }
*/
  }


  }

  packed_sw3 = (double *)tl_free(sizeof(double), num_doubles_sw3, packed_sw3);

  return;

}

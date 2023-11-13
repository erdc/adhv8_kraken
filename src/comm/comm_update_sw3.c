#include "global_header.h"

/* This updates the SW3 struct. It needs to be called when ghost nodes have been added after element renumbering */

void comm_update_sw3(SSW_3D *sw, SGRID *grid
          ) {

  int i,j,k; //counters
  int num_doubles_sw3;
  
  double *packed_sw3;

  num_doubles_sw3 = 25; // from ssw_3d.h no winds, waves, sediment, 2d or temp arrays;
  
#ifdef _SEDIMENT
  num_doubles_sw3 += 3;
  
#endif

  packed_sw3 = (double *)tl_alloc(sizeof(double), grid->nnodes * num_doubles_sw3); 

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
  }

  /*update the array */
  comm_update_double(packed_sw3, num_doubles_sw3, grid->smpi);

  /* update the sw3 struct (only ghost nodes to save time and be safe)*/
  for(j=0,i=0;i<grid->nnodes;i++){
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
  }

  packed_sw3 = (double *)tl_free(sizeof(double), grid->nnodes * num_doubles_sw3, packed_sw3);

  return;

}

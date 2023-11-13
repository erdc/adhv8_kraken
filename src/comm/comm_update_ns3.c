#include "global_header.h"

/* This updates the SW3 struct. It needs to be called when ghost nodes have been added after element renumbering */

void comm_update_ns3(SNS_3D *ns, SGRID *grid
          ) {

  int i,j,k; //counters
  int num_doubles_ns3;
  
  double *packed_ns3;

  num_doubles_ns3 = 25; // from sns_3d.h no winds, waves, sediment, 2d or temp arrays;
  
#ifdef _SEDIMENT
  num_doubles_ns3 += 3;
  
#endif

  packed_ns3 = (double *)tl_alloc(sizeof(double), grid->nnodes * num_doubles_ns3);

  /* pack the array */
  
  for(j=0,i=0;i<grid->nnodes;i++){
      packed_ns3[j++] = ns->displacement[i];
      packed_ns3[j++] = ns->old_displacement[i];
      packed_ns3[j++] = ns->older_displacement[i];
      packed_ns3[j++] = ns->grid_speed[i];
      packed_ns3[j++] = ns->old_grid_speed[i];
      packed_ns3[j++] = ns->density[i];
      packed_ns3[j++] = ns->dpl_perturbation[i];
      packed_ns3[j++] = ns->prs[i];
      packed_ns3[j++] = ns->old_prs[i];
      packed_ns3[j++] = ns->older_prs[i];
      packed_ns3[j++] = ns->error[i];
      packed_ns3[j++] = ns->vertical_node_flux[i];
      packed_ns3[j++] = ns->hyd_viscosity[i];
      packed_ns3[j++] = ns->trn_diffusivity[i];

      packed_ns3[j++] = ns->vel[i].x;
      packed_ns3[j++] = ns->vel[i].y;
      packed_ns3[j++] = ns->vel[i].z;
      packed_ns3[j++] = ns->old_vel[i].x;
      packed_ns3[j++] = ns->old_vel[i].y;
      packed_ns3[j++] = ns->old_vel[i].z;
      packed_ns3[j++] = ns->older_vel[i].x;
      packed_ns3[j++] = ns->older_vel[i].y;
      packed_ns3[j++] = ns->older_vel[i].z;
      packed_ns3[j++] = ns->tanvec[i].x;
      packed_ns3[j++] = ns->tanvec[i].y;

#ifdef _SEDIMENT
      packed_ns3[j++] = ns->bed_displacement[i];
      packed_ns3[j++] = ns->old_bed_displacement[i];
      packed_ns3[j++] = ns->older_bed_displacement[i];

#endif
  }

  /*update the array */
  comm_update_double(packed_ns3, num_doubles_ns3, grid->smpi);

  /* update the ns3 struct (only ghost nodes to save time and be safe)*/
  for(j=0,i=0;i<grid->nnodes;i++){
      ns->displacement[i] = packed_ns3[j++];
      ns->old_displacement[i] = packed_ns3[j++];
      ns->older_displacement[i] = packed_ns3[j++];
      ns->grid_speed[i] = packed_ns3[j++];
      ns->old_grid_speed[i] = packed_ns3[j++];
      ns->density[i] = packed_ns3[j++];
      ns->dpl_perturbation[i] = packed_ns3[j++];
      ns->prs[i] = packed_ns3[j++];
      ns->old_prs[i] = packed_ns3[j++];
      ns->older_prs[i] = packed_ns3[j++];
      ns->error[i] = packed_ns3[j++];
      ns->vertical_node_flux[i] = packed_ns3[j++];
      ns->hyd_viscosity[i] = packed_ns3[j++];
      ns->trn_diffusivity[i] = packed_ns3[j++];

      ns->vel[i].x = packed_ns3[j++];
      ns->vel[i].y = packed_ns3[j++];
      ns->vel[i].z = packed_ns3[j++];
      ns->old_vel[i].x = packed_ns3[j++];
      ns->old_vel[i].y = packed_ns3[j++];
      ns->old_vel[i].z = packed_ns3[j++];
      ns->older_vel[i].x = packed_ns3[j++];
      ns->older_vel[i].y = packed_ns3[j++];
      ns->older_vel[i].z = packed_ns3[j++];
      ns->tanvec[i].x = packed_ns3[j++];
      ns->tanvec[i].y = packed_ns3[j++];

#ifdef _SEDIMENT
      ns->bed_displacement[i] = packed_ns3[j++];
      ns->old_bed_displacement[i] = packed_ns3[j++];
      ns->older_bed_displacement[i] = packed_ns3[j++];

#endif
  }


  

  packed_ns3 = (double *)tl_free(sizeof(double), grid->nnodes * num_doubles_ns3, packed_ns3);

  return;

}

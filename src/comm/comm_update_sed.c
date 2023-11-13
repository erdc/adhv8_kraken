#include "global_header.h"

/* This updates the SED struct. It needs to be called when ghost nodes have been added after element renumbering */

void comm_update_sed(SSED *sed, SGRID *grid, int nconti, int ngrains, int nlayers) {

  int ind,i,j,jj,k; //counters
  int num_doubles_sed,num_doubles_sed_bed;
	int nnodes,nnodes_bed;


  double *packed_sed, *packed_sed_bed;

  nnodes=grid->nnodes;

  if(grid->ndim==2){
    num_doubles_sed = 2 + (11 * ngrains);
    nnodes_bed = grid->nnodes;
  }else{
    num_doubles_sed = 2 + (6 * ngrains);
    nnodes_bed = grid->nnodes_bed;
  }

  num_doubles_sed_bed = 22 + (17 * ngrains) + (5 * nconti) + (2* ngrains * nlayers) + (12 * nlayers);

	packed_sed = (double *)tl_alloc(sizeof(double), nnodes * num_doubles_sed);
	packed_sed_bed = (double *)tl_alloc(sizeof(double), nnodes_bed * num_doubles_sed_bed);

  /* pack the array */

  for(j=0,ind=0;ind<grid->nnodes;ind++){
    packed_sed[j++] = sed->susload_vector[ind].x;
    packed_sed[j++] = sed->susload_vector[ind].y;
    for(k=0;k<ngrains;k++){
      packed_sed[j++] = sed->susload[k].c[ind];
      packed_sed[j++] = sed->susload[k].old_c[ind];
      packed_sed[j++] = sed->susload[k].older_c[ind];
      packed_sed[j++] = sed->susload[k].source[ind];
      packed_sed[j++] = sed->susload[k].sink[ind];
      packed_sed[j++] = sed->susload[k].error[ind];
			if(grid->ndim == 2){
				packed_sed[j++] = sed->susload[i].mfcf[ind];
      	packed_sed[j++] = sed->susload[i].vcf[ind].x;
      	packed_sed[j++] = sed->susload[i].vcf[ind].y;
      	packed_sed[j++] = sed->susload[i].vor_vel[ind].x;
      	packed_sed[j++] = sed->susload[i].vor_vel[ind].y;
			}
    }
  }

  for(j=0,ind=0;ind<nnodes_bed;ind++){
    packed_sed_bed[j++] = sed->friction_coef[ind];
    packed_sed_bed[j++] = sed->as_ceiling[ind];
    packed_sed_bed[j++] = sed->old_as_ceiling[ind];
    packed_sed_bed[j++] = sed->bed_displacement[ind];
    packed_sed_bed[j++] = sed->old_bed_displacement[ind];
    packed_sed_bed[j++] = sed->bed_shear_stress[ind];
    packed_sed_bed[j++] = sed->bedload_vector[ind].x;
    packed_sed_bed[j++] = sed->bedload_vector[ind].y;
    packed_sed_bed[j++] = sed->bed_gradient[ind].x;
    packed_sed_bed[j++] = sed->bed_gradient[ind].y;
    for(k=0;k<5;k++){
      for(jj=0;jj<nconti;jj++){
        packed_sed_bed[j++] = sed->consolidation_bed_properties[k][jj][ind];
      }
    }
    for(i=0;i<ngrains;i++){
      packed_sed_bed[j++] = sed->susload[i].rouse_coef[ind];
      packed_sed_bed[j++] = sed->bedload[i].c[ind];
      packed_sed_bed[j++] = sed->bedload[i].old_c[ind];
      packed_sed_bed[j++] = sed->bedload[i].older_c[ind];
      packed_sed_bed[j++] = sed->bedload[i].thick[ind];
      packed_sed_bed[j++] = sed->bedload[i].old_thick[ind];
      packed_sed_bed[j++] = sed->bedload[i].older_thick[ind];
      packed_sed_bed[j++] = sed->bedload[i].source[ind];
      packed_sed_bed[j++] = sed->bedload[i].sink[ind];
      packed_sed_bed[j++] = sed->bedload[i].shear[ind];
      packed_sed_bed[j++] = sed->bedload[i].error[ind];
      packed_sed_bed[j++] = sed->bedload[i].v[ind].x;
      packed_sed_bed[j++] = sed->bedload[i].v[ind].y;
      packed_sed_bed[j++] = sed->bedload[i].flux[ind].x;
      packed_sed_bed[j++] = sed->bedload[i].flux[ind].y;
      packed_sed_bed[j++] = sed->active_layer[0].distribution[i][ind];
      packed_sed_bed[j++] = sed->old_active_layer[0].distribution[i][ind];
      for(jj=0;jj<nlayers;jj++){
        packed_sed_bed[j++] = sed->bed_layer[jj].distribution[i][ind];
        packed_sed_bed[j++] = sed->old_bed_layer[jj].distribution[i][ind];
      }
    }
    packed_sed_bed[j++] = sed->active_layer[0].bulk_density[ind];
    packed_sed_bed[j++] = sed->active_layer[0].porosity[ind];
    packed_sed_bed[j++] = sed->active_layer[0].critical_erosion_shear[ind];
    packed_sed_bed[j++] = sed->active_layer[0].erosion_rate_constant[ind];
    packed_sed_bed[j++] = sed->active_layer[0].erosion_rate_exponent[ind];
    packed_sed_bed[j++] = sed->active_layer[0].thickness[ind];
    //packed_sed_bed[j++] = sed->active_layer[0].sedflume.shear_stress[ind];
    //packed_sed_bed[j++] = sed->active_layer[0].sedflume.erosion_rate[ind];
    packed_sed_bed[j++] = sed->old_active_layer[0].bulk_density[ind];
    packed_sed_bed[j++] = sed->old_active_layer[0].porosity[ind];
    packed_sed_bed[j++] = sed->old_active_layer[0].critical_erosion_shear[ind];
    packed_sed_bed[j++] = sed->old_active_layer[0].erosion_rate_constant[ind];
    packed_sed_bed[j++] = sed->old_active_layer[0].erosion_rate_exponent[ind];
    packed_sed_bed[j++] = sed->old_active_layer[0].thickness[ind];
    //packed_sed_bed[j++] = sed->old_active_layer[0].sedflume.shear_stress[ind];
    //packed_sed_bed[j++] = sed->old_active_layer[0].sedflume.erosion_rate[ind];
    for (jj=0;jj<nlayers;jj++){
      packed_sed_bed[j++] = sed->bed_layer[jj].bulk_density[ind];
      packed_sed_bed[j++] = sed->bed_layer[jj].porosity[ind];
      packed_sed_bed[j++] = sed->bed_layer[jj].critical_erosion_shear[ind];
      packed_sed_bed[j++] = sed->bed_layer[jj].erosion_rate_constant[ind];
      packed_sed_bed[j++] = sed->bed_layer[jj].erosion_rate_exponent[ind];
      packed_sed_bed[j++] = sed->bed_layer[jj].thickness[ind];
      //packed_sed_bed[j++] = sed->bed_layer[j].sedflume.shear_stress[ind];
      //packed_sed_bed[j++] = sed->bed_layer[j].sedflume.erosion_rate[ind];
      packed_sed_bed[j++] = sed->old_bed_layer[jj].bulk_density[ind];
      packed_sed_bed[j++] = sed->old_bed_layer[jj].porosity[ind];
      packed_sed_bed[j++] = sed->old_bed_layer[jj].critical_erosion_shear[ind];
      packed_sed_bed[j++] = sed->old_bed_layer[jj].erosion_rate_constant[ind];
      packed_sed_bed[j++] = sed->old_bed_layer[jj].erosion_rate_exponent[ind];
      packed_sed_bed[j++] = sed->old_bed_layer[jj].thickness[ind];
    }
  }

  if(grid->ndim==2){
    comm_update_double(packed_sed_bed, num_doubles_sed_bed, grid->smpi);
    comm_update_double(packed_sed, num_doubles_sed, grid->smpi);
  }else{
    comm_update_double_surf(packed_sed_bed, num_doubles_sed_bed, grid->smpi);
    comm_update_double(packed_sed, num_doubles_sed, grid->smpi);
  }

  /* unpack the array */

  for(j=0,ind=0;ind<grid->nnodes;ind++){
    sed->susload_vector[ind].x = packed_sed[j++];
    sed->susload_vector[ind].y = packed_sed[j++];
    for(k=0;k<ngrains;k++){
      sed->susload[k].c[ind] = packed_sed[j++];
      sed->susload[k].old_c[ind] = packed_sed[j++];
      sed->susload[k].older_c[ind] = packed_sed[j++];
      sed->susload[k].source[ind] = packed_sed[j++];
      sed->susload[k].sink[ind] = packed_sed[j++];
      sed->susload[k].error[ind] = packed_sed[j++];
      if(grid->ndim == 2){
        sed->susload[i].mfcf[ind] = packed_sed[j++];
        sed->susload[i].vcf[ind].x = packed_sed[j++];
        sed->susload[i].vcf[ind].y = packed_sed[j++];
        sed->susload[i].vor_vel[ind].x = packed_sed[j++];
        sed->susload[i].vor_vel[ind].y = packed_sed[j++];
      }
    }
  }

  for(j=0,ind=0;ind<nnodes_bed;ind++){
    sed->friction_coef[ind] = packed_sed_bed[j++];
    sed->as_ceiling[ind] = packed_sed[j++];
    sed->old_as_ceiling[ind] = packed_sed[j++];
    sed->bed_displacement[ind] = packed_sed[j++];
    sed->old_bed_displacement[ind] = packed_sed[j++];
    sed->bed_shear_stress[ind] = packed_sed[j++];
    sed->bedload_vector[ind].x = packed_sed[j++];
    sed->bedload_vector[ind].y = packed_sed[j++];
    sed->bed_gradient[ind].x = packed_sed[j++];
    sed->bed_gradient[ind].y = packed_sed[j++];
    for(k=0;k<5;k++){
      for(jj=0;jj<nconti;jj++){
        sed->consolidation_bed_properties[k][jj][ind] = packed_sed[j++];
      }
    }
    for(i=0;i<ngrains;i++){
      sed->susload[i].rouse_coef[ind] = packed_sed[j++];
      sed->bedload[i].c[ind] = packed_sed[j++];
      sed->bedload[i].old_c[ind] = packed_sed[j++];
      sed->bedload[i].older_c[ind] = packed_sed[j++];
      sed->bedload[i].thick[ind] = packed_sed[j++];
      sed->bedload[i].old_thick[ind] = packed_sed[j++];
      sed->bedload[i].older_thick[ind] = packed_sed[j++];
      sed->bedload[i].source[ind] = packed_sed[j++];
      sed->bedload[i].sink[ind] = packed_sed[j++];
      sed->bedload[i].shear[ind] = packed_sed[j++];
      sed->bedload[i].error[ind] = packed_sed[j++];
      sed->bedload[i].v[ind].x = packed_sed[j++];
      sed->bedload[i].v[ind].y = packed_sed[j++];
      sed->bedload[i].flux[ind].x = packed_sed[j++];
      sed->bedload[i].flux[ind].y = packed_sed[j++];
      sed->active_layer[0].distribution[i][ind] = packed_sed[j++];
      sed->old_active_layer[0].distribution[i][ind] = packed_sed[j++];
      for(jj=0;jj<nlayers;jj++){
        sed->bed_layer[jj].distribution[i][ind] = packed_sed[j++];
        sed->old_bed_layer[jj].distribution[i][ind] = packed_sed[j++];
      }
    }
    sed->active_layer[0].bulk_density[ind] = packed_sed[j++];
    sed->active_layer[0].porosity[ind] = packed_sed[j++];
    sed->active_layer[0].critical_erosion_shear[ind] = packed_sed[j++];
    sed->active_layer[0].erosion_rate_constant[ind] = packed_sed[j++];
    sed->active_layer[0].erosion_rate_exponent[ind] = packed_sed[j++];
    sed->active_layer[0].thickness[ind] = packed_sed[j++];
    //sed->active_layer[0].sedflume.shear_stress[ind];
    //sed->active_layer[0].sedflume.erosion_rate[ind];
    sed->old_active_layer[0].bulk_density[ind] = packed_sed[j++];
    sed->old_active_layer[0].porosity[ind] = packed_sed[j++];
    sed->old_active_layer[0].critical_erosion_shear[ind] = packed_sed[j++];
    sed->old_active_layer[0].erosion_rate_constant[ind] = packed_sed[j++];
    sed->old_active_layer[0].erosion_rate_exponent[ind] = packed_sed[j++];
    sed->old_active_layer[0].thickness[ind] = packed_sed[j++];
    //sed->old_active_layer[0].sedflume.shear_stress[ind];
    //sed->old_active_layer[0].sedflume.erosion_rate[ind];
    for (jj=0;jj<nlayers;jj++){
      sed->bed_layer[jj].bulk_density[ind] = packed_sed[j++];
      sed->bed_layer[jj].porosity[ind] = packed_sed[j++];
      sed->bed_layer[jj].critical_erosion_shear[ind] = packed_sed[j++];
      sed->bed_layer[jj].erosion_rate_constant[ind] = packed_sed[j++];
      sed->bed_layer[jj].erosion_rate_exponent[ind] = packed_sed[j++];
      sed->bed_layer[jj].thickness[ind] = packed_sed[j++];
      //sed->bed_layer[j].sedflume.shear_stress[ind];
      //sed->bed_layer[j].sedflume.erosion_rate[ind];
      sed->old_bed_layer[jj].bulk_density[ind] = packed_sed[j++];
      sed->old_bed_layer[jj].porosity[ind] = packed_sed[j++];
      sed->old_bed_layer[jj].critical_erosion_shear[ind] = packed_sed[j++];
      sed->old_bed_layer[jj].erosion_rate_constant[ind] = packed_sed[j++];
      sed->old_bed_layer[jj].erosion_rate_exponent[ind] = packed_sed[j++];
      sed->old_bed_layer[jj].thickness[ind] = packed_sed[j++];
    }
  }

  packed_sed = (double *)tl_free(sizeof(double), nnodes * num_doubles_sed, packed_sed);
  packed_sed_bed = (double *)tl_free(sizeof(double), nnodes_bed * num_doubles_sed_bed, packed_sed_bed);

  return;
}


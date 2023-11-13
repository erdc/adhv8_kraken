/* this routine unpacks the double information for the given node in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void node_unpackd_bed(void *buffer, /* the void buffer to be unpacked */
                  int ind,       /* the node to be unpacked */
                  SSW *sw
#ifdef _SEDIMENT
                   , SSED *sed,
                   int nconti,
                   int ngrains,
                   int nlayers
#endif
    )
{
  int ibuff = 0;                /* current counter in the buffer */
  double *dbuffer;              /* the double version of the buffer */
  int i,j,k;
  /* cast the buffer */
  dbuffer = (double *) buffer;


  /* sets the buffer */
	sw->d3->bed_elevation[ind] = dbuffer[ibuff++];
	sw->d3->bottom_vel[ind].x = dbuffer[ibuff++];
	sw->d3->bottom_vel[ind].y = dbuffer[ibuff++];
#ifdef _SEDIMENT
    if (sed != NULL){
    sed->friction_coef[ind] = dbuffer[ibuff++];
    sed->as_ceiling[ind] = dbuffer[ibuff++];
    sed->old_as_ceiling[ind] = dbuffer[ibuff++];
    sed->bed_displacement[ind] = dbuffer[ibuff++];
    sed->old_bed_displacement[ind] = dbuffer[ibuff++];
    sed->bed_shear_stress[ind] = dbuffer[ibuff++];
    sed->bedload_vector[ind].x = dbuffer[ibuff++];
    sed->bedload_vector[ind].y = dbuffer[ibuff++];
    sed->bed_gradient[ind].x = dbuffer[ibuff++];
    sed->bed_gradient[ind].y = dbuffer[ibuff++];
    for(k=0;k<5;k++){
      for(j=0;j<nconti;j++){
        sed->consolidation_bed_properties[k][j][ind] = dbuffer[ibuff++];
      }
    }
    for(i=0;i<ngrains;i++){
      sed->susload[i].rouse_coef[ind] = dbuffer[ibuff++];
      sed->bedload[i].c[ind] = dbuffer[ibuff++];
      sed->bedload[i].old_c[ind] = dbuffer[ibuff++];
      sed->bedload[i].older_c[ind] = dbuffer[ibuff++];
      sed->bedload[i].thick[ind] = dbuffer[ibuff++];
      sed->bedload[i].old_thick[ind] = dbuffer[ibuff++];
      sed->bedload[i].older_thick[ind] = dbuffer[ibuff++];
      sed->bedload[i].source[ind] = dbuffer[ibuff++];
      sed->bedload[i].sink[ind] = dbuffer[ibuff++];
      sed->bedload[i].shear[ind] = dbuffer[ibuff++];
      sed->bedload[i].error[ind] = dbuffer[ibuff++];
      sed->bedload[i].v[ind].x = dbuffer[ibuff++];
      sed->bedload[i].v[ind].y = dbuffer[ibuff++];
      sed->bedload[i].flux[ind].x = dbuffer[ibuff++];
      sed->bedload[i].flux[ind].y = dbuffer[ibuff++];
      sed->active_layer[0].distribution[i][ind] = dbuffer[ibuff++];
      sed->old_active_layer[0].distribution[i][ind] = dbuffer[ibuff++];
      for(j=0;j<nlayers;j++){
        sed->bed_layer[j].distribution[i][ind] = dbuffer[ibuff++];
        sed->old_bed_layer[j].distribution[i][ind] = dbuffer[ibuff++];
      }
    }
    sed->active_layer[0].bulk_density[ind] = dbuffer[ibuff++];
    sed->active_layer[0].porosity[ind] = dbuffer[ibuff++];
    sed->active_layer[0].critical_erosion_shear[ind] = dbuffer[ibuff++];
    sed->active_layer[0].erosion_rate_constant[ind] = dbuffer[ibuff++];
    sed->active_layer[0].erosion_rate_exponent[ind] = dbuffer[ibuff++];
    sed->active_layer[0].thickness[ind] = dbuffer[ibuff++];
    //sed->active_layer[0].sedflume.shear_stress[ind];
    //sed->active_layer[0].sedflume.erosion_rate[ind];
    sed->old_active_layer[0].bulk_density[ind] = dbuffer[ibuff++];
    sed->old_active_layer[0].porosity[ind] = dbuffer[ibuff++];
    sed->old_active_layer[0].critical_erosion_shear[ind] = dbuffer[ibuff++];
    sed->old_active_layer[0].erosion_rate_constant[ind] = dbuffer[ibuff++];
    sed->old_active_layer[0].erosion_rate_exponent[ind] = dbuffer[ibuff++];
    sed->old_active_layer[0].thickness[ind] = dbuffer[ibuff++];
    //sed->old_active_layer[0].sedflume.shear_stress[ind];
    //sed->old_active_layer[0].sedflume.erosion_rate[ind];
    for (j=0;j<nlayers;j++){
      sed->bed_layer[j].bulk_density[ind] = dbuffer[ibuff++];
      sed->bed_layer[j].porosity[ind] = dbuffer[ibuff++];
      sed->bed_layer[j].critical_erosion_shear[ind] = dbuffer[ibuff++];
      sed->bed_layer[j].erosion_rate_constant[ind] = dbuffer[ibuff++];
      sed->bed_layer[j].erosion_rate_exponent[ind] = dbuffer[ibuff++];
      sed->bed_layer[j].thickness[ind] = dbuffer[ibuff++];
      //sed->bed_layer[j].sedflume.shear_stress[ind];
      //sed->bed_layer[j].sedflume.erosion_rate[ind];
      sed->old_bed_layer[j].bulk_density[ind] = dbuffer[ibuff++];
      sed->old_bed_layer[j].porosity[ind] = dbuffer[ibuff++];
      sed->old_bed_layer[j].critical_erosion_shear[ind] = dbuffer[ibuff++];
      sed->old_bed_layer[j].erosion_rate_constant[ind] = dbuffer[ibuff++];
      sed->old_bed_layer[j].erosion_rate_exponent[ind] = dbuffer[ibuff++];
      sed->old_bed_layer[j].thickness[ind] = dbuffer[ibuff++];
    }
    }
#endif
}
#else
void node_unpackd_bed(void)
{
}
#endif

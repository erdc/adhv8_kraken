/* this routine packs the double information for the given node in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void node_packd(void *buffer,   /* the void buffer to be packed */
                int ind,         /* the node to be packed */
                SSW *sw,
#ifdef _ADH_GROUNDWATER
                SGW *sgw,
#endif
                SGRID *grid,
                int ntransport,
                SCON *con
#ifdef _SEDIMENT
                , SSED *sed,
                int nconti,
                int ngrains,
                int nlayers
#endif

)
{
    int itrn,i,j,k;                     /* loop counter over the transported quantities */
    int ibuff = 0;                /* current counter in the buffer */
    double *dbuffer;              /* the double version of the buffer */
    
//    printf("myid: %d ind: %d nnode: %d\n",grid->smpi->myid,ind,grid->nnodes);
//    tag(MPI_COMM_WORLD);
//    MPI_Barrier(MPI_COMM_WORLD);
//
    /* cast the buffer */
    dbuffer = (double *) buffer;
    
    /* sets the buffer */
    dbuffer[ibuff++] = grid->node[ind].x;
    dbuffer[ibuff++] = grid->node[ind].y;
    dbuffer[ibuff++] = grid->node[ind].z;
    
    if (sw != NULL) {
        if (sw->d2 != NULL) {
            
            dbuffer[ibuff++] = sw->d2->vel[ind].x;
            dbuffer[ibuff++] = sw->d2->vel[ind].y;
            dbuffer[ibuff++] = sw->d2->old_vel[ind].x;
            dbuffer[ibuff++] = sw->d2->old_vel[ind].y;
            dbuffer[ibuff++] = sw->d2->older_vel[ind].x;
            dbuffer[ibuff++] = sw->d2->older_vel[ind].y;
            dbuffer[ibuff++] = sw->d2->head[ind];
            dbuffer[ibuff++] = sw->d2->old_head[ind];
            dbuffer[ibuff++] = sw->d2->older_head[ind];
            dbuffer[ibuff++] = sw->d2->bed_elevation[ind];
            dbuffer[ibuff++] = sw->d2->bed_displacement[ind];
            
            if (sw->d2->waves != NULL) {
                dbuffer[ibuff++] = sw->d2->waves[ind].rads.xx;
                dbuffer[ibuff++] = sw->d2->waves[ind].rads.xy;
                dbuffer[ibuff++] = sw->d2->waves[ind].rads.yy;
                dbuffer[ibuff++] = sw->d2->waves[ind].stress.x;
                dbuffer[ibuff++] = sw->d2->waves[ind].stress.y;
            }
            
            if (sw->d2->winds != NULL) {
                dbuffer[ibuff++] = sw->d2->winds[ind].stress.x;
                dbuffer[ibuff++] = sw->d2->winds[ind].stress.y;
            }
#ifdef _SEDIMENT
            if (sed != NULL){
                dbuffer[ibuff++] = sed->friction_coef[ind];
                dbuffer[ibuff++] = sed->as_ceiling[ind];
                dbuffer[ibuff++] = sed->old_as_ceiling[ind];
                dbuffer[ibuff++] = sed->bed_displacement[ind];
                dbuffer[ibuff++] = sed->old_bed_displacement[ind];
                dbuffer[ibuff++] = sed->bed_shear_stress[ind];
                dbuffer[ibuff++] = sed->bedload_vector[ind].x;
                dbuffer[ibuff++] = sed->bedload_vector[ind].y;
                dbuffer[ibuff++] = sed->susload_vector[ind].x;
                dbuffer[ibuff++] = sed->susload_vector[ind].y;
                dbuffer[ibuff++] = sed->bed_gradient[ind].x;
                dbuffer[ibuff++] = sed->bed_gradient[ind].y;
                for(k=0;k<5;k++){
                    for(j=0;j<nconti;j++){
                        dbuffer[ibuff++] = sed->consolidation_bed_properties[k][j][ind];
                    }
                }
                for(i=0;i<ngrains;i++){
                    dbuffer[ibuff++] = sed->susload[i].c[ind];
                    dbuffer[ibuff++] = sed->susload[i].old_c[ind];
                    dbuffer[ibuff++] = sed->susload[i].older_c[ind];
                    dbuffer[ibuff++] = sed->susload[i].source[ind];
                    dbuffer[ibuff++] = sed->susload[i].sink[ind];
                    dbuffer[ibuff++] = sed->susload[i].error[ind];
                    dbuffer[ibuff++] = sed->susload[i].rouse_coef[ind];
                    dbuffer[ibuff++] = sed->susload[i].mfcf[ind];
                    dbuffer[ibuff++] = sed->susload[i].vcf[ind].x;
                    dbuffer[ibuff++] = sed->susload[i].vcf[ind].y;
                    dbuffer[ibuff++] = sed->susload[i].vor_vel[ind].x;
                    dbuffer[ibuff++] = sed->susload[i].vor_vel[ind].y;
                    dbuffer[ibuff++] = sed->bedload[i].c[ind];
                    dbuffer[ibuff++] = sed->bedload[i].old_c[ind];
                    dbuffer[ibuff++] = sed->bedload[i].older_c[ind];
                    dbuffer[ibuff++] = sed->bedload[i].thick[ind];
                    dbuffer[ibuff++] = sed->bedload[i].old_thick[ind];
                    dbuffer[ibuff++] = sed->bedload[i].older_thick[ind];
                    dbuffer[ibuff++] = sed->bedload[i].source[ind];
                    dbuffer[ibuff++] = sed->bedload[i].sink[ind];
                    dbuffer[ibuff++] = sed->bedload[i].shear[ind];
                    dbuffer[ibuff++] = sed->bedload[i].error[ind];
                    dbuffer[ibuff++] = sed->bedload[i].v[ind].x;
                    dbuffer[ibuff++] = sed->bedload[i].v[ind].y;
                    dbuffer[ibuff++] = sed->bedload[i].flux[ind].x;
                    dbuffer[ibuff++] = sed->bedload[i].flux[ind].y;
                    dbuffer[ibuff++] = sed->active_layer[0].distribution[i][ind];
                    dbuffer[ibuff++] = sed->old_active_layer[0].distribution[i][ind];
                    for(j=0;j<nlayers;j++){
                        dbuffer[ibuff++] = sed->bed_layer[j].distribution[i][ind];
                        dbuffer[ibuff++] = sed->old_bed_layer[j].distribution[i][ind];
                    }
                }
                dbuffer[ibuff++] = sed->active_layer[0].bulk_density[ind];
                dbuffer[ibuff++] = sed->active_layer[0].porosity[ind];
                dbuffer[ibuff++] = sed->active_layer[0].critical_erosion_shear[ind];
                dbuffer[ibuff++] = sed->active_layer[0].erosion_rate_constant[ind];
                dbuffer[ibuff++] = sed->active_layer[0].erosion_rate_exponent[ind];
                dbuffer[ibuff++] = sed->active_layer[0].thickness[ind];
                //dbuffer[ibuff++] = sed->active_layer[0].sedflume.shear_stress[ind];
                //dbuffer[ibuff++] = sed->active_layer[0].sedflume.erosion_rate[ind];
                dbuffer[ibuff++] = sed->old_active_layer[0].bulk_density[ind];
                dbuffer[ibuff++] = sed->old_active_layer[0].porosity[ind];
                dbuffer[ibuff++] = sed->old_active_layer[0].critical_erosion_shear[ind];
                dbuffer[ibuff++] = sed->old_active_layer[0].erosion_rate_constant[ind];
                dbuffer[ibuff++] = sed->old_active_layer[0].erosion_rate_exponent[ind];
                dbuffer[ibuff++] = sed->old_active_layer[0].thickness[ind];
                //dbuffer[ibuff++] = sed->old_active_layer[0].sedflume.shear_stress[ind];
                //dbuffer[ibuff++] = sed->old_active_layer[0].sedflume.erosion_rate[ind];
                for (j=0;j<nlayers;j++){
                    dbuffer[ibuff++] = sed->bed_layer[j].bulk_density[ind];
                    dbuffer[ibuff++] = sed->bed_layer[j].porosity[ind];
                    dbuffer[ibuff++] = sed->bed_layer[j].critical_erosion_shear[ind];
                    dbuffer[ibuff++] = sed->bed_layer[j].erosion_rate_constant[ind];
                    dbuffer[ibuff++] = sed->bed_layer[j].erosion_rate_exponent[ind];
                    dbuffer[ibuff++] = sed->bed_layer[j].thickness[ind];
                    //dbuffer[ibuff++] = sed->bed_layer[j].sedflume.shear_stress[ind];
                    //dbuffer[ibuff++] = sed->bed_layer[j].sedflume.erosion_rate[ind];
                    dbuffer[ibuff++] = sed->old_bed_layer[j].bulk_density[ind];
                    dbuffer[ibuff++] = sed->old_bed_layer[j].porosity[ind];
                    dbuffer[ibuff++] = sed->old_bed_layer[j].critical_erosion_shear[ind];
                    dbuffer[ibuff++] = sed->old_bed_layer[j].erosion_rate_constant[ind];
                    dbuffer[ibuff++] = sed->old_bed_layer[j].erosion_rate_exponent[ind];
                    dbuffer[ibuff++] = sed->old_bed_layer[j].thickness[ind];
                    //dbuffer[ibuff++] = sed->old_bed_layer[j].sedflume.shear_stress[ind];
                    //dbuffer[ibuff++] = sed->old_bed_layer[j].sedflume.erosion_rate[ind];
                }
            }
#endif
        }
        
        if (sw->d3 != NULL) {

            dbuffer[ibuff++] = sw->d3->vel[ind].x;
            dbuffer[ibuff++] = sw->d3->vel[ind].y;
            dbuffer[ibuff++] = sw->d3->vel[ind].z;
            
            dbuffer[ibuff++] = sw->d3->old_vel[ind].x;
            dbuffer[ibuff++] = sw->d3->old_vel[ind].y;
            dbuffer[ibuff++] = sw->d3->old_vel[ind].z;
            dbuffer[ibuff++] = sw->d3->older_vel[ind].x;
            dbuffer[ibuff++] = sw->d3->older_vel[ind].y;
            dbuffer[ibuff++] = sw->d3->older_vel[ind].z;
            
            dbuffer[ibuff++] = sw->d3->prs[ind];
            dbuffer[ibuff++] = sw->d3->prs_plus[ind];
            dbuffer[ibuff++] = sw->d3->prs_minus[ind];
            dbuffer[ibuff++] = sw->d3->displacement[ind];
            dbuffer[ibuff++] = sw->d3->old_displacement[ind];
            dbuffer[ibuff++] = sw->d3->older_displacement[ind];
            dbuffer[ibuff++] = sw->d3->grid_speed[ind];
            dbuffer[ibuff++] = sw->d3->old_grid_speed[ind];
            // dbuffer[ibuff++] = sw->d3->elev_factor[ind];
            dbuffer[ibuff++] = sw->d3->density[ind];
            dbuffer[ibuff++] = sw->d3->vertical_node_flux[ind];
            dbuffer[ibuff++] = sw->d3->tanvec[ind].x;
            dbuffer[ibuff++] = sw->d3->tanvec[ind].y;
            dbuffer[ibuff++] = sw->d3->error[ind];
            dbuffer[ibuff++] = sw->d3->hyd_viscosity[ind];
            dbuffer[ibuff++] = sw->d3->trn_diffusivity[ind];
            dbuffer[ibuff++] = sw->d3->dpl_perturbation[ind];
#ifdef _SEDIMENT
            dbuffer[ibuff++] = sw->d3->bed_displacement[ind];
            dbuffer[ibuff++] = sw->d3->old_bed_displacement[ind];
            dbuffer[ibuff++] = sw->d3->older_bed_displacement[ind];
            if(sed!=NULL){
                dbuffer[ibuff++] = sed->susload_vector[ind].x;
                dbuffer[ibuff++] = sed->susload_vector[ind].y;
                for(i=0;i<ngrains;i++){
                    dbuffer[ibuff++] = sed->susload[i].c[ind];
                    dbuffer[ibuff++] = sed->susload[i].old_c[ind];
                    dbuffer[ibuff++] = sed->susload[i].older_c[ind];
                    dbuffer[ibuff++] = sed->susload[i].source[ind];
                    dbuffer[ibuff++] = sed->susload[i].sink[ind];
                    dbuffer[ibuff++] = sed->susload[i].error[ind];
                }
            }
#endif
            
        }
    }
    
#ifdef _ADH_GROUNDWATER
    if (sgw != NULL) {
        dbuffer[ibuff++] = sgw->gw_phead[ind];
        dbuffer[ibuff++] = sgw->old_gw_phead[ind];
        dbuffer[ibuff++] = sgw->older_gw_phead[ind];
        dbuffer[ibuff++] = sgw->predict_gw_phead[ind];
        dbuffer[ibuff++] = sgw->gw_density[ind];
        dbuffer[ibuff++] = sgw->error[ind];
    }
#endif

    if (con != NULL) {
        for (itrn = 0; itrn < ntransport; itrn++) {
            //dbuffer[ibuff++] = con[itrn].property[ind];
            dbuffer[ibuff++] = con[itrn].concentration[ind];
            dbuffer[ibuff++] = con[itrn].old_concentration[ind];
            dbuffer[ibuff++] = con[itrn].older_concentration[ind];
            dbuffer[ibuff++] = con[itrn].sink[ind];
            dbuffer[ibuff++] = con[itrn].source[ind];
            dbuffer[ibuff++] = con[itrn].nodal_decay_coefficient[ind];
            dbuffer[ibuff++] = con[itrn].error[ind];
            dbuffer[ibuff++] = con[itrn].mfcf[ind];
            dbuffer[ibuff++] = con[itrn].vcf[ind].x;
            dbuffer[ibuff++] = con[itrn].vcf[ind].y;
        }
    }
}
#else
void node_packd(void)
{
}
#endif

/* this routine reads the sed library arrays and sets up the link to the sed library */

#include "global_header.h"
#include "sedlib_header.h"

static int DEBUG = OFF;


void sedlib_link_bed(SMODEL *mod, int inode, int ndim, int depth_only, SVECT2D bed_gradient) {
    
    int i, j, k;	/* loop counters */
 
    // aliases
    SSED *sed = mod->sed;
    SSEDLIB *sedlib = mod->sed->sedlib; 
    SSW_2D *sw2d = mod->sw->d2;
    SSW_3D *sw3d = mod->sw->d3;
    int n_sed = mod->sed->nsed;
    int n_blay = mod->sed->nlayers;
//   if(sed->old_bed_displacement[inode]>=0)printf("MYID %d has dips lessthan 0 at node %d \n", mod->grid->smpi->myid, inode); 
    // initialize sedlib for this node
    ssedlib_init(sedlib);

    /** Sedlib variable = ADH variable **/
    /* Ints */
    sedlib->consolidate_flag = 0;
    if (sed->ncti > 0) sedlib->consolidate_flag = 1;
    sedlib->coh_settling_flag = sed->cohesive_settling_flag;
    sedlib->wave_current_flag = sed->wind_wave_flag;
    sedlib->noncoh_sedent_flag = sed->noncoh_ent_flag;
    sedlib->noncoh_bedload_flag = sed->noncoh_ble_flag;
    sedlib->hiding_factor_flag = sed->hid_fact_flag;
    
    // find nodal roughness type
    sedlib->roughness_type = sed->friction_flag[inode];
    for (i=0; i<n_sed; i++) {
        sedlib->grain[i].type = UNSET_INT;
        if (sed->grain[i].type == SND) {
            sedlib->grain[i].type = 0;
        } else if (sed->grain[i].type == CLA) {
            sedlib->grain[i].type = 1;
        }
    }

    for (i=0; i<n_blay; i++) {
        sedlib->layer[i].sedflume_flag = sed->bed_layer[i].sedflume_flag;
    }
   
    /* Doubles */
    sedlib->al_por = sed->old_active_layer->porosity[inode];
    sedlib->al_ces = sed->old_active_layer->critical_erosion_shear[inode];
    sedlib->al_erc = sed->old_active_layer->erosion_rate_constant[inode];
    sedlib->al_ere = sed->old_active_layer->erosion_rate_exponent[inode];
    sedlib->old_al_por = sed->old_active_layer->porosity[inode];
    sedlib->old_al_ces = sed->old_active_layer->critical_erosion_shear[inode];
    sedlib->old_al_erc = sed->old_active_layer->erosion_rate_constant[inode];
    sedlib->old_al_ere = sed->old_active_layer->erosion_rate_exponent[inode];
    sedlib->as_ceiling = sed->old_as_ceiling[inode];
    sedlib->old_as_ceiling = sed->old_as_ceiling[inode];
    sedlib->disp = sed->old_bed_displacement[inode];
    sedlib->old_disp =  sed->old_bed_displacement[inode];
    sedlib->al_thick = sedlib->disp - sedlib->as_ceiling;

    int id3d = UNSET_INT;
    if (ndim == 2) {
        id3d = inode; // no mapping required in 2d
        sedlib->depth_old = floorf(sw2d->head[inode]/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        sedlib->depth_new = sedlib->depth_old;
        sedlib->depth_old = MAX(sw2d->head[inode],.0001);
        sedlib->depth_new = MAX(sedlib->depth_new,mod->drying_lower_limit);
        sedlib->v_da.x = floorf(sw2d->vel[inode].x/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        sedlib->v_da.y = floorf(sw2d->vel[inode].y/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        sedlib->v_nb.x = floorf(sw2d->vel[inode].x/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        sedlib->v_nb.y = floorf(sw2d->vel[inode].y/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        
        SVECT2D vort; svect2d_init(&(vort));
        if (mod->flag.VORTICITY) {
            vort = sed_vorticity_velocity_components(sedlib->depth_new, sw2d->vel[inode], mod->con[mod->vorticity_id], inode);
            sedlib->v_vort.x = floorf(vort.x/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
            sedlib->v_vort.y = floorf(vort.y/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        }
        
        sedlib->depth_averaged_flag = ON;

    } else if (ndim == 3) {
        id3d = mod->grid->nodeID_2d_to_3d_bed[inode];
        //printf("node2d to 3d: %d\n",mod->grid->nodeID_2d_to_3d_bed[inode]);
        sedlib->depth_old = floorf((sw3d->prs[id3d] * mod->density / sw3d->density[id3d])/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        sedlib->depth_new = sedlib->depth_old;
        sedlib->v_da.x = floorf(sw3d->vel[id3d].x/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        sedlib->v_da.y = floorf(sw3d->vel[id3d].y/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        sedlib->v_nb.x = sedlib->v_da.x;
        sedlib->v_nb.y = sedlib->v_da.y;
        sedlib->v_vort.x = 0;
        sedlib->v_vort.y = 0;
        
        sedlib->depth_averaged_flag =  OFF;
    }
    
    //// winds (does this take stress or velocity?
    //v_wind.x = 0.;
    //v_wind.y = 0.;
    //if (mod->flag.WIND) {
    //  v_wind.x = winds[n].x;
    //  v_wind.y = winds[n].y;
    //}
    //// attenuate winds
    //v_wind.x *= wind_attn;
    //v_wind.y *= wind_attn;
    
    
    sedlib->bed_slope.x = floorf(bed_gradient.x/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
    sedlib->bed_slope.y = floorf(bed_gradient.y/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;

    if (depth_only) {
        sedlib->v_da.x = 0.0;
        sedlib->v_da.y = 0.0;
        sedlib->v_vort.x = 0.0;
        sedlib->v_vort.y = 0.0;
        sedlib->v_wind.x = 0.0;
        sedlib->v_wind.y = 0.0;
        sedlib->bed_slope.x = 0.0;
        sedlib->bed_slope.y = 0.0;
    }

    sedlib->kin_visc = mod->viscosity;
    sedlib->grav = mod->gravity;
    sedlib->rho = mod->density;
    
    sedlib->roughness_coefficient = sed->friction_coef[inode];;

    sedlib->ice_roughness = 0.0;
    //if (irh > 1.E-8) ice_roughness = irh;
    
    sedlib->time_step = mod->dt;
    
    /* (cjt) */
    sedlib->wave.height = 0.0;
    sedlib->wave.period = 0.0;
    sedlib->wave.angle  = 0.0;
    sedlib->wave.number = 0.0;
    sedlib->wave.ibreak = 0.0;
    sedlib->wave.speed  = 0.0;
    sedlib->wave.edb    = 0.0;
    sedlib->wave.edf    = 0.0;
    sedlib->wave.seb    = 0.0;
    sedlib->wave.sef    = 0.0;
    sedlib->wave.wcf    = 0.0;
    sedlib->wave.beta   = 0.0;
    sedlib->wave.gama   = 0.0;
    
    /* ADH flag for if short waves are used */
    //      WAVE_SED = WAVE;
    //      if (WAVE_SED) {   /* (cjt) */
    //
    //        if (noncoh_ent_flag == 3) {
    //                wave_seb  = noncoh_sed_ent_const[0];
    //                wave_sef  = noncoh_sed_ent_const[1];
    //        }
    //
    //        if (noncoh_ble_flag == 3) {
    //                wave_wcf  = noncoh_sed_ble_const[0];
    //                wave_beta = noncoh_sed_ble_const[1];
    //                wave_gama = noncoh_sed_ble_const[2];
    //        }
    //
    //		/* Hmo is sent in, however, Hrms is used by the processes, so convert */
    //        wave_height = waves[n].height / sqrt(2) ;
    //        wave_period = waves[n].period;
    //        wave_angle  = waves[n].angle;
    //        wave_number = waves[n].number;
    //        wave_ibreak = waves[n].breaking;
    //        wave_speed  = waves[n].speed;
    //        wave_edb    = waves[n].ediss_break;
    //        wave_edf    = waves[n].ediss_frict;
    //      }
    
    /** read the arrays into the sedlib arrays **/
    for(i = 0; i < n_blay; i++) {
        sedlib->layer[i].ces = sed->old_bed_layer[i].critical_erosion_shear[inode];
        sedlib->layer[i].erc = sed->old_bed_layer[i].erosion_rate_constant[inode];
        sedlib->layer[i].ere = sed->old_bed_layer[i].erosion_rate_exponent[inode];
        sedlib->layer[i].thick = sed->old_bed_layer[i].thickness[inode];
    }
    
    for(i = 0; i < n_sed; i++) {
        sedlib->grain[i].al_dist = sed->old_active_layer->distribution[i][inode];
        sedlib->grain[i].bl_con = floorf(sed->bedload[i].old_c[inode]/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        sedlib->grain[i].bl_con *= sed->grain[i].reference_c;
        sedlib->grain[i].bl_thick = sed->bedload[i].thick[inode];
        sedlib->grain[i].sl_con = floorf(sed->susload[i].old_c[id3d]/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        sedlib->grain[i].sl_con *= sed->grain[i].reference_c;
        sedlib->grain[i].sl_bcon = floorf(sed->susload[i].old_c[id3d]/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
        sedlib->grain[i].sl_bcon *= sed->grain[i].reference_c;
        sedlib->grain[i].sg = sed->grain[i].specific_gravity;
        sedlib->grain[i].diameter = sed->grain[i].diameter;
        sedlib->grain[i].por = sed->grain[i].porosity;
        sedlib->grain[i].old_al_dist = sed->old_active_layer->distribution[i][inode];
        if (sed->grain[i].type == CLA) {
            sedlib->grain[i].coh_free_settling_vel = sed->grain[i].clay.settling_velocity;
            sedlib->grain[i].ces = sed->grain[i].clay.tau_ce;
            sedlib->grain[i].cds = sed->grain[i].clay.tau_cd;
            sedlib->grain[i].erc = sed->grain[i].clay.erode_const;
            sedlib->grain[i].ere = 1.0;
            for(j = 0; j < 4; j++) {
                sedlib->coh_settling_constants[j] = sed->cohesive_settling_const[j];
            }
        }
    }
    sedlib->grain[n_sed].al_dist = 0.0;
    sedlib->grain[n_sed].bl_con = 0.0;
    sedlib->grain[n_sed].bl_thick = 1.0;
    sedlib->grain[n_sed].sl_con = 0.0;
    sedlib->grain[n_sed].sl_bcon = 0.0;
    sedlib->grain[n_sed].sg = sed->grain[n_sed-1].specific_gravity;
    sedlib->grain[n_sed].old_al_dist = 0.0;
    
    for(i = 0; i < n_blay; i++){
        for(j = 0; j < n_sed; j++) {
            sedlib->bl_dist[i][j] = sed->old_bed_layer[i].distribution[j][inode];
        }
    }
    /* for now
    if (n_sedfl > 0) {
        for(i = 0; i < n_blay; i++) {
            for(j = 0; j < n_sedfl; j++) {
                sedlib->sedfl_shear[i][j] = sed->bed_layer[i].sedflume[j].shear_stress[inode];
                sedlib->sedfl_erate[i][j] = sed->bed_layer[i].sedflume[j].erosion_rate[inode];
            }
        }
    }
    */
 

    if (sed->ncti > 0) {
        for(i = 0; i < sed->ncti; i++){
            for(j = 0; j < 5; j++)
                sedlib->consolidate_table[i][j] = sed->consolidation_bed_properties[j][i][inode];
        }
    }

    double temp;
    for(i = 0; i < n_blay; i++) {
        temp = 0.0;
        for(j = 0; j < n_sed; j++) {
            temp += sedlib->grain[j].sg * sedlib->bl_dist[i][j];
        }
        sedlib->layer[i].por = (mod->density * temp - sed->old_bed_layer[i].bulk_density[inode]) / (mod->density * temp - mod->density);
    }
    
    
    if (sedlib->depth_old <= 0.0001) {
        for(i = 0; i < n_sed; i++)
            sedlib->grain[i].sl_con = 0.0;
    }
    
    //*****************************************************************************************//
    //*****************************************************************************************//
    //*****************************************************************************************//
    if (DEBUG && id3d==10) {
        ssedlib_printScreen(*sedlib);
        exit(-1);
    }
    /** CALL SEDLIB **/
    sedlib_core_bed(sedlib);

    //*****************************************************************************************//
    //*****************************************************************************************//
    //*****************************************************************************************//
    
    /** Redefine ADH variables to match SED values */
    /** ADH variable = Sedlib variable**/
    /** Only necessary for values that the Sedlib recalculates **/
    int inode_prob, surf_id;
    if(solv_isnan(sedlib->disp) || solv_isinf(sedlib->disp)) {
      if(mod->grid->ndim==3){
        inode_prob = mod->grid->node[mod->grid->nodeID_2d_to_3d_bed[inode]].gid;
        surf_id = mod->grid->node[mod->grid->nodeID_2d_to_3d_bed[inode]].global_bed_id;
      }else{
        inode_prob = mod->grid->node[inode].gid;
        surf_id = mod->grid->node[inode].global_bed_id;
      }
        printf("SEDLIB BED FAIL MYID %d AT NODE %d GID %d GLOBAL_SURF %d :IGNORING SEDLIB RESULTS FOR THIS TIME STEP \n", mod->grid->smpi->myid, inode, inode_prob,surf_id);
        printf("displacement: %20.10f \n",sedlib->disp);
        //exit(-1);
    } else {
        /* Ints */
        for(i=0;i<n_blay;i++) {
            sed->bed_layer[i].sedflume_flag = sedlib->layer[i].sedflume_flag;
        }
            
        /* Doubles */
        sed->active_layer->porosity[inode] = sedlib->al_por;
        sed->active_layer->critical_erosion_shear[inode] = sedlib->al_ces;
        sed->active_layer->erosion_rate_constant[inode] = sedlib->al_erc;
        sed->active_layer->erosion_rate_exponent[inode] = sedlib->al_ere;
        sed->active_layer->thickness[inode] = sedlib->al_thick;
        sed->as_ceiling[inode] = sedlib->disp - sedlib->al_thick;
        sed->bed_displacement[inode] = sedlib->disp;
        sed->bed_shear_stress[inode] = sedlib->grain_shear_stress_mag;
        
        /** Read these arrays into the ADH variables **/
        
        for(i = 0; i < n_blay; i++) {
            temp = 0.0;
            for(j = 0; j < n_sed; j++) {
                temp += sed->grain[j].specific_gravity * sedlib->bl_dist[i][j];
            }
            sed->bed_layer[i].bulk_density[inode] = mod->density * (temp * (1. - sedlib->layer[i].por) + sedlib->layer[i].por);
        }
        
        for(i = 0; i < n_blay; i++) {
            sed->bed_layer[i].critical_erosion_shear[inode] = sedlib->layer[i].ces;
            sed->bed_layer[i].erosion_rate_constant[inode] = sedlib->layer[i].erc;
            sed->bed_layer[i].erosion_rate_exponent[inode] = sedlib->layer[i].ere;
            sed->bed_layer[i].thickness[inode] = sedlib->layer[i].thick;
            for(j = 0; j < n_sed; j++) {
                sed->bed_layer[i].distribution[j][inode] = MAX(sedlib->bl_dist[i][j],0.0);
            }
        }
        
        for(i = 0; i < n_sed; i++) {
            sed->active_layer->distribution[i][inode] = MAX(sedlib->grain[i].al_dist,0.0);
            sed->bedload[i].thick[inode] = sedlib->grain[i].bl_thick;
            sed->susload[i].rouse_coef[inode] = sedlib->grain[i].sl_nbf;
            
            sedlib->grain[i].bl_net_flux = floorf(sedlib->grain[i].bl_net_flux/SMALL)*SMALL;
            sedlib->grain[i].sl_net_flux = floorf(sedlib->grain[i].sl_net_flux/SMALL)*SMALL;
            
            sed->bedload[i].source[inode] = - sedlib->grain[i].bl_net_flux / sed->grain[i].reference_c;
            sed->susload[i].source[id3d] = - sedlib->grain[i].sl_net_flux / sed->grain[i].reference_c;
        }
    }
}



/* this routine reads the sed library arrays and sets up the link to the sed library */

#include "global_header.h"
#include "sedlib_header.h"

void sedlib_link_transport(SMODEL *mod, int inode, int ndim, int depth_only, SVECT2D bed_gradient) {
    
    int i, j, k;	/* loop counters */
    
    // aliases
    SSED *sed = mod->sed;
    SSEDLIB *sedlib = mod->sed->sedlib;
    SSW_2D *sw2d = mod->sw->d2;
    SSW_3D *sw3d = mod->sw->d3;
    int n_sed = mod->sed->nsed;
    int n_blay = mod->sed->nlayers;
    
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
    
    sedlib->bed_slope.x = floorf(bed_gradient.x/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
    sedlib->bed_slope.y = floorf(bed_gradient.y/NOT_QUITE_SMALL)*NOT_QUITE_SMALL;
    
    if (depth_only == 1) {
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
 
    //      /* ADH flag for if short waves are used */
    //      WAVE_SED = WAVE;
    //
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
    //        /* Hmo is sent in, however, Hrms is used by the processes, so convert */
    //        wave_height = waves[n].height / sqrt(2) ;
    //
    //        wave_period = waves[n].period;
    //        wave_angle  = waves[n].angle;
    //        wave_number = waves[n].number;
    //        wave_ibreak = waves[n].breaking;
    //        wave_speed  = waves[n].speed;
    //        wave_edb    = waves[n].ediss_break;
    //        wave_edf    = waves[n].ediss_frict;
    //      }
    
    
    /** read the arrays into the sedlib arrays **/
    for(i = 0; i < n_sed; i++) {
        sedlib->grain[i].al_dist = sed->old_active_layer->distribution[i][inode];
        sedlib->grain[i].old_al_dist = sed->old_active_layer->distribution[i][inode];
        sedlib->grain[i].sg = sed->grain[i].specific_gravity;
        sedlib->grain[i].diameter = sed->grain[i].diameter;
        sedlib->grain[i].sl_nbf = MAX(sed->susload[i].rouse_coef[inode],1.0);
    }
    sedlib->grain[n_sed].al_dist = 0.0;
    sedlib->grain[n_sed].old_al_dist = 0.0;
    sedlib->grain[n_sed].sg = sed->grain[n_sed-1].specific_gravity;    

    for(i = 0; i < n_blay; i++){
        for(j = 0; j < n_sed; j++) {
            sedlib->bl_dist[i][j] = sed->bed_layer[i].distribution[j][inode];
        }
    }
    
    
    /** Call Sedlib **/
    for(i = 0; i < n_sed; i++) {
        
        //******************************************************************************************//
        //******************************************************************************************//
        sedlib_core_bedload_transport(sedlib, i);
        //******************************************************************************************//
        //******************************************************************************************//
        
        sed->bedload[i].thick[inode] = sedlib->grain[i].bl_thick;
        sed->bedload[i].v[inode].x = sedlib->bl_vel.x;
        sed->bedload[i].v[inode].y = sedlib->bl_vel.y;
        
        // if depth averaged shallow water is on, get transport correction factors
        if (ndim == 2 ) {
            sedlib->sl_mfcf = 1.0;
            sedlib->sl_vcf.x = 0.0;
            sedlib->sl_vcf.y = 0.0;
            
            //******************************************************************************************//
            sedlib_core_suspended_transport(sedlib, i);
            //******************************************************************************************//
            sed->susload[i].mfcf[inode] = sedlib->sl_mfcf;
            sed->susload[i].vcf[inode].x = sedlib->sl_vcf.x;
            sed->susload[i].vcf[inode].y = sedlib->sl_vcf.y;
            if(solv_isnan(sedlib->sl_mfcf) || solv_isinf(sedlib->sl_mfcf)) {
                printf("SEDLIB TRANSPORT FAIL AT PROCESSOR %d NODE %d RESIDENT PE %d :IGNORING SEDLIB RESULTS FOR THIS TIME STEP \n", mod->grid->smpi->myid, mod->grid->node[inode].gid, mod->grid->node[inode].resident_pe);
            }
        }
        
    }
    
}



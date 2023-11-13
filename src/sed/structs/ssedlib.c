#include "sedlib_header.h"

//********************************************************************************//
//********************************************************************************//
void ssedlib_wave_init(SSEDLIB_WAVE *wave) {
    wave->ibreak = 0.;
    wave->height = 0.;
    wave->period = 0.;
    wave->angle = 0.;
    wave->number = 0.;
    wave->speed = 0.;
    wave->edb = 0.;
    wave->edf = 0.;
    wave->sef = 0.;
    wave->seb = 0.;
    wave->wcf = 0.;
    wave->beta = 0.;
    wave->gama = 0.;
} 

//********************************************************************************//
//********************************************************************************//
void ssedlib_bed_layer_init(SSEDLIB_BED_LAYER *layer, int nlayers) {
    int ilayer = 0;
    for (ilayer=0; ilayer<nlayers; ilayer++) {
        layer[ilayer].type = 0;
        layer[ilayer].sedflume_flag = 0;
        layer[ilayer].fully_eroded_layer = 0;
        layer[ilayer].por = 0.;
        layer[ilayer].ces = 0.;
        layer[ilayer].erc = 0.;
        layer[ilayer].ere = 0.;
        layer[ilayer].thick = 0.;
    }
}

//********************************************************************************//
//********************************************************************************//
void ssedlib_grain_init(SSEDLIB_GRAIN *grain, int ngrains) {
    int igrain = 0;
    for (igrain=0; igrain<ngrains; igrain++) {
        grain[igrain].type = UNSET_INT;
        grain[igrain].bed_mass_change = 0.;
        grain[igrain].net_mass_change = 0.;
        grain[igrain].diameter = 0.;
        grain[igrain].settling_vel = 0.;
        grain[igrain].coh_free_settling_vel = 0.;
        grain[igrain].por = 0.;
        grain[igrain].ces = 0.;
        grain[igrain].cds = 0.;
        grain[igrain].erc = 0.;
        grain[igrain].ere = 0.;
        grain[igrain].critical_shear_stress = 0.;
        grain[igrain].critical_shear_stress_plane_bed = 0.;
        grain[igrain].sg = 0.;
        grain[igrain].al_dist = 0.;
        grain[igrain].as_dist = 0.;
        grain[igrain].bl_con = 0.;
        grain[igrain].bl_thick = 0.;
        grain[igrain].sl_con = 0.;
        grain[igrain].sl_bcon = 0.;
        grain[igrain].sl_sg = 0.;
        grain[igrain].bl_erosion_flux = 0.;
        grain[igrain].bl_deposition_flux = 0.;
        grain[igrain].sl_erosion_flux = 0.;
        grain[igrain].sl_deposition_flux = 0.;
        grain[igrain].sl_erosion_flux_nbl = 0.;
        grain[igrain].bl_erosion_flux_nbl = 0.;
        grain[igrain].bl_deposition_flux_nbl = 0.;
        grain[igrain].bl_net_flux = 0.;
        grain[igrain].sl_net_flux = 0.;
        grain[igrain].sl_nbf = 0.;
        grain[igrain].old_al_dist = 0.;
    }
}

//********************************************************************************//
//********************************************************************************//

void ssedlib_alloc(SSEDLIB **sedlib, int n_blay, int n_sed, int n_sedfl, int n_conti) {
 
    int i=0, j=0;

    (*sedlib) = (SSEDLIB *) tl_alloc(sizeof(SSEDLIB), 1);
    SSEDLIB *sed = *sedlib;
    
    sed->n_blay = n_blay;
    sed->n_sed = n_sed;
    sed->n_sp1 = n_sed + 1;
    sed->n_sedfl = 0; //n_sedfl;
    sed->n_conti = n_conti;
    
    sed->layer = (SSEDLIB_BED_LAYER *) tl_alloc(sizeof(SSEDLIB_BED_LAYER), sed->n_blay);
    sed->grain = (SSEDLIB_GRAIN *) tl_alloc(sizeof(SSEDLIB_GRAIN), sed->n_sp1);
    
    // ?? noncoh_sedent_constants  = (double *) tl_alloc(sizeof(double), 1);
    sed->coh_settling_constants = (double *) tl_alloc(sizeof(double), 4);
    
    sed->al_coh_prop = (double *)  tl_alloc(sizeof(double), 5);
    
    if (sed->n_sedfl > 0) {
        sed->sedfl_shear = (double **) tl_alloc(sizeof(double *), sed->n_blay);
        for (i = 0; i < sed->n_blay; i++) {
            sed->sedfl_shear[i] = (double *) tl_alloc(sizeof(double), sed->n_sedfl);
        }
        
        sed->sedfl_erate = (double **) tl_alloc(sizeof(double *), sed->n_blay);
        for (i = 0; i < sed->n_blay; i++) {
            sed->sedfl_erate[i] = (double *) tl_alloc(sizeof(double), sed->n_sedfl);
        }
    }
    
    if (sed->n_conti > 0) {
        sed->consolidate_table = (double **) tl_alloc(sizeof(double *), sed->n_conti);
        for (i = 0; i < sed->n_conti; i++) {
            sed->consolidate_table[i] = (double *) tl_alloc (sizeof(double), 5);
        }
    }
    
    sed->bl_dist = (double **) tl_alloc(sizeof(double *), sed->n_blay);
    for (i = 0; i < sed->n_blay; i++) {
        sed->bl_dist[i] =(double *) tl_alloc(sizeof(double), sed->n_sed);
    }
    
    ssedlib_init(sed);
    
}

//********************************************************************************//
//********************************************************************************//

void ssedlib_init(SSEDLIB *sed) {
    
        int i=0, j=0;

    sed->WAVE_SED = OFF;
    sed->roughness_type = UNSET_INT;
    sed->oald_resolve_flag = 0;
    sed->layer_being_eroded = 0;
    sed->consolidate_flag = UNSET_INT;
    sed->coh_settling_flag = UNSET_INT;
    sed->wave_current_flag = UNSET_INT;
    sed->noncoh_sedent_flag = UNSET_INT;
    sed->noncoh_bedload_flag = UNSET_INT;
    sed->hiding_factor_flag = UNSET_INT;
    sed->depth_averaged_flag = UNSET_INT;
    
    sed->al_por = 0.; 
    sed->old_al_por = 0.; 
    sed->as_por = 0.; 
    sed->al_ces = 0.;
    sed->old_al_ces = 0.;
    sed->as_ces = 0.;
    sed->al_erc = 0.;
    sed->old_al_erc = 0.;
    sed->as_erc = 0.;
    sed->al_ere = 0.;
    sed->old_al_ere = 0.;
    sed->as_ere = 0.;
    sed->al_thick = 0.;
    sed->as_ceiling = 0.;
    sed->old_as_ceiling = 0.;
    sed->disp = 0.;
    sed->old_disp = 0.;
    sed->depth_old = 0.;
    sed->depth_new = 0.;
    sed->kin_visc = 0.;
    sed->grav = 0.;
    sed->rho = 0.;
    sed->ice_roughness = 0.;
    sed->roughness_coefficient = 0.;
    sed->total_erodable_thickness = 0.;
    sed->minimum_erodable_thickness = 0.;
    sed->strata_dep_coef = 0.;
    sed->total_shear_stress_mag = 0.;
    sed->grain_shear_stress_mag = 0.;
    sed->near_bed_velocity_mag = 0.;
    sed->time_step = 0.;
    sed->ero_limit_ts = 0.;
    sed->sl_mfcf = 0.;
    sed->esrh = 0.;
    sed->tsrh = 0.;
    
    svect2d_init(&(sed->grain_shear_stress));
    svect2d_init(&(sed->v_da));
    svect2d_init(&(sed->v_nb));
    svect2d_init(&(sed->v_wind));
    svect2d_init(&(sed->v_vort));
    svect2d_init(&(sed->bl_vel));
    svect2d_init(&(sed->bed_slope));
    svect2d_init(&(sed->sl_vcf));
    

    for (i=0; i<4; i++) {
        sed->coh_settling_constants[i] = 0.;
    }

    for (i=0; i<5; i++) {
        sed->al_coh_prop[i] = 0.;
    }

    for (i = 0; i < sed->n_blay; i++) {
          for (j = 0; j < sed->n_sed; j++) {
              sed->bl_dist[i][j] = 0.;
          }
    }


    ssedlib_bed_layer_init(sed->layer, sed->n_blay);
    ssedlib_grain_init(sed->grain, (sed->n_sed)+1);
    ssedlib_wave_init(&(sed->wave));
    
}

//**********************************************************************************//
//**********************************************************************************//
//**********************************************************************************//

void ssedlib_free(SSEDLIB *sed) {

        int i;

        sed->coh_settling_constants = (double *) tl_free(sizeof(double), 4, sed->coh_settling_constants);
        sed->al_coh_prop = (double *)  tl_free(sizeof(double), 5, sed->al_coh_prop);

        if (sed->n_sedfl > 0) {
            for (i = 0; i < sed->n_blay; i++) {
                sed->sedfl_shear[i] = (double *) tl_free(sizeof(double), sed->n_sedfl, sed->sedfl_shear[i]);
                sed->sedfl_erate[i] = (double *) tl_free(sizeof(double), sed->n_sedfl, sed->sedfl_erate[i]);
            }
            sed->sedfl_shear = (double **) tl_free(sizeof(double *), sed->n_blay, sed->sedfl_shear);
            sed->sedfl_erate = (double **) tl_free(sizeof(double *), sed->n_blay, sed->sedfl_erate);
        }

        for (i = 0; i < sed->n_blay; i++) {
            sed->bl_dist[i] =(double *) tl_free(sizeof(double), sed->n_sed, sed->bl_dist[i]);
        }
        sed->bl_dist = (double **) tl_free(sizeof(double *), sed->n_blay, sed->bl_dist);

        if (sed->n_conti > 0) {
        for (i = 0; i < sed->n_conti; i++) {
            sed->consolidate_table[i] = (double *) tl_free(sizeof(double), 5, sed->consolidate_table[i]);
        }
        sed->consolidate_table = (double **) tl_free(sizeof(double *), sed->n_conti, sed->consolidate_table);
        }

        sed->layer = (SSEDLIB_BED_LAYER *) tl_free(sizeof(SSEDLIB_BED_LAYER), sed->n_blay, sed->layer);
        sed->grain = (SSEDLIB_GRAIN *) tl_free(sizeof(SSEDLIB_GRAIN), sed->n_sp1, sed->grain);
        
        sed = (SSEDLIB *) tl_free(sizeof(SSEDLIB), 1, sed);
}

//**********************************************************************************//
//**********************************************************************************//
//**********************************************************************************//

void ssedlib_printScreen(SSEDLIB sed) {
    
    int i=0, j=0;

    printf("\n");
    printf("WAVE_SED: %d\n",sed.WAVE_SED);
    printf("n_blay: %d\n",sed.n_blay);
    printf("n_sed: %d\n",sed.n_sed);
    printf("n_sp1: %d\n",sed.n_sp1);
    printf("n_sedfl: %d\n",sed.n_sedfl);
    printf("n_conti: %d\n",sed.n_conti);
    
    printf("roughness_type: %d\n",sed.roughness_type);
    printf("oald_resolve_flag: %d\n",sed.oald_resolve_flag);
    printf("layer_being_eroded: %d\n",sed.layer_being_eroded);
    printf("consolidate_flag: %d\n",sed.consolidate_flag);
    printf("coh_settling_flag: %d\n",sed.coh_settling_flag);
    printf("wave_current_flag: %d\n",sed.wave_current_flag);
    printf("noncoh_sedent_flag: %d\n",sed.noncoh_sedent_flag);
    printf("noncoh_bedload_flag: %d\n",sed.noncoh_bedload_flag);
    printf("hiding_factor_flag: %d\n",sed.hiding_factor_flag);
    printf("depth_averaged_flag: %d\n",sed.depth_averaged_flag);

    printf("al_por: %20.10f \t\t ol_al_por: %20.10f \n",sed.al_por, sed.old_al_por);
    printf("al_ces: %20.10f \t\t ol_al_ces: %20.10f \n",sed.al_ces, sed.old_al_ces);
    printf("al_erc: %20.10f \t\t ol_al_erc: %20.10f \n",sed.al_erc, sed.old_al_erc);
    printf("al_ere: %20.10f \t\t ol_al_ere: %20.10f \n",sed.al_ere, sed.old_al_ere);
    printf("al_thick: %20.10f \n",sed.al_thick);
    
    printf("as_por: %20.10f \t\t as_ces: %20.10f \t\t as_erc: %20.10f \t\t as_ere: %20.10f \t\t as_ceiling: %20.10f \t\t old_as_ceiling: %20.10f \n",sed.as_por, sed.as_ces, sed.as_erc, sed.as_ere, sed.as_ceiling, sed.old_as_ceiling);
    
    printf("disp: %20.10f \t\t old_disp: %20.10f\n",sed.disp, sed.old_disp);
    printf("depth_old: %20.10f \t\t depth_new: %20.10f \n",sed.depth_old, sed.depth_new);
    printf("kin visc: %20.10f \t\t grav: %20.10f \t\t rho: %20.10f \n",sed.kin_visc, sed.grav, sed.rho);

    printf("ice_roughness: %20.10f \t\t roughness coefficient: %20.10f \n",sed.ice_roughness, sed.roughness_coefficient);
    printf("total_erodable_thickness: %20.10f \t\t minimum_erodable_thickness: %20.10f\n",sed.total_erodable_thickness, sed.minimum_erodable_thickness);
    printf("strata_dep_coef: %20.10f \n",sed.strata_dep_coef);
    printf("total_shear_stress_mag; %20.10f \n",sed.total_shear_stress_mag);
    printf("grain_shear_stress_mag: %20.10f \n",sed.grain_shear_stress_mag);
    printf("near_bed_velocity_mag: %20.10f \n",sed.near_bed_velocity_mag);
    printf("time_step: %20.10f \n",sed.time_step);
    printf("ero_limit_ts: %20.10f \n",sed.ero_limit_ts);
    printf("sl_mfcf: %20.10f \n",sed.sl_mfcf);
    printf("esrh: %20.10f \t\t  tsrh: %20.10f \n",sed.esrh, sed.tsrh);

    printf("grain_shear_stress: "); svect2d_printScreen(sed.grain_shear_stress);
    printf("v_da: "); svect2d_printScreen(sed.v_da);
    printf("v_nb: "); svect2d_printScreen(sed.v_nb);
    printf("v_wind: "); svect2d_printScreen(sed.v_wind);
    printf("v_vort: "); svect2d_printScreen(sed.v_vort);
    printf("bl_vel: "); svect2d_printScreen(sed.bl_vel);
    printf("bed_slope: "); svect2d_printScreen(sed.bed_slope);
    printf("sl_vcf: "); svect2d_printScreen(sed.sl_vcf);

    printf("coh_settling_constants: ");
    for (i=0; i<4; i++) {
        printf("%20.10f \t",sed.coh_settling_constants[i]);
    }
    printf("\n");

    printf("al_coh_prop: ");
    for (i=0; i<5; i++) {
        printf("%20.10f \t",sed.al_coh_prop[i]);
    }
    printf("\n");


    if (sed.n_sedfl > 0) {
        for (i = 0; i < sed.n_blay; i++) {
            for (j = 0; j < sed.n_sedfl; j++) {
                printf("sedfl_shear[%d][%d]: %20.10e \t\t sedfl_erate: %20.10e\n",i,j,sed.sedfl_shear[i][j], sed.sedfl_erate[i][j]);        
            }
        }
    }

    if (sed.n_conti > 0) {
        for (i = 0; i < sed.n_conti; i++) {
            for (j=0; j<5; j++) {
                printf("consolidate_table[%d][%d]: %20.10e\n",i,j,sed.consolidate_table[i][j]);
            }
        }
    }

    for (i = 0; i < sed.n_blay; i++) {
        for (j=0; j<sed.n_sed; j++) {
            printf("bl_dist[%d][%d]: %20.10f\n",i,j,sed.bl_dist[i][j]);
        }
    }

    for (i = 0; i < sed.n_blay; i++) {
        printf("bed_layer: %d --------------\n",i);
        ssedlib_bed_layer_printScreen(sed.layer[i]);
    }

    for (i=0; i<sed.n_sed+1; i++) {
        printf("grain: %d ------------------\n",i);
        ssedlib_grain_printScreen(sed.grain[i]);
    }

}

//**********************************************************************************//
//**********************************************************************************//
//**********************************************************************************//

void ssedlib_grain_printScreen(SSEDLIB_GRAIN grn) {

    printf("type: %d\n",grn.type);
    printf("bed_mass_change: %20.10f\n",grn.bed_mass_change);
    printf("net_mass_change: %20.10f\n",grn.net_mass_change);
    printf("diameter: %20.10f\n",grn.diameter);
    printf("settling_vel: %20.10f\n",grn.settling_vel);
    printf("coh_free_settling_vel: %20.10f\n",grn.coh_free_settling_vel);
    printf("por: %20.10f \t ces: %20.10f \t cds: %20.10f \t erc: %20.10f \t ere: %20.10f\n",grn.por, grn.ces, grn.cds, grn.erc, grn.ere);
    printf("critical_shear_stress: %20.10f\n",grn.critical_shear_stress);
    printf("critical_shear_stress_plane_bed: %20.10f\n",grn.critical_shear_stress_plane_bed);
    printf("sg: %20.10f\n",grn.sg);
    printf("al_dist: %20.10f\n",grn.al_dist);
    printf("as_dist: %20.10f\n",grn.as_dist);
    printf("bl_con: %20.10f\n",grn.bl_con);
    printf("bl_thick: %20.10f\n",grn.bl_thick);
    printf("sl_con: %20.10f\n",grn.sl_con);
    printf("sl_bcon: %20.10f\n",grn.sl_bcon);
    printf("sl_sg: %20.10f\n",grn.sl_sg);
    printf("bl_erosion_flux: %20.10f\n",grn.bl_erosion_flux);
    printf("bl_deposition_flux: %20.10f\n",grn.bl_deposition_flux);
    printf("sl_erosion_flux: %20.10f\n",grn.sl_erosion_flux);
    printf("sl_deposition_flux: %20.10f\n",grn.sl_deposition_flux);
    printf("sl_erosion_flux_nbl: %20.10f\n",grn.sl_erosion_flux_nbl);
    printf("bl_erosion_flux_nbl: %20.10f\n",grn.bl_erosion_flux_nbl);
    printf("bl_deposition_flux_nbl: %20.10f\n",grn.bl_deposition_flux_nbl);
    printf("bl_net_flux: %20.10f\n",grn.bl_net_flux);
    printf("sl_net_flux: %20.10f\n",grn.sl_net_flux);
    printf("sl_nbf: %20.10f\n",grn.sl_nbf);
    printf("old_al_dist: %20.10f\n",grn.old_al_dist);

}


//**********************************************************************************//
//**********************************************************************************//
//**********************************************************************************//

void ssedlib_bed_layer_printScreen(SSEDLIB_BED_LAYER layer) {

    printf("type: %d\n",layer.type);
    printf("sedflume_flag: %d\n",layer.sedflume_flag);
    printf("fully_eroded_layer: %d\n",layer.fully_eroded_layer);

    printf("por: %20.10f \n",layer.por);
    printf("ces: %20.10f \n",layer.ces);
    printf("erc: %20.10f \n",layer.erc);
    printf("ere: %20.10f \n",layer.ere);
    printf("thick: %20.10f \n",layer.thick);
}

//**********************************************************************************//
//**********************************************************************************//
//**********************************************************************************//

void ssedlib_wave_printScreen(SSEDLIB_WAVE wave) {
    printf("ibreak: %d\n",wave.ibreak);
    printf("height: %20.10f \n",wave.height);
    printf("period: %20.10f \n",wave.period);
    printf("angle: %20.10f \n",wave.angle);
    printf("number: %20.10f \n",wave.number);
    printf("speed: %20.10f \n",wave.speed);
    printf("edb: %20.10f \n",wave.edb);
    printf("edf: %20.10f \n",wave.edf);
    printf("sef: %20.10f \n",wave.sef);
    printf("seb: %20.10f \n",wave.seb);
    printf("wcf: %20.10f \n",wave.wcf);
    printf("beta: %20.10f \n",wave.beta);
    printf("gam: %20.10f \n",wave.gama);
}
    



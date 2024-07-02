#ifdef _SEDIMENT
#include "global_header.h"
/******************************************************************************/
/* PROTOTYPES FOR FILE SCOPE FUNCTIONS ****************************************/
void ssediment_prep(SMODEL *, SSED *, int);

/******************************************************************************/
/******************************************************************************/
/* allocate and initialize suspended load struct */
void ssusload_alloc(SSUSLOAD **susload, int ngrains, int nnodes_bed, int nnodes_sus) {
    assert(nnodes_bed>0 && nnodes_sus>0 && ngrains>0);

    // allocate suspended load grains array
    (*susload) = (SSUSLOAD *) tl_alloc(sizeof(SSUSLOAD), ngrains);
    SSUSLOAD *load = (*susload);
    
    int ised=0;
    for (ised=0; ised<ngrains; ised++) {
        load[ised].c = (double *) tl_alloc(sizeof(double), nnodes_sus);
        load[ised].old_c = (double *) tl_alloc(sizeof(double), nnodes_sus);
        load[ised].older_c = (double *) tl_alloc(sizeof(double), nnodes_sus);
        load[ised].source = (double *) tl_alloc(sizeof(double), nnodes_sus);
        load[ised].sink = (double *) tl_alloc(sizeof(double), nnodes_sus);
        load[ised].error = (double *) tl_alloc(sizeof(double), nnodes_sus);
        load[ised].rouse_coef = (double *) tl_alloc(sizeof(double), nnodes_bed); // just on the bed
    
        // just on bed and only needed in 2d
        load[ised].mfcf = NULL;
        load[ised].vcf = NULL;
        load[ised].vor_vel = NULL;
        if (nnodes_sus == nnodes_bed) { // then we are in 2d
            load[ised].mfcf = (double *) tl_alloc(sizeof(double), nnodes_bed);
            load[ised].vcf = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes_bed);
            load[ised].vor_vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes_bed);
        }
    }
}

/******************************************************************************/
/* re-allocate and initialize suspended load struct */
void ssusload_realloc_sus(SSUSLOAD *load, int ngrains, int nnodes_sus_new, int nnodes_sus_old) {
    assert( nnodes_sus_new>0 && ngrains>0);
    
    int ised=0;
    for (ised=0; ised<ngrains; ised++) {
        load[ised].c = (double *) tl_realloc(sizeof(double), nnodes_sus_new, nnodes_sus_old, load[ised].c);
        load[ised].old_c = (double *) tl_realloc(sizeof(double), nnodes_sus_new, nnodes_sus_old, load[ised].old_c);
        load[ised].older_c = (double *) tl_realloc(sizeof(double), nnodes_sus_new, nnodes_sus_old, load[ised].older_c);
        load[ised].source = (double *) tl_realloc(sizeof(double), nnodes_sus_new, nnodes_sus_old, load[ised].source);
        load[ised].sink = (double *) tl_realloc(sizeof(double), nnodes_sus_new, nnodes_sus_old, load[ised].sink);
        load[ised].error = (double *) tl_realloc(sizeof(double), nnodes_sus_new, nnodes_sus_old, load[ised].error);
        if (load[ised].mfcf != NULL) { // then we are in 2d
            load[ised].mfcf = (double *) tl_realloc(sizeof(double), nnodes_sus_new, nnodes_sus_old, load[ised].mfcf);
            load[ised].vcf = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_sus_new, nnodes_sus_old, load[ised].vcf);
            load[ised].vor_vel = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_sus_new, nnodes_sus_old, load[ised].vor_vel);
        }
    }
}

/******************************************************************************/
/* initialize suspended load struct */
void ssusload_init_range_sus(SSUSLOAD *load, int ngrains, int nnodes_sus1, int nnodes_sus2) {
    assert(nnodes_sus1>=0 && nnodes_sus2>0 && ngrains > 0);
    
    int ised=0;
    for (ised=0; ised<ngrains; ised++) {
        sarray_dbl_init_value_range(load[ised].c, 0., nnodes_sus1, nnodes_sus2);
        sarray_dbl_init_value_range(load[ised].old_c, 0., nnodes_sus1, nnodes_sus2);
        sarray_dbl_init_value_range(load[ised].older_c, 0., nnodes_sus1, nnodes_sus2);
        sarray_dbl_init_value_range(load[ised].source, 0., nnodes_sus1, nnodes_sus2);
        sarray_dbl_init_value_range(load[ised].sink, 0.,nnodes_sus1, nnodes_sus2);
        sarray_dbl_init_value_range(load[ised].error, 0., nnodes_sus1, nnodes_sus2);
            
        if (load[ised].mfcf != NULL) {  // just on bed and only needed in 2d
            sarray_dbl_init_value_range(load[ised].mfcf, 1., nnodes_sus1, nnodes_sus2);
            svect2d_init_array_value_range(load[ised].vcf, 0., 0., nnodes_sus1, nnodes_sus2);
            svect2d_init_array_value_range(load[ised].vor_vel, 0., 0., nnodes_sus1, nnodes_sus2);
        }
    }
}

/*----------------------------------------------------------------------------*/
/* avg suspended load node - for adaption (only for 2d right now!) */
void ssusload_node_avg(SSUSLOAD *load, int ngrains, int node_new, int node1, int node2) {
    assert(ngrains>0);

    int ised=0;
    for (ised=0; ised<ngrains; ised++) {
        load[ised].c[node_new] = 0.5*(load[ised].c[node1] + load[ised].c[node2]);
        load[ised].old_c[node_new] = 0.5*(load[ised].old_c[node1] + load[ised].old_c[node2]);
        load[ised].older_c[node_new] = 0.5*(load[ised].older_c[node1] + load[ised].older_c[node2]);
        load[ised].source[node_new] = 0.5*(load[ised].source[node1] + load[ised].source[node2]);
        load[ised].sink[node_new] = 0.5*(load[ised].sink[node1] + load[ised].sink[node2]);
        load[ised].error[node_new] = 0.5*(load[ised].error[node1] + load[ised].error[node2]);
        load[ised].rouse_coef[node_new] = 0.5*(load[ised].rouse_coef[node1] + load[ised].rouse_coef[node2]);
        load[ised].mfcf[node_new] = 0.5*(load[ised].mfcf[node1] + load[ised].mfcf[node2]);
        load[ised].vcf[node_new] = svect2d_avg(load[ised].vcf[node1], load[ised].vcf[node2]);
        load[ised].vor_vel[node_new] = svect2d_avg(load[ised].vor_vel[node1], load[ised].vor_vel[node2]);
    }
}

/*----------------------------------------------------------------------------*/
/* renumber suspended load nodes - for adaption (only for 2d right now!) */
void ssusload_renumber(SSUSLOAD *load, int ngrains, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp, SVECT2D *vtmp) {
    assert(ngrains>0 && max_nnode>0);

    int ised=0;
    for (ised=0; ised<ngrains; ised++) {
        node_renumber_double(max_nnode, load[ised].c, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].old_c, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].older_c, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].source, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].sink, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].error, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].rouse_coef, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].mfcf, dtmp, new_numbers, order_tmp);
        node_renumber_vect2d(max_nnode, load[ised].vcf, vtmp, new_numbers, order_tmp);
        node_renumber_vect2d(max_nnode, load[ised].vor_vel, vtmp, new_numbers, order_tmp);
    }
}


/******************************************************************************/
/* re-allocate and initialize suspended load struct */
void ssusload_realloc_bed(SSUSLOAD *load, int ngrains, int nnodes_bed_new, int nnodes_bed_old) {
    assert(nnodes_bed_new>0 && ngrains>0);

    int ised=0;
    for (ised=0; ised<ngrains; ised++) {
        load[ised].rouse_coef = (double *) tl_realloc(sizeof(double), nnodes_bed_new, nnodes_bed_old, load[ised].rouse_coef); // just on the bed
    }
}

/******************************************************************************/
/* initialize suspended load struct */
void ssusload_init_range_bed(SSUSLOAD *load, int ngrains, int nnodes_bed1, int nnodes_bed2) {
    assert(nnodes_bed1>=0 && nnodes_bed2>0 && ngrains > 0);

    int ised=0;
    for (ised=0; ised<ngrains; ised++) {
      sarray_dbl_init_value_range(load[ised].rouse_coef, 1., nnodes_bed1, nnodes_bed2); // just on the bed
    }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/* free suspended load struct */
void ssusload_free(SSUSLOAD *load, int nsed, int nnodes_bed, int nnodes_sus) {
    assert(nsed>0 && nnodes_bed>0 && nnodes_sus>0);

    int ised=0;
    for (ised=0; ised<nsed; ised++) {
        // just on bed and only needed in 2d
        if (load[ised].vor_vel != NULL) load[ised].vor_vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes_bed, load[ised].vor_vel);
        if (load[ised].vcf != NULL) load[ised].vcf = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes_bed, load[ised].vcf);
        if (load[ised].mfcf != NULL) load[ised].mfcf = (double *) tl_free(sizeof(double), nnodes_bed, load[ised].mfcf);
        load[ised].rouse_coef = (double *) tl_free(sizeof(double), nnodes_bed, load[ised].rouse_coef);
    
        // these variables are on the complete 2d/3d grid
        load[ised].error = (double *) tl_free(sizeof(double), nnodes_sus, load[ised].error);
        load[ised].source = (double *) tl_free(sizeof(double), nnodes_sus, load[ised].source);
        load[ised].sink = (double *) tl_free(sizeof(double), nnodes_sus, load[ised].sink);
        load[ised].older_c = (double *) tl_free(sizeof(double), nnodes_sus, load[ised].older_c);
        load[ised].old_c = (double *) tl_free(sizeof(double), nnodes_sus, load[ised].old_c);
        load[ised].c = (double *) tl_free(sizeof(double), nnodes_sus, load[ised].c);
    }

    load = (SSUSLOAD *) tl_free(sizeof(SSUSLOAD), nsed, load);
}

/******************************************************************************/
/******************************************************************************/
/* allocate and initialize bed load struct */
void sbedload_alloc(SBEDLOAD **bedload, int nsed, int nnodes) {
    assert(nsed>0 && nnodes>0);

    // allocate bed load grains array
    (*bedload) = (SBEDLOAD *) tl_alloc(sizeof(SBEDLOAD), nsed);
    SBEDLOAD *load = (*bedload);

    int ised=0;
    for (ised=0; ised<nsed; ised++) {
        load[ised].c = (double *) tl_alloc(sizeof(double), nnodes);
        load[ised].old_c = (double *) tl_alloc(sizeof(double), nnodes);
        load[ised].older_c = (double *) tl_alloc(sizeof(double), nnodes);
        load[ised].thick = (double *) tl_alloc(sizeof(double), nnodes);
        load[ised].old_thick = (double *) tl_alloc(sizeof(double), nnodes);
        load[ised].older_thick = (double *) tl_alloc(sizeof(double), nnodes);
        load[ised].source = (double *) tl_alloc(sizeof(double), nnodes);
        load[ised].sink = (double *) tl_alloc(sizeof(double), nnodes);
        load[ised].shear = (double *) tl_alloc(sizeof(double), nnodes);
        load[ised].error = (double *) tl_alloc(sizeof(double), nnodes);
        load[ised].v = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
        load[ised].flux = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    }
}

/******************************************************************************/
/* re-allocate and initialize bed load struct */
void sbedload_realloc(SBEDLOAD *load, int nsed, int nnodes_new, int nnodes_old) {
    assert(nsed>0 && nnodes_old>=0 && nnodes_new>0);

    int ised=0;
    for (ised=0; ised<nsed; ised++) {
        load[ised].c = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, load[ised].c);
        load[ised].old_c = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, load[ised].old_c);
        load[ised].older_c = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, load[ised].older_c);
        load[ised].thick = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, load[ised].thick);
        load[ised].old_thick = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, load[ised].old_thick);
        load[ised].older_thick = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, load[ised].older_thick);
        load[ised].shear = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, load[ised].shear);
        load[ised].source = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, load[ised].source);
        load[ised].sink = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, load[ised].sink);
        load[ised].error = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, load[ised].error);
        load[ised].v = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, load[ised].v);
        load[ised].flux = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, load[ised].flux);
    }
}

/******************************************************************************/
/* initialize bed load struct */
void sbedload_init_range(SBEDLOAD *load, int nsed, int nnodes1, int nnodes2) {
    assert(nsed>0 && nnodes1>=0 && nnodes2>0);

    int ised=0;
    for (ised=0; ised<nsed; ised++) {
        sarray_dbl_init_value_range(load[ised].c, 0., nnodes1, nnodes2);
        sarray_dbl_init_value_range(load[ised].old_c, 0., nnodes1, nnodes2);
        sarray_dbl_init_value_range(load[ised].older_c, 0., nnodes1, nnodes2);
        sarray_dbl_init_value_range(load[ised].thick, 1., nnodes1, nnodes2);
        sarray_dbl_init_value_range(load[ised].old_thick, 1., nnodes1, nnodes2);
        sarray_dbl_init_value_range(load[ised].older_thick, 0., nnodes1, nnodes2);
        sarray_dbl_init_value_range(load[ised].source, 0., nnodes1, nnodes2);
        sarray_dbl_init_value_range(load[ised].sink, 0.,nnodes1, nnodes2);
        sarray_dbl_init_value_range(load[ised].shear, 0.,nnodes1, nnodes2);
        sarray_dbl_init_value_range(load[ised].error, 0., nnodes1, nnodes2);
        svect2d_init_array_value_range(load[ised].v, 0., 0., nnodes1, nnodes2);
        svect2d_init_array_value_range(load[ised].flux, 0., 0., nnodes1, nnodes2);
    }
}

/*----------------------------------------------------------------------------*/
/* avg bed load node - for adaption (only for 2d right now!) */
void sbedload_node_avg(SBEDLOAD *load, int nsed, int node_new, int node1, int node2) {
    assert(nsed>0 && node1>=0 && node2>=0 && node_new>=0);
    int ised=0;
    for (ised=0; ised<nsed; ised++) { 
        load[ised].c[node_new] = 0.5*(load[ised].c[node1] + load[ised].c[node2]);
        load[ised].old_c[node_new] = 0.5*(load[ised].old_c[node1] + load[ised].old_c[node2]);
        load[ised].older_c[node_new] = 0.5*(load[ised].older_c[node1] + load[ised].older_c[node2]);
        load[ised].thick[node_new] = 0.5*(load[ised].thick[node1] + load[ised].thick[node2]);
        load[ised].old_thick[node_new] = 0.5*(load[ised].old_thick[node1] + load[ised].old_thick[node2]);
        load[ised].older_thick[node_new] = 0.5*(load[ised].older_thick[node1] + load[ised].older_thick[node2]);
        load[ised].source[node_new] = 0.5*(load[ised].source[node1] + load[ised].source[node2]);
        load[ised].sink[node_new] = 0.5*(load[ised].sink[node1] + load[ised].sink[node2]);
        load[ised].shear[node_new] = 0.5*(load[ised].shear[node1] + load[ised].shear[node2]);
        load[ised].error[node_new] = 0.5*(load[ised].error[node1] + load[ised].error[node2]);
        load[ised].v[node_new] = svect2d_avg(load[ised].v[node1], load[ised].v[node2]);
        load[ised].flux[node_new] = svect2d_avg(load[ised].flux[node1], load[ised].flux[node2]);
    }
}

/*----------------------------------------------------------------------------*/
/* renumber bed load nodes - for adaption (only for 2d right now!) */
void sbedload_renumber(SBEDLOAD *load, int nsed, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp, SVECT2D *vtmp) {
    assert(nsed>0 && max_nnode >=0);
    int ised=0;
    for (ised=0; ised<nsed; ised++) {
        node_renumber_double(max_nnode, load[ised].c, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].old_c, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].older_c, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].thick, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].old_thick, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].older_thick, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].source, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].sink, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].shear, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, load[ised].error, dtmp, new_numbers, order_tmp);
        node_renumber_vect2d(max_nnode, load[ised].v, vtmp, new_numbers, order_tmp);
        node_renumber_vect2d(max_nnode, load[ised].flux, vtmp, new_numbers, order_tmp);
    }
}

/*----------------------------------------------------------------------------*/
/* free bed load struct */
void sbedload_free(SBEDLOAD *load, int nsed, int nnodes) {
    assert(nsed>0 && nnodes>0);
    int ised=0;
    for (ised=0; ised<nsed; ised++) {
        load[ised].flux = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, load[ised].flux);
        load[ised].v = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, load[ised].v);
        load[ised].error = (double *) tl_free(sizeof(double), nnodes, load[ised].error);
        load[ised].source = (double *) tl_free(sizeof(double), nnodes, load[ised].source);
        load[ised].sink = (double *) tl_free(sizeof(double), nnodes, load[ised].sink); 
        load[ised].shear = (double *) tl_free(sizeof(double), nnodes, load[ised].shear);
        load[ised].older_thick = (double *) tl_free(sizeof(double), nnodes, load[ised].older_thick);
        load[ised].old_thick = (double *) tl_free(sizeof(double), nnodes, load[ised].old_thick);
        load[ised].thick = (double *) tl_free(sizeof(double), nnodes, load[ised].thick);
        load[ised].older_c = (double *) tl_free(sizeof(double), nnodes, load[ised].older_c);
        load[ised].old_c = (double *) tl_free(sizeof(double), nnodes, load[ised].old_c);
        load[ised].c = (double *) tl_free(sizeof(double), nnodes, load[ised].c);
    }
    load = (SBEDLOAD *) tl_free(sizeof(SBEDLOAD), nsed, load);
}

/******************************************************************************/
/******************************************************************************/
/* allocate and initialize bed layers */
/* initialize bed layer struct */
void sbed_layer_alloc(SBED_LAYER **bedlayer, int nlayers, int ngrains, int nnodes) {
    
    (*bedlayer) = (SBED_LAYER *) tl_alloc(sizeof(SBED_LAYER), nlayers);
    SBED_LAYER *layer = (*bedlayer);
    
    int ilayer;
    for (ilayer=0; ilayer<nlayers; ilayer++) {
        layer[ilayer].type = 0;
        layer[ilayer].cbed_flag = NO;
        layer[ilayer].sedflume_flag = 0;

        /* allocate and initialize sedflume */
        ssedflume_alloc_init(&(layer[ilayer].sedflume), nnodes);

        layer[ilayer].bulk_density = (double *) tl_alloc(sizeof(double), nnodes);
        layer[ilayer].porosity = (double *) tl_alloc(sizeof(double), nnodes);
        layer[ilayer].critical_erosion_shear = (double *) tl_alloc(sizeof(double), nnodes);
        layer[ilayer].erosion_rate_constant = (double *) tl_alloc(sizeof(double), nnodes);
        layer[ilayer].erosion_rate_exponent = (double *) tl_alloc(sizeof(double), nnodes);
        layer[ilayer].thickness = (double *) tl_alloc(sizeof(double), nnodes);
    
        layer[ilayer].distribution = (double **) tl_alloc(sizeof(double *), ngrains);
        int i;
        for (i=0; i<ngrains; i++) {
            layer[ilayer].distribution[i] = (double *) tl_alloc(sizeof(double), nnodes);
        }
    }
}

/******************************************************************************/
/* re-allocate and initialize bed load struct */
void sbed_layer_realloc(SBED_LAYER *layer, int nlayers, int ngrains, int nnodes_new, int nnodes_old) {
    assert(nlayers>0 && ngrains>0 && nnodes_new>0);
    int ilayer;
    for (ilayer=0; ilayer<nlayers; ilayer++) {
        ssedflume_realloc_init(layer[ilayer].sedflume, nnodes_new, nnodes_old);
        layer[ilayer].bulk_density = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, layer[ilayer].bulk_density);
        layer[ilayer].porosity = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, layer[ilayer].porosity);
        layer[ilayer].critical_erosion_shear = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, layer[ilayer].critical_erosion_shear);
        layer[ilayer].erosion_rate_constant = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, layer[ilayer].erosion_rate_constant);
        layer[ilayer].erosion_rate_exponent = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, layer[ilayer].erosion_rate_exponent);
        layer[ilayer].thickness = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, layer[ilayer].thickness);
        int i;
        for (i=0; i<ngrains; i++) {
            layer[ilayer].distribution[i] = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, layer[ilayer].distribution[i]);
        }
    }
}

/******************************************************************************/
/* initialize bed load struct */
void sbed_layer_init_range(SBED_LAYER *layer, int nlayers, int ngrains, int nnodes1, int nnodes2) {
    assert(nlayers>0 && ngrains>0 && nnodes1>=0 && nnodes2>0);
    
    int ilayer;
    for (ilayer=0; ilayer<nlayers; ilayer++) {
        ssedflume_init_range(layer[ilayer].sedflume, nnodes1, nnodes2);
        sarray_dbl_init_value_range(layer[ilayer].bulk_density, 0., nnodes1, nnodes2);
        sarray_dbl_init_value_range(layer[ilayer].porosity, 0., nnodes1, nnodes2);
        sarray_dbl_init_value_range(layer[ilayer].critical_erosion_shear, 0., nnodes1, nnodes2);
        sarray_dbl_init_value_range(layer[ilayer].erosion_rate_constant, 0., nnodes1, nnodes2);
        sarray_dbl_init_value_range(layer[ilayer].erosion_rate_exponent, 0., nnodes1, nnodes2);
        sarray_dbl_init_value_range(layer[ilayer].thickness, 0., nnodes1, nnodes2);
        int i;
        for (i=0; i<ngrains; i++) {
            sarray_dbl_init_value_range(layer[ilayer].distribution[i], 0., nnodes1, nnodes2);
        }
    }
}

/*----------------------------------------------------------------------------*/
/* avg bed layer node - for adaption (only for 2d right now!) */
void sbed_layer_node_avg(SBED_LAYER *layer, int nlayers, int ngrains, int node_new, int node1, int node2, double weight1, double weight2, int active_layer_flag) {
    assert(nlayers>0 && ngrains>0 && node_new>=0 && node1>=0 && node2>=0 && weight1>=0 && weight2>=0);

    int i=UNSET_INT, ilayer=UNSET_INT;
    double sum = 0., weight1_new = 0., weight2_new = 0.;

    for (ilayer=0; ilayer<nlayers; ilayer++) { // only calling this function one layer at a time, otherwise below won't work
        ssedflume_node_avg(layer[ilayer].sedflume, node_new, node1, node2);
        layer[ilayer].bulk_density[node_new] = weight1*layer[ilayer].bulk_density[node1] + weight2*layer[ilayer].bulk_density[node2];
        layer[ilayer].porosity[node_new] = weight1*layer[ilayer].porosity[node1] + weight2*layer[ilayer].porosity[node2];
        layer[ilayer].critical_erosion_shear[node_new] = weight1*layer[ilayer].critical_erosion_shear[node1] + weight2*layer[ilayer].critical_erosion_shear[node2];
        layer[ilayer].erosion_rate_constant[node_new] = weight1*layer[ilayer].erosion_rate_constant[node1] + weight2*layer[ilayer].erosion_rate_constant[node2];
        layer[ilayer].erosion_rate_exponent[node_new] = weight1*layer[ilayer].erosion_rate_exponent[node1] + weight2*layer[ilayer].erosion_rate_exponent[node2];
        layer[ilayer].thickness[node_new] = weight1*layer[ilayer].thickness[node1] + weight2*layer[ilayer].thickness[node2];
    
        weight1_new = weight1;
        weight2_new = weight2;
        if (active_layer_flag) {
            weight1_new *= (1. - layer[ilayer].porosity[node1]);
            weight2_new *= (1. - layer[ilayer].porosity[node2]);
        } else {
            weight1_new *= layer[ilayer].bulk_density[node1];
            weight2_new *= layer[ilayer].bulk_density[node2];
        }
        sum = weight1_new + weight2_new;
        weight1_new /= sum;
        weight2_new /= sum;
    
        for (i=0; i<ngrains; i++) {
            layer[ilayer].distribution[i][node_new] = weight1_new*layer[ilayer].distribution[i][node1] + weight2_new*layer[ilayer].distribution[i][node2];
        }
    }
}

/*----------------------------------------------------------------------------*/
/* renumber bed layer nodes - for adaption (only for 2d right now!) */
void sbed_layer_renumber(SBED_LAYER *layer, int nlayers, int ngrains, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp) {
    assert(nlayers>0 && ngrains>0 && max_nnode>=0);
    int ilayer;
    for (ilayer=0; ilayer<nlayers; ilayer++) {
        ssedflume_renumber(layer[ilayer].sedflume, max_nnode, new_numbers, order_tmp, dtmp);
        node_renumber_double(max_nnode, layer[ilayer].bulk_density, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, layer[ilayer].porosity, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, layer[ilayer].critical_erosion_shear, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, layer[ilayer].erosion_rate_constant, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, layer[ilayer].erosion_rate_exponent, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, layer[ilayer].thickness, dtmp, new_numbers, order_tmp);
        int i;
        for (i=0; i<ngrains; i++) {
            node_renumber_double(max_nnode, layer[ilayer].distribution[i], dtmp, new_numbers, order_tmp);
        }
    }
}

/*----------------------------------------------------------------------------*/
/* free bed layer struct */
void sbed_layer_free(SBED_LAYER *layer, int nlayers, int ngrains, int nnodes) {
    assert(nlayers>0 && ngrains>0 && nnodes>=0); 
    int i, ilayer;
    for (ilayer=0; ilayer<nlayers; ilayer++) {
        for (i=0; i<ngrains; i++) {
            layer[ilayer].distribution[i] = (double *) tl_free(sizeof(double), nnodes, layer[ilayer].distribution[i]);
        }
        layer[ilayer].distribution = (double **) tl_free(sizeof(double *), ngrains, layer[ilayer].distribution);
    
        layer[ilayer].bulk_density = (double *) tl_free(sizeof(double), nnodes, layer[ilayer].bulk_density);
        layer[ilayer].porosity = (double *) tl_free(sizeof(double), nnodes, layer[ilayer].porosity);
        layer[ilayer].critical_erosion_shear = (double *) tl_free(sizeof(double), nnodes, layer[ilayer].critical_erosion_shear);
        layer[ilayer].erosion_rate_constant = (double *) tl_free(sizeof(double), nnodes, layer[ilayer].erosion_rate_constant);
        layer[ilayer].erosion_rate_exponent = (double *) tl_free(sizeof(double), nnodes, layer[ilayer].erosion_rate_exponent);
        layer[ilayer].thickness = (double *) tl_free(sizeof(double), nnodes, layer[ilayer].thickness);
    
        ssedflume_free(layer[ilayer].sedflume,nnodes);
    }
    layer = (SBED_LAYER *) tl_free(sizeof(SBED_LAYER), nlayers, layer);
    
}

/*----------------------------------------------------------------------------*/
/* copy bed layer struct */
void sbed_layer_copy(SBED_LAYER *to, SBED_LAYER *from, int nlayers, int ngrains, int nnodes) {
    assert(nlayers>0 && ngrains>0 && nnodes>=0); 
    to->type = from->type;
    to->sedflume_flag = from->sedflume_flag;
    
    int ilayer = 0, igrain = 0, inode = 0;
    for (ilayer=0; ilayer<nlayers; ilayer++) {
        for (inode=0; inode<nnodes; inode++) {
            to[ilayer].bulk_density[inode] = from[ilayer].bulk_density[inode];
            to[ilayer].porosity[inode] = from[ilayer].porosity[inode];
            to[ilayer].critical_erosion_shear[inode] = from[ilayer].critical_erosion_shear[inode];
            to[ilayer].erosion_rate_constant[inode] = from[ilayer].erosion_rate_constant[inode];
            to[ilayer].erosion_rate_exponent[inode] = from[ilayer].erosion_rate_exponent[inode];
            to[ilayer].thickness[inode] = from[ilayer].thickness[inode];
        }
        for (igrain=0; igrain<ngrains; igrain++) {
            for (inode=0; inode<nnodes; inode++) {
                to[ilayer].distribution[igrain][inode] = from[ilayer].distribution[igrain][inode];
            }
        }
    }
}

/******************************************************************************/
/******************************************************************************/
/* initialize clay */
void sgrain_init(SGRAIN *grain) {
    grain->type = 0;
    grain->specific_gravity = 0.;
    grain->diameter = 0.;
    grain->porosity = 0.;
    sclay_init(&(grain->clay));
    ssand_init(&(grain->sand));
}

/*----------------------------------------------------------------------------*/
/* copy grain struct */
void sgrain_copy(SGRAIN *to, SGRAIN *from) {
    to->type = from->type;
    to->specific_gravity = from->specific_gravity;
    to->diameter = from->diameter;
    to->porosity = from->porosity;
    to->clay = from->clay;  // shallow copy is enough here
    to->sand = from->sand;  // shallow copy is enough here
}

/******************************************************************************/
/******************************************************************************/
/* initialize clay */
void sclay_init(SCLAY *clay) {
    clay->sed_to_clay = 0;
    clay->clay_to_sed = 0;
    clay->settling_velocity = 0.;
    clay->tau_ce = 0.;
    clay->tau_cd = 0.;
    clay->erode_const = 0;
}

// note :: if pointers are added to the above, shallow copies are not good enough

/******************************************************************************/
/******************************************************************************/
/* initialize sand */
void ssand_init(SSAND *sand) {
    sand->sed_to_snd = 0;
    sand->snd_to_sed = 0;
    sand->settling_velocity = 0.;
    sand->tau_ce = 0.;
    sand->tau_cd = 0.;
    sand->erode_const = 0;
}

// note :: if pointers are added to the above, shallow copies are not good enough

/******************************************************************************/
/******************************************************************************/
/* allocate and initialize sedflume arrays */
void ssedflume_alloc_init(SSEDFLUME **sedflume, int nnodes) {
    (*sedflume) = (SSEDFLUME *) tl_alloc(sizeof(SSEDFLUME), 1);
    SSEDFLUME *flume = (*sedflume);
    
    flume->shear_stress = (double *) tl_alloc(sizeof(double), nnodes);
    flume->erosion_rate = (double *) tl_alloc(sizeof(double), nnodes);
    ssedflume_init_range(flume, 0, nnodes);
}
/*----------------------------------------------------------------------------*/
void ssedflume_realloc_init(SSEDFLUME *flume, int nnodes_new, int nnodes_old) {
    flume->shear_stress = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, flume->shear_stress);
    flume->erosion_rate = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, flume->erosion_rate);
    ssedflume_init_range(flume, nnodes_old, nnodes_new);
}
/*----------------------------------------------------------------------------*/
void ssedflume_init_range(SSEDFLUME *flume, int nnodes_start, int nnodes_end) {
    sarray_dbl_init_value_range(flume->shear_stress, 0., nnodes_start, nnodes_end);
    sarray_dbl_init_value_range(flume->erosion_rate, 0., nnodes_start, nnodes_end);
}
/*----------------------------------------------------------------------------*/
void ssedflume_node_avg(SSEDFLUME *flume, int node_new, int node1, int node2) {
    flume->shear_stress[node_new] = 0.5 * (flume->shear_stress[node1] + flume->shear_stress[node2]);
    flume->erosion_rate[node_new] = 0.5 * (flume->erosion_rate[node1] + flume->erosion_rate[node2]);
}
/*----------------------------------------------------------------------------*/
void ssedflume_renumber(SSEDFLUME *flume, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp) {
    node_renumber_double(max_nnode, flume->shear_stress, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, flume->erosion_rate, dtmp, new_numbers, order_tmp);
}
/*----------------------------------------------------------------------------*/
/* free sedflume struct */
void ssedflume_free(SSEDFLUME *flume, int nnodes) {
    flume->shear_stress = (double *) tl_free(sizeof(double), nnodes, flume->shear_stress);
    flume->erosion_rate = (double *) tl_free(sizeof(double), nnodes, flume->erosion_rate);
    flume = (SSEDFLUME *) tl_free(sizeof(SSEDFLUME), 1, flume);
}

/******************************************************************************/
/******************************************************************************/
/* allocates ssediment struct */
void ssediment_alloc_init(int ndim, int nnodes_sus, int nnodes_bed, int nlayers, int nsand, int nclay, int nsilt, int ngrains, int nsfssi, int nconti, SSED **sediment) {

    int i=0, j=0, k=0;
    int icons = 0, inode = 0, ilayer = 0, igrain = 0;

	/* allocated sediment transport struct */
    (*sediment) = (SSED *) tl_alloc(sizeof(SSED), 1);
    SSED *sed = (*sediment);         // alias

    // default flags
    sed->wind_wave_flag = 0;
    sed->cohesive_settling_flag = 0;
    sed->noncoh_ent_flag = 0;
    sed->noncoh_ble_flag = 0;
    sed->hid_fact_flag = 0;
    sed->bed_type_flag = 0;

    // sediment counts
    sed->nnodes_sus = nnodes_sus;
    sed->nnodes_bed = nnodes_bed;
    sed->nse_params = 1;
    sed->nbe_params = 1;
    sed->nsand = nsand;
    sed->nclay = nclay;
    sed->nsilt = nsilt;
    sed->nsed = ngrains;
    sed->nlayers = nlayers;
    sed->nsfssi = nsfssi;
    sed->ncti = nconti;

    // arrays
	  sed->cohesive_settling_const = (double *) tl_alloc(sizeof(double), COHSET);
    sarray_dbl_init(sed->cohesive_settling_const, COHSET);

    sed->wind_wave_const = (double *) tl_alloc(sizeof(double), WINDWAVE);
    sarray_dbl_init(sed->wind_wave_const, WINDWAVE);

    sed->noncoh_ent_const = (double *) tl_alloc(sizeof(double), 1);
    *(sed->noncoh_ent_const) = 0.;

    sed->noncoh_ble_const = (double *) tl_alloc(sizeof(double), 1);
    *(sed->noncoh_ble_const) = 0.;

    sed->div_coef = (double *) tl_alloc(sizeof(double), DIVCOEF);
    sarray_dbl_init(sed->div_coef, DIVCOEF);
    
    sed->friction_flag = (int *) tl_alloc(sizeof(int), nnodes_bed);
    sed->call_flag = (int *) tl_alloc(sizeof(int), nnodes_bed);
    sed->no_displacement_flag = (int *) tl_alloc(sizeof(int), nnodes_bed);
    sed->friction_coef = (double *) tl_alloc(sizeof(double), nnodes_bed);
    sed->as_ceiling = (double *) tl_alloc(sizeof(double), nnodes_bed);
    sed->old_as_ceiling = (double *) tl_alloc(sizeof(double), nnodes_bed);
    sed->bed_displacement = (double *) tl_alloc(sizeof(double), nnodes_bed);
    sed->old_bed_displacement = (double *) tl_alloc(sizeof(double), nnodes_bed);
    sed->bed_shear_stress = (double *) tl_alloc(sizeof(double), nnodes_bed);

    sed->bedload_vector = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes_bed);
    sed->susload_vector = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes_sus);
    sed->bed_gradient = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes_bed);

    /* allocate consolidation tables */
    sed->consolidation_bed_properties = NULL;
    if (nconti > 0) {
        sed->consolidation_bed_properties = (double ***) tl_alloc(sizeof(double **), 5);
        for (k = 0; k < 5; k++) {
            sed->consolidation_bed_properties[k] = (double **) tl_alloc(sizeof(double *), nconti);
            for (j = 0; j < nconti; j++) {
                sed->consolidation_bed_properties[k][j] = (double *) tl_alloc(sizeof(double), nnodes_bed);
                for (i = 0; i < nnodes_bed; i++) {
                    sed->consolidation_bed_properties[k][j][i] = 0.;
                }
            }
        }
    }
    
    /* allocate and initialize grains */
    sed->grain = (SGRAIN *) tl_alloc(sizeof(SGRAIN), ngrains);
    for (igrain=0; igrain<ngrains; igrain++) {
        sgrain_init(&(sed->grain[igrain]));
    }
    
    /* allocate and initialize bed and suspended load structs */
    ssusload_alloc(&(sed->susload), ngrains, nnodes_bed, nnodes_sus);
    sbedload_alloc(&(sed->bedload), ngrains, nnodes_bed);

    // active and bed layers (sed->bed_layer[ilayer].thickness[inode]
    sbed_layer_alloc(&(sed->active_layer), ONE, ngrains, nnodes_bed);
    sbed_layer_alloc(&(sed->old_active_layer), ONE, ngrains, nnodes_bed);
    sbed_layer_alloc(&(sed->bed_layer), nlayers, ngrains, nnodes_bed);
    sbed_layer_alloc(&(sed->old_bed_layer), nlayers, ngrains, nnodes_bed);

    // allocate sedlib sediment variables
    int nsedfl = 0; // for now
    ssedlib_alloc(&(sed->sedlib), nlayers, ngrains, nsedfl, nconti);
    
    // intialize
    int old_nodes = 0;
    ssediment_init_range_bed(sed, old_nodes, nnodes_bed);
    ssediment_init_range_sus(sed, old_nodes, nnodes_sus);
}

/******************************************************************************/
/* re-allocate and initialize suspended sediment struct */
void ssediment_realloc_init_sus(SSED *sed, int nnodes_sus_old, int nnodes_sus_new) {
    
    assert(nnodes_sus_new>0);
    
    // re-allocate
    
    sed->susload_vector = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_sus_new, nnodes_sus_old, sed->susload_vector);
        
    ssusload_realloc_sus(sed->susload, sed->nsed, nnodes_sus_new, nnodes_sus_old);
    
    
    // initialize just new nodes
    ssediment_init_range_sus(sed, nnodes_sus_old, nnodes_sus_new);

    sed->nnodes_sus = nnodes_sus_new;

}

/******************************************************************************/

/******************************************************************************/
/* re-allocate and initialize bed sediment struct */
void ssediment_realloc_init_bed(SSED *sed, int nnodes_bed_old, int nnodes_bed_new) {

    assert(nnodes_bed_new>0);

    // re-allocate
    sed->call_flag = (int *) tl_realloc(sizeof(int), nnodes_bed_new, nnodes_bed_old, sed->call_flag);
    sed->friction_flag = (int*) tl_realloc(sizeof(int), nnodes_bed_new, nnodes_bed_old, sed->friction_flag);
    sed->no_displacement_flag = (int*) tl_realloc(sizeof(int), nnodes_bed_new, nnodes_bed_old, sed->no_displacement_flag);
    sed->as_ceiling = (double *) tl_realloc(sizeof(double), nnodes_bed_new, nnodes_bed_old, sed->as_ceiling);
    sed->old_as_ceiling = (double *) tl_realloc(sizeof(double), nnodes_bed_new, nnodes_bed_old, sed->old_as_ceiling);
    sed->bed_displacement = (double *) tl_realloc(sizeof(double), nnodes_bed_new, nnodes_bed_old, sed->bed_displacement);
    sed->old_bed_displacement = (double *) tl_realloc(sizeof(double), nnodes_bed_new, nnodes_bed_old, sed->old_bed_displacement);
    sed->bed_shear_stress = (double *) tl_realloc(sizeof(double), nnodes_bed_new, nnodes_bed_old, sed->bed_shear_stress);
    sed->friction_coef = (double *) tl_realloc(sizeof(double), nnodes_bed_new, nnodes_bed_old, sed->friction_coef);

    sed->bedload_vector = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_bed_new, nnodes_bed_old, sed->bedload_vector);
    sed->bed_gradient = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_bed_new, nnodes_bed_old, sed->bed_gradient);

    int j, k;
    if (sed->ncti > 0) {
        for (k = 0; k < 5; k++) {
            for (j = 0; j < sed->ncti; j++) {
                sed->consolidation_bed_properties[k][j] = (double *) tl_realloc(sizeof(double), nnodes_bed_new, nnodes_bed_old, sed->consolidation_bed_properties[k][j]);
            }
        }
    }

    ssusload_realloc_bed(sed->susload, sed->nsed, nnodes_bed_new, nnodes_bed_old);
    sbedload_realloc(sed->bedload, sed->nsed, nnodes_bed_new, nnodes_bed_old);

    sbed_layer_realloc(sed->active_layer, 1, sed->nsed, nnodes_bed_new, nnodes_bed_old);
    sbed_layer_realloc(sed->old_active_layer, 1, sed->nsed, nnodes_bed_new, nnodes_bed_old);
    sbed_layer_realloc(sed->bed_layer, sed->nlayers, sed->nsed, nnodes_bed_new, nnodes_bed_old);
    sbed_layer_realloc(sed->old_bed_layer, sed->nlayers, sed->nsed, nnodes_bed_new, nnodes_bed_old);

    // initialize just new nodes
    ssediment_init_range_bed(sed, nnodes_bed_old, nnodes_bed_new);

    sed->nnodes_bed = nnodes_bed_new;
}

/******************************************************************************/
/* initialize suspended load struct */
void ssediment_init_range_sus(SSED *sed, int nnodes_sus1, int nnodes_sus2) {
    
    assert(nnodes_sus1>=0);
    assert(nnodes_sus2>0);
    
    // these variables are on the complete 2d/3d grid
 
    svect2d_init_array_value_range(sed->susload_vector, 0., 0., nnodes_sus1, nnodes_sus2);
    
    ssusload_init_range_sus(sed->susload, sed->nsed, nnodes_sus1, nnodes_sus2);
 }

/*----------------------------------------------------------------------------*/

/******************************************************************************/
/* initialize suspended load struct */
void ssediment_init_range_bed(SSED *sed, int nnodes_bed1, int nnodes_bed2) {

    assert(nnodes_bed1>=0);
    assert(nnodes_bed2>0);

    // these variables are on the complete 2d/3d grid
    sarray_int_init_value_range(sed->call_flag, OFF, nnodes_bed1, nnodes_bed2); // ??
    sarray_int_init_value_range(sed->friction_flag, 0, nnodes_bed1, nnodes_bed2); // ??
    sarray_int_init_value_range(sed->no_displacement_flag, 0, nnodes_bed1, nnodes_bed2); // ??
    sarray_dbl_init_value_range(sed->as_ceiling, 0., nnodes_bed1, nnodes_bed2);
    sarray_dbl_init_value_range(sed->old_as_ceiling, 0.,nnodes_bed1, nnodes_bed2);
    sarray_dbl_init_value_range(sed->bed_displacement, 0., nnodes_bed1, nnodes_bed2);
    sarray_dbl_init_value_range(sed->old_bed_displacement, 0., nnodes_bed1, nnodes_bed2);
    sarray_dbl_init_value_range(sed->bed_shear_stress, 0., nnodes_bed1, nnodes_bed2);
    sarray_dbl_init_value_range(sed->friction_coef, 0., nnodes_bed1, nnodes_bed2);

    svect2d_init_array_value_range(sed->bedload_vector, 0., 0., nnodes_bed1, nnodes_bed2);
    svect2d_init_array_value_range(sed->bed_gradient, 0., 0., nnodes_bed1, nnodes_bed2);

    int i, j, k;
    if (sed->ncti > 0) {
        for (k = 0; k < 5; k++) {
            for (j = 0; j < sed->ncti; j++) {
                for (i = nnodes_bed1; i < nnodes_bed2; i++) {
                    sed->consolidation_bed_properties[k][j][i] = 0.;
                }
            }
        }
    }

    ssusload_init_range_bed(sed->susload, sed->nsed, nnodes_bed1, nnodes_bed2);
    sbedload_init_range(sed->bedload, sed->nsed, nnodes_bed1, nnodes_bed2);

    sbed_layer_init_range(sed->active_layer, 1, sed->nsed, nnodes_bed1, nnodes_bed2);
    sbed_layer_init_range(sed->old_active_layer, 1, sed->nsed, nnodes_bed1, nnodes_bed2);
    sbed_layer_init_range(sed->bed_layer, sed->nlayers, sed->nsed, nnodes_bed1, nnodes_bed2);
    sbed_layer_init_range(sed->old_bed_layer, sed->nlayers, sed->nsed, nnodes_bed1, nnodes_bed2);
}

/*----------------------------------------------------------------------------*/
/* avg sediment stuct node - for adaption (only for 2d right now!) */
void ssediment_node_avg(SSED *sed, int node_new, int node1, int node2, SGRID *grid) {
    int i, j, k, igrain, ilayer;
    int node_new_bed,node1_bed,node2_bed;

    if(grid->ndim==3){
      node_new_bed =  grid->nodeID_3d_to_2d_bed[node_new];
      node1_bed =  grid->nodeID_3d_to_2d_bed[node1];
      node2_bed =  grid->nodeID_3d_to_2d_bed[node2];
    }else{
      node_new_bed =  node_new;
      node1_bed =  node1;
      node2_bed = node2;
    }

  if((node_new_bed != UNSET_INT) && (node1_bed != UNSET_INT) && (node2_bed != UNSET_INT)){    
    sed->friction_flag[node_new_bed] = sed->friction_flag[node1_bed];
    sed->call_flag[node_new_bed] = OFF;
    sed->no_displacement_flag[node_new_bed] = OFF;
    if (sed->no_displacement_flag[node1_bed] == ON && sed->no_displacement_flag[node2_bed] == ON) {
        sed->no_displacement_flag[node_new_bed] = ON;
    }

    sed->as_ceiling[node_new_bed] = 0.5*(sed->as_ceiling[node1_bed] + sed->as_ceiling[node2_bed]);
    sed->old_as_ceiling[node_new_bed] = 0.5*(sed->old_as_ceiling[node1_bed] + sed->old_as_ceiling[node2_bed]);
    sed->bed_displacement[node_new_bed] = 0.5*(sed->bed_displacement[node1_bed] + sed->bed_displacement[node2_bed]);
    sed->old_bed_displacement[node_new_bed] = 0.5*(sed->old_bed_displacement[node1_bed] + sed->old_bed_displacement[node2_bed]);
    sed->bed_shear_stress[node_new_bed] = 0.5*(sed->bed_shear_stress[node1_bed] + sed->bed_shear_stress[node2_bed]);
    sed->friction_coef[node_new_bed] = 0.5*(sed->friction_coef[node1_bed] + sed->friction_coef[node2_bed]);
    
    
    sed->bedload_vector[node_new_bed] = svect2d_avg(sed->bedload_vector[node1_bed], sed->bedload_vector[node2_bed]);
    sed->bed_gradient[node_new_bed] = svect2d_avg(sed->bed_gradient[node1_bed], sed->bed_gradient[node2_bed]);
    
    /* allocate and initialize bed and suspended load structs */
    
    sbedload_node_avg(sed->bedload, sed->nsed, node_new_bed, node1_bed, node2_bed);
    
    // average layers ---------------------------------------------------- //
    double thick1 = 0., thick2 = 0., thick = 0., weight1 = 0.5, weight2 = 0.5;
    
    // new active layer
    weight1 = 0.5; weight2 = 0.5;
    thick1 = sed->bed_displacement[node1_bed] - sed->as_ceiling[node1_bed];
    thick2 = sed->bed_displacement[node2_bed] - sed->as_ceiling[node2_bed];
    thick = thick1 + thick2;
    if (thick > SMALL) {
        weight1 = thick1 / thick;
        weight2 = thick2 / thick;
    }
    sbed_layer_node_avg(sed->active_layer, 1, sed->nsed, node_new_bed, node1_bed, node2_bed, weight1, weight2, YES);
    
    // old active layer
    weight1 = 0.5; weight2 = 0.5;
    thick1 = sed->old_bed_displacement[node1_bed] - sed->old_as_ceiling[node1_bed];
    thick2 = sed->old_bed_displacement[node2_bed] - sed->old_as_ceiling[node2_bed];
    thick = thick1 + thick2;
    if (thick > SMALL) {
        weight1 = thick1 / thick;
        weight2 = thick2 / thick;
    }
    sbed_layer_node_avg(sed->old_active_layer, 1, sed->nsed, node_new_bed, node1_bed, node2_bed, weight1, weight2, YES);
    
    // new bed layers
    for (ilayer=0; ilayer<sed->nlayers; ilayer++) {
        weight1 = 0.5; weight2 = 0.5;
        thick1 = sed->bed_layer[ilayer].thickness[node1_bed];
        thick2 = sed->bed_layer[ilayer].thickness[node2_bed];
        thick = thick1 + thick2;
        if (thick > SMALL) {
            weight1 = thick1 / thick;
            weight2 = thick2 / thick;
        }
        sbed_layer_node_avg(&(sed->bed_layer[ilayer]), 1, sed->nsed, node_new_bed, node1_bed, node2_bed, weight1, weight2, NO);
    }
    
    // old bed layers
    for (ilayer=0; ilayer<sed->nlayers; ilayer++) {
        weight1 = 0.5; weight2 = 0.5;
        thick1 = sed->old_bed_layer[ilayer].thickness[node1_bed];
        thick2 = sed->old_bed_layer[ilayer].thickness[node2_bed];
        thick = thick1 + thick2;
        if (thick > SMALL) {
            weight1 = thick1 / thick;
            weight2 = thick2 / thick;
        }
        sbed_layer_node_avg(&(sed->old_bed_layer[ilayer]), 1, sed->nsed, node_new_bed, node1_bed, node2_bed, weight1, weight2, NO);
    }
    // ---------------------------------------------------------------------//
    
    /* allocate consolidation tables */
    for (k = 0; k < 5; k++) {
        for (j = 0; j < sed->ncti; j++) {
           sed->consolidation_bed_properties[k][j][node_new_bed]=0.5*(sed->consolidation_bed_properties[k][j][node1_bed] + sed->consolidation_bed_properties[k][j][node2_bed]);
         }
    }
  }

  sed->susload_vector[node_new] = svect2d_avg(sed->susload_vector[node1], sed->susload_vector[node2]);
  ssusload_node_avg(sed->susload, sed->nsed, node_new, node1, node2);
}

/*----------------------------------------------------------------------------*/
/* renumber sediment struct nodes - for adaption (only for 2d right now!) */
void ssediment_renumber(SSED *sed, int max_nnode, int *new_numbers, int *order_tmp, int *itmp, double *dtmp, SVECT2D *vtmp) {
    node_renumber_int(max_nnode, sed->friction_flag, itmp, new_numbers, order_tmp);
    node_renumber_int(max_nnode, sed->call_flag, itmp, new_numbers, order_tmp);
    node_renumber_int(max_nnode, sed->no_displacement_flag, itmp, new_numbers, order_tmp);
    
    node_renumber_double(max_nnode, sed->as_ceiling, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sed->old_as_ceiling, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sed->bed_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sed->old_bed_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sed->bed_shear_stress, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sed->friction_coef, dtmp, new_numbers, order_tmp);

    node_renumber_vect2d(max_nnode, sed->susload_vector, vtmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, sed->bedload_vector, vtmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, sed->bed_gradient, vtmp, new_numbers, order_tmp);
    
    /* allocate and initialize bed and suspended load structs */
    ssusload_renumber(sed->susload, sed->nsed, max_nnode, new_numbers, order_tmp, dtmp, vtmp);
    sbedload_renumber(sed->bedload, sed->nsed, max_nnode, new_numbers, order_tmp, dtmp, vtmp);
    
    // active and bed layers (sed->bed_layer[ilayer].thickness[inode]
    sbed_layer_renumber(sed->active_layer, 1, sed->nsed, max_nnode, new_numbers, order_tmp, dtmp);
    sbed_layer_renumber(sed->old_active_layer, 1, sed->nsed, max_nnode, new_numbers, order_tmp, dtmp);
    sbed_layer_renumber(sed->bed_layer, sed->nlayers, sed->nsed, max_nnode, new_numbers, order_tmp, dtmp);
    sbed_layer_renumber(sed->old_bed_layer, sed->nlayers, sed->nsed, max_nnode, new_numbers, order_tmp, dtmp);
    
    /* allocate consolidation tables */
    int j, k;
    for (k = 0; k < 5; k++) {
        for (j = 0; j < sed->ncti; j++) {
            node_renumber_double(max_nnode, sed->consolidation_bed_properties[k][j], dtmp, new_numbers, order_tmp);
        }
    }
}


/*----------------------------------------------------------------------------*/
/* renumber sediment struct nodes - for adaption (only for 2d right now!) */
void ssediment_renumber_3d_sus(SSED *sed, int max_nnode, int *new_numbers, int *order_tmp, int *itmp, double *dtmp, SVECT2D *vtmp) {
    
    node_renumber_vect2d(max_nnode, sed->susload_vector, vtmp, new_numbers, order_tmp);

    /* renumber suspended load struct arays allocated on full grid */
    int ised;
    
    for(ised=0;ised<sed->nsed;ised++){
        node_renumber_double(max_nnode, sed->susload[ised].c, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, sed->susload[ised].old_c, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, sed->susload[ised].older_c, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, sed->susload[ised].source, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, sed->susload[ised].sink, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, sed->susload[ised].error, dtmp, new_numbers, order_tmp);
    }
    
    
}


/*----------------------------------------------------------------------------*/
/* renumber sediment struct nodes - for adaption (only for 2d right now!) */
void ssediment_renumber_3d_bed(SSED *sed, int max_nnode, int *new_numbers, int *order_tmp, int *itmp, double *dtmp, SVECT2D *vtmp) {
    node_renumber_int(max_nnode, sed->friction_flag, itmp, new_numbers, order_tmp);
    node_renumber_int(max_nnode, sed->call_flag, itmp, new_numbers, order_tmp);
    node_renumber_int(max_nnode, sed->no_displacement_flag, itmp, new_numbers, order_tmp);

    node_renumber_double(max_nnode, sed->as_ceiling, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sed->old_as_ceiling, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sed->bed_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sed->old_bed_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sed->bed_shear_stress, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sed->friction_coef, dtmp, new_numbers, order_tmp);

    node_renumber_vect2d(max_nnode, sed->bedload_vector, vtmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, sed->bed_gradient, vtmp, new_numbers, order_tmp);

    /* renumber bed load struct and suspended load structs bed node arrays*/
    sbedload_renumber(sed->bedload, sed->nsed, max_nnode, new_numbers, order_tmp, dtmp, vtmp);
    int ised;
    for(ised=0;ised<sed->nsed;ised++) node_renumber_double(max_nnode, sed->susload[ised].rouse_coef, dtmp, new_numbers, order_tmp);

    // active and bed layers (sed->bed_layer[ilayer].thickness[inode]
    sbed_layer_renumber(sed->active_layer, 1, sed->nsed, max_nnode, new_numbers, order_tmp, dtmp);
    sbed_layer_renumber(sed->old_active_layer, 1, sed->nsed, max_nnode, new_numbers, order_tmp, dtmp);
    sbed_layer_renumber(sed->bed_layer, sed->nlayers, sed->nsed, max_nnode, new_numbers, order_tmp, dtmp);
    sbed_layer_renumber(sed->old_bed_layer, sed->nlayers, sed->nsed, max_nnode, new_numbers, order_tmp, dtmp);

    /* allocate consolidation tables */
    int j, k;
    for (k = 0; k < 5; k++) {
        for (j = 0; j < sed->ncti; j++) {
            node_renumber_double(max_nnode, sed->consolidation_bed_properties[k][j], dtmp, new_numbers, order_tmp);
        }
    }
}



/*----------------------------------------------------------------------------*/
/* free ssediment struct */
void ssediment_free(SSED *sed) {

    int i=0, j=0, k=0;
    int inode=0, ilayer=0;
#ifdef _DEBUG
    printf("... freeing sediment and sedlib arrays ...\n");
#endif
    if (sed->consolidation_bed_properties != NULL) {
        for (k = 0; k < 5; k++) {
            for (j = 0; j < sed->ncti; j++) {
                sed->consolidation_bed_properties[k][j] = (double *) tl_free(sizeof(double), sed->nnodes_bed, sed->consolidation_bed_properties[k][j]);
            }
            sed->consolidation_bed_properties[k] = (double **) tl_free(sizeof(double *), sed->ncti, sed->consolidation_bed_properties[k]);
        }
        sed->consolidation_bed_properties = (double ***) tl_free(sizeof(double **), 5, sed->consolidation_bed_properties);
    }

    sbed_layer_free(sed->bed_layer, sed->nlayers, sed->nsed, sed->nnodes_bed); 
    sbed_layer_free(sed->old_bed_layer, sed->nlayers, sed->nsed, sed->nnodes_bed);
    sbed_layer_free(sed->active_layer, ONE, sed->nsed, sed->nnodes_bed);
    sbed_layer_free(sed->old_active_layer, ONE, sed->nsed, sed->nnodes_bed);
    ssusload_free(sed->susload, sed->nsed, sed->nnodes_bed, sed->nnodes_sus);
    sbedload_free(sed->bedload, sed->nsed, sed->nnodes_bed);
    
    sed->grain = (SGRAIN *) tl_free(sizeof(SGRAIN), sed->nsed, sed->grain);
    sed->cohesive_settling_const = (double *) tl_free(sizeof(double), COHSET, sed->cohesive_settling_const);
    sed->wind_wave_const = (double *) tl_free(sizeof(double), WINDWAVE, sed->wind_wave_const);
    sed->noncoh_ent_const = (double *) tl_free(sizeof(double), 1, sed->noncoh_ent_const);
    sed->noncoh_ble_const = (double *) tl_free(sizeof(double), 1, sed->noncoh_ble_const);
    sed->div_coef = (double *) tl_free(sizeof(double), DIVCOEF, sed->div_coef);
    sed->call_flag = (int *) tl_free(sizeof(int), sed->nnodes_bed, sed->call_flag);
    sed->friction_flag = (int *) tl_free(sizeof(int), sed->nnodes_bed, sed->friction_flag);
    sed->friction_coef = (double *) tl_free(sizeof(double), sed->nnodes_bed, sed->friction_coef);
    sed->no_displacement_flag = (int *) tl_free(sizeof(int), sed->nnodes_bed, sed->no_displacement_flag);
    sed->as_ceiling = (double *) tl_free(sizeof(double), sed->nnodes_bed, sed->as_ceiling);
    sed->old_as_ceiling = (double *) tl_free(sizeof(double), sed->nnodes_bed, sed->old_as_ceiling);
    sed->bed_displacement = (double *) tl_free(sizeof(double), sed->nnodes_bed, sed->bed_displacement);
    sed->bed_shear_stress = (double *) tl_free(sizeof(double), sed->nnodes_bed, sed->bed_shear_stress);
    sed->old_bed_displacement = (double *) tl_free(sizeof(double), sed->nnodes_bed, sed->old_bed_displacement);
    sed->bedload_vector = (SVECT2D *) tl_free(sizeof(SVECT2D), sed->nnodes_bed, sed->bedload_vector);
    sed->susload_vector = (SVECT2D *) tl_free(sizeof(SVECT2D), sed->nnodes_sus, sed->susload_vector);
    sed->bed_gradient = (SVECT2D *) tl_free(sizeof(SVECT2D), sed->nnodes_bed, sed->bed_gradient);
    
    // free sedlib struct
    ssedlib_free(sed->sedlib);

    sed = (SSED *) tl_free(sizeof(SSED), 1, sed);

}

/******************************************************************************/
/******************************************************************************/
/* get sediment variables ready for model run */
/* note :: this happens after sediment hotstart read */
void ssediment_prep(SMODEL *mod, SSED *sed, int iald_flag) {
    
    int ie, ilayer=0, i=0, inode=0, ised=0, ielem=UNSET_INT;
    double blsg = 0.;

    double density = 0.;
    if (mod->density < 1.) {
        density = 1000.;
    }
    
    int nelems2d;
    if (mod->grid->ndim == 2) {
        nelems2d = mod->grid->nelems2d;
    } else if (mod->grid->ndim == 3) {
        nelems2d = mod->grid->nelems2d_bed;
    }
    
    // calculate and store nodal friction flags
    int *flag = (int *) tl_alloc(sizeof(int), nelems2d);
    for(ie = 0; ie < nelems2d; ie++) {
        ielem = ie;
        if (mod->grid->ndim==3) {
            ielem = mod->grid->elem2d_bed[ie];
        }
        if (mod->str_values[mod->grid->elem2d[ielem].string].fterms.mng_flag == YES) {
            flag[ie] = 1;
        } else {
            flag[ie] = 0;
        }
//        printf("istring: %d flag: %d\n",mod->grid->elem2d[ie].string, flag[ie]);
    }
    //exit(-1);
    if (mod->grid->ndim==3) {
        elem2dbed_to_node_int(mod->grid, sed->friction_flag, flag, 2, -1);
    } else {
        elem2d_to_node_int(mod->grid, sed->friction_flag, flag, 2, -1);
    }
    flag = (int *) tl_free(sizeof(int), nelems2d, flag);
    
//    for (inode=0; inode<mod->grid->nnodes_bed; inode++) {
//        printf("node: %d   flag: %d\n",inode,sed->friction_flag[inode]);
//    }
//    exit(-1);
    
    // calculate and store nodal friction coefficients
    double *coef = (double *) tl_alloc(sizeof(double), nelems2d);
    for(ie = 0; ie < nelems2d; ie++) {
        ielem = ie;
        if (mod->grid->ndim==3) {
            ielem = mod->grid->elem2d_bed[ie]; // in this case, ie is really nelems2d_bed
        }
        if (mod->str_values[mod->grid->elem2d[ielem].string].fterms.mng_flag == YES) {
            coef[ie] = mod->str_values[mod->grid->elem2d[ielem].string].fterms.manningsn;
        } else {
            coef[ie] = mod->str_values[mod->grid->elem2d[ielem].string].fterms.eqrheight;
        }
//        printf("istring: %d flag: %20.10f\n",mod->grid->elem2d[ielem].string, coef[ie]);
    }
    if (mod->grid->ndim==3) {
        elem2dbed_to_node_double(mod->grid, mod->str_values, sed->friction_coef, coef);
    } else {
        elem2d_to_node_double(mod->grid, mod->str_values, sed->friction_coef, coef);
    }
    coef = (double *) tl_free(sizeof(double), nelems2d, coef);
    
//    for (inode=0; inode<mod->grid->nnodes_bed; inode++) {
//        printf("node: %d   coef: %10.5f\n",inode,sed->friction_coef[inode]);
//    }
//    exit(-1);
    
    // initialize bed layers
    for (ilayer=0; ilayer<sed->nlayers; ilayer++) {
        if (sed->bed_layer[ilayer].cbed_flag == NO) {
            printf("No cohesive properties assigned to bed layer %d\n", i + 1);
            for (inode=0; inode<mod->grid->nnodes_bed; inode++) {
                sed->bed_layer[ilayer].bulk_density[inode] = 0.0;
                sed->bed_layer[ilayer].critical_erosion_shear[inode] = 0.0;
                sed->bed_layer[ilayer].erosion_rate_constant[inode] = 0.0;
                sed->bed_layer[ilayer].erosion_rate_exponent[inode] = 0.0;
                for (ised = 0; ised < sed->nsed; ised++) {
                    if (sed->grain[ised].type == SND) {
                        sed->bed_layer[ilayer].bulk_density[inode] += sed->bed_layer[ilayer].distribution[ised][inode] * mod->density * (sed->grain[ised].porosity + (1. - sed->grain[ised].porosity) * sed->grain[ised].specific_gravity);
                        sed->bed_layer[ilayer].critical_erosion_shear[inode] += SMALL;
                        sed->bed_layer[ilayer].erosion_rate_constant[inode] += SMALL;
                        sed->bed_layer[ilayer].erosion_rate_exponent[inode] += sed->bed_layer[ilayer].distribution[ised][inode] * 1.0;
                    }
                    if (sed->grain[ised].type == CLA) {
                        sed->bed_layer[ilayer].bulk_density[inode] += sed->bed_layer[ilayer].distribution[ised][inode] * mod->density * (sed->grain[ised].porosity + (1. - sed->grain[ised].porosity) * sed->grain[ised].specific_gravity);						
                        sed->bed_layer[ilayer].critical_erosion_shear[inode] += sed->bed_layer[ilayer].distribution[ised][inode] * sed->grain[ised].clay.tau_ce;
                        sed->bed_layer[ilayer].erosion_rate_constant[inode] += sed->bed_layer[ilayer].distribution[ised][inode] * sed->grain[ised].clay.erode_const;
                        sed->bed_layer[ilayer].erosion_rate_exponent[inode] += sed->bed_layer[ilayer].distribution[ised][inode] * 1.0;
                    }
                }
            }
        }
    }


    // set old bed layer to current
    sbed_layer_copy(sed->old_bed_layer, sed->bed_layer, sed->nlayers, sed->nsed, mod->grid->nnodes_bed);
    
    
    // initiate active layer as top bed layer
    int top_layer = sed->nlayers - 1;
    for (inode=0; inode<mod->grid->nnodes_bed; inode++) {
        
        // if the active layer distribution is not specified in the hotstart file, define here
        if (iald_flag == OFF) {
            for (ised=0; ised<sed->nsed; ised++) {
                sed->active_layer->distribution[ised][inode] = sed->bed_layer[top_layer].distribution[ised][inode];
                sed->old_active_layer->distribution[ised][inode] = sed->active_layer->distribution[ised][inode];
            }
        }
        
        blsg = 0.0;
        for (ised=0; ised<sed->nsed; ised++) {
            blsg += sed->active_layer->distribution[ised][inode] * sed->grain[ised].specific_gravity;
        }

        sed->active_layer->porosity[inode] = (sed->bed_layer[top_layer].bulk_density[inode] - blsg * mod->density) / (mod->density * (1. - blsg));
        sed->old_active_layer->porosity[inode] = sed->active_layer->porosity[inode];
        
        sed->active_layer->critical_erosion_shear[inode] = sed->bed_layer[top_layer].critical_erosion_shear[inode];
        sed->old_active_layer->critical_erosion_shear[inode] = sed->active_layer->critical_erosion_shear[inode];
        
        sed->active_layer->erosion_rate_constant[inode] = sed->bed_layer[top_layer].erosion_rate_constant[inode];
        sed->old_active_layer->erosion_rate_constant[inode] = sed->active_layer->erosion_rate_constant[inode];
        
        sed->active_layer->erosion_rate_exponent[inode] = sed->bed_layer[top_layer].erosion_rate_exponent[inode];
        sed->old_active_layer->erosion_rate_exponent[inode] = sed->active_layer->erosion_rate_exponent[inode];

        // initialize active stratum ceiling
        sed->as_ceiling[inode] = sed->bed_displacement[inode];
        sed->old_as_ceiling[inode] = sed->bed_displacement[inode];
        
    }
}


/******************************************************************************/
/******************************************************************************/
/* open sediment output files */
void ssediment_open_output(SMODEL *mod) {

    assert(mod->io); /* should be valid */
    int super = strlen(mod->io->sup.filename); /* whether super file exists */
   
    open_output_file( &(mod->io->fout_bed), "sediment bed file", super);
    print_header_sediment(mod, mod->io->fout_bed.fp, PS_FLAG_BED, UNSET_INT);
    
    open_output_file( &(mod->io->fout_bed_flux), "sediment bedload flux file", super);
    print_header_sediment(mod, mod->io->fout_bed_flux.fp, PS_FLAG_BED_FLUX, UNSET_INT);
    
    open_output_file( &(mod->io->fout_active_layer), "sediment active layer file", super);
    print_header_sediment(mod, mod->io->fout_active_layer.fp, PS_FLAG_ALAYER, UNSET_INT);

    int ilayer = 0;
    for (ilayer=0; ilayer<mod->sed->nlayers; ilayer++) {
        open_output_file( &(mod->io->fout_bed_layer[ilayer]), "sediment bed layer", super);
        print_header_sediment(mod, mod->io->fout_bed_layer[ilayer].fp, PS_FLAG_BED_LAYER, ilayer);
    }

    int ised = 0;
    for (ised=0; ised<mod->sed->nsed; ised++) {    
        open_output_file( &(mod->io->fout_sl_grain[ised]), "sediment suspended grain file", super);
        print_header_sediment(mod, mod->io->fout_sl_grain[ised].fp, PS_FLAG_SL_GRAIN, UNSET_INT);

        open_output_file( &(mod->io->fout_bl_grain[ised]), "sediment bedload grain file", super);
        print_header_sediment(mod, mod->io->fout_bl_grain[ised].fp, PS_FLAG_BL_GRAIN, UNSET_INT);
    }

}

/******************************************************************************/
/******************************************************************************/
/* print time-step */

/******************************************************************************/
/******************************************************************************/
/* print time-step MPI*/
void ssediment_print_ts(SGRID *grid, SIO *io, SSED *sed, double time, double outfact, int ndim, SFILE sup, char *proj_name, int it1, int it2, int **ndata, int my_nnode_max, int *my_nnode_ext, int **ndata_sur, int my_nnode_max_sur, int *my_nnode_ext_sur, int flag) {

    int ip, i;
    double *gdata, *edata;
    int ierr, gsize;
    int global_nnode;
#ifdef _MESSG
    MPI_Status *msg_status= grid->smpi->msg_status;
#endif
     // grains
    int ised=0, inode=0;
    int max = 0, global_max;
    int nnodes_sus_max;
    int nnodes_bed_max;
    int offset, index, max_node;
    int nnodes_bed=my_nnode_ext_sur[grid->smpi->myid];
    int nnodes_sus=my_nnode_ext[grid->smpi->myid];
    int macro_nelems = grid->macro_nelems2d;
    char *mesh_tag;
    int super = strlen(sup.filename);

    if (ndim==3){mesh_tag = "mesh3d";
    }else{mesh_tag = "mesh2d";}

    SFILE fout,fout1,fout2;
    
    if(max < (2 * my_nnode_max * sed->nsed)) max = 2 * my_nnode_max * sed->nsed;
    if(max < (my_nnode_max_sur * sed->nlayers * (5+sed->nsed))) max = my_nnode_max_sur * sed->nlayers * (5+sed->nsed);
    if (ndim==3)macro_nelems = grid->macro_nelems3d;
    /* communicate the suspended load */
    if (ndim == 2){offset = 3;}
    else{offset = 2;}

    if (grid->smpi->myid == 0){
      gdata = (double *) tl_alloc(sizeof(double), grid->macro_nnodes * offset);
      edata = (double *) tl_alloc(sizeof(double), my_nnode_max * offset);
      for (ised = 0; ised < sed->nsed; ised++) {
        // suspended grains and concentration error estimator
        
        for(inode=0; inode < nnodes_sus; inode++) {
          index = grid->node[inode].gid * offset;
          gdata[index] = sed->susload[ised].c[inode] * 1.E+6 * sed->grain[ised].reference_c;
          gdata[index + 1] = sed->susload[ised].error[inode];
          if (ndim == 2) gdata[index + 2] = sed->susload[ised].rouse_coef[inode];
        }
#ifdef _MESSG
        for (ip = 1; ip < grid->smpi->npes; ip++) {
          ierr = MPI_Recv(edata, my_nnode_ext[ip]*offset, MPI_DOUBLE, ip, 999, grid->smpi->ADH_COMM, msg_status);        
          if (ierr != MPI_SUCCESS)
            messg_err(ierr);
          inode =0; 
          for (i = 0; i < my_nnode_ext[ip]; i++){
            index = ndata[ip][i] * offset;
            gdata[index] = edata[inode++];
            gdata[index + 1] = edata[inode++];
            if (ndim  == 2) gdata[index + 2] = edata[inode++];
            
          }
        }
#endif
        for (i=0;i<=flag;i++){
          if(i==0){
            fout=io->fout_sl_grain[ised];
            max_node=grid->orig_macro_nnodes;
          }else{
            max_node=grid->macro_nnodes;
            init_adh_file(&(fout));     // initialize file
            build_filename3(fout.filename, MAXLINE, proj_name, "_sl_grain_", ised+1,".dat-", it1,".", it2); // build the new filename
            open_output_file(&(fout), "adapted suspended grain file", super); // open the file
            fprintf(fout.fp, "DATASET\n");
            pdata(proj_name, "", fout.fp, 0, "Suspended Grain", "", "BEGSCL", mesh_tag, grid->macro_nnodes, macro_nelems); // header
            tc_timeunits(fout.fp, outfact);
          }
          fprintf(fout.fp, "TS 0 %15.8e\n", time);
          for(inode=0; inode < max_node; inode++) {
            index = inode * offset; 
            fprintf(fout.fp, "%16.8e %16.8e", gdata[index], gdata[index +1]);
            if (ndim == 2) {
              fprintf(fout.fp, "%16.8e", gdata[index + 2]);
            }
            fprintf(fout.fp, "\n");
          }
          if(i>0) print_trailer(fout.fp);
        }
      }
      gdata = (double *) tl_free(sizeof(double), grid->orig_macro_nnodes * offset, gdata);
      edata = (double *) tl_free(sizeof(double), my_nnode_max * offset, edata);
    }else {
#ifdef _MESSG
      edata = (double *) tl_alloc(sizeof(double), my_nnode_ext[grid->smpi->myid] * offset);
      for (ised = 0; ised < sed->nsed; ised++) {
        for(inode=0; inode < my_nnode_ext[grid->smpi->myid]; inode++) {
          edata[(inode*offset)] = sed->susload[ised].c[inode] * 1.E+6 * sed->grain[ised].reference_c;
          edata[(inode*offset) + 1] = sed->susload[ised].error[inode];
          if (ndim == 2) edata[(inode*offset) + 2] = sed->susload[ised].rouse_coef[inode];
        }
        ierr = MPI_Send(edata, (my_nnode_ext[grid->smpi->myid] * offset), MPI_DOUBLE, 0, 999, grid->smpi->ADH_COMM);
        if (ierr != MPI_SUCCESS)
          messg_err(ierr);
      }
      edata = (double *) tl_free(sizeof(double), (my_nnode_ext[grid->smpi->myid] * offset), edata);
#endif
    }


    //bed grains (cjt :: this will give problems for now in 3d, since it loops over all 2d elements instead of just bed)
    offset = 2;
    if (grid->smpi->myid == 0){
      gdata = (double *) tl_alloc(sizeof(double), grid->macro_nnodes_bed * offset);
      edata = (double *) tl_alloc(sizeof(double), my_nnode_max_sur * offset);

      for (ised = 0; ised < sed->nsed; ised++) {
        
        for(inode=0; inode < nnodes_bed; inode++) {
          if(ndim==3){ 
            index = grid->node[grid->nodeID_2d_to_3d_bed[inode]].global_bed_id * offset;
          }
          else{
            index = grid->node[inode].gid * offset;
          }
          gdata[index] = sed->bedload[ised].c[inode] * 1.E+6 * sed->grain[ised].reference_c;
          gdata[index + 1] = sed->bedload[ised].error[inode];
        }
#ifdef _MESSG
        for (ip = 1; ip < grid->smpi->npes; ip++) {
          ierr = MPI_Recv(edata, my_nnode_ext_sur[ip] * offset, MPI_DOUBLE, ip, 999, grid->smpi->ADH_COMM, msg_status);
          if (ierr != MPI_SUCCESS)
            messg_err(ierr);
          inode=0;
          for (i = 0; i < my_nnode_ext_sur[ip]; i++){
              index = ndata_sur[ip][i] * offset;
              gdata[index] = edata[inode++];
              gdata[index + 1] = edata[inode++];
          }
        }
#endif
        for (i=0;i<=flag;i++){
          if(i==0){
            fout=io->fout_bl_grain[ised];
            max_node=grid->orig_macro_nnodes_bed;
          }else{
            max_node=grid->macro_nnodes_bed;
            init_adh_file(&(fout));     // initialize file
        		build_filename3(fout.filename, MAXLINE, proj_name, "_bl_grain_", ised+1,".dat-", it1,".", it2); // build the new filename
        		open_output_file(&(fout), "adapted bed grain file", super); // open the file
        		fprintf(fout.fp, "DATASET\n");
        		pdata(proj_name, "", fout.fp, 0, "Bed Grain", "", "BEGSCL", mesh_tag, grid->macro_nnodes, grid->macro_nelems2d); // header
        		tc_timeunits(fout.fp, outfact);
          }
          fprintf(fout.fp, "TS 0 %15.8e\n", time);
          for(inode=0; inode < max_node; inode++) {
            index = inode * offset;
            fprintf(io->fout_bl_grain[ised].fp, "%16.8e %16.8e \n", gdata[index], gdata[index +1]);
          }
          if(i>0) print_trailer(fout.fp);
        }
      }
      gdata = (double *) tl_free(sizeof(double), grid->macro_nnodes_bed * offset, gdata);
      edata = (double *) tl_free(sizeof(double), my_nnode_max_sur * offset, edata);
    }
    else{
#ifdef _MESSG
      edata = (double *) tl_alloc(sizeof(double), (my_nnode_ext_sur[grid->smpi->myid] * offset));
      for (ised = 0; ised < sed->nsed; ised++) {
        for(inode=0; inode < nnodes_bed; inode++) {
          edata[(inode*offset)] = sed->bedload[ised].c[inode] * 1.E+6 * sed->grain[ised].reference_c;
          edata[(inode*offset) + 1] = sed->bedload[ised].error[inode];
        }
        ierr = MPI_Send(edata, (my_nnode_ext_sur[grid->smpi->myid] * offset), MPI_DOUBLE, 0, 999, grid->smpi->ADH_COMM);
        if (ierr != MPI_SUCCESS)
          messg_err(ierr);
      }
      edata = (double *) tl_free(sizeof(double), (my_nnode_ext_sur[grid->smpi->myid] * offset), edata);
#endif
    }
    
    // bed layers
    int ilayer=0;
    offset = 5+sed->nsed;
    for (ilayer = 0; ilayer<sed->nlayers; ilayer++) {
      if (grid->smpi->myid == 0){
        gdata = (double *) tl_alloc(sizeof(double), grid->macro_nnodes_bed * offset);
        edata = (double *) tl_alloc(sizeof(double), my_nnode_max_sur * offset);
        for(inode=0; inode < nnodes_bed; inode++) {
          if(ndim==3){
            index = grid->node[grid->nodeID_2d_to_3d_bed[inode]].global_bed_id * offset;
          }
          else{
            index = grid->node[inode].gid * offset;
          }
        	gdata[index] = sed->bed_layer[ilayer].thickness[inode],
         	gdata[index + 1] = sed->bed_layer[ilayer].bulk_density[inode],
         	gdata[index + 2] = sed->bed_layer[ilayer].critical_erosion_shear[inode],
         	gdata[index + 3] = sed->bed_layer[ilayer].erosion_rate_constant[inode],
         	gdata[index + 4] = sed->bed_layer[ilayer].erosion_rate_exponent[inode];
         	for (ised = 0; ised < sed->nsed; ised++) {
           	gdata[index + 4 + ised + 1] = sed->bed_layer[ilayer].distribution[ised][inode];
         	}
        }
#ifdef _MESSG
        for (ip = 1; ip < grid->smpi->npes; ip++) {
          ierr = MPI_Recv(edata, my_nnode_ext_sur[ip] * offset, MPI_DOUBLE, ip, 999, grid->smpi->ADH_COMM, msg_status);
          if (ierr != MPI_SUCCESS)
            messg_err(ierr);
          inode=0;
          for (i = 0; i < my_nnode_ext_sur[ip]; i++){
            index = ndata[ip][i] * offset;
            gdata[index] = edata[inode++];
            gdata[index + 1] = edata[inode++];
            gdata[index + 2] = edata[inode++];
            gdata[index + 3] = edata[inode++];
            gdata[index + 4] = edata[inode++];
            for (ised = 0; ised < sed->nsed; ised++) {
              gdata[index + 4 + ised + 1] = edata[inode++];
            }
          }
        }
#endif
        for (i=0;i<=flag;i++){
          if(i==0){
            fout=io->fout_bed_layer[ilayer];
            max_node=grid->orig_macro_nnodes_bed;
          }else{
            max_node=grid->macro_nnodes_bed;
        		init_adh_file(&(fout));     // initialize file
        		build_filename3(fout.filename, MAXLINE, proj_name, "_bed_layer_", ilayer+1,".dat-", it1,".", it2); // build the new filename
        		open_output_file(&(fout), "adapted bed layer file", super); // open the file
        		fprintf(fout.fp, "DATASET\n");
        		pdata(proj_name, "", fout.fp, 0, "Bed Layer", "", "BEGSCL", mesh_tag, grid->macro_nnodes_bed, grid->macro_nelems2d); // header
          }
        	fprintf(fout.fp, "TS 0 %15.8e\n", time);
          for(inode=0; inode < max_node; inode++) {
            index = inode * offset;
            fprintf(fout.fp, "%16.8e %16.8e %16.8e %16.8e %16.8e", gdata[index], gdata[index +1], gdata[index+2],gdata[index+3],gdata[index+4]);
            for (ised = 0; ised < sed->nsed; ised++) {
              fprintf(fout.fp, "%16.8e", gdata[index + 4 + ised +1]);
            }
            fprintf(fout.fp, "\n");
          }
          if(i>0) print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer
        }
        gdata = (double *) tl_free(sizeof(double), grid->orig_macro_nnodes_bed * offset, gdata);
        edata = (double *) tl_free(sizeof(double), my_nnode_max_sur * offset, edata);
      }
      else{
#ifdef _MESSG
        edata = (double *) tl_alloc(sizeof(double), (my_nnode_ext_sur[grid->smpi->myid] * offset));
        for (ised = 0; ised < sed->nsed; ised++) {
          for(inode=0; inode < nnodes_bed; inode++) {
            edata[(inode*offset)] = sed->bed_layer[ilayer].thickness[inode],
            edata[(inode*offset) + 1] = sed->bed_layer[ilayer].bulk_density[inode],
            edata[(inode*offset) + 2] = sed->bed_layer[ilayer].critical_erosion_shear[inode],
            edata[(inode*offset) + 3] = sed->bed_layer[ilayer].erosion_rate_constant[inode],
            edata[(inode*offset) + 4] = sed->bed_layer[ilayer].erosion_rate_exponent[inode];
            for (ised = 0; ised < sed->nsed; ised++) {
              edata[(inode*offset) + 4 + ised + 1] = sed->bed_layer[ilayer].distribution[ised][inode];
            }
          }
          ierr = MPI_Send(edata, (my_nnode_ext_sur[grid->smpi->myid] * offset), MPI_DOUBLE, 0, 999, grid->smpi->ADH_COMM);
          if (ierr != MPI_SUCCESS)
            messg_err(ierr);
        }
        edata = (double *) tl_free(sizeof(double), (my_nnode_ext_sur[grid->smpi->myid] * offset), edata);
#endif
      }
    }

    // active layer

    // bed
   
    // bed load flux
       
      if (grid->smpi->myid == 0){
        gdata = (double *) tl_alloc(sizeof(double), grid->macro_nnodes_bed * offset);
        edata = (double *) tl_alloc(sizeof(double), my_nnode_max_sur * offset);
        fprintf(io->fout_active_layer.fp, "TS 0 %15.8e\n", time);
        fprintf(io->fout_bed.fp, "TS 0 %15.8e\n", time);
        fprintf(io->fout_bed_flux.fp, "TS 0 %15.8e\n", time);
        for(inode=0; inode < nnodes_bed; inode++) {
          if(ndim==3){
            index = grid->node[grid->nodeID_2d_to_3d_bed[inode]].global_bed_id * offset;
          }
          else{
            index = grid->node[inode].gid * offset;
          }
          gdata[index] = sed->bed_displacement[inode];
          gdata[index + 1] = sed->bed_shear_stress[inode];
          gdata[index + 2] = sed->bedload_vector[inode].x;
          gdata[index + 3] = sed->bedload_vector[inode].y;
          gdata[index + 4] = sed->active_layer->thickness[inode];
          for (ised = 0; ised < sed->nsed; ised++) {
            gdata[index + 4 + ised + 1] = sed->active_layer->distribution[ised][inode];
          }
        }
#ifdef _MESSG
        for (ip = 1; ip < grid->smpi->npes; ip++) {
          ierr = MPI_Recv(edata, my_nnode_ext_sur[ip] * offset, MPI_DOUBLE, ip, 999, grid->smpi->ADH_COMM, msg_status);
          if (ierr != MPI_SUCCESS)
            messg_err(ierr); 
          inode=0;       
          for (i = 0; i < my_nnode_ext_sur[ip]; i++){
            index = ndata[ip][i] * offset;
            gdata[index] = edata[inode++];
            gdata[index + 1] = edata[inode++];
            gdata[index + 2] = edata[inode++];
            gdata[index + 3] = edata[inode++];
            gdata[index + 4] = edata[inode++];
            for (ised = 0; ised < sed->nsed; ised++) {
              gdata[index + 4 + ised + 1] = edata[inode++];
            }
          }
        }
#endif
       for (i=0;i<=flag;i++){
         if(i==0){
           fout=io->fout_active_layer;
           fout1=io->fout_bed;
           fout2=io->fout_bed_flux;
           max_node=grid->orig_macro_nnodes_bed;
         }else{
           max_node=grid->macro_nnodes_bed;
           init_adh_file(&(fout));     // initialize file
    			 build_filename2(fout.filename, MAXLINE, proj_name, "_active_layer.dat-", it1,".", it2);
    			 open_output_file(&(fout), "adapted active layer file", super); // open the file
    	     fprintf(fout.fp, "DATASET\n");
           pdata(proj_name, "", fout.fp, 0, "Active Layer", "", "BEGSCL", mesh_tag, grid->nnodes_bed, grid->macro_nelems2d); // header
    			 init_adh_file(&(fout1));     // initialize file
    			 build_filename2(fout1.filename, MAXLINE, proj_name, "_bed.dat-", it1,".", it2);
    			 open_output_file(&(fout1), "adapted bed file", super); // open the file
    			 fprintf(fout1.fp, "DATASET\n");
    			 pdata(proj_name, "", fout1.fp, 0, "Bed Properties", "", "BEGSCL", mesh_tag, grid->nnodes_bed, grid->macro_nelems2d); // header
    			 init_adh_file(&(fout2));     // initialize file
    			 build_filename2(fout2.filename, MAXLINE, proj_name, "_bed_flux.dat-", it1,".", it2);
    			 open_output_file(&(fout2), "adapted bed flux file", super); // open the file
    			 fprintf(fout2.fp, "DATASET\n");
    			 pdata(proj_name, "", fout2.fp, 0, "Bed Load Flux", "", "BEGSCL", mesh_tag, grid->initial_nnodes, grid->macro_nelems2d); // header
         }
         fprintf(fout.fp, "TS 0 %15.8e\n", time);
         fprintf(fout1.fp, "TS 0 %15.8e\n", time);
         fprintf(fout2.fp, "TS 0 %15.8e\n", time);
         for(inode=0; inode < max_node; inode++) {
           index = inode * offset;
           fprintf(fout1.fp, "%16.8e %16.8e \n",gdata[index], gdata[index+1]);
           fprintf(fout2.fp, "%16.8e %16.8e \n",gdata[index+2],gdata[index+3]);
           fprintf(fout.fp, "%16.8e", gdata[index+4]);
           for (ised = 0; ised < sed->nsed; ised++) {
             fprintf(fout.fp, "%16.8e",gdata[index+4+ised+1]);
           }
           fprintf(fout.fp, "\n");
         }
         if(i>0){
           print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer
           print_trailer(fout1.fp); // add ENNDS, close the file and NULL the pointer
           print_trailer(fout2.fp); // add ENNDS, close the file and NULL the pointer
         }
       }
        gdata = (double *) tl_free(sizeof(double), grid->macro_nnodes_bed * offset, gdata);
        edata = (double *) tl_free(sizeof(double), my_nnode_max_sur * offset, edata);
      }
      else{
#ifdef _MESSG
        edata = (double *) tl_alloc(sizeof(double), (my_nnode_ext_sur[grid->smpi->myid] * offset));
        for(inode=0; inode < nnodes_bed; inode++) {
          edata[inode*offset] = sed->bed_displacement[inode];
          edata[(inode*offset) + 1] = sed->bed_shear_stress[inode];
          edata[(inode*offset) + 2] = sed->bedload_vector[inode].x;
          edata[(inode*offset) + 3] = sed->bedload_vector[inode].y;
          edata[(inode*offset) + 4] = sed->active_layer->thickness[inode];
          for (ised = 0; ised < sed->nsed; ised++) {
            edata[(inode*offset) + 4 + ised + 1] = sed->active_layer->distribution[ised][inode];
          }
        }
        ierr = MPI_Send(edata, (my_nnode_ext_sur[grid->smpi->myid] * offset), MPI_DOUBLE, 0, 999, grid->smpi->ADH_COMM);
        if (ierr != MPI_SUCCESS)
          messg_err(ierr);
        
        edata = (double *) tl_free(sizeof(double), (my_nnode_ext_sur[grid->smpi->myid] * offset), edata);
#endif
      }
     
}

/******************************************************************************/
/******************************************************************************/
/* print time-step */
void ssediment_print_adapted_ts(SGRID *grid, SSED *sed, double time, double outfact, int ndim, SFILE sup, char *proj_name, int it1, int it2) {
    
    int super = strlen(sup.filename);
    char *mesh_tag;
    
    int nnodes_sus = grid->initial_nnodes;      // only print the initial nnodes (not nodes added via adaption)
    int nnodes_bed = grid->initial_nnodes_bed;  // only print the initial nnodes (not nodes added via adaption)
    
    // convert time to specified output units
    double c_time = time * outfact;
    
    int nelems = 0;
    if (ndim == 2) {
        mesh_tag = "mesh2d";
        nelems = grid->nelems2d;
    } else if (ndim == 3) {
        mesh_tag = "mesh3d";
        nelems = grid->nelems3d;
    } else {
        return;
    }
    
    SFILE fout;
    
    // grains
    int ised=0, inode=0;
    for (ised = 0; ised < sed->nsed; ised++) {
        
        // suspended grains
        init_adh_file(&(fout));     // initialize file
        build_filename3(fout.filename, MAXLINE, proj_name, "_sl_grain_", ised+1,".dat-", it1,".", it2); // build the new filename
        open_output_file(&(fout), "adapted suspended grain file", super); // open the file
        fprintf(fout.fp, "DATASET\n");
        pdata(proj_name, "", fout.fp, 0, "Suspended Grain", "", "BEGSCL", mesh_tag, grid->initial_nnodes, grid->initial_nelems); // header
        tc_timeunits(fout.fp, outfact);
        fprintf(fout.fp, "TS 0 %15.8e\n", time);
        for(inode=0; inode < nnodes_sus; inode++) {
            //fprintf(fout.fp, "%16.8e %16.8e", sed->susload[ised].c[inode], sed->susload[ised].error[inode]);
            fprintf(fout.fp, "%16.8e %16.8e", sed->susload[ised].c[inode] * 1.E+6 * sed->grain[ised].reference_c, sed->susload[ised].error[inode]);
            if (ndim == 2) {
                fprintf(fout.fp, "%16.8e", sed->susload[ised].rouse_coef[inode]);
            }
            fprintf(fout.fp, "\n");
        }
        print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer
        
        //bed grains (cjt :: this will give problems for now in 3d, since it loops over all 2d elements instead of just bed)
        init_adh_file(&(fout));     // initialize file
        build_filename3(fout.filename, MAXLINE, proj_name, "_bl_grain_", ised+1,".dat-", it1,".", it2); // build the new filename
        open_output_file(&(fout), "adapted bed grain file", super); // open the file
        fprintf(fout.fp, "DATASET\n");
        pdata(proj_name, "", fout.fp, 0, "Bed Grain", "", "BEGSCL", mesh_tag, grid->initial_nnodes, grid->initial_nelems); // header
        tc_timeunits(fout.fp, outfact);
        fprintf(fout.fp, "TS 0 %15.8e\n", time);
        for(inode=0; inode < nnodes_bed; inode++) {
            //fprintf(fout.fp, "%16.8e %16.8e \n", sed->bedload[ised].c[inode], sed->bedload[ised].error[inode]);
            fprintf(fout.fp, "%16.8e %16.8e \n", sed->bedload[ised].c[inode] * 1.E+6 * sed->grain[ised].reference_c, sed->bedload[ised].error[inode]);
        }
        print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer
        
    }
    
    // bed layers
    int ilayer=0;
    for (ilayer = 0; ilayer<sed->nlayers; ilayer++) {
        init_adh_file(&(fout));     // initialize file
        build_filename3(fout.filename, MAXLINE, proj_name, "_bed_layer_", ilayer+1,".dat-", it1,".", it2); // build the new filename
        open_output_file(&(fout), "adapted bed layer file", super); // open the file
        fprintf(fout.fp, "DATASET\n");
        pdata(proj_name, "", fout.fp, 0, "Bed Layer", "", "BEGSCL", mesh_tag, grid->initial_nnodes, grid->initial_nelems); // header
        fprintf(fout.fp, "TS 0 %15.8e\n", time);
        for(inode=0; inode < nnodes_bed; inode++) {
            fprintf(fout.fp, "%16.8e %16.8e %16.8e %16.8e %16.8e",
                    sed->bed_layer[ilayer].thickness[inode],
                    sed->bed_layer[ilayer].bulk_density[inode],
                    sed->bed_layer[ilayer].critical_erosion_shear[inode],
                    sed->bed_layer[ilayer].erosion_rate_constant[inode],
                    sed->bed_layer[ilayer].erosion_rate_exponent[inode]
                    );
            for (ised = 0; ised < sed->nsed; ised++) {
                fprintf(fout.fp, "%16.8e", sed->bed_layer[ilayer].distribution[ised][inode]);
            }
            fprintf(fout.fp, "\n");
        }
        print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer
    }
    
    // active layer
    init_adh_file(&(fout));     // initialize file
    build_filename2(fout.filename, MAXLINE, proj_name, "_active_layer.dat-", it1,".", it2);
    open_output_file(&(fout), "adapted active layer file", super); // open the file
    fprintf(fout.fp, "DATASET\n");
    pdata(proj_name, "", fout.fp, 0, "Active Layer", "", "BEGSCL", mesh_tag, grid->initial_nnodes, grid->initial_nelems); // header
    fprintf(fout.fp, "TS 0 %15.8e\n", time);
    for(inode=0; inode < nnodes_bed; inode++) {
        fprintf(fout.fp, "%16.8e", sed->active_layer->thickness[inode]);
        for (ised = 0; ised < sed->nsed; ised++) {
            fprintf(fout.fp, "%16.8e", sed->active_layer->distribution[ised][inode]);
        }
        fprintf(fout.fp, "\n");
    }
    print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer
    
    // bed
    init_adh_file(&(fout));     // initialize file
    build_filename2(fout.filename, MAXLINE, proj_name, "_bed.dat-", it1,".", it2);
    open_output_file(&(fout), "adapted bed file", super); // open the file
    fprintf(fout.fp, "DATASET\n");
    pdata(proj_name, "", fout.fp, 0, "Bed Properties", "", "BEGSCL", mesh_tag, grid->initial_nnodes, grid->initial_nelems); // header
    fprintf(fout.fp, "TS 0 %15.8e\n", time);
    for(inode=0; inode < nnodes_bed; inode++) {
        fprintf(fout.fp, "%16.8e %16.8e \n",sed->bed_displacement[inode], sed->bed_shear_stress[inode]);
    }
    print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer
    
    // bed load flux
    init_adh_file(&(fout));     // initialize file
    build_filename2(fout.filename, MAXLINE, proj_name, "_bed_flux.dat-", it1,".", it2);
    open_output_file(&(fout), "adapted bed flux file", super); // open the file
    fprintf(fout.fp, "DATASET\n");
    pdata(proj_name, "", fout.fp, 0, "Bed Load Flux", "", "BEGSCL", mesh_tag, grid->initial_nnodes, grid->initial_nelems); // header
    fprintf(fout.fp, "TS 0 %15.8e\n", time);
    for(inode=0; inode < nnodes_bed; inode++) {
        fprintf(fout.fp, "%16.8e %16.8e \n", sed->bedload_vector[inode].x, sed->bedload_vector[inode].y);
    }
    print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer
    
}

/******************************************************************************/
/******************************************************************************/
/* open sediment hotstart files */
void ssediment_read_hot(SMODEL *mod) {

    char line[MAXLINE];           /* the input line */
    char name[MAXLINE];           /* the name of the current data set */
    char msg[MAXLINE];            /* for reporting to screen */
    char *data = NULL;
    int i, j, k, ised = 0, inode = 0, jnode = 0, ibl = 0;     /* loop counter */
    double *tmp_data_set;         /* temporary space to store the current data set in */
    
    int iald_flag = OFF;

    double c_inv = 1.0;

    // aliases
    SGRID *grid = mod->grid;
    SIO io = *(mod->io);
    SSED *sed = mod->sed;
    
    assert(mod);
    assert(&io);
    
    // these flag determine how to initialize arrays based on what is read in
    int hot_flag_sed[mod->sed->nsed]; sarray_int_init(hot_flag_sed, mod->sed->nsed);
    int hot_flag_psed[mod->sed->nsed]; sarray_int_init(hot_flag_psed, mod->sed->nsed);
    
    printf("-- reading shallow water hotstart file for sediment initiation: %s\n",io.hot.filename);
    
    //-------------------------------------------------------------------------------------------------//
    //-------------------------------------------------------------------------------------------------//
    // Read full grid variables first
    //-------------------------------------------------------------------------------------------------//
    /* allocates space to read the current data set */
    tmp_data_set = (double *) tl_alloc(sizeof(double), 3 * grid->macro_nnodes);
    
    /* loops over the lines in the input hotstart file looking for the file type */
    while ((fgets(line, MAXLINE, io.hot.fp) != NULL) && (strncmp(line, "DATASET", 7) != AGREE));
    
    /* reads the data sets while there are data sets */
    while (read_data_set(io, &io.hot, grid, tmp_data_set, line) == YES) {
        
        /* parse the name from the name line */
        io_save_line(&io, io.hot.fp, io.hot.filename, line);
        if (parse_card(line, &data) != CARD_NAME) {
            io_read_error(io, "Expected to read data set name.", TRUE);
        }
        read_text_field_custom(io, &data, name, MAXLINE, NULL, "data set name", 1, TRUE);
        
        /* convert name to uppercase to simplify comparisons */
        convert_to_uppercase(name);
        
        /* read sediment concentration initial conditions */
        if ((strncmp("ISED", name, 4) == AGREE)  || (strncmp("SUSPENDED SEDIMENT CONCENTRATION", name, 32) == AGREE)) {
            
            /* get the constituent number from the name line */
            ised = read_int_field_custom(io, &data, NULL, "sediment grain ID", 2, TRUE) - 1;
            if (ised < 0 || ised >= sed->nsed) {
                io_read_error(io, "Tried to read concentration for non existent sediment " "constituent.", TRUE);
            }
            hot_flag_sed[ised] = ON;
            
            sprintf(msg, " Reading initial sediment concentrations (grain " "ID: %d).", ised + 1);
            root_print(msg);
            
            /* if it is sediment it is read in in micromass per unit mass and we need it mass per unit mass */
            c_inv = 1.e-6 / sed->grain[ised].reference_c;
            for (i = 0; i < grid->nnodes; i++) {
                sed->susload[ised].c[i] = tmp_data_set[grid->node[i].gid] * c_inv;
            }
        }
        
        else if ((strncmp("ISEDP", name, 4) == AGREE)  || (strncmp("OLD SUSPENDED SEDIMENT CONCENTRATION", name, 36) == AGREE)) {
            
            /* get the constituent number from the name line */
            ised = read_int_field_custom(io, &data, NULL, "sediment grain ID", 2, TRUE) - 1;
            if (ised < 0 || ised >= sed->nsed) {
                io_read_error(io, "Tried to read concentration for non existent sediment " "constituent.", TRUE);
            }
            
            hot_flag_psed[ised] = ON;
            
            sprintf(msg, " Reading previous sediment concentrations (grain " "ID: %d).", ised + 1);
            root_print(msg);
            
            /* if it is sediment it is read in in micromass per unit mass and we need it mass per unit mass */
            c_inv = 1.e-6 / sed->grain[ised].reference_c;
            for (i = 0; i < grid->nnodes; i++) {
                sed->susload[ised].old_c[i] = tmp_data_set[grid->node[i].gid] * c_inv;
            }
        }
    }
    io_save_line(&io, NULL, "", "");
    rewind(io.hot.fp);

    //-------------------------------------------------------------------------------------------------//
    /* force the initial condition to meet Dirichlet boundary conditions */
    root_print("------ (WARNING) Sediment Dirichlet boundary conditions are being applied to nodes (superseding initial conditions).");
    int istring = -1, isers = -1;
    for (i = 0; i < grid->nnodes; i++) {
        istring = grid->node[i].string;
        if (istring > NORMAL) {
            for (ised = 0; ised < mod->sed->nsed; ised++) {
                if (mod->str_values[istring].sed[ised].bc_flag == BCT_DIR) {
                        isers = mod->str_values[istring].sed[ised].iu_0;
                        sed->susload[ised].c[i] = sseries_get_value(isers, mod->series_head,0) * 1.E-6;
                }
            }
        }
    }

    //-------------------------------------------------------------------------------------------------//
    /* initialize arrays */
    for (ised=0; ised<mod->sed->nsed; ised++) {
        if (hot_flag_psed[ised] == OFF) {
            sarray_copy_dbl(mod->sed->susload[ised].old_c, sed->susload[ised].c, grid->nnodes);
        } else {
        }
    }
    

    //-------------------------------------------------------------------------------------------------//
    //-------------------------------------------------------------------------------------------------//
    // Read bed (2d) variables
    //-------------------------------------------------------------------------------------------------//

    /* loops over the lines in the input hotstart file looking for the file type */
    while ((fgets(line, MAXLINE, io.hot.fp) != NULL) && (strncmp(line, "DATASET", 7) != AGREE));

    /* reads the data sets while there are data sets */
    while (read_data_set(io, &io.hot, grid, tmp_data_set, line) == YES) {
        
        /* parse the name from the name line */
        io_save_line(&io, io.hot.fp, io.hot.filename, line);
        if (parse_card(line, &data) != CARD_NAME) {
            io_read_error(io, "Expected to read data set name.", TRUE);
        }
        read_text_field_custom(io, &data, name, MAXLINE, NULL, "data set name", 1, TRUE);
        
        /* convert name to uppercase to simplify comparisons */
        convert_to_uppercase(name);
        
        if ((strncmp("IBD", name, 3) == AGREE) || (strncmp("BED DISPLACEMENT", name, 16) == AGREE)) {
            root_print(" Reading initial bed displacements.");
            for (i = 0; i < grid->nnodes_bed; i++) {
                if(grid->ndim==3){ j = grid->node[grid->nodeID_2d_to_3d_bed[i]].global_surf_id; }
                else { j = grid->node[i].global_surf_id; } 
                sed->bed_displacement[i] = tmp_data_set[j];
                sed->old_bed_displacement[i] = sed->bed_displacement[i];
            }
        }
        else if ((strncmp("IPBD", name, 4) == AGREE) || (strncmp("OLD BED DISPLACEMENT", name, 20) == AGREE)) {
            root_print(" Reading previous bed displacements.");
            for (i = 0; i < grid->nnodes_bed; i++) {
              if(grid->ndim==3){ j = grid->node[grid->nodeID_2d_to_3d_bed[i]].global_surf_id; }
              else { j = grid->node[i].global_surf_id; }
              sed->old_bed_displacement[i] = tmp_data_set[j];
            }
        }
        else if ((strncmp("IBLT", name, 4) == AGREE) || (strncmp("BED LAYER THICKNESS", name, 19) == AGREE)) {
            root_print(" Reading previous bed layer thickness.");
            ibl = read_int_field_custom(io, &data, NULL, "bed layer ID", 2, TRUE) - 1;
            if (ibl < 0 || ibl >= sed->nlayers) {
                io_read_error(io, "Tried to read thickness for non existent bed layer.", TRUE);
            }
            sprintf(msg, " Reading initial bed layer thicknesses (bed layer ID: %d).", ibl + 1);
            root_print(msg);
            for (i = 0; i < grid->nnodes_bed; i++) {
              if(grid->ndim==3){ j = grid->node[grid->nodeID_2d_to_3d_bed[i]].global_surf_id; }
              else { j = grid->node[i].global_surf_id; }
              sed->bed_layer[ibl].thickness[i] = tmp_data_set[j];
            }
        }
        else if ((strncmp("IBLD", name, 4) == AGREE)  || (strncmp("BED LAYER DISTRIBUTION", name, 22) == AGREE)) {
            root_print(" Reading previous bed layer distribution.");
            ibl = read_int_field_custom(io, &data, NULL, "bed layer ID", 2, TRUE) - 1;
            if (ibl < 0 || ibl >= sed->nlayers) {
                io_read_error(io, "Tried to read distribution for non existent bed layer.", TRUE);
            }
            sprintf(msg, " Reading initial bed layer distributions (bed layer ID: %d).", ibl + 1);
            root_print(msg);
            for (i = 0; i < grid->nnodes_bed; i++) {
                for (k = 0; k < sed->nsed; k++) {
                  if(grid->ndim==3){ j = grid->node[grid->nodeID_2d_to_3d_bed[i]].global_surf_id + k; }
                  else { j = grid->node[i].global_surf_id + k; }
                  sed->bed_layer[ibl].distribution[k][i] = tmp_data_set[j];
                }
            }
        }
        else if ((strncmp("ICBP", name, 4) == AGREE)  || (strncmp("COHESIVE BED PROPERTIES", name, 23) == AGREE)) {
            root_print(" Reading previous cohesive bed properties.");
            ibl = read_int_field_custom(io, &data, NULL, "bed layer ID", 2, TRUE) - 1;
            if (ibl < 0 || ibl >= sed->nlayers) {
                io_read_error(io, "Tried to read properties for non existent bed layer.", TRUE);
            }
            sprintf(msg, " Reading initial cohesive bed layer properties (bed layer " "ID: %d).", ibl + 1);
            root_print(msg);
            for (i = 0; i < grid->nnodes_bed; i++) {
              if(grid->ndim==3){ j = grid->node[grid->nodeID_2d_to_3d_bed[i]].global_surf_id; }
              else { j = grid->node[i].global_surf_id; }
              sed->bed_layer[ibl].bulk_density[i] = tmp_data_set[j];
              sed->bed_layer[ibl].critical_erosion_shear[i] = tmp_data_set[j + 1];
              sed->bed_layer[ibl].erosion_rate_constant[i] = tmp_data_set[j + 2];
              sed->bed_layer[ibl].erosion_rate_exponent[i] = tmp_data_set[j + 3];
            }
        }
        else if ((strncmp("IALT", name, 4) == AGREE) || (strncmp("ACTIVE LAYTER THICKNESS", name, 23) == AGREE)) {
            root_print(" Reading initial active layer thicknesses.");
            for (i = 0; i < grid->nnodes_bed; i++) {
              if(grid->ndim==3){ j = grid->node[grid->nodeID_2d_to_3d_bed[i]].global_surf_id; }
              else { j = grid->node[i].global_surf_id; } 
              sed->active_layer->thickness[i] = tmp_data_set[j];
            }
        }
        else if ((strncmp("IALD", name, 4) == AGREE)  || (strncmp("ACTIVE LAYER DISTRIBUTION", name, 25) == AGREE)) {
            root_print(" Reading initial active layer distributions.");
            iald_flag = ON;
            for (i = 0; i < grid->nnodes_bed; i++) {
                for (k = 0, j = 0; k < sed->nsed; k++, j++) {
                  if(grid->ndim==3){ j = grid->node[grid->nodeID_2d_to_3d_bed[i]].global_surf_id + k; }
                  else { j = grid->node[i].global_surf_id + k; }
                  sed->active_layer->distribution[k][i] = tmp_data_set[j];
                  sed->old_active_layer->distribution[k][i] = tmp_data_set[j];
                }
            }
        }
    }
    io_save_line(&io, NULL, "", "");
    rewind(io.hot.fp);
    
    /* clean up memory */
    tmp_data_set = (double *) tl_free(sizeof(double), 3 * grid->macro_nnodes, tmp_data_set);

    //***************************************************************************************************//
    //***************************************************************************************************//
    // prep all sediment variables
    ssediment_prep(mod, mod->sed, iald_flag);
    //***************************************************************************************************//
    //***************************************************************************************************//
}

/******************************************************************************/
/******************************************************************************/
/* calculate error estimator */
void ssediment_calculate_2d_error(SSED *sed, SSW_2D *sw, SGRID *grid, SMAT *mat, double dt) {

    int ised = UNSET_INT;
    int nelems2d = grid->nelems2d; // alias
    
    /* calculates the error indicator for the flow */
    int ie = UNSET_INT, i = UNSET_INT;
    SVECT2D elem_vel[NDONTRI]; svect2d_init_array(elem_vel, NDONTRI);
    double elem_head[NDONTRI]; sarray_dbl_init(elem_head, NDONTRI);
    double elem_c[NDONTRI]; sarray_dbl_init(elem_c, NDONTRI);
    double elem_old_c[NDONTRI]; sarray_dbl_init(elem_old_c, NDONTRI);
    double error = 0., dcdx = 0., dcdy = 0., dudx = 0., dvdy = 0.;
    double coef_0 = 0., coef_1 = 0., coef_2 = 0.;
    
    sarray_int_init(sw->iarray, grid->nnodes); // element per node count for averaging
    for (ised = 0; ised < sed->nsed; ised++) {
        sarray_dbl_init(sed->susload[ised].error, grid->nnodes); // suspended load nodal error
        sarray_dbl_init(sed->bedload[ised].error, grid->nnodes); // bed load nodal error
    }
    
    for (ie = 0; ie < grid->nelems2d; ie++) {
        
        ELEM2D_GET_LOCAL(sw->head, elem_head, grid->elem2d[ie].nodes);
        ELEM2D_GET_LOCAL_VECT2D(sw->vel, elem_vel, grid->elem2d[ie].nodes);
        dudx = 0.; dvdy = 0.;
        for (i=0; i<NDONTRI; i++) {
            sw->iarray[grid->elem2d[ie].nodes[i]] += 1;
            dudx += grid->elem2d[ie].grad_shp[i].x * elem_vel[i].x;
            dvdy += grid->elem2d[ie].grad_shp[i].y * elem_vel[i].y;
        }
        
        // sediment suspended load
        for (ised = 0; ised < sed->nsed; ised++) {
            error = 0;
            ELEM2D_GET_LOCAL(sed->susload[ised].c, elem_c, grid->elem2d[ie].nodes);
            ELEM2D_GET_LOCAL(sed->susload[ised].old_c, elem_old_c, grid->elem2d[ie].nodes);
            
            dcdx = 0.; dcdy = 0.;
            for (i=0; i<NDONTRI; i++) {
                dcdx += grid->elem2d[ie].grad_shp[i].x * elem_c[i];
                dcdy += grid->elem2d[ie].grad_shp[i].y * elem_c[i];
            }
            
            coef_0 = elem_head[0] * (elem_c[0] - elem_old_c[0]) / dt + elem_vel[0].x * elem_head[0] * dcdx + elem_head[0] * elem_c[0] * dudx + elem_head[0] * elem_vel[0].y * dcdy + elem_head[0] * elem_c[0] * dvdy;
            coef_1 = elem_head[1] * (elem_c[1] - elem_old_c[1]) / dt + elem_vel[1].x * elem_head[1] * dcdx + elem_head[1] * elem_c[1] * dudx + elem_head[1] * elem_vel[1].y * dcdy + elem_head[1] * elem_c[1] * dvdy;
            coef_2 = elem_head[2] * (elem_c[2] - elem_old_c[2]) / dt + elem_vel[2].x * elem_head[2] * dcdx + elem_head[2] * elem_c[2] * dudx + elem_head[2] * elem_vel[2].y * dcdy + elem_head[2] * elem_c[2] * dvdy;
            
            /* GLB 0213 changed to subtract out the mean.  This is to account for sources and sinks (such as sediment erosion and deposition)
             a finite mean should be a good approximaiton of the source/sink term (the mean should be near zero for a conservative solution)
             if we don't account for this, we will resolve on sources and sinks rather than errors in the solution */
            error = (1./3.) * (coef_0 + coef_1 + coef_2);
            coef_0 -= error;
            coef_1 -= error;
            coef_2 -= error;
            
            coef_0 = coef_0 * coef_0;
            coef_1 = coef_1 * coef_1;
            coef_2 = coef_2 * coef_2;
            error = sqrt(coef_0 + coef_1 + coef_2) * grid->elem2d[ie].djac;
            
            // nodal error calculation
            for (i=0; i<NDONTRI; i++) {
                sed->susload[ised].error[grid->elem2d[ie].nodes[i]] += error;
            }
            
            /* scales the tolerance */
            error /= mat[grid->elem2d[ie].mat].sed[ised].refine_tolerance;
            
            /* resets the error if larger */
            if (error > grid->elem_error[ie]) {
                grid->elem_error[ie] = error;
            }
        }
        
        // sediment bed load
        for (ised = 0; ised < sed->nsed; ised++) {
            error = 0;
            ELEM2D_GET_LOCAL(sed->bedload[ised].c, elem_c, grid->elem2d[ie].nodes);
            ELEM2D_GET_LOCAL(sed->bedload[ised].old_c, elem_old_c, grid->elem2d[ie].nodes);
            
            dcdx = 0.; dcdy = 0.;
            for (i=0; i<NDONTRI; i++) {
                dcdx += grid->elem2d[ie].grad_shp[i].x * elem_c[i];
                dcdy += grid->elem2d[ie].grad_shp[i].y * elem_c[i];
            }
            
            coef_0 = (elem_c[0] - elem_old_c[0]) / dt + elem_vel[0].x * dcdx + elem_vel[0].y * dcdy;
            coef_1 = (elem_c[1] - elem_old_c[1]) / dt + elem_vel[1].x * dcdx + elem_vel[1].y * dcdy;
            coef_2 = (elem_c[2] - elem_old_c[2]) / dt + elem_vel[2].x * dcdx + elem_vel[2].y * dcdy;
            
            /* GLB 0213 changed to subtract out the mean.  This is to account for sources and sinks (such as sediment erosion and deposition)
             a finite mean should be a good approximaiton of the source/sink term (the mean should be near zero for a conservative solution)
             if we don't account for this, we will resolve on sources and sinks rather than errors in the solution */
            error = (1./3.) * (coef_0 + coef_1 + coef_2);
            coef_0 -= error;
            coef_1 -= error;
            coef_2 -= error;
            
            coef_0 = coef_0 * coef_0;
            coef_1 = coef_1 * coef_1;
            coef_2 = coef_2 * coef_2;
            error = sqrt(coef_0 + coef_1 + coef_2) * grid->elem2d[ie].djac;
            
            // nodal error calculation
            for (i=0; i<NDONTRI; i++) {
                sed->bedload[ised].error[grid->elem2d[ie].nodes[i]] += error;
            }
            
            /* scales the tolerance */
            error /= mat[grid->elem2d[ie].mat].sed[ised].refine_tolerance;
            
            /* resets the error if larger */
            if (error > grid->elem_error[ie]) {
                grid->elem_error[ie] = error;
            }
        }
    }
    // average the elemental to nodal errors for output
    for (i=0; i<grid->nnodes; i++) {
        for (ised = 0; ised < sed->nsed; ised++) {
            sed->susload[ised].error[i] /= sw->iarray[i];
            sed->bedload[ised].error[i] /= sw->iarray[i];
            
        }
    }
    sarray_int_init(sw->iarray, grid->nnodes); // reset initialize utility array
    
}

/******************************************************************************/
/******************************************************************************/
/* calculate error estimator */
void ssediment_calculate_3d_error(SSED *sed, SSW_3D *sw, SGRID *grid, SMAT *mat, double dt)  {
    
    int ised = UNSET_INT;
    int nelems3d = grid->nelems3d; // alias
    
    /* calculates the error indicator for the flow */
    int ie = UNSET_INT, i = UNSET_INT;
    SVECT elem_vel[NDONTET]; svect_init_array(elem_vel, NDONTET);
    double elem_c[NDONTET]; sarray_dbl_init(elem_c, NDONTET);
    double elem_old_c[NDONTET]; sarray_dbl_init(elem_old_c, NDONTET);
    double error = 0., dcdx = 0., dcdy = 0., dcdz = 0.;
    double coef_0 = 0., coef_1 = 0., coef_2 = 0., coef_3 = 0.;
    
    // sediment suspended load
    sarray_int_init(sw->iarray, grid->nnodes); // element per node count for averaging
    for (ised = 0; ised < sed->nsed; ised++) {
        sarray_dbl_init(sed->susload[ised].error, grid->nnodes); // suspended load nodal error
    }
    for (ie = 0; ie < grid->nelems3d; ie++) {
        
        for (i=0; i<NDONTET; i++) {
            sw->iarray[grid->elem3d[ie].nodes[i]] += 1;
        }
        
        ELEM3D_GET_LOCAL_VECT(sw->vel, elem_vel, grid->elem3d[ie].nodes);
        
        for (ised = 0; ised < sed->nsed; ised++) {
            error = 0;
            ELEM3D_GET_LOCAL(sed->susload[ised].c, elem_c, grid->elem3d[ie].nodes);
            ELEM3D_GET_LOCAL(sed->susload[ised].old_c, elem_old_c, grid->elem3d[ie].nodes);
            
            dcdx = 0.; dcdy = 0.; dcdz = 0.;
            for (i=0; i<NDONTET; i++) {
                dcdx += grid->elem3d[ie].grad_shp[i].x * elem_c[i];
                dcdy += grid->elem3d[ie].grad_shp[i].y * elem_c[i];
                dcdz += grid->elem3d[ie].grad_shp[i].z * elem_c[i];
            }
            coef_0 = (elem_c[0] - elem_old_c[0]) / dt + elem_vel[0].x * dcdx + elem_vel[0].y * dcdy + elem_vel[0].z * dcdz;
            coef_1 = (elem_c[1] - elem_old_c[1]) / dt + elem_vel[1].x * dcdx + elem_vel[1].y * dcdy + elem_vel[1].z * dcdz;
            coef_2 = (elem_c[2] - elem_old_c[2]) / dt + elem_vel[2].x * dcdx + elem_vel[2].y * dcdy + elem_vel[2].z * dcdz;
            coef_3 = (elem_c[3] - elem_old_c[3]) / dt + elem_vel[3].x * dcdx + elem_vel[3].y * dcdy + elem_vel[3].z * dcdz;
            
            /* GLB 0213 changed to subtract out the mean.  This is to account for sources and sinks (such as sediment erosion and deposition)
             a finite mean should be a good approximaiton of the source/sink term (the mean should be near zero for a conservative solution)
             if we don't account for this, we will resolve on sources and sinks rather than errors in the solution */
            error = (1./4.) * (coef_0 + coef_1 + coef_2 + coef_3);
            coef_0 -= error;
            coef_1 -= error;
            coef_2 -= error;
            coef_3 -= error;
            
            coef_0 = coef_0 * coef_0;
            coef_1 = coef_1 * coef_1;
            coef_2 = coef_2 * coef_2;
            coef_3 = coef_3 * coef_3;
            error = sqrt(coef_0 + coef_1 + coef_2 + coef_3) * grid->elem3d[ie].djac;

            // nodal error calculation
            for (i=0; i<NDONTET; i++) {
                sed->susload[ised].error[grid->elem3d[ie].nodes[i]] += error;
            }
            
            /* scales the tolerance */
            error /= mat[grid->elem3d[ie].mat].sed[ised].refine_tolerance;
            
            /* resets the error if larger */
            if (error > grid->elem_error[ie]) {
                grid->elem_error[ie] = error;
            }
        }
    }
    // average the elemental to nodal errors for output
    for (i=0; i<grid->nnodes; i++) {
        for (ised = 0; ised < sed->nsed; ised++) {
            sed->susload[ised].error[i] /= sw->iarray[i];
            
        }
    }
    sarray_int_init(sw->iarray, grid->nnodes); // reset initialize utility array
    
//    // sediment bed load // this is bugged, it needs to be over just bed 2d elems (cjt)
//    for (ie = 0; ie < sed->nelems_bed; ie++) {
//        
//        ELEM2D_GET_LOCAL(sw->head, elem_head, grid->elem2d[ie].nodes);
//        ELEM2D_GET_LOCAL_VECT2D(sw->vel, elem_vel, grid->elem2d[ie].nodes);
//    
//        // sediment bed load
//        for (ised = 0; ised < sed->nsed; ised++) {
//            error = 0;
//            sed->bedload[ised].elem_error[ie] = 0;
//            ELEM2D_GET_LOCAL(sed->bedload[ised].c, elem_c, grid->elem2d[ie].nodes);
//            ELEM2D_GET_LOCAL(sed->bedload[ised].old_c, elem_old_c, grid->elem2d[ie].nodes);
//            
//            dcdx = 0.; dcdy = 0.;
//            for (i=0; i<NDONTRI; i++) {
//                dcdx += grid->elem2d[ie].grad_shp[i].x * elem_c[i];
//                dcdy += grid->elem2d[ie].grad_shp[i].y * elem_c[i];
//            }
//            
//            coef_0 = (elem_c[0] - elem_old_c[0]) / dt + elem_vel[0].x * dcdx + elem_vel[0].y * dcdy;
//            coef_1 = (elem_c[1] - elem_old_c[1]) / dt + elem_vel[1].x * dcdx + elem_vel[1].y * dcdy;
//            coef_2 = (elem_c[2] - elem_old_c[2]) / dt + elem_vel[2].x * dcdx + elem_vel[2].y * dcdy;
//            
//            /* GLB 0213 changed to subtract out the mean.  This is to account for sources and sinks (such as sediment erosion and deposition)
//             a finite mean should be a good approximaiton of the source/sink term (the mean should be near zero for a conservative solution)
//             if we don't account for this, we will resolve on sources and sinks rather than errors in the solution */
//            error = (1./3.) * (coef_0 + coef_1 + coef_2);
//            coef_0 -= error;
//            coef_1 -= error;
//            coef_2 -= error;
//            
//            coef_0 = coef_0 * coef_0;
//            coef_1 = coef_1 * coef_1;
//            coef_2 = coef_2 * coef_2;
//            error = sqrt(coef_0 + coef_1 + coef_2) * grid->elem2d[ie].djac;
//            
//            if (error > sed->bedload[ised].elem_error[ie]) {
//                sed->bedload[ised].elem_error[ie] = error;
//            }
//            
//            /* scales the tolerance */
//            error /= mat[grid->elem2d[ie].mat].sed[ised].refine_tolerance;
//            
//            /* resets the error if larger */
//            if (error > sw->elem_error[ie]) {
//                sw->elem_error[ie] = error;
//            }
//        }
//        
//    }
    
    // average the elemental to nodal errors for output
    for (i=0; i<grid->nnodes_bed; i++) {
        for (ised = 0; ised < sed->nsed; ised++) {
            sed->bedload[ised].error[i] /= sw->iarray[i];
            
        }
    }
    sarray_int_init(sw->iarray, grid->nnodes); // reset initialize utility array
    
}

#endif

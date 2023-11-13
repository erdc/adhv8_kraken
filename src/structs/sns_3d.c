#include "global_header.h"

//*******************************************************************
//*******************************************************************

static int DEBUG = OFF;

//*******************************************************************
//*******************************************************************
//*******************************************************************
// allocates and initializes the 3d arrays within the SNS_3D struct
void sns_3d_alloc_init(SNS_3D **ns3d, SGRID *grid, SIO *io, SFLAGS flags) {
    int i,j;
    
#ifdef _DEBUG
    printf("-- MYID %d 3D NS initializing\n", grid->smpi->myid);
#endif
    
    /* initialize 3d shallow water model */
    (*ns3d) = (SNS_3D *) tl_alloc(sizeof(SNS_3D), 1);
    SNS_3D *ns = (*ns3d);         // alias
    
    /* allocate ns3d 3d variables */
    int nnodes = grid->nnodes;
    ns->iarray = (int *) tl_alloc(sizeof(int), nnodes);
    ns->darray = (double *) tl_alloc(sizeof(double), nnodes);
    ns->displacement = (double *) tl_alloc(sizeof(double), nnodes);
    ns->old_displacement = (double *) tl_alloc(sizeof(double), nnodes);
    ns->older_displacement = (double *) tl_alloc(sizeof(double), nnodes);
    ns->dpl_perturbation = (double *) tl_alloc(sizeof(double), nnodes);
    ns->grid_speed = (double *) tl_alloc(sizeof(double), nnodes);
    ns->old_grid_speed = tl_alloc(sizeof(double), nnodes);
    ns->density = tl_alloc(sizeof(double), nnodes);
    ns->prs = tl_alloc(sizeof(double), nnodes);
    ns->old_prs = tl_alloc(sizeof(double), nnodes);
    ns->older_prs = tl_alloc(sizeof(double), nnodes);
    ns->error = tl_alloc(sizeof(double), nnodes);
    ns->vertical_node_flux = tl_alloc(sizeof(double), nnodes);
    ns->vel = (SVECT *) tl_alloc(sizeof(SVECT), nnodes);
    ns->old_vel = (SVECT *) tl_alloc(sizeof(SVECT), nnodes);
    ns->older_vel = (SVECT *) tl_alloc(sizeof(SVECT), nnodes);
    ns->tanvec = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    ns->hyd_viscosity = (double *) tl_alloc(sizeof(double), nnodes); /* GSAVANT */
    ns->trn_diffusivity = (double *) tl_alloc(sizeof(double), nnodes); /* GSAVANT */
    ns->depth = (double *) tl_alloc(sizeof(double), nnodes);
#ifdef _SEDIMENT
    ns->bed_displacement = (double *) tl_alloc(sizeof(double), nnodes);
    ns->old_bed_displacement = (double *) tl_alloc(sizeof(double), nnodes);
    ns->older_bed_displacement = (double *) tl_alloc(sizeof(double), nnodes);
#endif
    
    ns->elem_rhs_realloc = 0;
    
    /* winds and waves (optional) */
    ns->winds = NULL;
    ns->waves = NULL;
    if (flags.WIND) swind_alloc(&(ns->winds), nnodes);
    if (flags.WAVE) swave_alloc(&(ns->waves), nnodes);
    
    /* utility array size */
    ns->vwork = NULL;
    ns->vwork_size = 0;
    
    
    /* initialize the arrays */
    sns_3d_init(ns, grid);
    
}

/**************************************************************/
void sns_3d_check_nan(SNS_3D *ns, SGRID *grid){
    int i;
    for (i=0;i<grid->nnodes;i++){
        if(solv_isnan( ns->displacement[i]) || solv_isinf(ns->displacement[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->old_displacement[i]) || solv_isinf(ns->old_displacement[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->older_displacement[i]) || solv_isinf(ns->older_displacement[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->grid_speed[i]) || solv_isinf(ns->grid_speed[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->old_grid_speed[i]) || solv_isinf(ns->old_grid_speed[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->density[i]) || solv_isinf(ns->density[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->prs[i]) || solv_isinf(ns->prs[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->old_prs[i]) || solv_isinf(ns->old_prs[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->older_prs[i]) || solv_isinf(ns->older_prs[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->error[i]) || solv_isinf(ns->error[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->vertical_node_flux[i]) || solv_isinf(ns->vertical_node_flux[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->vel[i].x) || solv_isinf(ns->vel[i].x)) tl_error("Disp NaN");
        if(solv_isnan( ns->old_vel[i].x) || solv_isinf(ns->old_vel[i].x)) tl_error("Disp NaN");
        if(solv_isnan( ns->older_vel[i].x) || solv_isinf(ns->older_vel[i].x)) tl_error("Disp NaN");
        if(solv_isnan( ns->tanvec[i].x) || solv_isinf(ns->tanvec[i].x)) tl_error("Disp NaN");
        if(solv_isnan( ns->hyd_viscosity[i]) || solv_isinf(ns->hyd_viscosity[i])) tl_error("Disp NaN"); /* GSAVANT */
        if(solv_isnan( ns->trn_diffusivity[i]) || solv_isinf(ns->trn_diffusivity[i])) tl_error("Disp NaN"); /* GSAVANT */
        if(solv_isnan( ns->dpl_perturbation[i]) || solv_isinf(ns->dpl_perturbation[i])) tl_error("Disp NaN");
        if(solv_isnan( ns->depth[i]) || solv_isinf(ns->depth[i])) tl_error("Disp NaN");
    }
}

//*******************************************************************
//*******************************************************************
//*******************************************************************
// allocates and initializes the 3d arrays within the NS_3D struct
void sns_3d_init(SNS_3D *ns, SGRID *grid) {
    
    /* allocate ns3d 3d variables */
    int nnodes = grid->nnodes;
    sarray_init_int(ns->iarray, nnodes);
    sarray_init_dbl(ns->darray, nnodes);
    sarray_init_dbl(ns->displacement, nnodes);
    sarray_init_dbl(ns->old_displacement, nnodes);
    sarray_init_dbl(ns->older_displacement, nnodes);
    sarray_init_dbl(ns->dpl_perturbation, nnodes);
    sarray_init_dbl(ns->grid_speed, nnodes);
    sarray_init_dbl(ns->old_grid_speed, nnodes);
    sarray_init_dbl(ns->density, nnodes);
    sarray_init_dbl(ns->prs, nnodes);
    sarray_init_dbl(ns->old_prs, nnodes);
    sarray_init_dbl(ns->older_prs, nnodes);
    sarray_init_dbl(ns->error, nnodes);
    sarray_init_dbl(ns->vertical_node_flux, nnodes);
    svect_init_array(ns->vel, nnodes);
    svect_init_array(ns->old_vel, nnodes);
    svect_init_array(ns->older_vel, nnodes);
    svect2d_init_array(ns->tanvec, nnodes);
    sarray_init_dbl(ns->hyd_viscosity, nnodes); /* GSAVANT */
    sarray_init_dbl(ns->trn_diffusivity, nnodes); /* GSAVANT */
    sarray_init_dbl(ns->depth, nnodes); /* GSAVANT */
#ifdef _SEDIMENT
    sarray_init_dbl(ns->bed_displacement, nnodes);
    sarray_init_dbl(ns->old_bed_displacement, nnodes);
    sarray_init_dbl(ns->older_bed_displacement, nnodes);
#endif
    
    /* allocate ns3d 2d variables */
    nnodes = grid->max_nnodes_sur;
    /*
    sarray_init_dbl(ns->depth, nnodes);
    sarray_init_dbl(ns->old_depth, nnodes);
    svect2d_init_array(ns->depth_avg_vel, nnodes);
    svect2d_init_array(ns->old_depth_avg_vel, nnodes);
    svect2d_init_array(ns->surface_vel, nnodes);
    svect2d_init_array(ns->bottom_vel, nnodes);
    
    int inode;
    for (inode=0; inode<nnodes; inode++) {
        ns->bed_elevation[inode] = grid->node[grid->nodeID_2d_to_3d_bed[inode]].z;
    }
    */
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sns_3d_free(SNS_3D *ns, SGRID *grid, SFLAGS flags) {
    int i;
    
    if (ns != NULL) {
        assert(grid);
        int nnodes;
        
#ifdef _DEBUG
        if (DEBUG) {
            tl_check_all_pickets(__FILE__, __LINE__);
        }
#endif
        
        /* free 3d grid solution arrays */
#ifdef _DEBUG
        printf("... MYID %d freeing ns3d 3d grid solution arrays ...\n", grid->smpi->myid);
#endif
        nnodes = grid->max_nnodes;
        ns->iarray = (int *) tl_free(sizeof(int), nnodes, ns->iarray);
        ns->darray = (double *) tl_free(sizeof(double), nnodes, ns->darray);
        ns->displacement = (double *) tl_free(sizeof(double), nnodes, ns->displacement);
        ns->old_displacement = (double *) tl_free(sizeof(double), nnodes, ns->old_displacement);
        ns->older_displacement = (double *) tl_free(sizeof(double), nnodes, ns->older_displacement);
        ns->dpl_perturbation = (double *) tl_free(sizeof(double), nnodes, ns->dpl_perturbation);
        ns->grid_speed = (double *) tl_free(sizeof(double), nnodes, ns->grid_speed);
        ns->old_grid_speed = (double *) tl_free(sizeof(double), nnodes, ns->old_grid_speed);
        ns->density = (double *) tl_free(sizeof(double), nnodes, ns->density);
        ns->prs = (double *) tl_free(sizeof(double), nnodes, ns->prs);
        ns->old_prs = (double *) tl_free(sizeof(double), nnodes, ns->old_prs);
        ns->older_prs = (double *) tl_free(sizeof(double), nnodes, ns->older_prs);
        ns->error = (double *) tl_free(sizeof(double), nnodes, ns->error);
        ns->vertical_node_flux = (double *) tl_free(sizeof(double), nnodes, ns->vertical_node_flux);
        ns->vel = (SVECT *) tl_free(sizeof(SVECT), nnodes, ns->vel);
        ns->old_vel = (SVECT *) tl_free(sizeof(SVECT), nnodes, ns->old_vel);
        ns->older_vel = (SVECT *) tl_free(sizeof(SVECT), nnodes, ns->older_vel);
        ns->tanvec = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, ns->tanvec);
        ns->hyd_viscosity = (double *) tl_free(sizeof(double), nnodes, ns->hyd_viscosity); /* GSAVANT */
        ns->trn_diffusivity = (double *) tl_free(sizeof(double), nnodes, ns->trn_diffusivity); /* GSAVANT */
        ns->depth = (double *) tl_free(sizeof(double), nnodes, ns->depth); /* GSAVANT */
#ifdef _SEDIMENT
        ns->bed_displacement = (double *) tl_free(sizeof(double), nnodes, ns->bed_displacement);
        ns->old_bed_displacement = (double *) tl_free(sizeof(double), nnodes, ns->old_bed_displacement);
        ns->older_bed_displacement = (double *) tl_free(sizeof(double), nnodes, ns->older_bed_displacement);
#endif
        
        /* free 2d grid solutions arrays */
#ifdef _DEBUG
        printf("... MYID %d freeing ns3d 2d grid solution arrays ...\n", grid->smpi->myid);
#endif
        
        nnodes = grid->max_nnodes_sur;
        /*
        ns->depth = (double *) tl_free(sizeof(double), nnodes, ns->depth);
        ns->old_depth = (double *) tl_free(sizeof(double), nnodes, ns->old_depth);
        ns->depth_avg_vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, ns->depth_avg_vel);
        ns->old_depth_avg_vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, ns->old_depth_avg_vel);
        ns->surface_vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, ns->surface_vel);
        ns->bottom_vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, ns->bottom_vel);
        ns->bed_elevation = (double *) tl_free(sizeof(double), nnodes, ns->bed_elevation);
        */

        if (flags.WIND) swind_free(ns->winds, grid->max_nnodes_sur);
        if (flags.WAVE) swave_free(ns->waves, grid->max_nnodes_sur);
        
        /* free ns struct */
#ifdef _DEBUG
#endif
        ns = (SNS_3D *) tl_free(sizeof(SNS_3D), 1, ns);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// check ns3d solutions for NaNs and InFs
void sns_3d_checkall(SNS_3D *ns, int nnodes2d, int nnodes3d){
    sarray_integrity_check_dbl(ns->displacement, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(ns->old_displacement, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(ns->older_displacement, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(ns->grid_speed, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(ns->old_grid_speed, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(ns->density, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(ns->prs, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(ns->old_prs, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(ns->older_prs, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(ns->error, nnodes3d, __LINE__, __FILE__);
    svect_integrity_check_array(ns->vel, nnodes3d, __LINE__, __FILE__);
    svect_integrity_check_array(ns->old_vel, nnodes3d, __LINE__, __FILE__);
    svect_integrity_check_array(ns->older_vel, nnodes3d, __LINE__, __FILE__);
    svect2d_integrity_check_array(ns->tanvec, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(ns->hyd_viscosity, nnodes3d, __LINE__, __FILE__); /* GSAVANT */
    sarray_integrity_check_dbl(ns->trn_diffusivity, nnodes3d, __LINE__, __FILE__); /* GSAVANT */
    sarray_integrity_check_dbl(ns->depth, nnodes3d, __LINE__, __FILE__);
    
    /*
    sarray_integrity_check_dbl(ns->depth, nnodes2d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(ns->old_depth, nnodes2d, __LINE__, __FILE__);
    svect2d_integrity_check_array(ns->depth_avg_vel, nnodes2d, __LINE__, __FILE__);
    svect2d_integrity_check_array(ns->old_depth_avg_vel, nnodes2d, __LINE__, __FILE__);
    */
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// print a ns3d time-step
void sns_3d_print_ts(SNS_3D *ns, SIO *info, SGRID *grid, double time, int outfact, SFILE sup, char *proj_name, int it1, int it2, int **ndata, int my_nnode_max, int *my_nnode_ext, int **ndata_sur, int my_nnode_max_sur, int *my_nnode_ext_sur, int flag, SFILE_OUTPUT file_output) { 
    FILE *fout;
    // 3d grid prints
#ifdef _DEBUG
    assert(info->fout_ns3_displacement.fp);
#endif
    if (grid->smpi->myid <= 0) fprintf(info->fout_ns3_displacement.fp, "TS 0 %15.8e\n", time);
    sarray_print_dbl(grid, info->fout_ns3_displacement.fp, fout, ns->displacement, ndata, my_nnode_max, my_nnode_ext, 0);

    if (grid->smpi->myid <= 0) fprintf(info->fout_ns3_pressure.fp, "TS 0 %15.8e\n", time);
        sarray_print_dbl(grid, info->fout_ns3_pressure.fp, fout, ns->prs, ndata, my_nnode_max, my_nnode_ext, 0);

    if (grid->smpi->myid <= 0) fprintf(info->fout_ns3_vel.fp, "TS 0 %15.8e\n", time);
    svect_print_array_MPI(grid, info->fout_ns3_vel.fp, fout, ns->vel, ndata, my_nnode_max, my_nnode_ext, 0);

    // maximum nodal continuity errors (here we want to loop over all (initial + adapted) nnodes)
    // elem3d_to_node_double(grid, string, ns->darray, grid->elem_error); // jacobian weighted
    sarray_init_int(ns->iarray, grid->nnodes);
    sarray_init_dbl(ns->darray, grid->nnodes);
    int ie, inode;
    for (ie=0; ie<grid->nelems3d; ie++) {
        for (inode=0; inode<grid->elem3d[ie].nnodes; inode++) {
            ns->darray[grid->elem3d[ie].nodes[inode]] += grid->elem_error[ie];
            ns->iarray[grid->elem3d[ie].nodes[inode]] += 1;
        }
    }
    for (inode=0; inode<grid->nnodes; inode++) {
        ns->darray[inode] /= (double)ns->iarray[inode];
    }
    if (grid->smpi->myid <= 0) fprintf(info->fout_ns3_error.fp, "TS 0 %15.8e\n", time);
    sarray_print_dbl(grid, info->fout_ns3_error.fp, fout, ns->darray, ndata, my_nnode_max, my_nnode_ext, 0);
    
    // hydro nodal continuity errors
    if (grid->smpi->myid <= 0)fprintf(info->fout_ns3_error_hydro.fp, "TS 0 %15.8e\n", time);
    sarray_print_dbl(grid, info->fout_ns3_error_hydro.fp, fout, ns->error, ndata, my_nnode_max, my_nnode_ext, 0);
    
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sns_3d_open_output(SMODEL *mod) {
    
    assert(mod->io); /* should be valid */
    int super = strlen(mod->io->sup.filename); /* whether super file exists */
    
    // 3d grid prints
    
    open_output_file( &(mod->io->fout_ns3_displacement), "ns3d displacement file", super);
    print_header(mod, mod->io->fout_ns3_displacement.fp, PS_FLAG_DPL );

    open_output_file( &(mod->io->fout_ns3_vel), "ns3d velocity file", super);
    print_header(mod, mod->io->fout_ns3_vel.fp, PS_FLAG_VEL );
    
    open_output_file( &(mod->io->fout_ns3_pressure), "ns3d pressure file", super);
    print_header(mod, mod->io->fout_ns3_pressure.fp, PS_FLAG_PRS );

    if (mod->file_output.grid_speed == ON) {
        open_output_file( &(mod->io->fout_ns3_grid_speed), "ns3d vertical grid speed file", super);
        print_header(mod, mod->io->fout_ns3_grid_speed.fp, PS_FLAG_GSP );
    }
    
    open_output_file( &(mod->io->fout_ns3_error), "ns3d error file", super);
    print_header(mod, mod->io->fout_ns3_error.fp, PS_FLAG_ERR );
    
    open_output_file( &(mod->io->fout_ns3_error_hydro), "ns3d hydrodynamic error file", super);
    print_header(mod, mod->io->fout_ns3_error_hydro.fp, PS_FLAG_ERR_HYDRO );
 
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

/* Approximates elemental error for the 3d shallow water equations */
void sns_3d_calculate_elem_error(SNS_3D *ns, SGRID *grid, SMAT *mat, double dt) {
    
    /* calculates the error indicator for the flow */
    int ie = UNSET_INT, i=UNSET_INT;;
    SVECT elem_vel[MAX_NNODES_ON_ELEM3D]; svect_init_array(elem_vel, MAX_NNODES_ON_ELEM3D);
    double error = 0., continuity = 0.;
    
    sarray_init_int(ns->iarray, grid->nnodes); // element per node count for averaging
    sarray_init_dbl(ns->error, grid->nnodes);  // nodal error for printing
    
    int nelems3d = grid->nelems3d; // alias
    SVECT  elem_grad_u, elem_grad_v, elem_grad_w, elem_nds[NDONPRISM];
    double u[NDONPRISM], v[NDONPRISM], w[NDONPRISM], elem_dpl[NDONPRISM];
    double elem_volume;
    for (ie = 0; ie < nelems3d; ie++) {
        error = 0.;
        grid->elem_error[ie] = 0.;
        
        /* calculate the residual for the flow */
        ELEM3D_GET_LOCAL_VECT(ns->vel, elem_vel, grid->elem3d[ie].nodes);
        
        /* Gajanan gkc adding after finding out elem->grad_shp is NaN in current version. */
        int nnodes = grid->elem3d[ie].nnodes;
        SNODE elem_nodes[NDONPRISM];
        for (i=0; i<grid->elem3d[ie].nnodes; i++) {
            snode_copy(&(elem_nodes[i]), grid->node[grid->elem3d[ie].nodes[i]]);
        }
        //printf("\nElem %7i:", ie);
        if (nnodes == NDONTET){
            get_tet_linear_djac_gradPhi(&(grid->elem3d[ie]), elem_nodes, NULL);
            continuity = 0.;
            for (i=0; i<grid->elem3d[ie].nnodes; i++) {
                continuity += svect_dotp(grid->elem3d[ie].grad_shp[i], elem_vel[i]);
            }
            error = sqrt(continuity * continuity) * grid->elem3d[ie].djac;
        }
        else if (nnodes == NDONPRISM){
            /* Do the prism calculations here to find out the continuity error. */
            /* Gajanan gkc temp - quick trials of a simple error indicator. */
            /* What we need to calculate is sqrt(integral over volume of { div(u)*div(u) }) */
            /* What I am doing below is just (integral over volume of { div(u) }). */
            
            for (i=0; i<grid->elem3d[ie].nnodes; i++){
                elem_nds[i].x = elem_nodes[i].x;
                elem_nds[i].y = elem_nodes[i].y;
                elem_nds[i].z = elem_nodes[i].z + elem_dpl[i];
            }
            global_to_local_dbl(ns->displacement, elem_dpl, grid->elem3d[ie].nodes, grid->elem3d[ie].nnodes);
            global_to_local_svect(ns->vel, elem_vel, grid->elem3d[ie].nodes, grid->elem3d[ie].nnodes);
            dumpVector(elem_vel, nnodes, u, v, w);
            
            elem_grad_u = integrate_triPrism_df_full(elem_nds, 1.0, u);
            elem_grad_v = integrate_triPrism_df_full(elem_nds, 1.0, v);
            elem_grad_w = integrate_triPrism_df_full(elem_nds, 1.0, w);
            //printf("\n    elem_grad_u.x = %15.6e, elem_grad_v.y = %15.6e, elem_grad_w.z = %15.6e,", elem_grad_u.x, elem_grad_v.y, elem_grad_w.z);
            //continuity      = fabs(elem_grad_u.x) + fabs(elem_grad_v.y) + fabs(elem_grad_w.z);
            continuity      = fabs(elem_grad_u.x + elem_grad_v.y + elem_grad_w.z);
            error = continuity;
            
            /* For now, to force refinement everywhere. */
            //continuity = 100.0;
            //error = 100.0;
        }
        
        for (i=0; i<grid->elem3d[ie].nnodes; i++) {
            ns->error[grid->elem3d[ie].nodes[i]] += error;
            ns->iarray[grid->elem3d[ie].nodes[i]] += 1;
        }
        
        /* scales by the tolerance */
        error /= mat[grid->elem3d[ie].mat].ns->refine_tolerance;
        
        // sets as max error, since only one NS component can be run for now (con & sed can change)
        grid->elem_error[ie] = error;
    }
    
    // average the elemental to nodal errors for output
    for (i=0; i<grid->nnodes; i++) {
        ns->error[i] /= (double)ns->iarray[i];
    }
    sarray_init_int(ns->iarray, grid->nnodes); // reset initialize utility array
}


/***********************************************************/
/***********************************************************/
/***********************************************************/
// averages two nodes for a  new node
void sns_3d_renumber(SNS_3D *ns, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp, SVECT2D *v2tmp, SVECT *vtmp) {
    // note darray and iarray are just utility arrays, no need to worry about them here
    
    node_renumber_double(max_nnode, ns->displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->old_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->older_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_vect(max_nnode, ns->vel, vtmp, new_numbers, order_tmp);
    node_renumber_vect(max_nnode, ns->old_vel, vtmp, new_numbers, order_tmp);
    node_renumber_vect(max_nnode, ns->older_vel, vtmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, ns->tanvec, v2tmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->density, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->error, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->grid_speed, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->old_grid_speed, dtmp, new_numbers, order_tmp);
    //node_renumber_double(max_nnode, ns->elev_factor, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->prs, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->old_prs, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->older_prs, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->vertical_node_flux, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->depth, dtmp, new_numbers, order_tmp);
#ifdef _SEDIMENT
    node_renumber_double(max_nnode, ns->bed_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->old_bed_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, ns->older_bed_displacement, dtmp, new_numbers, order_tmp);
#endif
    
    
}

/************************************************************/
/************************************************************/
/************************************************************/
// reallocs 3d variables when a new node increment is added
void sns_3d_realloc_init(SNS_3D *ns, int nnodes_old, int nnodes_new) {
    
    /* re-allocate ns3d 3d variables */
    ns->iarray = (int *) tl_realloc(sizeof(int), nnodes_new, nnodes_old, ns->iarray);
    ns->darray = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->darray);
    ns->displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->displacement);
    ns->old_displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->old_displacement);
    ns->older_displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->older_displacement);
    ns->grid_speed = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->grid_speed);
    ns->old_grid_speed = tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->old_grid_speed);
    ns->density = tl_realloc(sizeof(double), nnodes_new, nnodes_old,ns->density);
    ns->prs = tl_realloc(sizeof(double), nnodes_new, nnodes_old,ns->prs);
    ns->old_prs = tl_realloc(sizeof(double), nnodes_new, nnodes_old,ns->old_prs);
    ns->older_prs = tl_realloc(sizeof(double), nnodes_new, nnodes_old,ns->older_prs);
    ns->error = tl_realloc(sizeof(double), nnodes_new, nnodes_old,ns->error);
    ns->vertical_node_flux = tl_realloc(sizeof(double), nnodes_new, nnodes_old,ns->vertical_node_flux);
    ns->vel = (SVECT *) tl_realloc(sizeof(SVECT), nnodes_new, nnodes_old,ns->vel);
    ns->old_vel = (SVECT *) tl_realloc(sizeof(SVECT), nnodes_new, nnodes_old, ns->old_vel);
    ns->older_vel = (SVECT *) tl_realloc(sizeof(SVECT), nnodes_new, nnodes_old, ns->older_vel);
    ns->tanvec = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, ns->tanvec);
    ns->hyd_viscosity = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->hyd_viscosity); /* GSAVANT */
    ns->trn_diffusivity = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->trn_diffusivity); /* GSAVANT */
    ns->depth = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->depth);
#ifdef _SEDIMENT
    ns->bed_displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->bed_displacement);
    ns->old_bed_displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->old_bed_displacement);
    ns->older_bed_displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->older_bed_displacement);
#endif
    
    ns->dpl_perturbation = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, ns->dpl_perturbation);
    
    int i;
    for (i=nnodes_old;i<nnodes_new;i++){
        ns->iarray[i] = UNSET_INT;
        ns->darray[i] = 0.;
        ns->displacement[i] = 0.;
        ns->old_displacement[i] = 0.;
        ns->older_displacement[i] = 0.;
        ns->grid_speed[i] = 0.;
        ns->old_grid_speed[i] = 0.;
        ns->density[i] = 0.;
        ns->prs[i] = 0.;
        ns->old_prs[i] = 0.;
        ns->older_prs[i] = 0.;
        ns->error[i] = 0.;
        ns->vertical_node_flux[i] = 0.;
        svect_init(&ns->vel[i]);
        svect_init(&ns->old_vel[i]);
        svect_init(&ns->older_vel[i]);
        svect2d_init(&ns->tanvec[i]);
        ns->hyd_viscosity[i] = 0.; /* GSAVANT */
        ns->trn_diffusivity[i] = 0.; /* GSAVANT */
        ns->depth[i] = 0.;
#ifdef _SEDIMENT
        ns->bed_displacement[i] = 0.;
        ns->old_bed_displacement[i] = 0. ;
        ns->older_bed_displacement[i] = 0.;
#endif
        ns->dpl_perturbation[i] = 0.;
    }
}
/***********************************************************/
/***********************************************************/
/***********************************************************/
// averages two nodes for a new node - cjt :: since 2d and 3d arrays are numbered different must send both until a map is made
void sns_3d_node_avg(SNS_3D *ns, int node_new, int node1, int node2, SGRID *grid) {
    
    // 3D :: (cjt) note: darray and iarray are just utility arrays, no need to worry about them here
    ns->displacement[node_new] = 0.5*(ns->displacement[node1] + ns->displacement[node2]);
    ns->old_displacement[node_new] = 0.5*(ns->old_displacement[node1] + ns->old_displacement[node2]);
    ns->older_displacement[node_new] = 0.5*(ns->older_displacement[node1] + ns->older_displacement[node2]);
    ns->grid_speed[node_new] = 0.5*(ns->grid_speed[node1] + ns->grid_speed[node2]);
    ns->old_grid_speed[node_new] = 0.5*(ns->old_grid_speed[node1] + ns->old_grid_speed[node2]);
    ns->density[node_new] = 0.5*(ns->density[node1] + ns->density[node2]);
    ns->prs[node_new] = 0.5*(ns->prs[node1] + ns->prs[node2]);
    ns->old_prs[node_new] = 0.5*(ns->old_prs[node1] + ns->old_prs[node2]);
    ns->older_prs[node_new] = 0.5*(ns->older_prs[node1] + ns->older_prs[node2]);
    ns->error[node_new] = 0.5*(ns->error[node1] + ns->error[node2]);
    ns->vertical_node_flux[node_new] = 0.5*(ns->vertical_node_flux[node1] + ns->vertical_node_flux[node2]);
    ns->vel[node_new] = svect_avg(ns->vel[node1], ns->vel[node2]);
    ns->old_vel[node_new] = svect_avg(ns->old_vel[node1], ns->old_vel[node2]);
    ns->older_vel[node_new] = svect_avg(ns->older_vel[node1], ns->older_vel[node2]);
    ns->tanvec[node_new] = svect2d_avg(ns->tanvec[node1], ns->tanvec[node2]);
    ns->hyd_viscosity[node_new] = 0.5*(ns->hyd_viscosity[node1] + ns->hyd_viscosity[node2]);
    ns->trn_diffusivity[node_new] = 0.5*(ns->trn_diffusivity[node1] + ns->trn_diffusivity[node2]);
    ns->depth[node_new] = 0.5*(ns->depth[node1] + ns->depth[node2]);
#ifdef _SEDIMENT
    ns->bed_displacement[node_new] = 0.5*(ns->bed_displacement[node1] + ns->bed_displacement[node2]);
    ns->old_bed_displacement[node_new] = 0.5*(ns->old_bed_displacement[node1] + ns->old_bed_displacement[node2]);
    ns->older_bed_displacement[node_new] = 0.5*(ns->older_bed_displacement[node1] + ns->older_bed_displacement[node2]);
#endif
    ns->dpl_perturbation[node_new] = 0.5*(ns->dpl_perturbation[node1] + ns->dpl_perturbation[node2]);
    
}


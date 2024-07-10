#include "global_header.h"

//*******************************************************************
//*******************************************************************

static int DEBUG = OFF;

//*******************************************************************
//*******************************************************************
//*******************************************************************
// allocates and initializes the 3d arrays within the SSW_3D struct
void ssw_3d_alloc_init(SSW_3D **sw3d, SGRID *grid, SIO *io, SFLAGS flags) {
    int i,j;
 
#ifdef _DEBUG
    int myid = 0, ierr = UNSET_INT;;
#ifdef _MESSG
    ierr = MPI_Comm_rank(cstorm_comm, &myid);
#endif
    printf("\n---- MYID %d :: MPI_COMM_WORLD_RANK: %d :: 3D SW INTIALIZATION\n", grid->smpi->myid,myid);
#endif
    /* initialize 3d shallow water model */
    (*sw3d) = (SSW_3D *) tl_alloc(sizeof(SSW_3D), 1);
    SSW_3D *sw = (*sw3d);         // alias
    
    /* allocate sw3d 3d variables */
    int nnodes = grid->nnodes;
    sw->iarray = (int *) tl_alloc(sizeof(int), nnodes);
    sw->darray = (double *) tl_alloc(sizeof(double), nnodes);
    sw->displacement = (double *) tl_alloc(sizeof(double), nnodes);
    sw->old_displacement = (double *) tl_alloc(sizeof(double), nnodes);
    sw->older_displacement = (double *) tl_alloc(sizeof(double), nnodes);
    sw->dpl_perturbation = (double *) tl_alloc(sizeof(double), nnodes);
    sw->grid_speed = (double *) tl_alloc(sizeof(double), nnodes);
    sw->old_grid_speed = tl_alloc(sizeof(double), nnodes);
    sw->density = tl_alloc(sizeof(double), nnodes);
    sw->prs = tl_alloc(sizeof(double), nnodes);
    sw->prs_plus = tl_alloc(sizeof(double), nnodes);
    sw->prs_minus = tl_alloc(sizeof(double), nnodes);
    sw->error = tl_alloc(sizeof(double), nnodes);
    sw->vertical_node_flux = tl_alloc(sizeof(double), nnodes);
    sw->vel = (SVECT *) tl_alloc(sizeof(SVECT), nnodes);
    sw->old_vel = (SVECT *) tl_alloc(sizeof(SVECT), nnodes);
    sw->older_vel = (SVECT *) tl_alloc(sizeof(SVECT), nnodes);
    sw->tanvec = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    sw->hyd_viscosity = (double *) tl_alloc(sizeof(double), nnodes); /* GSAVANT */
    sw->trn_diffusivity = (double *) tl_alloc(sizeof(double), nnodes); /* GSAVANT */
    sw->grad_bed = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
#ifdef _SEDIMENT
    sw->bed_displacement = (double *) tl_alloc(sizeof(double), nnodes);
    sw->old_bed_displacement = (double *) tl_alloc(sizeof(double), nnodes);
    sw->older_bed_displacement = (double *) tl_alloc(sizeof(double), nnodes);
#endif
    
    if (debug.no_hydro == OFF) {
        sw->elem_rhs_supg_dacont = (double **) tl_alloc(sizeof(double *), MAX_NNODES_ON_ELEM3D);
        sw->elem_rhs_supg_cont = (double **) tl_alloc(sizeof(double *), MAX_NNODES_ON_ELEM3D);
        for (i=0; i<MAX_NNODES_ON_ELEM3D; i++) {
            sw->elem_rhs_supg_dacont[i] = (double *) tl_alloc(sizeof(double), grid->nelems3d);
            sw->elem_rhs_supg_cont[i] = (double *) tl_alloc(sizeof(double), grid->nelems3d);
        }
    }
    
    sw->elem_rhs_realloc = 0;
    /* allocate sw3d 2d variables */
    nnodes = grid->max_nnodes_sur;
    sw->depth = (double *) tl_alloc(sizeof(double), nnodes);
    sw->old_depth = (double *) tl_alloc(sizeof(double), nnodes);
    sw->depth_avg_vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    sw->old_depth_avg_vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    sw->surface_vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    sw->bottom_vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    sw->bed_elevation = (double *) tl_alloc(sizeof(double), nnodes);
    
    /* winds and waves (optional) */
    sw->winds = NULL;
    sw->waves = NULL;
    if (flags.WIND) swind_alloc(&(sw->winds), nnodes);
    if (flags.WAVE) swave_alloc(&(sw->waves), nnodes);
    
    /* utility array size */
    sw->vwork = NULL;
    sw->vwork_size = 0;
    
    
    /* initialize the arrays */
    ssw_3d_init(sw, grid);
    
#ifdef _DEBUG
    printf("\n---- MYID %d :: MPI_COMM_WORLD_RANK: %d :: 3D SW FINALIZING INTIALIZATION\n", grid->smpi->myid,myid);
#endif
}

/**************************************************************/
void ssw_3d_check_nan(SSW_3D *sw, SGRID *grid){
    int i;
    for (i=0;i<grid->nnodes;i++){
        if(solv_isnan(sw->displacement[i]) || solv_isinf(sw->displacement[i])) tl_error("Disp NaN");
        if(solv_isnan( sw->old_displacement[i]) || solv_isinf(sw->old_displacement[i])) tl_error("Disp NaN");
        if(solv_isnan( sw->older_displacement[i]) || solv_isinf(sw->older_displacement[i])) tl_error("Disp NaN");
        if(solv_isnan( sw->grid_speed[i]) || solv_isinf(sw->grid_speed[i])) tl_error("Disp NaN");
        if(solv_isnan( sw->old_grid_speed[i]) || solv_isinf(sw->old_grid_speed[i])) tl_error("Disp NaN");
        if(solv_isnan( sw->density[i]) || solv_isinf(sw->density[i])) tl_error("Disp NaN");
        if(solv_isnan( sw->prs[i]) || solv_isinf(sw->prs[i])) tl_error("Disp NaN");
        if(solv_isnan( sw->prs_plus[i]) || solv_isinf(sw->prs_plus[i])) tl_error("Disp NaN");
        if(solv_isnan( sw->prs_minus[i]) || solv_isinf(sw->prs_minus[i])) tl_error("Disp NaN");
        if(solv_isnan( sw->error[i]) || solv_isinf(sw->error[i])) tl_error("Disp NaN");
        if(solv_isnan( sw->vertical_node_flux[i]) || solv_isinf(sw->vertical_node_flux[i])) tl_error("Disp NaN");
        if(solv_isnan( sw->vel[i].x) || solv_isinf(sw->vel[i].x)) tl_error("Disp NaN");
        if(solv_isnan( sw->old_vel[i].x) || solv_isinf(sw->old_vel[i].x)) tl_error("Disp NaN");
        if(solv_isnan( sw->older_vel[i].x) || solv_isinf(sw->older_vel[i].x)) tl_error("Disp NaN");
        if(solv_isnan( sw->tanvec[i].x) || solv_isinf(sw->tanvec[i].x)) tl_error("Disp NaN");
        if(solv_isnan( sw->hyd_viscosity[i]) || solv_isinf(sw->hyd_viscosity[i])) tl_error("Disp NaN"); /* GSAVANT */
        if(solv_isnan( sw->trn_diffusivity[i]) || solv_isinf(sw->trn_diffusivity[i])) tl_error("Disp NaN"); /* GSAVANT */
        if(solv_isnan( sw->dpl_perturbation[i]) || solv_isinf(sw->dpl_perturbation[i])) tl_error("Disp NaN");
    }
}

//*******************************************************************
//*******************************************************************
//*******************************************************************
// allocates and initializes the 3d arrays within the SSW_3D struct
void ssw_3d_init(SSW_3D *sw, SGRID *grid) {
    
    /* allocate sw3d 3d variables */
    int nnodes = grid->nnodes;
    sarray_init_int(sw->iarray, nnodes);
    sarray_init_dbl(sw->darray, nnodes);
    sarray_init_dbl(sw->displacement, nnodes);
    sarray_init_dbl(sw->old_displacement, nnodes);
    sarray_init_dbl(sw->older_displacement, nnodes);
    sarray_init_dbl(sw->dpl_perturbation, nnodes);
    sarray_init_dbl(sw->grid_speed, nnodes);
    sarray_init_dbl(sw->old_grid_speed, nnodes);
    sarray_init_dbl(sw->density, nnodes);
    sarray_init_dbl(sw->prs, nnodes);
    sarray_init_dbl(sw->prs_plus, nnodes);
    sarray_init_dbl(sw->prs_minus, nnodes);
    sarray_init_dbl(sw->error, nnodes);
    sarray_init_dbl(sw->vertical_node_flux, nnodes);
    svect_init_array(sw->vel, nnodes);
    svect_init_array(sw->old_vel, nnodes);
    svect_init_array(sw->older_vel, nnodes);
    svect2d_init_array(sw->tanvec, nnodes);
    sarray_init_dbl(sw->hyd_viscosity, nnodes); /* GSAVANT */
    sarray_init_dbl(sw->trn_diffusivity, nnodes); /* GSAVANT */
    svect2d_init_array(sw->grad_bed, nnodes);
#ifdef _SEDIMENT
    sarray_init_dbl(sw->bed_displacement, nnodes);
    sarray_init_dbl(sw->old_bed_displacement, nnodes);
    sarray_init_dbl(sw->older_bed_displacement, nnodes);
#endif
    
    int i;
    if (debug.no_hydro == OFF) {
        for (i=0; i<MAX_NNODES_ON_ELEM3D; i++) {
            sarray_init_dbl(sw->elem_rhs_supg_dacont[i], grid->nelems3d);
            sarray_init_dbl(sw->elem_rhs_supg_cont[i], grid->nelems3d);
        }
    }
    
    /* allocate sw3d 2d variables */
    nnodes = grid->max_nnodes_sur;
    sarray_init_dbl(sw->depth, nnodes);
    sarray_init_dbl(sw->old_depth, nnodes);
    svect2d_init_array(sw->depth_avg_vel, nnodes);
    svect2d_init_array(sw->old_depth_avg_vel, nnodes);
    svect2d_init_array(sw->surface_vel, nnodes);
    svect2d_init_array(sw->bottom_vel, nnodes);
    
    int inode;
    for (inode=0; inode<nnodes; inode++) {
        sw->bed_elevation[inode] = grid->node[grid->nodeID_2d_to_3d_bed[inode]].z;
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void ssw_3d_free(SSW_3D *sw, SGRID *grid, SFLAGS flags) {
    int i;
    
    if (sw != NULL) {
        assert(grid);
        int nnodes;
        
#ifdef _DEBUG
        if (DEBUG) {
            tl_check_all_pickets(__FILE__, __LINE__);
        }
#endif
        
        /* free 3d grid solution arrays */
#ifdef _DEBUG
        printf("... MYID %d freeing sw3d 3d grid solution arrays ...\n", grid->smpi->myid);
#endif
        nnodes = grid->max_nnodes;
        sw->iarray = (int *) tl_free(sizeof(int), nnodes, sw->iarray);
        sw->darray = (double *) tl_free(sizeof(double), nnodes, sw->darray);
        sw->displacement = (double *) tl_free(sizeof(double), nnodes, sw->displacement);
        sw->old_displacement = (double *) tl_free(sizeof(double), nnodes, sw->old_displacement);
        sw->older_displacement = (double *) tl_free(sizeof(double), nnodes, sw->older_displacement);
        sw->dpl_perturbation = (double *) tl_free(sizeof(double), nnodes, sw->dpl_perturbation);
        sw->grid_speed = (double *) tl_free(sizeof(double), nnodes, sw->grid_speed);
        sw->old_grid_speed = (double *) tl_free(sizeof(double), nnodes, sw->old_grid_speed);
        sw->density = (double *) tl_free(sizeof(double), nnodes, sw->density);
        sw->prs = (double *) tl_free(sizeof(double), nnodes, sw->prs);
        sw->prs_plus = (double *) tl_free(sizeof(double), nnodes, sw->prs_plus);
        sw->prs_minus = (double *) tl_free(sizeof(double), nnodes, sw->prs_minus);
        sw->error = (double *) tl_free(sizeof(double), nnodes, sw->error);
        sw->vertical_node_flux = (double *) tl_free(sizeof(double), nnodes, sw->vertical_node_flux);
        sw->vel = (SVECT *) tl_free(sizeof(SVECT), nnodes, sw->vel);
        sw->old_vel = (SVECT *) tl_free(sizeof(SVECT), nnodes, sw->old_vel);
        sw->older_vel = (SVECT *) tl_free(sizeof(SVECT), nnodes, sw->older_vel);
        sw->tanvec = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, sw->tanvec);
        sw->hyd_viscosity = (double *) tl_free(sizeof(double), nnodes, sw->hyd_viscosity); /* GSAVANT */
        sw->trn_diffusivity = (double *) tl_free(sizeof(double), nnodes, sw->trn_diffusivity); /* GSAVANT */
        sw->grad_bed = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, sw->grad_bed);
#ifdef _SEDIMENT
        sw->bed_displacement = (double *) tl_free(sizeof(double), nnodes, sw->bed_displacement);
        sw->old_bed_displacement = (double *) tl_free(sizeof(double), nnodes, sw->old_bed_displacement);
        sw->older_bed_displacement = (double *) tl_free(sizeof(double), nnodes, sw->older_bed_displacement);
#endif
        
        if (debug.no_hydro == OFF) {
            for (i=0; i<MAX_NNODES_ON_ELEM3D; i++) {
                if ((flags.ADAPTED_THE_GRID == YES) && (sw->elem_rhs_realloc == 0)) {
                    sw->elem_rhs_supg_dacont[i] = (double *) tl_free(sizeof(double), grid->nelems3d_old, sw->elem_rhs_supg_dacont[i]);
                    sw->elem_rhs_supg_cont[i] = (double *) tl_free(sizeof(double), grid->nelems3d_old, sw->elem_rhs_supg_cont[i]);
                } else {
                    sw->elem_rhs_supg_dacont[i] = (double *) tl_free(sizeof(double), grid->nelems3d, sw->elem_rhs_supg_dacont[i]);
                    sw->elem_rhs_supg_cont[i] = (double *) tl_free(sizeof(double), grid->nelems3d, sw->elem_rhs_supg_cont[i]);
                }
            }
            sw->elem_rhs_supg_dacont = (double **) tl_free(sizeof(double *), MAX_NNODES_ON_ELEM3D, sw->elem_rhs_supg_dacont);
            sw->elem_rhs_supg_cont = (double **) tl_free(sizeof(double *), MAX_NNODES_ON_ELEM3D, sw->elem_rhs_supg_cont);
        }
        
        /* free 2d grid solutions arrays */
#ifdef _DEBUG
        printf("... MYID %d freeing sw3d 2d grid solution arrays ...\n", grid->smpi->myid);
#endif
        
        nnodes = grid->max_nnodes_sur;
        sw->depth = (double *) tl_free(sizeof(double), nnodes, sw->depth);
        sw->old_depth = (double *) tl_free(sizeof(double), nnodes, sw->old_depth);
        sw->depth_avg_vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, sw->depth_avg_vel);
        sw->old_depth_avg_vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, sw->old_depth_avg_vel);
        sw->surface_vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, sw->surface_vel);
        sw->bottom_vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, sw->bottom_vel);
        sw->bed_elevation = (double *) tl_free(sizeof(double), nnodes, sw->bed_elevation);
        
        if (flags.WIND) swind_free(sw->winds, grid->max_nnodes_sur);
        if (flags.WAVE) swave_free(sw->waves, grid->max_nnodes_sur);
        
        /* free sw3 struct */
#ifdef _DEBUG
#endif
        sw = (SSW_3D *) tl_free(sizeof(SSW_3D), 1, sw);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// check sw3d solutions for NaNs and InFs
void ssw_3d_checkall(SSW_3D *sw, int nnodes2d, int nnodes3d){
    sarray_integrity_check_dbl(sw->displacement, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(sw->old_displacement, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(sw->older_displacement, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(sw->grid_speed, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(sw->old_grid_speed, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(sw->density, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(sw->prs, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(sw->prs_plus, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(sw->prs_minus, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(sw->error, nnodes3d, __LINE__, __FILE__);
    svect_integrity_check_array(sw->vel, nnodes3d, __LINE__, __FILE__);
    svect_integrity_check_array(sw->old_vel, nnodes3d, __LINE__, __FILE__);
    svect_integrity_check_array(sw->older_vel, nnodes3d, __LINE__, __FILE__);
    svect2d_integrity_check_array(sw->tanvec, nnodes3d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(sw->hyd_viscosity, nnodes3d, __LINE__, __FILE__); /* GSAVANT */
    sarray_integrity_check_dbl(sw->trn_diffusivity, nnodes3d, __LINE__, __FILE__); /* GSAVANT */
    
    sarray_integrity_check_dbl(sw->depth, nnodes2d, __LINE__, __FILE__);
    sarray_integrity_check_dbl(sw->old_depth, nnodes2d, __LINE__, __FILE__);
    svect2d_integrity_check_array(sw->depth_avg_vel, nnodes2d, __LINE__, __FILE__);
    svect2d_integrity_check_array(sw->old_depth_avg_vel, nnodes2d, __LINE__, __FILE__);
    
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// print a sw3d time-step
void ssw_3d_print_ts(SSW_3D *sw, SIO *info, SGRID *grid, double time, int outfact, SFILE sup, char *proj_name, int it1, int it2, int **ndata, int my_nnode_max, int *my_nnode_ext, int **ndata_sur, int my_nnode_max_sur, int *my_nnode_ext_sur, int flag, SFILE_OUTPUT file_output) {
    
    FILE *fout;
    // 3d grid prints
    if (grid->smpi->myid <= 0) fprintf(info->fout_sw3_displacement.fp, "TS 0 %15.8e\n", time);
    sarray_print_dbl(grid, info->fout_sw3_displacement.fp, fout, sw->displacement, ndata, my_nnode_max, my_nnode_ext, 0);
    if (grid->smpi->myid <= 0) fprintf(info->fout_sw3_vel.fp, "TS 0 %15.8e\n", time);
    svect_print_array_MPI(grid, info->fout_sw3_vel.fp, fout, sw->vel, ndata, my_nnode_max, my_nnode_ext, 0);
    if (file_output.grid_speed == ON) {
        if (grid->smpi->myid <= 0) fprintf(info->fout_sw3_grid_speed.fp, "TS 0 %15.8e\n", time);
        sarray_print_dbl(grid, info->fout_sw3_grid_speed.fp, fout, sw->grid_speed, ndata, my_nnode_max, my_nnode_ext, 0);
    }
    if (file_output.pressure == ON) {
        if (grid->smpi->myid <= 0) fprintf(info->fout_sw3_pressure.fp, "TS 0 %15.8e\n", time);
        sarray_print_dbl(grid, info->fout_sw3_pressure.fp, fout, sw->prs, ndata, my_nnode_max, my_nnode_ext, 0);
    }
    if (file_output.hyd_vis == ON) {
        if (grid->smpi->myid <= 0) fprintf(info->fout_sw3_hyd_viscosity.fp, "TS 0 %15.8e\n", time);
        sarray_print_dbl(grid, info->fout_sw3_hyd_viscosity.fp, fout, sw->hyd_viscosity, ndata, my_nnode_max, my_nnode_ext, 0);
    }
    if (file_output.trn_dif == ON) {
        if (grid->smpi->myid <= 0) fprintf(info->fout_sw3_trn_diffusivity.fp, "TS 0 %15.8e\n", time);
        sarray_print_dbl(grid, info->fout_sw3_trn_diffusivity.fp, fout, sw->trn_diffusivity, ndata, my_nnode_max, my_nnode_ext, 0);
    }
    
    // maximum nodal continuity errors (here we want to loop over all (initial + adapted) nnodes)
    // elem3d_to_node_double(grid, string, sw->darray, grid->elem_error); // jacobian weighted
    sarray_init_int(sw->iarray, grid->nnodes);
    sarray_init_dbl(sw->darray, grid->nnodes);
    int ie, inode;
    for (ie=0; ie<grid->nelems3d; ie++) {
        for (inode=0; inode<grid->elem3d[ie].nnodes; inode++) {
            sw->darray[grid->elem3d[ie].nodes[inode]] += grid->elem_error[ie];
            sw->iarray[grid->elem3d[ie].nodes[inode]] += 1;
        }
    }
    for (inode=0; inode<grid->nnodes; inode++) {
        sw->darray[inode] /= (double)sw->iarray[inode];
    }
    if (grid->smpi->myid <= 0) fprintf(info->fout_sw3_error.fp, "TS 0 %15.8e\n", time);
    sarray_print_dbl(grid, info->fout_sw3_error.fp, fout, sw->darray, ndata, my_nnode_max, my_nnode_ext, 0);
    
    // hydro nodal continuity errors
    if (grid->smpi->myid <= 0)fprintf(info->fout_sw3_error_hydro.fp, "TS 0 %15.8e\n", time);
    sarray_print_dbl(grid, info->fout_sw3_error_hydro.fp, fout, sw->error, ndata, my_nnode_max, my_nnode_ext, 0);
    
    // 2d grid prints
    if (grid->smpi->myid <= 0) fprintf(info->fout_sw3_depth.fp, "TS 0 %15.8e\n", time);
    sarray_print_dbl(grid, info->fout_sw3_depth.fp, fout, sw->depth, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur, 0);
    if (file_output.depth_avg_velocity == ON) {
        if (grid->smpi->myid <= 0) fprintf(info->fout_sw3_depth_avg_vel.fp, "TS 0 %15.8e\n", time);
        svect2d_print_array_MPI(grid, info->fout_sw3_depth_avg_vel.fp, fout, sw->depth_avg_vel, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur, 0);
    }
    if (file_output.surface_velocity == ON) {
        if (grid->smpi->myid <= 0) fprintf(info->fout_sw3_surface_vel.fp, "TS 0 %15.8e\n", time);
        svect2d_print_array_MPI(grid, info->fout_sw3_surface_vel.fp, fout, sw->surface_vel, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur, 0);
    }
    if (file_output.bed_velocity == ON) {
        if (grid->smpi->myid <= 0) fprintf(info->fout_sw3_bottom_vel.fp, "TS 0 %15.8e\n", time);
        svect2d_print_array_MPI(grid, info->fout_sw3_bottom_vel.fp, fout, sw->bottom_vel, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur, 0);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// print a sw3d time-step
/*void ssw_3d_print_ts_MPI(SSW_3D *sw, SIO *info, double time, SGRID *grid, STR_VALUE *string) {
 
 int nnodes3d = grid->orig_macro_nnodes;        // first nnodes are the original ones
 int nnodes2d = grid->macro_nnodes_bed;    // first nnodes are the original ones
 #ifdef _MESSG
 
 // print headers
 if (grid->smpi->myid <= 0 ){
 fprintf(info->fout_sw3_displacement.fp, "TS 0 %15.8e\n", time);
 fprintf(info->fout_sw3_vel.fp, "TS 0 %15.8e\n", time);
 if (file_output.grid_speed == ON) {
 fprintf(info->fout_sw3_grid_speed.fp, "TS 0 %15.8e\n", time);
 }
 if (file_output.pressure == ON) {
 fprintf(info->fout_sw3_pressure.fp, "TS 0 %15.8e\n", time);
 }
 fprintf(info->fout_sw3_error.fp, "TS 0 %15.8e\n", time);
 fprintf(info->fout_sw3_error_hydro.fp, "TS 0 %15.8e\n", time);
 fprintf(info->fout_sw3_depth.fp, "TS 0 %15.8e\n", time);
 if (file_output.depth_avg_velocity == ON) {
 fprintf(info->fout_sw3_depth_avg_vel.fp, "TS 0 %15.8e\n", time);
 }
 if (file_output.surface_velocity == ON) {
 fprintf(info->fout_sw3_surface_vel.fp, "TS 0 %15.8e\n", time);
 }
 if (file_output.bed_velocity == ON) {
 fprintf(info->fout_sw3_bottom_vel.fp, "TS 0 %15.8e\n", time);
 }
 }
 
 
 // 3d grid prints
 sarray_print_dbl_MPI(grid, info->fout_sw3_displacement.fp, sw->displacement, 0);
 svect_print_array_MPI(grid, info->fout_sw3_vel.fp, sw->vel);
 if (file_output.grid_speed == ON) {
 sarray_print_dbl_MPI(grid, info->fout_sw3_grid_speed.fp, sw->grid_speed, 0);
 }
 if (file_output.pressure == ON) {
 sarray_print_dbl_MPI(grid, info->fout_sw3_pressure.fp, sw->prs, 0);
 }
 
 // maximum nodal continuity errors (here we want to loop over all (initial + adapted) nnodes)
 // elem3d_to_node_double(grid, string, sw->darray, grid->elem_error); // jacobian weighted
 sarray_init_int(sw->iarray, grid->nnodes);
 sarray_init_dbl(sw->darray, grid->nnodes);
 int ie, inode;
 for (ie=0; ie<grid->nelems3d; ie++) {
 for (inode=0; inode<MAX_NNODES_ON_ELEM3D; inode++) {
 sw->darray[grid->elem3d[ie].nodes[inode]] += grid->elem_error[ie];
 sw->iarray[grid->elem3d[ie].nodes[inode]] += 1;
 }
 }
 for (inode=0; inode<grid->nnodes; inode++) {
 sw->darray[inode] /= (double)sw->iarray[inode];
 }
 
 sarray_print_dbl_MPI(grid, info->fout_sw3_error.fp, sw->darray, 0);
 
 // hydro nodal continuity errors
 sarray_print_dbl_MPI(grid, info->fout_sw3_error_hydro.fp, sw->error, 0);
 
 // 2d grid prints
 
 sarray_print_dbl_MPI(grid, info->fout_sw3_depth.fp, sw->depth, 1);
 if (file_output.depth_avg_velocity == ON) {
 svect2d_print_array_MPI(grid, info->fout_sw3_depth_avg_vel.fp, sw->depth_avg_vel, 1);
 }
 if (file_output.surface_velocity == ON) {
 svect2d_print_array_MPI(grid, info->fout_sw3_surface_vel.fp, sw->surface_vel, 1);
 }
 if (file_output.bed_velocity == ON) {
 svect2d_print_array_MPI(grid, info->fout_sw3_bottom_vel.fp, sw->bottom_vel, 1);
 }
 
 #endif
 
 }*/

/***********************************************************/
/***********************************************************/
/***********************************************************/

void ssw_3d_open_output(SMODEL *mod) {
    
    assert(mod->io); /* should be valid */
    int super = strlen(mod->io->sup.filename); /* whether super file exists */
    
    // 3d grid prints
    open_output_file( &(mod->io->fout_sw3_displacement), "sw3d displacement file", super);
    print_header(mod, mod->io->fout_sw3_displacement.fp, PS_FLAG_DPL );
    
    open_output_file( &(mod->io->fout_sw3_vel), "sw3d velocity file", super);
    print_header(mod, mod->io->fout_sw3_vel.fp, PS_FLAG_VEL );
    
    if (mod->file_output.grid_speed == ON) {
        open_output_file( &(mod->io->fout_sw3_grid_speed), "sw3d vertical grid speed file", super);
        print_header(mod, mod->io->fout_sw3_grid_speed.fp, PS_FLAG_GSP );
    }
    
    if (mod->file_output.pressure == ON) {
        open_output_file( &(mod->io->fout_sw3_pressure), "sw3d pressure file", super);
        print_header(mod, mod->io->fout_sw3_pressure.fp, PS_FLAG_PRS );
    }
    
    if (mod->file_output.hyd_vis == ON) {
        open_output_file( &(mod->io->fout_sw3_hyd_viscosity), "sw3d Hydro Viscosity file", super);
        print_header(mod, mod->io->fout_sw3_hyd_viscosity.fp, PS_FLAG_HYV );   /* GSAVANT */
    }
    
    if (mod->file_output.trn_dif == ON) {
        open_output_file( &(mod->io->fout_sw3_trn_diffusivity), "sw3d Hydro Diffusion file", super);
        print_header(mod, mod->io->fout_sw3_trn_diffusivity.fp, PS_FLAG_TRD );  /* GSAVANT */
    }
    //    open_output_file( &(mod->io->fout_sw3_old_displacement), "sw3d old displacement file", super);
    //    print_header(mod, mod->io->fout_sw3_old_displacement.fp, PS_FLAG_PDPL );
    //
    //    open_output_file( &(mod->io->fout_sw3_old_vel), "sw3d old velocity file", super);
    //    print_header(mod, mod->io->fout_sw3_old_vel.fp, PS_FLAG_PVEL );
    //
    open_output_file( &(mod->io->fout_sw3_error), "sw3d error file", super);
    print_header(mod, mod->io->fout_sw3_error.fp, PS_FLAG_ERR );
    
    open_output_file( &(mod->io->fout_sw3_error_hydro), "sw3d hydrodynamic error file", super);
    print_header(mod, mod->io->fout_sw3_error_hydro.fp, PS_FLAG_ERR_HYDRO );
    
    // 2d grid prints
    open_output_file( &(mod->io->fout_sw3_depth), "sw3d 2d depth file", super);
    print_header(mod, mod->io->fout_sw3_depth.fp, PS_FLAG_DEPTH);
    
    if (mod->file_output.depth_avg_velocity == ON) {
        open_output_file( &(mod->io->fout_sw3_depth_avg_vel), "sw3d 2d depth averaged velocity file", super);
        print_header(mod, mod->io->fout_sw3_depth_avg_vel.fp, PS_FLAG_DAVG_VEL );
    }
    
    if (mod->file_output.surface_velocity == ON) {
        open_output_file( &(mod->io->fout_sw3_surface_vel), "sw3d 2d surface velocity file", super);
        print_header(mod, mod->io->fout_sw3_surface_vel.fp, PS_FLAG_SURF_VEL );
    }
    
    if (mod->file_output.bed_velocity == ON) {
        open_output_file( &(mod->io->fout_sw3_bottom_vel), "sw3d 2d bedvelocity file", super);
        print_header(mod, mod->io->fout_sw3_bottom_vel.fp, PS_FLAG_BOTT_VEL );
    }
}


/***********************************************************/
/***********************************************************/
/***********************************************************/

void ssw_3d_open_input(SMODEL *mod) {
    
    //assert(mod->io); /* should be valid */
    //int super = strlen(mod->io->sup.filename); /* whether super file exists */
    
    /* 2d face file if SW3d is used */
    //open_input_file( &(mod->io->face), "boundary face file", super);
    
}


/***********************************************************/
/***********************************************************/
/***********************************************************/

/* Approximates elemental error for the 3d shallow water equations */
void ssw_3d_calculate_elem_error(SSW_3D *sw, SGRID *grid, SMAT *mat, double dt) {
    
    /* calculates the error indicator for the flow */
    int ie = UNSET_INT, i=UNSET_INT;;
    SVECT elem_vel[MAX_NNODES_ON_ELEM3D]; svect_init_array(elem_vel, MAX_NNODES_ON_ELEM3D);
    double error = 0., continuity = 0.;
    
    sarray_init_int(sw->iarray, grid->nnodes); // element per node count for averaging
    sarray_init_dbl(sw->error, grid->nnodes);  // nodal error for printing
    
    int nelems3d = grid->nelems3d; // alias
    SVECT  elem_grad_u, elem_grad_v, elem_grad_w, elem_nds[NDONPRISM];
    double u[NDONPRISM], v[NDONPRISM], w[NDONPRISM], elem_dpl[NDONPRISM];
    double elem_volume;
    for (ie = 0; ie < nelems3d; ie++) {
        error = 0.;
        grid->elem_error[ie] = 0.;
        
        /* calculate the residual for the flow */
        ELEM3D_GET_LOCAL_VECT(sw->vel, elem_vel, grid->elem3d[ie].nodes);
        
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
            global_to_local_dbl(sw->displacement, elem_dpl, grid->elem3d[ie].nodes, grid->elem3d[ie].nnodes);
            global_to_local_svect(sw->vel, elem_vel, grid->elem3d[ie].nodes, grid->elem3d[ie].nnodes);
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
            sw->error[grid->elem3d[ie].nodes[i]] += error;
            sw->iarray[grid->elem3d[ie].nodes[i]] += 1;
        }
        
        /* scales by the tolerance */
        error /= mat[grid->elem3d[ie].mat].sw->refine_tolerance;
        
        // sets as max error, since only one SW component can be run for now (con & sed can change)
        grid->elem_error[ie] = error;
    }
    
    // average the elemental to nodal errors for output
    for (i=0; i<grid->nnodes; i++) {
        sw->error[i] /= (double)sw->iarray[i];
    }
    sarray_init_int(sw->iarray, grid->nnodes); // reset initialize utility array
}


/***********************************************************/
/***********************************************************/
/***********************************************************/
// averages two nodes for a  new node
void ssw_3d_renumber(SSW_3D *sw, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp, SVECT2D *v2tmp, SVECT *vtmp) {
    // note darray and iarray are just utility arrays, no need to worry about them here
    
    node_renumber_double(max_nnode, sw->displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->old_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->older_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_vect(max_nnode, sw->vel, vtmp, new_numbers, order_tmp);
    node_renumber_vect(max_nnode, sw->old_vel, vtmp, new_numbers, order_tmp);
    node_renumber_vect(max_nnode, sw->older_vel, vtmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, sw->tanvec, v2tmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->density, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->error, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->grid_speed, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->old_grid_speed, dtmp, new_numbers, order_tmp);
    //node_renumber_double(max_nnode, sw->elev_factor, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->prs, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->prs_plus, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->prs_minus, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->vertical_node_flux, dtmp, new_numbers, order_tmp);
#ifdef _SEDIMENT
    node_renumber_double(max_nnode, sw->bed_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->old_bed_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->older_bed_displacement, dtmp, new_numbers, order_tmp);
#endif
    
    
}

/************************************************************/
/************************************************************/
/************************************************************/
// reallocs 3d variables when a new node increment is added
void ssw_3d_realloc_init(SSW_3D *sw, int nnodes_old, int nnodes_new) {
    
    /* re-allocate sw3d 3d variables */
    sw->iarray = (int *) tl_realloc(sizeof(int), nnodes_new, nnodes_old, sw->iarray);
    sw->darray = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->darray);
    sw->displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->displacement);
    sw->old_displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->old_displacement);
    sw->older_displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->older_displacement);
    sw->grid_speed = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->grid_speed);
    sw->old_grid_speed = tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->old_grid_speed);
    sw->density = tl_realloc(sizeof(double), nnodes_new, nnodes_old,sw->density);
    sw->prs = tl_realloc(sizeof(double), nnodes_new, nnodes_old,sw->prs);
    sw->prs_plus = tl_realloc(sizeof(double), nnodes_new, nnodes_old,sw->prs_plus);
    sw->prs_minus = tl_realloc(sizeof(double), nnodes_new, nnodes_old,sw->prs_minus);
    sw->error = tl_realloc(sizeof(double), nnodes_new, nnodes_old,sw->error);
    sw->vertical_node_flux = tl_realloc(sizeof(double), nnodes_new, nnodes_old,sw->vertical_node_flux);
    sw->vel = (SVECT *) tl_realloc(sizeof(SVECT), nnodes_new, nnodes_old,sw->vel);
    sw->old_vel = (SVECT *) tl_realloc(sizeof(SVECT), nnodes_new, nnodes_old, sw->old_vel);
    sw->older_vel = (SVECT *) tl_realloc(sizeof(SVECT), nnodes_new, nnodes_old, sw->older_vel);
    sw->tanvec = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, sw->tanvec);
    sw->hyd_viscosity = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->hyd_viscosity); /* GSAVANT */
    sw->trn_diffusivity = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->trn_diffusivity); /* GSAVANT */
    sw->grad_bed = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, sw->grad_bed);
#ifdef _SEDIMENT
    sw->bed_displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->bed_displacement);
    sw->old_bed_displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->old_bed_displacement);
    sw->older_bed_displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->older_bed_displacement);
#endif
    
    sw->dpl_perturbation = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->dpl_perturbation);
    
    int i;
    for (i=nnodes_old;i<nnodes_new;i++){
        sw->iarray[i] = UNSET_INT;
        sw->darray[i] = 0.;
        sw->displacement[i] = 0.;
        sw->old_displacement[i] = 0.;
        sw->older_displacement[i] = 0.;
        sw->grid_speed[i] = 0.;
        sw->old_grid_speed[i] = 0.;
        sw->density[i] = 0.;
        sw->prs[i] = 0.;
        sw->prs_plus[i] = 0.;
        sw->prs_minus[i] = 0.;
        sw->error[i] = 0.;
        sw->vertical_node_flux[i] = 0.;
        svect_init(&sw->vel[i]);
        svect_init(&sw->old_vel[i]);
        svect_init(&sw->older_vel[i]);
        svect2d_init(&sw->tanvec[i]);
        sw->hyd_viscosity[i] = 0.; /* GSAVANT */
        sw->trn_diffusivity[i] = 0.; /* GSAVANT */
        sw->grad_bed[i].x = 0.;
        sw->grad_bed[i].y = 0.;
#ifdef _SEDIMENT
        sw->bed_displacement[i] = 0.;
        sw->old_bed_displacement[i] = 0. ;
        sw->older_bed_displacement[i] = 0.;
#endif
        sw->dpl_perturbation[i] = 0.;
    }
}
/***********************************************************/
/***********************************************************/
/***********************************************************/
// averages two nodes for a new node - cjt :: since 2d and 3d arrays are numbered different must send both until a map is made
void ssw_3d_node_avg(SSW_3D *sw, int node_new, int node1, int node2, SGRID *grid) {
    
    // 3D :: (cjt) note: darray and iarray are just utility arrays, no need to worry about them here
    sw->displacement[node_new] = 0.5*(sw->displacement[node1] + sw->displacement[node2]);
    sw->old_displacement[node_new] = 0.5*(sw->old_displacement[node1] + sw->old_displacement[node2]);
    sw->older_displacement[node_new] = 0.5*(sw->older_displacement[node1] + sw->older_displacement[node2]);
    sw->grid_speed[node_new] = 0.5*(sw->grid_speed[node1] + sw->grid_speed[node2]);
    sw->old_grid_speed[node_new] = 0.5*(sw->old_grid_speed[node1] + sw->old_grid_speed[node2]);
    sw->density[node_new] = 0.5*(sw->density[node1] + sw->density[node2]);
    sw->prs[node_new] = 0.5*(sw->prs[node1] + sw->prs[node2]);
    sw->prs_plus[node_new] = 0.5*(sw->prs_plus[node1] + sw->prs_plus[node2]);
    sw->prs_minus[node_new] = 0.5*(sw->prs_minus[node1] + sw->prs_minus[node2]);
    sw->error[node_new] = 0.5*(sw->error[node1] + sw->error[node2]);
    sw->vertical_node_flux[node_new] = 0.5*(sw->vertical_node_flux[node1] + sw->vertical_node_flux[node2]);
    sw->vel[node_new] = svect_avg(sw->vel[node1], sw->vel[node2]);
    sw->old_vel[node_new] = svect_avg(sw->old_vel[node1], sw->old_vel[node2]);
    sw->older_vel[node_new] = svect_avg(sw->older_vel[node1], sw->older_vel[node2]);
    sw->tanvec[node_new] = svect2d_avg(sw->tanvec[node1], sw->tanvec[node2]);
    sw->hyd_viscosity[node_new] = 0.5*(sw->hyd_viscosity[node1] + sw->hyd_viscosity[node2]);
    sw->trn_diffusivity[node_new] = 0.5*(sw->trn_diffusivity[node1] + sw->trn_diffusivity[node2]);
#ifdef _SEDIMENT
    sw->bed_displacement[node_new] = 0.5*(sw->bed_displacement[node1] + sw->bed_displacement[node2]);
    sw->old_bed_displacement[node_new] = 0.5*(sw->old_bed_displacement[node1] + sw->old_bed_displacement[node2]);
    sw->older_bed_displacement[node_new] = 0.5*(sw->older_bed_displacement[node1] + sw->older_bed_displacement[node2]);
#endif
    sw->dpl_perturbation[node_new] = 0.5*(sw->dpl_perturbation[node1] + sw->dpl_perturbation[node2]);
    
}

/**********************************************************/
void ssw_3d_node_avg_sur(SSW_3D *sw, int node_new, int node1, int node2, SGRID *grid) {
    
    sw->depth[node_new] = 0.5*(sw->depth[node1] + sw->depth[node2]);
    sw->old_depth[node_new] = 0.5*(sw->old_depth[node1] + sw->old_depth[node2]);
    sw->depth_avg_vel[node_new] = svect2d_avg(sw->depth_avg_vel[node1], sw->depth_avg_vel[node2]);
    sw->old_depth_avg_vel[node_new] = svect2d_avg(sw->old_depth_avg_vel[node1], sw->old_depth_avg_vel[node2]);
    sw->surface_vel[node_new] = svect2d_avg(sw->surface_vel[node1], sw->surface_vel[node2]);
    if(sw->waves!=NULL) {
        sw->waves[node_new].stress.x = 0.5*(sw->waves[node1].stress.x + sw->waves[node2].stress.x);
        sw->waves[node_new].stress.y = 0.5*(sw->waves[node1].stress.y + sw->waves[node2].stress.y);
        sw->waves[node_new].rads.xx = 0.5*(sw->waves[node1].rads.xx + sw->waves[node2].rads.xx);
        sw->waves[node_new].rads.xy = 0.5*(sw->waves[node1].rads.xy + sw->waves[node2].rads.xy);
        sw->waves[node_new].rads.yy = 0.5*(sw->waves[node1].rads.yy + sw->waves[node2].rads.yy);
    }
    if(sw->winds!=NULL) {
        sw->winds[node_new].stress.x = 0.5*(sw->winds[node1].stress.x + sw->winds[node2].stress.x);
        sw->winds[node_new].stress.y = 0.5*(sw->winds[node1].stress.y + sw->winds[node2].stress.y);
    }
    
}
/*********************************************************/

void ssw_3d_node_avg_bed(SSW_3D *sw, int node_new, int node1, int node2, SGRID *grid) {
    sw->bed_elevation[node_new] = 0.5*(sw->bed_elevation[node1] + sw->bed_elevation[node2]);
    sw->bottom_vel[node_new] = svect2d_avg(sw->bottom_vel[node1], sw->bottom_vel[node2]);
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
// reallocs 3d surface/bed variables when a new node increment is added
void ssw_3d_realloc_init_surface(SSW_3D *sw, int nnodes_old, int nnodes_new) {
    sw->depth = tl_realloc(sizeof(double), nnodes_new, nnodes_old,sw->depth);
    sw->old_depth = tl_realloc(sizeof(double), nnodes_new, nnodes_old,sw->old_depth);
    sw->depth_avg_vel = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, sw->depth_avg_vel);
    sw->old_depth_avg_vel = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, sw->old_depth_avg_vel);
    sw->surface_vel = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, sw->surface_vel);
    
    if(sw->winds != NULL) swind_realloc(&sw->winds, nnodes_old, nnodes_new);
    if(sw->waves != NULL) swave_realloc(&sw->waves, nnodes_old, nnodes_new);
    int i;
    for (i=nnodes_old;i<nnodes_new;i++){
        sw->depth[i] = 0.;
        sw->old_depth[i] = 0.;
        svect2d_init(&(sw->depth_avg_vel[i]));
        svect2d_init(&(sw->old_depth_avg_vel[i]));
        svect2d_init(&(sw->surface_vel[i]));
        if(sw->winds != NULL) {
            sw->winds[i].stress.x = 0.;
            sw->winds[i].stress.y = 0.;
        }
        if(sw->waves != NULL) {
            sw->waves[i].stress.x = 0.;
            sw->waves[i].stress.y = 0.;
            sw->waves[i].rads.xx = 0.;
            sw->waves[i].rads.xy = 0.;
            sw->waves[i].rads.yy = 0.;
        }
    }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
// reallocs 3d surface/bed variables when a new node increment is added
void ssw_3d_realloc_init_bed(SSW_3D *sw, int nnodes_old, int nnodes_new) {
    sw->bed_elevation = tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->bed_elevation);
    sw->bottom_vel = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, sw->bottom_vel);
    int i;
    for (i=nnodes_old;i<nnodes_new;i++){
        sw->bed_elevation[i] = 0.;
        sw->bottom_vel[i].x = 0.;
        sw->bottom_vel[i].y = 0.;
    }
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
// renumbers 3d surface/bed variables after full 3d struct is renumbered and columns rebuilt
void ssw_3d_renumber_surface(SSW_3D *sw, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp, SVECT2D *v2tmp) {
    node_renumber_double(max_nnode, sw->depth, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->old_depth, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->bed_elevation, dtmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, sw->depth_avg_vel, v2tmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, sw->old_depth_avg_vel, v2tmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, sw->surface_vel, v2tmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, sw->bottom_vel, v2tmp, new_numbers, order_tmp);
    if(sw->waves != NULL) swave_renumber(sw->waves, max_nnode, new_numbers, order_tmp);
    if(sw->winds != NULL) swind_renumber(sw->winds, max_nnode, new_numbers, order_tmp, v2tmp);
    
}


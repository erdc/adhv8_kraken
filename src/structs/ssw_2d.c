#include "global_header.h"
/***********************************************************/
/***********************************************************/
/***********************************************************/

// note: need nstations here to allocate node_contrib and wind/waves
void ssw_2d_alloc_init(SSW_2D **sw2d, SGRID *grid, SIO *io, SFLAGS flags) {
    int i,inode;
 
#ifdef _DEBUG
    int myid = 0, ierr = UNSET_INT;;
#ifdef _MESSG
    ierr = MPI_Comm_rank(cstorm_comm, &myid);
#endif
    printf("\n---- MYID %d :: MPI_COMM_WORLD_RANK: %d :: 2D SW INTIALIZATION\n", grid->smpi->myid,myid);
#endif

    /* initialize 2d shallow model */
    (*sw2d) = (SSW_2D *) tl_alloc(sizeof(SSW_2D), 1);
    SSW_2D *sw = *sw2d;       // alias

    /* intialize optional variables to NULL */
    sw->winds = NULL;
    sw->waves = NULL;

    /* allocate sw2d variables */
    int nnodes = grid->nnodes;
    sw->head = (double *) tl_alloc(sizeof(double), nnodes);
    sw->old_head = (double *) tl_alloc(sizeof(double), nnodes);
    sw->older_head = (double *) tl_alloc(sizeof(double), nnodes);
    sw->vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    sw->old_vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    sw->older_vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    sw->density = (double *) tl_alloc(sizeof(double), nnodes);
    sw->error = (double *) tl_alloc(sizeof(double), nnodes);
    sw->bed_displacement = (double *) tl_alloc(sizeof(double), nnodes); 
    sw->bed_elevation = (double *) tl_alloc(sizeof(double), nnodes);
    sw->darray = (double *) tl_alloc(sizeof(double), nnodes);
    sw->iarray = (int *) tl_alloc(sizeof(int), nnodes);
    sw->dacontResid = (double *) tl_alloc(sizeof(double), nnodes);

    /* initialize */
    ssw_2d_init(sw, 0, nnodes);

    for (inode=0; inode<nnodes; inode++) {
        sw->bed_elevation[inode] = grid->node[inode].z;
    }

    if (flags.WIND) swind_alloc(&(sw->winds), nnodes); 
    if (flags.WAVE) swave_alloc(&(sw->waves), nnodes);

    /* utility array size */
    sw->vwork = NULL;
    sw->vwork_size = 0;

    sw->elem_rhs_realloc = 0;
    //printf("allocated elem_rhs_dacont_extra_terms : %d\n",grid->nelems2d);
    sw->elem_rhs_dacont_extra_terms = (double **) tl_alloc(sizeof(double *), MAX_NNODES_ON_ELEM2D);
    for (i=0; i<MAX_NNODES_ON_ELEM2D; i++) {
        sw->elem_rhs_dacont_extra_terms[i] = (double *) tl_alloc(sizeof(double), grid->nelems2d);
    }
#ifdef _DEBUG
  printf("\n---- MYID %d :: MPI_COMM_WORLD_RANK: %d :: 2D SW FINALIZING INTIALIZATION\n", grid->smpi->myid,myid);
#endif
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
void ssw_2d_init(SSW_2D *sw, int nnodes_start, int nnodes_end) {
    int inode = UNSET_INT;
    for (inode = nnodes_start; inode < nnodes_end; inode++) {
        sw->head[inode] = 0.;
        sw->old_head[inode] = 0.;
        sw->older_head[inode] = 0.;
        svect2d_init(&(sw->vel[inode]));
        svect2d_init(&(sw->old_vel[inode]));
        svect2d_init(&(sw->older_vel[inode]));
        sw->density[inode] = 0.;
        sw->error[inode] = 0.;
        sw->bed_displacement[inode] = 0.;
        sw->bed_elevation[inode] = 0.;
        sw->darray[inode] = 0.;
        sw->iarray[inode] = 0;
        sw->dacontResid[inode] = 0.;
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// reallocs when a new node increment is added
void ssw_2d_realloc_init(SSW_2D *sw, int nnodes_old, int nnodes_new) {
     
    /* re-allocate sw2d variables */
    sw->head = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->head);
    sw->old_head = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->old_head);
    sw->older_head = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->older_head);
    sw->vel = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, sw->vel);
    sw->old_vel = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, sw->old_vel);
    sw->older_vel = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, sw->older_vel);
    sw->density = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->density);
    sw->error = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->error);
    sw->bed_displacement = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->bed_displacement);
    sw->bed_elevation = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->bed_elevation);
    sw->darray = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->darray);
    sw->iarray = (int *) tl_realloc(sizeof(int), nnodes_new, nnodes_old, sw->iarray);
    sw->dacontResid = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, sw->dacontResid);
    if(sw->winds != NULL) swind_realloc(&sw->winds, nnodes_old, nnodes_new);
    if(sw->waves != NULL) swave_realloc(&sw->waves, nnodes_old, nnodes_new);

    /* initialize new nodes */
    ssw_2d_init(sw, nnodes_old, nnodes_new);
    
    int i;
    if(sw->winds!=NULL || sw->waves!=NULL) {
    	for (i=nnodes_old;i<nnodes_new;i++){
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
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// averages two nodes for a new node
void ssw_2d_node_avg(SSW_2D *sw, int node_new, int node1, int node2) {
    // note darray and iarray are just utility arrays, no need to worry about them here
    sw->head[node_new] = 0.5*(sw->head[node1] + sw->head[node2]);
    sw->old_head[node_new] = 0.5*(sw->old_head[node1] + sw->old_head[node2]);
    sw->older_head[node_new] = 0.5*(sw->older_head[node1] + sw->older_head[node2]);
    sw->vel[node_new] = svect2d_avg(sw->vel[node1], sw->vel[node2]);
    sw->old_vel[node_new] = svect2d_avg(sw->old_vel[node1], sw->old_vel[node2]);
    sw->older_vel[node_new] = svect2d_avg(sw->older_vel[node1], sw->older_vel[node2]);
    sw->density[node_new] = 0.5*(sw->density[node1] + sw->density[node2]);
    sw->error[node_new] = 0.5*(sw->error[node1] + sw->error[node2]);
    sw->bed_displacement[node_new] = 0.5*(sw->bed_displacement[node1] + sw->bed_displacement[node2]);
    sw->bed_elevation[node_new] = 0.5*(sw->bed_elevation[node1] + sw->bed_elevation[node2]);
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

/***********************************************************/
/***********************************************************/
/***********************************************************/
// averages two nodes for a  new node
void ssw_2d_renumber(SSW_2D *sw, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp, SVECT2D *vtmp) {
    // note darray and iarray are just utility arrays, no need to worry about them here
    
    node_renumber_double(max_nnode, sw->head, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->old_head, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->older_head, dtmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, sw->vel, vtmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, sw->old_vel, vtmp, new_numbers, order_tmp);
    node_renumber_vect2d(max_nnode, sw->older_vel, vtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->density, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->error, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->bed_displacement, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, sw->bed_elevation, dtmp, new_numbers, order_tmp);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// check all arrays for INF or NAN
// NOTE :: darray and iarray are just utility arrays, no need to worry about them here
void ssw_2d_checkall(SSW_2D sw, int nnodes) {
    int inode = UNSET_INT;
    for (inode = 0; inode < nnodes; inode++) {
        Is_Double_Inf_or_NaN(sw.head[inode], __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.old_head[inode], __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.older_head[inode], __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.vel[inode].x, __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.vel[inode].y, __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.old_vel[inode].x, __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.old_vel[inode].y, __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.older_vel[inode].x, __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.older_vel[inode].y, __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.density[inode], __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.error[inode], __FILE__, __LINE__);
    }
    if (sw.winds != NULL) {
        Is_Double_Inf_or_NaN(sw.winds[inode].stress.x, __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.winds[inode].stress.y, __FILE__, __LINE__);
    }    
    if (sw.waves != NULL) {
        Is_Double_Inf_or_NaN(sw.waves[inode].stress.x, __FILE__, __LINE__);
        Is_Double_Inf_or_NaN(sw.waves[inode].stress.y, __FILE__, __LINE__);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

/*void ssw_2d_print_ts(SSW_2D *sw, SIO *info, double time, SGRID *grid, STR_VALUE *string) {
    int ie = 0, inode = 0;
    int nnodes = grid->initial_nnodes; // alias - first nnodes are the original ones
    
    // basic hydro
    fprintf(info->fout_sw2_head.fp, "TS 0 %15.8e\n", time);
    sarray_print_dbl(info->fout_sw2_head.fp, sw->head, nnodes);
    fprintf(info->fout_sw2_old_head.fp, "TS 0 %15.8e\n", time);
    sarray_print_dbl(info->fout_sw2_old_head.fp, sw->old_head, nnodes);
    fprintf(info->fout_sw2_vel.fp, "TS 0 %15.8e\n", time);
    svect2d_print_array(info->fout_sw2_vel.fp, sw->vel, nnodes);
    fprintf(info->fout_sw2_old_vel.fp, "TS 0 %15.8e\n", time);
    svect2d_print_array(info->fout_sw2_old_vel.fp, sw->old_vel, nnodes);
    
    // density
    //print_double_array(fout, sw->density, nnodes, time);
    
    // maximum nodal continuity errors (here we want to loop over all (initial + adapted) nnodes)
    // elem2d_to_node_double(grid, string, sw->darray, grid->elem_error); // jacobian weighted
    sarray_init_int(sw->iarray, grid->nnodes);
    sarray_init_dbl(sw->darray, grid->nnodes);
    for (ie=0; ie<grid->nelems2d; ie++) {
        for (inode=0; inode<NDONTRI; inode++) {
            sw->darray[grid->elem2d[ie].nodes[inode]] += grid->elem_error[ie];
            sw->iarray[grid->elem2d[ie].nodes[inode]] += 1;
        }
    }
    for (inode=0; inode<grid->nnodes; inode++) {
        sw->darray[inode] /= (double)sw->iarray[inode];
    }
    fprintf(info->fout_sw2_error.fp, "TS 0 %15.8e\n", time);
    sarray_print_dbl(info->fout_sw2_error.fp, sw->darray, nnodes);
    
    // hydro nodal continuity errors
    fprintf(info->fout_sw2_error_hydro.fp, "TS 0 %15.8e\n", time);
    sarray_print_dbl(info->fout_sw2_error_hydro.fp, sw->error, nnodes);

}
*/
/***********************************************************/
/***********************************************************/
/***********************************************************/


  
/*---------------------------------------------------------*/
/*---------------------------------------------------------*/

void ssw_2d_print_ts(SSW_2D *sw, SIO *info, SGRID *grid, double time, int outfact, SFILE sup, char *proj_name, int it1, int it2, int **ndata, int my_nnode_max, int *my_nnode_ext, int flag, SFILE_OUTPUT file_output) {
    int super = strlen(sup.filename);
    int ie = 0, inode = 0;
    
    SFILE fout;

    if (grid->smpi->myid <= 0) {
      fprintf(info->fout_sw2_head.fp, "TS 0 %15.8e\n", time);
      fprintf(info->fout_sw2_old_head.fp, "TS 0 %15.8e\n", time);
      fprintf(info->fout_sw2_vel.fp, "TS 0 %15.8e\n", time);
      fprintf(info->fout_sw2_old_vel.fp, "TS 0 %15.8e\n", time);
      fprintf(info->fout_sw2_error.fp, "TS 0 %15.8e\n", time);
      fprintf(info->fout_sw2_error_hydro.fp, "TS 0 %15.8e\n", time);
    	init_adh_file(&(fout));     // initialize file
    	if(flag){
    		build_filename2(fout.filename, MAXLINE, proj_name, "_dep.dat-", it1,".", it2); // build the new filename
    		open_output_file(&(fout), "adapted head file", super); // open the file
    		fprintf(fout.fp, "DATASET\n");
    		pdata(proj_name, "", fout.fp, 0, "Depth", "", "BEGSCL", "mesh2d", grid->macro_nnodes, grid->macro_nelems2d); // header
    		tc_timeunits(fout.fp, outfact);
    		fprintf(fout.fp, "TS 0 %15.8e\n", time);    // print the time0-header
      }
    }
    
    sarray_print_dbl(grid, info->fout_sw2_head.fp, fout.fp, sw->head, ndata, my_nnode_max, my_nnode_ext, flag);    // print the variables
    sarray_print_dbl(grid, info->fout_sw2_old_head.fp, fout.fp, sw->old_head, ndata, my_nnode_max, my_nnode_ext, 0);

    if ((grid->smpi->myid <= 0) && flag) print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer
    
    if ((grid->smpi->myid <= 0) && flag) {
      init_adh_file(&(fout));
      build_filename2(fout.filename, MAXLINE, proj_name, "_vel.dat-", it1,".", it2);
      open_output_file(&(fout), "adapted velocity file", super);
      fprintf(fout.fp, "DATASET\n");
      pdata(proj_name, "", fout.fp, 0, "Depth-Averaged Velocity", "", "BEGVEC", "mesh2d", grid->macro_nnodes, grid->macro_nelems2d);
      tc_timeunits(fout.fp, outfact);
      fprintf(fout.fp, "TS 0 %15.8e\n", time);
    }

    svect2d_print_array_MPI(grid, info->fout_sw2_vel.fp, fout.fp, sw->vel, ndata, my_nnode_max, my_nnode_ext, flag);
    svect2d_print_array_MPI(grid, info->fout_sw2_old_vel.fp, fout.fp, sw->old_vel, ndata, my_nnode_max, my_nnode_ext, 0);

    if ((grid->smpi->myid <= 0) && flag) print_trailer(fout.fp);
    
    if ((grid->smpi->myid <= 0) && flag) {
      init_adh_file(&(fout));
      build_filename2(fout.filename, MAXLINE, proj_name, "_error.dat-", it1,".", it2);
      open_output_file(&(fout), "adapted error file", super);
      fprintf(fout.fp, "DATASET\n");
      pdata(proj_name, "", fout.fp, 0, "Error", "", "BEGSCL", "mesh2d", grid->macro_nnodes, grid->macro_nelems2d);
      tc_timeunits(fout.fp, outfact);
      fprintf(fout.fp, "TS 0 %15.8e\n", time);
    }
    // maximum nodal continuity errors (here we want to loop over all (initial + adapted) nnodes)
    // elem2d_to_node_double(grid, string, sw->darray, grid->elem_error); // jacobian weighted
    sarray_init_int(sw->iarray, grid->nnodes);
    sarray_init_dbl(sw->darray, grid->nnodes);
    for (ie=0; ie<grid->nelems2d; ie++) {
        for (inode=0; inode<NDONTRI; inode++) {
            sw->darray[grid->elem2d[ie].nodes[inode]] += grid->elem_error[ie];
            sw->iarray[grid->elem2d[ie].nodes[inode]] += 1;
        }
    }
    for (inode=0; inode<grid->nnodes; inode++) {
        sw->darray[inode] /= (double)sw->iarray[inode];
    }
    
    sarray_print_dbl(grid, info->fout_sw2_error.fp, fout.fp, sw->darray, ndata, my_nnode_max, my_nnode_ext, flag);

    if ((grid->smpi->myid <= 0) && flag) print_trailer(fout.fp);
    
    if ((grid->smpi->myid <= 0) && flag) {
      init_adh_file(&(fout));
      build_filename2(fout.filename, MAXLINE, proj_name, "_error_hydro.dat-", it1,".", it2);
      open_output_file(&(fout), "adapted error hydro file", super);
      fprintf(fout.fp, "DATASET\n");
      pdata(proj_name, "", fout.fp, 0, "Hydro Error", "", "BEGSCL", "mesh2d", grid->macro_nnodes, grid->macro_nelems2d);
      tc_timeunits(fout.fp, outfact);
      fprintf(fout.fp, "TS 0 %15.8e\n", time);
    }
    sarray_print_dbl(grid, info->fout_sw2_error_hydro.fp, fout.fp, sw->error, ndata, my_nnode_max, my_nnode_ext, flag);
    if ((grid->smpi->myid <= 0) && flag)print_trailer(fout.fp);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void ssw_2d_open_output(SMODEL * mod) {

    assert(mod->io);    /* should be valid */
    int super = strlen(mod->io->sup.filename);  /* whether super file exists */

    open_output_file(&(mod->io->fout_sw2_old_head), "sw2d old depth file", super);
    print_header(mod, mod->io->fout_sw2_old_head.fp, PS_FLAG_PDEPTH);

    open_output_file(&(mod->io->fout_sw2_vel), "sw2d velocity file", super);
    print_header(mod, mod->io->fout_sw2_vel.fp, PS_FLAG_OLVEL);

    open_output_file(&(mod->io->fout_sw2_old_vel), "sw2d old velocity file", super);
    print_header(mod, mod->io->fout_sw2_old_vel.fp, PS_FLAG_POLVEL);

    open_output_file(&(mod->io->fout_sw2_error), "sw2d error file", super);
    print_header(mod, mod->io->fout_sw2_error.fp, PS_FLAG_ERR);

    open_output_file(&(mod->io->fout_sw2_error_hydro), "sw2d hydrodynamic error file", super);
    print_header(mod, mod->io->fout_sw2_error_hydro.fp, PS_FLAG_ERR_HYDRO);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void ssw_2d_free(SSW_2D *sw, SGRID *grid, SFLAGS flags) {

    if (sw != NULL) {
        int nnodes = grid->max_nnodes;
        

        /* free sw2 arrays */
#ifdef _DEBUG
        printf("... MYID %d freeing sw2d solution arrays ...\n", grid->smpi->myid);
#endif
        sw->head = (double *) tl_free(sizeof(double), nnodes, sw->head);
        sw->old_head = (double *) tl_free(sizeof(double), nnodes, sw->old_head);
        sw->older_head = (double *) tl_free(sizeof(double), nnodes, sw->older_head);
        sw->density = (double *) tl_free(sizeof(double), nnodes, sw->density);
        sw->error = (double *) tl_free(sizeof(double), nnodes, sw->error);
        sw->vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, sw->vel);
        sw->old_vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, sw->old_vel);
        sw->older_vel = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, sw->older_vel);
        sw->bed_displacement = (double *) tl_free(sizeof(double), nnodes, sw->bed_displacement);
        sw->bed_elevation = (double *) tl_free(sizeof(double), nnodes, sw->bed_elevation);
        sw->darray = (double *) tl_free(sizeof(double), nnodes, sw->darray);
        sw->iarray = (int *) tl_free(sizeof(int), nnodes, sw->iarray);
        sw->dacontResid = (double *) tl_free(sizeof(double), nnodes, sw->dacontResid);
        if (sw->vwork != NULL) sw->vwork = (void *) tl_free(sizeof(char), sw->vwork_size, sw->vwork);
        
        if (flags.WIND) swind_free(sw->winds, nnodes);
        if (flags.WAVE) swave_free(sw->waves, nnodes);


		int i;
		for (i=0; i<NDONQUAD; i++) {		
			sw->elem_rhs_dacont_extra_terms[i] = (double *) tl_free(sizeof(double), grid->nelems2d, sw->elem_rhs_dacont_extra_terms[i]);
		}
		sw->elem_rhs_dacont_extra_terms = (double **) tl_free(sizeof(double *), NDONQUAD,sw->elem_rhs_dacont_extra_terms);

        /* free sw2 struct */
#ifdef _DEBUG
        printf("... MYID %d freeing ssw_2d struct ...\n", grid->smpi->myid);
#endif
        sw = (SSW_2D *) tl_free(sizeof(SSW_2D), 1, sw);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

/* Approximates elemental error for the 2d shallow water equations */
void ssw_2d_calculate_elem_error(SSW_2D *sw, SGRID *grid, SMAT *mat, double dt) {
    
    /* calculates the error indicator for the flow */
    int ie = UNSET_INT, i=UNSET_INT;;
    SVECT2D elem_vel[NDONTRI]; svect2d_init_array(elem_vel, NDONTRI);
    double elem_head[NDONTRI]; sarray_init_dbl(elem_head, NDONTRI);
    double elem_old_head[NDONTRI]; sarray_init_dbl(elem_old_head, NDONTRI);
    double error = 0., dhdx = 0., dhdy = 0., dudx = 0., dvdy = 0.;
    double coef_0 = 0., coef_1 = 0., coef_2 = 0.;
    
    sarray_init_int(sw->iarray, grid->nnodes); // element per node count for averaging
    sarray_init_dbl(sw->error, grid->nnodes);  // nodal error for printing

    int nelems2d = grid->nelems2d; // alias
    for (ie = 0; ie < nelems2d; ie++) {
        error = 0.;
        grid->elem_error[ie] = 0.;

        ELEM2D_GET_LOCAL_VECT2D(sw->vel, elem_vel, grid->elem2d[ie].nodes);
        ELEM2D_GET_LOCAL(sw->head, elem_head, grid->elem2d[ie].nodes);
        ELEM2D_GET_LOCAL(sw->old_head, elem_old_head, grid->elem2d[ie].nodes);
        
        dhdx = 0.; dhdy = 0.;
        dudx = 0.; dvdy = 0.;
        for (i=0; i<NDONTRI; i++) {
            dhdx += grid->elem2d[ie].grad_shp[i].x * elem_head[i];
            dhdy += grid->elem2d[ie].grad_shp[i].y * elem_head[i];
            dudx += grid->elem2d[ie].grad_shp[i].x * elem_vel[i].x;
            dvdy += grid->elem2d[ie].grad_shp[i].y * elem_vel[i].y;
        }
        
        /* steady-state continuity */
        coef_0 = (elem_head[0] - elem_old_head[0]) / dt;
        coef_0 += elem_vel[0].x * dhdx;
        coef_0 += elem_vel[0].y * dhdy;
        coef_0 += elem_head[0] * (dudx + dvdy);
        coef_0 = coef_0 * coef_0;
        
        coef_1 = (elem_head[1] - elem_old_head[1]) / dt;
        coef_1 += elem_vel[1].x * dhdx;
        coef_1 += elem_vel[1].y * dhdy;
        coef_1 += elem_head[1] * (dudx + dvdy);
        coef_1 = coef_1 * coef_1;
        
        coef_2 = (elem_head[2] - elem_old_head[2]) / dt;
        coef_2 += elem_vel[2].x * dhdx;
        coef_2 += elem_vel[2].y * dhdy;
        coef_2 += elem_head[2] * (dudx + dvdy);
        coef_2 = coef_2 * coef_2;
        
        error = sqrt(coef_0 + coef_1 + coef_2) * grid->elem2d[ie].djac;
        for (i=0; i<NDONTRI; i++) {
            sw->error[grid->elem2d[ie].nodes[i]] += error;
            sw->iarray[grid->elem2d[ie].nodes[i]] += 1;
        }

        /* scales by the tolerance */
        error /= mat[grid->elem2d[ie].mat].sw->refine_tolerance;
        
        /* sets as max error, since only one SW component can be run for now (con & sed can change) */
        grid->elem_error[ie] = error;

        /* don't refine or relax any element that is partially wet */
        if (grid->wd_flag[ie] > 0) {
            grid->elem_error[ie] = 0.5;
        }
        
    }

    // average the elemental to nodal errors for output
    for (i=0; i<grid->nnodes; i++) {
        sw->error[i] /= (double)sw->iarray[i];
    }
    sarray_init_int(sw->iarray, grid->nnodes); // reset initialize utility array

}


/* Approximates elemental error for the 2d diffusive water equations */
void ssw_diffusive_calculate_elem_error(SSW_2D *sw, SGRID *grid, SMAT *mat, double dt) {
    
    /* calculates the error indicator for the flow */
    int ie = UNSET_INT, i=UNSET_INT;;
    SVECT2D elem_vel[NDONTRI]; svect2d_init_array(elem_vel, NDONTRI);
    double elem_head[NDONTRI]; sarray_init_dbl(elem_head, NDONTRI);
    double elem_old_head[NDONTRI]; sarray_init_dbl(elem_old_head, NDONTRI);
    double error = 0., dhdx = 0., dhdy = 0., dudx = 0., dvdy = 0.;
	double dzdx, dzdy;
    double coef_0 = 0., coef_1 = 0., coef_2 = 0.;
    
    sarray_init_int(sw->iarray, grid->nnodes); // element per node count for averaging
    sarray_init_dbl(sw->error, grid->nnodes);  // nodal error for printing

    int nelems2d = grid->nelems2d; // alias
    for (ie = 0; ie < nelems2d; ie++) {
        error = 0.;
        grid->elem_error[ie] = 0.;

        ELEM2D_GET_LOCAL_VECT2D(sw->vel, elem_vel, grid->elem2d[ie].nodes);
        ELEM2D_GET_LOCAL(sw->head, elem_head, grid->elem2d[ie].nodes);
        ELEM2D_GET_LOCAL(sw->old_head, elem_old_head, grid->elem2d[ie].nodes);
        
        dhdx = 0.; dhdy = 0.;
        dudx = 0.; dvdy = 0.;
		dzdx = 0.; dzdy = 0.;
        for (i=0; i<NDONTRI; i++) {
            dhdx += grid->elem2d[ie].grad_shp[i].x * elem_head[i];
            dhdy += grid->elem2d[ie].grad_shp[i].y * elem_head[i];
            dudx += grid->elem2d[ie].grad_shp[i].x * elem_vel[i].x;
            dvdy += grid->elem2d[ie].grad_shp[i].y * elem_vel[i].y;
			dzdx += grid->elem2d[ie].grad_shp[i].x * elem_head[i] + grid->elem2d[ie].grad_shp[i].x * grid->node[grid->elem2d[ie].nodes[i]].z * grid->elem2d[ie].grad_shp[i].x;
			dzdy += grid->elem2d[ie].grad_shp[i].y * elem_head[i] + grid->elem2d[ie].grad_shp[i].y * grid->node[grid->elem2d[ie].nodes[i]].z * grid->elem2d[ie].grad_shp[i].y;
			elem_vel[0].x = dzdx * elem_head[0];
			elem_vel[1].x = dzdx * elem_head[0];
			elem_vel[2].x = dzdx * elem_head[0];
			elem_vel[0].y = dzdy * elem_head[0];
			elem_vel[1].y = dzdy * elem_head[0];
			elem_vel[2].y = dzdy * elem_head[0];
        }
        
        /* steady-state continuity */
        coef_0 = (elem_head[0] - elem_old_head[0]) / dt;
        coef_0 += elem_vel[0].x * dhdx;
        coef_0 += elem_vel[0].y * dhdy;
        coef_0 += elem_head[0] * (dudx + dvdy);
        coef_0 = coef_0 * coef_0;
        
        coef_1 = (elem_head[1] - elem_old_head[1]) / dt;
        coef_1 += elem_vel[1].x * dhdx;
        coef_1 += elem_vel[1].y * dhdy;
        coef_1 += elem_head[1] * (dudx + dvdy);
        coef_1 = coef_1 * coef_1;
        
        coef_2 = (elem_head[2] - elem_old_head[2]) / dt;
        coef_2 += elem_vel[2].x * dhdx;
        coef_2 += elem_vel[2].y * dhdy;
        coef_2 += elem_head[2] * (dudx + dvdy);
        coef_2 = coef_2 * coef_2;
        
        error = sqrt(coef_0 + coef_1 + coef_2) * grid->elem2d[ie].djac;
        for (i=0; i<NDONTRI; i++) {
            sw->error[grid->elem2d[ie].nodes[i]] += error;
            sw->iarray[grid->elem2d[ie].nodes[i]] += 1;
        }

        /* scales by the tolerance */
        error /= mat[grid->elem2d[ie].mat].sw->refine_tolerance;
        
        /* sets as max error, since only one SW component can be run for now (con & sed can change) */
        grid->elem_error[ie] = error;


        
    }

    // average the elemental to nodal errors for output
    for (i=0; i<grid->nnodes; i++) {
        sw->error[i] /= (double)sw->iarray[i];
    }
    sarray_init_int(sw->iarray, grid->nnodes); // reset initialize utility array

}

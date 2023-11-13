#include "global_header.h"
static int ncon_properties = 4; // # of properties per constituent

/******************************************************************************/
/******************************************************************************/
// allocate and initialize all variables in SCON struct
void scon_alloc_init(int nnodes, int nelems, int ncon, SCON **constituent) {
    
    /* allocate the sconstinuent struct */
    (*constituent) = (SCON *) tl_alloc(sizeof(SCON), ncon);
    SCON *con = (*constituent);         // alias
    
    int icon;
    for (icon=0; icon<ncon; icon++) {
        con[icon].concentration = (double *) tl_alloc(sizeof(double), nnodes);
        con[icon].old_concentration = (double *) tl_alloc(sizeof(double), nnodes);
        con[icon].older_concentration = (double *) tl_alloc(sizeof(double), nnodes);
        con[icon].sink = (double *) tl_alloc(sizeof(double), nnodes);
        con[icon].source = (double *) tl_alloc(sizeof(double), nnodes);
        con[icon].nodal_decay_coefficient = (double *) tl_alloc(sizeof(double), nnodes);
        con[icon].error = (double *) tl_alloc(sizeof(double), nnodes);
        con[icon].mfcf = (double *) tl_alloc(sizeof(double), nnodes);
        con[icon].vcf = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);

        con[icon].property = (double *) tl_alloc(sizeof(double), ncon_properties);
        sarray_init_dbl(con[icon].property, ncon_properties);
    }

    scon_init(con, ncon, 0, nnodes);
}

/******************************************************************************/
/******************************************************************************/
// initialize
void scon_init(SCON *con, int ncon, int nnodes_start, int nnodes_end) {

    int icon = UNSET_INT;
    int inode = UNSET_INT;
    for (icon=0; icon<ncon; icon++) {
        for (inode = nnodes_start; inode < nnodes_end; inode++) {
            con[icon].concentration[inode] = 0.;
            con[icon].old_concentration[inode] = 0.;
            con[icon].older_concentration[inode] = 0.;
            con[icon].sink[inode] = 0.;
            con[icon].source[inode] = 0.;
            con[icon].nodal_decay_coefficient[inode] = 0.;
            con[icon].error[inode] = 0.;
            con[icon].mfcf[inode] = 1.;
            svect2d_init(&(con[icon].vcf[inode]));
        }
    }
}
/******************************************************************************/
/******************************************************************************/
// allocate and initialize all variables in SCON struct
void scon_realloc_init(SCON *con, int ncon, int nnodes_old, int nnodes_new) {

    int icon;
    for (icon=0; icon<ncon; icon++) {
        con[icon].concentration = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, con[icon].concentration);
        con[icon].old_concentration = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, con[icon].old_concentration);
        con[icon].older_concentration = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, con[icon].older_concentration);
        con[icon].sink = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, con[icon].sink);
        con[icon].source = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, con[icon].source);
        con[icon].nodal_decay_coefficient = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, con[icon].nodal_decay_coefficient);
        con[icon].error = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, con[icon].error);
        con[icon].mfcf = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, con[icon].mfcf);
        con[icon].vcf = (SVECT2D *) tl_realloc(sizeof(SVECT2D), nnodes_new, nnodes_old, con[icon].vcf);
    }
    scon_init(con, ncon, nnodes_old, nnodes_new);
}

/******************************************************************************/
/******************************************************************************/

// averages two nodes for a  new node
void scon_node_avg(SCON *con, int ncon, int node_new, int node1, int node2) {
    int icon = UNSET_INT;
    for (icon=0; icon<ncon; icon++) {
        con[icon].concentration[node_new] = 0.5*(con[icon].concentration[node1] + con[icon].concentration[node2]);
        con[icon].old_concentration[node_new] = 0.5*(con[icon].old_concentration[node1] + con[icon].old_concentration[node2]);
        con[icon].older_concentration[node_new] = 0.5*(con[icon].older_concentration[node1] + con[icon].older_concentration[node2]);        
        con[icon].sink[node_new] = 0.5*(con[icon].sink[node1] + con[icon].sink[node2]);
        con[icon].source[node_new] = 0.5*(con[icon].source[node1] + con[icon].source[node2]);
        con[icon].nodal_decay_coefficient[node_new] = 0.5*(con[icon].nodal_decay_coefficient[node1] + con[icon].nodal_decay_coefficient[node2]);
        con[icon].error[node_new] = 0.5*(con[icon].error[node1] + con[icon].error[node2]);
        con[icon].mfcf[node_new] = 0.5*(con[icon].mfcf[node1] + con[icon].mfcf[node2]);
        con[icon].vcf[node_new] = svect2d_avg(con[icon].vcf[node1],con[icon].vcf[node2]);
    }
}

/******************************************************************************/
/******************************************************************************/
// renumbers after a node is added
void scon_renumber(SCON *con, int ncon, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp, SVECT2D *vtmp) {
    int icon = UNSET_INT;
    for (icon=0; icon<ncon; icon++) {
        node_renumber_double(max_nnode, con[icon].concentration, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, con[icon].old_concentration, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, con[icon].older_concentration, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, con[icon].sink, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, con[icon].source, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, con[icon].nodal_decay_coefficient, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, con[icon].error, dtmp, new_numbers, order_tmp);
        node_renumber_double(max_nnode, con[icon].mfcf, dtmp, new_numbers, order_tmp);
        node_renumber_vect2d(max_nnode, con[icon].vcf, vtmp, new_numbers, order_tmp);
    }
}

/******************************************************************************/
/******************************************************************************/
void scon_free(int nnodes, int nelems, int ncon, SCON *con) {
    int icon;
    for (icon=0; icon<ncon; icon++) {
        con[icon].property = (double *) tl_free(sizeof(double), ncon_properties, con[icon].property);
        con[icon].concentration = (double *) tl_free(sizeof(double), nnodes, con[icon].concentration);
        con[icon].old_concentration = (double *) tl_free(sizeof(double), nnodes, con[icon].old_concentration);
        con[icon].older_concentration = (double *) tl_free(sizeof(double), nnodes, con[icon].older_concentration);
        con[icon].sink = (double *) tl_free(sizeof(double), nnodes, con[icon].sink);
        con[icon].source = (double *) tl_free(sizeof(double), nnodes, con[icon].source);
        con[icon].nodal_decay_coefficient = (double *) tl_free(sizeof(double), nnodes, con[icon].nodal_decay_coefficient);
        con[icon].error = (double *) tl_free(sizeof(double), nnodes, con[icon].error);
        con[icon].mfcf = (double *) tl_free(sizeof(double), nnodes, con[icon].mfcf);
        con[icon].vcf = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, con[icon].vcf);
    }
    con = (SCON *) tl_free(sizeof(SCON), ncon, con);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/


/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
/***********************************************************/
/***********************************************************/
/***********************************************************/
/*void scon_print_ts_MPI(int ntransport, double outfact, SIO *info, SCON *con, SGRID *grid, double time, SFLAGS flag, STR_VALUE *string) {

     double c_time = time * outfact;
#ifdef _MESSG
    int nnodes = 0;
    int nelems = 0;
    nnodes = grid->initial_nnodes;  // first nnodes are the original ones
    if (flag.SW2_FLOW) {
        nelems = grid->nelems2d;
    } else if (flag.SW3_FLOW) {
        nelems = grid->nelems3d;
    } else {
        return;
    }

    int itrns;
    double *darray = (double *) tl_alloc(sizeof(double), nnodes);
    for (itrns = 0; itrns < ntransport; itrns++) {

        // current concentrations
        sarray_scale(darray, con[itrns].concentration, con[itrns].property[0], nnodes);
        if (grid->smpi->myid <= 0){
          fprintf(info->fout_con[itrns].fp, "TS 0 %15.8e\n", time);
          fprintf(info->fout_error_con[itrns].fp, "TS 0 %15.8e\n", time);
        }
        sarray_print_dbl_MPI(grid, info->fout_con[itrns].fp, darray, 0);

        // previous concentrations
        //sarray_scale(pconc_dim, pconc, con[itrns].property[0], nnodes);
        //print_double_array(info->fout_pcon[itrns].fp, pconc_dim, nnodes, time);

        // concentration error estimator
        sarray_print_dbl_MPI(grid, info->fout_error_con[itrns].fp, con[itrns].error, 0);
    }

    trns = 0; itrns < ntransport; itrns++) {

        // current concentrations
        sarray_scale(darray, con[itrns].concentration, con[itrns].property[0], nnodes);

        fprintf(info->fout_con[itrns].fp, "TS 0 %15.8e\n", time);
        sarray_print_dbl(info->fout_con[itrns].fp, darray, nnodes);

        // previous concentrations
        //sarray_scale(pconc_dim, pconc, con[itrns].property[0], nnodes);
        //print_double_array(info->fout_pcon[itrns].fp, pconc_dim, nnodes, time);

        // concentration error estimator
        fprintf(info->fout_error_con[itrns].fp, "TS 0 %15.8e\n", time);
        sarray_print_dbl(info->fout_error_con[itrns].fp, con[itrns].error, nnodes);
    }

    if (ntransport > 0) {
        darray = (double *) tl_free(sizeof(double), nnodes, darray);
    }if (ntransport > 0) {
        darray = (double *) tl_free(sizeof(double), nnodes, darray);
    }
#endif
}
*/
/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
void scon_print_ts(int ntransport, SCON *con, SIO *info, SGRID *grid, double time, int outfact, SFLAGS flag, SFILE sup, char *proj_name, int it1, int it2, int **ndata, int my_nnode_max, int *my_nnode_ext, int flag1) {
    
    int super = strlen(sup.filename);
    char *mesh_tag;
    int nnodes = 0;
    int nelems = 0;
    nnodes = grid->nnodes;  // first nnodes are the original ones
    if (flag.SW2_FLOW || flag.NS2_FLOW) {
        mesh_tag = "mesh2d";
        nelems = grid->macro_nelems2d;
    } else if (flag.SW3_FLOW || flag.NS3_FLOW) {
        mesh_tag = "mesh3d";
        nelems = grid->macro_nelems3d;
    } 
    SFILE fout;
    
    int itrns;
    double *darray = (double *) tl_alloc(sizeof(double), nnodes);
    for (itrns = 0; itrns < ntransport; itrns++) {
      // current concentrations
      sarray_scale_dbl(darray, con[itrns].concentration, con[itrns].property[0], nnodes);
      if (grid->smpi->myid <= 0) {
        
				fprintf(info->fout_con[itrns].fp, "TS 0 %15.8e\n", time);
        if(flag1){
				  init_adh_file(&(fout));     // initialize file
         build_filename3(fout.filename, MAXLINE, proj_name, "_con", itrns+1,".dat-", it1,".", it2); // build the new filename
          open_output_file(&(fout), "adapted concentration file", super); // open the file
          fprintf(fout.fp, "DATASET\n");
          pdata(proj_name, "", fout.fp, 0, "Concentration", "", "BEGSCL", mesh_tag, grid->macro_nnodes, nelems); // header
          tc_timeunits(fout.fp, outfact);
          fprintf(fout.fp, "TS 0 %15.8e\n", time);    // print the time0-header
        }
      }
      sarray_print_dbl(grid, info->fout_con[itrns].fp, fout.fp, darray, ndata, my_nnode_max, my_nnode_ext, flag1);
      if ((grid->smpi->myid <= 0)&&flag1) print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer

      // concentration error estimator
      if (grid->smpi->myid <= 0) {
        fprintf(info->fout_error_con[itrns].fp, "TS 0 %15.8e\n", time);
        if(flag1){
          init_adh_file(&(fout));     // initialize file
          build_filename3(fout.filename, MAXLINE, proj_name, "_con", itrns+1,"_error.dat-", it1,".", it2); // build the new filename
          open_output_file(&(fout), "adapted concentration error file", super); // open the file
          fprintf(fout.fp, "DATASET\n");
          pdata(proj_name, "", fout.fp, 0, "Concentration Error", "", "BEGSCL", mesh_tag, grid->macro_nnodes, nelems); // header
          tc_timeunits(fout.fp, outfact);
          fprintf(fout.fp, "TS 0 %15.8e\n", time);    // print the time0-header
        }
      }
      sarray_print_dbl(grid, info->fout_error_con[itrns].fp, fout.fp, con[itrns].error, ndata, my_nnode_max, my_nnode_ext, flag1);
      if ((grid->smpi->myid <= 0)&&flag1) print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer
    }

    darray = (double *) tl_free(sizeof(double), nnodes, darray);    
}


/***********************************************************/
/***********************************************************/
/***********************************************************/

void scon_open_output(SMODEL *mod) {
    
    assert(mod->io); /* should be valid */
    int super = strlen(mod->io->sup.filename); /* whether super file exists */
    
    for (mod->itrns=0; mod->itrns<mod->ntransport; mod->itrns++) {
        open_output_file( &(mod->io->fout_con[mod->itrns]), "transport constituent file", super);
        print_header(mod, mod->io->fout_con[mod->itrns].fp, PS_FLAG_CON);

        open_output_file( &(mod->io->fout_error_con[mod->itrns]), "transport constituent error file", super);
        print_header(mod, mod->io->fout_error_con[mod->itrns].fp, PS_FLAG_ERR_CON);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

/* error indicator for the transport equations */
void scon_calculate_elem2d_error(int ntransport, SCON *con, SVECT2D *vel, double *head, SGRID *grid, SMAT *mat, double dt) {
    
    int itrn = UNSET_INT;
    int nelems2d = grid->nelems2d; // alias
    
    /* calculates the error indicator for the flow */
    int ie = UNSET_INT, i = UNSET_INT;
    SVECT2D elem_vel[NDONTRI]; svect2d_init_array(elem_vel, NDONTRI);
    double elem_head[NDONTRI]; sarray_init_dbl(elem_head, NDONTRI);
    double elem_c[NDONTRI]; sarray_init_dbl(elem_c, NDONTRI);
    double elem_old_c[NDONTRI]; sarray_init_dbl(elem_old_c, NDONTRI);
    double error = 0., dcdx = 0., dcdy = 0., dudx = 0., dvdy = 0.;
    double coef_0 = 0., coef_1 = 0., coef_2 = 0.;
    
    int iarray[grid->nnodes]; // element per node count for averaging
    sarray_init_int(iarray, grid->nnodes);
    for (itrn = 0; itrn < ntransport; itrn++) {
        sarray_init_dbl(con[itrn].error, grid->nnodes); // nodal error
    }
    
    for (ie = 0; ie < nelems2d; ie++) {
        ELEM2D_GET_LOCAL(head, elem_head, grid->elem2d[ie].nodes);
        ELEM2D_GET_LOCAL_VECT2D(vel, elem_vel, grid->elem2d[ie].nodes);
        dudx = 0.; dvdy = 0.;
        for (i=0; i<NDONTRI; i++) {
           iarray[grid->elem2d[ie].nodes[i]] += 1;
           dudx += grid->elem2d[ie].grad_shp[i].x * elem_vel[i].x;
           dvdy += grid->elem2d[ie].grad_shp[i].y * elem_vel[i].y;
        } 

        for (itrn = 0; itrn < ntransport; itrn++) {
            error = 0;
            ELEM2D_GET_LOCAL(con[itrn].concentration, elem_c, grid->elem2d[ie].nodes);
            ELEM2D_GET_LOCAL(con[itrn].old_concentration, elem_old_c, grid->elem2d[ie].nodes);
            
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
                con[itrn].error[grid->elem2d[ie].nodes[i]] += error;
            }
            
            /* scales the tolerance */
            error /= mat[grid->elem2d[ie].mat].trn[itrn].refine_tolerance;
            
            /* resets the error if larger */
            if (error > grid->elem_error[ie]) {
                grid->elem_error[ie] = error;
            }
        }
    }
    
    // average the elemental to nodal errors for output
    for (i=0; i<grid->nnodes; i++) {
        for (itrn = 0; itrn < ntransport; itrn++) {
            con[itrn].error[i] /= iarray[i];
        }
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

/* error indicator for the transport equations */
void scon_calculate_elem3d_error(int ntransport, SCON *con, SVECT *vel, double *dpl, double *old_dpl, double *older_dpl, SGRID *grid, SMAT *mat, double dt) {
    
    int itrn = UNSET_INT;
    int nelems3d = grid->nelems3d; // alias
    
    /* calculates the error indicator for the flow */
    int ie = UNSET_INT, i = UNSET_INT;
    double error = 0.;
    
    int iarray[grid->nnodes]; // element per node count for averaging
    sarray_init_int(iarray, grid->nnodes);
    for (itrn = 0; itrn < ntransport; itrn++) {
        sarray_init_dbl(con[itrn].error, grid->nnodes); // nodal error
    }
    
    
    SQUAD *quad;
    int isElementTetrahedron = TRUE, nnodes = 0, nnodes_quad = 0;
    double elem_volume = 0., elem_volume_old = 0., elem_volume_older = 0., elem_djac = 0., elem_djac_old = 0., elem_djac_older = 0.;
    double elem_avg_strong_resid = 0.;

    for (ie = 0; ie < nelems3d; ie++) {
        
        SELEM_3D *elem3d = &(grid->elem3d[ie]);
        
        nnodes =elem3d->nnodes;
        nnodes_quad =elem3d->nnodes_quad;
        if (nnodes == NDONPRISM) isElementTetrahedron = FALSE;
        SNODE elem_nodes[nnodes];
        for (i=0; i<nnodes; i++) {snode_copy(&(elem_nodes[i]), grid->node[elem3d->nodes[i]]);}
        
        // displacements
        double elem_dpl[nnodes], elem_dpl_old[nnodes], elem_dpl_older[nnodes];
        global_to_local_dbl(dpl,       elem_dpl,      elem3d->nodes, nnodes);
        global_to_local_dbl(old_dpl,   elem_dpl_old,  elem3d->nodes, nnodes);
        global_to_local_dbl(older_dpl, elem_dpl_older,elem3d->nodes, nnodes);
        
        // nodal vectors with displacements
        SVECT elem_nds[nnodes_quad], elem_nds_old[nnodes_quad], elem_nds_older[nnodes_quad];
        elem_get_midpt_locations2(elem_nodes, elem_dpl,       elem_nds,       nnodes, nnodes_quad, elem3d->edges);
        elem_get_midpt_locations2(elem_nodes, elem_dpl_old,   elem_nds_old,   nnodes, nnodes_quad, elem3d->edges);
        elem_get_midpt_locations2(elem_nodes, elem_dpl_older, elem_nds_older, nnodes, nnodes_quad, elem3d->edges);
        
        SVECT grad_phi[nnodes], grad_phi_old[nnodes], grad_phi_older[nnodes];
        if (isElementTetrahedron == TRUE) { // djacs cancel here (djac = volume);
            elem_djac = get_tet_linear_djac_gradPhi2(NULL, elem_nds, grad_phi);
            elem_djac_old = get_tet_linear_djac_gradPhi2(NULL, elem_nds_old, grad_phi_old);
            elem_djac_older = get_tet_linear_djac_gradPhi2(NULL, elem_nds_older, grad_phi_older);
            elem_volume = elem_djac;
            
#ifdef _DEBUG
            for (i=0;i<nnodes;i++) {
                Is_Double_Inf_or_NaN(elem_nds[i].x,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(elem_nds[i].y,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(elem_nds[i].z,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(grad_phi[i].x,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(grad_phi[i].y,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(grad_phi[i].z,__FILE__,__LINE__);
            }
            
            Is_Double_Inf_or_NaN(elem_djac,__FILE__,__LINE__);
            Is_Double_Inf_or_NaN(elem_djac_old,__FILE__,__LINE__);
            Is_Double_Inf_or_NaN(elem_djac_older,__FILE__,__LINE__);
            Is_Double_Inf_or_NaN(elem_volume,__FILE__,__LINE__);
#endif
        } else {
            elem_volume = get_triprism_volume(elem_nds);
            elem_volume_old = get_triprism_volume(elem_nds_old);
            elem_volume_older = get_triprism_volume(elem_nds_older);
        }
        
        
        SVECT elem_vel[nnodes];
        ELEM3D_GET_LOCAL_VECT(vel, elem_vel,elem3d->nodes);
        
        // grid and relative velocities
        SVECT elem_rel_vel[nnodes], elem_grid_vel[nnodes];
        svect_init_array(elem_grid_vel, nnodes);
        for (i=0; i<nnodes; i++) {elem_grid_vel[i].z = (elem_dpl[i] - elem_dpl_old[i])/dt;}
        svect_subtract_array2(elem_rel_vel, elem_vel, elem_grid_vel, nnodes);
        double elem_ur[nnodes], elem_vr[nnodes], elem_wr[nnodes];
        dumpVector(elem_rel_vel, nnodes, elem_ur, elem_vr, elem_wr);
        
        // weights
        for (i=0; i<nnodes; i++) {
            iarray[grid->elem3d[ie].nodes[i]] += 1;
        }
        
        for (itrn = 0; itrn < ntransport; itrn++) {
            
            // constituent concentrations
            double elem_c[nnodes], elem_c_old[nnodes], elem_c_older[nnodes];
            global_to_local_dbl(con[itrn].concentration, elem_c, elem3d->nodes, nnodes);
            global_to_local_dbl(con[itrn].old_concentration, elem_c_old, elem3d->nodes, nnodes);
            global_to_local_dbl(con[itrn].older_concentration, elem_c_older, elem3d->nodes, nnodes);
            
            // calculate tetrehedral jacobians and elemental averages
            SVECT elem_grad_u, elem_grad_v, elem_grad_w, elem_grad_c, elem_grad_density;
            double elem_avg_c = 0., elem_avg_c_old = 0., elem_avg_c_older = 0.;
            double elem_avg_ur = 0., elem_avg_vr = 0., elem_avg_wr = 0.;
            SVECT elem_avg_grad_c;
            if (isElementTetrahedron == TRUE) { // djacs cancel here (djac = volume);
                quad = grid->quad_tet;
                
                grad_phi_dot_v(grad_phi, elem_vel, &elem_grad_u, &elem_grad_v, &elem_grad_w, nnodes);
                grad_phi_f(grad_phi, elem_c, &elem_grad_c, nnodes);
                
                elem_avg_c = integrate_tetrahedron_f(1.,1.,elem_c);
                elem_avg_c_old = integrate_tetrahedron_f(1.,1.,elem_c_old);
                elem_avg_c_older = integrate_tetrahedron_f(1.,1.,elem_c_older);
                elem_avg_ur = integrate_tetrahedron_f(1.,1.,elem_ur);
                elem_avg_vr = integrate_tetrahedron_f(1.,1.,elem_vr);
                elem_avg_wr = integrate_tetrahedron_f(1.,1.,elem_wr);
                elem_avg_grad_c = elem_grad_c;
#ifdef _DEBUG
                Is_Double_Inf_or_NaN(elem_avg_c,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(elem_avg_c_old,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(elem_avg_c_older,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(elem_avg_ur,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(elem_avg_vr,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(elem_avg_wr,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(elem_avg_grad_c.x,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(elem_avg_grad_c.y,__FILE__,__LINE__);
                Is_Double_Inf_or_NaN(elem_avg_grad_c.z,__FILE__,__LINE__);
#endif
            } else {
                quad = grid->quad_prism;
                
                elem_avg_c =             integrate_triPrism_f(elem_nds, 1./elem_volume, elem_c);
                elem_avg_c_old =         integrate_triPrism_f(elem_nds_old, 1./elem_volume_old, elem_c_old);
                elem_avg_c_older =       integrate_triPrism_f(elem_nds_older, 1./elem_volume_older, elem_c_older);
                elem_avg_ur =            integrate_triPrism_f(elem_nds, 1./elem_volume, elem_ur);
                elem_avg_vr =            integrate_triPrism_f(elem_nds, 1./elem_volume, elem_vr);
                elem_avg_wr =            integrate_triPrism_f(elem_nds, 1./elem_volume, elem_wr);
                elem_avg_grad_c =        integrate_triPrism_df_full(elem_nds, 1./elem_volume, elem_c);
            }
            
            SVECT elem_avg_vel_rel;
            elem_avg_vel_rel.x = elem_avg_ur;
            elem_avg_vel_rel.y = elem_avg_vr;
            elem_avg_vel_rel.z = elem_avg_wr;
            double elem_avg_c_32dt = get_second_order(elem_avg_c, elem_avg_c_old);
            double elem_avg_c_12dt = get_second_order(elem_avg_c_old, elem_avg_c_older);
            elem_avg_strong_resid = (elem_avg_c_32dt - elem_avg_c_12dt)/dt + elem_avg_grad_c.x * elem_avg_ur + elem_avg_grad_c.y * elem_avg_vr + elem_avg_grad_c.z * elem_avg_wr;
            error = fabs(elem_avg_strong_resid * elem_volume);
            
            // OVER-RIDE!
            //error = 1000;
            
            // nodal error calculation
            for (i=0; i<nnodes; i++) {
                con[itrn].error[grid->elem3d[ie].nodes[i]] += error;
            }
            
            /* scales the tolerance */
            error /= mat[grid->elem3d[ie].mat].trn[itrn].refine_tolerance;
            
            //printf("error: %20.10f\n",elem_volume);
            
            /* resets the error if larger */
            if (error >elem3d->error) {
               elem3d->error = error;
            }
            
            /* if no hydro is being run, use transport error for adaption */
            if (debug.no_hydro == ON) {
                grid->elem_error[ie] = error;
            }
        }
    }
    
    // average the elemental to nodal errors for output
    for (i=0; i<grid->nnodes; i++) {
        for (itrn = 0; itrn < ntransport; itrn++) {
            con[itrn].error[i] /= iarray[i];
        }
    }
}



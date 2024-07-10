#include "global_header.h"

/***********************************************************/
/***********************************************************/
/***********************************************************/
#ifdef _ADH_GROUNDWATER

void sgw_alloc_init(SGW **sgw, SGRID *grid, SIO *io, SFLAGS flags) {
    int i,inode;
    
#ifdef _DEBUG
    printf("-- MYID %d GW initialization\n", grid->smpi->myid);
#endif
    
    /* initialize 2d shallow model */
    (*sgw) = (SGW *) tl_alloc(sizeof(SGW), 1);
    SGW *gw = *sgw;       // alias
    
    /* intialize optional variables to NULL */
    
    /* allocate gw variables */
    int nnodes = grid->nnodes;
    gw->gw_phead = (double *) tl_alloc(sizeof(double), nnodes);
    gw->old_gw_phead = (double *) tl_alloc(sizeof(double), nnodes);
    gw->older_gw_phead = (double *) tl_alloc(sizeof(double), nnodes);
    gw->predict_gw_phead = (double *) tl_alloc(sizeof(double), nnodes);
    gw->gw_density = (double *) tl_alloc(sizeof(double), nnodes);
    gw->elem_3d_data = (ELEMENT_3D_DATA *) tl_alloc(sizeof(ELEMENT_3D_DATA),grid->nelems3d);
    gw->elem_gw_flux= (SVECT *) tl_alloc(sizeof(SVECT),grid->nelems3d);
    gw->error = (double *) tl_alloc(sizeof(double), nnodes);
    gw->darray = (double *) tl_alloc(sizeof(double), nnodes);
    gw->iarray = (int *) tl_alloc(sizeof(int), nnodes);
    
    /* initialize */
    sgw_init(gw, 0, nnodes, 0, grid->nelems3d);
    
    
    /* utility array size */
    gw->vwork = NULL;
    gw->vwork_size = 0;
}
/***********************************************************/
/***********************************************************/
/***********************************************************/
void sgw_init(SGW *gw, int nnodes_start, int nnodes_end, int nelems_start, int nelems_end) {
    int inode = UNSET_INT;
    int ielem = UNSET_INT, ii= UNSET_INT;
    for (inode = nnodes_start; inode < nnodes_end; inode++) {
        gw->gw_phead[inode] = 0.;
        gw->old_gw_phead[inode] = 0.;
        gw->older_gw_phead[inode] = 0.;
        gw->predict_gw_phead[inode] = 0.;
        gw->gw_density[inode] = 0.;
        gw->error[inode] = 0.;
        gw->darray[inode] = 0.;
        gw->iarray[inode] = 0;
    }
    for (ielem = nelems_start; ielem < nelems_end; ielem++) {
        gw->elem_gw_flux[ielem].x = 0.;
        gw->elem_gw_flux[ielem].y = 0.;
        gw->elem_gw_flux[ielem].z = 0.;
        
        for (ii=0; ii < MAX_NNODES_ON_ELEM3D; ii++) {
            gw->elem_3d_data[ielem].saturation[ii] = 0.;
            gw->elem_3d_data[ielem].old_saturation[ii] = 0.;
        }
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sgw_free(SGW *gw, SGRID *grid, SFLAGS flags) {
    int DEBUG=0;
    if (gw != NULL) {
        assert(grid);
        int nnodes;
#ifdef _DEBUG
        if (DEBUG) {
            tl_check_all_pickets(__FILE__, __LINE__);
        }
#endif
        /* free 3d grid solution arrays */
#ifdef _DEBUG
        printf("... MYID %d freeing GW grid solution arrays ...\n", grid->smpi->myid);
#endif
        nnodes = grid->max_nnodes;
        gw->gw_phead = (double *) tl_free(sizeof(double), nnodes, gw->gw_phead);
        gw->old_gw_phead = (double *) tl_free(sizeof(double), nnodes, gw->old_gw_phead);
        gw->older_gw_phead = (double *) tl_free(sizeof(double), nnodes, gw->older_gw_phead);
        gw->predict_gw_phead = (double *) tl_free(sizeof(double), nnodes, gw->predict_gw_phead);
        gw->gw_density = (double *) tl_free(sizeof(double), nnodes, gw->gw_density);
        gw->error = (double *) tl_free(sizeof(double), nnodes, gw->error);
        gw->iarray = (int *) tl_free(sizeof(int), nnodes, gw->iarray);
        gw->darray = (double *) tl_free(sizeof(double), nnodes, gw->darray);
        gw->elem_3d_data = (ELEMENT_3D_DATA *) tl_free(sizeof(ELEMENT_3D_DATA), grid->nelems3d, gw->elem_3d_data);
        gw->elem_gw_flux = (SVECT *) tl_free(sizeof(SVECT), grid->nelems3d, gw->elem_gw_flux);
        gw = (SGW *) tl_free(sizeof(SGW), 1, gw);
        
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sgw_open_output(SMODEL *mod) {
    
    assert(mod->io); /* should be valid */
    int super = strlen(mod->io->sup.filename); /* whether super file exists */
    
    //
    open_output_file( &(mod->io->fout_gw_phead), "pressure head file", super);
    print_header(mod, mod->io->fout_gw_phead.fp, PS_FLAG_PHEAD );
    
    open_output_file( &(mod->io->fout_gw_thead), "total head file", super);
    print_header(mod, mod->io->fout_gw_thead.fp, PS_FLAG_THEAD );
    
    open_output_file( &(mod->io->fout_gw_density), "density file", super);
    print_header(mod, mod->io->fout_gw_density.fp, PS_FLAG_DENS );
    
    open_output_file( &(mod->io->fout_gw_sat), "saturation file", super);
    print_header(mod, mod->io->fout_gw_sat.fp, PS_FLAG_SAT );
    
    open_output_file( &(mod->io->fout_gw_flx), "nodal flux file", super);
    print_header(mod, mod->io->fout_gw_flx.fp, PS_FLAG_FLX );
    
    open_output_file( &(mod->io->fout_gw_vel), "groundwater velocity file", super);
    print_header(mod, mod->io->fout_gw_vel.fp, PS_FLAG_VEL );
    
    
    open_output_file( &(mod->io->fout_gw_error), "gw error file", super);
    print_header(mod, mod->io->fout_gw_error.fp, PS_FLAG_ERR );
    
}
void sgw_3d_open_input(SMODEL *mod) {
    
}


void sgw_read_hot(SMODEL *mod) {
    char line[MAXLINE];           /* the input line */
    char name[MAXLINE];           /* the name of the current data set */
    char msg[MAXLINE];            /* for reporting to screen */
    char *data = NULL;
    int i, ie, itrns = 0, inode = 0, jnode = 0;     /* loop counter */
    double *tmp_data_set;         /* temporary space to store the current data set in */
    
    
    assert(mod);
    assert(mod->io);
    assert(mod->sgw);
    // these flag determine how to initialize arrays based on what is read in
    int hot_flag_con[mod->ntransport]; sarray_init_int(hot_flag_con, mod->ntransport);
    
    //aliases for convenience
    SIO io = *(mod->io);
    SGRID *grid = mod->grid;
    SGW *sgw = mod->sgw;
    int nnodes = grid->nnodes;
    int macro_nnodes = grid->macro_nnodes;
    int flag_found_phead = FALSE;
    int flag_found_thead = FALSE;
    
#ifdef _DEBUG
    printf("-- MYID %d reading the GW hotstart file: %s\n",grid->smpi->myid,io.hot.filename);
#endif
    
    /* allocates space to read the current data set */
    tmp_data_set = (double *) tl_alloc(sizeof(double), macro_nnodes);
    
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
        
        if (strncmp("IPH", name, 3) == AGREE) {
#ifdef _DEBUG
            root_print("------ reading groundwater pressure head.");
#endif
            flag_found_phead = TRUE;
            for (inode = 0; inode < nnodes; inode++)
                sgw->gw_phead[inode] = tmp_data_set[grid->node[inode].gid];
        }
        else if (strncmp("ITH", name, 3) == AGREE) {
#ifdef _DEBUG
            root_print("------ reading groundwater total freshwater head.");
#endif
            flag_found_thead = TRUE;
            for (inode = 0; inode < nnodes; inode++) {
                sgw->gw_phead[inode] = tmp_data_set[grid->node[inode].gid];
                sgw->gw_phead[inode]-= grid->node[inode].z;
            }
        }
        else if (strncmp("ICON", name, 4) == AGREE) {
            /* get the constituent number from the name line */
            itrns = read_int_field_custom(io, &data, NULL, "constituent ID", 2, TRUE) - 1;
            if (itrns < 0 || itrns >= mod->ntransport) {
                io_read_error(io, "Tried to read concentration for non existent transport " "constituent.", TRUE);
            }
            hot_flag_con[itrns] = ON;
#ifdef _DEBUG
            sprintf(msg, "------ reading transport concentrations (constituent " "ID: %d).", itrns + 1);
            root_print(msg);
#endif
            for (i = 0; i < nnodes; i++) {
                mod->con[itrns].concentration[i] = tmp_data_set[grid->node[i].gid];
            }
        }
    }
    io_save_line(&io, NULL, "", "");
    rewind(io.hot.fp);
    
    assert (flag_found_phead || flag_found_thead);
    for (itrns=0; itrns < mod->ntransport; itrns++) {
        assert(hot_flag_con[itrns] == ON);
    }
    
    /* force the initial condition to meet Dirichlet boundary conditions */
#ifdef _DEBUG
    root_print("------ (WARNING) Dirichlet boundary conditions are being applied to nodes (superseding initial conditions).");
#endif
    int istring = -1, isers = -1;
    for (i = 0; i < nnodes; i++) {
        istring = grid->node[i].string;
        if (istring > NORMAL) {
            /*mwg debug
             printf("sgw_read_hot node %d istring %d flow.bc_flag= %d\n",i,istring,mod->str_values[istring].flow.bc_flag);
             */
            if (mod->str_values[istring].flow.bc_flag == BCT_DIR) {
                //tl_error("\nhere1\n");
                isers = mod->str_values[istring].flow.iu_0;
                sgw->gw_phead[i] = sseries_get_value(isers, mod->series_head,0);
                sgw->gw_phead[i]-= grid->node[i].z;
                /*mwf debug
                 printf("sgw_read_hot BCT_DIR value= %g \n",sgw->gw_phead[i]);
                 */
                
            } else if (mod->str_values[istring].flow.bc_flag == BCT_PRS_DIR) {
                //tl_error("\nhere12\n");
                //printf("myid: %d global node[%d]: %d\n",grid->smpi->myid,i+1,grid->node[i].gid+1);
                //mwf this is not in the original AdH
                //assert(fabs(grid->node[i].z) < 1e-5 || fabs(grid->node[i].z + 1) < 1e-5);
                isers = mod->str_values[istring].flow.iu_0;
                sgw->gw_phead[i] = sseries_get_value(isers, mod->series_head,0);
                /*mwf debug
                 printf("sgw_read_hot BCT_PSI_DIR value= %g \n",sgw->gw_phead[i]);
                 */
            }
            
        }
    }
    
    istring = -1, isers = -1;
    for (i = 0; i < nnodes; i++) {
        istring = grid->node[i].string;
        if (istring > NORMAL) {
            for (itrns = 0; itrns < mod->ntransport; itrns++) {
                if (mod->str_values[istring].trans[itrns].bc_flag == BCT_DIR) {
                    mod->con[itrns].concentration[i] = sseries_get_value(isers, mod->series_head,0);
                }
            }
        }
    }
    
    
    /* assume for now that previous values are the same as hot-start values */
    sarray_copy_dbl(sgw->old_gw_phead, sgw->gw_phead, nnodes);
    sarray_copy_dbl(sgw->predict_gw_phead, sgw->gw_phead, nnodes);
    sarray_copy_dbl(sgw->older_gw_phead, sgw->old_gw_phead, nnodes);
    
    for (itrns=0; itrns<mod->ntransport; itrns++) {
        sarray_copy_dbl(mod->con[itrns].old_concentration, mod->con[itrns].concentration, nnodes);
        sarray_copy_dbl(mod->con[itrns].older_concentration, mod->con[itrns].old_concentration, nnodes);
    }
    
    /* default initialization for auxiliary variables */
    for (i = 0; i < nnodes; i++) {
        sgw->gw_density[i] = UNSET_FLT;
    }
    /* go head and try to calculate now
     TODO for SALINITY and HEAT set actual salinity or heat variables
     */
    tl_density_calculator_metric(mod->density, NULL, 1.0, NULL, 1.0, mod->grid->nnodes, mod->sgw->gw_density, 3);

    /* clean up memory */
    tmp_data_set = (double *) tl_free(sizeof(double), macro_nnodes, tmp_data_set);
    
}
/**
 calculate the saturations
 **/
void sgw_evaluate_element_saturations(SMODEL *mod) {
    SGRID *grid = mod->grid;
    SGW *sgw = mod->sgw;
    SMAT_GW * gwm;
    int ie,imat,II,i;
    double sat,psi,slope;
    SSERIES * psk;
    for (ie = 0; ie < grid->nelems3d; ie++) {
        imat = grid->elem3d[ie].mat;
        assert(imat >= 0); assert(imat < mod->nmat);
        gwm  = mod->mat[imat].gw;
        assert(gwm);
        assert(gwm->isat >= 0); assert(gwm->isat < 2*mod->nmat);
        psk = sseries_search(gwm->isat,NULL,mod->series_gw_psk_head);
        assert(psk);
        /*mwf debug
         printf("sgw_read_hot elem %d material %d\n",ie,imat);
         */
        for (i=0; i < grid->elem3d[ie].nnodes; i++) {
            II = grid->elem3d[ie].nodes[i];
            psi= sgw->gw_phead[II];
            /*mwf debug
             printf("sgw_read_hot local node %d --> %d val %g\n",i,II,psi);
             */
            sat = sgw_eval_sat(psk,psi,&sat,&slope);
            sgw->elem_3d_data[ie].saturation[i] = sat;
        }
    }
}

double sgw_eval_sat(SSERIES* psk_pc, double psi, double *sat, double *slope) {
    assert(psk_pc);
    double psi_eval=psi;
    if (psi_eval > psk_pc->entry[psk_pc->size-1].time) {
        *sat=1.;
        *slope=0.;
        return *sat;
    }
    if (psi < psk_pc->entry[0].time) {
        psi_eval=psk_pc->entry[1].time;
    }
    int interval = sseries_get_interval(*psk_pc,psi_eval);
    *sat = tc_eval_series(*psk_pc,interval,psi_eval,0);
    *slope=tc_eval_series_slope(*psk_pc,interval,psi_eval,0);
    return *sat;
}
double sgw_eval_kr(SSERIES* psk_kr, double psi, double *kr, double *slope) {
    assert(psk_kr);
    double psi_eval=psi;
    if (psi_eval > psk_kr->entry[psk_kr->size-1].time) {
        *kr=1.;
        *slope=0.;
        return *kr;
    }
    if (psi < psk_kr->entry[0].time) {
        psi_eval=psk_kr->entry[1].time;
    }
    
    int interval = sseries_get_interval(*psk_kr,psi_eval);
    *kr = tc_eval_series(*psk_kr,interval,psi_eval,0);
    *slope=tc_eval_series_slope(*psk_kr,interval,psi_eval,0);
    return *kr;
}
double sgw_eval_kr_elem(SSERIES* psk_kr, int nnodes, double *elem_gw_phead) {
    double alpha = 0.90;          /* weighting factor for relative permeability */
    double elem_kr = 1.0;         /* the elemental relative conductivity/permeability */
    double p_max = -1.0e6;        /* pressures for finding kr */
    double p_min = 1.0e6;         /* pressures for finding kr */
    double p_avg = DZERO;         /* pressures for finding kr */
    double k_max = DZERO;         /* Relative Conductivity */
    double k_min = DZERO;         /* Relative Conductivity */
    double k_avg = DZERO;         /* Relative Conductivity */
    double p_eval= DZERO;
    
    int ii = 0;                   /* loop counter */
    double slope = DZERO;
    
    /* use the arithmetic average of the max and min k_r
     needed to get influx into dry elements */
    for (ii = 0; ii < nnodes; ii++)
    {
        p_max = MAX(p_max, elem_gw_phead[ii]);
        p_min = MIN(p_min, elem_gw_phead[ii]);
    }
    p_avg = alpha * p_max + (1.0 - alpha) * p_min;
    k_min = sgw_eval_kr(psk_kr, p_min, &k_min, &slope);
    k_max = sgw_eval_kr(psk_kr, p_max, &k_max, &slope);
    elem_kr = alpha * k_max + (1.0 - alpha) * k_min;
    
    /** mwf try just straight average? /
     elem_kr = 0.; alpha = 1./nnodes;
     for (ii = 0; ii < nnodes; ii++)
     {
     p_eval = elem_gw_phead[ii];
     k_avg = sgw_eval_kr(psk_kr, p_eval, &k_avg, &slope);
     elem_kr += k_avg*alpha;
     }
     */
    assert (0. <= elem_kr);
    assert (elem_kr <= 1.0);
    return (elem_kr);
    
}

void sgw_evaluate_element_fluxes(SMODEL *mod){
    SELEM_3D *elem3d;
    SGW *gw = mod->sgw;
    SMAT_GW *gwm;
    SVECT *grad_shp;
    SSERIES *psk_kr;
    double ref_density = mod->density;
    int i, imat, nnodes;
    
    for (i=0; i<mod->grid->nelems3d; i++){
        elem3d = &(mod->grid->elem3d[i]);
        nnodes = elem3d->nnodes;
        double elem_phead[nnodes], elem_density[nnodes];
        imat = elem3d->mat;
        grad_shp = elem3d->grad_shp;
        gwm = mod->mat[imat].gw;
        
        psk_kr = sseries_search(gwm->ikr,NULL,mod->series_gw_psk_head);
        global_to_local_dbl(gw->gw_density,elem_density, elem3d->nodes,nnodes);
        global_to_local_dbl(gw->gw_phead, elem_phead, elem3d->nodes, nnodes);
        
        sgw_eval_flux(nnodes, psk_kr, gwm->k, ref_density, elem_density, elem_phead, grad_shp, &(gw->elem_gw_flux[i]));
    }
}

void sgw_eval_flux(
                   int nnodes,
                   SSERIES *psk_kr,
                   const STENSOR k, /* DO NOT modify this variable! */
                   double ref_density,
                   double* elem_density,
                   double *elem_phead,
                   SVECT *grad_shp,
                   SVECT *darcy_flux
                   ) {
    
    int i;
    SVECT grad_phead, gradient_z;
    SVECT k_grad_phead,k_grad_z;
    STENSOR k_total;
    double avg_density=0.0;
    double elem_kr;
    
    /* Calculate average density */
    for (i=0;i<nnodes;i++){
        avg_density+=elem_density[i];
    }
    avg_density/=(double)nnodes;
    
    VT_3D_TCOPY(k, k_total);
    
    assert(psk_kr);
    elem_kr = sgw_eval_kr_elem(psk_kr,nnodes,elem_phead);
    
    /* scales by the relative permeability */
    VT_3D_TSCALE(k_total, elem_kr);
    
    /* Set grad_z = {0,0,1} */
    gradient_z.x=0.; gradient_z.y=0.; gradient_z.z=1.;
    
    /* Initialize grad_phead */
    grad_phead.x=0.; grad_phead.y=0.; grad_phead.z=0.;
    for (i=0; i < nnodes; i++) {
        grad_phead.x += elem_phead[i]*grad_shp[i].x;
        grad_phead.y += elem_phead[i]*grad_shp[i].y;
        grad_phead.z += elem_phead[i]*grad_shp[i].z;
    }
    
    VT_3D_TENS_VECT_PROD(k_grad_phead,k_total,grad_phead);
    VT_3D_TENS_VECT_PROD(k_grad_z,k_total,gradient_z);
    
    /*mwf not allowing density dependence for now */
    assert(fabs(avg_density/ref_density-1.0) <= SMALL);
    darcy_flux->x = -(k_grad_phead.x + avg_density/ref_density*k_grad_z.x);
    darcy_flux->y = -(k_grad_phead.y + avg_density/ref_density*k_grad_z.y);
    darcy_flux->z = -(k_grad_phead.z + avg_density/ref_density*k_grad_z.z);
}

void sgw_print_ts(SMODEL *mod, SGW *gw, SIO *info, SGRID *grid, double time, int outfact, SFLAGS flag, SFILE sup, char *proj_name, int it1, int it2, int **ndata, int my_nnode_max, int *my_nnode_ext, int **ndata_sur, int my_nnode_max_sur, int *my_nnode_ext_sur, int flag1) {
    int super = strlen(sup.filename);
    int ie = 0, inode = 0;
    
    SFILE fout;
    int my_nnode = grid->my_nnodes;
    my_nnode = grid->nnodes; // cjt
    
    /* allocate temporary space for nodally-averaged variables to be printed */
    double * gw_print_head = (double *) tl_alloc(sizeof(double), my_nnode);
    double * nodal_sat = (double *) tl_alloc(sizeof(double), my_nnode);
    int * nodal_cnt = (int *) tl_alloc(sizeof(int), my_nnode);
    SVECT * nodal_gw_vel = (SVECT *) tl_alloc(sizeof(SVECT), my_nnode);
    
    if (grid->smpi->myid <= 0) {
        /*mwf debug
         printf("In sgw_print_ts time=%g \n",time);
         */
        
        fprintf(info->fout_gw_phead.fp, "TS 0 %15.8e\n", time);
        fprintf(info->fout_gw_thead.fp, "TS 0 %15.8e\n", time);
        fprintf(info->fout_gw_density.fp, "TS 0 %15.8e\n", time);
        fprintf(info->fout_gw_sat.fp, "TS 0 %15.8e\n", time);
        fprintf(info->fout_gw_flx.fp, "TS 0 %15.8e\n", time);
        fprintf(info->fout_gw_vel.fp, "TS 0 %15.8e\n", time);
        fprintf(info->fout_gw_error.fp, "TS 0 %15.8e\n", time);
        
        init_adh_file(&(fout));
        if (flag1) {
            build_filename2(fout.filename, MAXLINE, proj_name, "_phd.dat-", it1,".", it2); // build the new filename
            open_output_file(&(fout), "adapted pressure head file", super); // open the file
            fprintf(fout.fp, "DATASET\n");
            pdata(proj_name, "", fout.fp, 0, "Pressure Head", "", "BEGSCL", "mesh3d", grid->macro_nnodes, grid->macro_nelems3d); // header
            tc_timeunits(fout.fp, outfact);
            fprintf(fout.fp, "TS 0 %15.8e\n", time);    // print the time0-header
        }
        
    }/* end myid <= 0*/
    sarray_print_dbl(grid, info->fout_gw_phead.fp, fout.fp, gw->gw_phead, ndata, my_nnode_max, my_nnode_ext, flag1);    // print the variables
    
    if ((grid->smpi->myid <= 0) && flag1) print_trailer(fout.fp); // add ENNDS, close the file and NULL the pointer
    
    
    if ((grid->smpi->myid <= 0) && flag1) {
        init_adh_file(&(fout));     // initialize file
        build_filename2(fout.filename, MAXLINE, proj_name, "_thd.dat-", it1,".", it2);
        open_output_file(&(fout), "adapted total head file", super);
        fprintf(fout.fp, "DATASET\n");
        pdata(proj_name, "", fout.fp, 0, "Error", "", "BEGSCL", "mesh3d", grid->macro_nnodes, grid->macro_nelems3d);
        tc_timeunits(fout.fp, outfact);
        fprintf(fout.fp, "TS 0 %15.8e\n", time);
    }
    for (inode = 0; inode < my_nnode; inode++) {
        /* need correct version for density dependent flow as well */
        gw_print_head[inode] = gw->gw_phead[inode] + grid->node[inode].z;
    }
    sarray_print_dbl(grid, info->fout_gw_thead.fp, fout.fp, gw_print_head, ndata, my_nnode_max, my_nnode_ext, flag1);    // print the variables
    if ((grid->smpi->myid <= 0) && flag1) print_trailer(fout.fp);

    
    /** TODO:
     add printing of groundwater flux
     
     **/
    if ((grid->smpi->myid <= 0) && flag1) {
        init_adh_file(&(fout));     // initialize file
        build_filename2(fout.filename, MAXLINE, proj_name, "_den.dat-", it1,".", it2);
        open_output_file(&(fout), "adapted density file", super);
        fprintf(fout.fp, "DATASET\n");
        pdata(proj_name, "", fout.fp, 0, "Density", "", "BEGSCL", "mesh3d", grid->macro_nnodes, grid->macro_nelems3d);
        tc_timeunits(fout.fp, outfact);
        fprintf(fout.fp, "TS 0 %15.8e\n", time);
    }
    /** density -- assuming already calculated **/
    sarray_print_dbl(grid, info->fout_gw_density.fp, fout.fp, gw->gw_density, ndata, my_nnode_max, my_nnode_ext, flag1);    // print the variables
    if ((grid->smpi->myid <= 0) && flag1) print_trailer(fout.fp);
    
#ifdef _DEBUG
    int DEBUG=OFF;
    if (DEBUG) tl_check_all_pickets(__FILE__,__LINE__);
#endif

    if ((grid->smpi->myid <= 0) && flag1) {
        init_adh_file(&(fout));
        build_filename2(fout.filename, MAXLINE, proj_name, "_sat.dat-", it1,".", it2);
        open_output_file(&(fout), "adapted saturation file", super);
        fprintf(fout.fp, "DATASET\n");
        pdata(proj_name, "", fout.fp, 0, "Saturation", "", "BEGSCL", "mesh3d", grid->macro_nnodes, grid->macro_nelems3d);
        tc_timeunits(fout.fp, outfact);
        fprintf(fout.fp, "TS 0 %15.8e\n", time);
    }
    
    /** average saturations at nodes **/
    sgw_nodal_sat_avg(gw,grid,my_nnode, nodal_sat, nodal_cnt);
    sarray_print_dbl(grid, info->fout_gw_sat.fp, fout.fp, nodal_sat, ndata, my_nnode_max, my_nnode_ext, flag1);    // print the variables
    if ((grid->smpi->myid <= 0) && flag1) print_trailer(fout.fp);

    
    /** error -- assuming already calculated **/
    if ((grid->smpi->myid <= 0) && flag1) {
        init_adh_file(&(fout));
        build_filename2(fout.filename, MAXLINE, proj_name, "_error.dat-", it1,".", it2);
        open_output_file(&(fout), "adapted error file", super);
        fprintf(fout.fp, "DATASET\n");
        pdata(proj_name, "", fout.fp, 0, "Error", "", "BEGSCL", "mesh3d", grid->macro_nnodes, grid->macro_nelems3d);
        tc_timeunits(fout.fp, outfact);
        fprintf(fout.fp, "TS 0 %15.8e\n", time);
    }
    sarray_print_dbl(grid, info->fout_gw_error.fp, fout.fp, gw->error, ndata, my_nnode_max, my_nnode_ext, flag1);    // print the variables
    if ((grid->smpi->myid <= 0) && flag1) print_trailer(fout.fp);
    
    /* Print velocities */
    if ((grid->smpi->myid <= 0) && flag1) {
        init_adh_file(&(fout));
        build_filename2(fout.filename, MAXLINE, proj_name, "_vel.dat-", it1,".", it2);
        open_output_file(&(fout), "adapted velocity file", super);
        fprintf(fout.fp, "DATASET\n");
        pdata(proj_name, "", fout.fp, 0, "Velocity", "", "BEGVEC", "mesh3d", grid->macro_nnodes, grid->macro_nelems3d);
        tc_timeunits(fout.fp, outfact);
        fprintf(fout.fp, "TS 0 %15.8e\n", time);
    }
    
    /** average elemental darcy velocity to nodes **/
    SVECT * nodal_vel = (SVECT *) tl_alloc(sizeof(SVECT), grid->nnodes);
    sgw_elem3d_flux_to_nodal_vel(grid, mod->mat, gw, nodal_vel);
    svect_print_array_MPI(grid, info->fout_gw_vel.fp, fout.fp, nodal_vel, ndata, my_nnode_max, my_nnode_ext, 0);
    if ((grid->smpi->myid <= 0) && flag1) print_trailer(fout.fp);
    nodal_vel = (SVECT *) tl_free(sizeof(SVECT), grid->nnodes, nodal_vel);
    
    /* deallocate temporary space for nodally-averaged variables to be printed */
    gw_print_head = (double *) tl_free(sizeof(double), my_nnode, gw_print_head);
    nodal_sat = (double *) tl_free(sizeof(double), my_nnode, nodal_sat);
    nodal_cnt = (int *) tl_free(sizeof(int), my_nnode, nodal_cnt);
    nodal_gw_vel = (SVECT *) tl_free(sizeof(SVECT), my_nnode, nodal_gw_vel);
    
}

void sgw_nodal_sat_avg(SGW* gw, SGRID* grid, int max_num_nodes, double *sat_avg, int* cnt_arr) {
    int ii,ie,II;
    double wgt;
    assert(gw);
    assert(grid);
    assert(sat_avg);
    assert(cnt_arr);
    
    for (ii=0; ii < max_num_nodes; ii++) {
        sat_avg[ii] = 0.0;
        grid->node[ii].area = 0.;
        cnt_arr[ii] = 0.;
    }
        //tl_check_all_pickets(__FILE__,__LINE__);
    
    for (ie=0; ie < grid->nelems3d; ie++) {
        for (ii=0; ii < grid->elem3d[ie].nnodes; ii++) {
            II = grid->elem3d[ie].nodes[ii];
            //printf("II: %d my_nnodes: %d \n",II,max_num_nodes);
            wgt = 0.25*grid->elem3d[ie].djac;
            sat_avg[II] += gw->elem_3d_data[ie].saturation[ii]*wgt;
            grid->node[II].area += wgt;
            cnt_arr[II] += 1;
        }
    }
    
        //tl_check_all_pickets(__FILE__,__LINE__);
    
    for (ii=0; ii < max_num_nodes; ii++) {
        if (cnt_arr[ii] > 0) {
            assert(grid->node[ii].area > 0.);
            sat_avg[ii] /= grid->node[ii].area;
        }
    }
    /*mwf debug
     for (ii=0; ii < max_num_nodes; ii++) {
     printf("sgw_nodal_sat_avg node ii= %d cnt=%d area= %g\n",ii,cnt_arr[ii],grid->node[ii].area);
     }
     */
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Calculate nodally averaged velocities from element-centered darcy flux and saturations.
 * \author    Gajanan Choudhary, PhD.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in]     grid    (SGRID *) a pointer to the grid
 * @param[in]     mat     (SMAT *) a pointer to the materials
 * @param[in]     gw      (SGW *) a pointer to the groundwater solution variables
 * @param[in,out] vel_avg (SVECT *) a pointer to node-averaged velocity, assumed allocated
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sgw_elem3d_flux_to_nodal_vel(SGRID* grid, SMAT *mat, SGW* gw, SVECT *vel_avg) {
    int ii,ie,II, imat;
    double wgt, porosity;
    
    /* Allocate temporary variable */
    int *nodal_cnt = (int *) tl_alloc(sizeof(int), grid->nnodes);
    double *nodal_wgt = (double *) tl_alloc(sizeof(double), grid->nnodes);
    
    assert(mat);
    assert(gw);
    assert(grid);
    assert(vel_avg);
    assert(nodal_cnt);
    
    for (ii=0; ii < grid->nnodes; ii++) {
        vel_avg[ii].x = 0.0;
        vel_avg[ii].y = 0.0;
        vel_avg[ii].z = 0.0;
        nodal_wgt[ii] = 0.;
        nodal_cnt[ii] = 0.;
    }
    
    for (ie=0; ie < grid->nelems3d; ie++) {
        wgt = 0.25*grid->elem3d[ie].djac;
        //printf("\nelem[%3i].djac=%f",ie,4.0*wgt);
        imat = grid->elem3d[ie].mat;
        porosity = mat[imat].gw->porosity;
        for (ii=0; ii < grid->elem3d[ie].nnodes; ii++) {
            II = grid->elem3d[ie].nodes[ii];
            nodal_wgt[II] += wgt;
            vel_avg[II].x += wgt * gw->elem_gw_flux[ie].x/(porosity*gw->elem_3d_data[ie].saturation[ii] + SMALL);
            vel_avg[II].y += wgt * gw->elem_gw_flux[ie].y/(porosity*gw->elem_3d_data[ie].saturation[ii] + SMALL);
            vel_avg[II].z += wgt * gw->elem_gw_flux[ie].z/(porosity*gw->elem_3d_data[ie].saturation[ii] + SMALL);
            nodal_cnt[II] += 1;
        }
    }
#ifdef _MESSG
    comm_update_double((double*) vel_avg, 3, grid->smpi);
    comm_update_double(nodal_wgt, 1, grid->smpi);
    comm_update_int(nodal_cnt, 1, grid->smpi);
#endif
    for (ii=0; ii < grid->nnodes; ii++) {
        if (nodal_cnt[ii] > 0) {
            assert(nodal_wgt[ii] > 0.);
            vel_avg[ii].x /= nodal_wgt[ii];
            vel_avg[ii].y /= nodal_wgt[ii];
            vel_avg[ii].z /= nodal_wgt[ii];
        }
    }
    /* Free the temporary variable */
    nodal_cnt = (int *)    tl_free(sizeof(int), grid->nnodes, nodal_cnt);
    nodal_wgt = (double *) tl_free(sizeof(double), grid->nnodes, nodal_wgt);
}


void sgw_set_psk_series(SMODEL* mod) {
    /* loop through the materials and create relative perm and capillary pressure for each one*/
    int imat,isers,size,nnodes,ientry;
    int npsk=0;
    SSERIES *series_kr,*series_pc;
    SMAT_GW* gwm;
    int generate_psk= 0;
    
    assert(mod);
    assert(mod->flag.GW_FLOW);
    for (imat=0; imat < mod->nmat; imat++) {
        gwm = mod->mat[imat].gw;
        assert (gwm);
        double alpha = gwm->vangen_alpha;
        double n_exponent = gwm->vangen_n;
        double m_exponent = 1. - 1./n_exponent;
        double lambda= gwm->brooks_lambda;
        double pd    = gwm->brooks_pd;
        double residual_saturation = gwm->residual_sat;
        
        generate_psk= (gwm->isat == UNSET_INT) || (gwm->ikr == UNSET_INT) ;
        
        if (gwm->isat != UNSET_INT) {
            assert (gwm->ikr != UNSET_INT);
        }
        if (gwm->ikr != UNSET_INT) {
            assert (gwm->isat != UNSET_INT);
        }
        if (!generate_psk) {
            break;
        }
        if (alpha > SMALL || lambda > SMALL) {
            nnodes=1;
            if (alpha > SMALL) {
                size = gwm->vangen_num_xy;
            } else {
                size = gwm->brooks_num_xy;
            }
            /*capillary pressure */
            
            sseries_alloc(&series_pc,size,CONSTITUITIVE_SERIES,nnodes);
            series_pc->id=npsk++;
            series_pc->nnodes=nnodes;
            gwm->isat = series_pc->id;
            
            /* relative perm */
            sseries_alloc(&series_kr,size,CONSTITUITIVE_SERIES,nnodes);
            series_kr->id=npsk++;
            series_kr->nnodes=nnodes;
            gwm->ikr = series_kr->id;
            /*mwf debug */
            printf("Building psk series\n");
            
            if (alpha > SMALL) {
                for (ientry=0; ientry < series_kr->size; ientry++) {
                    double cap_press_head = gwm->vangen_max_cp * (1. - log10(ientry + 1) / log10(gwm->vangen_num_xy));
                    double eff_sat = pow(1.0 + pow(alpha * cap_press_head, n_exponent), -m_exponent);
                    double kr = sqrt(eff_sat) * pow((1.0 - pow(1.0 - pow(eff_sat, 1.0 / m_exponent), m_exponent)), 2);
                    series_kr->entry[ientry].time = -cap_press_head;
                    series_kr->entry[ientry].value[0]=kr;
                    
                    series_pc->entry[ientry].time = -cap_press_head;
                    series_pc->entry[ientry].value[0]=residual_saturation + (1.0 - residual_saturation) * eff_sat;
                }
            } else {
                printf("Brooks Corey psk series not implemented yet\n");
                assert(0);
            }
            /*compute the slopes */
            for (ientry=1; ientry < series_kr->size; ientry++) {
                series_kr->entry[ientry-1].slope[0] = (series_kr->entry[ientry].value - series_kr->entry[ientry-1].value) / (series_kr->entry[ientry].time - series_kr->entry[ientry-1].time);
            }
            for (ientry=1; ientry < series_pc->size; ientry++) {
                series_pc->entry[ientry-1].slope[0] = (series_pc->entry[ientry].value - series_pc->entry[ientry-1].value) / (series_pc->entry[ientry].time - series_pc->entry[ientry-1].time);
            }
            /*compute the areas */
            for (ientry=0; ientry < series_kr->size-1; ientry++) {
                series_kr->entry[ientry].area[0] = tc_trap_area(series_kr->entry[ientry].value[0],series_kr->entry[ientry+1].value[0],
                                                                series_kr->entry[ientry].time, series_kr->entry[ientry+1].time);
            }
            for (ientry=0; ientry < series_pc->size-1; ientry++) {
                series_pc->entry[ientry].area[0] = tc_trap_area(series_pc->entry[ientry].value[0],series_pc->entry[ientry+1].value[0],
                                                                series_pc->entry[ientry].time, series_pc->entry[ientry+1].time);
            }
            
            /*mwf debug
             printf("DEBUG sgw_set_psk_series mat %d pc \n",imat);
             sseries_printScreen(*series_pc,1);
             */
            /*mwf end debug*/
            sseries_add(series_pc,&(mod->series_gw_psk_head),&(mod->series_gw_psk_curr), TRUE);
            sseries_free(series_pc);
            /*
             printf("DEBUG sgw_set_psk_series mat %d pc curr \n",imat);
             sseries_printScreen(*mod->series_gw_psk_curr,1);
             */
            /*
             printf("DEBUG sgw_set_psk_series mat %d kr \n",imat);
             sseries_printScreen(*series_kr,1);
             */
            sseries_add(series_kr,&(mod->series_gw_psk_head),&(mod->series_gw_psk_curr), TRUE);
            sseries_free(series_kr);
            /*
             printf("DEBUG sgw_set_psk_series mat %d kr curr \n",imat);
             sseries_printScreen(*mod->series_gw_psk_curr,1);
             */
            
        }
    }
    
}

/***********************************************************/
/* Functions related to mesh adaption. */
/***********************************************************/

/***********************************************************/
/***********************************************************/
/***********************************************************/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Calculate elemental errors for mesh (un)refinement.
 * \author    Gajanan Choudhary, PhD.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] gw      (SGW *) a pointer to the groundwater solution variables
 * @param[in]     grid    (SGRID *) a pointer to the grid
 * @param[in]     mat     (SMAT *) a pointer to the materials
 * @param[in]     dt      (double *) time step size
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* Approximates elemental error for the 3d richards equations */
void sgw_3d_calculate_elem_error(SGW *gw, SGRID *grid, SMAT *mat, double dt) {
    
    /* calculates the error indicator for the flow */
    int ie = UNSET_INT, i=UNSET_INT;;
    SVECT elem_vel[MAX_NNODES_ON_ELEM3D]; svect_init_array(elem_vel, MAX_NNODES_ON_ELEM3D);
    double elem_head[MAX_NNODES_ON_ELEM3D]; sarray_init_dbl(elem_head, MAX_NNODES_ON_ELEM3D);
    double elem_z[MAX_NNODES_ON_ELEM3D]; sarray_init_dbl(elem_z, MAX_NNODES_ON_ELEM3D);
    double error = 0., continuity = 0.;
    double dhdx, dhdy, dhdz;
    
    sarray_init_dbl(gw->darray, grid->nnodes); // element per node count for averaging
    sarray_init_dbl(gw->error, grid->nnodes);  // nodal error for printing
    
    /* calculate the velocities for the flow */
    SVECT * nodal_vel = (SVECT *) tl_alloc(sizeof(SVECT), grid->nnodes);
    sgw_elem3d_flux_to_nodal_vel(grid, mat, gw, nodal_vel);
    
    int nelems3d = grid->nelems3d; // alias
    for (ie = 0; ie < nelems3d; ie++) {
        error = 0.;
        grid->elem_error[ie] = 0.;
        
        ELEM3D_GET_LOCAL(gw->gw_phead, elem_head, grid->elem3d[ie].nodes);
        ELEM3D_GET_LOCAL_VECT(nodal_vel, elem_vel, grid->elem3d[ie].nodes);
        
        /* Gajanan gkc adding after finding out elem->grad_shp is NaN in current version. */
        int nnodes = grid->elem3d[ie].nnodes;
        //SNODE elem_nodes[NDONPRISM];
        for (i=0; i<grid->elem3d[ie].nnodes; i++) {
            //snode_copy(&(elem_nodes[i]), grid->node[grid->elem3d[ie].nodes[i]]);
            elem_z[i] = grid->node[grid->elem3d[ie].nodes[i]].z;
        }
        //printf("\nElem %7i:", ie);
        if (nnodes == NDONTET){
            // continuity = 0.;
            // for (i=0; i<grid->elem3d[ie].nnodes; i++) {
            //     continuity += svect_dotp(grid->elem3d[ie].grad_shp[i], elem_vel[i]);
            // }
            // error = sqrt(continuity * continuity) * grid->elem3d[ie].djac;
            dhdx = 0.; dhdy = 0.; dhdz = 0.;
            for (i=0; i<grid->elem3d[ie].nnodes; i++) {
                dhdx += grid->elem3d[ie].grad_shp[i].x * (elem_head[i] + elem_z[i]);
                dhdy += grid->elem3d[ie].grad_shp[i].y * (elem_head[i] + elem_z[i]);
                dhdz += grid->elem3d[ie].grad_shp[i].z * (elem_head[i] + elem_z[i]);
            }
            error = grid->elem3d[ie].djac * fabs(dhdx+dhdy+dhdz) * sqrt(1/24.);
            
        }
        else if (nnodes == NDONPRISM){
            /* Gajanan gkc: prisms not supported for now. */
        }
        
        for (i=0; i<grid->elem3d[ie].nnodes; i++) {
            gw->error[grid->elem3d[ie].nodes[i]] += error * grid->elem3d[ie].djac;
            gw->darray[grid->elem3d[ie].nodes[i]] +=  grid->elem3d[ie].djac;
        }
        
        /* scales by the tolerance */
        error /= mat[grid->elem3d[ie].mat].gw->refine_tolerance;
        
        // sets as max error, since only one gw component can be run for now (con & sed can change)
        grid->elem_error[ie] = error;
    }
    
    // average the elemental to nodal errors for output
    for (i=0; i<grid->nnodes; i++) {
        gw->error[i] /= (double)gw->darray[i];
    }
    sarray_init_dbl(gw->darray, grid->nnodes); // reset initialize utility array
    nodal_vel = (SVECT *) tl_free(sizeof(SVECT), grid->nnodes, nodal_vel);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// averages two nodes for a  new node
void sgw_3d_renumber(SGW *gw, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp, SVECT2D *v2tmp, SVECT *vtmp) {
    // note darray and iarray are just utility arrays, no need to worry about them here
    
    node_renumber_double(max_nnode, gw->gw_phead, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, gw->old_gw_phead, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, gw->older_gw_phead, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, gw->predict_gw_phead, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, gw->gw_density, dtmp, new_numbers, order_tmp);
    node_renumber_double(max_nnode, gw->error, dtmp, new_numbers, order_tmp);
    //node_renumber_vect(max_nnode, gw->vel, vtmp, new_numbers, order_tmp);
    //node_renumber_vect(max_nnode, gw->old_vel, vtmp, new_numbers, order_tmp);
    //node_renumber_vect(max_nnode, gw->older_vel, vtmp, new_numbers, order_tmp);
    
    /* Gajanan gkc: Must also deal with elem_gw_flux and elem_3d_data.
     * This isn't being done just yet*/
}

/************************************************************/
/************************************************************/
/************************************************************/
// reallocs 3d variables when a new node increment is added
void sgw_3d_realloc_init(SGW *gw, int nnodes_old, int nnodes_new) {
    int i,j;
    
    gw->iarray = (int *) tl_realloc(sizeof(int), nnodes_new, nnodes_old, gw->iarray);
    gw->darray = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, gw->darray);
    gw->gw_phead = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, gw->gw_phead);
    gw->old_gw_phead = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, gw->old_gw_phead);
    gw->older_gw_phead = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, gw->older_gw_phead);
    gw->predict_gw_phead = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, gw->predict_gw_phead);
    gw->gw_density = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, gw->gw_density);
    gw->error = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, gw->error);
    //gw->vel = (SVECT *) tl_realloc(sizeof(SVECT), nnodes_new, nnodes_old,gw->vel);
    //gw->old_vel = (SVECT *) tl_realloc(sizeof(SVECT), nnodes_new, nnodes_old, gw->old_vel);
    //gw->older_vel = (SVECT *) tl_realloc(sizeof(SVECT), nnodes_new, nnodes_old, gw->older_vel);
    
    for (i=nnodes_old;i<nnodes_new;i++){
        gw->iarray[i] = UNSET_INT;
        gw->darray[i] = 0.;
        gw->gw_phead[i] = 0.;
        gw->gw_phead[i] = 0.;
        gw->old_gw_phead[i] = 0.;
        gw->older_gw_phead[i] = 0.;
        gw->predict_gw_phead[i] = 0.;
        gw->gw_density[i] = 0.;
        gw->error[i] = 0.;
    }
}

/************************************************************/
/************************************************************/
/************************************************************/
// reallocs 3d element variables when a new node increment is added
void sgw_3d_elem_realloc_init(SGW *gw, int nelems3d_new, int nelems3d_old, int nalloc_inc) {
    int i,j;
    
    if(nelems3d_new == nelems3d_old) {
        nelems3d_new += nalloc_inc;
        gw->elem_3d_data = (ELEMENT_3D_DATA *) tl_realloc(sizeof(ELEMENT_3D_DATA), nelems3d_new, nelems3d_old, gw->elem_3d_data);
        gw->elem_gw_flux = (SVECT *) tl_realloc(sizeof(SVECT), nelems3d_new, nelems3d_old, gw->elem_gw_flux);
        
        for (i=nelems3d_old;i<nelems3d_new;i++){
            for(j=0; j<MAX_NNODES_ON_ELEM3D; j++){
                gw->elem_3d_data[i].saturation[j] = 0.;
                gw->elem_3d_data[i].old_saturation[j] = 0.;
            }
            gw->elem_gw_flux[i].x = 0.;
            gw->elem_gw_flux[i].y = 0.;
            gw->elem_gw_flux[i].z = 0.;
        }
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// averages two nodes for a new node
void sgw_3d_node_avg(SGW *gw, int node_new, int node1, int node2, SGRID *grid) {
    
    // 3D :: (cjt) note: darray and iarray are just utility arrays, no need to worry about them here
    gw->gw_phead[node_new] = 0.5*(gw->gw_phead[node1] + gw->gw_phead[node2]);
    gw->old_gw_phead[node_new] = 0.5*(gw->old_gw_phead[node1] + gw->old_gw_phead[node2]);
    gw->older_gw_phead[node_new] = 0.5*(gw->older_gw_phead[node1] + gw->older_gw_phead[node2]);
    gw->predict_gw_phead[node_new] = 0.5*(gw->predict_gw_phead[node1] + gw->predict_gw_phead[node2]);
    gw->gw_density[node_new] = 0.5*(gw->gw_density[node1] + gw->gw_density[node2]);
    gw->error[node_new] = 0.5*(gw->error[node1] + gw->error[node2]);
    //gw->vel[node_new] = svect_avg(gw->vel[node1], gw->vel[node2]);
    //gw->old_vel[node_new] = svect_avg(gw->old_vel[node1], gw->old_vel[node2]);
    //gw->older_vel[node_new] = svect_avg(gw->older_vel[node1], gw->older_vel[node2]);
    
    /* Gajanan gkc: Must also deal with elem_gw_flux and elem_3d_data.
     * This isn't being done just yet*/
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// averages two elems for a new elem
void sgw_3d_elem_avg(SGW *gw, int orig_elem, int new_elem, int new_node_in_orig_elem, int new_node_in_new_elem, SGRID *grid) {
    
    int i;
    int nd1 = new_node_in_orig_elem;
    int nd2 = new_node_in_new_elem;
    
    // 3D :: (cjt) note: darray and iarray are just utility arrays, no need to worry about them here
    gw->elem_gw_flux[new_elem].x = gw->elem_gw_flux[orig_elem].x;
    gw->elem_gw_flux[new_elem].y = gw->elem_gw_flux[orig_elem].y;
    gw->elem_gw_flux[new_elem].z = gw->elem_gw_flux[orig_elem].z;
    for (i=0; i<MAX_NNODES_ON_ELEM3D; i++){
        gw->elem_3d_data[new_elem].saturation[i] = gw->elem_3d_data[orig_elem].saturation[i];
        gw->elem_3d_data[new_elem].old_saturation[i] = gw->elem_3d_data[orig_elem].old_saturation[i];
    }
    
    gw->elem_3d_data[orig_elem].saturation[new_node_in_orig_elem] =
    0.5 * ( gw->elem_3d_data[orig_elem].saturation[nd1] +
           gw->elem_3d_data[orig_elem].saturation[nd2] );
    gw->elem_3d_data[new_elem].old_saturation[new_node_in_new_elem] =
    0.5 * ( gw->elem_3d_data[new_elem].old_saturation[nd1] +
           gw->elem_3d_data[new_elem].old_saturation[nd2] );
    
    /* Gajanan gkc: Must also deal with elem_gw_flux and elem_3d_data.
     * This isn't being done just yet*/
}
/***********************************************************/
/***********************************************************/
/***********************************************************/
#endif

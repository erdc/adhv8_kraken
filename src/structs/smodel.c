/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates, initializes and reads in finite element geometry files.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mod           (SMODEL *)  a pointer to an AdH model struture
 * @param[in]  adh_root   (char *) a pointer to a root AdH filename
 *
 * \note CJT \:: Order is VERY important here!
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

// file scope variables
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smodel_init(SMODEL *mod, const char *adh_root) {
    
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) {
        DEBUG = ON;
    }
#endif
    
    int myid = 0, ierr = UNSET_INT;
    
    /* model base name */
    strcpy(mod->name,adh_root); //printf("mod->name: %s\n",mod->name);
    
    /* initialize debug options */
    sdebug_init(&(debug));
    
    /* initialize simulation screen output options */
    sscreen_output_init();
    
    /* initialize optional output files */
    sfile_output_init(&(mod->file_output));
    
    /* initialize test case parameters */
    //stestcase_init();  cjt :: move this to adh_init
    
    /* set model default parameters */
    int old_flag=0;
    if(mod->proc_flag==1)old_flag=mod->proc_flag;
    smodel_defaults(mod);
    mod->proc_flag=old_flag;
#ifdef _MESSG
    ierr = MPI_Comm_rank(cstorm_comm, &myid);
    printf("\n-- CSTORM COMM PE: %d :: ADH model initializing for project: %s\n",myid,adh_root);
#else
    printf("\n-- ADH model initializing for project: %s\n", adh_root);
#endif
    /* sets up the linked lists */
    tl_list_setup();
    
    /* initialize i/o */
    mod->io = (SIO *) tl_alloc(sizeof(SIO), 1);
    sio_init(mod->io, adh_root);
    if(mod->myid == 0){
        print_build_info();
        print_runtime_info(*(mod->io));
    }
    
    /* opens the basic input files */
    smodel_open_input(mod);
    
    /* read physics */
    read_physics(mod);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /* count the number of series and strings */
    read_bc_prep(mod);
    #ifdef _DEBUG
        if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    #endif
    
    //*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* allocate and initialize grids */
#ifdef _MESSG
    sgrid_alloc_init(&(mod->grid), mod->io, mod->proc_flag, mod->file_output, mod->model_comm);
#else
    sgrid_alloc_init(&(mod->grid), mod->io, mod->proc_flag, mod->file_output);
#endif
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    //*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    
#ifdef _MESSG
    if(mod->grid->smpi->npes>1) mod->grid->nmat = messg_imax(mod->grid->nmat, mod->grid->smpi->ADH_COMM);
#endif
    mod->nmat = mod->grid->nmat;

    /* initialize transport model if used */
    if (mod->ntransport > 0) {
        int nelems = 0;
        if (mod->grid->ndim == 2) {
            nelems = mod->grid->nelems2d;
        } else if (mod->grid->ndim == 3) {
            nelems = mod->grid->nelems3d;
        } else {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> grid dimension not supported in AdH.");
        }
        scon_alloc_init(mod->grid->nnodes, nelems, mod->ntransport, &(mod->con));
    }
    
    /* initialize materials */
    smat_init(mod, &(mod->mat), mod->nmat);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /* initialize strings */
    sstr_value_init(&(mod->str_values), mod->nstring, mod->ntransport, mod->nsed);
#ifdef _ADH_STRUCTURES
    if(mod->myid==0) printf("Number of weirs %d Number of Flap Gates %d \n",mod->nweir, mod->nflap);
    if (mod->nweir > 0)
        sstructures_init_weirs(&(mod->weir), mod->nstring, mod->nweir);
    if (mod->nflap > 0)
        sstructures_init_flaps(&(mod->flap), mod->nstring, mod->nflap);
    if (mod->nsluice > 0)
        sstructures_init_sluices(&(mod->sluice), mod->nstring, mod->nsluice);
#endif
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /* initialize solver to defaults */
    solv_initialize(&(mod->solver_info));
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /* initialize maximum number of equations */
    if (mod->flag.SW_FLOW == ON && mod->flag.DIFFUSIVE_WAVE == OFF) {
        mod->max_nsys = 3;
        mod->max_nsys_sq = 9;
    } else if (mod->flag.NS3_FLOW == ON) {
        mod->max_nsys = 4;
        mod->max_nsys_sq = 16;
    } else {
        mod->max_nsys = 1;
        mod->max_nsys_sq = 1;
    }

    //*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    read_bc(mod);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    /* initialize model physics */
    if (mod->flag.SW_FLOW && !mod->flag.GW_FLOW) {
        printf("\n\nCalling ssw_alloc_init");
        ssw_alloc_init(&(mod->sw), mod->grid, mod->io, mod->flag);
    } else if (mod->flag.NS_FLOW){
        sns_alloc_init(&(mod->ns), mod->grid, mod->io, mod->flag);
    }
#ifdef _ADH_GROUNDWATER
    else if (mod->flag.GW_FLOW){
#ifdef _DWGW_COUPLING
        if (mod->flag.SW_FLOW){ //GW-DW coupling since SW_FLOW is ALSO on
            printf("\n\nCalling ssw_alloc_init");
            printf("\nHacking into mod->grid->ndim. Setting it to 2; Don't forget to reset to 3 again!\n");
            mod->grid->ndim=2;
            ssw_alloc_init(&(mod->sw), mod->grid, mod->io, mod->flag);
            mod->grid->ndim=3;
            printf("\n----mod->grid->ndim reset to 3!\n");
        }
#endif
        printf("\n\nCalling sgw_alloc_init");
        sgw_alloc_init(&(mod->sgw), mod->grid, mod->io, mod->flag);
    }
#endif
    else {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> model physics not supported!");
    }
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /* define UNITS flag */
    tl_model_units(mod);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /* update the values for the time series */
    sseries_setall_ts_values(mod->series_head, mod->t_prev);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /* set string physics flags */
    sstr_value_set_physics(mod->str_values, mod->nstring, mod->flag);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /* set string roughness */
    sstr_value_set_roughness(mod->str_values, mod->nstring, mod->flag.SW2_FLOW, mod->flag.SW3_FLOW, mod->gravity, mod->manning_units_constant);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /* set strings total area MPI_REDUCE SUM */
    if (mod->grid->ndim == 3) {
        sstr_value_set_total_area(mod->str_values, mod->grid);
    }
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
#ifdef _ADH_GROUNDWATER
    if (mod->flag.GW_FLOW) {
      sgw_set_psk_series(mod);
    }
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    smodel_read_hot(mod);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#ifdef _ADH_GROUNDWATER
    if (mod->flag.GW_FLOW) {
        //sgw_evaluate_element_saturations(mod); /* Done in sgw_read_hot already. */
        sgw_evaluate_element_fluxes(mod);
    }
#endif

    /* open additional output files */
    if (mod->flag.WIND) {
        if (mod->file_output.wind) {
            sio_init_winds(mod->io);
        }
        if (mod->flag.WIND_STATION) {
            sseries_set_meteor_stations(mod->series_wind_head, mod->grid, WIND_SERIES);
        }
    }
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    if (mod->flag.WAVE) {
        if (mod->file_output.wave) {
            sio_init_waves(mod->io);
        }
        if (mod->flag.WAVE_STATION) {
            sseries_set_meteor_stations(mod->series_wave_head, mod->grid, WAVE_SERIES);
        }
    }
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    if (mod->flag.TRANSPORT) {
        sio_init_transport(mod->io, mod->ntransport);
    }
#ifdef _SEDIMENT
    if (mod->flag.SEDIMENT) {
        sio_init_sediment(mod->io, mod->nlayers, mod->nsed);
    }
#endif
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /* initialize dt */
    tc_init(mod);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
#ifdef _MESSG
    node_renumber(mod, 0);
    if (mod->grid->nelems3d>0) elem3d_renumber(mod->grid);
    elem2d_renumber(mod->grid);
    elem1d_renumber(mod->grid);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    if (mod->grid->type == COLUMNAR) build_columns(mod->grid, 1);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    if (mod->grid->type == COLUMNAR) node_renumber_surface(mod);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    if(mod->proc_flag == 1){
        comm_set_keys(mod->grid);
        partition_main(mod, 0);
#ifdef _DEBUG

        if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);

        /*
        if (mod->grid->ndim == 2) {
            sgrid_printScreen(mod->grid, mod->io->geo2d.filename);
        } else if (mod->grid->ndim == 3) {
            sgrid_printScreen(mod->grid, mod->io->geo3d.filename);
        }
        //tl_error("stop here\n");
        */
#endif
    }
#endif
    
    /* initialize total grid mass */
    if (screen_output.grid_mass_error == ON && (mod->proc_flag==1)) {
        if (mod->flag.SW2_FLOW) {
            mod->initial_grid_mass = tl_find_grid_mass_elem2d(mod->density, mod->str_values, mod->series_head, mod->sw->d2->head, mod->grid, mod->flag);
        } else if (mod->flag.SW3_FLOW) {
            mod->initial_grid_mass = tl_find_grid_mass_elem3d(mod->density, mod->grid, mod->sw->d3->displacement);
        }
#ifdef _MESSG
        mod->initial_grid_mass = messg_dsum(mod->initial_grid_mass, mod->grid->smpi->ADH_COMM);
#endif
    }
#ifdef _DEBUG
   if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif

    /* if sw3d, then calculate depth-averaged values for writing t=0 files*/
    if (mod->flag.SW3_FLOW == ON) {
        tl_calculate_depthavgvel(mod->grid, mod->sw->d3);
    }
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif

    // CJT :: if NS, we initialize pressure using a hydrostatic column
    /*
     if (mod->flag.NS3_FLOW == ON) {
     FILE *fpt = fopen("adh_out.hot","w");
     fprintf(fpt,"DATASET \n");
     fprintf(fpt,"OBJTYPE \"mesh3d\" \n");
     fprintf(fpt,"BEGSCL\n");
     fprintf(fpt,"ND %d\n",mod->grid->nnodes);
     fprintf(fpt,"NC %d\n",mod->grid->nelems3d);
     fprintf(fpt,"NAME IP\n");
     fprintf(fpt,"TS 0 0\n");
     int i, j, ncol = 6;
     for (i=0; i<mod->grid->nnodes; i++) {
     mod->ns->d3->prs[i] = 0.;
     printf("prs[%d]: %12.10f \n",i,mod->ns->d3->prs[i]);
     fprintf(fpt,"%20.10e\n",mod->ns->d3->prs[i]);
     for (j=1; j<ncol; j++) {
     mod->ns->d3->prs[i+j] = -mod->density * mod->gravity * (mod->grid->node[i+j].z - mod->grid->node[i].z);
     printf("prs[%d]: %12.10f \n",i+j,mod->ns->d3->prs[i+j]);
     fprintf(fpt,"%20.10e\n",mod->ns->d3->prs[i+j]);
     }
     i+=(ncol-1);
     }
     fprintf(fpt,"ENDDS\n");
     fclose(fpt);
     exit(-1);
     }
     */
    
    //print_grid_to_file(mod->grid,"INIT_DONE");
    //sgrid_check(mod->grid,__FILE__,__LINE__);
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
    //                                                                                              //
    // IMPORTANT NOTE: Add these lines AFTER node/element renumbering and partition_main(mod)       //
    //                 for MPI runs. Otherwise, the code will give errors / incorrect results.      //
    //                                                                                              //
#ifdef WINDLIB
    if (mod->flag.WIND_LIBRARY == ON){
        swindlib_init(mod->windlib, mod->grid->nnodes_sur);
        swindlib_firstread(mod, mod->windlib, mod->grid->nnodes_sur);
    }
#endif
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    // ABOVE LINES ADDED BY GAJANAN                                                                 //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //if(mod->grid->ndim ==3) ssw_3d_check_nan(mod->sw->d3, mod->grid);
    
    
#ifdef _ADH_GROUNDWATER
    // calculate initial groundwater saturations
    if (mod->flag.GW_FLOW == ON) {
        int i,ie;
        for (ie = 0; ie < mod->grid->nelems3d; ie++) {
            for (i=0; i < MAX_NNODES_ON_ELEM3D; i++) {
                mod->sgw->elem_3d_data[ie].saturation[i] = UNSET_FLT;
            }
        }
        sgw_evaluate_element_saturations(mod);
        for (ie = 0; ie < mod->grid->nelems3d; ie++) {
            for (i=0; i < mod->grid->elem3d[ie].nnodes; i++) {
                double sat = mod->sgw->elem_3d_data[ie].saturation[i];
                mod->sgw->elem_3d_data[ie].old_saturation[i] = sat;
            }
        }
    }
#endif

    smodel_check(mod);
    smodel_close_input(mod);
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif

    
//    MPI_Barrier(MPI_COMM_WORLD);
//    printf("\n");
//    int i;
//    for(i=0;i<mod->grid->nnodes;i++) {
//        if (fabs(mod->grid->node[i].z) < 1e-5) {
//            assert(fabs(mod->sgw->gw_phead[i] + 0.75)<1e-5);
//            assert(fabs(mod->sgw->old_gw_phead[i] + 0.75)<1e-5);
//            assert(fabs(mod->sgw->older_gw_phead[i] + 0.75)<1e-5);
//        } else {
//            assert(fabs(mod->sgw->gw_phead[i] + 10)<1e-5);
//            assert(fabs(mod->sgw->old_gw_phead[i] + 10)<1e-5);
//            assert(fabs(mod->sgw->older_gw_phead[i] + 10)<1e-5);
//        }
//    }
//    printf("myid: %d gw_phead[%d]: %f\n",mod->grid->smpi->myid,mod->grid->node[i].gid+1,mod->sgw->gw_phead[i]);
//    MPI_Barrier(MPI_COMM_WORLD);
//    tl_error("temp");
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smodel_read_hot(SMODEL *mod) {

#ifdef _DEBUG
    int myid = 0, ierr = UNSET_INT;;
#ifdef _MESSG
        ierr = MPI_Comm_rank(cstorm_comm, &myid);
#endif
      printf("\n---- MYID %d :: MPI_COMM_WORLD_RANK: %d MODEL READING HOTSTART %s\n", mod->grid->smpi->myid,myid,mod->io->hot.filename);
#endif
    
    if (mod->flag.SW_FLOW && !mod->flag.GW_FLOW) {
        ssw_read_hot(mod);
    } else if (mod->flag.NS_FLOW) {
        sns_read_hot(mod);
    }
#ifdef _ADH_GROUNDWATER
    else if (mod->flag.GW_FLOW){
#ifdef _DWGW_COUPLING
        if (mod->flag.SW_FLOW){ //GW-DW Coupling is SW_FLOW is ALSO on
            printf("\nHacking into mod->grid->ndim. Setting it to 2; Don't forget to reset to 3 again!\n");
            mod->grid->ndim=2;
            ssw_read_hot(mod);
            mod->grid->ndim=3;
            printf("\n----mod->grid->ndim reset to 3!\n");
    
            sgw_read_hot(mod);
        }
#endif 
      sgw_read_hot(mod);
    }
#endif 
    else {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> model physics not supported!");
    }
    
#ifdef _SEDIMENT
    if (mod->flag.SEDIMENT) {
        ssediment_read_hot(mod);
    }
#endif

#ifdef _DEBUG
        printf("\n---- MYID %d :: MPI_COMM_WORLD_RANK: %d MODEL FINISHED READING HOTSTART %s\n", mod->grid->smpi->myid,myid,mod->io->hot.filename);
#endif
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smodel_print(SMODEL *mod)
{
    
    printf("\n****************************************************\n");
    if (mod->flag.SW2_FLOW) {
        printf("MODEL: 2D - SHALLOW WATER");
        if (mod->flag.SW2_TRANSPORT)
            printf(" - with transport\n");
    }
    else if (mod->flag.SW3_FLOW) {
        printf("MODEL: 3D - SHALLOW WATER\n");
        if (mod->flag.SW3_TRANSPORT)
            printf(" - with transport\n");
    }
#ifdef _ADH_GROUNDWATER
    else if (mod->flag.GW_FLOW) {
        printf("MODEL: 3D - GROUNDWATER\n");
        if (mod->flag.GW_TRANSPORT)
            printf(" - with transport\n");
    }
#endif 
    else
        printf("MODEL: UNRECOGNIZED\n");
    
    printf("\nModel Flags ----------------------------------------\n");
    printf("moving_grid: %d\n", mod->flag.MG);
    printf("adaption: %d\n", mod->flag.GRID_ADAPTION);
    printf("icm: %d\n", mod->flag.ICM);
    printf("nsm: %d\n", mod->flag.NSM);
    printf("wave: %d\n", mod->flag.WAVE);
    printf("wind: %d\n", mod->flag.WIND);
    printf("coriolis: %d\n", mod->flag.CORIOLIS);
    printf("muc: %d \n", mod->flag.MUC);
    printf("conveyance: %d \n", mod->flag.CONVEYANCE);
    printf("ice: %d \n", mod->flag.ICE);
    printf("tidal_flag: %d\n", mod->flag.TIDE);
    printf("o_flag: %d \n", mod->o_flag);
    
    printf("\nSolver Variables  -----------------------------------\n");
    printf("inc_nonlin; %f\n", mod->inc_nonlin);
    printf("tol_nonlin: %f\n", mod->tol_nonlin);
    printf("max_nonlin_it: %d\n", mod->max_nonlin_it);
    printf("nalloc_inc: %d\n", mod->nalloc_inc);
    printf("nblock: %d\n", mod->nblock);
    
    printf("\nWetting/Drying Variables  ----------------------------\n");
    printf("drying_lower_limit: %f\n", mod->drying_lower_limit);
    printf("drying_upper_limit: %f\n", mod->drying_upper_limit);
    printf("wd_lower_tol: %f\n", mod->wd_lower_tol);
    printf("wd_upper_tol: %f\n", mod->wd_upper_tol);
    printf("wd_rate_lower_tol: %f\n", mod->wd_rate_lower_tol);
    printf("wd_rate_upper_tol: %f\n", mod->wd_rate_upper_tol);
    
    printf("\nModel Time Variables  -------------------------------\n");
    printf("t_init: %f\n", mod->t_init);
    printf("t_prev: %f\n", mod->t_prev);
    printf("t_final: %f\n", mod->t_final);
    printf("tau_temporal: %f\n", mod->tau_temporal);
    printf("out_level: %d\n", mod->out_level);
    
    printf("\nTime-Series: %d  -------------------------------------\n", mod->nseries);
    sseries_printScreen(*(mod->series_dt), 1);
    sseries_printScreen(*(mod->series_out), 1);
    if (mod->series_wind_head != NULL) sseries_printScreen_list(1, mod->series_wind_head);
    if (mod->series_wave_head != NULL) sseries_printScreen_list(1, mod->series_wave_head);
    if (mod->series_head != NULL) sseries_printScreen_list(1, mod->series_head);
    
    printf("\nStrings  ---------------------------------------------\n");
    printf("nstring: %d\n", mod->nstring);
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smodel_open_input(SMODEL * mod) {
    
    char msg[MAXLINE];
    assert(mod->io);    /* should be valid */
    
    root_print("\n ---------\n Run Files\n ---------");
    
    /* super file */
    int super = strlen(mod->io->sup.filename);  /* whether super file exists */
    if (super) {
        sprintf(msg, " Super file: %s - FOUND", mod->io->sup.filename);
        root_print(msg);
    }
    
    /* fundamental AdH input files */
    open_input_file(&(mod->io->bc), "boundary condition file", super);
    open_input_file(&(mod->io->hot), "hotstart file", super);
    
    // check which geo file exists, if both 2D and 3D exists, open 3D
    if (doesFileExist(mod->io->geo2d.filename) && doesFileExist(mod->io->geo3d.filename)) {
        open_input_file(&(mod->io->geo3d), "3d geometry file", super);
        open_input_file( &(mod->io->face), "3d boundary face file", super);
        return;
    } else if (doesFileExist(mod->io->geo2d.filename)) {
        open_input_file(&(mod->io->geo2d), "2d geometry file", super);
        return;
    } else if (doesFileExist(mod->io->geo3d.filename)) {
        open_input_file(&(mod->io->geo3d), "3d geometry file", super);
        open_input_file( &(mod->io->face), "3d boundary face file", super);
        return;
    } else {
        tl_error("2d or 3d geo file needed\n");
    }
}

void smodel_close_input(SMODEL * mod) {

    fclose(mod->io->bc.fp);
    fclose(mod->io->hot.fp);
    
    if (doesFileExist(mod->io->geo2d.filename) && doesFileExist(mod->io->geo3d.filename)) {
        fclose(mod->io->geo3d.fp);
        fclose(mod->io->face.fp);
        return;
    } else if (doesFileExist(mod->io->geo2d.filename)) {
        fclose(mod->io->geo2d.fp);
        return;
    } else if (doesFileExist(mod->io->geo3d.filename)) {
        fclose(mod->io->geo3d.fp);
        fclose(mod->io->face.fp);
        return;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smodel_open_output(SMODEL * mod)
{
    
    assert(mod->io);
    if (mod->flag.SW_FLOW && !mod->flag.GW_FLOW) { // includes diffusive wave
        ssw_open_output(mod);
    }
    
    if (mod->flag.NS_FLOW) {
        sns_open_output(mod);
    }
    
    if (mod->flag.TRANSPORT) {
        scon_open_output(mod);
    }
    
#ifdef _SEDIMENT
    if (mod->flag.SEDIMENT) {
        ssediment_open_output(mod);
    }
#endif
    
#ifdef _ADH_GROUNDWATER
    if (mod->flag.GW_FLOW){
#ifdef _DWGW_COUPLING
        if (mod->flag.SW_FLOW){ //GW-DW Coupling since SW_FLOW is ALSO on
            printf("\nHacking into mod->grid->ndim. Setting it to 2; Don't forget to reset to 3 again!\n");
            mod->grid->ndim=2;
            ssw_open_output(mod);
            mod->grid->ndim=3;
            printf("\n----mod->grid->ndim reset to 3!\n");
        }
#endif
        sgw_open_output(mod);
    }

#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smodel_print_ts(SMODEL *mod, double time) {
    
    int i, j;                        /* loop counter */
    int **ndata, **ndata_sur;
    int ierr;
    int *my_nnode_ext, *my_nnode_ext_sur;             /* external processor nodal information */
    int my_nnode_max;
    int my_nnode_max_sur;
    int my_nnode;
    int my_nnode_sur;
    
    SGRID *grid = mod->grid;
    ndata = (int **) tl_alloc(sizeof(int *), mod->grid->smpi->npes);
    ndata_sur = (int **) tl_alloc(sizeof(int *), mod->grid->smpi->npes);
    
    my_nnode = grid->my_nnodes;
    my_nnode_sur =  grid->my_nnodes_sur;
    
    my_nnode_ext = (int *) tl_alloc(sizeof(int), grid->smpi->npes);
    my_nnode_ext_sur = (int *) tl_alloc(sizeof(int), grid->smpi->npes);
    
    /*set global node ID array for printing */
#ifdef _MESSG
    MPI_Status *msg_status= grid->smpi->msg_status;
    
    ierr = MPI_Reduce(&my_nnode, &my_nnode_max, 1, MPI_INT, MPI_MAX, 0, mod->grid->smpi->ADH_COMM);
    if (ierr != MPI_SUCCESS)
        messg_err(ierr);
    ierr = MPI_Reduce(&my_nnode_sur, &my_nnode_max_sur, 1, MPI_INT, MPI_MAX, 0, mod->grid->smpi->ADH_COMM);
    if (ierr != MPI_SUCCESS)
        messg_err(ierr);
    if(mod->grid->smpi->myid != 0) {
        my_nnode_max=my_nnode;
        my_nnode_max_sur=my_nnode_sur;
    }
#else
    my_nnode_max = my_nnode;
    my_nnode_max_sur = my_nnode_sur;
#endif
    
#ifdef _MESSG
    ierr = MPI_Allgather(&my_nnode, 1, MPI_INT, my_nnode_ext, 1, MPI_INT, mod->grid->smpi->ADH_COMM);
    if (ierr != MPI_SUCCESS) {
        messg_err(ierr);
    }
    ierr = MPI_Allgather(&my_nnode_sur, 1, MPI_INT, my_nnode_ext_sur, 1, MPI_INT, mod->grid->smpi->ADH_COMM);
    if (ierr != MPI_SUCCESS) {
        messg_err(ierr);
    }
#else
    my_nnode_ext[0] = my_nnode;
    my_nnode_ext_sur[0] = my_nnode_sur;
    
#endif
    
    if (grid->smpi->myid <= 0) {
        
        for (i=0;i<mod->grid->smpi->npes;i++){
            ndata[i] = (int *) tl_alloc(sizeof(int), my_nnode_max);
        }
        
        for (i=0;i<mod->grid->smpi->npes;i++){
            ndata_sur[i] = (int *) tl_alloc(sizeof(int), my_nnode_max_sur);
        }
        
        j=0;
        for (i=0;i<my_nnode;i++){
            ndata[0][i]=grid->node[i].gid;
            if (grid->type == COLUMNAR) {
                if(grid->ndim==3){
                    if(grid->nodeID_3d_to_2d_sur[i] != UNSET_INT) ndata_sur[0][j++]=grid->node[i].global_surf_id;
                }else{
                    ndata_sur[0][j++]=grid->node[i].global_surf_id;
                }
            }
        }
        
#ifdef _MESSG
        for (i=1;i<mod->grid->smpi->npes;i++){
            ierr = MPI_Recv(ndata[i], my_nnode_max, MPI_INT, i, 998, mod->grid->smpi->ADH_COMM, msg_status);
            if (ierr != MPI_SUCCESS)
                messg_err(ierr);
            ierr = MPI_Recv(ndata_sur[i], my_nnode_max_sur, MPI_INT, i, 999, mod->grid->smpi->ADH_COMM, msg_status);
            if (ierr != MPI_SUCCESS)
                messg_err(ierr);
        }
#endif
    }else{
#ifdef _MESSG

        ndata[0] = (int *) tl_alloc(sizeof(int), my_nnode);
        ndata_sur[0] = (int *) tl_alloc(sizeof(int), my_nnode_sur);
        
        j=0;
        
        for (i = 0; i < my_nnode; i++) {
            ndata[0][i] = grid->node[i].gid;
            if (grid->type == COLUMNAR) {
                if(grid->ndim==3){
                    if(grid->nodeID_3d_to_2d_sur[i] != UNSET_INT) ndata_sur[0][j++]=grid->node[i].global_surf_id;
                }else{
                    ndata_sur[0][j++]=grid->node[i].global_surf_id;
                }
            }
        }
        
        ierr = MPI_Send(ndata[0], my_nnode, MPI_INT, 0, 998, mod->grid->smpi->ADH_COMM);
        if (ierr != MPI_SUCCESS)
            messg_err(ierr);
        ierr = MPI_Send(ndata_sur[0], my_nnode_sur, MPI_INT, 0, 999, mod->grid->smpi->ADH_COMM);
        if (ierr != MPI_SUCCESS)
            messg_err(ierr);
#endif
    }
    
    double tuse = mod->t_prev + mod->dt;
    int it1 = (int) (tuse);
    int it2 = (int) (1000.0 * (tuse - (double) (it1)));
    
#ifdef _ADH_HDF5
    if (mod->flag.PC_FILE_XDMF != ON){
#endif
        if (mod->proc_flag == 1) {
        if (mod->file_output.adapt_grid) {
            sgrid_print_adapted_ts(mod->grid, mod->io->sup, mod->io->proj_name, it1, it2, ndata, my_nnode_max, my_nnode_ext);
        }
        
        if (mod->flag.SW_FLOW && !mod->flag.GW_FLOW) {
            ssw_print_ts(mod->sw, mod->io, mod->grid, time, mod->series_out->outfact, mod->flag, mod->io->sup, mod->io->proj_name, it1, it2, ndata, my_nnode_max, my_nnode_ext, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur, mod->file_output.adapt_sw,mod->file_output);
        }
        
        if (mod->flag.NS_FLOW) {
            sns_print_ts(mod->ns, mod->io, mod->grid, time, mod->series_out->outfact, mod->flag, mod->io->sup, mod->io->proj_name, it1, it2, ndata, my_nnode_max, my_nnode_ext, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur, mod->file_output.adapt_ns,mod->file_output);
        }
        
        if (mod->flag.TRANSPORT) {
            scon_print_ts(mod->ntransport, mod->con, mod->io, mod->grid, time, mod->series_out->outfact, mod->flag, mod->io->sup, mod->io->proj_name, it1, it2, ndata, my_nnode_max, my_nnode_ext, mod->file_output.adapt_con);
        }
        
#ifdef _SEDIMENT
        if (mod->flag.SEDIMENT) {
            ssediment_print_ts(mod->grid, mod->io, mod->sed, time, mod->series_out->outfact, mod->grid->ndim, mod->io->sup, mod->io->proj_name, it1, it2, ndata, my_nnode_max, my_nnode_ext, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur, mod->file_output.adapt_sed);
        }
#endif
        }
        
#ifdef _ADH_GROUNDWATER
        if (mod->flag.GW_FLOW) {
#ifdef _DWGW_COUPLING
            if (mod->flag.SW2_FLOW) { // DW-GW Coupling since SW2_FLOW is ALSO on
                //printf("\nHacking into mod->grid->ndim. Setting it to 2; Don't forget to reset to 3 again!\n");
                mod->grid->ndim=2;
                ssw_print_ts(mod->sw, mod->io, mod->grid, time, mod->series_out->outfact, mod->flag, mod->io->sup, mod->io->proj_name, it1, it2, ndata, my_nnode_max, my_nnode_ext, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur, mod->file_output.adapt_sw,mod->file_output);
                mod->grid->ndim=3;
                //printf("\n----mod->grid->ndim reset to 3!\n");
            }
#endif
            sgw_print_ts(mod, mod->sgw, mod->io, mod->grid, time, mod->series_out->outfact, mod->flag, mod->io->sup, mod->io->proj_name, it1, it2, ndata, my_nnode_max, my_nnode_ext, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur, mod->file_output.adapt_gw);

        }
#endif 

#ifdef _ADH_HDF5
    }
    else {
        //printf("\n Time Outfact = %f",mod->series_out->outfact);
        if (mod->proc_flag == 1) {
            ps_print_xdmf(mod, time, mod->series_out->outfact);
        }
    }
#endif
    
    if(grid->smpi->myid <= 0) {
        for (i=0;i<grid->smpi->npes;i++){
            ndata[i] = (int *) tl_free(sizeof(int ), my_nnode_max, ndata[i]);
            ndata_sur[i] = (int *) tl_free(sizeof(int ), my_nnode_max_sur, ndata_sur[i]);
        }
    }else{
        ndata[0] = (int *) tl_free(sizeof(int), my_nnode, ndata[0]);
        ndata_sur[0] = (int *) tl_free(sizeof(int), my_nnode_sur, ndata_sur[0]);
    }
    
    ndata = (int **) tl_free(sizeof(int *), mod->grid->smpi->npes, ndata);
    ndata_sur = (int **) tl_free(sizeof(int *), mod->grid->smpi->npes, ndata_sur);
    my_nnode_ext = (int *) tl_free(sizeof(int), grid->smpi->npes, my_nnode_ext);
    my_nnode_ext_sur = (int *) tl_free(sizeof(int), grid->smpi->npes, my_nnode_ext_sur);
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smodel_check(SMODEL *mod)
{
    // make sure that bottom friction is nonzero for all 2d bed faces when sedlib is turned on
    int istring = UNSET_INT;
    if (mod->flag.SEDLIB == ON) {
        for (istring=0; istring<mod->nstring; istring++) {
            printf("istring %d nstring %d flag %d BCT_BED %d roughness %f \n", istring, mod->nstring, mod->str_values[istring].flow.bc_flag, BCT_BED, mod->str_values[istring].roughness );
            if ((mod->str_values[istring].flow.bc_flag == BCT_BED) && (fabs(mod->str_values[istring].roughness)<1e-12)) {
                tl_error(">> SEDLIB requires a non-zero friction on all bed faces.");
            }
        }
    }
    
    if (mod->flag.SW_FLOW + mod->flag.NS_FLOW + mod->flag.GW_FLOW == 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> must define a type of model physics via the OP card!");
    }
    
    if (mod->flag.SW2_FLOW + mod->flag.SW3_FLOW> 1) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> only one type of physics can be given for each bc file!");
    }
    
    if (mod->flag.SW2_FLOW + mod->flag.SW3_FLOW + mod->flag.GW_FLOW> 2) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> Cannot use multiple physics in the same bc file, except for OP DIF and OP GW!");
    }
    
    if ((mod->inc_nonlin) < 0 && (mod->tol_nonlin < 0)) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> no tolerance set!");
    }
    
    if (mod->flag.OUTPUT != ON) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> No output print control specified.  You must use an OS or OC card.");
    }
    
    if (mod->flag.TRANSPORT && mod->ntransport < 1) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> transport flag on, but there are 0 constituents\n");
    }
    
    /* check materials */
    smat_checkall(mod);
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smodel_free(SMODEL *mod) {
    
    int i, k, imod, nnodes = 0;
    int nelems=0;
    
#ifdef _DEBUG
    printf("\n");
#endif
    
    if (mod != NULL) {
        if (mod->sw != NULL) {
            nnodes = mod->grid->nnodes;
        }
    }
    
    if (mod->bc_mask != NULL) {
        mod->bc_mask = (int *) tl_free(sizeof(int), mod->grid->nnodes_matrix * mod->max_nsys, mod->bc_mask);
    }
    
    /* free solver info node block (should be done in Solver_Info (free) ) */
    // CJT :: GJ COMMENTED THIS OUT, NOT SURE IF IT SHOULD BE
    // gkc :: 2020-04-16: I commented this out because this is now being done in ssuperModel.c.
    //    if (mod->sw != NULL) {
    //        if ((mod->sw->d2 != NULL || mod->sw->d3 != NULL) && (mod->solver_info.node_block != NULL)) {
    //            mod->solver_info.node_block = (int *) tl_free(sizeof(int), mod->grid->nnodes_matrix, mod->solver_info.node_block);
    //        }
    //    }
    //
    //    solv_bcgstab_clean();
    //    solv_blk_set_clean(&(mod->solver_info.profile));
    // gkc :: 2020-04-16: However, choosing to keep the groundwater free here for now. This might also have to be
    //                    commented out and  added to ssuperModel.c. Will check some day and update this comment I
    //                    hope.
#ifdef _ADH_GROUNDWATER
        if (mod->sgw != NULL) {
            if (mod->solver_info.node_block != NULL) {
                mod->solver_info.node_block = (int *) tl_free(sizeof(int), mod->grid->nnodes_matrix, mod->solver_info.node_block);
            }
        }
#endif
    
#ifdef _UMFPACK
    solv_blk_free_sparse();
#endif
    
    /* free file struct */
    if(mod->io != NULL) sio_free(mod->io, mod->ntransport, mod->nlayers, mod->nsed);
    
    /* free materials */
    if (mod->mat != NULL) smat_free(mod, mod->mat, mod->nmat);
    
    /* free strings */
    if(mod->str_values != NULL) sstr_value_free(mod->str_values, mod->nstring, mod->ntransport, mod->nsed);
    
    /* free structures */
#ifdef _ADH_STRUCTURES
    if (mod->nweir > 0)
        sstructures_free_weirs(mod->weir, mod->nstring, mod->nweir);
    if (mod->nflap > 0)
        sstructures_free_flaps(mod->flap, mod->nstring, mod->nflap);
    if (mod->nsluice > 0)
        sstructures_free_sluices(mod->sluice, mod->nstring, mod->nsluice);
#endif
    
    /* free series linked list */
    if(mod->series_wind_head != NULL) sseries_free_list(&(mod->series_wind_head));
    if(mod->series_wave_head != NULL) sseries_free_list(&(mod->series_wave_head));
    if(mod->series_head != NULL)      sseries_free_list(&(mod->series_head));
#ifdef _ADH_GROUNDWATER
    if(mod->series_gw_psk_head != NULL) sseries_free_list(&(mod->series_gw_psk_head));
#endif
    
    /* free individual series */
    if(mod->series_dt != NULL)  sseries_free(mod->series_dt);
    if(mod->series_out != NULL) sseries_free(mod->series_out);
    
    /* free physics stuct */
    if (mod->flag.SW_FLOW == ON && mod->sw != NULL) {
        ssw_free(mod->sw, mod->grid, mod->flag);
    }
    if (mod->flag.NS_FLOW == ON && mod->ns != NULL) {
        sns_free(mod->ns, mod->grid, mod->flag);
    }
#ifdef _ADH_GROUNDWATER
    if (mod->flag.GW_FLOW == ON && mod->sgw != NULL) {
        sgw_free(mod->sgw, mod->grid, mod->flag);
    }
#endif

    /* free constituents */
    if (mod->flag.TRANSPORT == ON && mod->con != NULL) {
        
        if (mod->flag.SW2_FLOW || mod->flag.NS2_FLOW) {
            nelems = mod->grid->max_nelems2d;
        } else {
            nelems = mod->grid->max_nelems3d;
        }
        scon_free(mod->grid->max_nnodes, nelems, mod->ntransport, mod->con);
    }
    
#ifdef _SEDIMENT
    /* free sediment constituents */
    if (mod->flag.SEDIMENT == ON && mod->sed != NULL) {
        ssediment_free(mod->sed);
    }
#endif

#ifdef WINDLIB
    if (mod->flag.WIND_LIBRARY==ON) {
        swindlib_free(mod->windlib, mod->grid->nnodes_sur);
    }
#endif
    
    mod->fmap = NULL;        // alias reset to NULL; freed in ssuperModel.c
    mod->fmap_wvel = NULL;   // alias reset to NULL; freed in ssuperModel.c
    
    /* free grid */
    if (mod->grid != NULL) {
        
#ifdef _MESSG
        mod->grid->supersmpi = NULL;   // alias reset to NULL; freed in ssuperModel.c
        mod->grid->supersmpi_wvel = NULL;   // alias reset to NULL; freed in ssuperModel.c
#endif
        
#ifdef _DEBUG
        printf("... MYID %d freeing grid ...\n", mod->grid->smpi->myid);
#endif
        sgrid_free(mod->grid);
        
#ifdef _DEBUG
        printf("... freeing grid ...\n");
#endif
    }
}

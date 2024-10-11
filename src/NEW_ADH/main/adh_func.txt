#include "adh.h"

// -----------------------------------------------------
// CSTORM ----------------------------------------------
#ifdef _MESSG
extern MPI_Comm cstorm_comm;
#endif
extern SMODEL *mod_cstorm;
static int CSTORM_FLAG = 0;
static int DEBUG = OFF;

// -----------------------------------------------------
// FOR MODEL MONOLITHIC AN FLUX COUPLING ---------------
static int *my_nelem_ext;
static int my_nelem_max, my_nelem;
static int **ndata;
void update_dt(SSUPER_MODEL *sm);
void gather_fluxes(SMODEL *, int, int);
void set_flux_messaging(SMODEL *);

// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------

int adh_init_func_(SMODEL_DESIGN dmod, int *argc, char ***argv
#ifdef _MESSG
                   ,MPI_Comm *comm
#endif
) {
    int i,j;
    
    // Command line arguments
    int input_test = OFF;
    if ( (*argc) == 3) {
        if (strcmp((*argv)[2],"-s") != AGREE) { // sanity check option
            fprintf(stderr, "COMMAND LINE ERROR: Third command line argument must be -s.\n");
            exit(-1);
        } else {
            printf("\n");
            printf("********************************************\n");
            printf("Running AdH only for file input checking ...\n");
            printf("********************************************\n");
            input_test = ON;
        }
    }
    
    // Get AdH root filename
    char adh_root[MAXLINE];
    if (*argc == 0) {
        FILE *fp_cstorm = io_fopen("adh.sup", "r", FALSE);
        if (fp_cstorm == NULL) {
            tl_error("Error :: The -cstorm input flag is used with binary, please provide adh.sup file\n");
        }
        fscanf(fp_cstorm, "%s", adh_root);
        fclose(fp_cstorm);
        CSTORM_FLAG = 1;
    } else {
        strcpy(adh_root, (*argv)[1]);
    }
    
    // Initialize AdH MPI
    int myid = 0, npes = 1, ierr_code = UNSET_INT;
#ifdef _MESSG
    int flag;
    ierr_code = MPI_Initialized(&flag); // has it been initialized?
    if(ierr_code != MPI_SUCCESS) messg_err(ierr_code);
    if (!flag) {
        ierr_code = MPI_Init(argc, argv);
        if(ierr_code != MPI_SUCCESS) messg_err(ierr_code);
        MPI_Comm_dup(MPI_COMM_WORLD, &cstorm_comm);
    } else {
        MPI_Comm_dup(*comm, &cstorm_comm);
    }
    if(ierr_code != MPI_SUCCESS) {
        messg_err(ierr_code);
    }
    ierr_code = MPI_Comm_rank(cstorm_comm, &myid);
    ierr_code = MPI_Comm_size(cstorm_comm, &npes);
#endif

    // Initialize PETSc
#ifdef _PETSC
    static char help[] = "PETSc for AdH.\n";
    ierr = PetscInitialize(argc,argv,(char*)0,help);if (ierr) return ierr;
#endif
    
    // Initialize Debug
#ifdef _MESSG
    debug_initialize(cstorm_comm);
#else
    debug_initialize();
#endif
    sdebug_init(&(debug));
    
    // Initialize test case parameters
    stestcase_init();
    
    // initialize simulation screen output options
    sscreen_output_init();
    
    // sets up the linked lists
    tl_list_setup();
    
    // Allocate and read Designer Model grid
    void sgrid_read(&(dmod->grid),adh_root
#ifdef _MESSG
                    ,cstorm_comm
#endif
    );
    
    // Initialize Designer Model Super Models
    smodel_super_alloc_init(&(dmod->superModel),dmod->nSuperModels);
    
    
    // Stop code if only file checking
    if (input_test == YES) {
        tl_error("Exiting after file check!");
    }
    
    // Open all output files
    mod = NULL;
    for (i=0; i<dmod->nsupermodels; i++) {
        if(dmod->grid->smpi->myid <= 0) smodel_super_open_output(&(dmod->superModel[i]));
    }
    
    // Prep test case solutions if used
    if (dmod->test_case_flag.on) {
        for (i=0; i<dmod->nsupermodels; i++) {
            testcase_prep(&(dmod->superModel[i]));
        }
    }
    
    // Write the initial grid and conditions
    for (i=0; i<dmod->nsupermodels; i++) {
#ifdef _ADH_HDF5
        if (dmod->superModel[i]->flag.PC_FILE_XDMF) {
            xdmf_init(dmod->superModel[i], dmod->grid->smpi->npes, dmod->grid->smpi->myid);
        }
#endif
        // write initial conditions (either zero or read via hot
        smodel_print_ts(dmod->superModel[i], dmod->superModel[i]->t_prev);
    }
    
    // cjt :: for now, cstorm only couples to one model runs
    if (CSTORM_FLAG == 1) {
        if (dmod->nsupermodels > 1) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> Only 1 supermodel can be ran with CSTORM at this time.");
        }
        mod_cstorm = &(dmod->superModel[0]);
    }
    
    return ierr;
}

// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------

int adh_finalize_func_(SSUPER_MODEL **superModel_ptr, SSUPER_INTERFACE **superInterface_ptr, int nsupermodels, int nsuperinterfaces) {
    
    int i,j,k,ip,ierr = 0;
    SSUPER_MODEL *sm=(*superModel_ptr); /* alias */
    SMODEL *mod;
    
    /* cjt :: write total linear and nonlinear iterations */
    printf("\n \n");
    for (i=0; i<nsupermodels; i++) {
        if(sm[i].proc_flag == 1){
            for (j=0; j<sm[i].nsubmodels; j++) {
                if(sm[i].submodel[j].proc_flag==1){
                    mod = &(sm[i].submodel[j]);
                    
                    if(ndata != NULL){
                        if (mod->grid->smpi->myid==0) {
                            for(ip=0;ip<mod->grid->smpi->npes;ip++)
                                ndata[ip] = (int *) tl_free(sizeof(int), my_nelem_max, ndata[i]);
                        } else {
                            ndata[0] = (int *) tl_free(sizeof(int), my_nelem, ndata[0]);
                        }
                        ndata = (int **) tl_free(sizeof(int *), mod->grid->smpi->npes, ndata);
                        my_nelem_ext = (int *) tl_free(sizeof(int), mod->grid->smpi->npes, my_nelem_ext);
                    }
#ifdef _ADH_HDF5
                    if (mod->proc_flag) {
                    if (mod->flag.PC_FILE_XDMF) {
                        xdmf_finalize(&(mod->hdf5), mod->io, mod->grid->smpi->npes, mod->grid->smpi->myid,
                                (mod->flag.SW3_FLOW
#ifdef _ADH_GROUNDWATER
                                 || mod->flag.GW_FLOW
#endif
                                )
                                );
                    }
                    }
                    else {
                        if (mod->grid->smpi->myid <= 0) print_xms_trailers(mod);
                    }
#else
                    if(mod->grid->smpi->myid <= 0) print_xms_trailers(mod);
#endif
                    
                    // cjt :: write total linear and nonlinear iterations
                    if (mod->grid->smpi->myid == 0) {
                        if (sm[i].nsubmodels == 1) {
                            if (sm[i].submodel[j].flag.SW2_FLOW) {
                                printf("model: %s :: total number of HYD iterations: %d\n",sm[i].submodel[j].name, mod->nonlinear_it_total);
                            } else if (sm[i].submodel[j].flag.SW3_FLOW) {
                                printf("model: %s :: total number of HVEL nonlinear iterations: %d\n",sm[i].submodel[j].name, mod->nonlinear_it_total_hvel);
                                printf("model: %s :: total number of WVEL nonlinear iterations: %d\n",sm[i].submodel[j].name, mod->nonlinear_it_total_wvel);
                            }
                        } else if (sm[i].nsubmodels > 1) {
                            if (sm[i].submodel[j].flag.SW2_FLOW) {
                                printf("model: %s :: total number of HYBRID nonlinear iterations: %d\n",sm[i].submodel[j].name, sm->nonlinear_it_total);
                            } else if (sm[i].submodel[j].flag.SW3_FLOW) {
                                printf("model: %s :: total number of HYBRID nonlinear iterations: %d\n",sm[i].submodel[j].name, sm->nonlinear_it_total_hvel);
                                printf("model: %s :: total number of WVEL nonlinear iterations: %d\n",sm[i].submodel[j].name, sm->nonlinear_it_total_wvel);
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    /* free test case memory if used */
    for (i=0; i<nsupermodels; i++) {
        if (test_case_flag.on) {
            testcase_clean(&(sm[i]));
        }
    }
    
    /* free models */
    sm = NULL;
    
    ssuperModel_free(superModel_ptr, superInterface_ptr, nsupermodels, nsuperinterfaces);
    
    /* these are globally allocated at the moment */
    solv_bcgstab_clean();
#ifdef _UMFPACK
    solv_blk_free_sparse();
#endif
    
    /* finalize debugging */
    debug_finalize();
    
    return ierr;
}

// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------

int adh_run_func_(SSUPER_MODEL *superModel, SSUPER_INTERFACE *superInterface,
                  int nsupermodels, int nsuperinterfaces, double *run_time) {
    
    int i,j,k, ie;
    int isuper_model = UNSET_INT, isuperinterface = UNSET_INT, ierr=0, isers=0;
    int sum_finished_models=0, isubmodel=UNSET_INT, was_grid_adapted=0, myid = 0, npes = 1;
    int dt_store[nsupermodels];
    int model_finish[nsupermodels];
    SSUPER_MODEL *sm = NULL;
    SSUPER_INTERFACE *si = NULL;
    SMODEL *mod = NULL;
    
    // initialize all models to unfinished
    sarray_init_int(model_finish, nsupermodels);
    
#ifdef _MESSG
    MPI_Comm_rank(cstorm_comm, &myid);
    MPI_Comm_size(cstorm_comm, &npes);
    int was_grid_adapted_max=0;
#endif
    // cjt :: i orginally added below, but sent back to the above
    //mod = &(superModel[0].submodel[0]);
    //myid = mod->grid->smpi->myid;
    
    // obtain final run times
    for (isuper_model=0; isuper_model<nsupermodels; isuper_model++) {
        assert(&(superModel[isuper_model]) != NULL);
        sm = &(superModel[isuper_model]);
        for(isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++){
            mod = &(sm->submodel[isubmodel]);
            if (CSTORM_FLAG == 1) {
                mod->t_final_store = mod->t_final;
                mod->t_final = mod->t_prev + *run_time;
                mod->flag.CSTORM_WSID = 100; // anything but off for now
            } else {
                // check to make sure all submodels have same final time
                if (isubmodel > 1 && fabs(mod->t_final-mod->t_final_store) > 1e-6){
                    tl_error("ERROR :: All submodel simulation end times (TC FINAL) must be consistent!");
                }
                mod->t_final_store = mod->t_final;
            }
        }
    }
    
    
    if (myid == 0) {
        printf("\n**************************************************************************************************");
        printf("\n************************************** ADH SIMULATION START **************************************\n");
    }
    
    
    for (isuper_model=0; isuper_model<nsupermodels; isuper_model++) { sm->dt = 100000000.; }

    
    // ---------------------------------------------------
    // time loop -----------------------------------------
    do {
        
        sum_finished_models = 0;          // keep track of how many models have finished
        int fluxcouplingconverged = NO;   // interative flux control
        
        // ---------------------------------------------------
        // flux convergence loop -----------------------------
        do {
            assert (nsupermodels<3);
            
            // ---------------------------------------------------
            // supermodel loop -----------------------------------
            for (isuper_model=0; isuper_model<nsupermodels; isuper_model++) {
                sm = &(superModel[isuper_model]);
                
                mod = &(sm->submodel[0]); // non-fully coupled models only
                if (CSTORM_FLAG == 1) {
                    dt_store[isuper_model] = mod->dt;
                    if (*run_time < mod->dt) { // shorten dt for CSTORM coupling if run-time is less
                        mod->dt = *run_time;
                    }
                }
                mod = NULL;
                
                if (model_finish[isuper_model] == 1) continue;
                if (sm->proc_flag==0) model_finish[isuper_model] = 1; //processor is not assigned to this supermodel
                
                
                
                // Solve one finite element time-step implicitly
                do {
                    update_dt(sm);
                    // Loop over all submodels in the current supermodel
                    for (isubmodel=0; isubmodel<sm->nsubmodels; isubmodel++){
                        mod = &(sm->submodel[isubmodel]);
                        if(mod->proc_flag==1){
                            // update time-series values
                            sseries_setall_ts_valuesAVG(mod->series_head, mod->t_prev, mod->t_prev + mod->dt);
                            
                            // update meterologic data
                            meteor_update(mod, mod->t_prev);
                        }
                    }
                    
                } while (fe_main(sm) == NO);
                
            }
            
            fluxcouplingconverged = YES; /* For now, till flux coupling is set-up. */
            if (nsuperinterfaces > 0 && nsupermodels > 1) {
                tl_error("Flux coupling not yet supported.");
            }
            
        } while (fluxcouplingconverged == NO);
        
        // ---------------------------------------------------
        // supermodel loop -----------------------------------
        for (isuper_model=0; isuper_model<nsupermodels; isuper_model++) {
            sm = &(superModel[isuper_model]);
            
            if (model_finish[isuper_model] == 1) continue;
            
            for (isubmodel=0; isubmodel<sm->nsubmodels; isubmodel++){
                mod = &(sm->submodel[isubmodel]);
                if(mod->proc_flag == 1){
                    if (mod->series_out->entry[mod->ientry_out].time == mod->t_prev) {
                        mod->ientry_out++;
                    }
                }
            }
            
            for (isubmodel=0; isubmodel<sm->nsubmodels; isubmodel++){
                if(sm->submodel[isubmodel].proc_flag==1){
                    mod = &(sm->submodel[isubmodel]);
                    while ((mod->series_out->size > mod->ientry_out) &&
                           (mod->series_out->entry[mod->ientry_out].time <= mod->t_prev + mod->dt) &&
                           (mod->series_out->entry[mod->ientry_out].time <= mod->t_final)) {
                        
                        if (fabs(mod->t_prev + mod->dt - mod->t_final) < SMALL8) {
                            smodel_print_ts(mod, mod->t_final);
                            while (mod->series_out->entry[mod->ientry_out].time <= mod->t_prev + mod->dt &&
                                   mod->series_out->size > mod->ientry_out) {
                                mod->ientry_out++;
                            }
                        } else {
                            smodel_print_ts(mod, mod->t_prev + mod->dt);
                            mod->ientry_out++;
                            while (mod->series_out->entry[mod->ientry_out].time <= mod->t_prev + mod->dt) {
                                mod->ientry_out++;
                            }
                        }
                        
                    }
                }
            }
            
            // Adapt mesh if required - Loop over all submodels in the current supermodel*/
            for (isubmodel=0; isubmodel<sm->nsubmodels; isubmodel++){
                if(sm->submodel[isubmodel].proc_flag==1){
                    mod = &(sm->submodel[isubmodel]);
                    if (mod->flag.GRID_ADAPTION) {
                        adpt_main(mod);
                        if(mod->flag.ADAPTED_THE_GRID == YES) was_grid_adapted = 1;
                    }
                }
            }
            
#ifdef _MESSG
            ierr = MPI_Allreduce(&was_grid_adapted,&was_grid_adapted_max,1,MPI_INT,MPI_MAX,sm->supersmpi->ADH_COMM);
            was_grid_adapted=was_grid_adapted_max;
#endif
            
            if (was_grid_adapted==1){
                ssuperModel_forward_map(sm);
                if(sm->supersmpi_wvel!=NULL)ssuperModel_forward_map_wvel(sm);
            }
            was_grid_adapted=0;
            
            // Add this CSTORM variable here in case of any reallocation (such as adaption), pointer must be re-defined.
            if(coupled_normal_flux != NULL){
                for (isubmodel=0; isubmodel<sm->nsubmodels; isubmodel++){
                    mod = &(sm->submodel[isubmodel]);
                    if(mod->proc_flag==1){
                        calculate_supermodel_fluxes(isuper_model, isubmodel, mod);
#ifdef _MESSG
                        gather_fluxes(mod, isuper_model, isubmodel);
#endif
                    }
#ifdef _MESSG
                    if(mod->grid->smpi->myid==0) {
                        if(mod->grid->ndim == 2) {
                            ie=mod->grid->macro_nelems1d;
                        }else{
                            ie=mod->grid->macro_nelems2d;
                        }
                    }else{
                        ie=0;
                    }
                    MPI_Bcast(&ie, 1, MPI_INT,root_ids[isuper_model][isubmodel],cstorm_comm);
                    MPI_Bcast(&coupled_normal_flux[isuper_model][isubmodel],ie,MPI_DOUBLE,root_ids[isuper_model][isubmodel], cstorm_comm);
#endif
                }
            }
            mod_cstorm = &(sm[0].submodel[0]);
            
#ifdef _DEBUG
            if (DEBUG) {
#ifdef _MESSG
                tag(cstorm_comm);
#else
                tag();
#endif
            }
#endif
            
            // test case write
            if (test_case_flag.on) {
                write_testcase_error(sm);
            }
            
#ifdef _DEBUG
            if (DEBUG) {
#ifdef _MESSG
                tag(cstorm_comm);
                messg_barrier(cstorm_comm);
#else
                tag();
#endif
            }
#endif
            
            // is this model finished?
            int tc_end_mod_count = 0;
            for (isubmodel=0; isubmodel<superModel->nsubmodels; isubmodel++){
                mod = &(sm->submodel[isubmodel]);
                if(mod->proc_flag==1){
                    if (tc_end(mod) == YES) tc_end_mod_count++;
                    
                } else {
                    mod->dt *= DT_ENLARGE_FACTOR; // enlarge DT. It will be minimized if necessary in update_dt
                    tc_end_mod_count++;
                }
            }
            
#ifdef _DEBUG
            if (DEBUG) {
#ifdef _MESSG
                tag(cstorm_comm);
                messg_barrier(cstorm_comm);
#else
                tag();
#endif
            }
#endif
            
            if (tc_end_mod_count == sm->nsubmodels){
                model_finish[isuper_model] = 1; // Meaning this supermodel has reached its final time
            }
            else {
                update_dt(sm);
            }
            sum_finished_models = sarray_sum_int(model_finish,nsupermodels);
            
        } // End super model loop
        
    } while (sum_finished_models != nsupermodels);
    
    
    // supermodel loop -----------------------------------
    // ---------------------------------------------------
    
    // cjt :: put back stored values in case CSTORM is applied
    for (isuper_model=0; isuper_model<nsupermodels; isuper_model++) {
        mod = &(sm[isuper_model].submodel[0]); // non-fully coupled models only
        mod->dt = dt_store[isuper_model];
        mod->t_final = mod->t_final_store;
    }
    
    
    /* --------------------------------------------------- */
    /* --------------------------------------------------- */
    
    /* Post Process Chop */
    for (isuper_model=0; isuper_model<nsupermodels; isuper_model++) {
        sm = &(superModel[isuper_model]);
        for (isubmodel=0; isubmodel<sm->nsubmodels; isubmodel++){
            mod = &(sm->submodel[isubmodel]);
            if (mod->file_output.chop == ON){
                if (mod->grid->smpi->myid <= 0){
                    tl_3d_to_2d();
                }
            }
        }
    }
    
    return ierr;
}

// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------
// Update time-step for one superModel
void update_dt(SSUPER_MODEL *sm){
    SMODEL *mod;
    int i;
    double min_dt = UNSET_FLT; /* Not required, but purposely using min_dt to keep track */
    double new_dt;
    
    //fflush(stdout); for (i=0;i<100000;i++);messg_barrier(cstorm_comm); tag();
    
    for (i=0; i<sm->nsubmodels; i++){
#ifdef _MESSG
        if(sm->submodel[i].proc_flag==1){
#endif
            mod = &(sm->submodel[i]);
            if ((min_dt > mod->dt)||(min_dt < 0.0)){ min_dt = mod->dt; }
#ifdef _MESSG
        }
#endif
    }
    
#ifdef _MESSG
    new_dt = messg_dmin(min_dt, sm->supersmpi->ADH_COMM);
    sm->dt = new_dt;
#else
    new_dt = min_dt;
    sm->dt = min_dt;
#endif
    
    for (i=0; i<sm->nsubmodels; i++){
#ifdef _MESSG
        if(sm->submodel[i].proc_flag==1){
#endif
            mod = &(sm->submodel[i]);
            mod->dt = new_dt;
#ifdef _MESSG
        }
#endif
        // printf("\tsubmodel[%i].dt = %.6e\n", i, mod->dt);
    }
    
    //fflush(stdout); for (i=0;i<100000;i++);messg_barrier(cstorm_comm); tag();
}

// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------

/* Flux array messaging functions  */

void set_flux_messaging(SMODEL *mod){
    int i,j;
    int ierr;
#ifdef _MESSG
    MPI_Status *msg_status;
    ndata = (int **) tl_alloc(sizeof(int *), mod->grid->smpi->npes);
    my_nelem_ext = (int *) tl_alloc(sizeof(int), mod->grid->smpi->npes);
    
    
    if(mod->grid->ndim==2){
        my_nelem = mod->grid->nelems1d;
        ierr = MPI_Reduce(&my_nelem, &my_nelem_max, 1, MPI_INT, MPI_MAX, 0, mod->grid->smpi->ADH_COMM);
        ierr = MPI_Allgather(&my_nelem, 1, MPI_INT, my_nelem_ext, 1, MPI_INT, mod->grid->smpi->ADH_COMM);
        if(mod->grid->smpi->myid != 0) {
            my_nelem_max=my_nelem;
        }
        if (mod->grid->smpi->myid <= 0){
            for (i=0;i<mod->grid->smpi->npes;i++){
                ndata[i] = (int *) tl_alloc(sizeof(int), my_nelem_max);
            }
            for (i=1;i<mod->grid->smpi->npes;i++){
                ierr = MPI_Recv(ndata[i], my_nelem_max, MPI_INT, i, 998, mod->grid->smpi->ADH_COMM, msg_status);
            }
        }else{
            ndata[0]=(int *) tl_alloc(sizeof(int), my_nelem);
            for(i=0;i<mod->grid->nelems1d;i++){
                ndata[0][i] = mod->grid->elem1d[i].gid;
            }
            ierr = MPI_Send(ndata[0], my_nelem, MPI_INT, 0, 998, mod->grid->smpi->ADH_COMM);
        }
    }else{
        my_nelem = mod->grid->nelems2d;
        ierr = MPI_Reduce(&my_nelem, &my_nelem_max, 1, MPI_INT, MPI_MAX, 0, mod->grid->smpi->ADH_COMM);
        ierr = MPI_Allgather(&my_nelem, 1, MPI_INT, my_nelem_ext, 1, MPI_INT, mod->grid->smpi->ADH_COMM);
        if(mod->grid->smpi->myid != 0) {
            my_nelem_max=my_nelem;
        }
        if (mod->grid->smpi->myid <= 0){
            for (i=0;i<mod->grid->smpi->npes;i++){
                ndata[i] = (int *) tl_alloc(sizeof(int), my_nelem_max);
            }
            for (i=1;i<mod->grid->smpi->npes;i++){
                ierr = MPI_Recv(ndata[i], my_nelem_max, MPI_INT, i, 998, mod->grid->smpi->ADH_COMM, msg_status);
            }
        }else{
            ndata[0]=(int *) tl_alloc(sizeof(int), my_nelem);
            for(i=0;i<mod->grid->nelems2d;i++){
                ndata[0][i] = mod->grid->elem2d[i].gid;
            }
            ierr = MPI_Send(ndata[0], my_nelem, MPI_INT, 0, 998, mod->grid->smpi->ADH_COMM);
        }
    }
#endif
    return;
}

// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------

void gather_fluxes(SMODEL *mod, int isup, int isub){
    double *gdata, *edata;
    int i, ip;
    int ierr;
#ifdef _MESSG
    MPI_Status *msg_status;
    edata = (double *) tl_alloc(sizeof(double), my_nelem_max);
    if (mod->grid->smpi->myid <= 0) {
        if(mod->grid->ndim ==2){
            gdata = (double *) tl_alloc(sizeof(double), mod->grid->macro_nelems1d);
            for (i=0;i<mod->grid->nelems1d;i++){
                gdata[mod->grid->elem1d[i].gid] = coupled_normal_flux[isup][isub][mod->grid->elem1d[i].gid];
            }
            for (ip = 1; ip < mod->grid->smpi->npes; ip++) {
                ierr = MPI_Recv(edata, my_nelem_ext[ip], MPI_DOUBLE, ip, 999, mod->grid->smpi->ADH_COMM, msg_status);
                for (i = 0; i < my_nelem_ext[ip]; i++){
                    gdata[ndata[ip][i]] = edata[i];
                }
            }
            for (i=0;i<mod->grid->macro_nelems1d;i++){
                coupled_normal_flux[isup][isub][i] = gdata[i];
            }
        }else{
            gdata = (double *) tl_alloc(sizeof(double), mod->grid->macro_nelems2d);
            for (i=0;i<mod->grid->nelems2d;i++){
                gdata[mod->grid->elem1d[i].gid] = coupled_normal_flux[isup][isub][mod->grid->elem2d[i].gid];
            }
            for (ip = 1; ip < mod->grid->smpi->npes; ip++) {
                ierr = MPI_Recv(edata, my_nelem_ext[ip], MPI_DOUBLE, ip, 999, mod->grid->smpi->ADH_COMM, msg_status);
                for (i = 0; i < my_nelem_ext[ip]; i++){
                    gdata[ndata[ip][i]] = edata[i];
                }
            }
            for (i=0;i<mod->grid->macro_nelems2d;i++){
                coupled_normal_flux[isup][isub][i] = gdata[i];
            }
        }
    }else{
        if(mod->grid->ndim ==2){
            for (i=0;i<mod->grid->nelems1d;i++){
                edata[mod->grid->elem1d[i].gid] = coupled_normal_flux[isup][isub][mod->grid->elem1d[i].gid];
            }
            ierr = MPI_Send(edata, my_nelem_ext[mod->grid->smpi->myid], MPI_DOUBLE, 0, 999, mod->grid->smpi->ADH_COMM);
        }else{
            for (i=0;i<mod->grid->nelems2d;i++){
                edata[mod->grid->elem2d[i].gid] = coupled_normal_flux[isup][isub][mod->grid->elem2d[i].gid];
            }
            ierr = MPI_Send(edata, my_nelem_ext[mod->grid->smpi->myid], MPI_DOUBLE, 0, 999, mod->grid->smpi->ADH_COMM);
        }
    }
#endif
    return;
}

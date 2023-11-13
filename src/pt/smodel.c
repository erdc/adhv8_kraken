#include "global_header.h"

void smodel_read_opt(SMODEL * mod);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// constructor for user input file read
int smodel_init(SMODEL *mod) {
    int ip;
    char *line = NULL, card[100], card2[100];
    size_t len = 0;
    ssize_t readfile = -1;
    
    RUN_VERBOSE = TRUE;
    
    int npes=1,myid=0;
#ifdef _MPI
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    mod->screen_print_timestep = true;
    
    mod->np = 0;
    mod->fp_dpl = NULL;
    mod->fp_vel = NULL;
    mod->fp_wq = NULL;
    mod->fp_salinity = NULL;
    mod->fp_doxygen  = NULL;
    mod->fp_exposure = NULL;
    mod->fp_out = NULL;
    mod->fp_error = NULL;
    mod->get_analytic_p = NULL;
    
    mod->dpl = NULL;
    mod->dpl_tL = NULL;
    mod->dpl_tR = NULL;
    
    mod->oxygen = NULL;
    mod->oxygen_tL = NULL;
    mod->oxygen_tR = NULL;
    
    mod->salinity = NULL;
    mod->salinity_tL = NULL;
    mod->salinity_tR = NULL;
    
    mod->sunlight = NULL;
    mod->sunlight_tL = NULL;
    mod->sunlight_tR = NULL;
    
    // default time integrator is RK2
    mod->time_update_flag = 1;
    
    mod->dt = 0.;
    mod->t0 = 0.;
    mod->tf = 0.;
    mod->time = 0.;
    mod->tL = 0.;
    mod->tR = 0.;
    mod->tL_dpl = 0.;
    mod->tR_dpl = 0.;
    mod->tL_wq = 0.;
    mod->tR_wq = 0.;
    mod->max_t_abs_error_mag = 0.;
    mod->max_t_rel_error_mag = 0.;

    mod->normalize = false;
    mod->grid = NULL;
    mod->p = NULL;

    // allocate and read grid
    mod->grid = (SGRID *)malloc(sizeof(SGRID));
    sgrid_read_adh(&mod->grid,mod->file_base);

    // Normalize, if requested in the model options.
    smodel_read_opt(mod);
    sgrid_normalize(mod->normalize, mod->grid);
    
    // Read the rest of the model input file
    smodel_read_input(mod);
    mod->time = mod->t0;
    assert(mod->np > 0);

    // allocate particle initial position arrays and displacements
    mod->vel =    (SVECT *)malloc(sizeof(SVECT) * mod->grid->nnodes);
    mod->vel_tR = (SVECT *)malloc(sizeof(SVECT) * mod->grid->nnodes);
    mod->vel_tL = (SVECT *)malloc(sizeof(SVECT) * mod->grid->nnodes);
    svect_init_array(mod->pdpl, mod->np); // particle displacements, not grid
    svect_init_array(mod->vel,    mod->grid->nnodes);
    svect_init_array(mod->vel_tL, mod->grid->nnodes);
    svect_init_array(mod->vel_tR, mod->grid->nnodes);
    
    // default timeunits are seconds  // 0(s), 1(m), 2(h), 3(d)
    mod->time_units_dep = 1;
    mod->time_units_vel = 1;
    mod->time_units_wq = 1;
    mod->time_units_salinity = 1;
    mod->time_units_deoxygen = 1;
    mod->time_units_expossure = 1;
    
    // open velocity file
    strcpy(mod->vel_filename,mod->file_base);
    strcat(mod->vel_filename, "_vel.dat");
    if (RUN_VERBOSE) printf("OPENING VELOCITY FILE: %s & READING FIRST 2 SNAPS \n",mod->vel_filename);
    mod->fp_vel = fopen(mod->vel_filename,"r");
    //mod->tL = read_adh_velocity_frame(mod->fp_vel, mod->grid, mod->vel_tL);
    //mod->tR = read_adh_velocity_frame(mod->fp_vel, mod->grid, mod->vel_tR);
    //svect_copy_array(mod->vel, mod->vel_tL, mod->grid->nnodes); // assume for now that the start of the firt snap is simulation t0
    //if (RUN_VERBOSE) printf(":: tL: %f tR: %f \n",mod->tL,mod->tR);
    // get units
    while(1) {
        readfile = getline(&line, &len, mod->fp_vel); // TS 0 time
        if (readfile == -1) {
            printf("ERROR :: velocity file end reached without finding timeunits!\n");
            exit(-1);
        }
        //printf("%s\b",line);
        sscanf(line,"%s %s",card,card2);
        if (strcmp(card, "TIMEUNITS") == 0) {
            if (strcmp(card2, "SECONDS") == 0) mod->time_units_vel = 1;
            if (strcmp(card2, "MINUTES") == 0) mod->time_units_vel = 60;
            if (strcmp(card2, "HOURS") == 0)   mod->time_units_vel = 3600;
            if (strcmp(card2, "DAYS") == 0)    mod->time_units_vel = 86400;
            break;
        }
    }
    if (RUN_VERBOSE) printf("FINISHED OPENING VELOCITY FILE: %s :: timeunits: %s :: time_units_vel: %d\n",mod->vel_filename,card2,mod->time_units_vel);
    
    // open grid displacement file
    strcpy(mod->dpl_filename,mod->file_base);
    strcat(mod->dpl_filename, "_dpl.dat");
    mod->fp_dpl = fopen(mod->dpl_filename,"r");
    if (mod->fp_dpl != NULL) {
        if (RUN_VERBOSE) printf("OPENING GRID DISPLACEMENT FILE: %s & READING FIRST 2 SNAPS ",mod->dpl_filename);
        mod->dpl =      (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->dpl_tR =   (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->dpl_tL =   (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->tL_dpl = read_adh_dpl_frame(mod->fp_dpl, mod->grid, mod->dpl_tL, mod->time_units_dep);
        mod->tR_dpl = read_adh_dpl_frame(mod->fp_dpl, mod->grid, mod->dpl_tR, mod->time_units_dep);
//        while (mod->time + mod->dt - 1e-6 > mod->tR_dpl) {
//            if (mod->tR_dpl < mod->time) {
//                mod->tL_dpl = mod->tR_dpl;
//                sarray_copy_dbl(mod->dpl_tL,mod->dpl_tR,mod->grid->nnodes);
//            }
//            mod->tR_dpl = read_adh_dpl_frame(mod->fp_dpl, mod->grid, mod->dpl_tR);
//            if (RUN_VERBOSE) printf("** updated grid displacement snap :: tL: %f tR: %f ** \n",mod->tL_dpl,mod->tR_dpl);
//        }
//        if (RUN_VERBOSE) printf(":: tL: %f tR: %f \n",mod->tL_dpl,mod->tR_dpl);
        if (RUN_VERBOSE) printf("FINISHED OPENING GRID DISPLACEMENT FILE: %s\n",mod->dpl_filename);
    } else {
        mod->dpl =      NULL;
        mod->dpl_tR =   NULL;
        mod->dpl_tL =   NULL;
    }
    
    // open water quality file
    strcpy(mod->wq_filename,mod->file_base);
    strcat(mod->wq_filename, "_wq.dat");
    mod->fp_wq = fopen(mod->wq_filename,"r");
    if (mod->fp_wq != NULL) {
        if (RUN_VERBOSE) printf("OPENING WATER QUALITY FILE FILE: %s & READING FIRST 2 SNAPS ",mod->wq_filename);
        mod->oxygen =      (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->oxygen_tR =   (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->oxygen_tL =   (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->salinity =    (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->salinity_tR = (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->salinity_tL = (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->sunlight =    (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->sunlight_tR = (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->sunlight_tL = (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->tL_wq = read_adh_wq_frame(mod->fp_wq, mod->grid, mod->oxygen_tL, mod->salinity_tL, mod->sunlight_tL, mod->time_units_wq);
        mod->tR_wq = read_adh_wq_frame(mod->fp_wq, mod->grid, mod->oxygen_tR, mod->salinity_tR, mod->sunlight_tR, mod->time_units_wq);
        sarray_copy_dbl(mod->oxygen, mod->oxygen_tL, mod->grid->nnodes);     // assume for now that the start of the first snap is simulation t0
        sarray_copy_dbl(mod->salinity, mod->salinity_tL, mod->grid->nnodes); // assume for now that the start of the first snap is simulation t0
        sarray_copy_dbl(mod->sunlight, mod->sunlight_tL, mod->grid->nnodes); // assume for now that the start of the first snap is simulation t0
        if (RUN_VERBOSE) printf(":: tL: %f tR: %f \n",mod->tL_wq,mod->tR_wq);
        if (RUN_VERBOSE) printf("FINISHED OPENING WATER QUALITY FILE: %s\n",mod->wq_filename);
    } else {
        mod->oxygen =      NULL;
        mod->oxygen_tR =   NULL;
        mod->oxygen_tL =   NULL;
        mod->salinity =    NULL;
        mod->salinity_tR = NULL;
        mod->salinity_tL = NULL;
        mod->sunlight =    NULL;
        mod->sunlight_tR = NULL;
        mod->sunlight_tL = NULL;
    }
    
    // open output file
    char scatter_file[100];
    strcpy(scatter_file,mod->file_base); strcat(scatter_file, "_scatter_out.csv");
    if (RUN_VERBOSE) printf("OPENING PARTICLE OUTPUT FILE: %s\n",scatter_file);
    mod->fp_out = fopen(scatter_file,"w");
    fprintf(mod->fp_out,"particle ID, t, x, y, z\n");
    if (RUN_VERBOSE) printf("FINISHED OPENING PARTICLE OUTPUT FILE: scatter.csv\n");
    
    // open error file output if used
    if (mod->get_analytic_p != NULL) {
        char error_file[100];
        strcpy(error_file,mod->file_base); strcat(error_file, "_error_out.csv");
        mod->fp_error = fopen(error_file,"w");
        fprintf(mod->fp_error,"particle id, t, x, y, z, xa, ya, za, abs_error, rel_error, max_t_abs, max_t_rel\n");
    }
    
#ifdef _ADH_HDF5
    if (RUN_VERBOSE) printf("Initializing XDMF grid output files\n"); fflush(stdout);
    if (myid ==0) xdmf_init(mod->grid, 1 /*npes*/, 0 /*myid*/,mod->file_base); // all pe's have same grid
    if (RUN_VERBOSE) printf("Initializing XDMF particle output files\n");
    xdmf_init_particles(&mod->hdf5,mod->np,mod->p,npes,myid,mod->file_base);
#endif
    
    return 0;
}

// constructor for input variable call
int smodel_init_noRead(SMODEL *mod, SGRID *grid, char *file_base, double t0, double tf, double dt, int np, int *node_ID, int *elem_ID, int spool, int ptype) {
    int ip;
    int npes=1,myid=0;
    ssize_t readfile = -1; // do still need to read in time units (seconds/minutes/etc.) from velocity file
    char *line=NULL, card[100], card2[100];
    size_t len = 0;
#ifdef _MPI
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
    
    mod->screen_print_timestep = true;
    
    mod->np = np;
    mod->fp_vel = NULL;
    mod->fp_wq = NULL;
    mod->fp_salinity = NULL;
    mod->fp_doxygen  = NULL;
    mod->fp_exposure = NULL;
    mod->fp_out = NULL;
    mod->fp_error = NULL;
    
    mod->dt = dt;
    mod->t0 = t0;
    mod->tf = tf;
    mod->time = t0;
    mod->spool = spool;
    mod->tL = 0.;
    mod->tR = 0.;
    mod->tL_wq = 0.;
    mod->tR_wq = 0.;
    mod->max_t_abs_error_mag = 0.;
    mod->max_t_rel_error_mag = 0.;
    
    // default time integrator is RK2
    mod->time_update_flag = 1;
    
    mod->normalize = false;
    mod->grid = grid;
    mod->p = NULL;
    
    // allocate particle initial position arrays and displacements
    mod->p0 =     (SVECT *)malloc(sizeof(SVECT) * mod->np);
    mod->pdpl =   (SVECT *)malloc(sizeof(SVECT) * mod->np);
    mod->vel =    (SVECT *)malloc(sizeof(SVECT) * mod->grid->nnodes);
    mod->vel_tR = (SVECT *)malloc(sizeof(SVECT) * mod->grid->nnodes);
    mod->vel_tL = (SVECT *)malloc(sizeof(SVECT) * mod->grid->nnodes);
    svect_init_array(mod->pdpl,   mod->np);
    svect_init_array(mod->vel,    mod->grid->nnodes);
    svect_init_array(mod->vel_tL, mod->grid->nnodes);
    svect_init_array(mod->vel_tR, mod->grid->nnodes);
    
    // default timeunits are seconds  // 0(s), 1(m), 2(h), 3(d)
    mod->time_units_dep = 1;
    mod->time_units_vel = 1;
    mod->time_units_wq = 1;
    mod->time_units_salinity = 1;
    mod->time_units_deoxygen = 1;
    mod->time_units_expossure = 1;

    // allocate and initialize particles
    mod->p = (SPARTICLE *)malloc(sizeof(SPARTICLE) * mod->np);
    sparticle_init(mod->p,mod->np);
    for (ip=0;ip<mod->np;ip++) {
        mod->p[ip].ptype = ptype;
        mod->p[ip].isActive = YES;
        mod->p[ip].nodeID = node_ID[ip];
        mod->p[ip].elemID = elem_ID[ip];
        mod->p[ip].r.x = mod->grid->node[node_ID[ip]].x;
        mod->p[ip].r.y = mod->grid->node[node_ID[ip]].y;
        mod->p[ip].r.z = mod->grid->node[node_ID[ip]].z;
        mod->p0[ip] = mod->p[ip].r;
    }
    
    // open velocity file
    strcpy(mod->vel_filename,mod->file_base);
    strcat(mod->vel_filename, "_vel.dat");
    if (RUN_VERBOSE) printf("OPENING VELOCITY FILE: %s & READING FIRST 2 SNAPS ",mod->vel_filename);
    mod->fp_vel = fopen(mod->vel_filename,"r");
//    mod->tL = read_adh_velocity_frame(mod->fp_vel, mod->grid, mod->vel_tL);
//    mod->tR = read_adh_velocity_frame(mod->fp_vel, mod->grid, mod->vel_tR);
//    if (RUN_VERBOSE) printf(":: tL: %f tR: %f \n",mod->tL,mod->tR);
    mod->time_units_vel = 1;
    while(1) {
        readfile = getline(&line, &len, mod->fp_vel); // TS 0 time
        if (readfile == -1) {
            printf("ERROR :: velocity file end reached without finding timeunits!\n");
            exit(-1);
        }
        //printf("%s\b",line);
        sscanf(line,"%s %s",card,card2);
        if (strcmp(card, "TIMEUNITS") == 0) {
            if (strcmp(card2, "SECONDS") == 0) mod->time_units_vel = 1;
            if (strcmp(card2, "MINUTES") == 0) mod->time_units_vel = 60;
            if (strcmp(card2, "HOURS") == 0)   mod->time_units_vel = 3600;
            if (strcmp(card2, "DAYS") == 0)    mod->time_units_vel = 86400;
            break;
        }
    }
    if (RUN_VERBOSE) printf("FINISHED OPENING VELOCITY FILE: %s :: timeunits: %s :: time_units_vel: %d\n",mod->vel_filename,card2,mod->time_units_vel);
    
    // open grid displacement file
    strcpy(mod->dpl_filename,mod->file_base);
    strcat(mod->dpl_filename, "_dpl.dat");
    mod->fp_dpl = fopen(mod->dpl_filename,"r");
    if (mod->fp_dpl != NULL) {
        if (RUN_VERBOSE) printf("OPENING GRID DISPLACEMENT FILE: %s & READING FIRST 2 SNAPS ",mod->dpl_filename);
        mod->dpl =      (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->dpl_tR =   (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->dpl_tL =   (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->tL_dpl = read_adh_dpl_frame(mod->fp_dpl, mod->grid, mod->dpl_tL, mod->time_units_dep);
        mod->tR_dpl = read_adh_dpl_frame(mod->fp_dpl, mod->grid, mod->dpl_tR, mod->time_units_dep);
//        while (mod->time + mod->dt - 1e-6 > mod->tR_dpl) {
//            if (mod->tR_dpl < mod->time) {
//                mod->tL_dpl = mod->tR_dpl;
//                sarray_copy_dbl(mod->dpl_tL,mod->dpl_tR,mod->grid->nnodes);
//            }
//            mod->tR_dpl = read_adh_dpl_frame(mod->fp_dpl, mod->grid, mod->dpl_tR);
//            if (RUN_VERBOSE) printf("** updated grid displacement snap :: tL: %f tR: %f ** \n",mod->tL_dpl,mod->tR_dpl);
//        }
//        if (RUN_VERBOSE) printf(":: tL: %f tR: %f \n",mod->tL_dpl,mod->tR_dpl);
        if (RUN_VERBOSE) printf("FINISHED OPENING GRID DISPLACEMENT FILE: %s\n",mod->dpl_filename);
    } else {
        mod->dpl =      NULL;
        mod->dpl_tR =   NULL;
        mod->dpl_tL =   NULL;
    }
    
    // open water quality file
    if (mod->get_analytic_wq != NULL) {
        strcpy(mod->wq_filename,mod->file_base);
        strcat(mod->wq_filename, "_wq.dat");
        mod->fp_wq = fopen(mod->wq_filename,"r");
        mod->oxygen =      (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->oxygen_tR =   (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->oxygen_tL =   (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->salinity =    (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->salinity_tR = (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->salinity_tL = (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->sunlight =    (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->sunlight_tR = (double *)calloc(mod->grid->nnodes,sizeof(double));
        mod->sunlight_tL = (double *)calloc(mod->grid->nnodes,sizeof(double));
        if (RUN_VERBOSE) printf("OPENING WATER QUALITY FILE FILE: %s & READING FIRST 2 SNAPS ",mod->wq_filename);
        mod->tL_wq = read_adh_wq_frame(mod->fp_wq, mod->grid, mod->oxygen_tL, mod->salinity_tL, mod->sunlight_tL, mod->time_units_wq);
        mod->tR_wq = read_adh_wq_frame(mod->fp_wq, mod->grid, mod->oxygen_tR, mod->salinity_tR, mod->sunlight_tR, mod->time_units_wq);
        sarray_copy_dbl(mod->oxygen, mod->oxygen_tL, mod->grid->nnodes);     // assume for now that the start of the first snap is simulation t0
        sarray_copy_dbl(mod->salinity, mod->salinity_tL, mod->grid->nnodes); // assume for now that the start of the first snap is simulation t0
        sarray_copy_dbl(mod->sunlight, mod->sunlight_tL, mod->grid->nnodes); // assume for now that the start of the first snap is simulation t0
        if (RUN_VERBOSE) printf(":: tL: %f tR: %f \n",mod->tL_wq,mod->tR_wq);
        if (RUN_VERBOSE) printf("FINISHED OPENING VELOCITY FILE: %s\n",mod->wq_filename);
    } else {
        mod->oxygen =      NULL;
        mod->oxygen_tR =   NULL;
        mod->oxygen_tL =   NULL;
        mod->salinity =    NULL;
        mod->salinity_tR = NULL;
        mod->salinity_tL = NULL;
        mod->sunlight =    NULL;
        mod->sunlight_tR = NULL;
        mod->sunlight_tL = NULL;
    }
    
    // open output file
    if (RUN_VERBOSE) printf("OPENING PARTICLE OUTPUT FILE: scatter.csv\n");
    mod->fp_out = fopen("scatter.csv","w");
    fprintf(mod->fp_out,"particle ID, t, x, y, z\n");
    if (RUN_VERBOSE) printf("FINISHED OPENING PARTICLE OUTPUT FILE: scatter.csv\n");
    
    // open error file output if used
    if (mod->get_analytic_p != NULL) {
        char error_file[100];
        strcpy(error_file,mod->file_base); strcat(error_file, "_error_out.csv");
        mod->fp_error = fopen(error_file,"w");
        fprintf(mod->fp_error,"particle id, t, x, y, z, xa, ya, za, abs_error, rel_error, max_t_abs, max_t_rel\n");
    }
    
#ifdef _ADH_HDF5
    if (RUN_VERBOSE) printf("Initializing XDMF grid output files\n");
    if (myid == 0) xdmf_init(mod->grid, 1 /*npes*/, 0 /*myid*/, mod->file_base);  // all pe's have same grid
    if (RUN_VERBOSE) printf("Initializing XDMF particle output files\n");
    xdmf_init_particles(&mod->hdf5,mod->np,mod->p, 1 /*npes*/, 0 /*myid*/, mod->file_base);
#endif
    
    return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int smodel_finalize(SMODEL *mod) {
    
    int npes=1,myid=0;
#ifdef _MPI
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
    
#ifdef _ADH_HDF5
    // CJT :: note - not set up for MPI yet.  Need a particle finalize.
    //if (RUN_VERBOSE) printf("Finalizing XDMF grid output files\n");
    //xdmf_finalize(&(mod->grid->hdf5), mod->file_base, 1 /*npes*/, 0 /*myid*/, 0 /*flag*/); // all PEs have same grid
    if (RUN_VERBOSE) printf("Finalizing XDMF particle output files\n");
    xdmf_particle_finalize(&(mod->hdf5), mod->file_base,npes,myid, 0 /*flag*/);
#endif
    
    if (RUN_VERBOSE) printf("Freeing model arrays\n");
    if (mod->p0 != NULL)     free(mod->p0);
    if (mod->pdpl != NULL)   free(mod->pdpl);
    if (mod->vel != NULL)    free(mod->vel);
    if (mod->vel_tL != NULL) free(mod->vel_tL);
    if (mod->vel_tR != NULL) free(mod->vel_tR);
    if (mod->dpl    != NULL) free(mod->dpl);
    if (mod->dpl_tL != NULL) free(mod->dpl_tL);
    if (mod->dpl_tR != NULL) free(mod->dpl_tR);
    
    if (RUN_VERBOSE) printf("Freeing water quality arrays\n");
    if (mod->oxygen    != NULL) free(mod->oxygen);
    if (mod->oxygen_tL != NULL) free(mod->oxygen_tL);
    if (mod->oxygen_tR != NULL) free(mod->oxygen_tR);
    if (mod->salinity    != NULL) free(mod->salinity);
    if (mod->salinity_tL != NULL) free(mod->salinity_tL);
    if (mod->salinity_tR != NULL) free(mod->salinity_tR);
    if (mod->sunlight    != NULL) free(mod->sunlight);
    if (mod->sunlight_tL != NULL) free(mod->sunlight_tL);
    if (mod->sunlight_tR != NULL) free(mod->sunlight_tR);
    
    // deallocate grid
    if (RUN_VERBOSE) printf("Freeing model grid\n");
    if (mod->grid != NULL) sgrid_free(mod->grid);
    
    // close data file for particle trajectory output
    if (RUN_VERBOSE) printf("Freeing input and output files\n");
    if (mod->fp_out != NULL) fclose(mod->fp_out);
    if (mod->fp_vel != NULL) fclose(mod->fp_vel);
    if (mod->fp_wq != NULL) fclose(mod->fp_wq);
    if (mod->fp_salinity != NULL) fclose(mod->fp_salinity);
    if (mod->fp_doxygen  != NULL) fclose(mod->fp_doxygen);
    if (mod->fp_exposure != NULL) fclose(mod->fp_exposure);
    if (mod->fp_error != NULL)    fclose(mod->fp_error);
    
    
    free(mod);
    
    return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int smodel_run_lockstep(SMODEL *mod) {
    
    int ip, inode, flag, it=0;
    SVECT p_new;
    
    // calculate errors and print initial conditions
    //smodel_print_ts(mod);
    
    // time-loop
    while(1) {
        it++;
        
        if (mod->time > mod->tf + 1e-6) break;

        // update velocity time-snap extents
        if (mod->time < mod->tf) { // do not read on the last time-step
            bool break_time_loop = false;
            while(1)  {
                if (mod->time + mod->dt - 1e-6 < mod->tR) break;
                if (mod->tR - 1e-6  < mod->time + mod->dt) {
                    mod->tL = mod->tR;
                    svect_copy_array(mod->vel_tL, mod->vel_tR, mod->grid->nnodes);
                }
                mod->tR = read_adh_velocity_frame(mod->fp_vel, mod->grid, mod->vel_tR, mod->time_units_vel, mod->normalize);
                if (mod->tR == -1) {
                    printf("** in smodel_run_lockstep :: reached end of velocity file before last time step. Stopping. \n");
                    break_time_loop = true;
                }
                if (RUN_VERBOSE) printf("** smodel_run_lockstep :: updated velocity snap :: t: %f t+dt: %f || tL: %f tR: %f ** \n",mod->time,mod->time+mod->dt,mod->tL,mod->tR);
            }
            if (break_time_loop) break;
            
            // update grid displacemment time-snap extents
            if (mod->fp_dpl != NULL) {
                while(1)  {
                    if (mod->time + mod->dt - 1e-6 < mod->tR_dpl) break;
                    if (mod->tR_dpl - 1e-6  < mod->time + mod->dt) {
                        mod->tL_dpl = mod->tR_dpl;
                        sarray_copy_dbl(mod->dpl_tL,mod->dpl_tR,mod->grid->nnodes);
                    }
                    mod->tR_dpl = read_adh_dpl_frame(mod->fp_dpl, mod->grid, mod->dpl_tR, mod->time_units_dep);
                    if (RUN_VERBOSE) printf("** smodel_run_lockstep :: updated grid displacement snap :: tL: %f tR: %f ** \n",mod->tL_dpl,mod->tR_dpl);
                }
            }
            
            // update water quality time-snap extents :: NOTE: snap size must be bigger than dt
            if (mod->fp_wq != NULL) {
                while(1)  {
                    if (mod->time + mod->dt - 1e-6 < mod->tR_wq) break;
                    if (mod->tR_wq - 1e-6  < mod->time + mod->dt) {
                        mod->tL_wq = mod->tR_wq;
                        sarray_copy_dbl(mod->oxygen_tL,  mod->oxygen_tR,  mod->grid->nnodes);
                        sarray_copy_dbl(mod->salinity_tL,mod->salinity_tR,mod->grid->nnodes);
                        sarray_copy_dbl(mod->sunlight_tL,mod->sunlight_tR,mod->grid->nnodes);
                        mod->tR_wq = read_adh_wq_frame(mod->fp_wq, mod->grid, mod->oxygen_tR, mod->salinity_tR, mod->sunlight_tR, mod->time_units_wq);
                        if (RUN_VERBOSE) printf("** smodel_run_lockstep :: updated water quality snap :: tL: %f tR: %f ** \n",mod->tL_wq,mod->tR_wq);
                    }
                }
            }
        }

        // calculate displacement for paraview simulation
        for (ip=0; ip<mod->np; ip++) mod->pdpl[ip] = svect_subtract(mod->p[ip].r,mod->p0[ip]);

        smodel_print_ts(mod,it);
        
        // update particle positions from hydrodynamic advection
        time_update(mod,2);
    }
    return 0;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int smodel_read_input(SMODEL *mod) {
    
    int i,node_np,nd,ie,npes=1,myid=0;
    int ptype;
    
#ifdef _MPI
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    // open model parameter file
    char pfile[100];
    strcpy(pfile,mod->file_base);
    strcat(pfile, ".in");
    if (myid == 0) {fflush(stdout); printf("OPENING & READING MODEL INPUT FILE: %s\n",pfile);}
    FILE *fp_poly;
    FILE *fp = fopen(pfile,"r");
    
    
    // read parameters
    size_t buffer_size = 80;
    char *buffer = malloc(buffer_size * sizeof(char));
    char delim[] = " ";  // line delimiter
    char *ptr = NULL;    // to read the card on line
    
    // read each line and print it to the screen
    mod->np_global = 0;
    mod->npoly = 0;
    while(getline(&buffer, &buffer_size, fp) != -1) {
        ptr = strtok(buffer, delim);
        //printf("%s ", ptr);
        if (strcmp(ptr,"DT") == 0) {
            ptr = strtok(NULL, delim);
            //printf("%s\n", ptr);
            mod->dt = atof(ptr);
        }
        if (strcmp(ptr,"T0") == 0) {
            ptr = strtok(NULL, delim);
            //printf("%s\n", ptr);
            mod->t0 = atof(ptr);
        }
        if (strcmp(ptr,"TF") == 0) {
            ptr = strtok(NULL, delim);
            //printf("%s\n", ptr);
            mod->tf = atof(ptr);
        }
        if (strcmp(ptr,"PNODE") == 0) {
            ptr = strtok(NULL, delim); // the total # of particles started at this node
            //printf("%s\n", ptr);
            mod->np_global += atoi(ptr);
        }
        if (strcmp(ptr,"PBARY") == 0) {
            ptr = strtok(NULL, delim); // the total # of particles started at the barycenter of this element
            //printf("%s\n", ptr);
            mod->np_global += atoi(ptr);
        }
        if (strcmp(ptr,"NWRITE") == 0) {
            ptr = strtok(NULL, delim);
            //printf("%s\n", ptr);
            mod->spool = atoi(ptr);
        }
        if (strcmp(ptr,"REEF") == 0) {
            mod->npoly++;
        }
        if (strcmp(ptr,"NORMALIZE") == 0) {
            ptr = strtok(NULL, delim);
            mod->normalize = (bool)atoi(ptr);
        }
    }
    rewind(fp);

    // enforce compatibility between dt, t0, & tf
    double nt = floor((mod->tf - mod->t0)/mod->dt);
    double new_tf = mod->t0 + nt*mod->dt;
    if (fabs(new_tf - mod->tf) > 1.e-6) {
        printf("NOTE :: smodel_read_input() :: setting tf = %f so that # of timesteps (:= (tf - t0)/dt) is an integer.\n", new_tf);
        mod->tf = new_tf;
    }
    
    // read polyon/reef data
    mod->poly = (SPOLYGON *)malloc(sizeof(SPOLYGON) * mod->npoly);
    mod->npoly = 0;
    while(getline(&buffer, &buffer_size, fp) != -1) {
        ptr = strtok(buffer, delim); //printf("%s ", ptr);
        if (strcmp(ptr,"REEF") == 0) {
            ptr = strtok(NULL, delim); // polygon file name
            spolygon_read_file(&mod->poly[mod->npoly], ptr);
            mod->npoly++;
        }
    }
    rewind(fp);
    
    
    // split particles among processors in distributed run
    int particle_pe_start[npes];
    int particle_pe_end[npes];
    particle_pe_start[0] = 0;
#ifdef _MPI
    double dpe = floor(mod->np_global / (double) npes);
    particle_pe_end[0] = dpe-1;
    for (i=1; i<npes; i++) {
        particle_pe_start[i] = particle_pe_end[i-1] + 1;
        particle_pe_end[i] = particle_pe_start[i] + dpe;
    }
#endif
    particle_pe_end[npes-1] = mod->np_global-1;
    mod->np = particle_pe_end[myid] - particle_pe_start[myid] + 1;
//    for (i=0; i<npes; i++) {
//        printf("myid: %d || np: %d || pe: %d || start: %d || end: %d \n",myid,mod->np,i,particle_pe_start[i],particle_pe_end[i]);
//    }
//    exit(-1);
    
    // allocate particle storage
    mod->p = (SPARTICLE *)malloc(sizeof(SPARTICLE) * mod->np);
    mod->p0 = (SVECT *)malloc(sizeof(SVECT) * mod->np);
    mod->pdpl = (SVECT *)malloc(sizeof(SVECT) * mod->np);
    
    // initialize particles
    sparticle_init(mod->p,mod->np);
    
    // now read to store particle data
    mod->np_global = 0, mod->np = 0;
    while(getline(&buffer, &buffer_size, fp) != -1) {
        ptr = strtok(buffer, delim);
        //printf("%s ", ptr);
        bool pnode = false;
        bool pbary = false;
        if (strcmp(ptr,"PNODE") == 0) pnode = true;
        if (strcmp(ptr,"PBARY") == 0) pbary = true;
        if (pnode || pbary) {
            assert((pnode && pbary) == false); // they can't both be true
            ptr = strtok(NULL, delim); // the total # of particles started at this node/element
            //printf("%s ", ptr);
            node_np = atoi(ptr);
            ptr = strtok(NULL, delim); // the type of particles created
            //printf("%s ", ptr);
            if (strcmp(ptr,"GENERAL") == 0) ptype = GENERAL;
            if (strcmp(ptr,"OYSTER") == 0) ptype = OYSTER;
            ptr = strtok(NULL, delim); // read node/elem # on hydraulic grid
            //printf("%s ", ptr);
            if (pnode) {
                nd = atoi(ptr) - 1;
                ie = mod->grid->nc_elems[nd][0]; // choose 1st element from hydraulic grid connectivity table containing the given node
            } else { // pbary == true
                ie = atoi(ptr) - 1;
            }
            for (i=0;i<node_np; i++) {
                if (mod->np_global >= particle_pe_start[myid] && mod->np_global <= particle_pe_end[myid]) {
                    mod->p[mod->np].gID = mod->np_global;
                    mod->p[mod->np].lID = mod->np;
                    mod->p[mod->np].age = 0.0;
                    mod->p[mod->np].ptype = ptype;
                    mod->p[mod->np].isActive = YES;
                    if (pnode) {
                        mod->p[mod->np].nodeID = nd;
                        mod->p[mod->np].r.x = mod->grid->node[nd].x;
                        mod->p[mod->np].r.y = mod->grid->node[nd].y;
                        mod->p[mod->np].r.z = mod->grid->node[nd].z;
                    } else { // pbary == true
                        double bx = 0, by = 0, bz = 0;
                        for (int it=0; it<3; it++) {
                            int nd_ID = mod->grid->elem2d[ie].nodes[it];
                            bx += mod->grid->node[nd_ID].x;
                            by += mod->grid->node[nd_ID].y;
                            bz += mod->grid->node[nd_ID].z;
                        }
                        mod->p[mod->np].r.x = bx / 3.0;
                        mod->p[mod->np].r.y = by / 3.0;
                        mod->p[mod->np].r.z = bz / 3.0;

                    }
                    mod->p[mod->np].elemID = ie;
                    mod->p0[mod->np] = mod->p[mod->np].r;
                    
                    // make sure point is in element
                    double weights[3];
                    int lnd = UNSET_INT, ledge = UNSET_INT, flag = UNSET_INT;
                    flag = compute_interpolate_weights_2D_triangle(mod->grid,ie,mod->p[mod->np].r.x,mod->p[mod->np].r.y,weights,&lnd,&ledge);
                    assert(flag == 1);
                    
                    
                    mod->np++;
                    //printf("id: %d r: %f %f %f\n",mod->np,mod->p[mod->np].r.x,mod->p[mod->np].r.y,mod->p[mod->np].r.z);
                }
                mod->np_global++;
            }
        }
    }
    
    if (RUN_VERBOSE) {
        printf("myid: %d || Total # of particles: %d\n",myid,mod->np);
        printf("nelems connected: %d\n",mod->grid->nc_nelems[mod->p[0].nodeID]);
        for (i=0; i<mod->grid->nc_nelems[mod->p[0].nodeID]; i++) {
            printf("elem %d connected\n",mod->grid->nc_elems[mod->p[0].nodeID][i]);
        }
        for (i=0; i<mod->np; i++) printf("particle %d || nd: %d || ie: %d || {%20.15e, %20.15e, %20.15e}\n",i,mod->p[i].nodeID,mod->p[i].elemID,mod->p[i].r.x,mod->p[i].r.y,mod->p[i].r.z);
        printf("nd: %d || {%20.15e, %20.15e, %20.15e}\n",mod->p[0].nodeID,mod->grid->node[mod->p[0].nodeID].x,mod->grid->node[mod->p[0].nodeID].y,mod->grid->node[mod->p[0].nodeID].z);
    }
    
    free(buffer);
    fclose(fp);
#ifdef _MPI
     MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (myid == 0) {fflush(stdout); printf("FINISHED OPENING MODEL INPUT FILE: %s\n",pfile); }
        
    return 0;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void smodel_print_ts(SMODEL *mod, int it) {
    int ip,ipe;
    SVECT abs_error, rel_error, pa;
    double abs_error_mag, rel_error_mag, max_t_abs_error_mag, max_t_rel_error_mag;
    int npes=1,myid=0;
#ifdef _MPI
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    // print particle data to screen
    if (mod->screen_print_timestep) sparticle_printScreen(mod->p,mod->np,mod->time);
    if ((it - 1) % mod->spool != 0) return;
    
    // print particle data to file
    for (ip=0; ip<mod->np; ip++) {
        SVECT r0  = svect_denormalize(mod->normalize, mod->grid, mod->p0[ip]);
        SVECT rIP = svect_denormalize(mod->normalize, mod->grid, mod->p[ip].r);
        fprintf(mod->fp_out,"%d, %10.5e, %10.5e, %10.5e, %10.5e\n",mod->p[ip].gID+1,mod->time,rIP.x,rIP.y,rIP.z);
        if (mod->fp_error != NULL)  {
            assert(mod->get_analytic_p != NULL);
            //pa = get_errors_test(mod->p[ip].r, mod->p0[ip], mod->time, mod->get_analytic_p, &abs_error, &rel_error);
            pa = get_errors_test(rIP, r0, mod->p[ip].age, mod->get_analytic_p, &abs_error, &rel_error);
            abs_error_mag = svect_mag(abs_error);
            rel_error_mag = svect_mag(rel_error);
            if (abs_error_mag > mod->max_t_abs_error_mag) mod->max_t_abs_error_mag = abs_error_mag;
            if (rel_error_mag > mod->max_t_rel_error_mag) mod->max_t_rel_error_mag = rel_error_mag;
            // write errors to file :: particle id || t || x || y || z || xa || ya || za || abs_error || rel_error || max_t_abs || max_t_rel
            fprintf(mod->fp_error,"%d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
                    mod->p[ip].gID+1,
                    mod->time,
                    rIP.x,
                    rIP.y,
                    rIP.z,
                    pa.x,pa.y,pa.z,
                    abs_error_mag,
                    rel_error_mag,
                    mod->max_t_abs_error_mag,
                    mod->max_t_rel_error_mag
                    );
            
//            printf("%d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
//                    mod->p[ip].gID+1,
//                    mod->time,
//                    mod->p[ip].r.x,
//                    mod->p[ip].r.y,
//                    mod->p[ip].r.z,
//                    pa.x,pa.y,pa.z,
//                    abs_error_mag,
//                    rel_error_mag,
//                    mod->max_t_abs_error_mag,
//                    mod->max_t_rel_error_mag
//                    );
        }
    }
    
    
#ifdef _ADH_HDF5
    if (RUN_VERBOSE) printf("smodel_print_ts || Writing XDMF grid data output files\n");
    // NOTE :: all pe's have same grid
    // NOTE :: really want the velocity at the current time here, this is at time = mod->time - dt
    if (myid == 0) ps_print_xdmf(mod->grid, mod->vel, mod->dpl, mod->time);
    if (RUN_VERBOSE) printf("smodel_print_ts || Writing XDMF particle data output files\n");
    ps_print_particle_xdmf(&mod->hdf5, mod->pdpl, mod->np, mod->time);
#endif
    
}

void smodel_read_opt(SMODEL * mod) {
    // open model parameter file
    char pfile[100];
    strcpy(pfile,mod->file_base);
    strcat(pfile, ".in");
    FILE *fp = fopen(pfile,"r");

    // read parameters
    size_t buffer_size = 80;
    char *buffer = malloc(buffer_size * sizeof(char));
    char delim[] = " ";  // line delimiter
    char *ptr = NULL;    // to read the card on line

    // read each line and print it to the screen
    while(getline(&buffer, &buffer_size, fp) != -1) {
        ptr = strtok(buffer, delim);
        if (strcmp(ptr,"NORMALIZE") == 0) {
            ptr = strtok(NULL, delim);
            mod->normalize = (bool)atoi(ptr);
        }
    }

    fclose(fp);
}

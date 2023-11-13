/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     The FE engine driver for a given super model
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gary Brown
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gajanan Choudhary, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  \returns YES for a good calculation and NO for a bad calculation
 *
 * @param[in] SSUPER_MODEL *sm :: ONE  Super Model
 * \note CJT \:: matrix size will grow with refinement. It never shrinks with unrefinement!
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
extern double DEBUG_TIME;
static double total_time_mass_flux = 0.;


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int fe_main(SSUPER_MODEL *sm) {
    int i, inode = 0, block = 0, isys = 0, k, SUCCESS = NO, isubmodel = 0, itrns = 0;
    SMODEL *mod = NULL;
    
    // CJT :: Hardwire for testing initial residuals with finite solutions
    //    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++){
    //        mod = &(sm->submodel[isubmodel]);
    //        if (mod->flag.SW3_FLOW == ON) {
    //            for (i=0; i<mod->grid->nnodes; i++) {
    //                mod->sw->d3->vel[i].x = 0.1;
    //                mod->sw->d3->vel[i].y = 0.05;
    //                mod->sw->d3->vel[i].z = 0.001;
    //
    //                mod->sw->d3->old_vel[i].x = 0.011;
    //                mod->sw->d3->old_vel[i].y = 0.0051;
    //                mod->sw->d3->old_vel[i].z = 0.00011;
    //
    //                mod->sw->d3->older_vel[i].x = 0.012;
    //                mod->sw->d3->older_vel[i].y = 0.0052;
    //                mod->sw->d3->older_vel[i].z = 0.00012;
    //
    //                mod->sw->d3->displacement[i] += .01;
    //                mod->sw->d3->old_displacement[i] += .011;
    //                mod->sw->d3->older_displacement[i] += .012;
    //            }
    //        } else if (mod->flag.SW2_FLOW == ON) {
    //            for (i=0; i<mod->grid->nnodes; i++) {
    //                mod->sw->d2->vel[i].x = 0.01;
    //                mod->sw->d2->vel[i].y = 0.005;
    //
    //                mod->sw->d2->old_vel[i].x = 0.011;
    //                mod->sw->d2->old_vel[i].y = 0.0051;
    //
    //                mod->sw->d2->older_vel[i].x = 0.012;
    //                mod->sw->d2->older_vel[i].y = 0.0052;
    //
    //                mod->sw->d2->head[i] += .001;
    //                mod->sw->d2->old_head[i] += .0011;
    //                mod->sw->d2->older_head[i] += .0012;
    //            }
    //        }
    //    }
    
    // start debug time after
    int DEBUG = OFF;
    DEBUG_TIME = 15e8;
    
    // if time-step fails, decrease ts, otherwise enlarge to user dt
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++){
        sm->submodel[isubmodel].flag.TIME_ADAPT_FAIL = NO;
    }
    
    // Is this a monolithic coupling run?
    int monolithic_run = NO;
    if (sm->nsubmodels > 1) monolithic_run = YES;
    int flag_DW = OFF, flag_SW_2D = OFF, flag_SW_3D = OFF, flag_NS_3D = OFF, flag_GW_3D = OFF;
    
    
    //****************************************************************************
    //****************************************************************************
    // allocates FE matrix memory if necessary -----------------------------------
    
#ifdef _PETSC
    int n, ierr;
    int resident_id, resident_pe;
    // Check if solver exists and create if not
    // TODO: See if we should move this elsewhere
    if(sm->ksp == PETSC_NULL){
        ierr = KSPCreate(PETSC_COMM_WORLD, &(sm->ksp));
        ierr = KSPSetFromOptions(sm->ksp);
    }
#endif

    // TODO: Remove the hard-coding for single submodel runs
    // in the local to global map and for the grid refinement check.
    
    //****************************************************************************
    // PETSC CODE HERE
#ifdef _PETSC
        
    // Currently hard-coded for a single model - SAM
    mod = &(sm->submodel[0]);

    // Determine if the grid was updated
    // (refined or unrefined) on any pe
#ifdef _MESSG
    int REF_FLAG_RECV=0, UNREF_FLAG_RECV=0;
    int REF_FLAG = mod->flag.GRID_REFINED;
    int UNREF_FLAG = mod->flag.GRID_UNREFINED;
    MPI_Allreduce(&REF_FLAG, &REF_FLAG_RECV, 1, MPI_INT, MPI_MAX, mod->grid->smpi->ADH_COMM);
    MPI_Allreduce(&UNREF_FLAG, &UNREF_FLAG_RECV, 1, MPI_INT, MPI_MAX, mod->grid->smpi->ADH_COMM);
#else
    REF_FLAG_RECV = mod->flag.GRID_REFINED;
    UNREF_FLAG_RECV = mod->flag.GRID_UNREFINED;
#endif
    //if (sm->total_my_nnodes > sm->total_nnodes_matrix) {
    if (REF_FLAG_RECV==YES || UNREF_FLAG_RECV==YES || sm->A == PETSC_NULL){

        printf("REF = %i\n",REF_FLAG_RECV);
        printf("UNREF = %i\n",UNREF_FLAG_RECV);
        n = sm->my_nnodes*sm->max_nsys;

        // Check if Jacobian, sol, and residual have already been created.
        // If so, destroy each of them before creating new PETSc objects.
        if(sm->A != PETSC_NULL){
            printf("\n\nDestroying old matrix\n\n");
            ierr = MatDestroy(&(sm->A));
            ierr = KSPReset(sm->ksp);
            ierr = KSPSetFromOptions(sm->ksp);
        }
        if(sm->P != PETSC_NULL){
            ierr = MatDestroy(&(sm->P));
        }

        if(sm->sol != PETSC_NULL){
            ierr = VecDestroy(&(sm->sol));
        }
        if(sm->residual != PETSC_NULL){
            ierr = VecDestroy(&(sm->residual));
        }

        // Create Jacobian matrix
        ierr = MatCreate(PETSC_COMM_WORLD, &(sm->A));
        ierr = MatSetSizes(sm->A,n,n,PETSC_DETERMINE,PETSC_DETERMINE);
        ierr = MatSetFromOptions(sm->A);
        ierr = MatSetBlockSize(sm->A, sm->max_nsys);
        ierr = MatSetUp(sm->A);

        // Create preallocator matrix
        if(sm->PREALLOC_FLAG == ON){
            ierr = MatCreate(PETSC_COMM_WORLD, &(sm->P));
            ierr = MatSetType(sm->P, MATPREALLOCATOR);
            ierr = MatSetSizes(sm->P,n,n,PETSC_DETERMINE,PETSC_DETERMINE);
            ierr = MatSetBlockSize(sm->P, sm->max_nsys);
            ierr = MatSetUp(sm->P);
        }
        
        //// Currently hard-coded for a single model - SAM
        //mod = &(sm->submodel[0]);

        // Get ownership ranges of each proc
        ierr = MatGetOwnershipRange(sm->A, &(sm->Istart), &(sm->Iend));
        ierr = MatGetOwnershipRanges(sm->A, &(sm->ownership_range));
        
        // Create the Newtonian residual array
        ierr = VecCreate(PETSC_COMM_WORLD, &(sm->residual));
        ierr = VecSetSizes(sm->residual, n, PETSC_DETERMINE);
        ierr = VecSetFromOptions(sm->residual);
        ierr = VecSetBlockSize(sm->residual, sm->max_nsys);
        ierr = VecSetUp(sm->residual);
        
        // Create the Newton solution (array of solution increments)
        ierr = VecCreate(PETSC_COMM_WORLD, &(sm->sol));
        ierr = VecSetSizes(sm->sol, n, PETSC_DETERMINE);
        ierr = VecSetFromOptions(sm->sol);
        ierr = VecSetUp(sm->sol);

        // Local to global mapping object for convenient indexing when
        // setting values in residual/load/assemble functions
        ISLocalToGlobalMapping ltog_map;
       
        // TODO: Change the local to global mapping to work for
        // multiple models - SAM
        // Create local to global mapping for PETSc objects
        int petsc_global_indices[sm->nnodes];

        // Loop over the residential nodes first and add offset
        for(inode = 0; inode < sm->my_nnodes; inode++){
            petsc_global_indices[inode] = inode + (sm->Istart/sm->max_nsys);
        }
        // Loop over the ghost nodes and map them to the PETSc
        // index using information from the residential processors
        for(inode = sm->my_nnodes; inode < sm->nnodes; inode++){
            // Get the PETSc-global index of the non-residential node
            resident_id = mod->grid->node[inode].resident_id;
            resident_pe = mod->grid->node[inode].resident_pe;
            
            petsc_global_indices[inode] =
                resident_id + (sm->ownership_range[resident_pe]/sm->max_nsys);
        }

        // Create the local to global PETSc map
        ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, sm->max_nsys,
                sm->nnodes, petsc_global_indices, PETSC_COPY_VALUES, &(ltog_map));

        // Apply the mapping to A and residual
        ierr = MatSetLocalToGlobalMapping(sm->A, ltog_map, ltog_map);
        if(sm->PREALLOC_FLAG == ON){
            ierr = MatSetLocalToGlobalMapping(sm->P, ltog_map, ltog_map);
        }
        ierr = VecSetLocalToGlobalMapping(sm->residual, ltog_map);

        // Destroy mapping object now that it has been applied
        ierr = ISLocalToGlobalMappingDestroy(&(ltog_map));
        

        if(sm->PREALLOC_FLAG == ON){
            // Preallocate space for the matrix
            fe_preallocate(sm, 0);
            //ierr = MatView(sm->A,PETSC_VIEWER_STDOUT_WORLD);
            
            // Temporarily turn off allocation errors
            ierr = MatSetOption(sm->A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
        }
        
        // Below needs to be converted into PETSc vectors
        // TODO: Ask Corey about some of these
        //sm->solver_info.node_block = (int *) tl_realloc(sizeof(int), sm->total_nnodes, sm->total_nnodes_matrix, sm->solver_info.node_block);
        //sm->bc_mask = (int *) tl_realloc(sizeof(int), sm->total_nnodes * sm->max_nsys, sm->old_bc_mask_size, sm->bc_mask); // CJT :: prob not needed!
        sm->bc_mask = (int *) tl_realloc(sizeof(int), sm->total_nnodes * sm->max_nsys, sm->old_bc_mask_size, sm->bc_mask); // CJT :: prob not needed!
        sm->old_bc_mask_size = sm->total_nnodes * sm->max_nsys;
//#ifdef _MESSG
//        for (inode = sm->total_nnodes_matrix; inode < sm->total_nnodes; inode++) {
//            sm->solver_info.node_block[inode] = sm->supersmpi->myid;
//        }
//#endif
        // No longer use this for adaption condition above.
        // Could get rid of this line.
        sm->total_nnodes_matrix = sm->total_my_nnodes;
    }
#else
    //****************************************************************************
    // UMFPACK CODE HERE
    if (sm->total_nnodes > sm->total_nnodes_matrix) {
        sm->solver_info.node_block = (int *) tl_realloc(sizeof(int), sm->total_nnodes, sm->total_nnodes_matrix, sm->solver_info.node_block);
        sm->bc_mask = (int *) tl_realloc(sizeof(int), sm->total_nnodes * sm->max_nsys, sm->total_nnodes_matrix * sm->max_nsys, sm->bc_mask); // CJT :: prob not needed!
        sm->residual = (double *) tl_realloc(sizeof(double), sm->total_nnodes * sm->max_nsys, sm->total_nnodes_matrix * sm->max_nsys, sm->residual);
        sm->sol = (double *) tl_realloc(sizeof(double), sm->total_nnodes * sm->max_nsys, sm->total_nnodes_matrix * sm->max_nsys,sm->sol);
        sm->scale_vect = (double *) tl_realloc(sizeof(double), sm->total_nnodes * sm->max_nsys,sm->total_nnodes_matrix * sm->max_nsys, sm->scale_vect);
        sm->diagonal = (double *) tl_realloc(sizeof(double), sm->total_nnodes * sm->max_nsys_sq, sm->total_nnodes_matrix * sm->max_nsys_sq, sm->diagonal);
        sm->matrix = (SPARSE_VECT *) tl_realloc(sizeof(SPARSE_VECT), sm->total_nnodes, sm->total_nnodes_matrix, sm->matrix);
        
        // initialize
        sarray_init_dbl(sm->sol, sm->max_nsys*sm->total_nnodes);
        sarray_init_value_dbl(sm->scale_vect, sm->max_nsys*sm->total_nnodes, 1.0);
        sarray_init_dbl(sm->diagonal, sm->max_nsys_sq*sm->total_nnodes);
        for (inode = sm->total_nnodes_matrix; inode < sm->total_nnodes; inode++) {
#ifdef _MESSG
            sm->solver_info.node_block[inode] = sm->supersmpi->myid;
#endif
            sm->matrix[inode].value           = (double *) tl_alloc(sizeof(double), SPV_BLOCK * sm->max_nsys_sq);
            sm->matrix[inode].index           = (int *)    tl_alloc(sizeof(int),    SPV_BLOCK);
            sm->matrix[inode].max_size        = SPV_BLOCK;
            sm->matrix[inode].size            = 0;
        }
        sm->total_nnodes_matrix = sm->total_nnodes;
    }
#ifdef _MESSG
    comm_update_int(sm->solver_info.node_block, 1, sm->supersmpi);
#endif
#endif
    
    // Creating a separate bc_mask for all submodels.
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++){
        if(sm->submodel[isubmodel].proc_flag==1){
            mod = &(sm->submodel[isubmodel]);
            if (mod->grid->nnodes > mod->grid->nnodes_matrix) {
                mod->bc_mask = (int *) tl_realloc(sizeof(int), mod->grid->nnodes * mod->max_nsys, mod->grid->nnodes_matrix * mod->max_nsys, mod->bc_mask);
                mod->grid->nnodes_matrix = mod->grid->nnodes;
            }
        }
    }
   
#ifndef _PETSC 
    // Loop over all the nodes here, since block can completely change
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++){
        if(sm->submodel[isubmodel].proc_flag==1){
            mod = &(sm->submodel[isubmodel]);
            for (inode = 0; inode < mod->grid->nnodes; inode++) {
                block = mod->grid->node[inode].block;
#ifdef _MESSG
                sm->solver_info.node_block[mod->fmap[inode]] = sm->supersmpi->myid;
#else
                sm->solver_info.node_block[mod->fmap[inode]] = block;
#endif
            }
        }
    }
#ifdef _MESSG
    comm_update_int(sm->solver_info.node_block, 1, sm->supersmpi);
#endif
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // Initialize hydrodynamic submodels
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
        if(sm->submodel[isubmodel].proc_flag==1){
            mod = &(sm->submodel[isubmodel]); // alias
            
            // Diffusive Wave
            if (mod->flag.DIFFUSIVE_WAVE && mod->flag.GW_FLOW==OFF){
                flag_DW = ON;
                mod->nsys = 1;
                mod->nsys_sq = 1;
                fe_diffusive(mod);
            }
            
#ifdef _ADH_GROUNDWATER
            // Ground water
            if (mod->flag.GW_FLOW){
#ifdef _DWGW_COUPLING
                // Diffusive wave - Ground water coupling
                if (mod->flag.DIFFUSIVE_WAVE){
                    flag_DW = ON;
                    mod->nsys = 1;
                    mod->nsys_sq = 1;
                    fe_diffusive(mod);
                }
#endif
                flag_GW_3D = ON;
                mod->nsys = 1;
                mod->nsys_sq = 1;
                fe_gw(mod);
            }
#endif
            
            // 2D Shallow Water
            else if (mod->flag.SW2_FLOW){
                flag_SW_2D = ON;
                mod->nsys = 3;
                mod->nsys_sq = 9;
                if (debug.no_hydro == OFF) {
                    fe_sw2(mod);
                } else {
                    printf("\nAssigning Predetermined hydro values\n");
                    for (inode = 0; inode < mod->grid->nnodes; inode++) {
                        mod->sw->d2->vel[inode].x = debug.u_vel;
                        mod->sw->d2->vel[inode].y = debug.v_vel;
                        mod->sw->d2->head[inode] = debug.hard_disp;
                        mod->sw->d2->old_vel[inode].x = debug.u_vel;
                        mod->sw->d2->old_vel[inode].y = debug.v_vel;
                        mod->sw->d2->old_head[inode] = debug.hard_disp;
                        mod->sw->d2->older_vel[inode].x = debug.u_vel;
                        mod->sw->d2->older_vel[inode].y = debug.v_vel;
                        mod->sw->d2->older_head[inode] = debug.hard_disp;
                    }
                }
            }
            
            // 3D Shallow Water
            else if (mod->flag.SW3_FLOW){
                flag_SW_3D = ON;
                mod->nsys = 3;
                mod->nsys_sq = 9;
                if (debug.no_hydro == OFF) {
                    fe_sw3(mod);
                } else {
                    printf("\nAssigning Predetermined hydro values\n");
                    double scale = 1.;
                    for (inode = 0; inode < mod->grid->nnodes; inode++) {
                        mod->sw->d3->vel[inode].x = scale * debug.u_vel;
                        mod->sw->d3->vel[inode].y = scale * debug.v_vel;
                        mod->sw->d3->vel[inode].z = scale * debug.w_vel;
                        mod->sw->d3->displacement[inode] = debug.hard_disp;
                        mod->sw->d3->old_vel[inode].x = scale * debug.u_vel;
                        mod->sw->d3->old_vel[inode].y = scale * debug.v_vel;
                        mod->sw->d3->old_vel[inode].z = scale * debug.w_vel;
                        mod->sw->d3->old_displacement[inode] = debug.hard_disp;
                        mod->sw->d3->older_vel[inode].x = scale * debug.u_vel;
                        mod->sw->d3->older_vel[inode].y = scale * debug.v_vel;
                        mod->sw->d3->older_vel[inode].z = scale * debug.w_vel;
                        mod->sw->d3->older_displacement[inode] = debug.hard_disp;
                    }
                }
            }
            
            // Navier Stokes Flow
            else if (mod->flag.NS_FLOW) {
                flag_NS_3D = ON;
                mod->nsys = 4;
                mod->nsys_sq = 16;
                if (debug.no_hydro == OFF) {
                    fe_ns3(mod);
                } else {
                    printf("\nAssigning Predetermined hydro values\n");
                    double scale = 1.;
                    for (inode = 0; inode < mod->grid->nnodes; inode++) {
                        mod->ns->d3->vel[inode].x = scale * debug.u_vel;
                        mod->ns->d3->vel[inode].y = scale * debug.v_vel;
                        mod->ns->d3->vel[inode].z = scale * debug.w_vel;
                        mod->ns->d3->displacement[inode] = debug.hard_disp;
                        mod->ns->d3->old_vel[inode].x = scale * debug.u_vel;
                        mod->ns->d3->old_vel[inode].y = scale * debug.v_vel;
                        mod->ns->d3->old_vel[inode].z = scale * debug.w_vel;
                        mod->ns->d3->old_displacement[inode] = debug.hard_disp;
                        mod->ns->d3->older_vel[inode].x = scale * debug.u_vel;
                        mod->ns->d3->older_vel[inode].y = scale * debug.v_vel;
                        mod->ns->d3->older_vel[inode].z = scale * debug.w_vel;
                        mod->ns->d3->older_displacement[inode] = debug.hard_disp;
                    }
                }
            }
        }
    }
    mod = NULL; /* To prevent silly mistakes in this file later. */
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // HYDRODYNAMICS
    
    if (debug.no_hydro == OFF) {
        if (sm->nsubmodels > 1) {
            
            /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            /* MONOLITHIC MULTI-MODEL COUPLING */
            /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

            // 2D/3D Shallow Water
            sm->nsys = 3;
            sm->nsys_sq = 9;
            SUCCESS = fe_sw_hybrid_solve(sm, 0);
#ifdef _MESSG
            // Make sure everyone is on the same page with the Newton solve success.
            // Sometimes, the SUCESS of hvel pes, which works on all PES will leak over
            // into wvel.
            MPI_Allreduce(&SUCCESS, &SUCCESS, 1, MPI_INT, MPI_MIN, cstorm_comm);
#endif
            if (SUCCESS == NO) return NO;
            
            /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            /* SINGLE MODEL RUNS */
            /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            
        } else {
            
            mod = &(sm->submodel[0]);
            sm->nsys = mod->nsys;
            sm->nsys_sq = mod->nsys_sq;
            
            // Diffusive Wave
            if (mod->flag.DIFFUSIVE_WAVE == ON && mod->flag.GW_FLOW == OFF) {
                SUCCESS = fe_diffusive_solve(sm,0);
                if (SUCCESS == NO) return NO;
            }
            
#ifdef _ADH_GROUNDWATER
            // Ground water
            if (mod->flag.GW_FLOW == ON){
#ifdef _DWGW_COUPLING
                if (mod->flag.DIFFUSIVE_WAVE == ON){ // GW-DW Coupling
                    SUCCESS = fe_dwgw_hybrid_solve(sm,0);
                    if (SUCCESS == NO) return NO;
                }
                else { // Solo GW
                    SUCCESS = fe_gw_solve(sm,0);
                    if (SUCCESS == NO) return NO;
                }
#else
                SUCCESS = fe_gw_solve(sm,0);
                if (SUCCESS == NO) return NO;
#endif
            }
#endif

            // 2D Shallow Water
            if (mod->flag.SW2_FLOW == ON && mod->flag.DIFFUSIVE_WAVE == OFF) {
                SUCCESS = fe_sw2_solve(sm,0);
                if (SUCCESS == NO) return NO;
            }
            
            // 3D Shallow Water
            if (mod->flag.SW3_FLOW == ON) {
                SUCCESS = fe_sw3_solve(sm,0);
                if (SUCCESS == NO) return NO;
            }
            
            // 3D Navier Stokes
            if (mod->flag.NS_FLOW == ON) {
                SUCCESS = fe_ns3_solve(sm,0);
                if (SUCCESS == NO) return NO;
            }
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GENERAL CONSTITUENT AND SEDIMENT TRANSPORT
    
    if (sm->nsubmodels > 1) {
        
        /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /* MONOLITHIC MULTI-MODEL COUPLING */
        /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
        // 2D/3D Coupled Transport
        SUCCESS = fe_transport_hybrid(sm);
        if (SUCCESS == NO) return NO;
        
        /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /* SINGLE MODEL RUNS */
        /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
    } else {
        
        mod = &(sm->submodel[0]);
        
        // 2D/3D General Constituent Transport
        if (mod->flag.TRANSPORT) {
            if (mod->flag.SW2_TRANSPORT || mod->flag.SW3_TRANSPORT || mod->flag.NS2_TRANSPORT || mod->flag.NS3_TRANSPORT) {
                SUCCESS = fe_transport(sm,0);
                if (SUCCESS == NO) return NO;
            }
        }
        
        // Sediment Transport
#ifdef _SEDIMENT
        if (mod->flag.SEDIMENT) {
            mod->is_sediment_running = TRUE;
            int ndim = 0;
            if (mod->flag.SW2_FLOW) {ndim = 2;}
            else if (mod->flag.SW3_FLOW) {ndim = 3;}
            SUCCESS = fe_sediment_transport(mod, ndim);
            if (SUCCESS == NO){
                mod->is_sediment_running = FALSE; /* We failed so set the flag back to FALSE */
                return NO;
            }
            mod->is_sediment_running = FALSE;
        }
#endif
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Calculate total mass errors
    if(sm->supersmpi->myid==0) printf("\n");
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
        if(sm->submodel[isubmodel].proc_flag==1){
            mod = &(sm->submodel[isubmodel]);
            if (screen_output.grid_mass_error == ON) {
                if (mod->grid->ndim == 2){
                    mod->grid_mass_error = tl_find_grid_mass_error_elem2d(mod->density, mod->sw->d2->head, mod->sw->d2->vel, mod->grid, mod->flag, mod->initial_grid_mass, mod->series_head, mod->str_values, mod->dt, &total_time_mass_flux);
                    
                    // transport constituent mass :: CJT :: need to store consituent initial mass here!
                    //for (itrns=0; itrns<sm->ntransport; itrns++){
                    //    mod->grid_mass_error = tl_find_grid_mass_error_elem2d(1.0, mod->con[itrns].concentration, mod->grid, mod->flag, mod->initial_grid_mass, mod->series_head, mod->str_values, mod->dt);
                    //}
                }
                else if (mod->grid->ndim == 3){
                    
                    SVECT *vel;
                    double *dpl, *old_dpl, *older_dpl, new_grid_mass = 0.;
                    if (flag_SW_3D == ON) {
                        vel = mod->sw->d3->vel;
                        dpl = mod->sw->d3->displacement;
                        old_dpl = mod->sw->d3->old_displacement;
                        older_dpl = mod->sw->d3->older_displacement;
                    } else if (flag_NS_3D == ON) {
                        vel = mod->ns->d3->vel;
                        dpl = mod->ns->d3->displacement;
                        old_dpl = mod->ns->d3->old_displacement;
                        older_dpl = mod->ns->d3->older_displacement;
                    } else if (flag_GW_3D == ON) {
                        // Gajanan gkc: Aren't all dpl's supposed to be 0 for GW? Only vel needs to be assigned I guess?
                        //vel = mod->sgw->vel;
                        //dpl = mod->sgw->displacement???;
                        //old_dpl = mod->sgw->old_displacement???;
                        //older_dpl = mod->sgw->older_displacement???;
                    }
                    
                    if (flag_GW_3D == ON){
                        //mod->grid_mass_error = tl_find_grid_mass_error_tet_gw(mod->initial_grid_mass,
                        //        mod->grid, gw->gw_phead, gw->elem_3d_data, mod->series_head,
                        //        mod->str_values, mod->dt);
                    }
                    else{
                        mod->grid_mass_error = tl_find_3d_grid_mass_error(mod->str_values,
                                mod->series_head, mod->initial_grid_mass, mod->density,
                                mod->grid, vel, dpl, old_dpl, older_dpl, mod->dt,
                                &new_grid_mass, &total_time_mass_flux);
                    }
                    
                    // transport constituent mass :: CJT :: need to store consituent initial mass here!
                    //for (itrns=0; itrns<sm->ntransport; itrns++){
                    //    mod->grid_mass_error = tl_find_grid_error_elem3d(mod->initial_grid_mass, mod->grid, mod->con[itrns].concentration, mod->series_head, mod->str_values, mod->dt);
                    //}
                }
                //printf("  grid_mass_error submodel[%i]: %10.5e\n", isubmodel, mod->grid_mass_error);
            }
        }
    }
   
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Calculate approximate continuity error over elements for adaption
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
        if(sm->submodel[isubmodel].proc_flag==1) {
            mod = &(sm->submodel[isubmodel]);
            if ((mod->flag.SW2_FLOW || mod->flag.DIFFUSIVE_WAVE)
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
                    && (mod->flag.GW_FLOW==OFF)
#endif
#endif
                    ){
                ssw_2d_calculate_elem_error(mod->sw->d2, mod->grid, mod->mat, mod->dt);
                
                if (mod->flag.SW2_TRANSPORT) {
                    scon_calculate_elem2d_error(mod->ntransport, mod->con, mod->sw->d2->vel, mod->sw->d2->head, mod->grid, mod->mat, mod->dt);
                }
            } else if (mod->flag.SW3_FLOW) {
                ssw_3d_calculate_elem_error(mod->sw->d3, mod->grid, mod->mat, mod->dt);
                
                if (mod->flag.SW3_TRANSPORT) {
                    scon_calculate_elem3d_error(mod->ntransport, mod->con, mod->sw->d3->vel, mod->sw->d3->displacement,
                                                mod->sw->d3->old_displacement, mod->sw->d3->older_displacement, mod->grid, mod->mat, mod->dt);
                }
            } else if (mod->flag.NS3_FLOW) {
                sns_3d_calculate_elem_error(mod->ns->d3, mod->grid, mod->mat, mod->dt);
                
                if (mod->flag.NS3_TRANSPORT) {
                    scon_calculate_elem3d_error(mod->ntransport, mod->con, mod->ns->d3->vel, mod->ns->d3->displacement,
                                                mod->ns->d3->old_displacement, mod->ns->d3->older_displacement, mod->grid, mod->mat, mod->dt);
                }
#ifdef _ADH_GROUNDWATER
            } else if (mod->flag.GW_FLOW) {
                sgw_3d_calculate_elem_error(mod->sgw, mod->grid, mod->mat, mod->dt);
                
                //if (mod->flag.GW_TRANSPORT) {
                //    scon_calculate_elem3d_error(mod->ntransport, mod->con, mod->sgw->vel???, mod->sgw->displacement???,
                //                                mod->sgw->old_displacement???, mod->sgw->older_displacement???, mod->grid, mod->mat, mod->dt);
                //}
#endif
            }
        }
    }

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* gkc adding following line in July 2016. */
    /* Need 3D model depth to be updated before writing the output to a file. */
    /* Since fe_sw3_hvel_inc.c isn't doing this, we'll do it here. */
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
        if(sm->submodel[isubmodel].proc_flag==1){
            mod = &(sm->submodel[isubmodel]); // alias
            if (mod->flag.SW3_FLOW){
                tl_calculate_depthavgvel(mod->grid, mod->sw->d3);
                
                // Convert elemental eddy viscosity to nodal viscosity
                if (mod->file_output.hyd_vis == ON)
                    elem3d_to_node_double(mod->grid, mod->str_values, mod->sw->d3->hyd_viscosity , mod->grid->hyd_eddy); //GSAVANT
                if (mod->file_output.trn_dif == ON)
                    elem3d_to_node_double(mod->grid, mod->str_values, mod->sw->d3->trn_diffusivity , mod->grid->trn_diff); //GSAVANT
            }
        }
    }
    
#ifdef _DEBUG
    if (DEBUG==ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif
    
    return (SUCCESS);
}


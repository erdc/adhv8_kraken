/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file fe_newton.c This file collections functions responsible for
 *          the Netwon solve.                                                               */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static double *x0;         /* solution = u+u0 - u0 is the shift */
static int isize = 0;
static double PETSC_RTOL = 1e-12;
static double PETSC_ATOL = 1e-11;
static double PETSC_DTOL = 1e5;
static int PETSC_MAXIT = 10000;
#ifdef _PETSC
    static const char *PETSC_PRECON = PCLU;// PCILU;//PCLU;
    static const char *PETSC_KSP = KSPPREONLY;//KSPGMRES;//KSPPREONLY; KSPBICGSTAB
#endif
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     The Newton solver using FE method and central finite difference approximation
 *  for Jacobian (aka a Secant method). The function takes in a SMODEL_SUPER and solves the
 *  nonlinear problem, for time dependent problems this effectively marches the solution
 *  forward one step.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 * 
 * \returns TRUE for a good calculation and FALSE for a bad calculation
 *
 * @param[in,out] sm (SMODEL_SUPER *) - a pointer to an AdH super model
 * @param[in] isuperModel (int) - the supermodel index
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int fe_newton(SMODEL_SUPER* sm)
{

    //pointer to the grid
    SGRID *grid = sm->grid;
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG OPTIONS
#ifdef _DEBUG
    int DEBUG_FULL = OFF;       // prints everything but matrix
    int DEBUG_MATRIX = OFF;     // prints the matrix
    int DEBUG_PICKETS = OFF;    // check memory
    int DEBUG_INIT_RESID = OFF;             // print the initial residual only
    int DEBUG_STOP_AFTER_INIT_RESID = OFF;  // stop of this intitial residual print
    //if (init_fnctn == fe_sw_hybrid_wvel_init) DEBUG_FULL = ON;
    //if (init_fnctn == fe_transport_init) DEBUG_STOP_AFTER_INIT_RESID = ON;
#endif
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//figure out MPI later
#ifdef _MESSG
        MPI *smpi = grid->smpi
#endif

    
    //uses model to access all elemental routines and builds/solves linear systems
    int status;
    int check = NO;		/* newton solver check */
    int keep_chugging = NO;     /* newton solver check */
    //double next_dt = 0;         /* temp variable for the time step for STEADY STATE */
    //double ratio = 0;           /* STEADY STATE VARIABLES */
    int i;                      /* loop counter */
    //int idof;                   /* location of the maximum residual */
    int imax_dof = 0;          /* location of the maximum residual */
    int iinc_dof = 0;          /* location of the maximum ..increment? */
    int it;                     /* the nonlinear iteration counter */
    int linesearch_cuts;        /* the number of line searches */
    int solv_flag = YES;              /* flag that tells if the linear solver converged */
    int iend;                   /* loop limit */
    double resid_max_norm = 0.0;    /* the max norm of the residual */
    double resid_l2_norm = 0.0; /* the l2 norm of the residual */
    double old_resid_norm;      /* the previous norm of the residual - used for the line search */
    //double partial_max_norm;    /* the partial residual norm prior to summing over processors */
    //double partial_l2_norm;     /* the partial residual norm prior to summing over processors */
    double inc_max_norm = 0.0;  /* the max increment change */
    //static double initial_residual; /* Initial_Residual */
    //double spatial_residual_proc;   /* Individual Processor resid */
    //double lolh;
    //static int count_std;
    //int mn_node = 0, im_node = 0;
    //int tot_nnode;
    int UMFail_max = NO;
    int solv_flag_min = YES;
    int myid = 0, npes = 1, proc_count = 0;
    //convenient aliasing
    SLIN_SYS *lin_sys = sm->lin_sys;
    double *dsol = lin_sys->dsol;
    double *residual = lin_sys->residual;

    int ndofs, my_ndofs, macro_ndofs;
    ndofs = *(sm->ndofs);
    my_ndofs = *(sm->my_ndofs);
    macro_ndofs = *(sm->macro_ndofs);
#ifdef _PETSC
    int ierr = 0; // SAM HERE
    int its;
    int isize_prev;
    MatInfo info;
    if(isize<ndofs){
        isize_prev = isize;
        isize = ndofs;
        x0 = (double *) tl_realloc(sizeof(double), isize, isize_prev, x0);
    }
#ifdef _DEBUG
    PetscViewer viewer;
#endif
#endif


    
#ifdef _DEBUG
    int nnodes = grid->nnodes;
    int my_nnodes = grid->my_nnodes;
    int macro_nnodes = grid->macro_nnodes;
    if (DEBUG_FULL) {
            printf("fe_newton :: myid: %d ndim: %d :: grid->my_nnodes: %d :: grid->nnodes: %d  :: macro_nnodes: %d\n",
                   grid->smpi->myid,
                   grid->ndim,
                   grid->my_nnodes,
                   grid->nnodes,
                   grid->macro_nnodes);
    }
#endif
    //In an ideal world the residual for SW models with h <0 = 0 but this
    // must not be case, need to dig into this later
    
    // Is there a 2D SW model inside the superModel?
    //int is_2d_hydro = ON;

    //old way of checking, for now just assume this is always on
//    for (i=0; i<sm->nsubmodels; i++){
//        if (sm->submodel[i].flag.DIFFUSIVE_WAVE==ON && sm->submodel[i].flag.GW_FLOW==ON){
//            is_2d_hydro = OFF;
//            break;
//        }
//        else if (sm->submodel[i].flag.SW2_FLOW || sm->submodel[i].flag.DIFFUSIVE_WAVE) {
//            is_2d_hydro = ON;
//            break;
//        }
//    }
    

    // this array determines what nodes will be including in the norm calculations
    //is it better to allocate this and store or just compute on fly?
    //why do we include ghost nodes??
    int *include_dof = NULL;
    //if (init_fnctn != fe_sw_hybrid_wvel_init) {
    //include_node_norm = (int *) tl_alloc(sizeof(int), nnodes);
    //sarray_init_value_int(include_node_norm, nnodes, TRUE);
    //}
    

    //Mark, currently ignore
    //not part of super model
    //Solver_Info *solver = &(sm->solver_info);
    //solver->UMFail = FALSE;
    
    // aliases

    //Deal with sediment differently now
    //need to come back to this
    /*    
    int index;
#ifdef _SEDIMENT
    if (init_fnctn == fe_sw3_transport_init || init_fnctn == fe_sw2_transport_init || init_fnctn == fe_sw2_bedload_init) {
        if (sm->submodel[0].is_sediment_running) {
            index = sm->submodel[0].ised+1;
        } else {
            index = sm->submodel[0].itrns+1;
        }
    }
#else
    if (init_fnctn == fe_transport_init) {
        index = sm->submodel[0].itrns+1;
    }
#endif
    */


    //Not sure what this is for 
    /*
    SVECT mn_coord = { 0.0, 0.0, 0.0 };
    SVECT im_coord = { 0.0, 0.0, 0.0 };

#ifdef _ADH_GROUNDWATER
    char prn_head[14][5] = { "HYD", "TRN", "BLT", "SLT", "GRD", "HVEL", "WVEL", "SALT", "CLAY", "SAND", "TRUE", "DIFF", "NS", "GW" };
#else
    char prn_head[13][5] = { "HYD", "TRN", "BLT", "SLT", "GRD", "HVEL", "WVEL", "SALT", "CLAY", "SAND", "TRUE", "DIFF", "NS" };
#endif
    */

//what is our mpi struct going to look like, leave for now
#ifdef _MESSG
    myid = smpi->myid;
    npes = smpi->npes;
#endif
    
#ifdef _DEBUG
    if (DEBUG_PICKETS == ON) {
        tl_check_all_pickets(__FILE__, __LINE__);
    }
#endif
    
    /* sets my_ndof */
    //includes ghosts
    //why? this is confusing
    //dont see this being used so comment this out
    //my_ndof = sm->ndofs;
    
    /* initial setup */
    //(*init_fnctn) (sm,isuperModel);
    //split up into 2 steps
    //printf("attempting initial\n");
    initialize_system(sm);
    //need to think about this bit
    initialize_dirichlet_bc(sm);


    //it = 0;
    //(*update_fnctn) (sm,isuperModel);
    it=0;
    //updates any other quantities FALSET in solution? Needed for solution with MPI too?
    update_function(sm);
    //printf("attempting update\n");
    //*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // get initial residual
    //(*residual_fnctn) (sm,isuperModel);
    //elemental loops are now part of call
    assemble_residual(sm,sm->grid);
//    for (int i=0;i<sm->ndofs;i++){
//        printf("resid before assembly[%d] = %f\n",i,sm->residual[i]);
//    }
    // Should we keep as is or calculate on fly?
    // this is for surface water calculations to ignore dry nodes

    // this is only on 2d top surface
    //Mark, ignore for now, need to reincorporate later
//    if (is_2d_hydro == FALSE && include_node_norm != NULL) {
//        //Mark to do, clean up once structs are more clear
//        include_node_in_residual(grid->nnodes, mod->sw->d2->head, mod->sw->d2->old_head, mod->fmap, include_node_norm);
//    }



    //Old call
    /*
    get_residual_norms(my_nnodes, nnodes, macro_nnodes, sm->residual, sm->sol,
#ifdef _PETSC
            sm->max_nsys,
#endif
            sm->nsys, &resid_max_norm, &resid_l2_norm, &inc_max_norm, &imax_node, &iinc_node, include_node_norm
#ifdef _MESSG
                       ,supersmpi
#endif
                       );
    */
    //New call, returns dof # instead of node #. Can back out nodes later if desired
    get_residual_norms(&resid_max_norm, &resid_l2_norm, &inc_max_norm,
        &imax_dof, &iinc_dof, include_dof,
        my_ndofs, ndofs, macro_ndofs, residual, dsol, sm->bc_mask
        );

#ifdef _MESSG
    resid_max_norm = messg_dmax(resid_max_norm, smpi->ADH_COMM);
#endif
    
#ifdef _DEBUG
    if (DEBUG_FULL)
        printScreen_dble_array("residual before solving", residual, ndofs, __LINE__, __FILE__);
    if (DEBUG_FULL || DEBUG_INIT_RESID || DEBUG_STOP_AFTER_INIT_RESID) {
        printf("\n **Initial resid_max_norm = %18.9e \n", resid_max_norm);
        if (DEBUG_STOP_AFTER_INIT_RESID) tl_error("Stopping after initial residual");
    }
#endif
    
#ifdef _DEBUG
    if (DEBUG_FULL) {
#ifdef _MESSG
        tag(smpi->ADH_COMM);
#else
        tag();
#endif
    }
#endif
    
    //printf("\n **Initial resid_max_norm = %18.9e \n", resid_max_norm);
    char prn_head[14][5] = { "HYD", "TRN", "BLT", "SLT", "GRD", "HVEL", "WVEL", "SALT", "CLAY", "SAND", "TRUE", "DIFF", "NS", "GW" };
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // loops over the nonlinear iterations
    do {
        //Mark, talk with Corey about this solver info stuff, for now just say it is false?
        //if (sm->solver_info.LINEAR_PROBLEM ==TRUE) {
        //    if(myid==0)printf("\n%s_%d TIME: %7.5e DT: %7.5e Progress: %3.2f%% | NIT: %2d | ", prn_head[0],index, sm->t_prev, sm->dt, (100.0 * (sm->t_prev - sm->t_init) / (sm->t_final - sm->t_init)), it + 1);} /* gkc warning come back later and fix use of submodel[0]. */
        //else{
            if(myid==0){
                printf("\n%s TIME: %7.5e DT: %7.5e Progress: %3.2f%% | NIT: %2d | ", prn_head[0], *(sm->t_prev), *(sm->dt), (100.0 * (*(sm->t_prev) - *(sm->t_init)) / (*(sm->t_final) - *(sm->t_init))), it + 1);
            }
        //}
        
        
        sm->nonlinear_it_total++;
        //printf("nonlinear it: %d\n",sm->nonlinear_it_total);

        //Old code. dont think  this applies anymore
        /*
        if (init_fnctn == fe_hvel_init || init_fnctn == fe_sw_hybrid_init) {
            sm->nonlinear_it_total_hvel++;
            sm->nonlinear_it_total++;
        }  else if (init_fnctn == fe_wvel_init) {
            sm->nonlinear_it_total_wvel++;
            sm->nonlinear_it_total++;
        } else if (init_fnctn == fe_sw2_init) {
            sm->nonlinear_it_total++;
        } else if (init_fnctn == fe_diffusive_init) {
            sm->nonlinear_it_total++;
        }
#ifdef _ADH_GROUNDWATER	
        else if (init_fnctn == fe_gw_init) {
            sm->nonlinear_it_total++;
        }	
#endif
        */

#ifdef _DEBUG
        /* forms the matrix */
        if (DEBUG_FULL) printf("\nLoading from fe_newton before another newton iteration ...\n");
#endif
        //loads global sparse system of equations
        //(*load_fnctn) (sm,isuperModel);
        assemble_jacobian(sm);


        
//        double temp1 = l_infty_norm(sm->nnz_diag, sm->vals_diag);
//        double temp2 = l2_norm(sm->vals_diag, sm->nnz_diag);
//        double temp3 = l2_norm(sm->residual,sm->ndofs);
//        double temp4 = l_infty_norm(sm->ndofs,sm->residual);
//        printf("CSR norms: %.17e, %.17e, %.17e, %.17e\n",temp1,temp2,temp3,temp4);
        //Screen_print_CSR(sm->indptr_diag, sm->cols_diag, sm->vals_diag, sm->ndofs);
        //check diagonal for near 0 entries and add to dirichlet data, ensure 0 increment
        check_diag(sm);
        apply_Dirichlet_BC(sm);

        //sarray_printScreen_int(sm->bc_mask, ndofs,"bc mask");
        //printf("Dirichlet applied\n");
        //lin_sys_CSR_printScreen(lin_sys);
//        get_residual_norms(&resid_max_norm, &resid_l2_norm, &inc_max_norm,
//        &imax_dof, &iinc_dof, include_dof,
//        sm->my_ndofs, sm->ndofs, sm->macro_ndofs, sm->residual, sm->dsol, sm->bc_mask
//        );
//        printf("residual norms computed after Dirichlet: %.17e, %.17e, %.17e\n",resid_max_norm,resid_l2_norm,inc_max_norm);
//        temp1 = l_infty_norm(sm->nnz_diag, sm->vals_diag);
//        temp2 = l2_norm(sm->vals_diag, sm->nnz_diag);
//        temp3 = l2_norm(sm->residual,sm->ndofs);
//        temp4 = l_infty_norm(sm->ndofs,sm->residual);
//        printf("CSR norms after Dirichlet: %.17e, %.17e, %.17e, %.17e\n",temp1,temp2,temp3,temp4);
//        printf("RHS\n");
//        for (int i=0;i<sm->ndofs;i++){
//            printf("Before solve resid[%d] = %f\n",i,sm->residual[i]);
//        }
        /* Set initial guess */
        //maybe redundant?
        sarray_init_dbl(dsol,ndofs);
        

#ifdef _DEBUG
        if (DEBUG_FULL) {
#ifdef _MESSG
            for(proc_count=0;proc_count<supersmpi->npes;proc_count++){
                if(proc_count==supersmpi->myid){
                    printf("*********** myid %d nnodes: %d nsys: nsys %d :: ",supersmpi->myid,nnodes,nsys);
#endif
                    printf("BEFORE SOLVE!\n");
                    printScreen_dble_array("sol before solving", sm->sol, ndofs, __LINE__, __FILE__);
                    printScreen_dble_array("residual before solving", residual, ndofs, __LINE__, __FILE__);
                    //printScreen_int_array("bc_mask before solving",sm->bc_mask, nnodes * nsys, __LINE__, __FILE__);
                    //printScreen_dble_array("scale_vect before solving",sm->scale_vect, nnodes * nsys, __LINE__, __FILE__);
#ifdef _MESSG
                }
            }
#endif
        }
#endif


        //Set up system + solve
        //Either use proprietary solver
        //or pipe to PETSC 

        //PETSC option
#ifdef _PETSC
        //do linear scaling, idk if it makes difference?
        scale_linear_system(lin_sys->indptr_diag, lin_sys->cols_diag, lin_sys->vals_diag, lin_sys->indptr_off_diag, lin_sys->cols_off_diag,
        lin_sys->vals_off_diag, lin_sys->residual, lin_sys->dsol, lin_sys->scale_vect, *(lin_sys->local_size), *(lin_sys->size), myid, lin_sys->ghosts, lin_sys->nghost);
        
        //store temporary vector of inital guess for linear system
        sarray_copy_dbl(x0, dsol,  ndofs);

        //printf("Using PETSC solver\n");
        // Solver
        //be sure values get updated since we used set with split arrays ...
        KSPSetOperators(lin_sys->ksp,lin_sys->A,lin_sys->A);
        
        //Want user to specify this but hard code for now
        //precon and stuff
        KSPSetType(lin_sys->ksp,PETSC_KSP); //GMRES iteratice solver
        //KSPSetType(lin_sys->ksp,KSPPREONLY); //Direct solve
        PC          pc;      /* preconditioner context */
        KSPGetPC(lin_sys->ksp, &pc);
        PCSetType(pc,PETSC_PRECON);//ILU only works in serial
        //some other options, seems to break stuff?
        //PetscOptionsSetValue(NULL, "-pc_factor_levels", "0");
        //PetscOptionsSetValue(NULL, "-ksp_gmres_restart", "400");
        
        //PetscCall(KSPSetTolerances(sm->ksp, 1.e-11, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
        // TODO: can I call this once in fe_main instead?
        // Does the matrix operator need to be specified before SetFromOptions each time?
        //KSPSetFromOptions(lin_sys->ksp);
        // TODO: I don't think sol needs to be initialized to zero but should double-check
        // set resid and solution into PETSdc objects should be fine on initialization
        // what about ghost values?
        //solution goes in X
        //ierr = MatView(sm->A, PETSC_VIEWER_STDOUT_WORLD);
        //ierr = VecView(sm->B, PETSC_VIEWER_STDOUT_WORLD);
        //HARD CODED FOR NOW, NEED TO CHANGE
        KSPSetTolerances(lin_sys->ksp, PETSC_RTOL, PETSC_ATOL, PETSC_DTOL, PETSC_MAXIT);
        //KSPSetTolerances(lin_sys->ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
        //KSPSetInitialGuessNonzero(lin_sys->ksp,PETSC_TRUE);
        //KSPSetFromOptions(lin_sys->ksp);
        //PCSetFromOptions(pc);
        //PCSetUp(pc);
        //KSPSetUp(lin_sys->ksp);
        //KSPView(lin_sys->ksp, PETSC_VIEWER_DEFAULT );
        //solve
        KSPSolve(lin_sys->ksp,lin_sys->B,lin_sys->X);
        //scatter forward appears to update array as we need
        //forward sends owned dofs -> ghosts
        //VecGhostUpdateBegin(sm->X,INSERT_VALUES,SCATTER_FORWARD);
        //VecGhostUpdateEnd(sm->X,INSERT_VALUES,SCATTER_FORWARD);
        KSPGetIterationNumber(lin_sys->ksp, &its);
        //printf("KSP iterations: %i\n",its);
        //unscale, idk if helps or not
        unscale_linear_system(dsol,x0,lin_sys->scale_vect,*(lin_sys->local_size));
        //printf("SOLVER STATUS   %d\n\n",status);
//    for(int j =0;j<ndofs;j++){
//            printf("dsol[%d] = %.17e \n",j,lin_sys->dsol[j]);
//        }
        //solv_flag = TRUE; // TODO: Does this need to change - SAM
        status = 0;
       
#ifdef _DEBUG
        // Viewer stuff
        if(DEBUG_MATRIX){
            // View Matrix
            PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Amat.m",&(viewer));
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
            MatView(lin_sys->A,viewer);
            PetscViewerPopFormat(viewer);
        }
        if(DEBUG_FULL){
            // View sol
            PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Svec.m",&(viewer));
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
            VecView(lin_sys->X,viewer);
            PetscViewerPopFormat(viewer);
            // View residual
            PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Rvec.m",&(viewer));
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
            VecView(lin_sys->B,viewer);
            PetscViewerPopFormat(viewer);
            // Solver
            PetscViewerASCIIOpen(PETSC_COMM_WORLD,"KSP.m",&(viewer));
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
            KSPView(lin_sys->ksp,viewer);
            PetscViewerPopFormat(viewer);
        }
#endif
#else        
        //Screen_print_CSR(sm->indptr_diag, sm->cols_diag, sm->vals_diag, sm->ndofs);
        //Non-PETSc option
        //does a scaling 

        //Mark commenting
        scale_linear_system(lin_sys->indptr_diag, lin_sys->cols_diag, lin_sys->vals_diag, lin_sys->indptr_off_diag, lin_sys->cols_off_diag,
        lin_sys->vals_off_diag, lin_sys->residual, lin_sys->dsol, lin_sys->scale_vect, *(lin_sys->local_size), *(lin_sys->size), myid, lin_sys->ghosts, lin_sys->nghost);
        //printf("linear system scaled\n");

//        temp1 = l_infty_norm(sm->nnz_diag, sm->vals_diag);
//        temp2 = l2_norm(sm->vals_diag, sm->nnz_diag);
//        temp3 = l2_norm(sm->residual,sm->ndofs);
//        temp4 = l_infty_norm(sm->ndofs,sm->residual);
//        printf("CSR norms after scaling: %.17e, %.17e, %.17e, %.17e\n",temp1,temp2,temp3,temp4);
        //factor the matrix preconditioner
        status = prep_umfpack(lin_sys->indptr_diag,lin_sys->cols_diag,lin_sys->vals_diag, *(lin_sys->local_size));
        
        

        //print before
        //printf("Just before solve\n");
        //Screen_print_CSR(sm->indptr_diag, sm->cols_diag, sm->vals_diag, sm->ndofs);
       // printf("rhs\n");
       // for(int j =0;j<sm->ndofs;j++){
       //     printf("bc mask [%d] = %d, resid[%d] = %.17e \n",j,sm->bc_mask[j],j,sm->residual[j]);
       // }
        //direct solve, need to turn of scaling if you want to try this by itself
        //status = solve_umfpack(sm->dsol, sm->indptr_diag, sm->cols_diag, sm->vals_diag, sm->residual, sm->local_size);
        



        //actually solve linear system, returns solution
        status = solve_linear_sys_bcgstab(dsol, lin_sys->indptr_diag, lin_sys->cols_diag, lin_sys->vals_diag, lin_sys->indptr_off_diag, lin_sys->cols_off_diag,
        lin_sys->vals_off_diag, residual, lin_sys->scale_vect, *(lin_sys->local_size), *(lin_sys->size) ,myid, lin_sys->ghosts, lin_sys->nghost);
        //printf("BCGSTAB completed\n");
        //printf("SOLVER STATUS   %d\n\n",status);
//        for(int j =0;j<ndofs;j++){
//            printf("dsol[%d] = %.17e \n",j,dsol[j]);
//        }
//        for (int i=0;i<ndofs;i++){
//          printf("Dsol increment[%d] = %f\n",i,dsol[i]);
//        }
#endif        
        
        /* adds the increment to the solution */
        //(*inc_fnctn) (sm,isuperModel);
        increment_function(sm);
        //printf("Solution incremented\n");

        //Mark needs to check, why is this
        //flip sign after update?
        //while not steady state??
//        if (sm.flag.STEADY_STATE == OFF) {
            //Commenting for now/
            for (i = 0, iend = ndofs; i < iend; i++){
                dsol[i] = -dsol[i];
            }

//        }
        
        /* calculates the residual after the increment?*/
        old_resid_norm = resid_l2_norm;
        update_function(sm);
        //for (int i=0;i<sm->ndofs;i++){
        //    printf("solution after increment, before new residual call[%d] = %f\n",i,sm->sol[i]);
        //}
        assemble_residual(sm,grid);
//                for (int i=0;i<sm->ndofs;i++){
//        printf("solution after increment, after new residual call[%d] = %f\n",i,sm->sol[i]);
//        }
//        for (int i=0;i<sm->ndofs;i++){
//        printf("residual after increment[%d] = %f\n",i,sm->residual[i]);
//        }
        //(*update_fnctn) (sm,isuperModel);
        //(*residual_fnctn) (sm,isuperModel);
        
#ifdef _DEBUG
#ifndef _PETSC
        if (DEBUG_FULL) {
#ifdef _MESSG
            for(proc_count=0;proc_count<supersmpi->npes;proc_count++){
                if(proc_count==supersmpi->myid){
                    printf("***********myid %d ",supersmpi->myid);
#endif
                    printf("AFTER SOLVE!\n");
                    printScreen_dble_array("sol after solving", sm->sol, ndofs, __LINE__, __FILE__);
                    printScreen_dble_array("residual before solving", residual, ndofs, __LINE__, __FILE__);
                    //printScreen_int_array("bc_mask after solving",sm->bc_mask, nnodes * nsys, __LINE__, __FILE__);
                    //printScreen_dble_array("scale_vect after solving",sm->scale_vect, nnodes * nsys, __LINE__, __FILE__);
                    //if (DEBUG_MATRIX) printScreen_matrix("matrix after solving", sm->diagonal, sm->matrix, nnodes, sm->max_nsys_sq,__LINE__, __FILE__);
#ifdef _MESSG
                }
            }
#endif
        }
#endif
#endif
        
        /* calculates the norms of the residual */
        //ignore for now
//        if (is_2d_hydro && include_node_norm != NULL) {
//            //Mark todo
//            include_node_in_residual(mod->grid->nnodes, mod->sw->d2->head, mod->sw->d2->old_head, mod->fmap, include_node_norm);
//        }
        //old call
//        get_residual_norms(my_nnodes, nnodes, macro_nnodes, sm->residual, sm->sol,
//#ifdef _PETSC
//                sm->max_nsys,
//#endif
//                sm->nsys, &resid_max_norm, &resid_l2_norm, &inc_max_norm, &imax_node, &iinc_node, include_node_norm
//#ifdef _MESSG
//                       ,supersmpi
//#endif
//                              );
        //new call
        get_residual_norms(&resid_max_norm, &resid_l2_norm, &inc_max_norm,
        &imax_dof, &iinc_dof, include_dof,
        my_ndofs, ndofs, macro_ndofs, residual, dsol, sm->bc_mask
        );
        //printf("residual norms computed: %.17e, %.17e, %.17e\n",resid_max_norm,resid_l2_norm,inc_max_norm);
#ifdef _DEBUG
        if (DEBUG_FULL) {
            printf("\n pe: %d resid_max_norm before line search = %18.9e my_nnodes: %d nnodes: %d macro_nnodes: %d \n", myid,resid_max_norm,my_nnodes, nnodes, macro_nnodes);
            tl_error("wha");
        }
#endif

        if (solv_isnan(resid_l2_norm) || solv_isinf(resid_l2_norm)) {
            solv_flag = NO;
            //what is solve_atf
//            for (i=0; i<sm->nsubmodels; i++) {
//                if(sm->submodel[i].proc_flag==1){
//                    sm->submodel[i].flag.SOLVE_ATF = FALSE;
//                }
//            }
        }
//Mark, do we really need UMFAIL too? Commenting out for now       
#ifdef _MESSG
        //UMFail_max = messg_imax(solver->UMFail, supersmpi->ADH_COMM);
        solv_flag_min = messg_imin(solv_flag, supersmpi->ADH_COMM);
#else
        //UMFail_max = messg_imax(solver->UMFail);
        solv_flag_min = messg_imin(solv_flag);
#endif     
        //Mark todo, what is this?
//        if (solver->PRN_NEWTON != OFF) {
//            // CJT :: right now, no way to tell which model this node resides on (need an inverse fmap!
//            // Gajanan gkc :: 2020.04.16 :: Cannot have a reliable inverse fmap since it is not invertible (multiple nodes may be mapped to a single row).
//            //                              Therefore, the following lins can cause a segmentation fault in coupled runs!
//            if (sm->nsubmodels == 1) {
//                mn_node    = sm->submodel[0].grid->node[imax_node].original_id;
//                mn_coord.x = sm->submodel[0].grid->node[imax_node].x;
//                mn_coord.y = sm->submodel[0].grid->node[imax_node].y;
//                mn_coord.z = sm->submodel[0].grid->node[imax_node].z;
//                
//                im_node    = sm->submodel[0].grid->node[iinc_node].original_id;
//                im_coord.x = sm->submodel[0].grid->node[iinc_node].x;
//                im_coord.y = sm->submodel[0].grid->node[iinc_node].y;
//                im_coord.z = sm->submodel[0].grid->node[iinc_node].z;
//            }
//#ifdef _MESSG
//            //MPI_Allreduce(&nnodes, &tot_nnode, 1, MPI_INT, MPI_SUM, sm->supersmpi->ADH_COMM);
//            tot_nnode = macro_nnodes;
//            messg_max_norm_loc(&inc_max_norm, &im_node, &im_coord, supersmpi->ADH_COMM, supersmpi->myid);
//            messg_max_norm_loc(&resid_max_norm, &mn_node, &mn_coord, supersmpi->ADH_COMM, supersmpi->myid);
//#else
//            tot_nnode = macro_nnodes;
//#endif
//        } else {
#ifdef _MESSG
            inc_max_norm = messg_dmax(inc_max_norm, supersmpi->ADH_COMM);
            resid_max_norm = messg_dmax(resid_max_norm, supersmpi->ADH_COMM);
#endif
//        }
        
        /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /* Perform A Line Search For Non-Linear Step If Initial (Full) Step Fails to: */
        /* (1) Infty Norm of Non-Linear Residual Does Not Satisfy Absolute Tolerance */
        /* (2) No Simple Decrease in l2 norm of Residual */
        /* (3) Not Enough Line Search Cuts Yet */
        linesearch_cuts = 0;

        while ((resid_max_norm > sm->tol_nonlin) && (resid_l2_norm > old_resid_norm) && (linesearch_cuts < sm->max_nonlin_linesearch_cuts)) {
            //printf("IN LINE SEARCH\n");
#ifdef _DEBUG
            if (DEBUG_FULL) {
                if (myid <= 0 && linesearch_cuts == 0) {
                    printf("\n\n IN LINE SEARCH: tol_nonlin: %30.20f \t\t resid_max_norm: %30.20f \t\t resid_l2_norm: %30.20f \t\t old_resid_norm: %30.20f", sm->tol_nonlin, resid_max_norm, resid_l2_norm, old_resid_norm);
                }
            }
#endif
            
            /* reduces the increment */
            for (i = 0, iend = ndofs; i < iend; i++){
                dsol[i] *= 0.5;
            }

            /* calculates the residual */
            increment_function(sm);
            update_function(sm);
            assemble_residual(sm,grid);

            //(*inc_fnctn) (sm,isuperModel);
            //(*update_fnctn) (sm,isuperModel);
            //(*residual_fnctn) (sm,isuperModel);
            
            /* calculates the norms of the residual */
//            if (is_2d_hydro && include_node_norm != NULL) {
//                //Mark, need to fix this once structs are finalized
//                include_node_in_residual(mod->grid->nnodes, mod->sw->d2->head, mod->sw->d2->old_head, mod->fmap, include_node_norm);
//            }
            //old call
//            get_residual_norms(my_nnodes, nnodes, macro_nnodes, sm->residual, sm->sol,
//#ifdef _PETSC
//                    sm->max_nsys,
//#endif
//                    sm->nsys, &resid_max_norm, &resid_l2_norm, &inc_max_norm, &imax_node, &iinc_node, include_node_norm
//#ifdef _MESSG
//                               ,supersmpi
//#endif
//                               );
            //new call
            get_residual_norms(&resid_max_norm, &resid_l2_norm, &inc_max_norm,
            &imax_dof, &iinc_dof, include_dof,
            my_ndofs, ndofs, macro_ndofs, residual, dsol, sm->bc_mask
            );
            
#ifdef _MESSG
            resid_max_norm = messg_dmax(resid_max_norm, supersmpi->ADH_COMM);
            inc_max_norm = messg_dmax(inc_max_norm, supersmpi->ADH_COMM);
#endif
            /* increments the line search counter */
            linesearch_cuts++;
            
#ifdef _DEBUG
            if (DEBUG_FULL) {
                printf("pe %d of %d :: new norms :: resid_max_norm: %30.20f \t resid_l2_norm: %30.20f \t old_resid_norm: %30.20f \t linesearch_cuts: %d \t max_nonlin_linesearch_cuts: %d \n",myid,npes,resid_max_norm,resid_l2_norm,old_resid_norm,linesearch_cuts,sm->max_nonlin_linesearch_cuts);
            }
#endif
        }
        /* End Line Search */
        /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
        if(myid==0) {
            //Mark, need to figure out screen outpu later
            //if (screen_output.residuals) {
                printf(" MAX_NL: %6.4e MAX_INC: %6.4e",resid_max_norm,inc_max_norm);
            //}
            //if (screen_output.worse_node_nonlinear) {
                //printf(" NTL_FALSEDE: %d %8.6e %8.6e %8.4e", mn_node + 1, mn_coord.x, mn_coord.y, mn_coord.z);
            //}
            //if (screen_output.worse_node_linear) {
                //printf(" ITL_FALSEDE: %d %8.6e %8.6e %8.4e", im_node + 1, im_coord.x, im_coord.y, im_coord.z);
//#ifdef _ADH_GROUNDWATER
                /*mwf debug */
                //printf(" --- inc[%d]= %g val= %g val prev= %g --- \n",im_node+1,-sm->sol[im_node],
                //       sm->submodel[0].sgw->gw_phead[im_node],sm->submodel[0].sgw->gw_phead[im_node]+sm->sol[im_node]);
//#endif
            //}
            printf(" NNODES: %2d",grid->nnodes);
            if (linesearch_cuts > 0)  printf(" LS");
            printf("\n");
        }
        
        /* increment the iteration counter */
        it++;
        
	#ifdef _DEBUG
        if (DEBUG_FULL) {
#ifdef _MESSG
            for(proc_count=0;proc_count<supersmpi->npes;proc_count++){
                if(proc_count==supersmpi->myid){
                    printf("***********myid %d ",supersmpi->myid);
#endif
                    printf("\nsolver variables:\n");
                    printf("inc_max_norm: %30.20f \t inc_nonlin: %30.20f \n",inc_max_norm,sm->inc_nonlin);
                    printf("resid_max_norm: %30.20f \t tol_nonlin: %30.20f \n",resid_max_norm,sm->tol_nonlin);
                    printf("it: %d \t sm->max_nonlin_it: %d\n",it, sm->max_nonlin_it);
                    printf("solv_flag_min: %d UMFail_max: %d\n",solv_flag_min,UMFail_max);
                    printf("LINEAR_PROBLEM: %d \n",sm->LINEAR_PROBLEM);
                    printf("inc_max_norm: %30.20f \t inc_nonlin: %30.20f \n",inc_max_norm,sm->inc_nonlin);
#ifdef _MESSG
                }
            }
#endif
        }
	#endif
        
        /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        // CONVERGENCE CHECK
        check = YES;
        if (solv_flag_min == YES) {
            check *= YES;
        } else {
            check = OFF;
        }
        //printf("CHECK VAL %d\n",check);
        if (it < sm->max_nonlin_it) {
            check *= YES;
        } else {
            check = OFF;
        }
        //printf("CHECK VAL %d\n",check);
        if (resid_max_norm < 1.0E+20) {
            check *= YES;
        } else check = OFF;
        //printf("CHECK VAL %d\n",check);
        if (UMFail_max == NO) {
            check *= YES;
        } else {
	       //printf("fe_newton UMFail_max != FALSE\n");
	       check = OFF;
        }
        //printf("CHECK VAL %d\n",check);
        //Mark, need to figure this out
        //if (solver->LINEAR_PROBLEM == TRUE) {
        //    if (sm->tol_nonlin > SMALL && resid_max_norm > sm->tol_nonlin) {
        //        check *= TRUE;
        //    } else check = OFF;
        //} else {
        if (resid_max_norm > sm->tol_nonlin && inc_max_norm > sm->inc_nonlin) {
                check *= YES;
        } else {
            check = OFF;
        }
        //printf("CHECK VAL %d\n",check);
        //check=OFF;


#ifdef _MESSG  // make sure that if one subdomain fails, they all do (cjt)
        MPI_Allreduce(&check, &keep_chugging, 1, MPI_INT, MPI_MAX, supersmpi->ADH_COMM);
#else
        keep_chugging = check;
#endif
    //printf("CHECK VAL %d\n",check);
    //printf("CHECK VAL %d\n",check);
    }while (keep_chugging);

#ifdef _PETSC
#ifdef _DEBUG
    // Destroy viewer
    if(DEBUG_MATRIX || DEBUG_FULL){
        ierr = PetscViewerDestroy(&(viewer));
    }
#endif
#endif
    
    if (include_dof!= NULL) include_dof = (int *) tl_free(sizeof(int), ndofs, include_dof);
    
    //this will need to change, ask corey
    //for (i=0;i<sm->nsubmodels; i++){
//        if(sm->submodel[i].proc_flag==1){
//            sm->submodel[i].flag.SOLVE_ATF = FALSE;
//        }
//    }
    sm->it_count_nonlin += it;
    
    /* checks the residual - if failure then reduce the time step and try again */

    //Mark need to revisit failure criteria, modifying for now
    //original here
//    if (solv_flag_min == FALSE || UMFail_max == TRUE
//        || ((sm->LINEAR_PROBLEM == TRUE && (sm->tol_nonlin > SMALL && resid_max_norm > sm->tol_nonlin))
//            || (sm->LINEAR_PROBLEM == FALSE && ((resid_max_norm > sm->tol_nonlin) && (inc_max_norm > sm->inc_nonlin))))
//        || (resid_max_norm >= 1.0E+20)) 

//modified
    if (UMFail_max == YES
        || ((sm->LINEAR_PROBLEM == YES && (sm->tol_nonlin > SMALL && resid_max_norm > sm->tol_nonlin))
            || (sm->LINEAR_PROBLEM == NO && ((resid_max_norm > sm->tol_nonlin) && (inc_max_norm > sm->inc_nonlin))))
        || (resid_max_norm >= 1.0E+20)){

        
        //Mark commenting this for now
        //how should we handle this?
//        for (i=0;i<sm->nsubmodels; i++){
//            if(sm->submodel[i].proc_flag==1){
//                sm->submodel[i].flag.SOLVE_ATF = FALSE;
//                sm->submodel[i].flag.UMFPACK_FAIL = FALSE;
//          }
//        }
        
        /* we are forcing it to accept the result at the max nonlinear iteration */
        if (sm->force_nonlin_it == YES) return (TRUE);
        
        tc_scale(sm->dt
#ifdef _MESSG
                 , supersmpi->ADH_COMM
#endif
                 );
        //i dont think this exists anymore, check with Corey
//        for (i=0;i<sm->nsubmodels; i++){
//            sm->submodel[i].dt=sm->dt;
//            if(sm->submodel[i].proc_flag==1){
//                sm->submodel[i].flag.TIME_ADAPT_FAIL = TRUE;
//            }
//        }
        
        /* if failed to converge, then reset to old values */
        printf("Failed to converge, reinitiliaing!!\n");
        initialize_system(sm);
        //need to think about this bit
        initialize_dirichlet_bc(sm);
        //(*init_fnctn) (sm, isuperModel);
        
        if (myid <= 0 && solv_flag_min == YES)
#ifdef _DEBUG
            printf(" #");  // nonlinear failure
#else
        printf(" #");
#endif
        sm->it_count_nonlin_failed += it;
        return (FALSE);
    }

    //Ask Corey what this is for
    
    /* if it makes it to here then the calculations were good */
//    for (i=0;i<sm->nsubmodels; i++){
//        if(sm->submodel[i].proc_flag==1){
//            sm->submodel[i].flag.SOLVE_ATF = TRUE;
//        }
//    }
    return (TRUE);
}





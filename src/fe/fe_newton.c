/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_sw2_2d_elem_resid.c This file collections functions responsible for
 *          the Netwon solve.                                                               */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns a series of residual norms.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] resid_max_norm (double *) the max (L-infinity) residual error norm for local processor (myid)
 * @param[out] resid_l2_norm  (double *) the l2 residual error norm maxed over all processors
 * @param[out] inc_max_norm   (double *) the incremental max (L-infinity) norm for local processor (myid)
 * @param[out] imax_node      (int *) the node ID at which the residual max occurs
 * @param[out] iinc_node      (int *) the node ID at which the incremental max occurs
 * @param[in]  grid           (SGRID *) a pointer to a grid
 * @param[in]  residual       (double *) the Newton residual
 * @param[in]  sol            (double *) the Newton solution (array of solution increments)
 * @param[in]  nsys           (int) the number of equations in the system
 * @param[in]  include_node   (int *) an array of nodal flags to determine whether the nodes are including in the norm calculations
 *
 * \note CJT\::
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void get_residual_norms(int my_nnodes, int nnodes, int macro_nnodes
#ifdef _PETSC
        , Vec residual, Vec sol, int max_nsys
#else
        , double *residual, double *sol
#endif
        , int nsys, double *resid_max_norm, double *resid_l2_norm, double *inc_max_norm, int *imax_node, int *iinc_node, int *include_node
#ifdef _MESSG
                        , SMPI *supersmpi
#endif
) {
    
    int i, j, idof, include_node_flag;
    double lolh = 1.0;
    double partial_l2_norm = 0.;
    double partial_max_norm = 0.;
    double partial_inc_max_norm = 0.;
    int ierr = 0;
    
#ifdef _MESSG
    int myid = supersmpi->myid; /* gkc temp warning come back later and check. */
#else
    int myid = 0;
#endif
    
    *resid_max_norm = 0.;
    *resid_l2_norm = 0.;
    *inc_max_norm = 0.;
    *imax_node = 0;
    *iinc_node = 0;
    
#ifdef _PETSC
    const double *read_sol, *read_residual;

    // Get read-only pointer to PETSc vectors
    ierr = VecGetArrayRead(sol, &(read_sol));
    ierr = VecGetArrayRead(residual, &(read_residual));
#endif
    
    for (i = 0; i < my_nnodes; i++) {
        for (j = 0; j<nsys; j++) {
#ifdef _PETSC
            // Here, i gets multiplied by max_nsys because
            // PETSC Vecs (and Mats) are sized based on max_nsys
            idof = i*max_nsys + j;
#else
            idof = i*nsys + j;
#endif
            
            include_node_flag = 1;
            if (include_node != NULL) {
                if (include_node[i] != YES) include_node_flag = 0;
            }
           
#ifdef _PETSC
            // Use pointers to PETSc residual and sol
            if ( (fabs(read_sol[idof]) > partial_inc_max_norm) && include_node_flag == 1) {
                partial_inc_max_norm = fabs(read_sol[idof]);
                *iinc_node = i;
            }
            if (solv_isnan(read_residual[idof])){
                printf("\nERROR :: myid %d  || residual[%d] = %f || equation: %d\n", myid, idof, read_residual[idof], j);
                printf("ERROR :: fe_newton.c :: get_residual_norms :: residual[%d] = %20.10f\n",idof,read_residual[idof]);
                ierr = 1;
                exit(-1);
            }
            
            if ((fabs(read_residual[idof]) > partial_max_norm) && include_node_flag == 1) {
                partial_max_norm = fabs(read_residual[idof]);
                *imax_node = i;
            }
            
            partial_l2_norm += read_residual[idof] * read_residual[idof];
#else
            // Directly use residual and sol double arrays
            if ( (fabs(sol[idof]) > partial_inc_max_norm) && include_node_flag == 1) {
                partial_inc_max_norm = fabs(sol[idof]);
                *iinc_node = i;
            }
            if (solv_isnan(residual[idof])){
                printf("\nERROR :: myid %d  || residual[%d] = %f || equation: %d\n", myid, idof, residual[idof], j);
                printf("ERROR :: fe_newton.c :: get_residual_norms :: residual[%d] = %20.10f\n",idof,residual[idof]);
                ierr = 1;
                exit(-1);
            }
            
            if ((fabs(residual[idof]) > partial_max_norm) && include_node_flag == 1) {
                partial_max_norm = fabs(residual[idof]);
                *imax_node = i;
            }
            
            partial_l2_norm += residual[idof] * residual[idof];
#endif
        }
    }
    *resid_max_norm = partial_max_norm;
    
    // Clean up the read array pointers to PETSc vectors
#ifdef _PETSC
    ierr = VecRestoreArrayRead(sol, &(read_sol));
    ierr = VecRestoreArrayRead(residual, &(read_residual));
#endif
    
    // CJT ::  below should only include include_node == YES, but is shouldn't make a huge difference
#ifdef _MESSG
    *resid_l2_norm = sqrt(messg_dsum(partial_l2_norm, supersmpi->ADH_COMM) / (macro_nnodes * nsys));
#else
    *resid_l2_norm = sqrt(messg_dsum(partial_l2_norm) / (nnodes * nsys));
#endif
    *inc_max_norm = partial_inc_max_norm;
    if((ierr > 0) && (myid == 0)) printf("\n +++++++ WARNING NaN generated!!! ++++++++ \n");
}

void print_concentration(SSUPER_MODEL *sm, int isuperModel, int myid, int npes);
void print_concentration_file(SSUPER_MODEL *sm, int isuperModel, int myid, int npes);
void print_sw2_file(SSUPER_MODEL *sm, int isuperModel, int myid, int npes);
void print_sw3_file(SSUPER_MODEL *sm, int isuperModel, int myid, int npes);
void print_dw_file(SSUPER_MODEL *sm, int isuperModel, int myid, int npes);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns an array of flags determining whether to include the node in the error norm
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] include_node_norm (int *) an integer array of nodal flags
 * @param[in]  is_2d_hydro       (int)   a flag to indicate if this is a 2D SW model run
 * @param[in]  head              (double *) the current depth
 * @param[in]  old_head          (double *) the old depth
 *
 * \note CJT\:: Excludes nodes that are currently dry or were dry last time-step
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void include_node_in_residual(int nnodes, double *head, double *old_head, int *fmap, int *include_node_norm) {
    int i;
    sarray_init_value_int(include_node_norm, nnodes, YES);
    for (i=0; i<nnodes; i++) {
        if (MIN(head[i],old_head[i]) > 0) {
            include_node_norm[fmap[i]] = YES;
        } else {
            include_node_norm[fmap[i]] = NO;
        }
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     The FE Newton driver
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in]  mod             (SMODEL *) a pointer to an AdH model
 * @param[in]  grid            (SGRID *) a pointer to a grid in the AdH model
 * @param[in]  init_fnctn      (void *) a pointer to a function for initializing the AdH solutions
 * @param[in]  update_fnctn    (void *) a pointer to a function for updating the solution for parallel runs
 * @param[in]  residual_fnctn  (void *) a pointer to a function for calculating the nonlinear residual
 * @param[in]  load_fnctn      (void *) a pointer to a function for building the Jacobi matrix
 * @param[in]  inc_fnctn       (void *) a pointer to a function for incrementing the solution
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int fe_newton(SSUPER_MODEL *sm,                           /* input supermodel */
              int isuperModel,                            /* the supermodel id */
              int my_nnodes, int nnodes, int macro_nnodes,/* Resident, Resident + ghost, gathered resident */
#ifdef _MESSG
              SMPI *supersmpi,                            /* either sm->supersmpi or sm->supersmpi_wvel */
#endif
              void (*init_fnctn) (SSUPER_MODEL *, int),        /* a function to set the initial iterate */
              void (*update_fnctn) (SSUPER_MODEL *, int),      /* a function to update the solution for parallel runs */
              void (*residual_fnctn) (SSUPER_MODEL *, int),    /* a function to calculate the nonlinear residual */
              void (*load_fnctn) (SSUPER_MODEL *, int),        /* a function to load the matrix */
              void (*inc_fnctn) (SSUPER_MODEL *, int)          /* a function to increment the solution */
)
{
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG OPTIONS
    int DEBUG_FULL = OFF;       // prints everything but matrix
    int DEBUG_MATRIX = OFF;     // prints the matrix
    int DEBUG_PICKETS = OFF;    // check memory
    
    int DEBUG_INIT_RESID = OFF;             // print the initial residual only
    int DEBUG_STOP_AFTER_INIT_RESID = OFF;  // stop of this intitial residual print
    
    //if (init_fnctn == fe_sw_hybrid_wvel_init) DEBUG_FULL = ON;
    //if (init_fnctn == fe_transport_init) DEBUG_STOP_AFTER_INIT_RESID = ON;
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    
    SMODEL *mod;
    
    int check = NO;		/* newton solver check */
    int keep_chugging = NO;     /* newton solver check */
    double next_dt = 0;         /* temp variable for the time step for STEADY STATE */
    double ratio = 0;           /* STEADY STATE VARIABLES */
    int i, j, k;                /* loop counter */
    int idof;                   /* location of the maximum residual */
    int imax_node = 0;          /* location of the maximum residual */
    int iinc_node = 0;          /* location of the maximum ..increment? */
    int it;                     /* the nonlinear iteration counter */
    int linesearch_cuts;        /* the number of line searches */
    int solv_flag;              /* flag that tells if the linear solver converged */
    int iend;                   /* loop limit */
    int my_ndof;                /* the number of degrees of freedom I am responsible for */
    double resid_max_norm = 0.0;    /* the max norm of the residual */
    double resid_l2_norm = 0.0; /* the l2 norm of the residual */
    double old_resid_norm;      /* the previous norm of the residual - used for the line search */
    double partial_max_norm;    /* the partial residual norm prior to summing over processors */
    double partial_l2_norm;     /* the partial residual norm prior to summing over processors */
    double inc_max_norm = 0.0;  /* the max increment change */
    static double initial_residual; /* Initial_Residual */
    double spatial_residual_proc;   /* Individual Processor resid */
    double lolh;
    static int count_std;
    int mn_node = 0, im_node = 0;
    int tot_nnode;
    int UMFail_max = NO;
    int solv_flag_min = YES;
    int myid = 0, npes = 1, proc_count = 0;
#ifdef _PETSC
    int ierr = 0; // SAM HERE
    int its;
    MatInfo info;
#ifdef _DEBUG
    PetscViewer viewer;
#endif
#endif
    
#ifdef _DEBUG
    if (DEBUG_FULL) {
        for (i=0; i<sm->nsubmodels; i++){
            printf("fe_newton :: myid: %d ndim: %d :: my_nnodes: %d grid->my_nnodes: %d :: nnodes: %d grid->nnodes: %d  :: macro_nnodes: %d grid->macro_nnodes: %d\n",
                   sm->submodel[i].grid->smpi->myid,
                   sm->submodel[i].grid->ndim,
                   my_nnodes,sm->submodel[i].grid->my_nnodes,
                   nnodes,sm->submodel[i].grid->nnodes,
                   macro_nnodes,sm->submodel[i].grid->macro_nnodes);
        }
    }
#endif
    
    // Is there a 2D SW model inside the superModel?
    int is_2d_hydro = OFF;
    for (i=0; i<sm->nsubmodels; i++){
        if (sm->submodel[i].flag.DIFFUSIVE_WAVE==ON && sm->submodel[i].flag.GW_FLOW==ON){
            is_2d_hydro = OFF;
            break;
        }
        else if (sm->submodel[i].flag.SW2_FLOW || sm->submodel[i].flag.DIFFUSIVE_WAVE) {
            is_2d_hydro = ON;
            break;
        }
    }
    
    // this array determines what nodes will be including in the norm calculations
    int *include_node_norm = NULL;
    if (init_fnctn != fe_sw_hybrid_wvel_init) {
        include_node_norm = (int *) tl_alloc(sizeof(int), nnodes);
        sarray_init_value_int(include_node_norm, nnodes, YES);
    }
    
    Solver_Info *solver = &(sm->solver_info);
    solver->UMFail = NO;
    
    // aliases
    int nsys = sm->nsys;
    int nsys_sq = nsys * nsys;
    
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
    
    SVECT mn_coord = { 0.0, 0.0, 0.0 };
    SVECT im_coord = { 0.0, 0.0, 0.0 };

#ifdef _ADH_GROUNDWATER
    char prn_head[14][5] = { "HYD", "TRN", "BLT", "SLT", "GRD", "HVEL", "WVEL", "SALT", "CLAY", "SAND", "TEMP", "DIFF", "NS", "GW" };
#else
    char prn_head[13][5] = { "HYD", "TRN", "BLT", "SLT", "GRD", "HVEL", "WVEL", "SALT", "CLAY", "SAND", "TEMP", "DIFF", "NS" };
#endif
    
#ifdef _MESSG
    myid = supersmpi->myid;
    npes = supersmpi->npes;
#endif
    
#ifdef _DEBUG
    if (debug.newton == ON) {
        DEBUG_FULL = ON;
    }
    if ((init_fnctn == fe_hvel_init) && (debug.matrix_hvel == ON)) {
        DEBUG_MATRIX = ON;
    }
    if (DEBUG_PICKETS == ON) {
        tl_check_all_pickets(__FILE__, __LINE__);
    }
#endif
    
    /* sets my_ndof */
    my_ndof = nnodes * nsys;
    
    /* initial setup */
    (*init_fnctn) (sm,isuperModel);
    
    it = 0;
    (*update_fnctn) (sm,isuperModel);
    
    //*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // get initial residual
    (*residual_fnctn) (sm,isuperModel);
    
    // Do not require this if WVEL is being called, since no 2D is used
    if (is_2d_hydro == TRUE && include_node_norm != NULL) {
        for (i=0; i<sm->nsubmodels; i++) {
            if(sm->submodel[i].proc_flag==1 && sm->submodel[i].flag.SW2_FLOW){
                mod = &(sm->submodel[i]);
                include_node_in_residual(mod->grid->nnodes, mod->sw->d2->head, mod->sw->d2->old_head, mod->fmap, include_node_norm);
            }
        }
    }
    mod=NULL;
   
    get_residual_norms(my_nnodes, nnodes, macro_nnodes, sm->residual, sm->sol,
#ifdef _PETSC
            sm->max_nsys,
#endif
            sm->nsys, &resid_max_norm, &resid_l2_norm, &inc_max_norm, &imax_node, &iinc_node, include_node_norm
#ifdef _MESSG
                       ,supersmpi
#endif
                       );
#ifdef _MESSG
    resid_max_norm = messg_dmax(resid_max_norm, supersmpi->ADH_COMM);
#endif
    
#ifdef _DEBUG
#ifndef _PETSC
    if (DEBUG_FULL)
        printScreen_resid("residual before solving", sm->residual, nnodes, nsys, __LINE__, __FILE__);
    if (DEBUG_FULL || DEBUG_INIT_RESID || DEBUG_STOP_AFTER_INIT_RESID) {
        printf("\n **Initial resid_max_norm = %18.9e \n", resid_max_norm);
        if (DEBUG_STOP_AFTER_INIT_RESID) tl_error("Stopping after initial residual");
    }
#endif
#endif
    
#ifdef _DEBUG
    if (DEBUG_FULL) {
#ifdef _MESSG
        tag(supersmpi->ADH_COMM);
#else
        tag();
#endif
    }
#endif
    
    //printf("\n **Initial resid_max_norm = %18.9e \n", resid_max_norm);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // loops over the nonlinear iterations
    do {
        if (sm->solver_info.LINEAR_PROBLEM ==YES) {
            if(myid==0)printf("\n%s_%d TIME: %7.5e DT: %7.5e Progress: %3.2f%% | NIT: %2d | ", prn_head[solver->PRN_NEWTON - 1],index, sm->submodel[0].t_prev, sm->dt, (100.0 * (sm->submodel[0].t_prev - sm->submodel[0].t_init) / (sm->submodel[0].t_final - sm->submodel[0].t_init)), it + 1);} /* gkc warning come back later and fix use of submodel[0]. */
        else{
            if(myid==0)printf("\n%s TIME: %7.5e DT: %7.5e Progress: %3.2f%% | NIT: %2d | ", prn_head[solver->PRN_NEWTON - 1], sm->submodel[0].t_prev, sm->dt, (100.0 * (sm->submodel[0].t_prev - sm->submodel[0].t_init) / (sm->submodel[0].t_final - sm->submodel[0].t_init)), it + 1);
        }
        
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
        
#ifdef _DEBUG
        /* forms the matrix */
        if (DEBUG_FULL) printf("\nLoading from fe_newton before another newton iteration ...\n");
#endif
        (*load_fnctn) (sm,isuperModel);
#ifdef _PETSC
        // print preallocation statistics if using PETSc
        ierr = MatGetInfo(sm->A,MAT_LOCAL,&info);

        printf("Mallocs = %f\n",info.mallocs);
        printf("nz_allocated = %f\n",info.nz_allocated);
        printf("nz_used = %f\n",info.nz_used);
        printf("nz_unneeded = %f\n",info.nz_unneeded);
#endif
//        printf("myid: %d nsys: %d nnodes: %d\n",myid,nsys,nnodes);
//        tl_check_all_pickets(__FILE__,__LINE__);
//        MPI_Barrier(MPI_COMM_WORLD);
        
        /* Set initial guess */
#ifndef _PETSC
        solv_init_dbl(nnodes * nsys, sm->sol);
#endif
        
#ifdef _DEBUG
//        printf("myid: %d nsys: %d nnodes: %d\n",myid,nsys,nnodes);
//        tl_check_all_pickets(__FILE__,__LINE__);
//        MPI_Barrier(MPI_COMM_WORLD); tag(MPI_COMM_WORLD); tl_error("test");

	    /*mwf debug  
	    printScreen_matrix("matrix before solving", sm->diagonal, sm->matrix, nnodes, sm->max_nsys_sq,__LINE__, __FILE__);
	    printScreen_resid("residual before solving", sm->residual, grid->nnodes, nsys, __LINE__, __FILE__);
	    */
#endif

#ifdef _DEBUG
#ifndef _PETSC
        if (DEBUG_FULL) {
#ifdef _MESSG
            for(proc_count=0;proc_count<supersmpi->npes;proc_count++){
                if(proc_count==supersmpi->myid){
                    printf("*********** myid %d nnodes: %d nsys: nsys %d :: ",supersmpi->myid,nnodes,nsys);
#endif
                    printf("BEFORE SOLVE!\n");
                    printScreen_dble_array("sol before solving", sm->sol, nnodes * nsys, __LINE__, __FILE__);
                    printScreen_resid("residual before solving", sm->residual, nnodes, nsys, __LINE__, __FILE__);
                    //printScreen_int_array("bc_mask before solving",sm->bc_mask, nnodes * nsys, __LINE__, __FILE__);
                    //printScreen_dble_array("scale_vect before solving",sm->scale_vect, nnodes * nsys, __LINE__, __FILE__);
                    if (DEBUG_MATRIX) printScreen_matrix("matrix before solving", sm->diagonal, sm->matrix, nnodes, sm->max_nsys_sq,__LINE__, __FILE__);
#ifdef _MESSG
                }
            }
#endif
        }
#endif
#endif

#ifdef _PETSC
        // Solver
        KSPSetOperators(sm->ksp,sm->A,sm->A);
        // TODO: can I call this once in fe_main instead?
        // Does the matrix operator need to be specified before SetFromOptions each time?
        //KSPSetFromOptions(sm->ksp);
        // TODO: I don't think sol needs to be initialized to zero but should double-check
        KSPSolve(sm->ksp,sm->residual,sm->sol);
        KSPGetIterationNumber(sm->ksp, &its);
        printf("KSP iterations: %i\n",its);
        solv_flag = YES; // TODO: Does this need to change - SAM
       
#ifdef _DEBUG
        // Viewer stuff
        if(DEBUG_MATRIX){
            // View Matrix
            PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Amat.m",&(viewer));
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
            MatView(sm->A,viewer);
            PetscViewerPopFormat(viewer);
        }
        if(DEBUG_FULL){
            // View sol
            PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Svec.m",&(viewer));
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
            VecView(sm->sol,viewer);
            PetscViewerPopFormat(viewer);
            // View residual
            PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Rvec.m",&(viewer));
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
            VecView(sm->residual,viewer);
            PetscViewerPopFormat(viewer);
            // Solver
            PetscViewerASCIIOpen(PETSC_COMM_WORLD,"KSP.m",&(viewer));
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
            KSPView(sm->ksp,viewer);
            PetscViewerPopFormat(viewer);
        }
#endif
#else
        /* solves the matrix */
        solv_linear_sys_setup(solver, sm->bc_mask, sm->matrix, sm->diagonal, sm->residual, sm->sol, sm->scale_vect, my_nnodes, nnodes, nsys
#ifdef _MESSG
                              , supersmpi
#endif
                              );
        
        solv_flag = solv_linear_sys_solve(solver, sm->bc_mask, sm->matrix, sm->diagonal, sm->residual, sm->sol, sm->scale_vect, my_nnodes, nnodes, nsys
#ifdef _MESSG
                                          , supersmpi
#endif
                                          );
#endif
        
        /* adds the increment to the solution */
        (*inc_fnctn) (sm,isuperModel);
#ifdef _DEBUG
        if(init_fnctn == fe_transport_init){
            print_concentration_file(sm, isuperModel, myid, npes);
        }
        if(init_fnctn == fe_sw2_init){
            print_sw2_file(sm, isuperModel, myid, npes);
        }
        if(init_fnctn == fe_wvel_init){
            print_sw3_file(sm,isuperModel,myid,npes);
        }
        if(init_fnctn == fe_diffusive_init){
            printf("Printing diffusive wave solution to file\n");
            print_dw_file(sm,isuperModel,myid,npes);
        }
#endif
        
        if (sm->submodel[0].flag.STEADY_STATE == OFF) {
#ifdef _PETSC
            ierr = VecScale(sm->sol,-1.);
#else
            for (i = 0, iend = nnodes * nsys; i < iend; i++)
                sm->sol[i] = -sm->sol[i];
#endif
        }
        
        /* calculates the residual */
        old_resid_norm = resid_l2_norm;
        (*update_fnctn) (sm,isuperModel);
        (*residual_fnctn) (sm,isuperModel);
        
#ifdef _DEBUG
#ifndef _PETSC
        if (DEBUG_FULL) {
#ifdef _MESSG
            for(proc_count=0;proc_count<supersmpi->npes;proc_count++){
                if(proc_count==supersmpi->myid){
                    printf("***********myid %d ",supersmpi->myid);
#endif
                    printf("AFTER SOLVE!\n");
                    printScreen_dble_array("sol after solving", sm->sol, nnodes * nsys, __LINE__, __FILE__);
                    printScreen_resid("residual after solving", sm->residual, nnodes, nsys, __LINE__, __FILE__);
                    //printScreen_int_array("bc_mask after solving",sm->bc_mask, nnodes * nsys, __LINE__, __FILE__);
                    //printScreen_dble_array("scale_vect after solving",sm->scale_vect, nnodes * nsys, __LINE__, __FILE__);
                    if (DEBUG_MATRIX) printScreen_matrix("matrix after solving", sm->diagonal, sm->matrix, nnodes, sm->max_nsys_sq,__LINE__, __FILE__);
#ifdef _MESSG
                }
            }
#endif
        }
#endif
#endif
        
        /* calculates the norms of the residual */
        if (is_2d_hydro && include_node_norm != NULL) {
            for (i=0; i<sm->nsubmodels; i++) {          // gkc commenting out - need a work around for this
                if(sm->submodel[i].proc_flag==1 && sm->submodel[i].flag.SW2_FLOW){
                    mod = &(sm->submodel[i]);
                    include_node_in_residual(mod->grid->nnodes, mod->sw->d2->head, mod->sw->d2->old_head, mod->fmap, include_node_norm);
                }
            }
            mod=NULL;
        }
        get_residual_norms(my_nnodes, nnodes, macro_nnodes, sm->residual, sm->sol,
#ifdef _PETSC
                sm->max_nsys,
#endif
                sm->nsys, &resid_max_norm, &resid_l2_norm, &inc_max_norm, &imax_node, &iinc_node, include_node_norm
#ifdef _MESSG
                       ,supersmpi
#endif
                              );
#ifdef _DEBUG
        if (DEBUG_FULL) {
            printf("\n pe: %d resid_max_norm before line search = %18.9e my_nnodes: %d nnodes: %d macro_nnodes: %d \n", myid,resid_max_norm,my_nnodes, nnodes, macro_nnodes);
            tl_error("wha");
        }
#endif
        if (solv_isnan(resid_l2_norm) || solv_isinf(resid_l2_norm)) {
            solv_flag = NO;
            for (i=0; i<sm->nsubmodels; i++) {
                if(sm->submodel[i].proc_flag==1){
                    sm->submodel[i].flag.SOLVE_ATF = NO;
                }
            }
        }
        
#ifdef _MESSG
        UMFail_max = messg_imax(solver->UMFail, supersmpi->ADH_COMM);
        solv_flag_min = messg_imin(solv_flag, supersmpi->ADH_COMM);
#else
        UMFail_max = messg_imax(solver->UMFail);
        solv_flag_min = messg_imin(solv_flag);
#endif
        
        if (solver->PRN_NEWTON != OFF) {
            // CJT :: right now, no way to tell which model this node resides on (need an inverse fmap!
            // Gajanan gkc :: 2020.04.16 :: Cannot have a reliable inverse fmap since it is not invertible (multiple nodes may be mapped to a single row).
            //                              Therefore, the following lins can cause a segmentation fault in coupled runs!
            if (sm->nsubmodels == 1) {
                mn_node    = sm->submodel[0].grid->node[imax_node].original_id;
                mn_coord.x = sm->submodel[0].grid->node[imax_node].x;
                mn_coord.y = sm->submodel[0].grid->node[imax_node].y;
                mn_coord.z = sm->submodel[0].grid->node[imax_node].z;
                
                im_node    = sm->submodel[0].grid->node[iinc_node].original_id;
                im_coord.x = sm->submodel[0].grid->node[iinc_node].x;
                im_coord.y = sm->submodel[0].grid->node[iinc_node].y;
                im_coord.z = sm->submodel[0].grid->node[iinc_node].z;
            }
#ifdef _MESSG
            //MPI_Allreduce(&nnodes, &tot_nnode, 1, MPI_INT, MPI_SUM, sm->supersmpi->ADH_COMM);
            tot_nnode = macro_nnodes;
            messg_max_norm_loc(&inc_max_norm, &im_node, &im_coord, supersmpi->ADH_COMM, supersmpi->myid);
            messg_max_norm_loc(&resid_max_norm, &mn_node, &mn_coord, supersmpi->ADH_COMM, supersmpi->myid);
#else
            tot_nnode = macro_nnodes;
#endif
        } else {
#ifdef _MESSG
            inc_max_norm = messg_dmax(inc_max_norm, supersmpi->ADH_COMM);
            resid_max_norm = messg_dmax(resid_max_norm, supersmpi->ADH_COMM);
#endif
        }
        
        /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /* Perform A Line Search For Non-Linear Step If Initial (Full) Step Fails to: */
        /* (1) Infty Norm of Non-Linear Residual Does Not Satisfy Absolute Tolerance */
        /* (2) No Simple Decrease in l2 norm of Residual */
        /* (3) Not Enough Line Search Cuts Yet */
        linesearch_cuts = 0;
        while ((resid_max_norm > sm->tol_nonlin) && (resid_l2_norm > old_resid_norm) && (linesearch_cuts < solver->max_nonlin_linesearch_cuts)) {
            
#ifdef _DEBUG
            if (DEBUG_FULL) {
                if (myid <= 0 && linesearch_cuts == 0) {
                    printf("\n\n IN LINE SEARCH: tol_nonlin: %30.20f \t\t resid_max_norm: %30.20f \t\t resid_l2_norm: %30.20f \t\t old_resid_norm: %30.20f", sm->tol_nonlin, resid_max_norm, resid_l2_norm, old_resid_norm);
                }
            }
#endif
            
            /* reduces the increment */
#ifdef _PETSC
            ierr = VecScale(sm->sol,0.5);
#else
            for (i = 0, iend = nnodes * nsys; i < iend; i++)
                sm->sol[i] *= 0.5;
#endif
            
            /* calculates the residual */
            (*inc_fnctn) (sm,isuperModel);
            (*update_fnctn) (sm,isuperModel);
            (*residual_fnctn) (sm,isuperModel);
            
            /* calculates the norms of the residual */
            if (is_2d_hydro && include_node_norm != NULL) {
                for (i=0; i<sm->nsubmodels; i++) {
                    if(sm->submodel[i].proc_flag==1 && sm->submodel[i].flag.SW2_FLOW){
                        mod = &(sm->submodel[i]);
                        include_node_in_residual(mod->grid->nnodes, mod->sw->d2->head, mod->sw->d2->old_head, mod->fmap, include_node_norm);
                    }
                }
                mod=NULL;
            }
            get_residual_norms(my_nnodes, nnodes, macro_nnodes, sm->residual, sm->sol,
#ifdef _PETSC
                    sm->max_nsys,
#endif
                    sm->nsys, &resid_max_norm, &resid_l2_norm, &inc_max_norm, &imax_node, &iinc_node, include_node_norm
#ifdef _MESSG
                               ,supersmpi
#endif
                               );
            
#ifdef _MESSG
            resid_max_norm = messg_dmax(resid_max_norm, supersmpi->ADH_COMM);
            inc_max_norm = messg_dmax(inc_max_norm, supersmpi->ADH_COMM);
#endif
            /* increments the line search counter */
            linesearch_cuts++;
            
#ifdef _DEBUG
            if (DEBUG_FULL) {
                printf("pe %d of %d :: new norms :: resid_max_norm: %30.20f \t resid_l2_norm: %30.20f \t old_resid_norm: %30.20f \t linesearch_cuts: %d \t max_nonlin_linesearch_cuts: %d \n",myid,npes,resid_max_norm,resid_l2_norm,old_resid_norm,linesearch_cuts,solver->max_nonlin_linesearch_cuts);
            }
#endif
        }
        /* End Line Search */
        /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
        if(myid==0) {
            if (screen_output.residuals) {
                printf(" MAX_NL: %6.4e MAX_INC: %6.4e",resid_max_norm,inc_max_norm);
            }
            if (screen_output.worse_node_nonlinear) {
                printf(" NTL_NODE: %d %8.6e %8.6e %8.4e", mn_node + 1, mn_coord.x, mn_coord.y, mn_coord.z);
            }
            if (screen_output.worse_node_linear) {
                printf(" ITL_NODE: %d %8.6e %8.6e %8.4e", im_node + 1, im_coord.x, im_coord.y, im_coord.z);
#ifdef _ADH_GROUNDWATER
                /*mwf debug */
                //printf(" --- inc[%d]= %g val= %g val prev= %g --- \n",im_node+1,-sm->sol[im_node],
                //       sm->submodel[0].sgw->gw_phead[im_node],sm->submodel[0].sgw->gw_phead[im_node]+sm->sol[im_node]);
#endif
            }
            printf(" NNODES: %2d",tot_nnode);
            if (linesearch_cuts > 0)  printf(" LS");
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
                    printf("LINEAR_PROBLEM: %d \n",solver->LINEAR_PROBLEM);
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
        } else check = OFF;
        
        if (it < sm->max_nonlin_it) {
            check *= YES;
        } else check = OFF;
        
        if (resid_max_norm < 1.0E+20) {
            check *= YES;
        } else check = OFF;
        
        if (UMFail_max == NO) {
            check *= YES;
        } else {
	  /*mwf debug 
	  printf("fe_newton UMFail_max != NO\n");
	  */
	  check = OFF;
        }
        if (solver->LINEAR_PROBLEM == YES) {
            if (sm->tol_nonlin > SMALL && resid_max_norm > sm->tol_nonlin) {
                check *= YES;
            } else check = OFF;
        } else {
            if (resid_max_norm > sm->tol_nonlin && inc_max_norm > sm->inc_nonlin) {
                check *= YES;
            } else check = OFF;
        }
        
#ifdef _MESSG  // make sure that if one subdomain fails, they all do (cjt)
        MPI_Allreduce(&check, &keep_chugging, 1, MPI_INT, MPI_MAX, supersmpi->ADH_COMM);
#else
        keep_chugging = check;
#endif
        
    } while (keep_chugging);

#ifdef _PETSC
#ifdef _DEBUG
    // Destroy viewer
    if(DEBUG_MATRIX || DEBUG_FULL){
        ierr = PetscViewerDestroy(&(viewer));
    }
#endif
#endif
    
    if (include_node_norm != NULL) include_node_norm = (int *) tl_free(sizeof(int), nnodes, include_node_norm);
    
    for (i=0;i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            sm->submodel[i].flag.SOLVE_ATF = NO;
        }
    }
    solver->it_count_nonlin += it;
    
    /* checks the residual - if failure then reduce the time step and try again */
    if (solv_flag_min == NO || UMFail_max == YES
        || ((solver->LINEAR_PROBLEM == YES && (sm->tol_nonlin > SMALL && resid_max_norm > sm->tol_nonlin))
            || (solver->LINEAR_PROBLEM == NO && ((resid_max_norm > sm->tol_nonlin) && (inc_max_norm > sm->inc_nonlin))))
        || (resid_max_norm >= 1.0E+20)) {
        
        for (i=0;i<sm->nsubmodels; i++){
            if(sm->submodel[i].proc_flag==1){
                sm->submodel[i].flag.SOLVE_ATF = NO;
                sm->submodel[i].flag.UMFPACK_FAIL = NO;
            }
        }
        
        /* we are forcing it to accept the result at the max nonlinear iteration */
        if (solver->force_nonlin_it == YES) return (YES);
        
        tc_scale(&sm->dt
#ifdef _MESSG
                 , supersmpi->ADH_COMM
#endif
                 );
        
        for (i=0;i<sm->nsubmodels; i++){
            sm->submodel[i].dt=sm->dt;
            if(sm->submodel[i].proc_flag==1){
                sm->submodel[i].flag.TIME_ADAPT_FAIL = YES;
            }
        }
        
        /* if failed to converge, then reset to old values */
        (*init_fnctn) (sm, isuperModel);
        
        if (myid <= 0 && solv_flag_min == YES)
#ifdef _DEBUG
            printf(" #");  // nonlinear failure
#else
        printf(" #");
#endif
        solver->it_count_nonlin_failed += it;
        return (NO);
    }
    
    /* if it makes it to here then the calculations were good */
    for (i=0;i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            sm->submodel[i].flag.SOLVE_ATF = YES;
        }
    }
    
#ifdef _DEBUG
#ifndef _PETSC
    if (DEBUG_FULL){
        printScreen_dble_array("sol before leaving fe_newton_hybrid with YES", sm->sol, nnodes * nsys, __LINE__, __FILE__);
        printScreen_resid("residual before_leaving fe_newton_hybrid with YES", sm->residual, nnodes, nsys, __LINE__, __FILE__);
    }
#endif
#endif
    
    return (YES);
}

void print_concentration(SSUPER_MODEL *sm, int isuperModel, int myid, int npes){
    SMODEL *mod;
    int i, ierr;
    double *c;
    int ifmap, gid;

#ifdef _MESSG
#ifdef _PETSC
    ierr = MPI_Barrier(PETSC_COMM_WORLD);
#else
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    mod = &(sm->submodel[isuperModel]);
    
    if (mod->is_sediment_running) {
#ifdef _SEDIMENT
        c = mod->sed->susload[mod->ised].c;
#endif
    } else {
        c = mod->con[mod->itrns].concentration;
    }
    
    if(myid == 0){
        printf("\nCONCENTRATION:\n");
    }
    for(int pid = 0; pid < npes; pid++){
        if(myid == pid){
            for (i = 0; i < mod->grid->nnodes; i++) {
                ifmap = mod->fmap[i];
                gid = mod->grid->node[ifmap].gid;
                printf("%i: %.16e",gid,c[i]);
                if(i >= mod->grid->my_nnodes){
                    printf(" - ghost\n");
                }
                else{
                    printf("\n");
                }
            }
        }
#ifdef _MESSG
#ifdef _PETSC
        ierr = MPI_Barrier(PETSC_COMM_WORLD);
#else
        MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    }
}

void print_concentration_file(SSUPER_MODEL *sm, int isuperModel, int myid, int npes){
    SMODEL *mod;
    int i, ierr;
    double *c;
    int ifmap, gid;
    FILE *fp;

    // Open file
    fp = fopen("concentrations.txt","a");

#ifdef _MESSG
#ifdef _PETSC
    ierr = MPI_Barrier(PETSC_COMM_WORLD);
#else
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    mod = &(sm->submodel[isuperModel]);
    
    if (mod->is_sediment_running) {
#ifdef _SEDIMENT
        c = mod->sed->susload[mod->ised].c;
#endif
    } else {
        c = mod->con[mod->itrns].concentration;
    }
    
    if(myid == 0){
        fprintf(fp,"CONCENTRATIONS:\n");
    }
    for(int pid = 0; pid < npes; pid++){
        if(myid == pid){
            for (i = 0; i < mod->grid->my_nnodes; i++) {
                ifmap = mod->fmap[i];
                gid = mod->grid->node[ifmap].gid;
                fprintf(fp,"%4i %.16e\n",gid,c[i]);
            }
            fclose(fp);
        }
#ifdef _MESSG
#ifdef _PETSC
        ierr = MPI_Barrier(PETSC_COMM_WORLD);
#else
        MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    }
}

void print_sw2_file(SSUPER_MODEL *sm, int isuperModel, int myid, int npes){
    SMODEL *mod;
    int i, ierr;
    //double *c;
    int ifmap, gid;
    FILE *fp;

    // Open file
    fp = fopen("sw2_solution.txt","a");

#ifdef _MESSG
#ifdef _PETSC
    ierr = MPI_Barrier(PETSC_COMM_WORLD);
#else
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    mod = &(sm->submodel[isuperModel]);
    
//    if (mod->is_sediment_running) {
//#ifdef _SEDIMENT
//        c = mod->sed->susload[mod->ised].c;
//#endif
//    } else {
//        c = mod->con[mod->itrns].concentration;
//    }
//
    double x_vel, y_vel, head;
    
    if(myid == 0){
        fprintf(fp,"SW2 SOLUTION:\n");
    }
    for(int pid = 0; pid < npes; pid++){
        if(myid == pid){
            for (i = 0; i < mod->grid->my_nnodes; i++) {
                ifmap = mod->fmap[i];
                x_vel = mod->sw->d2->vel[i].x;
                y_vel = mod->sw->d2->vel[i].y;
                head = mod->sw->d2->head[i];
                gid = mod->grid->node[ifmap].gid;
                fprintf(fp,"%4i ",gid);
                fprintf(fp,"%.16e ",x_vel);
                fprintf(fp,"%.16e ",y_vel);
                fprintf(fp,"%.16e\n",head);
            }
            fclose(fp);
        }
#ifdef _MESSG
#ifdef _PETSC
        ierr = MPI_Barrier(PETSC_COMM_WORLD);
#else
        MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    }
}

void print_sw3_file(SSUPER_MODEL *sm, int isuperModel, int myid, int npes){
    SMODEL *mod;
    int i, ierr;
    int ifmap, gid;
    FILE *fp;

    // Open file
    fp = fopen("sw3_solution.txt","a");

#ifdef _MESSG
#ifdef _PETSC
    ierr = MPI_Barrier(PETSC_COMM_WORLD);
#else
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    mod = &(sm->submodel[isuperModel]);
    SSW_3D *sw3 = mod->sw->d3;
    ID_LIST_ITEM *ptr;
    int nd = 0;
    
    double x_vel, y_vel, z_vel, disp;
    
    if(myid == 0){
        fprintf(fp,"SW3 SOLUTION:\n");
    }
    for(int pid = 0; pid < npes; pid++){
        if(myid == pid){
            for (i = 0; i < mod->grid->my_nnodes_sur; i++){
                ptr = mod->grid->vertical_list[i];
                nd = ptr->id;
                
                ifmap = mod->fmap[nd];
                gid = mod->grid->node[ifmap].gid;
                disp = sw3->displacement[nd];
                fprintf(fp,"Surface node\n");
                fprintf(fp,"%4i ",gid);
                fprintf(fp,"%.16e\n",disp);

                while(ptr->next != NULL){
                    nd = ptr->id;
                    ifmap = mod->fmap[nd];
                    gid = mod->grid->node[ifmap].gid;
                    x_vel = sw3->vel[nd].x;
                    y_vel = sw3->vel[nd].y;
                    z_vel = sw3->vel[nd].z;
                    fprintf(fp,"%4i ",gid);
                    fprintf(fp,"%.16e ",x_vel);
                    fprintf(fp,"%.16e ",y_vel);
                    fprintf(fp,"%.16e\n",z_vel);
                    ptr = ptr->next;
                }
            }

            //for (i = 0; i < mod->grid->my_nnodes; i++) {
            //    ifmap = mod->fmap[i];
            //    x_vel = mod->sw->d2->vel[i].x;
            //    y_vel = mod->sw->d2->vel[i].y;
            //    head = mod->sw->d2->head[i];
            //    gid = mod->grid->node[ifmap].gid;
            //    fprintf(fp,"%4i ",gid);
            //    fprintf(fp,"%.16e ",x_vel);
            //    fprintf(fp,"%.16e ",y_vel);
            //    fprintf(fp,"%.16e\n",head);
            //}
            fclose(fp);
        }
#ifdef _MESSG
#ifdef _PETSC
        ierr = MPI_Barrier(PETSC_COMM_WORLD);
#else
        MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    }
}

void print_dw_file(SSUPER_MODEL *sm, int isuperModel, int myid, int npes){
    SMODEL *mod;
    int i, ierr;
    double *dw_head;
    int ifmap, gid;
    FILE *fp;

    // Open file
    fp = fopen("dw_solution.txt","a");

#ifdef _MESSG
#ifdef _PETSC
    ierr = MPI_Barrier(PETSC_COMM_WORLD);
#else
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    mod = &(sm->submodel[isuperModel]);
   
    dw_head = mod->sw->d2->head;
    
    if(myid == 0){
        fprintf(fp,"DW HEAD:\n");
    }
    for(int pid = 0; pid < npes; pid++){
        if(myid == pid){
            for (i = 0; i < mod->grid->my_nnodes; i++) {
                ifmap = mod->fmap[i];
                gid = mod->grid->node[ifmap].gid;
                fprintf(fp,"%4i %.16e\n",gid,dw_head[i]);
            }
            fclose(fp);
        }
#ifdef _MESSG
#ifdef _PETSC
        ierr = MPI_Barrier(PETSC_COMM_WORLD);
#else
        MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    }
}

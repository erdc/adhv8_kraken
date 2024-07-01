/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file fe_newton_tools.c This file  contains all routines for initialization, allocation, etc of the Jacobian and Residual */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
static int FILE_DEBUG = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file initalizes a SuperModel Newton residual
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct
 *  \returns ierr = error flag for initializing residual
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int fe_initialize_supermodel_residual(SSUPER_MODEL *sm) {
    ierr = 0;
#ifdef _PETSC
    ierr = VecZeroEntries(sm->residual);
#else
    sarray_init_dbl(sm->residual, sm->matrix_size);
#endif
    return ierr;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file assembles a local 3 degree-of-freedom element into the global SuperModel Newton residual
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in] SGRID *grid - the grid over which the monolithic residual resides
 *  @param[in] double *fmap - a map from the specific model dof to the supermodel dof
 *  @param[in] int *GnodeIDs - an array of length nnodes with the global node #'s of the local element nodes
 *  @param[in] int nnodes - the number of local nodes on the element
 *  @param[in] DOF_3 *elem_rhs - the local element right-hand-side
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct - adjusts the residual
 * \note elem_rhs[0] = x_eq, elem_rhs[1] = y_eq, elem_rhs[2] = c_eq,
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_elem_assemble(SSUPER_MODEL *sm, SGRID *grid, double *fmap, int *GnodeIDs, int nnodes, DOF_3 *elem_rhs) {
    int i,j;
    
#ifdef _PETSC
    PetscScalar values[3];
    for (i=0; i<nnodes; i++) {
        // No explicitly programmed ghost nodes so only add values from residential nodes
        if(grid->node[GNodeIDs[i]].resident_pe == grid->smpi->myid){
            j = fmap[GNodeIDs[i]];
            if (sm->max_ndof == 1) { // cjt :: why do we still use a 3dof matrix here?
                values[0] = -elem_rhs[i].c_eq;
                values[1] = 0;
                values[2] = 0;
            } else if (sm->max_ndof == 3) {
                values[0] = -elem_rhs[i].x_eq;
                values[1] = -elem_rhs[i].y_eq;
                values[2] = -elem_rhs[i].c_eq;
            } else {
                tl_error("max_ndof should be 1 or 3\n");
            }
            VecSetValuesBlockedLocal(sm->residual,1,&j,values,ADD_VALUES);
        }
    }
#else
    for (i=0; i<nnodes; i++) {
        if (sm->max_ndof == 1) {
            j = 1 * fmap[GNodeIDs[i]];
            sm->residual[j] = -elem_rhs[i].c_eq;
        } else if (sm->max_ndof == 3) {
            j = 3 * fmap[GNodeIDs[i]];
            sm->residual[j]     -= elem_rhs[i].x_eq;
            sm->residual[j + 1] -= elem_rhs[i].y_eq;
            sm->residual[j + 2] -= elem_rhs[i].c_eq;
        }
        else {
            tl_error("max_ndof should be 1 or 3\n");
        }
    }
#endif
    
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file assembles a local 3 degree-of-freedom element into the global SuperModel Newton residual
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in] SGRID *grid - the grid over which the monolithic residual resides
 *  @param[in] double *fmap - a map from the specific model dof to the supermodel dof
 *  @param[in] int *GnodeIDs - an array of length nnodes with the global node #'s of the local element nodes
 *  @param[in] int nnodes - the number of local nodes on the element
 *  @param[in] DOF_3 *elem_rhs - the local element right-hand-side
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct - adjusts the residual
 * \note elem_rhs[0] = x_eq, elem_rhs[1] = y_eq, elem_rhs[2] = c_eq,
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void assemble_residual(SSUPER_MODEL *sm, SGRID *grid, SMAT *mat) {
    int j,k;

    //seems like the easiest way?
    //maybe think about this
    double fmap3d = sm->fmap3d;
    double fmap2d = sm->fmap2d;
    double fmap1d = sm->fmap1d;
    //zero out stuff
    #ifdef _PETSC
        ierr = VecZeroEntries(sm->residual);
    #else
        sarray_init_dbl(sm->residual, sm->ndofs);
    #endif


    double temp[max_elem_dofs]
    //loop through all nelem3d
    for (j=0;j<grid->nelem3d;j++){
        for (k=0;k<sm.nSubMods3d[j];k++){

            //would this work for vector functions?
            //would it be better to put global residual outside inner loop if possible?
            sm.elem3d_physics[j][k]->fe_resid(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

            for (l=0;l<elem3d_physics[j][k].ndof;l++){        
                sm->residual[fmap3d[j,k,l]] +=  temp[l];
            }
        }

    }

    //loop through all nelem2d
    for (j=0;j<grid->nelem2d;j++){
        for (k=0;k<sm.nSubMods2d[j];k++){
            //would this work for vector functions?
            sm->residual[fmap2d[j,k]] += sm.elem2d_physics[j][k]->fe_resid(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

            for (l=0;l<elem3d_physics[j][k].ndof;l++){        
                sm->residual[fmap3d[j,k,l]] +=  temp[l];
            }
        }
    }


    //loop through all nelem1d
    for (j=0;j<grid->nelem1d;j++){
        for (k=0;k<sm.nSubMods1d[j];k++){
            //would this work for vector functions?
            sm->residual[fmap1d[j,k]] += sm.elem1d_physics[j][k]->fe_resid(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

            for (l=0;l<elem3d_physics[j][k].ndof;l++){        
                sm->residual[fmap3d[j,k,l]] +=  temp[l];
            }
        }
    }    
}



/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file assembles a local 3 element into the global SuperModel Newton residual
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in] SGRID *grid - the grid over which the monolithic residual resides
 *  @param[in] double *fmap - a map from the specific model dof to the supermodel dof
 *  @param[in] int *GnodeIDs - an array of length nnodes with the global node #'s of the local element nodes
 *  @param[in] int nnodes - the number of local nodes on the element
 *  @param[in] DOF_3 *elem_rhs - the local element right-hand-side
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct - adjusts the residual
 * \note elem_rhs[0] = x_eq, elem_rhs[1] = y_eq, elem_rhs[2] = c_eq,
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void assemble_matrix(SSUPER_MODEL *sm, SGRID *grid, SMAT *mat) {
    int j,k;

    //seems like the easiest way?
    //maybe think about this
    double fmap3d = sm->fmap3d;
    double fmap2d = sm->fmap2d;
    double fmap1d = sm->fmap1d;
    //zero out stuff
    #ifdef _PETSC
        ierr = VecZeroEntries(sm->residual);
    #else
        sarray_init_dbl(sm->residual, sm->ndofs);
    #endif
    
    //loop through all nelem3d
    for (j=0;j<grid->nelem3d;j++){
        for (k=0;k<sm.nSubMods3d[j];k++){
            //would this work for vector functions?
            //would it be better to put global residual outside inner loop if possible?
            sm->residual[fmap3d[j,k]] += sm.elem3d_physics[j][k].fe_resid(sm,grid,mat);
        }
    }

    //loop through all nelem2d
    for (j=0;j<grid->nelem2d;j++){
        for (k=0;k<sm.nSubMods2d[j];k++){
            //would this work for vector functions?
            sm->residual[fmap2d[j,k]] += sm.elem2d_physics[j][k].fe_resid(sm,grid,mat);
        }
    }


    //loop through all nelem1d
    for (j=0;j<grid->nelem1d;j++){
        for (k=0;k<sm.nSubMods1d[j];k++){
            //would this work for vector functions?
            sm->residual[fmap1d[j,k]] += sm.elem1d_physics[j][k].fe_resid(sm,grid,mat);
        }
    }    
}

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


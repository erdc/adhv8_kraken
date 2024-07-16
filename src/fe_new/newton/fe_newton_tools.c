/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file fe_newton_tools.c This file  contains miscellanious routines used in Newton solve */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
static int FILE_DEBUG = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function initializes the SuperModel Newton solution
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct - adjusts the residual
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void initialize_system(SSUPER_MODEL *sm) {
    int j,k;

    //seems like the easiest way?
    //maybe think about this
    //this will update sm solution structs that depend on time
    //is there a better way so if we add a variable we dont have to add bunches of things
    //everywhere?
    for (i = 0; i < sm->nhead ; i ++){
        sm->head[i] = sm->old_head[i]; 
    }


    for (i = 0; i < sm->nvel2d; i++) {
        sm->vel2d[i].x = sm->old_vel2d[i].x ;
        sm->vel2d[i].y = sm->old_vel2d[i].y;
    }


    for (i = 0; i < sm->nvel3d; i++) {
        sm->vel3d[i].x = sm->old_vel3d[i].x ;
        sm->vel3d[i].y = sm->old_vel3d[i].y;
        sm->vel3d[i].x = sm->old_vel3d[i].z ;
    }

    //need to add loop over number of constituents
    //need some sort of convention
    for (i = 0; i < sm->nconcentration ; i ++){
        sm->concentration[i] = sm->old_concentration[i];
    }

    for (i=0;i>sm->ndisplacement;i++){
        sm->displacement[i] = sm->old_displacement[i];
    }


    for (i=0;i>sm->nprs;i++){
        sm->prs[i] = sm->old_prs[i];
    }    


    for (i=0;i>sm->nc;i++){
        sm->c[i] = sm->old_c[i];
    }    


}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function initializes the SuperModel bc mask for any dirichlet conditions
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct - adjusts the residual
 *  @parma[in] grid (*SGRID) pointer to the grid from a design model
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void initialize_dirichlet_bc(SSUPER_MODEL *sm, SGRID *grid) {

    //Dirichlet condition handling go here
    int istart = 0, isers = 0, istr = 0;
    for (i = 0; i < grid->nnodes; i++) {

        //re-establish bc mask
        set_bc_mask_no(sm.bc_mask[i]);
        
        
        //what about different sting for supermodel?
        //maybe grid->node->string is array of length n sm?
        if (grid->node[i].string > NORMAL){
            istr = grid->node[i].string;

            //check any possible available boundary conditions

            //SW2 Model here, add more later
            //maybe more modular way to do this?
            if (sm->str_vals[istr] == BCT_VEL_DIR) {

                sm->bc_mask[i].u = YES;
                sm->bc_mask[i].v = YES;
                isers = sm->str_values[istr].ol_flow.ivx;
                sm->vel2d[i].x = sseries_get_value(isers, sm->series_head,0);
                isers = sm->str_values[istr].ol_flow.ivy;
                sm->vel2d[i].y = sseries_get_value(isers, sm->series_head,0);
            }
            else if (str_values[istr] == BCT_PRS_DIR) {
                sm->bc_mask[i].head = YES;
                isers = sm->str_values[istr].ol_flow.iu_0;
                sm->head[i] = sseries_get_value(isers, sm->series_head,0);
            }
            else if (str_values[istr] == BCT_VEL_PRS_DIR) {
                sm->bc_mask[i].head = YES;
                sm->bc_mask[i].u = YES;
                sm->bc_mask[i].v = YES;
                isers = sm->str_values[istr].ol_flow.ivx;
                sm->vel2d[i].x = sseries_get_value(isers, sm->series_head,0);
                isers = sm->str_values[istr].ol_flow.ivy;
                sm->vel2d[i].y = sseries_get_value(isers, sm->series_head,0);
                isers = sm->str_values[istr].ol_flow.iu_0;
                sm->head[i] = sseries_get_value(isers, sm->series_head,0);
            }

        }


     }
    
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Updates any model specific non-solution quantities using fe_update routines
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  SSUPER_MODEL *sm pointer to an instant of the SuperModel struct - adjusts the residual
 *  @param[in] SGRID *grid - the grid
 *  @param[in] SMAT *mat - the set of materials
 * \note elem_rhs[0] = x_eq, elem_rhs[1] = y_eq, elem_rhs[2] = c_eq,
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void update_function(SSUPER_MODEL *sm,SGRID *grid,SMAT *mat){
    int j,k;


    //loop through all nelem3d
    for (j=0;j<grid->nelem3d;j++){
        for (k=0;k<sm.nSubMods3d[j];k++){

            //would this work for vector functions?
            //would it be better to put global residual outside inner loop if possible?
            sm.elem3d_physics[j][k]->fe_update(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

        }

    }

    //loop through all nelem2d
    for (j=0;j<grid->nelem2d;j++){
        for (k=0;k<sm.nSubMods2d[j];k++){
            //would this work for vector functions?
            sm.elem2d_physics[j][k]->fe_update(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

        }
    }


    //loop through all nelem1d
    for (j=0;j<grid->nelem1d;j++){
        for (k=0;k<sm.nSubMods1d[j];k++){
            //would this work for vector functions?
            sm.elem1d_physics[j][k]->fe_update(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);

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

void get_residual_norms(SGRID *grid, int nnodes, int *ndof, int macro_ndof
#ifdef _PETSC
        , Vec residual, Vec sol
#else
        , double *residual, double *sol
#endif
        , double *resid_max_norm, double *resid_l2_norm, double *inc_max_norm, int *imax_node, int *iinc_node, int *include_node
        
) {
    
    int i, j, idof, include_node_flag;
    idof=0;
    double lolh = 1.0;
    double partial_l2_norm = 0.;
    double partial_max_norm = 0.;
    double partial_inc_max_norm = 0.;
    int ierr = 0;
    int ctr = 0;
    
#ifdef _MESSG
    int myid = grid->smpi->myid; /* gkc temp warning come back later and check. */
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
    
    for (i = 0; i < nnodes; i++) {

            //WARNING, NEED TO ONLY FIND NDOF that are residential!!!!!
            //need to fix this
            //i is different probably? this only works if ndof is locally owned on nodes 0 to my_nnode
            //check here
            if (grid->node[i].resident_pe == grid->smpi->myid){  
                for(j=0; j<ndof[i],j++){
                    idof+=1;
                    //think this is for multiple grids
                    //include_node_flag = 1;
                    //if (include_node != NULL) {
                    //    if (include_node[i] != YES) include_node_flag = 0;
                    //}
           
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
    }
    *resid_max_norm = partial_max_norm;
    
    // Clean up the read array pointers to PETSc vectors
#ifdef _PETSC
    ierr = VecRestoreArrayRead(sol, &(read_sol));
    ierr = VecRestoreArrayRead(residual, &(read_residual));
#endif
    
    // CJT ::  below should only include include_node == YES, but is shouldn't make a huge difference
#ifdef _MESSG
    *resid_l2_norm = sqrt(messg_dsum(partial_l2_norm, smpi->ADH_COMM) / (macro_ndof));
#else
    *resid_l2_norm = sqrt(messg_dsum(partial_l2_norm) / (ndof));
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


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file newton_tools.c This file  contains miscellanious routines used in Newton solve */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
//static int FILE_DEBUG = OFF;
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
 *  @param[in,out] sm (SMODEL_SUPER*) - pointer to an instant of the SuperModel struct
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void initialize_system(SMODEL_SUPER *sm) {
    int i,j,nmods,inc_index;

    //seems like the easiest way?
    //maybe think about this
    //this will update sm solution structs that depend on time
    //is there a better way so if we add a variable we dont have to add bunches of things
    //everywhere?

    //For now, update sol with old_sol. we can pull out specific things later
    //using dof maps
    for(i=0;i< *(sm->ndofs);i++){
        sm->sol[i] = sm->sol_old[i];
    }
    //then loop through for any initialization routines
    for(i=0;i<sm->nphysics_mat_1d;i++){
        nmods = sm->elem1d_physics_mat[i].nSubmodels;
        for (j=0;j<nmods;j++){
            inc_index = sm->elem1d_physics_mat[i].model[j].fe_init;
            //call wrapper for init function
            //as long as it is not unset
            if (inc_index!=UNSET_INT){
                fe_init[inc_index](sm);
            }
        }
    }
    //then loop through for any initialization routines
    for(i=0;i<sm->nphysics_mat_2d;i++){
        nmods = sm->elem2d_physics_mat[i].nSubmodels;
        for (j=0;j<nmods;j++){
            inc_index = sm->elem2d_physics_mat[i].model[j].fe_init;
            //call wrapper for init function
            //as long as it is not unset
            if (inc_index!=UNSET_INT){
                fe_init[inc_index](sm);
            }
        }
    }
    //then loop through for any initialization routines
    for(i=0;i<sm->nphysics_mat_3d;i++){
        nmods = sm->elem3d_physics_mat[i].nSubmodels;
        for (j=0;j<nmods;j++){
            inc_index = sm->elem3d_physics_mat[i].model[j].fe_init;
            //call wrapper for init function
            //as long as it is not unset
            if (inc_index!=UNSET_INT){
                fe_init[inc_index](sm);
            }
        }
    }
    //for (i = 0; i < sm->nhead ; i ++){
//        sm->head[i] = sm->old_head[i]; 
//    }//
//

//    for (i = 0; i < sm->nvel2d; i++) {
//        sm->vel2d[i].x = sm->old_vel2d[i].x ;
//        sm->vel2d[i].y = sm->old_vel2d[i].y;
//    }//
//

//    for (i = 0; i < sm->nvel3d; i++) {
//        sm->vel3d[i].x = sm->old_vel3d[i].x ;
//        sm->vel3d[i].y = sm->old_vel3d[i].y;
//        sm->vel3d[i].x = sm->old_vel3d[i].z ;
//    }//

//    //need to add loop over number of constituents
//    //need some sort of convention
//    for (i = 0; i < sm->nconcentration ; i ++){
//        sm->concentration[i] = sm->old_concentration[i];
//    }//

//    for (i=0;i>sm->ndisplacement;i++){
//        sm->displacement[i] = sm->old_displacement[i];
//    }//
//

//    for (i=0;i>sm->nprs;i++){
//        sm->prs[i] = sm->old_prs[i];
//    }    


//    for (i=0;i>sm->nc;i++){
//        sm->c[i] = sm->old_c[i];
//    }    


}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function initializes the SMODEL_SUPER bc mask for any dirichlet conditions
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out]  sm (SMODEL_SUPER*) - pointer to an instant of the SuperModel struct - adjusts the residual
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void initialize_dirichlet_bc(SMODEL_SUPER *sm) {

    //Dirichlet condition handling go here
    //maybe there is list of indices you want?
    //need to discuss

    //printf("Initializing Dirichlet B.C. \n");
    int idof;
    int ndofs = *(sm->ndofs);
    
    for(idof = 0; idof<ndofs; idof++){
        //need to use sseries strings to set dirichlet_data
        if(sm->bc_mask[idof] == YES){
            sm->sol[idof] = sm->dirichlet_data[idof];
        }
    }
//    int istart = 0, isers = 0, istr = 0;
//    for (i = 0; i < grid->nnodes; i++) {//

//        //re-establish bc mask
//        set_bc_mask_no(sm.bc_mask[i]);
//        
//        
//        //what about different sting for supermodel?
//        //maybe grid->node->string is array of length n sm?
//        if (grid->node[i].string > NORMAL){
//            istr = grid->node[i].string;//

//            //check any possible available boundary conditions//

//            //SW2 Model here, add more later
//            //maybe more modular way to do this?
//            if (sm->str_vals[istr] == BCT_VEL_DIR) {//

//                sm->bc_mask[i].u = YES;
//                sm->bc_mask[i].v = YES;
//                isers = sm->str_values[istr].ol_flow.ivx;
//                sm->vel2d[i].x = sseries_get_value(isers, sm->series_head,0);
//                isers = sm->str_values[istr].ol_flow.ivy;
//                sm->vel2d[i].y = sseries_get_value(isers, sm->series_head,0);
//            }
//            else if (str_values[istr] == BCT_PRS_DIR) {
//                sm->bc_mask[i].head = YES;
//                isers = sm->str_values[istr].ol_flow.iu_0;
//                sm->head[i] = sseries_get_value(isers, sm->series_head,0);
//            }
//            else if (str_values[istr] == BCT_VEL_PRS_DIR) {
//                sm->bc_mask[i].head = YES;
//                sm->bc_mask[i].u = YES;
//                sm->bc_mask[i].v = YES;
//                isers = sm->str_values[istr].ol_flow.ivx;
//                sm->vel2d[i].x = sseries_get_value(isers, sm->series_head,0);
//                isers = sm->str_values[istr].ol_flow.ivy;
//                sm->vel2d[i].y = sseries_get_value(isers, sm->series_head,0);
//                isers = sm->str_values[istr].ol_flow.iu_0;
//                sm->head[i] = sseries_get_value(isers, sm->series_head,0);
//            }//

//        }//
//

//     }
    
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
 *  @param[in,out] sm (SMODEL_SUPER*) - pointer to an instant of the SMODEL_SUPER struct
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void update_function(SMODEL_SUPER *sm){
    //int j,k;

    //is this just comm update double on sol?
    //need to check if this is necessary
    //printf("calling update function\n");
    return;

    ////loop through all nelem3d
//    for (j=0;j<grid->nelem3d;j++){
//        for (k=0;k<sm.nSubMods3d[j];k++){//

//            //would this work for vector functions?
//            //would it be better to put global residual outside inner loop if possible?
//            sm.elem3d_physics[j][k]->fe_update(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);//

//        }//

//    }//

//    //loop through all nelem2d
//    for (j=0;j<grid->nelem2d;j++){
//        for (k=0;k<sm.nSubMods2d[j];k++){
//            //would this work for vector functions?
//            sm.elem2d_physics[j][k]->fe_update(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);//

//        }
//    }//
//

//    //loop through all nelem1d
//    for (j=0;j<grid->nelem1d;j++){
//        for (k=0;k<sm.nSubMods1d[j];k++){
//            //would this work for vector functions?
//            sm.elem1d_physics[j][k]->fe_update(sm,temp,grid,mat,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);//

//        }
//    }    


}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Increments any model specific solution-dependent quantities after sucessful Newton iteration
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] sm (SMODEL_SUPER*) - pointer to an instant of the SMODEL_SUPER struct
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void increment_function(SMODEL_SUPER *sm){
    int i;


    for (i = 0; i < *(sm->ndofs); i++) {
        sm->sol[i] += sm->lin_sys->dsol[i];
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
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] resid_max_norm (double *) the max (L-infinity) residual error norm for local processor (myid)
 * @param[out] resid_l2_norm  (double *) the l2 residual error norm maxed over all processors
 * @param[out] inc_max_norm   (double *) the incremental max (L-infinity) norm for local processor (myid)
 * @param[out] imax_dof       (int *) the dof at which the residual max occurs
 * @param[out] iinc_dof       (int *) the dof at which the incremental max occurs
 * @param[in]  include_dof    (int *) an array of flags to determine whether the dofs are including in the norm calculations
 * @param[in]  my_ndofs       (int) number of locally owned dof
 * @param[in]  ndofs          (int) number of locally owned dof + number of ghosts
 * @param[in]  macro_ndofs    (int) number of all dofs across all processes
 * @param[in]  residual       (double *) the Newton residual
 * @param[in]  dsol           (double *) the Newton solution (array of solution increments)
 * @param[in]  bc_mask        (int*) array of flags determining what dofs are Dirichlet boundary conditions
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void get_residual_norms(double *resid_max_norm, double *resid_l2_norm, double *inc_max_norm,
         int *imax_dof, int *iinc_dof, int *include_dof,
         int my_ndofs, int ndofs, int macro_ndofs, double *residual, double *dsol, int *bc_mask){
    
    int i, include_dof_flag;
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
    *imax_dof = 0;
    *iinc_dof = 0;
    
    for (i = 0; i < my_ndofs; i++) {       
        include_dof_flag = 1;
        if (include_dof != NULL) {
            if (include_dof[i] != YES) include_dof_flag = 0;
        }
        //dirichlet bc shouldnt be included either
        if (bc_mask !=NULL){
            if (bc_mask[i] == YES) include_dof_flag = 0;
        }   

        // Directly use residual and sol double arrays
        if ( (fabs(dsol[i]) > partial_inc_max_norm) && include_dof_flag == 1) {
            partial_inc_max_norm = fabs(dsol[i]);
            *iinc_dof = i;
        }
        if (solv_isnan(residual[i])){
            printf("\nERROR :: myid %d  || residual[%d] = %f\n", myid, i, residual[i]);
            printf("ERROR :: fe_newton.c :: get_residual_norms :: residual[%d] = %20.10f\n",i,residual[i]);
            ierr = 1;
            exit(-1);
        }
            
        if ((fabs(residual[i]) > partial_max_norm) && include_dof_flag == 1) {
            partial_max_norm = fabs(residual[i]);
            *imax_dof = i;
        }
        
        //Mark added criterion ro include_dof_flag
        //printf("NODE %d, include dof flag %d bc_mask %d\n",i, include_dof_flag,bc_mask[i]);
        if (include_dof_flag == 1){
            //printf("NODE %d INCLUDED\n",i);
            partial_l2_norm += residual[i] * residual[i];
        }
    }
    *resid_max_norm = partial_max_norm;
    
    // CJT ::  below should only include include_node == YES, but is shouldn't make a huge difference
#ifdef _MESSG
    *resid_l2_norm = sqrt(messg_dsum(partial_l2_norm, supersmpi->ADH_COMM) / (macro_ndofs);
#else
    *resid_l2_norm = sqrt(partial_l2_norm)/ndofs;//sqrt(messg_dsum(partial_l2_norm) / (nnodes * nsys));
#endif
    *inc_max_norm = partial_inc_max_norm;
    if((ierr > 0) && (myid == 0)) printf("\n +++++++ WARNING NaN generated!!! ++++++++ \n");
}
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
//void include_node_in_residual(int nnodes, double *head, double *old_head, int *fmap, int *include_node_norm) {
//    int i;
//    sarray_init_value_int(include_node_norm, nnodes, YES);
//    for (i=0; i<nnodes; i++) {
//        if (MIN(head[i],old_head[i]) > 0) {
//            include_node_norm[fmap[i]] = YES;
//        } else {
//            include_node_norm[fmap[i]] = NO;
//        }
//    }
//}
//void print_concentration(SMODEL_SUPER *sm, int isuperModel, int myid, int npes){
//    SMODEL *mod;
//    int i, ierr;
//    double *c;
//    int ifmap, gid;//

//#ifdef _MESSG
//#ifdef _PETSC
//    ierr = MPI_Barrier(PETSC_COMM_WORLD);
//#else
//    MPI_Barrier(MPI_COMM_WORLD);
//#endif
//#endif
//    mod = &(sm->submodel[isuperModel]);
//    
//    if (mod->is_sediment_running) {
//#ifdef _SEDIMENT
//        c = mod->sed->susload[mod->ised].c;
//#endif
//    } else {
//        c = mod->con[mod->itrns].concentration;
//    }
//    
//    if(myid == 0){
//        printf("\nCONCENTRATION:\n");
//    }
//    for(int pid = 0; pid < npes; pid++){
//        if(myid == pid){
//            for (i = 0; i < mod->grid->nnodes; i++) {
//                ifmap = mod->fmap[i];
//                gid = mod->grid->node[ifmap].gid;
//                printf("%i: %.16e",gid,c[i]);
//                if(i >= mod->grid->my_nnodes){
//                    printf(" - ghost\n");
//                }
//                else{
//                    printf("\n");
//                }
//            }
//        }
//#ifdef _MESSG
//#ifdef _PETSC
//        ierr = MPI_Barrier(PETSC_COMM_WORLD);
//#else
//        MPI_Barrier(MPI_COMM_WORLD);
//#endif
//#endif
//    }
//}//
//
//
//

//void print_concentration_file(SMODEL_SUPER *sm, int isuperModel, int myid, int npes){
//    SMODEL *mod;
//    int i, ierr;
//    double *c;
//    int ifmap, gid;
//    FILE *fp;//

//    // Open file
//    fp = fopen("concentrations.txt","a");//

//#ifdef _MESSG
//#ifdef _PETSC
//    ierr = MPI_Barrier(PETSC_COMM_WORLD);
//#else
//    MPI_Barrier(MPI_COMM_WORLD);
//#endif
//#endif
//    mod = &(sm->submodel[isuperModel]);
//    
//    if (mod->is_sediment_running) {
//#ifdef _SEDIMENT
//        c = mod->sed->susload[mod->ised].c;
//#endif
//    } else {
//        c = mod->con[mod->itrns].concentration;
//    }
//    
//    if(myid == 0){
//        fprintf(fp,"CONCENTRATIONS:\n");
//    }
//    for(int pid = 0; pid < npes; pid++){
//        if(myid == pid){
//            for (i = 0; i < mod->grid->my_nnodes; i++) {
//                ifmap = mod->fmap[i];
//                gid = mod->grid->node[ifmap].gid;
//                fprintf(fp,"%4i %.16e\n",gid,c[i]);
//            }
//            fclose(fp);
//        }
//#ifdef _MESSG
//#ifdef _PETSC
//        ierr = MPI_Barrier(PETSC_COMM_WORLD);
//#else
//        MPI_Barrier(MPI_COMM_WORLD);
//#endif
//#endif
//    }
//}//

//void print_sw2_file(SMODEL_SUPER *sm, int isuperModel, int myid, int npes){
//    SMODEL *mod;
//    int i, ierr;
//    //double *c;
//    int ifmap, gid;
//    FILE *fp;//

//    // Open file
//    fp = fopen("sw2_solution.txt","a");//

//#ifdef _MESSG
//#ifdef _PETSC
//    ierr = MPI_Barrier(PETSC_COMM_WORLD);
//#else
//    MPI_Barrier(MPI_COMM_WORLD);
//#endif
//#endif
//    mod = &(sm->submodel[isuperModel]);
//    
////    if (mod->is_sediment_running) {
////#ifdef _SEDIMENT
////        c = mod->sed->susload[mod->ised].c;
////#endif
////    } else {
////        c = mod->con[mod->itrns].concentration;
////    }
////
//    double x_vel, y_vel, head;
//    
//    if(myid == 0){
//        fprintf(fp,"SW2 SOLUTION:\n");
//    }
//    for(int pid = 0; pid < npes; pid++){
//        if(myid == pid){
//            for (i = 0; i < mod->grid->my_nnodes; i++) {
//                ifmap = mod->fmap[i];
//                x_vel = mod->sw->d2->vel[i].x;
//                y_vel = mod->sw->d2->vel[i].y;
//                head = mod->sw->d2->head[i];
//                gid = mod->grid->node[ifmap].gid;
//                fprintf(fp,"%4i ",gid);
//                fprintf(fp,"%.16e ",x_vel);
//                fprintf(fp,"%.16e ",y_vel);
//                fprintf(fp,"%.16e\n",head);
//            }
//            fclose(fp);
//        }
//#ifdef _MESSG
//#ifdef _PETSC
//        ierr = MPI_Barrier(PETSC_COMM_WORLD);
//#else
//        MPI_Barrier(MPI_COMM_WORLD);
//#endif
//#endif
//    }
//}//

//void print_sw3_file(SMODEL_SUPER *sm, int isuperModel, int myid, int npes){
//    SMODEL *mod;
//    int i, ierr;
//    int ifmap, gid;
//    FILE *fp;//

//    // Open file
//    fp = fopen("sw3_solution.txt","a");//

//#ifdef _MESSG
//#ifdef _PETSC
//    ierr = MPI_Barrier(PETSC_COMM_WORLD);
//#else
//    MPI_Barrier(MPI_COMM_WORLD);
//#endif
//#endif
//    mod = &(sm->submodel[isuperModel]);
//    SSW_3D *sw3 = mod->sw->d3;
//    ID_LIST_ITEM *ptr;
//    int nd = 0;
//    
//    double x_vel, y_vel, z_vel, disp;
//    
//    if(myid == 0){
//        fprintf(fp,"SW3 SOLUTION:\n");
//    }
//    for(int pid = 0; pid < npes; pid++){
//        if(myid == pid){
//            for (i = 0; i < mod->grid->my_nnodes_sur; i++){
//                ptr = mod->grid->vertical_list[i];
//                nd = ptr->id;
//                
//                ifmap = mod->fmap[nd];
//                gid = mod->grid->node[ifmap].gid;
//                disp = sw3->displacement[nd];
//                fprintf(fp,"Surface node\n");
//                fprintf(fp,"%4i ",gid);
//                fprintf(fp,"%.16e\n",disp);//

//                while(ptr->next != NULL){
//                    nd = ptr->id;
//                    ifmap = mod->fmap[nd];
//                    gid = mod->grid->node[ifmap].gid;
//                    x_vel = sw3->vel[nd].x;
//                    y_vel = sw3->vel[nd].y;
//                    z_vel = sw3->vel[nd].z;
//                    fprintf(fp,"%4i ",gid);
//                    fprintf(fp,"%.16e ",x_vel);
//                    fprintf(fp,"%.16e ",y_vel);
//                    fprintf(fp,"%.16e\n",z_vel);
//                    ptr = ptr->next;
//                }
//            }//

//            //for (i = 0; i < mod->grid->my_nnodes; i++) {
//            //    ifmap = mod->fmap[i];
//            //    x_vel = mod->sw->d2->vel[i].x;
//            //    y_vel = mod->sw->d2->vel[i].y;
//            //    head = mod->sw->d2->head[i];
//            //    gid = mod->grid->node[ifmap].gid;
//            //    fprintf(fp,"%4i ",gid);
//            //    fprintf(fp,"%.16e ",x_vel);
//            //    fprintf(fp,"%.16e ",y_vel);
//            //    fprintf(fp,"%.16e\n",head);
//            //}
//            fclose(fp);
//        }
//#ifdef _MESSG
//#ifdef _PETSC
//        ierr = MPI_Barrier(PETSC_COMM_WORLD);
//#else
//        MPI_Barrier(MPI_COMM_WORLD);
//#endif
//#endif
//    }
//}//

//void print_dw_file(SMODEL_SUPER *sm, int isuperModel, int myid, int npes){
//    SMODEL *mod;
//    int i, ierr;
//    double *dw_head;
//    int ifmap, gid;
//    FILE *fp;//

//    // Open file
//    fp = fopen("dw_solution.txt","a");//

//#ifdef _MESSG
//#ifdef _PETSC
//    ierr = MPI_Barrier(PETSC_COMM_WORLD);
//#else
//    MPI_Barrier(MPI_COMM_WORLD);
//#endif
//#endif
//    mod = &(sm->submodel[isuperModel]);
//   
//    dw_head = mod->sw->d2->head;
//    
//    if(myid == 0){
//        fprintf(fp,"DW HEAD:\n");
//    }
//    for(int pid = 0; pid < npes; pid++){
//        if(myid == pid){
//            for (i = 0; i < mod->grid->my_nnodes; i++) {
//                ifmap = mod->fmap[i];
//                gid = mod->grid->node[ifmap].gid;
//                fprintf(fp,"%4i %.16e\n",gid,dw_head[i]);
//            }
//            fclose(fp);
//        }
//#ifdef _MESSG
//#ifdef _PETSC
//        ierr = MPI_Barrier(PETSC_COMM_WORLD);
//#else
//        MPI_Barrier(MPI_COMM_WORLD);
//#endif
//#endif
//    }
//}


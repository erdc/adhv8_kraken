/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file fe_newton_build.c This file  contains the routines for building an overlapping  SuperModel Jacobian and residual     */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
static int FILE_DEBUG = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file intializes a SuperModel's submodels
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] SSUPER_MODEL *sm :: ONE  Super Model
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_supermodel_init(SSUPER_MODEL *sm) {
    
    int isubmodel, inode;
    SMODEL *mod = NULL;
    
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
        if(sm->submodel[isubmodel].proc_flag==1){
            mod = &(sm->submodel[isubmodel]); // alias
            
            // Diffusive Wave
            if (mod->flag.DIFFUSIVE_WAVE){
                mod->nsys = 1;
                mod->nsys_sq = 1;
                fe_diffusive(mod);
            }
            
            // Ground water
            if (mod->flag.GW_FLOW){
                mod->nsys = 1;
                mod->nsys_sq = 1;
                fe_gw(mod);
            }
            
            // 2D Shallow Water
            else if (mod->flag.SW2_FLOW){
                mod->nsys = 3;
                mod->nsys_sq = 9;
                if (debug.no_hydro == OFF) {
                    fe_sw2(mod);
                } else {
                    printf("\nAssigning Predetermined hydro values\n");
                    for (inode = 0; inode < sm->grid->nnodes; inode++) {
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
}
    
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file builds a SuperModel Newton residual
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] SSUPER_MODEL *sm :: ONE  Super Model
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_build_supermodel_residual(SSUPER_MODEL *sm) {
    
    int isubmodel;
    
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // zero global residual vector
    int ierr = fe_initialize_supermodel_residual(sm);
    
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Loop over *ONLY* monolithically solved submodels
    for (isubmodel = 0; isubmodel<sm->nsubmodels; isubmodel++) {
        
        // 3D HydroModels
        fe_hvel_resid(sm,isubmodel);
        fe_ns_resid(sm,isubmodel);
        fe_gw_resid(sm,isubmodel);
        
        // 2D HydroModels
        fe_sw2_resid(sm,isubmodel);
        fe_diffusive_resid(sm,isubmodel);
        
        // Transport
        fe_transport_resid(sm,isubmodel);
    
#ifdef _DEBUG
    if (FILE_DEBUG==ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif
    
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file builds a SuperModel  Newton Jacobian
 *  \author    Count Corey J. Trahan
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  \returns   YES for a good calculation and NO for a bad calculation
 *
 * @param[in] SSUPER_MODEL *sm :: ONE  Super Model
 * \note CJT \:: matrix size will grow with refinement. It never shrinks with unrefinement!
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_build_supermodel_matrix(SSUPER_MODEL *sm) {
    
        int isubmodel;
        
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
        // zero global residual vector
        int ierr = fe_initialize_supermodel_matrix(sm);
        
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        // Loop over *ONLY* monolithically solved submodels
        for (isubmodel = 0; isubmodel<sm->nsubmodels; isubmodel++) {
            
            // 3D HydroModels
            fe_hvel_load(sm,isubmodel);
            fe_ns_load(sm,isubmodel);
            fe_gw_load(sm,isubmodel);
            
            // 2D HydroModels
            fe_sw2_load(sm,isubmodel);
            fe_diffusive_load(sm,isubmodel);
            
            // Transport
            fe_transport_load(sm,isubmodel);
        
    #ifdef _DEBUG
        if (FILE_DEBUG==ON) tl_check_all_pickets(__FILE__, __LINE__);
    #endif
        
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
}

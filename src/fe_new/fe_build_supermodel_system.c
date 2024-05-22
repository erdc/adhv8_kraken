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
    
    void fe_build_supermodel_residual(SSUPER_MODEL *sm) {
        
        int isubmodel;
        
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
        // zero global residual vector
        int ierr = fe_initialize_supermodel_residual(sm);
        
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
        STRUCT SDESIGN_MODEL {
            SSUPER_MODEL *sm; // each superModel is flux coupled
            SGRID *grid;
            FLUX_INTERFACE *... // NULL for overlapy coupling (same grid portion)
        }
        
        DESIGN_MODEL.sm[0].intialize(desMod->grid)
        DESIGN_MODEL.sm[1].intialize(desMod->grid)
    
        STRUCT SSUPER_MODEL {
            // one monolithic solve
            SGRID *grid /// passed over from DESIGN MODEL
            MONO_INTERFACE *.... //(already in code) // NULL for overlay coupling (same grid portion)
        }
    
        STRUCT SMODEL {
            //char *physic_flag
            void *fe_body_resid //(for all models)
            void *fe_body_load
        }
        
        // NOTE ::
        // for example, 3D SW, 2D TRNS on shared 2D surface elements
        // sm->smod2d[ielem].fe_resid[0] = fe_sw3_boundary_resid
        // sm->smod2d[ielem].fe_resid[0] = fe_sw2_trns_body_resid
        
        
        // Sw3D, GW
        int nsmod3d = 2 /// another SW3D, GW, 5 Contituents, nsmod3d = 7
        SMODEL smod3d[nelems3d]
        SMODEL smod2d[nelems2d]
        for (ielem=0; ieleme<nelems3d; ieleme++) {
            for (isubmod=0; isubmod<nsmod3d; isubmod++) {
                sm->smod3d[ielem].fe_body_resid[isubmod]
            }
        }
        for (ielem=0; ieleme<nelems2d; ieleme++) {
            for (isubmod=0; isubmod<nsmod2d; isubmod++) {
                sm->smod2d[ielem].fe_body_resid[isubmod]
            }
        }
        
        
        
        
    
        
        // intialization
        if sw2 by reading this 0011 at t=0
            void *fe_sw_body_resid[ie] = fe_sw2_body_resid
            :q
            superModel:
                SGRID *grid;
                SMODEL *mod_sw;
                SMODEL *mod_gw;
                SELEM_PHYSICS residuals3d[sm->grid->nelem3d]
                    fe_sw_body_resid[nelems];
                    fe_sw_body_resid[nelems]; ...
        
    SuperModel
        MOD_SURFACE_WATER *mod_sw;
        
        
    MOD_SURFACE_WATER
            void *
            void *body_resid[nelemns]
            void *body_load[nelemns]
            tau_temporal
            ....
        
        
        
        for (ie=0; ie>sm->grid->nelems3d; ie++) {
    
            
            // Surface water
            sm->mod_sw->body_resid[ie](sm->mod_sw, sm->grid, ie, ...)
            
            
            
            
        }
        
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

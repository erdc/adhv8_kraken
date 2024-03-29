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
            

            //Mark comment: We may want to rethink initializations
            // what exactly goes into the model class now that the
            // grid is independent of it
            // the initialization should definitely build up
            // the dof map

            // we probably want to make physics modules only elemental routines?


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

       //Mark Comments
    //assuming the supermodel has a graph of 3d elem, 2d elem and 1d elem

    //One hydro loop conditional
    //in smodel have one flag for Hydro (SW2,DW or SW3)
    //cut down on if conditions when no hydro is active
    //actually, find specific number of elements per physics
    // and loop this way

    //keep working on this!!
    //model
    int i;
    for (i = 0; i<sm->grid->nelem3dgw; i++){
        //pull data for that element
        //each model will have a list of elements!!!
        ie = sm->grid->elem3dgw[i]
        current_element = sm->grid->elem3d[i];
        global_dofs = current_element->dofs;
        active_physics= current_element->physics ; sm->model->?
        //call any possible marked physics
        //Think about exactly what residual routine needs
        //for now assuming all relevant info is in the element structure

        //Alternative idea, have the nsubmodels as elemental attribute and call those
        if(active_physics == GW_FLOW){
            elem_resid = fe_gw_body_resid(current_element);
        }
        if(active_physics == NS3){
            elem_resid += fe_ns_body_resid(current_element);
        }
        if(active_physics == SW3){
            elem_resid += fe_hvel_body_resid(current_element);
        }

        //load elem_resid into global vector using global_dofs
        sm->residual[global_dofs] = elem_resid

    }
    // now 2d elements
    for (i=0; i<sm->grid->nelem2d; i++){
        //pull data for that element
        current_element = sm->grid->elem2d[i];
        global_dofs = current_element->dofs;
        active_physics = current_element->physics;
        //we need to incorporate logic here of what is ok/not ok combo of physics
        if(active_physics == SW2_flow){
            elem_resid = fe_sw2_body_resid(current_element);
        }
        if(active_physics == DW_FLOW){
            elem_resid = fe_diffusive_body_resid(current_element);
        }
        if(active_physics == TRANSPORT){
            elem_resid = fe_transport_body_resid(current_element);
        }

        //maybe look at using bflag
        //also call any boundary resids from 3d models
        if(active_physics == GW_FLOW){
            elem_resid = fe_gw_boundary_resid(current_element);
        }
        if(active_physics == NS3){
            elem_resid = fe_ns_boundary_resid(current_element);
        }
        if(active_physics == SW3){
            elem_resid = fe_hvel_boundary_resid(current_element);
        }

        //load elem_resid into global vector using global dofs
        sm->residual[global_dofs] = elem_resid

    }
    //now 1d elements, these will only be the elements that are marked
    for (i=0;i<sm->grid->nelem1d;i++){
        //pull data for that element
        current_element = sm->grid->elem2d[i];
        global_dofs = current_element->dofs;
        active_physics = current_element->physics;
        //we need to incorporate logic here of what is ok/not ok combo of physics
        if(active_physics == SW2_flow){
            elem_resid = fe_sw2_boundary_resid(current_element);
        }
        if(active_physics == DW_FLOW){
            elem_resid = fe_diffusive_boundary_resid(current_element);
        }
        if(active_physics == TRANSPORT){
            elem_resid = fe_transport_boundary_resid(current_element);
        }


        //load elem_resid into global vector using global dofs
        global_vec[global_dofs] = elem_resid
    }
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

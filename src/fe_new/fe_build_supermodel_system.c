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
    


    //Mark Comments
    //assuming the supermodel has a graph of 3d elem, 2d elem and 1d elem

    //One hydro loop conditional
    //in smodel have one flag for Hydro (SW2,DW or SW3)
    //cut down on if conditions when no hydro is active

    //option 1 loop over models and then elements within each model
    // this requires no logic but multiple calls to global matrix
    int i;
    //loop through each possibly activated physics
    //for (isubmodel = 0; isubmodel<sm->nsubmodels; isubmodel++) {
    //GW first

    //maybe have a hardcoded order of models?


    //Subsurface models
    ///////////////
    //GW
    isubmodel=0;
    mod = &(sm->submodel[isubmodel]);
    for(i=0; i<sm->grid->nelem3dgw;i++){
        current_element = sm->grid->elem3dgw[i];
        global_dofs = current_element->dofs;
        elem_resid = fe_gw_body_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }

    for (i=0; i<sm->grid->nelem2dgw; i++){
        elem_resid = fe_gw_boundary_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }
    //////////////////

    //Surface dynamics here
    ////////////////
    //SW3
    isubmodel+=1;
    mod = &(sm->submodel[isubmodel]);
    for(i=0; i<sm->grid->nelem3dsw3;i++){
        current_element = sm->grid->elem3dsw3[i];
        global_dofs = current_element->dofs;
        elem_resid = fe_hvel_body_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }
    for (i=0; i<sm->grid->nelem2dsw3; i++){
        elem_resid = fe_hvel_boundary_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }
    //NS
    isubmodel+=1;
    mod = &(sm->submodel[isubmodel]);
    for(i=0; i<sm->grid->nelem3dns;i++){
        current_element = sm->grid->elem3dns[i];
        global_dofs = current_element->dofs;
        elem_resid = fe_ns_body_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }
    for (i=0; i<sm->grid->nelem2dns; i++){
        elem_resid = fe_ns_boundary_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }
    //SW2
    isubmodel+=1;
    mod = &(sm->submodel[isubmodel]);
    for(i=0; i<sm->grid->nelem2dsw2;i++){
        current_element = sm->grid->elem2dsw2[i];
        global_dofs = current_element->dofs;
        elem_resid = fe_sw2_body_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }
    for (i=0; i<sm->grid->nelem1dsw2; i++){
        elem_resid = fe_sw2_boundary_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }
    //DW
    isubmodel+=1;
    mod = &(sm->submodel[isubmodel]);
    for(i=0; i<sm->grid->nelem2ddw;i++){
        current_element = sm->grid->elem2ddw[i];
        global_dofs = current_element->dofs;
        elem_resid = fe_diffusive_body_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }
    for (i=0; i<sm->grid->nelem1ddw; i++){
        elem_resid = fe_diffusive_boundary_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }
    //////////////////

    //Transport dynamics
    ///////////////
    //passive scalar transport
    isubmodel+=1;
    mod = &(sm->submodel[isubmodel]);
    for(i=0; i<sm->grid->nelem2ddw;i++){
        current_element = sm->grid->elem2dtrans[i];
        global_dofs = current_element->dofs;
        elem_resid = fe_transport_body_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }
    for (i=0; i<sm->grid->nelem1ddw; i++){
        elem_resid = fe_transport_boundary_resid(mod);
        sm->residual[global_dofs] = elem_resid
    }

    //any other kind of transport goes here
    ///////////////////////

    ////////////////////
    //Waves

    ////////////////////

    // option 2
    //or is this faster?? this requires conditionals I think 
    //but only 1 assembly at each global element
    // maybe there is way to eliminate logic?

    //Idea:
    //partition elements into types, where they can only be a single kind of physics on to avoid extra logic
    for (i = 0; i<sm->grid->nelem3d; i++){
        //pull data for that element
        //each model will have a list of elements!!!
        ie = sm->grid->elem3d[i]
        current_element = sm->grid->elem3d[i];
        global_dofs = current_element->dofs;
        active_physics= current_element->physics ;
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

#ifdef _DEBUG
    if (FILE_DEBUG==ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif
}
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


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

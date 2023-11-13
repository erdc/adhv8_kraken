/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Calculates the 2D shallow water residual.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to the diffusive wave model struct
 *
 * \details   Solves the following weak, discrete equation: \n
 * \f{eqnarray*}{
 * \weakSWDAcont{e}{i}{h} \\
 * \weakSWMxDD{e}{i}{h} \\
 * \weakSWMxDD{e}{i}{h}
 * \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_sw2_resid(SSUPER_MODEL *sm, int imod) {

    SMODEL *mod = &(sm->submodel[imod]);
    
#ifdef _DEBUG
    assert(mod->flag.SW2_FLOW == ON);
    assert(mod->nsys == 3);
#endif
    
    int DEBUG = OFF;
    int DEBUG_PICKETS = OFF;
#ifdef _DEBUG
    if (mod->t_prev > DEBUG_TIME) debug.residual = ON;
    if (debug.residual == ON) {
        DEBUG = ON;
    }
    time_t time1;  time(&time1);
#endif
    
    int i = 0, j = 0, ie = UNSET_INT, GlobalNodeID = UNSET_INT, node_string = UNSET_INT, nnodes = 0;
    DOF_3 elem_rhs[NDONQUAD]; // allocated to max elemental nodes

#ifdef _PETSC
    // Declare PETSc-specific variables
    PetscInt ierr;
    PetscScalar values[3];
#endif
    
    // alias
    SSW_2D *sw2 = mod->sw->d2;
    SGRID *grid = mod->grid;
    
    /* initializes the residual */
    if (mod->amICoupled == NO) {
#ifdef _PETSC
        ierr = VecZeroEntries(sm->residual);
#else
        sarray_init_dbl(sm->residual, grid->nnodes * mod->nsys);
#endif
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    
    for (ie = 0; ie < grid->nelems2d; ie++) {
        if (grid->wd_flag[ie] == 0 || grid->wd_flag[ie] == 1 || grid->wd_flag[ie] == 2) {
            if (mod->str_values[grid->elem2d[ie].string].phys_flag == SW2_FLAG) {
                
                nnodes = grid->elem2d[ie].nnodes;
                
                fe_sw2_body_resid(mod,elem_rhs,ie,0,i,PERTURB_NONE,0,DEBUG);
#ifdef _DEBUG
                if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
#endif
                // dirchlet enforcement of zero momentum flux at dry nodes
                for (i=0; i<nnodes; i++) {
                    GlobalNodeID = mod->grid->elem2d[ie].nodes[i];
                    if (mod->sw->d2->old_head[GlobalNodeID] <= 0.) {
                        elem_rhs[i].x_eq = 0.0;
                        elem_rhs[i].y_eq = 0.0;
                        sw2->vel[grid->elem2d[ie].nodes[i]].x = 0.;
                        sw2->vel[grid->elem2d[ie].nodes[i]].y = 0.;
                    }
                }
                
                /// assembles global residual
                for (i=0; i<nnodes; i++) {
#ifdef _PETSC
                    GlobalNodeID = mod->grid->elem2d[ie].nodes[i];
                    // No explicitly programmed ghost nodes so only add values from
                    // residential nodes
                    if(grid->node[GlobalNodeID].resident_pe == grid->smpi->myid){
                        j = mod->fmap[GlobalNodeID];
                        values[0] = -elem_rhs[i].x_eq;
                        values[1] = -elem_rhs[i].y_eq;
                        values[2] = -elem_rhs[i].c_eq;
                        VecSetValuesBlockedLocal(sm->residual,1,&j,values,ADD_VALUES);

                        /* Old code commented out below. Keeping temporarily
                         * until current code fully tested. Probably can be
                         * deleted soon - SAM
                         */
                        //// Determine position in global PETSc vector
                        //j = (3 * mod->fmap[GlobalNodeID]) + sm->Istart;
                        ////// Debugging (temporary)
                        ////if((j-sm->Istart) >= (sm->my_nnodes*sm->nsys - 3)){
                        ////    printf("LOCAL: %i  %i  %i  %i\n",j-sm->Istart,GlobalNodeID,mod->fmap[GlobalNodeID],grid->smpi->myid);
                        ////}
                        ////if(j >= (sm->macro_nnodes*sm->nsys - 3)){
                        ////    printf("MACRO: %i>%i  %i  %i  %i\n",j,sm->macro_nnodes*sm->nsys,GlobalNodeID,mod->fmap[GlobalNodeID],grid->smpi->myid);
                        ////}
                        //// TODO: Block this
                        //VecSetValue(sm->residual, j,   (-elem_rhs[i].x_eq), ADD_VALUES);
                        //VecSetValue(sm->residual, j+1, (-elem_rhs[i].y_eq), ADD_VALUES);
                        //VecSetValue(sm->residual, j+2, (-elem_rhs[i].c_eq), ADD_VALUES);
                    }
#else
                    j = 3 * mod->fmap[mod->grid->elem2d[ie].nodes[i]];
                    sm->residual[j]     -= elem_rhs[i].x_eq;
                    sm->residual[j + 1] -= elem_rhs[i].y_eq;
                    sm->residual[j + 2] -= elem_rhs[i].c_eq;
#endif
                }
            }
        }
    }
    
    //printScreen_resid("residual before solving", sm->residual, grid->nnodes, sm->nsys, __LINE__, __FILE__);
    //tl_error("temp");
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 1D ELEMENT CONTRIBUTIONS
    
    for (ie = 0; ie < grid->nelems1d; ie++) {
        if (mod->str_values[ grid->elem1d[ie].string ].phys_flag == SW2_FLAG) {
            
            nnodes = grid->elem1d[ie].nnodes;
            
            fe_sw2_boundary_resid(mod,elem_rhs,ie,0,i,PERTURB_NONE,0,DEBUG);
#ifdef _DEBUG
            if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
#endif
            // dirchlet enforcement of zero momentum flux at dry nodes
            for (i=0; i<nnodes; i++) {
                GlobalNodeID = mod->grid->elem1d[ie].nodes[i];
                if (mod->sw->d2->old_head[GlobalNodeID] <= 0.) {
                    elem_rhs[i].x_eq = 0.0;
                    elem_rhs[i].y_eq = 0.0;
                    sw2->vel[grid->elem1d[ie].nodes[i]].x = 0.;
                    sw2->vel[grid->elem1d[ie].nodes[i]].y = 0.;
                }
            }
            
            /// assembles global residual
            for (i=0; i<nnodes; i++) {
#ifdef _PETSC
                GlobalNodeID = mod->grid->elem1d[ie].nodes[i];
                // No explicitly programmed ghost nodes so only add values from
                // residential nodes
                if(grid->node[GlobalNodeID].resident_pe == grid->smpi->myid){
                    // Determine position in global PETSc vector
                    j = mod->fmap[GlobalNodeID];
                    values[0] = -elem_rhs[i].x_eq;
                    values[1] = -elem_rhs[i].y_eq;
                    values[2] = -elem_rhs[i].c_eq;
                    VecSetValuesBlockedLocal(sm->residual,1,&j,values,ADD_VALUES);
                    
                    //j = (3 * mod->fmap[GlobalNodeID]) + sm->Istart;
                    //VecSetValue(sm->residual, j,   (-elem_rhs[i].x_eq), ADD_VALUES);
                    //VecSetValue(sm->residual, j+1, (-elem_rhs[i].y_eq), ADD_VALUES);
                    //VecSetValue(sm->residual, j+2, (-elem_rhs[i].c_eq), ADD_VALUES);
                }
#else
                j = 3 * mod->fmap[mod->grid->elem1d[ie].nodes[i]];
                sm->residual[j]     -= elem_rhs[i].x_eq;
                sm->residual[j + 1] -= elem_rhs[i].y_eq;
                sm->residual[j + 2] -= elem_rhs[i].c_eq;
#endif
            }
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
   
#ifdef _PETSC
    ierr = VecAssemblyBegin(sm->residual);
    ierr = VecAssemblyEnd(sm->residual);
#else 
    if (mod->amICoupled == NO) {
#ifdef _DEBUG
        if (DEBUG == ON
#ifdef _MESSG
            && grid->smpi->myid == 0
#endif
            ) {
            printScreen_resid("sw2 residual", sm->residual, grid->nnodes, sm->nsys, __LINE__, __FILE__);
            //exit(-1);
            
        }
        time_t time2;  time(&time2);
        TIME_IN_2D_SW_RESID += difftime(time2,time1);
#endif
    }
#endif
}

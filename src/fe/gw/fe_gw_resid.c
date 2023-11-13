/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates the diffusive wave residual for groundwater flow.
 *  @param[in,out]  model a pointer to an AdH model where Newton residual should be updated
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_gw_resid(SSUPER_MODEL *sm, int imod) {
#ifndef _PETSC

    SMODEL *mod = &(sm->submodel[imod]);

#ifdef _DEBUG
    time_t time1;  time(&time1);
    assert(mod->flag.GW_FLOW == ON);
    assert(mod->nsys == 1);
#endif

    int DEBUG = OFF;
    int DEBUG_PICKETS = OFF;
#ifdef _DEBUG
    if (mod->t_prev > DEBUG_TIME) debug.residual = ON;
    if (debug.residual == ON) {
        DEBUG = ON;
    }
#endif

    int ie, i;
    double elem_rhs[MAX_NNODES_ON_ELEM3D];     /* the sm->residual for each equation */
    // aliases
    SGW *gw = mod->sgw;
    SGRID *grid = mod->grid;

    /* initializes the residual */
    if (mod->amICoupled == NO) {
        sarray_init_dbl(sm->residual, grid->nnodes * mod->nsys);
    }

#ifdef _DEBUG
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* loops over the valid elements and loads the element contributions */
    for (ie = 0; ie < grid->nelems3d; ie++) {
        // forms the 3D element residual
        fe_gw_body_resid(mod, elem_rhs, ie, 0, UNSET_INT, UNSET_INT, 1, DEBUG);
#ifdef _DEBUG
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif
        // assembles the elemental residual contributions in global
        for (i = 0; i < grid->elem3d[ie].nnodes; i++) {
            sm->residual[mod->fmap[grid->elem3d[ie].nodes[i]]] -= elem_rhs[i];
        }          
    }

    /**********************************************************************/
    /* loops over the 2D elements and forms their contributions for boundary conditions*/
    for (ie = 0; ie < grid->nelems2d; ie++) {
        fe_gw_boundary_resid(mod, elem_rhs, ie, 0, UNSET_INT, UNSET_INT, 1, DEBUG);
#ifdef _DEBUG
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif
        /* assembles the element contributions */
        for (i = 0; i < grid->elem2d[ie].nnodes; i++) {
            sm->residual[mod->fmap[grid->elem2d[ie].nodes[i]]] -= elem_rhs[i];
        }
    }
    
    //tag(MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // enforce the Dirichlet boundary conditions

    int istart = 0, isers = 0, istr = 0;
    int I = 0;
    for (I = 0; I < mod->grid->nnodes; I++) {

        istart = mod->nsys * I;
        if (mod->grid->node[I].string > NORMAL) {
            istr = mod->grid->node[I].string;
            if (mod->str_values[istr].flow.bc_flag == BCT_PRS_DIR) {
                /*mwf debug
                printf("fe_gw_resid node %d BCT_PRS value= %g \n",i,gw->gw_phead[i]);
                */
                assert (mod->bc_mask[istart] == YES);
                sm->residual[mod->fmap[I]] =0.0;
            } else if (mod->str_values[istr].flow.bc_flag == BCT_DIR) {
                /*mwf debug
                printf("fe_gw_resid node %d BCT_DIR value= %g \n",i,gw->gw_phead[i]);
                */
                assert (mod->bc_mask[istart] == YES);
                sm->residual[mod->fmap[I]] =0.0;

            }
        }
    }
    
    //printf("nelems2d: %d\n",grid->nelems2d);
    //tl_error("test");

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // evaluate the saturations and global constitutive relations now
    sgw_evaluate_element_saturations(mod);
    sgw_evaluate_element_fluxes(mod); /* Gajanan gkc - Do we need this? Someone please check! */

     /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
     /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
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
        TIME_IN_GW_RESID += difftime(time2,time1);
#endif
    }
#endif
}/* end function */

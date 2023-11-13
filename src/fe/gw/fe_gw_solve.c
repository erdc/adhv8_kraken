/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     The FE engine for computing the next time-step of the 3D GW model. Returns
 *            NO if the nonlinear step fails, YES if successful.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \author    Gajanan Choudhary, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to a 3D SW wave model struct
 *
 * \details   This function was created by copying over and modifying fe_ns3_solve.c
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_gw_solve(SSUPER_MODEL *sm, int isubModel) {
  
    int isuperModel = 0;
    SMODEL *mod = &(sm->submodel[isubModel]);

    assert(mod->sgw);
    assert(mod->grid);
    
    // Monolithic Newton solve ++++
    mod->nsys = 1; sm->nsys = 1;
    mod->nsys_sq = 1; sm->nsys_sq = 1;
    sm->solver_info.refresh = YES;
    sm->solver_info.PRN_NEWTON = GW;
    sm->solver_info.LINEAR_PROBLEM = NO;
    mod->solver_info.refresh = YES;
    mod->solver_info.PRN_NEWTON = GW;
    mod->solver_info.LINEAR_PROBLEM = NO;
    
    if (fe_newton(sm, isuperModel, mod->grid->my_nnodes, mod->grid->nnodes, mod->grid->macro_nnodes,
#ifdef _MESSG
                  sm->supersmpi,
#endif
                  fe_gw_init, fe_gw_update, fe_gw_resid, fe_gw_load, fe_gw_inc) == NO) {
        return (NO);
    }
    
    sgw_evaluate_element_fluxes(mod);

    // reset to defaults
    sm->solver_info.PRN_NEWTON = OFF;
    mod->solver_info.PRN_NEWTON = OFF;
    
//    int i;
//    printf("\n");
//    for (i=0; i<mod->grid->nnodes; i++) {
//        if (mod->grid->node[i].gid == 49)  printf("myid: %d gnode: 49 lnode: %d after solve: %f \n",mod->grid->smpi->myid,mod->grid->node[i].id,mod->sgw->gw_phead[49]);
//        if (mod->grid->node[i].gid == 151) printf("myid: %d gnode: 151 lnode: %d after solve: %f \n",mod->grid->smpi->myid,mod->grid->node[i].id,mod->sgw->gw_phead[151]);
//        if (mod->grid->node[i].gid == 252) printf("myid: %d gnode: 252 lnode: %d after solve: %f \n",mod->grid->smpi->myid,mod->grid->node[i].id,mod->sgw->gw_phead[252]);
//        if (mod->grid->node[i].gid == 352) printf("myid: %d gnode: 352 lnode: %d after solve: %f \n",mod->grid->smpi->myid,mod->grid->node[i].id,mod->sgw->gw_phead[352]);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    tl_error("test");
    
    return (YES);
}

/* This routine calls the coupled DW-GW fe_newton solve */
#include "global_header.h"

int fe_dwgw_hybrid_solve(SSUPER_MODEL *sm, int isuperModel) {
    
    assert(isuperModel == 0);
    
    int isubmodel;
    
    /* newton solve (for DW-GW) */
    sm->nsys = 1;
    sm->nsys_sq = 1;
    sm->solver_info.refresh = YES;
    sm->solver_info.PRN_NEWTON = GW;
    sm->solver_info.LINEAR_PROBLEM = NO;
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
        sm->submodel[isubmodel].solver_info.refresh = YES;
        sm->submodel[isubmodel].solver_info.PRN_NEWTON = GW;
        sm->submodel[isubmodel].solver_info.LINEAR_PROBLEM = NO;
    }
    
    if (fe_newton(sm, isuperModel, sm->my_nnodes, sm->nnodes, sm->macro_nnodes,
#ifdef _MESSG
                  sm->supersmpi,
#endif
                  fe_dwgw_hybrid_init, fe_dwgw_hybrid_update, fe_dwgw_hybrid_resid, fe_dwgw_hybrid_load, fe_dwgw_hybrid_inc) == NO) {
        return (NO);
    }
    

    /* reset to defaults*/
    sm->solver_info.PRN_NEWTON = OFF;
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {

        sgw_evaluate_element_fluxes(&(sm->submodel[isubmodel]));

        sm->submodel[isubmodel].solver_info.PRN_NEWTON = OFF;
    }
}

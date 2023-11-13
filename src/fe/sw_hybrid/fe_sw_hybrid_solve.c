#include "global_header.h"

static int skip_wvel_calc = OFF;

int fe_sw_hybrid_solve(SSUPER_MODEL *sm, int isuperModel) {
    
    assert(isuperModel == 0); // CJT :: fe_main only called for one superModel at a time
    
    int isubmodel,inode,block;
    SMODEL *mod;
    
    /* newton solve (for SW 3, solve depth-average continuity and x,y momentum equations) */
    sm->nsys = 3;
    sm->nsys_sq = 9;
    sm->solver_info.refresh = YES;
    sm->solver_info.PRN_NEWTON = HYD; // for now... need HYBRID option
    sm->solver_info.LINEAR_PROBLEM = NO;
    
    if (fe_newton(sm, isuperModel, sm->my_nnodes, sm->nnodes, sm->macro_nnodes,
#ifdef _MESSG
                  sm->supersmpi,
#endif
                  fe_sw_hybrid_init, fe_sw_hybrid_update, fe_sw_hybrid_resid, fe_sw_hybrid_load, fe_sw_hybrid_inc) == NO) {
        return (NO);
    }
    
    
    sm->save_wvel_firstcall = 1; /* relevant only for 2D-3D coupled models */
    int count=0;
    
    // Following few lines check if there is any 3D submodel present in the supermodel :: cjt :: counts the number of *2D MODELS*
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
        if (!(sm->submodel[isubmodel].flag.SW3_FLOW)){
            count++;
        }
    }
    
    if ((count == sm->nsubmodels) || ((sm->supersmpi_wvel==NULL)&&(sm->supersmpi->npes>1)) ){
        skip_wvel_calc = ON;
    }
    else{
        skip_wvel_calc = OFF;
    }
    
    if (skip_wvel_calc == OFF){
        
        sm->nsys = 1;
        sm->nsys_sq = 1;
        sm->solver_info.refresh = YES;
        sm->solver_info.PRN_NEWTON = WVEL; // for now... need HYBRID option
        sm->solver_info.LINEAR_PROBLEM = NO;
        for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
            mod = &(sm->submodel[isubmodel]); // alias
            mod->nsys        = 1;
            mod->nsys_sq     = 1;
        }
        
        /* must loop over all the nodes here, since block can completely change */
        for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++){
            if(sm->submodel[isubmodel].proc_flag==1){
                mod = &(sm->submodel[isubmodel]);
                if (mod->flag.SW3_FLOW){
                    for (inode = 0; inode < mod->grid->nnodes; inode++) {
                        block = mod->grid->node[inode].block;
#ifdef _MESSG
                        sm->solver_info.node_block[mod->fmap_wvel[inode]] = sm->supersmpi_wvel->myid;
#else
                        sm->solver_info.node_block[mod->fmap_wvel[inode]] = block;
#endif
                    }
                }
            }
        }
#ifdef _MESSG
        comm_update_int(sm->solver_info.node_block, 1, sm->supersmpi_wvel);
#endif
        if (fe_newton(sm, isuperModel, sm->wvel_my_nnodes, sm->wvel_nnodes, sm->wvel_macro_nnodes,
#ifdef _MESSG
                      sm->supersmpi_wvel,
#endif
                      fe_sw_hybrid_wvel_init, fe_sw_hybrid_wvel_update, fe_sw_hybrid_wvel_resid, fe_sw_hybrid_wvel_load, fe_sw_hybrid_wvel_inc) == NO) {
            return (NO);
        }
        
        
    } else {
        
        for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
            if(sm->submodel[isubmodel].proc_flag==1){
                mod = &(sm->submodel[isubmodel]);
                if (mod->flag.SW3_FLOW){
                    for (inode=0; inode<mod->grid->nnodes; inode++){
                        mod->sw->d3->vel[inode].z = 0.0;
                    }
                }
            }
        }
    }
    
    /* reset to defaults*/
    sm->solver_info.PRN_NEWTON = OFF;
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
        mod = &(sm->submodel[isubmodel]); // alias
        mod->solver_info.PRN_NEWTON = OFF;
    }
    
    return YES;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     The FE engine for computing the next time-step of the 3D SW model. Returns
 *            NO if the nonlinear step fails, YES if successful.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to a 3D SW wave model struct
 *
 * \details   Solves the following weak, discrete SW 3D equations:
 * \nSTEP 1: - (h,u,v): \n
 * \f{eqnarray*}{
 *   \weakSWDaContReducedKinematic{i} \\
 *   \weakSWMxDDD{e}{i}{h} \\
 *   \weakSWMxDDD{e}{i}{h}
 * \f}
 * \nSTEP2 - Either of the following (w): \n
 * \f{eqnarray*}{
 *    \weakSWContSurDDD{e}{i}{h} \\
 *    \weakSWContBedDDD{e}{i}{h}
 * \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_sw3_solve(SSUPER_MODEL *sm, int isubModel) {
    
    int isuperModel = 0;
    SMODEL *mod = &(sm->submodel[isubModel]);

    assert(mod->sw->d3);
    assert(mod->grid);
    
    // alias
    SSW_3D *sw3d = mod->sw->d3;

    // Newton solve for depth-average continuity and x,y momentum equations
    mod->nsys = 3;      sm->nsys = 3;
    mod->nsys_sq = 9;   sm->nsys_sq = 9;
    sm->solver_info.refresh = YES;
    sm->solver_info.PRN_NEWTON = HVEL;
    sm->solver_info.LINEAR_PROBLEM = NO;
    if (fe_newton(sm, isuperModel, mod->grid->my_nnodes, mod->grid->nnodes, mod->grid->macro_nnodes, 
#ifdef _MESSG
                  sm->supersmpi,
#endif
                  fe_hvel_init, fe_hvel_update, fe_hvel_resid, fe_hvel_load, fe_hvel_inc) == NO) {
        return (NO);
    }

    // Newton solve for full continuity for w back-out
    ///*
    mod->nsys = 1;      sm->nsys = 1;
    mod->nsys_sq = 1;   sm->nsys_sq = 1;
    sm->solver_info.refresh = YES;
    sm->solver_info.PRN_NEWTON = WVEL;
    sm->solver_info.LINEAR_PROBLEM = NO; /// YES??
    if (fe_newton(sm, isuperModel, mod->grid->my_nnodes, mod->grid->nnodes, mod->grid->macro_nnodes, 
#ifdef _MESSG
                  sm->supersmpi,
#endif
                  fe_wvel_init, fe_wvel_update, fe_wvel_resid, fe_wvel_load, fe_wvel_inc) == NO) {
        return (NO);
    }
    //*/
    
    // reset to defaults
    sm->solver_info.PRN_NEWTON = OFF;
    
    return (YES);
}

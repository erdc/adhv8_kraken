/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     The FE engine for computing the next time-step of the 3D NS model. Returns
 *            NO if the nonlinear step fails, YES if successful.
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

int fe_ns3_solve(SSUPER_MODEL *sm, int isubModel) {
  
    int isuperModel = 0;
    SMODEL *mod = &(sm->submodel[isubModel]);

    assert(mod->ns->d3);
    assert(mod->grid);
    
    // alias
    SNS_3D *ns3d = mod->ns->d3;
    
    // Monolithich Newton solve for pressure and x,y,z momentum equations ++++
    mod->nsys = 4; sm->nsys = 4;
    mod->nsys_sq = 16; sm->nsys_sq = 16;
    mod->solver_info.refresh = YES;
    mod->solver_info.PRN_NEWTON = NS;
    mod->solver_info.LINEAR_PROBLEM = NO;
    if (fe_newton(sm, isuperModel, mod->grid->my_nnodes, mod->grid->nnodes, mod->grid->macro_nnodes, 
#ifdef _MESSG
                  sm->supersmpi,
#endif
                  fe_ns3_init, fe_ns3_update, fe_ns3_resid, fe_ns3_load, fe_ns3_inc) == NO) {
        return (NO);
    }
    
    // Split Operator Newton Solve +++++++++++++++++++++++++++++++++++++++++++
    // Newton solve for depth-average continuity and x,y momentum equations
    //    mod->nsys = 3;
    //    mod->nsys_sq = 9;
    //    mod->solver_info.refresh = YES;
    //    mod->solver_info.PRN_NEWTON = NS;
    //    mod->solver_info.LINEAR_PROBLEM = NO;
    //    if (fe_newton(mod, mod->grid, fe_ns_init, fe_ns_update, fe_ns_resid, fe_ns_load, fe_ns_inc) == NO) {
    //        return (NO);
    //    }
    
    // Newton solve for pressure
    //    mod->nsys = 1;
    //    mod->nsys_sq = 1;
    //    mod->solver_info.refresh = YES;
    //    mod->solver_info.PRN_NEWTON = PRS;
    //    mod->solver_info.LINEAR_PROBLEM = NO; /// YES??
    //    if (fe_newton(mod, mod->grid, fe_prs_init, fe_prs_update, fe_prs_resid, fe_prs_load, fe_prs_inc) == NO) {
    //        return (NO);
    //    }
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    // reset to defaults
    mod->solver_info.PRN_NEWTON = OFF;
    
    return (YES);
}

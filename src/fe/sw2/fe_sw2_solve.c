/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     The FE engine for computing the next time-step of the 2D SW model. Returns
 *            NO if the nonlinear step fails, YES if successful.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to a 2D SW wave model struct
 *
 * \details   Solves the following weak, discrete SW 2D equations: \n
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

int fe_sw2_solve(SSUPER_MODEL *sm, int isubModel) {
 
    int isuperModel = 0;
    SMODEL *mod = &(sm->submodel[isubModel]);

    assert(mod->sw->d2);
    assert(mod->grid);
    
    SSW_2D *sw2d = mod->sw->d2;
    
    // Newton solve
    sm->solver_info.refresh = YES;
    sm->solver_info.PRN_NEWTON = HYD;
    sm->solver_info.LINEAR_PROBLEM = NO;
 
    if (fe_newton(sm, isuperModel, mod->grid->my_nnodes, mod->grid->nnodes, mod->grid->macro_nnodes, 
#ifdef _MESSG
                  sm->supersmpi,
#endif
                  fe_sw2_init, fe_sw2_update, fe_sw2_resid, fe_sw2_load, fe_sw2_inc) == NO) {
        return (NO);
    }

    // reset to defaults
    sm->solver_info.PRN_NEWTON = OFF;
    
    return (YES);
}

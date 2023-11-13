/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     The FE engine for computing the next time-step of the Diffusive Wave model. Returns
 *            NO if the nonlinear step fails, YES if successful.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to the diffusive wave model struct
 *
 * \details   Solves the following weak, discrete diffusive wave equation: \n
 * \f{eqnarray*}{
 *  \resid{i}{dw}{} =\bodyTime{\,2d}{e}{\phidd{i}}{\depth{h}} \,-
 *           \bodyConv{\,2d}{e}{\phidd{i}}{ (\depth{h} \, k \, \grad \elev{h}) } \,-
 *           \bodySource{\,2d}{e}{\phidd{i}}{\srcGen{dw}{h}{}} \,+
 *           \bcConv{}{e}{\phidd{i}}{(\depth{h} \, k \, \grad \elev{h})} \,+
 *           \bodySUPG{\supg{i}{dw}{e}}
 * \f}
 * which can also be written as \n
 * \f{eqnarray*}{
 *  \resid{i}{dw}{} =\bodyTime{\,2d}{e}{\phidd{i}}{\depth{h}} \,-
 *           \bodyConv{\,2d}{e}{\phidd{i}}{ (\depth{h} \, \velb{h}} \,-
 *           \bodySource{\,2d}{e}{\phidd{i}}{\srcGen{dw}{h}{}} \,+
 *           \bcConv{}{e}{\phidd{i}}{(\velb{h} \,\depth{h})} \,+
 *           \bodySUPG{\supg{i}{dw}{e}}
 * \f}
 * where \n
 * \f{eqnarray*}{
 * u &= k \, \deriv{\elev{h}}{x} \\
 * v &= k \, \deriv{\elev{h}}{y} \\
 * k &= \frac{c^2 H^{4/3}}{n^2 \, \|\velh\|}
 * \f}
 * where c is the conversion factor and n is mannings coefficient.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_diffusive_solve(SSUPER_MODEL *sm, int isubModel) {
    
    int isuperModel = 0; // only 1 here
    SMODEL *mod = &(sm->submodel[isubModel]);

    SSW_2D *diff = mod->sw->d2; // alias

    /* newton solve */
    mod->nsys = 1;      sm->nsys = 1;
    mod->nsys_sq = 1;   sm->nsys_sq = 1;
    sm->solver_info.refresh = YES;
    sm->solver_info.PRN_NEWTON = DIFF;
    sm->solver_info.LINEAR_PROBLEM = NO;

    if (fe_newton(sm, isuperModel, mod->grid->my_nnodes, mod->grid->nnodes, mod->grid->macro_nnodes, 
#ifdef _MESSG
                  sm->supersmpi,
#endif
                  fe_diffusive_init, fe_diffusive_update, fe_diffusive_resid, fe_diffusive_load, fe_diffusive_inc) == NO)
    {return (NO);}

    /* reset to defaults */
    sm->solver_info.PRN_NEWTON = OFF;
    
    return (YES);
}

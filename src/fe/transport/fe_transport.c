/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     The FE engine for computing the next time-step of the 2D or 3D transport model.
 *            Returns NO if the nonlinear step fails, YES if successful.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to the diffusive wave model struct
 *
 * \details   Solves the following weak, discrete 2D or 3D transport equation: \n
 * \f{eqnarray*}{
 * \weakTrnsDA{2d}{e}{h}{\phidd{i}}{\cjhDA \, \depth{h}}{i}{t}{j} \\
 * \weakTrns{3d}{e}{h}{\phiddd{i}}{\cjh}{i}{t}{j} \\
 * \f}
 * \n where \f$\cjh\f$ is constituent j's concentration, \f$\cjhDA\f$ is depth-averaged constituent j's concentration \n
 * \f$ \depth{h} \f$ is the water depth and \f$ \velc{h}{t} \f$ and \f$ \velcDA{h}{t} \f$ are water velocities and  depth-averaged velocities.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_transport(SSUPER_MODEL *sm, int isubModel) {
    
    int isuperModel = 0; // only 1 superModel here
    SMODEL *mod = &(sm->submodel[isubModel]);

    int DEBUG = OFF;
    
#ifdef _DEBUG
    assert(mod->grid);
    assert(mod->grid->ndim == 2 || mod->grid->ndim == 3);
    if (DEBUG) tl_check_all_pickets(__FILE__, __LINE__);
#endif
    
    int nnodes = mod->grid->nnodes; // alias
    
    /* newton solve */
    mod->nsys = 1; sm->nsys = 1;
    mod->nsys_sq = 1; sm->nsys_sq = 1;
    sm->solver_info.refresh = YES;
    sm->solver_info.LINEAR_PROBLEM = YES;
    
    SCON *con = NULL;
    for (mod->itrns = 0; mod->itrns < mod->ntransport; mod->itrns++) {
        
        con = &mod->con[mod->itrns]; // alias
        
        if (con->type == SAL) {
            sm->solver_info.PRN_NEWTON = SALT;
        } else if (con->type == TMP) {
            sm->solver_info.PRN_NEWTON = TEMP;
        } else {
            sm->solver_info.PRN_NEWTON = TRN;
        }
        
        // prep solutions
        sarray_copy_dbl(con->older_concentration, con->old_concentration, nnodes);
        sarray_copy_dbl(con->old_concentration, con->concentration, nnodes);
        
        // call Newton Solver
        if (fe_newton(sm, isuperModel, mod->grid->my_nnodes, mod->grid->nnodes, mod->grid->macro_nnodes, 
#ifdef _MESSG
                  sm->supersmpi,
#endif
                  fe_transport_init, fe_transport_update, fe_transport_resid, fe_transport_load, fe_transport_inc) == NO) {
            return (NO);
        }
        
    }
    
#ifdef _NSM
    if (mod->flag.NSM == ON) nsmwq(mod);
#endif
    
    /* reset to defaults*/
    sm->solver_info.PRN_NEWTON = OFF;
    
    /* calculate approximate continuity error over elements for adaption */
    SVECT *vel3d = NULL;
    SVECT2D *vel2d = NULL;
    double *dpl, *old_dpl, *older_dpl;
    if (mod->grid->ndim == 2) {
        if (mod->flag.SW_FLOW == ON) {
            vel2d = mod->sw->d2->vel;
        } else if (mod->flag.NS_FLOW == ON) {
            vel2d = mod->ns->d2->vel;
        }
        scon_calculate_elem2d_error(mod->ntransport, mod->con, vel2d, mod->sw->d2->head, mod->grid, mod->mat, mod->dt);
    } else if (mod->grid->ndim == 3) {
        if (mod->flag.SW_FLOW == ON) {
            vel3d = mod->sw->d3->vel;
            dpl = mod->sw->d3->displacement;
            old_dpl = mod->sw->d3->old_displacement;
            older_dpl = mod->sw->d3->older_displacement;
        } else if (mod->flag.NS_FLOW == ON) {
            vel3d = mod->ns->d3->vel;
            dpl = mod->ns->d3->displacement;
            old_dpl = mod->ns->d3->old_displacement;
            older_dpl = mod->ns->d3->older_displacement;
        }
        scon_calculate_elem3d_error(mod->ntransport, mod->con, vel3d, dpl, old_dpl, older_dpl, mod->grid, mod->mat, mod->dt);
    }
    return (YES);
}

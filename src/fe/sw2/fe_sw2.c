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

int fe_sw2(SMODEL *mod) {
    
    assert(mod->sw->d2);
    assert(mod->grid);
    
    int i,j;

    // alias
    SSW_2D *sw2d = mod->sw->d2;
    int nnodes = mod->grid->nnodes;
    
    // prep solutions
    int inode = 0;
    for (inode = 0; inode < nnodes; inode++) {
        sw2d->older_head[inode]  = sw2d->old_head[inode];
        sw2d->older_vel[inode].x = sw2d->old_vel[inode].x;
        sw2d->older_vel[inode].y = sw2d->old_vel[inode].y;
        sw2d->old_head[inode]  = sw2d->head[inode];
        sw2d->old_vel[inode].x = sw2d->vel[inode].x;
        sw2d->old_vel[inode].y = sw2d->vel[inode].y;
    }
    
    // over-ride for testing purposes
    /*
    for (inode = 0; inode < nnodes; inode++) {
        sw2d->older_head[inode]  = 1.1 + inode*0.001;
        sw2d->older_vel[inode].x = 0.1 + inode*0.001;
        sw2d->older_vel[inode].y = 0.5 + inode*0.001;
        sw2d->old_head[inode]  = 2.1   + inode*0.001;
        sw2d->old_vel[inode].x = 0.2   + inode*0.001;
        sw2d->old_vel[inode].y = 0.25  + inode*0.001;
        sw2d->head[inode]  = 0.1       + inode*0.001;
        sw2d->vel[inode].x = 0.02      + inode*0.001;
        sw2d->vel[inode].y = 0.025     + inode*0.001;
    }
    */
    
    //if (mod->flag.CONVEYANCE == YES) fe_sw2_1d_conveyance(mod);
    
    // reallocated SW 2D elemental residuals stored for transport
    if ((mod->flag.ADAPTED_THE_GRID == YES) || (sw2d->elem_rhs_realloc==0)) {
        for (i=0; i<MAX_NNODES_ON_ELEM2D; i++) {
            sw2d->elem_rhs_dacont_extra_terms[i] = (double *) tl_realloc(sizeof(double), mod->grid->nelems2d, mod->grid->nelems2d_old,sw2d->elem_rhs_dacont_extra_terms[i]);
            for (j=0; j<mod->grid->nelems2d; j++) {
                sw2d->elem_rhs_dacont_extra_terms[i][j] = 0.;
            }
        }
        sw2d->elem_rhs_realloc = 1;
    } else {
        for (i=0; i<MAX_NNODES_ON_ELEM2D; i++) {
            for (j=0; j<mod->grid->nelems2d; j++) {
                sw2d->elem_rhs_dacont_extra_terms[i][j] = 0.;
            }
        }
    }
    
    // baroclinic sw2d density
    switch(mod->flag.BAROCLINIC) {
        case 1:
            tl_density_calculator_metric(mod->density, NULL, 1., mod->con[mod->salinity_id].concentration, mod->con[mod->salinity_id].property[0],  nnodes, sw2d->density, mod->flag.EOS);
            break;
        case 10:
            tl_density_calculator_metric(mod->density, mod->con[mod->temperature_id].concentration, 1., NULL, mod->con[mod->temperature_id].property[0], nnodes, sw2d->density, mod->flag.EOS);
            break;
        case 11:
            tl_density_calculator_metric(mod->density, mod->con[mod->temperature_id].concentration, 1., mod->con[mod->salinity_id].concentration, mod->con[mod->salinity_id].property[0], nnodes, sw2d->density, mod->flag.EOS);
            break;
    }
    
    if(mod->flag.BAROCLINIC > 0){
        for (inode = 0; inode < nnodes; inode++) {
            sw2d->density[inode] /= mod->density;
            sw2d->density[inode] -= 1.;
        }
    }
    
#ifdef _ADH_STRUCTURES
    if (mod->nweir > 0)
        weir_flow(mod->grid->node, nnodes, mod->nweir, sw2d->old_head, sw2d->older_vel,
                  mod->str_values, mod->weir, mod->flag, mod->gravity, mod->grid);
    if (mod->nflap > 0)
        flap_flow(mod->grid->node, nnodes, mod->nflap, sw2d->old_head, sw2d->older_vel,
                  mod->str_values, mod->flap, mod->flag, mod->gravity, mod->grid);
    if (mod->nsluice > 0)
        sluice_flow(mod->grid->node, nnodes, mod->nsluice, sw2d->old_head, sw2d->older_vel,
                    mod->str_values, mod->series_head, mod->sluice, mod->flag, mod->gravity, mod->grid);
#endif
    
    return (YES);
}

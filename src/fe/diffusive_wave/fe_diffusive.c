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

int fe_diffusive(SMODEL *mod) {
    
    assert(mod->sw->d2);
    assert(mod->grid);
    assert(mod->flag.DIFFUSIVE_WAVE == ON);

    SSW_2D *diff = mod->sw->d2; // alias

    /* prep solutions */
    int inode = 0;
    for (inode = 0; inode < mod->grid->nnodes; inode++) {
      diff->older_head[inode]  = diff->old_head[inode];
      diff->old_head[inode]    = diff->head[inode];
    }

    // Gajanan gkc: This is to avoid freeing error for elem_rhs_dacont_extra_terms 
    // when using mesh adaptivity in GW-DW coupling (and possibly also in DW-only
    // simulations). This variable is not realloc'd in DW models, but when being
    // freed, the latest nelems2d size is used, which causes freeing error in debug mode.

    int i,j;
    // reallocated SW 2D elemental residuals stored for transport
    if ((mod->flag.ADAPTED_THE_GRID == YES) && (diff->elem_rhs_realloc==0)) {
        for (i=0; i<MAX_NNODES_ON_ELEM2D; i++) {
            diff->elem_rhs_dacont_extra_terms[i] = (double *) tl_realloc(sizeof(double), mod->grid->nelems2d, mod->grid->nelems2d_old, diff->elem_rhs_dacont_extra_terms[i]);
            for (j=0; j<mod->grid->nelems2d; j++) {
                diff->elem_rhs_dacont_extra_terms[i][j] = 0.;
            }
        }
        diff->elem_rhs_realloc = 1;
    } else {
        for (i=0; i<MAX_NNODES_ON_ELEM2D; i++) {
            for (j=0; j<mod->grid->nelems2d; j++) {
                diff->elem_rhs_dacont_extra_terms[i][j] = 0.;
            }
        }
    }

    return (YES);
}

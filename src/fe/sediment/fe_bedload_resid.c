/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Calculates the 2D bed load transport residual.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to the diffusive wave model struct
 *
 * \details   Solves the following weak, discrete diffusive wave equation: \n
 * \f{eqnarray*}{ \weakBedTransport \f}
 * \n where \f$ c^{\,j,\,h}_{bl} \f$ is the bed concentration of constituent j and \f$ \delta^{\,h\,} \f$ is the bed load thickness.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_bedload_resid(SMODEL *mod) {
    
#ifdef _DEBUG
    assert(mod->flag.SEDIMENT == ON);
    assert(mod->nsys == 1);
#endif
    
    int DEBUG = OFF;
    int DEBUG_PICKETS = OFF;
#ifdef _DEBUG
    if (mod->t_prev > DEBUG_TIME) debug.residual = ON;
    if (debug.residual == ON) {
        DEBUG = ON;
    }
#endif

    // alias
    int ised = mod->ised;
    SSW_2D *sw2 = mod->sw->d2;
    SGRID *grid = mod->grid;
    
    double *depth = NULL;
    double *old_depth = NULL;
    double *older_depth = NULL;
    if (mod->flag.SW2_FLOW) {
        depth = mod->sw->d2->head;
        old_depth = mod->sw->d2->old_head;
        older_depth = mod->sw->d2->older_head;
    } else {
        depth = mod->sw->d3->depth;
        old_depth = mod->sw->d3->old_depth;
        older_depth = mod->sw->d3->old_depth;
    }
    
    int i = 0, ie = 0, idof = 0;
    double elem_rhs[NDONQUAD];     // allocate for maximum node 2D element

    /* note :: in 3d the matrix and residual are still dimensioned by nnodes */
    /* sets the number of degrees of freedom */
    int ndof = grid->nnodes * mod->nsys;

    /* initializes the residual */
    sarray_dbl_init(mod->residual, ndof);

#ifdef _DEBUG
    if (DEBUG) {
        
    }
#endif

    /*******************************************************************************/
    /*******************************************************************************/
    /* loops over the valid elements and loads the element contributions */
    for (ie = 0; ie < grid->nelems2d; ie++) {
        
        // if 3d, only loop over bottom faces
        if (mod->str_values[grid->elem2d[ie].string].phys_flag == SW3_FLAG && mod->str_values[grid->elem2d[ie].string].flow.bc_flag != BCT_BED) continue;

        // calculate the element residual
        fe_bedload_2d_elem_resid(mod,elem_rhs,ie,0,i,PERTURB_NONE,0,DEBUG);

#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
#endif
        
        // assembles the element contributions :: this needs work this is not a multi-dof assembly
        for (i = 0; i < grid->elem2d[ie].nnodes; i++) {
            mod->residual[grid->elem2d[ie].nodes[i]] -= elem_rhs[i];
        }
    }

    /*******************************************************************************/
    /*******************************************************************************/
    /* loops over the valid elements and loads the element contributions */
    for (ie = 0; ie < grid->nelems1d; ie++) {
        
        // for now assume all 1d elements in a 3d run are bedload related ...
        
        // calculate the element residual
        fe_bedload_1d_elem_resid(mod,elem_rhs,ie,0,i,PERTURB_NONE,0,DEBUG);
        
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
#endif
        
        // assembles the element contributions :: this needs work this is not a multi-dof assembly
        for (i = 0; i < grid->elem1d[ie].nnodes; i++) {
            mod->residual[grid->elem1d[ie].nodes[i]] -= elem_rhs[i];
        }
    }

    /*******************************************************************************/
    /*******************************************************************************/

    // checks the residual for very small entries
    for (i = 0; i < grid->nnodes; i++) {
        if (fabs(mod->residual[i]) < SMALL) mod->residual[i] = 0.0;
    }

#ifdef _DEBUG
   if (debug.residual == 1
#ifdef _MESSG
       && grid->smpi->myid == 0
#endif
     ) {
        printScreen_resid("bedload residual", mod->residual, grid->nnodes, mod->nsys, __LINE__, __FILE__);
        exit(-1);
    }
#endif

}

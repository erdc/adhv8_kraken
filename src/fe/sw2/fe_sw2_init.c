/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Initializes SW 2D solutions and boundary conditions for a Newton iterate.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Gary Brown, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout]  mod (SMODEL *) pointer to the model struct
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_sw2_init(SSUPER_MODEL *sm, int isubModel) {
    
    SMODEL *mod = &(sm->submodel[isubModel]);
    int DEBUG = OFF;

    // alias
    SSW_2D *sw = mod->sw->d2;

    // copy old solution to the current solution
    int i;
    for (i = 0; i < mod->grid->nnodes; i++) {
        sw->vel[i].x = sw->old_vel[i].x;
        sw->vel[i].y = sw->old_vel[i].y;
        sw->head[i]  = sw->old_head[i];
    }

    /* enforce the Dirichlet boundary conditions */
    int istart = 0, isers = 0, istr = 0, ndof = 3;
    for (i = 0; i < mod->grid->nnodes; i++) {
        
        istart = ndof * i;
        mod->bc_mask[istart] = NO;
        mod->bc_mask[istart + 1] = NO;
        mod->bc_mask[istart + 2] = NO;
        
        if (mod->grid->node[i].string > NORMAL) {
            istr = mod->grid->node[i].string;
            
            if (mod->str_values[istr].ol_flow.bc_flag == BCT_VEL_DIR) {
                mod->bc_mask[istart] = YES;
                mod->bc_mask[istart + 1] = YES;
                isers = mod->str_values[istr].ol_flow.ivx;
                sw->vel[i].x = sseries_get_value(isers, mod->series_head,0);
                isers = mod->str_values[istr].ol_flow.ivy;
                sw->vel[i].y = sseries_get_value(isers, mod->series_head,0);
            }
            else if (mod->str_values[istr].ol_flow.bc_flag == BCT_PRS_DIR) {
                mod->bc_mask[istart + 2] = YES;
                isers = mod->str_values[istr].ol_flow.iu_0;
                sw->head[i] = sseries_get_value(isers, mod->series_head,0);
            }
            else if (mod->str_values[istr].ol_flow.bc_flag == BCT_VEL_PRS_DIR) {
                mod->bc_mask[istart] = YES;
                mod->bc_mask[istart + 1] = YES;
                mod->bc_mask[istart + 2] = YES;
                isers = mod->str_values[istr].ol_flow.ivx;
                sw->vel[i].x = sseries_get_value(isers, mod->series_head,0);
                isers = mod->str_values[istr].ol_flow.ivy;
                sw->vel[i].y = sseries_get_value(isers, mod->series_head,0);
                isers = mod->str_values[istr].ol_flow.iu_0;
                sw->head[i] = sseries_get_value(isers, mod->series_head,0);
            }

        }
        
    }
    
    // get wet/dry flags
    if (mod->proc_flag == 1) {
        fe_sw2_wdflag(mod->sw->d2, mod->grid);
    }
    
#ifdef _DEBUG
    if (DEBUG) {
        tl_check_all_pickets(__FILE__, __LINE__);
    }
#endif
    
}

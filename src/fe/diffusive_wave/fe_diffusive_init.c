/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Initializes diffusive wave solutions and boundary conditions for a Newton iterate.
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

void fe_diffusive_init(SSUPER_MODEL *sm, int imod) {
 
    SMODEL *mod = &(sm->submodel[imod]);

    int DEBUG = OFF;

    int i, ie;                  /* loop counter */
    int istr;                   /* the string number for evaluating boundary conditions */
    int isers;                  /* the series number for evaluating boundary conditions */
    int istart;                 /* sets the starting location in the dof for each node */
    int ndof = 1;               /* sets the number of degrees of freedom per node */

    // alias
    SSW_2D *sw_diff = mod->sw->d2;

    /* copies the old solution to the current solution */
    for (i = 0; i < mod->grid->nnodes; i++) {
        sw_diff->head[i]  = sw_diff->old_head[i];
		istart = ndof * i;
        //mod->bc_mask[istart] = NO;
        //mod->bc_mask[istart + 1] = NO;
        //mod->bc_mask[istart + 2] = NO;
		if (mod->grid->node[i].string > NORMAL) {
            istr = mod->grid->node[i].string;
			if (mod->str_values[istr].pressure.bc_flag == BCT_PRS_DIR) {
                isers = mod->str_values[istr].pressure.iu_0;
                sw_diff->head[i] = sseries_get_value(isers, mod->series_head,0);
            }

		}
    }

#ifdef _DEBUG
    if (DEBUG) {
        tl_check_all_pickets(__FILE__, __LINE__);
    }
#endif
    
}

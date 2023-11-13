/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Initializes SW 3D WVEL solutions and boundary conditions for a Newton iterate.
 * \author    Charlie Berger, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Gary Brown, Ph.D.
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

void fe_wvel_init(SSUPER_MODEL *sm, int imod) {
 
    SMODEL *mod = &(sm->submodel[imod]);

    // aliases
    SGRID *grid = mod->grid;
    SSW_3D *sw3 = mod->sw->d3;
    
    // copies the old solution to the current solution
    int i;
    for(i = 0; i < grid->nnodes; i++) {
        sw3->vel[i].z = sw3->old_vel[i].z;
        mod->bc_mask[i]=NO;
    }
    
#ifdef _MESSG
    comm_update_int(mod->bc_mask, 1, grid->smpi);
#endif
    
    // enforce the Dirichlet boundary conditions
    for(i = 0; i < grid->nnodes; i++) {
        if(grid->node[i].string > NORMAL) {
            if(mod->str_values[grid->node[i].string].flow.bc_flag == BCT_VEL_DIR) {
                mod->bc_mask[i] = YES;
                sw3->vel[i].z = 0.;
            }
        }
    }
#ifdef _MESSG
    comm_update_int(mod->bc_mask, 1, grid->smpi);
#endif
    
}

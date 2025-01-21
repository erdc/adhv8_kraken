/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Assigns a wet/dry flag to all element.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] grid (SGRID *) a grid.  On return, assigned grid->wd_flag[ie]
 * @param[in]     sw2d (SSW_2D *) a shallow water 2D model
 *
 * \note CJT\:: Label wetting and drying AND neighboring elements as wet/dry elements
 * \note CJT\:: flag = 0 :: fully wet || flag = 1 :: some nodes are dry || flag = 2 :: all nodes dry
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//needs to be called inside an init routine
int fe_sw2_init(SMODEL_SUPER *sm) {
    fe_sw2_wdflag_legacy(sm);
    // update Dirchlet condition for nodes with depth less than zero
    //how are mappings going to look??
    //this should work for now
    int i,temp,tempu,tempv;
    //only want to loop over active nodes here not all nodes!!!! Will lead to bug
    for (i=0; i<sm->grid->nnodes; i++) {
        temp = get_cg_dof(PERTURB_H, i, sm->dof_map_local, sm->node_physics_mat);
        if (sm->sol_old[temp] <= 0.) {
            //get the two corresponding u and v entries
            tempu = get_cg_dof(PERTURB_U, i, sm->dof_map_local, sm->node_physics_mat);
            tempv = get_cg_dof(PERTURB_V, i, sm->dof_map_local, sm->node_physics_mat);
            sm->bc_mask[tempu] = YES;
            sm->bc_mask[tempv] = YES;
        }
    }

    
    return 0;
}

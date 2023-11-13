/* creates the pressure at the vertices over each vertical line */
// cjt :: actually returns = (1/g) *(P/rho_o) || So for non-zero surface pressures, the surface node will need to include gravity
// cjt :: g is multiplied later in fe_sw3_hvel_ resids, so must account somehow for surface pressure
#include "global_header.h"

void tl_calculate_pressure(SGRID *grid, SSW_3D *sw3, double density, double perturbation) {
    
    ID_LIST_ITEM *ptr;
    int nd;
    int inode;
    double prs_top, prs_top_plus, prs_top_minus, ele_top, ele_top_plus, ele_top_minus, den_top, perturb_scaled = 0., dum=0.;
    
    /* Loop over the surface nodes */
    for(inode = 0; inode < grid->nnodes_sur; inode++) {
        
        ptr = grid->vertical_list[inode];
        nd = ptr->id;
        
        /* the displacement for all the nodes in the vertical list are identical and the same as the surface */
        ele_top = grid->node[nd].z + sw3->displacement[nd];
        ele_top_plus  = grid->node[nd].z + (sw3->displacement[nd] + sw3->dpl_perturbation[nd]);
        ele_top_minus = grid->node[nd].z + (sw3->displacement[nd] - sw3->dpl_perturbation[nd]);

        den_top = sw3->density[nd];
        prs_top = 0.;
        prs_top_plus = 0.;
        prs_top_minus = 0.;
        
        sw3->prs[nd] = prs_top;
        sw3->prs_plus[nd] = prs_top_plus;
        sw3->prs_minus[nd] = prs_top_minus;
        
        while(ptr->next != NULL) {
            nd = ptr->id;
            //printf("tl_calculate_pressure :: perturbation :: %30.20e\n",sw3->dpl_perturbation[nd]);
            sw3->prs[nd] = prs_top + (ele_top - (sw3->displacement[nd] + grid->node[nd].z)) * 0.5 * (den_top + sw3->density[nd]) / density;
            sw3->prs_plus[nd] = prs_top_plus + (ele_top_plus - ((sw3->displacement[nd] + sw3->dpl_perturbation[nd]) + grid->node[nd].z)) * 0.5 * (den_top + sw3->density[nd]) / density;
            sw3->prs_minus[nd] = prs_top_minus + (ele_top_minus - ((sw3->displacement[nd] - sw3->dpl_perturbation[nd]) + grid->node[nd].z)) * 0.5 * (den_top + sw3->density[nd]) / density;

            ele_top = sw3->displacement[nd] + grid->node[nd].z;
            ele_top_plus  = (sw3->displacement[nd] + sw3->dpl_perturbation[nd]) + grid->node[nd].z;
            ele_top_minus = (sw3->displacement[nd] - sw3->dpl_perturbation[nd]) + grid->node[nd].z;

            den_top = sw3->density[nd];
            prs_top = sw3->prs[nd];
            prs_top_plus = sw3->prs_plus[nd];
            prs_top_minus = sw3->prs_minus[nd];

            ptr = ptr->next;         
   
        }
    }
    
}

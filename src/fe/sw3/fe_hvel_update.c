/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     MPI updates the SW 3D HVEL solutions after a Newton iterate.
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

void fe_hvel_update(SSUPER_MODEL *sm, int imod) {
    
    SMODEL *mod = &(sm->submodel[imod]);
    int DEBUG = OFF;
    
#ifdef _MESSG
    comm_update_double(mod->sw->d3->displacement, 1, mod->grid->smpi);
    comm_update_VECT(mod->sw->d3->vel, mod->grid->smpi);
#endif
    
    /* November 2017 - gkc */
    /* Gajanan adding updates below - This is a proposed bug-fix, particularly required for parallel     */
    /* calculations and correct XDMF output visualization. Otherwise, the 2D plots (depth, d.a.vel, and  */
    /* surface el.) will turn out to be wrong - you will see subdomain boundaries in parallel XDMF vis.  */
    /* which are a result of ghost node dpl and d.a.vel. not getting updated.                            */
    /* Also note: tl_calculate_depthavgvel had to be modified as well for this, to include ghost calcs.  */
    /* Someone pls check whether the following lines duplicated from fe_sw3_hvel_init should be removed. */
    
    
    SGRID *grid = mod->grid;
    SSW_3D *sw3 = mod->sw->d3;
    
    // update subsurface displacement
    tl_vertical_adapt(grid, sw3->displacement);
    
    // update displacement perturbations since surface dpl has changed
    tl_get_dpl_perturbation(grid, mod->perturbation, sw3->displacement, sw3->dpl_perturbation);
    
    // update depth averaged velocities :: cjt should be done, but makes solutions worse for some reason
    tl_calculate_depthavgvel(grid, sw3);
    
    // update pressures at vertices
    tl_calculate_pressure(mod->grid, sw3, mod->density, mod->perturbation);
    
    // update pressures at edge midpoints
    tl_find_edge_mdpt_pressure(mod->grid, sw3, mod->density, mod->perturbation);
    
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<mod->grid->smpi->npes;proc_count++){
            if(proc_count==mod->grid->smpi->myid){
                printf("***********myid %d",mod->grid->smpi->myid);
#endif
                svect2d_printScreen_array("Depth Averaged Velocity", sw3->depth_avg_vel, mod->grid->nnodes_sur, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
        // exit(-1);
#endif
    }
#endif
    
#ifdef _DEBUG
    
    // make sure pressures are (+) or write out debug info and exit run
    int i,j;
    for (i=0; i<grid->nnodes; i++) {
        if (sw3->prs[i] < 0 || sw3->prs_plus[i] < 0 || sw3->prs_minus[i] < 0) DEBUG = ON;
    }
    for (i=0; i<grid->num_midpts; i++) {
        for (j=0; j<5; j++) {
            if (grid->midpt_list[i]->value[j] < 0) DEBUG=ON;
        }
    }
    
    
    if (DEBUG) {
        printScreen_dble_array("hvel_update :: displacement", sw3->displacement, grid->nnodes, __LINE__, __FILE__);   printf("\n");
        printScreen_dble_array("hvel_update :: displacement perturbation", sw3->dpl_perturbation, grid->nnodes, __LINE__, __FILE__); printf("\n");
        svect_printScreen_array("hvel_update",sw3->vel,"velocity",grid->nnodes,__LINE__,__FILE__); printf("\n");
        printScreen_dble_array("hvel_update :: node pressure", sw3->prs, grid->nnodes, __LINE__, __FILE__); printf("\n");
        printScreen_dble_array("hvel_update :: node pressure +", sw3->prs_plus, grid->nnodes, __LINE__, __FILE__); printf("\n");
        printScreen_dble_array("hvel_update :: node pressure -", sw3->prs_minus, grid->nnodes, __LINE__, __FILE__); printf("\n");
        
        printf("hvel_update :: mdpt pressures\n");
        for (i=0; i<grid->num_midpts; i++) {
            for (j=0; j<5; j++) {
                printf("%20.10e",grid->midpt_list[i]->value[j]);
            }
            printf("\n");
        }
        tl_check_all_pickets(__FILE__, __LINE__);
        tl_error("STOP\n");
    }
#endif
}

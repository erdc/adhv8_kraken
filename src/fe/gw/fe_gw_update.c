/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     MPI updates the SW 2D solutions after a Newton iterate.
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

void fe_gw_update(SSUPER_MODEL *sm, int imod) {
    SMODEL *mod = &(sm->submodel[imod]);
    
//    int i;
//    for (i=0; i<mod->grid->nnodes; i++) {
//        if (mod->grid->node[i].gid == 49) printf("Before update: %f \n",mod->sgw->gw_phead[49]);
//    }

    
#ifdef _MESSG
    comm_update_double(mod->sgw->gw_phead, 1, mod->grid->smpi);
#endif
    
//    for (i=0; i<mod->grid->nnodes; i++) {
//        if (mod->grid->node[i].gid == 49) printf("After update: %f \n",mod->sgw->gw_phead[49]);
//    }
    
    //MPI_Barrier(MPI_COMM_WORLD);
    //tl_error("test");
    
    
    //sgw_evaluate_element_saturations(mod);
    //sgw_evaluate_element_fluxes(mod);
}

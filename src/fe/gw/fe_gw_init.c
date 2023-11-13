/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Initializes SW 2D solutions and boundary conditions for a Newton iterate.
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

void fe_gw_init(SSUPER_MODEL *sm, int imod) {

    SMODEL *mod = &(sm->submodel[imod]);

    int DEBUG = OFF;

    // alias
    SGW *gw = mod->sgw;

    // copy predictied solution to the current solution
    int i;
    assert(mod->old_dt > 0.);
    for (i = 0; i < mod->grid->nnodes; i++) {
        /** mwf linear predictor 
        gw->gw_phead[i]  = gw->old_gw_phead[i] + (mod->dt/mod->old_dt)*(gw->old_gw_phead[i]-gw->older_gw_phead[i]);
        */
        gw->gw_phead[i]  = gw->old_gw_phead[i];
        /*gw->gw_phead[i]  = gw->old_gw_phead[i] + 0.95*(mod->dt/mod->old_dt)*(gw->old_gw_phead[i]-gw->older_gw_phead[i]);*/
    }
    // enforce the Dirichlet boundary conditions
    int istart = 0, isers = 0, istr = 0, ndof = 1;
    for (i = 0; i < mod->grid->nnodes; i++) {

        istart = ndof * i;
        mod->bc_mask[istart] = NO;
        if (mod->grid->node[i].string > NORMAL) {
            istr = mod->grid->node[i].string;
            if (mod->str_values[istr].flow.bc_flag == BCT_PRS_DIR) {
                mod->bc_mask[istart] = YES;
                isers = mod->str_values[istr].flow.iu_0;
                gw->gw_phead[i] = sseries_get_value(isers, mod->series_head,0);
                //mwf debug 
#ifdef _DEBUG
                if (DEBUG==ON) printf("fe_gw_init node %d BCT_PRS_DIR value= %g \n",i,gw->gw_phead[i]);
#endif

            } else if (mod->str_values[istr].flow.bc_flag == BCT_DIR) {
                mod->bc_mask[istart] = YES;
                isers = mod->str_values[istr].flow.iu_0;
                gw->gw_phead[i] = sseries_get_value(isers, mod->series_head,0);
                gw->gw_phead[i]-= mod->grid->node[i].z;
#ifdef _DEBUG
                if (DEBUG==ON) printf("fe_gw_init node %d BCT_DIR value= %g \n",i,gw->gw_phead[i]);
#endif
            }
        } 
    }
    
    
//#ifdef _MESSG
//    comm_update_double(gw->predict_gw_phead, 1, mod->grid->smpi);
//#endif
    
    // copy predictied solution to the current solution
    for (i = 0; i < mod->grid->nnodes; i++) {
        gw->predict_gw_phead[i] = gw->gw_phead[i];
    }

//#ifdef _MESSG
//    comm_update_double(gw->gw_phead, 1, mod->grid->smpi);
//#endif

    sgw_evaluate_element_saturations(mod);
    //sgw_evaluate_element_fluxes(mod);
    
//#ifdef _MESSG
//    int ie;
//    for (ie = 0; ie < mod->grid->nelems3d; ie++) {
//        comm_update_double(gw->elem_3d_data[ie].saturation, 1, mod->grid->smpi);
//    }
//#endif

#ifdef _DEBUG
    if (DEBUG) {
        tl_check_all_pickets(__FILE__, __LINE__);
    }
#endif

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     The FE engine for computing the next time-step of the 2D/3D MONOLITHIC transport model.
 *            Returns NO if the nonlinear step fails, YES if successful.
 * \author    Corey Trahan, Ph.D.
 * \author    Gajanan Choudhary, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to the diffusive wave model struct
 *
 * \details   Solves the following weak, discrete 2D or 3D transport equation: \n
 * \f{eqnarray*}{
 * \weakTrnsDA{2d}{e}{h}{\phidd{i}}{\cjhDA \, \depth{h}}{i}{t}{j} \\
 * \weakTrns{3d}{e}{h}{\phiddd{i}}{\cjh}{i}{t}{j} \\
 * \f}
 * \n where \f$\cjh\f$ is constituent j's concentration, \f$\cjhDA\f$ is depth-averaged constituent j's concentration \n
 * \f$ \depth{h} \f$ is the water depth and \f$ \velc{h}{t} \f$ and \f$ \velcDA{h}{t} \f$ are water velocities and  depth-averaged velocities.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_transport_hybrid(SSUPER_MODEL *sm) {
 
    int isuperModel = 0;
    int DEBUG = OFF;
    
    int inode = 0;
    int isubmodel, itrns;
    int nnodes, block;
    int node_count = 0;
    SMODEL *mod;
    
#ifdef _DEBUG
    if (DEBUG) {
        tl_check_all_pickets(__FILE__, __LINE__);
    }
#endif
    
    ////////////////////////////////////////////////////////////////////
    // GAJANAN gkc Combined fe_sw2_transport.c and fe_sw3_transport.c //
    ////////////////////////////////////////////////////////////////////
    
    /* must loop over all the nodes here, since block can completely change */
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++){
        if(sm->submodel[isubmodel].proc_flag>0){
            mod = &(sm->submodel[isubmodel]);
            for (inode = 0; inode < mod->grid->nnodes; inode++) {
                block = mod->grid->node[inode].block;
#ifdef _MESSG
                sm->solver_info.node_block[mod->fmap[inode]] = sm->supersmpi->myid;
#else
                sm->solver_info.node_block[mod->fmap[inode]] = block;
#endif
            }
        }
    }
#ifdef _MESSG
    comm_update_int(sm->solver_info.node_block, 1, sm->supersmpi);
#endif
    
    sm->nsys = 1;
    sm->nsys_sq = 1;
    sm->solver_info.refresh = YES;
    sm->solver_info.LINEAR_PROBLEM = YES;
    
    for (itrns=0; itrns<sm->ntransport; itrns++){
        
        if      (sm->con_type[itrns] == SAL)  sm->solver_info.PRN_NEWTON = SALT;
        else if (sm->con_type[itrns] == TMP)  sm->solver_info.PRN_NEWTON = TEMP;
        else                                  sm->solver_info.PRN_NEWTON = TRN;
        // sm->itrns = itrns;
        
        for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
            if(sm->submodel[isubmodel].proc_flag>0){
                mod = &(sm->submodel[isubmodel]); // alias
                nnodes = mod->grid->nnodes;       /* gkc November 2016, important update! */
                mod->itrns = itrns;
                mod->nsys = 1;
                mod->nsys_sq = 1;
                mod->solver_info = sm->solver_info;
                
                /* prep solutions */
                sarray_copy_dbl(mod->con[itrns].older_concentration, mod->con[itrns].old_concentration, nnodes);
                sarray_copy_dbl(mod->con[itrns].old_concentration, mod->con[itrns].concentration, nnodes);
            }
        }
        
        if (fe_newton(sm, isuperModel, sm->my_nnodes, sm->nnodes, sm->macro_nnodes,
#ifdef _MESSG
                      sm->supersmpi,
#endif
                      fe_transport_hybrid_init, fe_transport_hybrid_update, fe_transport_hybrid_resid,
                      fe_transport_hybrid_load, fe_transport_hybrid_inc) == NO) {
            return (NO);
        }
    }
    
#ifdef _NSM
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
        if(sm->submodel[isubmodel].proc_flag>0){
            mod = &(sm->submodel[isubmodel]); // alias
            if (mod->flag.NSM == ON) {
                nsmwq(mod);
            }
        }
    }
#endif
    
    /* reset to defaults*/
    sm->solver_info.PRN_NEWTON = OFF;
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++) {
        mod = &(sm->submodel[isubmodel]); // alias
        mod->solver_info.PRN_NEWTON = OFF;
    }
    
    // printf("Exiting fe_sw_transport_hybrid with YES\n");
    return (YES);
}

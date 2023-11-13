/* This routine sets the initial iterate for the 3D shallow water calculations (RCB) */
#include "global_header.h"

/***********************************************************/
/***********************************************************/

static int DEBUG = OFF;

/***********************************************************/
/***********************************************************/

void fe_sw_hybrid_init(SSUPER_MODEL *sm, int idum) {
    
    int i, j, k, l;
    int imodel_id;
    
    // aliases
    SMODEL *mod;
    for (i=0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
        if(mod->proc_flag==1){
            if (mod->flag.SW2_FLOW){
                fe_sw2_init(sm,i);
            }
            else if (mod->flag.SW3_FLOW){
                fe_hvel_init(sm,i);
            }
        }
    }
    
    // Set bc_mask = NO at interface nodes if bc_mask = 2.
    for (j=0; j<sm->NumInterfaces; j++){
        
        if ((sm->submodel[sm->interface[j].model_id[0]].flag.SW3_FLOW && sm->submodel[sm->interface[j].model_id[1]].flag.SW3_FLOW) /* i.e. 3D-3D coupling */ ||
            (sm->submodel[sm->interface[j].model_id[0]].flag.SW2_FLOW && sm->submodel[sm->interface[j].model_id[1]].flag.SW2_FLOW) /* i.e. 2D-2D coupling */ ){
            continue; /* Do nothing for this 2D-2D or 3D-3D interface */
        }
        
        // What follows is only executed for 2D-3D coupling.
        for (imodel_id=0; imodel_id<2; imodel_id++){
            mod = &(sm->submodel[sm->interface[j].model_id[imodel_id]]);
            if(mod->proc_flag==1){
                for (k=0; k<sm->interface[j].NumNodeColumns; k++){
                    if (sm->interface[j].nodelist[k].size[imodel_id] > 0){
                        if (mod->bc_mask[sm->nsys * sm->interface[j].nodelist[k].couplednodes[imodel_id][0]] == 2){
                            for (l=0; l<sm->interface[j].nodelist[k].size[imodel_id]; l++){
                                mod->bc_mask[sm->nsys * sm->interface[j].nodelist[k].couplednodes[imodel_id][l]] = NO;
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    /* Create a supermodel bc_mask */
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            mod = &(sm->submodel[i]);
            int j_sub, j_sup;
            for (j=0; j<mod->grid->nnodes; j++){
                j_sub = mod->nsys * j;
                j_sup = mod->nsys * mod->fmap[j];
                sm->bc_mask[j_sup]   = mod->bc_mask[j_sub];
                sm->bc_mask[j_sup+1] = mod->bc_mask[j_sub+1];
                sm->bc_mask[j_sup+2] = mod->bc_mask[j_sub+2];
            }
        }
    }
    
#ifdef _MESSG
    comm_update_int(sm->bc_mask, sm->nsys, sm->supersmpi);
#endif
    
    
    
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                printf("***********myid %d",sm->supersmpi->myid);
#endif
                printScreen_int_array("bc_mask", sm->bc_mask, sm->nnodes*sm->nsys, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
        //exit(-1);
#endif
    }
#endif
}

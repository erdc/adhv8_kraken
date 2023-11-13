
#include "global_header.h"

/***********************************************************/
/***********************************************************/
static int DEBUG = OFF;
/***********************************************************/
/***********************************************************/

void fe_sw_hybrid_wvel_init(SSUPER_MODEL *sm, int idum) {
    int i, j;
    SMODEL *mod;
    
    for (i=0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
        if (mod->flag.SW3_FLOW){
            fe_wvel_init(sm,i);
        }
    }
    
    /* Create a supermodel bc_mask */
    for (i=0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
        if (mod->flag.SW3_FLOW){
            int j_sub, j_sup;
            for (j=0; j<mod->grid->nnodes; j++){
                j_sub = mod->nsys * j;
                j_sup = mod->nsys * mod->fmap_wvel[j];
                sm->bc_mask[j_sup]   = mod->bc_mask[j_sub];
                sm->bc_mask[j_sup+1] = mod->bc_mask[j_sub+1];
                sm->bc_mask[j_sup+2] = mod->bc_mask[j_sub+2];
            }
        }
    }
#ifdef _MESSG
    comm_update_int(sm->bc_mask, sm->nsys, sm->supersmpi_wvel);
#endif
    
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi_wvel->npes;proc_count++){
            if(proc_count==sm->supersmpi_wvel->myid){
                printf("***********myid %d",sm->supersmpi_wvel->myid);
#endif
                printScreen_int_array("bc_mask", sm->bc_mask, sm->wvel_nnodes*sm->nsys, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
        //exit(-1);
#endif
    }
#endif
}

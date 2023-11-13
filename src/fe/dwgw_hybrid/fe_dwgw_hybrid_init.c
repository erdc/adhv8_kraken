/* This routine sets the initial iterate for DW GW coupling */
#include "global_header.h"

/***********************************************************/
/***********************************************************/

void fe_dwgw_hybrid_init(SSUPER_MODEL *sm, int idum) {

    int DEBUG_LOCAL = OFF;

    assert(sm->submodel[0].flag.SW2_FLOW==ON);
    assert(sm->submodel[0].flag.DIFFUSIVE_WAVE==ON);
    
    /* Process the groundwater part */
    fe_gw_init(sm,0);
    
    /* Process the diffusive wave part */
    fe_diffusive_init(sm,0);
    
#ifdef _DEBUG
    if (DEBUG_LOCAL) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                printf("***********myid %d at line %i in file %s",sm->supersmpi->myid, __LINE__, __FILE__);
#endif
                //printScreen_int_array("bc_mask", sm->bc_mask, sm->nnodes*sm->nsys, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
        //exit(-1);
#endif
    }
#endif
}

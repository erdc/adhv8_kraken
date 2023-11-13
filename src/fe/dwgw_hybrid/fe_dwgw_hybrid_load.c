/* This routine calculates the Jacobian of the DW-GW coupled matrix */
#include "global_header.h"

/***********************************************************/
/***********************************************************/
/***********************************************************/
void fe_dwgw_hybrid_load(SSUPER_MODEL *sm, int idum) {
    
    int DEBUG_MATRIX = OFF;
    if (debug.matrix_hvel == ON) {
        DEBUG_MATRIX = ON;
    }
    
    assert(sm->nsys == 1);
    assert(sm->nsys_sq == 1);

    // initialize matrix
    init_adh_matrix(sm->total_nnodes, sm->max_nsys_sq, sm->matrix, sm->diagonal);

    /* Process the groundwater part */
    fe_gw_load(sm,0);

#ifdef _DEBUG
    if (DEBUG_MATRIX) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                printf("***********myid %d",sm->supersmpi->myid);
#endif
                printScreen_matrix("Hybrid DW-GW matrix", sm->diagonal, sm->matrix, sm->nnodes, sm->nsys_sq, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
#endif
        //exit(-1);
    }
#endif

    /* Process the diffusive wave part */
    sm->submodel[0].grid->ndim=2;    /* Hacked to 2 */
    fe_diffusive_load(sm,0);
    sm->submodel[0].grid->ndim=3;    /* Reset to 3 */

#ifdef _DEBUG
    if (DEBUG_MATRIX) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                printf("***********myid %d",sm->supersmpi->myid);
#endif
                printScreen_matrix("Hybrid DW-GW matrix", sm->diagonal, sm->matrix, sm->nnodes, sm->nsys_sq, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
#endif
        //exit(-1);
    }
#endif
}

/************************************************************************************/
/************************************************************************************/
/************************************************************************************/

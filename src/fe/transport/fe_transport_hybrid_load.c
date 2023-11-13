/* This routine loops over the submodels to calculate the Jacobian for the transport problem */
#include "global_header.h"

static int DEBUG_MATRIX = OFF;
/***********************************************************/
/***********************************************************/
/***********************************************************/
void fe_transport_hybrid_load(SSUPER_MODEL *sm, int dum) {
    
    int i, j, iend;
    
    if (debug.matrix_hvel == ON) {
        DEBUG_MATRIX = ON;
    }
    
    /* initializes the diagonal arrays */
    for (i = 0, iend = sm->total_nnodes * sm->max_nsys_sq; i < iend; i++)
        sm->diagonal[i] = 0.0;
    
    /* intialize matrix */
    for (i = 0; i < sm->total_nnodes; i++)
        spv_init(sm->matrix[i], sm->max_nsys_sq);
    for (i = 0; i < sm->total_nnodes; i++){
        spv_reset(sm->matrix+i);
        sm->matrix[i].size = 0;
    }
    
    /* Loop over models */
    for (i = 0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            fe_transport_load(sm,i);
        }
    }
    
#ifdef _MESSG
    //comm_update_double(sm->diagonal, sm->nsys_sq, sm->supersmpi);
#endif
    
    for(i=0;i<sm->nnodes;i++){ /* sm->nsys = 1 here! */
        if(sm->diagonal[i] < NOT_QUITE_SMALL)
            sm->diagonal[i] = NOT_QUITE_SMALL;
    }
    
#ifdef _DEBUG
    if (DEBUG_MATRIX) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                printf("***********myid %d",sm->supersmpi->myid);
#endif
                printScreen_matrix("Hybrid TRANSPORT matrix", sm->diagonal, sm->matrix, sm->nnodes, sm->max_nsys_sq, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
        //exit(-1);
#endif
    }
#endif
}

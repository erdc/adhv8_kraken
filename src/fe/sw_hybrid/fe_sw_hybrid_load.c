/* This routine loops over the submodels to calculate the Jacobian of the shallow water supermodel */
#include "global_header.h"

static int DEBUG_MATRIX = OFF;
/***********************************************************/
/***********************************************************/
/***********************************************************/
void fe_sw_hybrid_load(SSUPER_MODEL *sm, int idum) {
    
    int i, j, iend;
    //int node_count = 0;
    
    // aliases
    SMODEL *mod;
    
    if (debug.matrix_hvel == ON) {
        DEBUG_MATRIX = ON;
    }
    
    assert(sm->nsys == 3);
    assert(sm->nsys_sq == 9);

    // initialize matrix
    init_adh_matrix(sm->total_nnodes, sm->max_nsys_sq, sm->matrix, sm->diagonal);

    // printf("Entering individual fe_sw_loads\n");
    // GAJANAN gkc: loop over all models
    for (i = 0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
        if(mod->proc_flag==1){
            if(mod->flag.SW2_FLOW){
                // printf("Entering fe_sw2_load()\n");
                fe_sw2_load(sm,i);
                // printScreen_matrix("sw hybrid mod->matrix", sm->diagonal, sm->matrix, sm->nnodes, sm->max_nsys_sq, __LINE__, __FILE__);
                // exit(-1);
            }
            else if(mod->flag.SW3_FLOW){
                // printf("Entering fe_sw3_hvel_load()\n");
                fe_hvel_load(sm,i);
                // printScreen_matrix("sw hybrid mod->matrix", sm->diagonal, sm->matrix, sm->nnodes, sm->max_nsys_sq, __LINE__, __FILE__);
                // exit(-1);
            }
        }
        mod = NULL;
    }

    
#ifdef _MESSG
    //    comm_update_double(sm->diagonal, sm->nsys_sq, sm->supersmpi);
#endif
    
    int maxNonzeroCols=0;
    for (i = 0; i < sm->nnodes; i++) {
        for (j=(i*sm->nsys_sq); j<(i+1)*sm->nsys_sq; j+=(sm->nsys+1)){
            if (fabs(sm->diagonal[j]) < NOT_QUITE_SMALL) {
                sm->diagonal[j] = NOT_QUITE_SMALL;
            }
        }
        maxNonzeroCols=MAX(maxNonzeroCols, sm->matrix[i].size);
    }
    
    if (maxNonzeroCols>SPV_BLOCK){
        printf("\ninclude/define.h: SPV_BLOCK = %i", SPV_BLOCK);
        printf("\nMaximum number of non-zero blocks in HVEL matrix = %i", maxNonzeroCols);
        printf("\nRecompile AdH with a higher value of SPV_BLOCK (at least %i) in include/define.h.", maxNonzeroCols);
        exit(-1);
    }
    
#ifdef _DEBUG
    if (DEBUG_MATRIX) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                printf("***********myid %d",sm->supersmpi->myid);
#endif
                printScreen_matrix("Hybrid HVEL matrix", sm->diagonal, sm->matrix, sm->nnodes, sm->nsys_sq, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
        //exit(-1);
#endif
    }
#endif
}

/************************************************************************************/
/************************************************************************************/
/************************************************************************************/
